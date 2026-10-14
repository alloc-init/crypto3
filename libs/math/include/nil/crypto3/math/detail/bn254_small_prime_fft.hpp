//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_MATH_DETAIL_BN254_SMALL_PRIME_FFT_HPP
#define CRYPTO3_MATH_DETAIL_BN254_SMALL_PRIME_FFT_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/math/detail/small_prime_fft_arithmetic.hpp>
#include <nil/crypto3/math/detail/small_prime_fft_params.hpp>
#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>

namespace nil::crypto3::math::detail {

    /**
     * Convert canonical BN254 Fp12 coordinates to a prefix of the small-prime catalogue.
     * Coordinate order is Crypto3's [1,u,v,uv,v^2,uv^2,w,uw,vw,uvw,v^2w,uv^2w].
     * Each output cell is an exact integer in [0,q), congruent to the canonical base-field
     * coordinate modulo q. In particular, a base-field -1 is lifted to p-1 before reduction.
     *
     * Precomputed powers are reusable across operands and sizes. Conversion extracts each
     * coordinate once for all primes, using fixed-size scratch and caller-owned output storage.
     * Kernel selection occurs once per call; the calling thread's rounding mode must remain
     * unchanged until it returns. This conversion has no constant-time guarantee.
     */
    class bn254_small_prime_fft_converter {
        static constexpr std::size_t chunk_bits = 52;
        static constexpr std::size_t chunk_count = 5;
        static constexpr std::size_t block_size = 32;
        using wide_type = boost::multiprecision::uint128_modular_t;

        struct prime_context {
            small_prime_fft_arithmetic arithmetic;
            std::array<double, chunk_count - 1> powers;

            explicit prime_context(std::uint64_t modulus) : arithmetic(modulus) {
                std::uint64_t power = 1;
                for (auto &weight : powers) {
                    power = ((wide_type(power) << chunk_bits) % modulus).convert_to<std::uint64_t>();
                    // Center the weights so that |weight| < q/2. These conversions are exact
                    // under every rounding mode, and setup products fit in 102 bits.
                    const auto centered = static_cast<std::int64_t>(power) -
                                          (power > modulus / 2 ? static_cast<std::int64_t>(modulus) : 0);
                    weight = static_cast<double>(centered);
                }
            }
        };

    public:
        using field_type = algebra::fields::fp12_2over3over2<algebra::fields::alt_bn128_base_field<254>>;
        using value_type = field_type::value_type;

        explicit bn254_small_prime_fft_converter(std::size_t prime_count) {
            if (prime_count == 0 || prime_count > small_prime_fft_primes.size()) {
                throw std::invalid_argument("bn254_small_prime_fft_converter: invalid catalogue prefix");
            }
            primes_.reserve(prime_count);
            for (std::size_t i = 0; i < prime_count; ++i) {
                primes_.emplace_back(small_prime_fft_primes[i].modulus);
            }
        }

        std::size_t prime_count() const {
            return primes_.size();
        }

        /**
         * Write output[(prime*12 + coordinate)*stride + coefficient], padding each channel
         * with zeros from coefficients.size() to stride. Pass a subspan to convert a prefix.
         * stride need not be a power of two; the transform planner supplies that restriction.
         *
         * Require coefficients.size() <= stride and output.size() == prime_count()*12*stride.
         * Check dimensions without overflow before writing any output; throw invalid_argument
         * on mismatch. Empty input writes zeros, including when reusing a previous output buffer.
         */
        void convert(std::span<double> output, std::span<const value_type> coefficients, std::size_t stride) const {
            const std::size_t channels = field_type::arity * primes_.size();
            if (channels == 0 || coefficients.size() > stride || stride > output.size() / channels ||
                output.size() != channels * stride) {
                throw std::invalid_argument("bn254_small_prime_fft_converter: invalid buffer dimensions");
            }
            if (output.empty()) {
                return;
            }
            const bool use_fma = std::all_of(primes_.begin(), primes_.end(), [](const prime_context &prime) {
                return prime.arithmetic.fma_supported();
            });
            if (use_fma) {
                convert_impl<true>(output, coefficients, stride);
            } else {
                convert_impl<false>(output, coefficients, stride);
            }
        }

    private:
        template<bool UseFma>
        void convert_impl(std::span<double> output, std::span<const value_type> coefficients,
                          std::size_t stride) const {
            for (std::size_t channel = 0; channel < field_type::arity * primes_.size(); ++channel) {
                auto padding = output.subspan(channel * stride + coefficients.size(), stride - coefficients.size());
                std::fill(padding.begin(), padding.end(), 0.0);
            }

            std::array<std::array<double, chunk_count>, block_size> chunks;
            for (std::size_t coordinate = 0; coordinate < field_type::arity; ++coordinate) {
                for (std::size_t offset = 0; offset < coefficients.size();) {
                    const std::size_t count = std::min(block_size, coefficients.size() - offset);
                    for (std::size_t i = 0; i < count; ++i) {
                        std::array<std::uint64_t, chunk_count> words {};
                        // Use the canonical value, not the field's internal Montgomery representation.
                        const auto canonical = coefficients[offset + i].coordinate(coordinate).to_integral();
                        boost::multiprecision::export_bits(canonical, words.begin(), chunk_bits, false);
                        for (std::size_t j = 0; j < chunk_count; ++j) {
                            chunks[i][j] = static_cast<double>(words[j]);
                        }
                    }
                    for (std::size_t prime = 0; prime < primes_.size(); ++prime) {
                        const auto &context = primes_[prime];
                        const auto &arithmetic = context.arithmetic;
                        const double modulus = static_cast<double>(arithmetic.modulus());
                        double *destination =
                            output.data() + (prime * field_type::arity + coordinate) * stride + offset;
                        for (std::size_t i = 0; i < count; ++i) {
                            double sum;
                            if constexpr (UseFma) {
                                sum = arithmetic.reduce_fma(chunks[i][0]);
                            } else {
                                sum = arithmetic.reduce_integer(chunks[i][0]);
                            }
                            // Each chunk is < 2^52 and each weight has magnitude < q/2, so the
                            // products are < 4*q^2 (q > 2^49). The initial remainder is < q;
                            // adding four terms of magnitude < 1.5*q keeps every partial sum
                            // inside (-7*q,7*q), hence below 2^53 and exactly representable.
                            for (std::size_t j = 1; j < chunk_count; ++j) {
                                if constexpr (UseFma) {
                                    sum += arithmetic.multiply_fma(chunks[i][j], context.powers[j - 1]);
                                } else {
                                    sum += arithmetic.multiply_integer(chunks[i][j], context.powers[j - 1]);
                                }
                            }
                            if constexpr (UseFma) {
                                sum = arithmetic.reduce_fma(sum);
                            } else {
                                sum = arithmetic.reduce_integer(sum);
                            }
                            destination[i] = sum < 0 ? sum + modulus : sum;
                        }
                    }
                    offset += count;
                }
            }
        }

        std::vector<prime_context> primes_;
    };

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_DETAIL_BN254_SMALL_PRIME_FFT_HPP
