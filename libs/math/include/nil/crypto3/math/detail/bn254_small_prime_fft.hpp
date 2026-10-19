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
#include <bit>
#include <cassert>
#include <cmath>
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
#include <nil/crypto3/multiprecision/inverse.hpp>

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

    /**
     * Reconstruct signed integer coordinates from a prefix of the small-prime catalogue, then
     * reduce them into BN254 Fp12. The supplied bound C must cover the original integer expression,
     * with M > 2*C for the prime product M. Centering must precede reduction into the target field:
     * the CRT representative M-1 denotes -1, not the field element (M-1) mod p.
     *
     * Cofactors and their inverses are prepared once. Reconstruction uses fixed-size integer scratch
     * and caller-owned output, with no arbitrary-precision setup work in the coefficient loop.
     * This operation has no constant-time guarantee.
     */
    class bn254_small_prime_fft_reconstructor {
        using base_field_type = algebra::fields::alt_bn128_base_field<254>;
        static constexpr std::size_t block_size = 32;
        static constexpr std::size_t correction_count = std::bit_width(small_prime_fft_primes.size() - 1);
        // Each CRT term is < M, so their sum is < prime_count*M < 12*2^600 < 2^604.
        static constexpr unsigned crt_bits = 50 * small_prime_fft_primes.size() + correction_count;
        using integer_type = boost::multiprecision::number<boost::multiprecision::cpp_int_modular_backend<crt_bits>>;
        static constexpr unsigned export_split = 384;
        using export_low_type =
            boost::multiprecision::number<boost::multiprecision::cpp_int_modular_backend<export_split>>;
        using export_type = boost::multiprecision::uint512_modular_t;
        static_assert(50 * small_prime_fft_primes.size() <= 600, "CRT export folding requires M < 2^600");

        struct prime_context {
            small_prime_fft_arithmetic arithmetic;
            integer_type cofactor;
            double inverse;

            prime_context(std::uint64_t modulus, const integer_type &cofactor) :
                arithmetic(modulus), cofactor(cofactor) {
                using wide_type = boost::multiprecision::uint128_modular_t;
                const auto residue = (cofactor % modulus).convert_to<std::uint64_t>();
                inverse = static_cast<double>(boost::multiprecision::inverse_mod(wide_type(residue), wide_type(modulus))
                                                  .convert_to<std::uint64_t>());
            }
        };

    public:
        using field_type = bn254_small_prime_fft_converter::field_type;
        using value_type = field_type::value_type;

        bn254_small_prime_fft_reconstructor(const small_prime_fft_crt_params &params,
                                            const boost::multiprecision::cpp_int &coefficient_bound) {
            if (params.prime_count == 0 || params.prime_count > small_prime_fft_primes.size() ||
                coefficient_bound < 0) {
                throw std::invalid_argument("bn254_small_prime_fft_reconstructor: invalid CRT parameters");
            }
            boost::multiprecision::cpp_int product = 1;
            for (std::size_t i = 0; i < params.prime_count; ++i) {
                product *= small_prime_fft_primes[i].modulus;
            }
            if (params.modulus != product || 2 * coefficient_bound >= product) {
                throw std::invalid_argument(
                    "bn254_small_prime_fft_reconstructor: inconsistent modulus or insufficient capacity");
            }

            // The checks above establish positivity and capacity before conversion to fixed precision.
            using setup_backend = integer_type::backend_type::cpp_int_type;
            multiples_[0].backend().from_cpp_int(setup_backend(product.backend()));
            bound_.backend().from_cpp_int(setup_backend(coefficient_bound.backend()));
            half_modulus_ = multiples_[0] >> 1;
            const boost::multiprecision::cpp_int sum_bound = params.prime_count * product - 1;
            constexpr unsigned limb_bits = integer_type::backend_type::limb_bits;
            accumulator_limbs_ = (boost::multiprecision::msb(sum_bound) + 1 + limb_bits - 1) / limb_bits;
            export_radix_ = (export_type(1) << export_split) % export_type(base_field_type::modulus);
            for (std::size_t bit = 1; bit < correction_count; ++bit) {
                multiples_[bit] = multiples_[bit - 1] << 1;
            }
            primes_.reserve(params.prime_count);
            for (std::size_t i = 0; i < params.prime_count; ++i) {
                const auto modulus = small_prime_fft_primes[i].modulus;
                primes_.emplace_back(modulus, multiples_[0] / modulus);
            }
        }

        std::size_t prime_count() const {
            return primes_.size();
        }

        /**
         * Read residues[(prime*12 + coordinate)*stride + coefficient] for each requested output.
         * Require output.size() <= stride and residues.size() == prime_count()*12*stride.
         * Check dimensions without overflow before writing any output; throw invalid_argument on
         * mismatch. An empty output does no work; entries after output.size() in each channel are
         * not read. Input residues must already include any inverse-transform normalization.
         *
         * @pre Each consumed residue is an integer-valued double of magnitude < 2^53. The thread's
         * rounding mode must remain unchanged during the call; kernel selection occurs once.
         *
         * Throw overflow_error if a centered coordinate exceeds the declared bound; earlier output
         * coordinates may already have been written. This check cannot detect an incorrect original
         * bound when the integer expression has wrapped by M, so the caller must justify that bound.
         */
        void reconstruct(std::span<value_type> output, std::span<const double> residues, std::size_t stride) const {
            const std::size_t channels = field_type::arity * primes_.size();
            if (channels == 0 || output.size() > stride || stride > residues.size() / channels ||
                residues.size() != channels * stride) {
                throw std::invalid_argument("bn254_small_prime_fft_reconstructor: invalid buffer dimensions");
            }
            if (output.empty()) {
                return;
            }
            const bool use_fma = std::all_of(primes_.begin(), primes_.end(), [](const prime_context &prime) {
                return prime.arithmetic.fma_supported();
            });
            if (use_fma) {
                reconstruct_impl<true>(output, residues, stride);
            } else {
                reconstruct_impl<false>(output, residues, stride);
            }
        }

    private:
        void add_cofactor_product(integer_type &sum, const integer_type &cofactor, std::uint64_t digit) const {
            using limb_type = boost::multiprecision::limb_type;
            using wide_type = boost::multiprecision::double_limb_type;
            constexpr unsigned limb_bits = integer_type::backend_type::limb_bits;
            if constexpr (limb_bits >= 50) {
                auto *destination = sum.backend().limbs();
                const auto *source = cofactor.backend().limbs();
                wide_type carry = 0;
                // A digit is < 2^50. Each product plus the accumulator limb and carry fits
                // in a double limb. The complete sum is < prime_count*M, so its final carry
                // fits in the active width prepared by the constructor; higher limbs stay zero.
                for (std::size_t limb = 0; limb < accumulator_limbs_; ++limb) {
                    const wide_type value = wide_type(source[limb]) * digit + destination[limb] + carry;
                    destination[limb] = static_cast<limb_type>(value);
                    carry = value >> limb_bits;
                }
                assert(carry == 0);
                assert((destination[sum.backend().size() - 1] & ~integer_type::backend_type::upper_limb_mask) == 0);
            } else {
                // A 50-bit digit spans more than one limb on a 32-bit limb backend.
                sum += cofactor * digit;
            }
        }

        base_field_type::value_type export_magnitude(const integer_type &magnitude) const {
            // Centering and the coefficient-bound check have already established magnitude < M < 2^600.
            // Write magnitude = low + high*2^384 and replace 2^384 by its cached residue modulo p.
            // The fold is < 2^384 + 2^216*2^254 < 2^471, within the existing BN254 Barrett reducer's
            // range. This avoids its general-division fallback while retaining the field's normal
            // Montgomery conversion. The fixed-width low conversion intentionally discards high bits.
            const export_type high(magnitude >> export_split);
            export_type folded {export_low_type(magnitude)};
            if (high != 0) {
                folded += high * export_radix_;
            }
            return base_field_type::value_type(folded);
        }

        template<bool UseFma>
        void reconstruct_impl(std::span<value_type> output, std::span<const double> residues,
                              std::size_t stride) const {
            for (std::size_t coordinate = 0; coordinate < field_type::arity; ++coordinate) {
                for (std::size_t offset = 0; offset < output.size();) {
                    const std::size_t count = std::min(block_size, output.size() - offset);
                    std::array<integer_type, block_size> sums {};
                    for (std::size_t prime = 0; prime < primes_.size(); ++prime) {
                        const auto &context = primes_[prime];
                        const auto &arithmetic = context.arithmetic;
                        const double modulus = static_cast<double>(arithmetic.modulus());
                        const double *source =
                            residues.data() + (prime * field_type::arity + coordinate) * stride + offset;
                        std::size_t i = 0;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                        if constexpr (UseFma) {
                            if (count >= 4) {
                                const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic);
                                const __m256d inverse = _mm256_set1_pd(context.inverse);
                                const __m256d q = _mm256_set1_pd(modulus);
                                for (; count - i >= 4; i += 4) {
                                    __m256d digits =
                                        packed.multiply_fma(packed.reduce_fma(_mm256_loadu_pd(source + i)), inverse);
                                    const __m256d negative = _mm256_cmp_pd(digits, _mm256_setzero_pd(), _CMP_LT_OQ);
                                    digits = _mm256_add_pd(digits, _mm256_and_pd(negative, q));
                                    std::array<double, 4> canonical;
                                    _mm256_storeu_pd(canonical.data(), digits);
                                    for (std::size_t lane = 0; lane < canonical.size(); ++lane) {
                                        add_cofactor_product(sums[i + lane], context.cofactor,
                                                             static_cast<std::uint64_t>(canonical[lane]));
                                    }
                                }
                            }
                        }
#endif
                        for (; i < count; ++i) {
                            double digit;
                            if constexpr (UseFma) {
                                digit = arithmetic.multiply_fma(arithmetic.reduce_fma(source[i]), context.inverse);
                            } else {
                                digit =
                                    arithmetic.multiply_integer(arithmetic.reduce_integer(source[i]), context.inverse);
                            }
                            // Both operands have magnitude < q, so the product remainder is in (-q,q).
                            // Canonical digits give 0 <= digit*(M/q) < M and bound the complete sum.
                            const auto canonical = static_cast<std::uint64_t>(digit < 0 ? digit + modulus : digit);
                            add_cofactor_product(sums[i], context.cofactor, canonical);
                        }
                    }
                    for (std::size_t i = 0; i < count; ++i) {
                        auto &magnitude = sums[i];
                        // The sum is < prime_count*M. Subtract binary multiples of M from high to low
                        // to obtain [0,M), using at most four comparisons/subtractions and no division.
                        for (std::size_t bit = std::bit_width(primes_.size() - 1); bit-- > 0;) {
                            if (magnitude >= multiples_[bit]) {
                                magnitude -= multiples_[bit];
                            }
                        }
                        const bool negative = magnitude > half_modulus_;
                        if (negative) {
                            magnitude = multiples_[0] - magnitude;
                        }
                        if (magnitude > bound_) {
                            throw std::overflow_error(
                                "bn254_small_prime_fft_reconstructor: coefficient exceeds CRT bound");
                        }
                        const auto value = export_magnitude(magnitude);
                        output[offset + i].coordinate(coordinate) = negative ? -value : value;
                    }
                    offset += count;
                }
            }
        }

        std::vector<prime_context> primes_;
        std::array<integer_type, correction_count> multiples_;
        integer_type half_modulus_;
        integer_type bound_;
        std::size_t accumulator_limbs_;
        export_type export_radix_;
    };

    /**
     * In-place, power-of-two transforms for one prime from the catalogue. The forward transform
     * maps natural-order coefficients to evaluations at omega^k in bit-reversed index order.
     * The inverse consumes that same order, so convolution needs no permutation passes.
     *
     * A reusable table of max(1,N/2) centered twiddles serves both directions. Radix-4 recursion
     * fuses pairs of layers, with radix-2 leaves for odd depths and special cases for unit twiddles.
     * AVX2/FMA builds process four independent butterflies at a time; short layers stay scalar.
     * Transform calls allocate no buffers; recursion uses O(log N) stack space.
     *
     * @pre Inputs are integer-valued doubles in (-3*q,3*q). The calling thread's rounding mode
     * must remain unchanged during a call. Each call selects FMA or integer arithmetic once.
     * The forward output is in (-3*q,3*q); the normalized inverse output is in [0,q).
     * This operation has no constant-time guarantee.
     */
    class small_prime_fft_plan {
    public:
        small_prime_fft_plan(std::size_t prime_index, std::size_t size) :
            arithmetic_(checked_prime(prime_index, size).modulus), size_(size),
            twiddles_(std::max(std::size_t(1), size / 2), 1.0) {
            using wide_type = boost::multiprecision::uint128_modular_t;
            const auto &prime = small_prime_fft_primes[prime_index];
            wide_type root = prime.root;
            const unsigned depth = std::countr_zero(size_);
            for (unsigned level = prime.two_adicity; level > depth; --level) {
                root = root * root % prime.modulus;
            }
            omega_ = root.convert_to<std::uint64_t>();
            size_inverse_ = centered(
                static_cast<double>(boost::multiprecision::inverse_mod(wide_type(size_), wide_type(prime.modulus))
                                        .convert_to<std::uint64_t>()));

            // At a layer with 2^d blocks, block j uses omega_(2^(d+1))^bitreverse(j,d).
            // Each new table segment is the preceding prefix times a primitive root of twice
            // the previous order: [1], [i], [w8,i*w8], [w16,i*w16,w8*w16,i*w8*w16], ... .
            for (std::size_t half = size_ / 4; half != 0; half >>= 1) {
                twiddles_[half] = centered(static_cast<double>(root.convert_to<std::uint64_t>()));
                root = root * root % prime.modulus;
            }
            for (std::size_t half = 1; half < size_ / 2; half <<= 1) {
                const double step = twiddles_[half];
                for (std::size_t i = 1; i < half; ++i) {
                    twiddles_[half + i] = centered(arithmetic_.multiply_integer(twiddles_[i], step));
                }
            }
        }

        std::size_t size() const {
            return size_;
        }

        std::uint64_t modulus() const {
            return arithmetic_.modulus();
        }

        std::uint64_t omega() const {
            return omega_;
        }

        /** Exact centered representative of 1/N modulo q, also usable to fuse scaling into pointwise work. */
        double size_inverse() const {
            return size_inverse_;
        }

        /** Require exactly N entries; invalid dimensions throw before any writes. N=1 is the identity. */
        void fft(std::span<double> values) const {
            check_buffer(values);
            if (arithmetic_.fma_supported()) {
                forward_recursive<true>(values.data(), size_, 0);
            } else {
                forward_recursive<false>(values.data(), size_, 0);
            }
        }

        /**
         * Require exactly N entries in the forward spectrum order. By default, scale by 1/N and
         * return canonical residues. With normalize=false, return N times the inverse in (-2*q,2*q),
         * allowing the caller to fuse the scaling elsewhere. Invalid dimensions throw before writes.
         */
        void inverse_fft(std::span<double> values, bool normalize = true) const {
            check_buffer(values);
            if (arithmetic_.fma_supported()) {
                inverse_impl<true>(values, normalize);
            } else {
                inverse_impl<false>(values, normalize);
            }
        }

    private:
        static const small_prime_fft_prime &checked_prime(std::size_t prime_index, std::size_t size) {
            if (prime_index >= small_prime_fft_primes.size() || !std::has_single_bit(size) ||
                std::countr_zero(size) > small_prime_fft_primes[prime_index].two_adicity) {
                throw std::invalid_argument("small_prime_fft_plan: invalid prime index or transform size");
            }
            return small_prime_fft_primes[prime_index];
        }

        void check_buffer(std::span<const double> values) const {
            if (values.size() != size_ || twiddles_.empty()) {
                throw std::invalid_argument("small_prime_fft_plan: invalid buffer dimensions or moved-from plan");
            }
#ifndef NDEBUG
            const double limit = 3 * static_cast<double>(modulus());
            for (double value : values) {
                assert(std::isfinite(value) && std::trunc(value) == value && std::abs(value) < limit);
            }
#endif
        }

        double centered(double value) const {
            const double q = static_cast<double>(modulus());
            return value > q / 2 ? value - q : (value < -q / 2 ? value + q : value);
        }

        double inverse_twiddle(std::size_t index) const {
            if (index == 0) {
                return 1.0;
            }
            // Within [b,2*b), reflected exponents add to half the root order, hence w^-1=-w_reflected.
            const std::size_t base = std::bit_floor(index);
            return -twiddles_[base + (2 * base - 1 - index)];
        }

        // Keep arithmetic selection outside the butterfly loops, with compile-time dispatch inside them.
        template<bool UseFma, typename Arithmetic, typename Value>
        static Value reduce(const Arithmetic &arithmetic, Value value) {
            if constexpr (UseFma) {
                return arithmetic.reduce_fma(value);
            } else {
                return arithmetic.reduce_integer(value);
            }
        }

        template<bool UseFma, typename Arithmetic, typename Value>
        static Value multiply(const Arithmetic &arithmetic, Value left, Value right) {
            if constexpr (UseFma) {
                return arithmetic.multiply_fma(left, right);
            } else {
                return arithmetic.multiply_integer(left, right);
            }
        }

        template<bool UseFma, bool UnitTwiddle, typename Arithmetic, typename Value>
        static void forward_butterfly(const Arithmetic &arithmetic, Value &v0, Value &v1, Value &v2, Value &v3,
                                      Value w2, Value w, Value iw) {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
            const Value a0 = reduce<UseFma>(arithmetic, v0);
            const Value a1 = v1;
            const Value a2 = UnitTwiddle ? reduce<UseFma>(arithmetic, v2) : multiply<UseFma>(arithmetic, v2, w2);
            const Value a3 = UnitTwiddle ? reduce<UseFma>(arithmetic, v3) : multiply<UseFma>(arithmetic, v3, w2);
            const Value even_sum = a0 + a2;
            const Value even_difference = a0 - a2;
            const Value odd_sum =
                UnitTwiddle ? reduce<UseFma>(arithmetic, a1 + a3) : multiply<UseFma>(arithmetic, a1 + a3, w);
            const Value odd_difference = multiply<UseFma>(arithmetic, a1 - a3, iw);
            // Input magnitudes < 3*q and weights < q/2 imply |a2|,|a3| < q and
            // |a1 +/- a3| < 4*q. Every product is < 2*q^2, so the odd terms are < q.
            // Thus all outputs are < 3*q, with every integer sum exact below 2^53.
            v0 = even_sum + odd_sum;
            v1 = even_sum - odd_sum;
            v2 = even_difference + odd_difference;
            v3 = even_difference - odd_difference;
        }

        template<bool UseFma, bool UnitTwiddle>
        void forward_radix4(double *values, std::size_t quarter, std::size_t node) const {
            const double w2 = twiddles_[node];
            const double w = twiddles_[2 * node];
            const double iw = twiddles_[2 * node + 1];
            std::size_t i = 0;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
                if (quarter >= 4) {
                    const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
                    const __m256d packed_w2 = _mm256_set1_pd(w2);
                    const __m256d packed_w = _mm256_set1_pd(w);
                    const __m256d packed_iw = _mm256_set1_pd(iw);
                    for (; quarter - i >= 4; i += 4) {
                        __m256d v0 = _mm256_loadu_pd(values + i);
                        __m256d v1 = _mm256_loadu_pd(values + i + quarter);
                        __m256d v2 = _mm256_loadu_pd(values + i + 2 * quarter);
                        __m256d v3 = _mm256_loadu_pd(values + i + 3 * quarter);
                        forward_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, packed_w2, packed_w, packed_iw);
                        _mm256_storeu_pd(values + i, v0);
                        _mm256_storeu_pd(values + i + quarter, v1);
                        _mm256_storeu_pd(values + i + 2 * quarter, v2);
                        _mm256_storeu_pd(values + i + 3 * quarter, v3);
                    }
                }
            }
#endif
            for (; i < quarter; ++i) {
                forward_butterfly<UseFma, UnitTwiddle>(arithmetic_, values[i], values[i + quarter],
                                                       values[i + 2 * quarter], values[i + 3 * quarter], w2, w, iw);
            }
        }

        template<bool UseFma>
        void forward_recursive(double *values, std::size_t length, std::size_t node) const {
            if (length == 1) {
                return;
            }
            if (length == 2) {
                const double left = reduce<UseFma>(arithmetic_, values[0]);
                const double right = node == 0 ? reduce<UseFma>(arithmetic_, values[1]) :
                                                 multiply<UseFma>(arithmetic_, values[1], twiddles_[node]);
                values[0] = left + right;
                values[1] = left - right;
                return;
            }
            const std::size_t quarter = length / 4;
            if (node == 0) {
                forward_radix4<UseFma, true>(values, quarter, node);
            } else {
                forward_radix4<UseFma, false>(values, quarter, node);
            }
            if (quarter > 1) {
                for (std::size_t child = 0; child < 4; ++child) {
                    forward_recursive<UseFma>(values + child * quarter, quarter, 4 * node + child);
                }
            }
        }

        template<bool UseFma, bool UnitTwiddle, typename Arithmetic, typename Value>
        static void inverse_butterfly(const Arithmetic &arithmetic, Value &v0, Value &v1, Value &v2, Value &v3,
                                      Value w2, Value w, Value iw, bool leaf) {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
            // Leaf inputs may be forward residues of magnitude < 3*q. Reduce them before
            // adding four entries; non-leaf inputs already satisfy the inverse bound < 2*q.
            const Value b0 = leaf ? reduce<UseFma>(arithmetic, v0) : v0;
            const Value b1 = leaf ? reduce<UseFma>(arithmetic, v1) : v1;
            const Value b2 = leaf ? reduce<UseFma>(arithmetic, v2) : v2;
            const Value b3 = leaf ? reduce<UseFma>(arithmetic, v3) : v3;
            const Value even_sum = b0 + b1;
            const Value odd_sum = b2 + b3;
            const Value left =
                UnitTwiddle ? reduce<UseFma>(arithmetic, b0 - b1) : multiply<UseFma>(arithmetic, b0 - b1, w);
            const Value right = multiply<UseFma>(arithmetic, b2 - b3, iw);
            // Four inputs < 2*q give sums < 8*q < 2^53. Times a centered weight this is
            // < 4*q^2, so the a2 remainder is < 1.5*q. Left and right are each < q;
            // their sum is < 2*q and their weighted difference is < q. The invariant closes.
            v0 = reduce<UseFma>(arithmetic, even_sum + odd_sum);
            v1 = left + right;
            v2 = UnitTwiddle ? reduce<UseFma>(arithmetic, even_sum - odd_sum) :
                               multiply<UseFma>(arithmetic, even_sum - odd_sum, w2);
            v3 =
                UnitTwiddle ? reduce<UseFma>(arithmetic, left - right) : multiply<UseFma>(arithmetic, left - right, w2);
        }

        template<bool UseFma, bool UnitTwiddle>
        void inverse_radix4(double *values, std::size_t quarter, std::size_t node) const {
            const double w2 = inverse_twiddle(node);
            const double w = inverse_twiddle(2 * node);
            const double iw = inverse_twiddle(2 * node + 1);
            std::size_t i = 0;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
                if (quarter >= 4) {
                    const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
                    const __m256d packed_w2 = _mm256_set1_pd(w2);
                    const __m256d packed_w = _mm256_set1_pd(w);
                    const __m256d packed_iw = _mm256_set1_pd(iw);
                    for (; quarter - i >= 4; i += 4) {
                        __m256d v0 = _mm256_loadu_pd(values + i);
                        __m256d v1 = _mm256_loadu_pd(values + i + quarter);
                        __m256d v2 = _mm256_loadu_pd(values + i + 2 * quarter);
                        __m256d v3 = _mm256_loadu_pd(values + i + 3 * quarter);
                        inverse_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, packed_w2, packed_w, packed_iw,
                                                             false);
                        _mm256_storeu_pd(values + i, v0);
                        _mm256_storeu_pd(values + i + quarter, v1);
                        _mm256_storeu_pd(values + i + 2 * quarter, v2);
                        _mm256_storeu_pd(values + i + 3 * quarter, v3);
                    }
                }
            }
#endif
            for (; i < quarter; ++i) {
                inverse_butterfly<UseFma, UnitTwiddle>(arithmetic_, values[i], values[i + quarter],
                                                       values[i + 2 * quarter], values[i + 3 * quarter], w2, w, iw,
                                                       quarter == 1);
            }
        }

        template<bool UseFma>
        void inverse_recursive(double *values, std::size_t length, std::size_t node) const {
            if (length == 1) {
                values[0] = reduce<UseFma>(arithmetic_, values[0]);
                return;
            }
            if (length == 2) {
                const double left = reduce<UseFma>(arithmetic_, values[0]);
                const double right = reduce<UseFma>(arithmetic_, values[1]);
                values[0] = reduce<UseFma>(arithmetic_, left + right);
                values[1] = node == 0 ? reduce<UseFma>(arithmetic_, left - right) :
                                        multiply<UseFma>(arithmetic_, left - right, inverse_twiddle(node));
                return;
            }
            const std::size_t quarter = length / 4;
            if (quarter > 1) {
                for (std::size_t child = 0; child < 4; ++child) {
                    inverse_recursive<UseFma>(values + child * quarter, quarter, 4 * node + child);
                }
            }
            if (node == 0) {
                inverse_radix4<UseFma, true>(values, quarter, node);
            } else {
                inverse_radix4<UseFma, false>(values, quarter, node);
            }
        }

        template<bool UseFma>
        void inverse_impl(std::span<double> values, bool normalize) const {
            inverse_recursive<UseFma>(values.data(), size_, 0);
            if (normalize) {
                const double q = static_cast<double>(modulus());
                std::size_t i = 0;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                if constexpr (UseFma) {
                    if (values.size() >= 4) {
                        const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
                        const __m256d scale = _mm256_set1_pd(size_inverse_);
                        const __m256d packed_q = _mm256_set1_pd(q);
                        for (; values.size() - i >= 4; i += 4) {
                            __m256d value = packed.multiply_fma(_mm256_loadu_pd(values.data() + i), scale);
                            const __m256d negative = _mm256_cmp_pd(value, _mm256_setzero_pd(), _CMP_LT_OQ);
                            value = _mm256_add_pd(value, _mm256_and_pd(negative, packed_q));
                            _mm256_storeu_pd(values.data() + i, value);
                        }
                    }
                }
#endif
                for (; i < values.size(); ++i) {
                    double &value = values[i];
                    value = multiply<UseFma>(arithmetic_, value, size_inverse_);
                    value = value < 0 ? value + q : value;
                }
            }
        }

        small_prime_fft_arithmetic arithmetic_;
        std::size_t size_;
        std::uint64_t omega_;
        double size_inverse_;
        std::vector<double> twiddles_;
    };

    /**
     * Pointwise products in the BN254 tower u^2=-1, v^3=9+u, w^2=v, over one small modulus q.
     * This is a quotient ring; irreducibility of the tower polynomials modulo q is not required.
     * Coordinate channels have the same order as the converter, with one prime per call.
     *
     * Reuse the field library's Karatsuba tower identities and quadratic-extension squaring
     * identity with scalar or four-lane residue kernels. The existing field elements bind their base
     * arithmetic to the target modulus and cannot represent these small-prime intermediates.
     *
     * @pre Input cells are integer-valued doubles of magnitude < 2^53. scale is an integer with
     * |scale| < q/2, normally the FFT plan's centered 1/N. Buffers are disjoint or identical.
     * The thread's rounding mode stays unchanged during the call. Output cells are in (-q,q).
     * No buffers are allocated and no constant-time guarantee is provided.
     */
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
// The packed kernel uses __m256d as the coordinate type of the existing tower arrays.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
    class bn254_small_prime_fft_pointwise {

        template<bool UseFma, typename Value, typename Arithmetic>
        struct kernel {
            using quadratic_type = std::array<Value, 2>;
            using cubic_type = std::array<quadratic_type, 3>;
            using tower_type = std::array<cubic_type, 2>;

            const Arithmetic &arithmetic;

            Value reduce(Value value) const {
                if constexpr (UseFma) {
                    return arithmetic.reduce_fma(value);
                } else {
                    return arithmetic.reduce_integer(value);
                }
            }

            Value multiply(Value a, Value b) const {
                if constexpr (UseFma) {
                    return arithmetic.multiply_fma(a, b);
                } else {
                    return arithmetic.multiply_integer(a, b);
                }
            }

            Value add(Value a, Value b) const {
                return reduce(a + b);
            }

            Value subtract(Value a, Value b) const {
                return reduce(a - b);
            }

            template<typename T, std::size_t N>
            std::array<T, N> add(const std::array<T, N> &a, const std::array<T, N> &b) const {
                std::array<T, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                    result[i] = add(a[i], b[i]);
                }
                return result;
            }

            template<typename T, std::size_t N>
            std::array<T, N> subtract(const std::array<T, N> &a, const std::array<T, N> &b) const {
                std::array<T, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                    result[i] = subtract(a[i], b[i]);
                }
                return result;
            }

            // All tower inputs and returned coordinates are in (-q,q). Pair sums are < 2*q;
            // the Karatsuba product is < 4*q^2 and its remainder < 1.5*q. Correcting by the
            // two diagonal products gives a sum < 3.5*q, hence exact in binary64 before reduction.
            quadratic_type multiply(const quadratic_type &a, const quadratic_type &b) const {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
                const Value diagonal0 = multiply(a[0], b[0]);
                const Value diagonal1 = multiply(a[1], b[1]);
                const Value cross = multiply(a[0] + a[1], b[0] + b[1]);
                return {reduce(diagonal0 - diagonal1), reduce(cross - diagonal0 - diagonal1)};
            }

            quadratic_type multiply_xi(const quadratic_type &a) const {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
                // (9+u)*(a0+a1*u) = (9*a0-a1) + (a0+9*a1)*u. Reduce 8*a first:
                // |8*a| < 8*q < 2^53, whereas forming 9*a directly could lose an integer bit.
                return {reduce(reduce(8 * a[0]) + a[0] - a[1]), reduce(reduce(8 * a[1]) + a[1] + a[0])};
            }

            cubic_type multiply_v(const cubic_type &a) const {
                return {multiply_xi(a[2]), a[0], a[1]};
            }

            cubic_type multiply(const cubic_type &a, const cubic_type &b) const {
                const auto diagonal0 = multiply(a[0], b[0]);
                const auto diagonal1 = multiply(a[1], b[1]);
                const auto diagonal2 = multiply(a[2], b[2]);
                const auto cross12 =
                    subtract(subtract(multiply(add(a[1], a[2]), add(b[1], b[2])), diagonal1), diagonal2);
                const auto cross01 =
                    subtract(subtract(multiply(add(a[0], a[1]), add(b[0], b[1])), diagonal0), diagonal1);
                const auto cross02 =
                    subtract(subtract(multiply(add(a[0], a[2]), add(b[0], b[2])), diagonal0), diagonal2);
                return {add(diagonal0, multiply_xi(cross12)), add(cross01, multiply_xi(diagonal2)),
                        add(cross02, diagonal1)};
            }

            tower_type multiply(const tower_type &a, const tower_type &b) const {
                const auto diagonal0 = multiply(a[0], b[0]);
                const auto diagonal1 = multiply(a[1], b[1]);
                const auto cross = subtract(subtract(multiply(add(a[0], a[1]), add(b[0], b[1])), diagonal0), diagonal1);
                return {add(diagonal0, multiply_v(diagonal1)), cross};
            }

            tower_type square(const tower_type &a) const {
                // The quadratic-extension square needs two cubic products instead of three:
                // c0 = (a0+a1)*(a0+v*a1) - a0*a1 - v*a0*a1; c1 = 2*a0*a1.
                const auto cross = multiply(a[0], a[1]);
                const auto combined = multiply(add(a[0], a[1]), add(a[0], multiply_v(a[1])));
                return {subtract(subtract(combined, cross), multiply_v(cross)), add(cross, cross)};
            }
        };

    public:
        explicit bn254_small_prime_fft_pointwise(std::uint64_t modulus) : arithmetic_(modulus) {
        }

        /** Overwrite left with scale*left*right. Both buffers must contain exactly 12*stride cells. */
        void multiply(std::span<double> left, std::span<const double> right, std::size_t stride, double scale) const {
            check_buffers(left, right, stride, scale);
            if (arithmetic_.fma_supported()) {
                apply<true, false>(left, left, right, scale);
            } else {
                apply<false, false>(left, left, right, scale);
            }
        }

        /** Overwrite values with scale*values^2, loading each tower element once. */
        void square(std::span<double> values, std::size_t stride, double scale) const {
            check_buffers(values, values, stride, scale);
            if (arithmetic_.fma_supported()) {
                apply<true, true>(values, values, values, scale);
            } else {
                apply<false, true>(values, values, values, scale);
            }
        }

        /** Compute scale*(a*b+c*d). Every buffer has 12*stride cells; output may equal any input. */
        void multiply_sum(std::span<double> output, std::span<const double> a, std::span<const double> b,
                          std::span<const double> c, std::span<const double> d, std::size_t stride,
                          double scale) const {
            check_buffers(output, a, stride, scale);
            check_buffers(output, b, stride, scale);
            check_buffers(output, c, stride, scale);
            check_buffers(output, d, stride, scale);
            if (arithmetic_.fma_supported()) {
                apply<true, false, true>(output, a, b, scale, c, d);
            } else {
                apply<false, false, true>(output, a, b, scale, c, d);
            }
        }

    private:
        void check_buffers(std::span<const double> left, std::span<const double> right, std::size_t stride,
                           double scale) const {
            if (left.size() % 12 != 0 || left.size() / 12 != stride || left.size() != right.size()) {
                throw std::invalid_argument("bn254_small_prime_fft_pointwise: invalid buffer dimensions");
            }
            assert(std::isfinite(scale) && std::trunc(scale) == scale &&
                   std::abs(scale) < static_cast<double>(arithmetic_.modulus()) / 2);
        }

        template<bool UseFma, bool Square, bool Sum = false>
        void apply(std::span<double> output, std::span<const double> left, std::span<const double> right, double scale,
                   std::span<const double> extra_left = {}, std::span<const double> extra_right = {}) const {
            // Derive the indexing bound from the buffer extent checked by the public entry points.
            const std::size_t stride = output.size() / 12;
            std::size_t offset = 0;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
                offset = stride - stride % 4;
                if (offset != 0) {
                    const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
                    apply_range<true, Square, Sum>(packed, output, left, right, extra_left, extra_right, 0, offset,
                                                   _mm256_set1_pd(scale));
                }
            }
#endif
            apply_range<UseFma, Square, Sum>(arithmetic_, output, left, right, extra_left, extra_right, offset, stride,
                                             scale);
        }

        template<bool UseFma, bool Square, bool Sum, typename Arithmetic, typename Value>
        static void apply_range(const Arithmetic &arithmetic, std::span<double> output, std::span<const double> left,
                                std::span<const double> right, std::span<const double> extra_left,
                                std::span<const double> extra_right, std::size_t begin, std::size_t end, Value scale) {
            constexpr std::size_t width = sizeof(Value) / sizeof(double);
            const std::size_t stride = output.size() / 12;
            const kernel<UseFma, Value, Arithmetic> operations {arithmetic};
            using tower_type = typename kernel<UseFma, Value, Arithmetic>::tower_type;
            for (std::size_t i = begin; i < end; i += width) {
                tower_type a, b, c, d;
                for (std::size_t w = 0; w < 2; ++w) {
                    for (std::size_t v = 0; v < 3; ++v) {
                        for (std::size_t u = 0; u < 2; ++u) {
                            const std::size_t offset = (6 * w + 2 * v + u) * stride + i;
                            if constexpr (width == 1) {
                                a[w][v][u] = operations.reduce(left[offset]);
                                if constexpr (!Square) {
                                    b[w][v][u] = operations.reduce(right[offset]);
                                }
                                if constexpr (Sum) {
                                    c[w][v][u] = operations.reduce(extra_left[offset]);
                                    d[w][v][u] = operations.reduce(extra_right[offset]);
                                }
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                            } else {
                                a[w][v][u] = operations.reduce(_mm256_loadu_pd(left.data() + offset));
                                if constexpr (!Square) {
                                    b[w][v][u] = operations.reduce(_mm256_loadu_pd(right.data() + offset));
                                }
                                if constexpr (Sum) {
                                    c[w][v][u] = operations.reduce(_mm256_loadu_pd(extra_left.data() + offset));
                                    d[w][v][u] = operations.reduce(_mm256_loadu_pd(extra_right.data() + offset));
                                }
#endif
                            }
                        }
                    }
                }
                const tower_type result = [&]() {
                    if constexpr (Sum) {
                        // Each product coordinate is in (-q,q); adding two is exact below 2*q.
                        // Reduce that sum before the usual inverse-transform normalization.
                        return operations.add(operations.multiply(a, b), operations.multiply(c, d));
                    } else if constexpr (Square) {
                        return operations.square(a);
                    } else {
                        return operations.multiply(a, b);
                    }
                }();
                for (std::size_t w = 0; w < 2; ++w) {
                    for (std::size_t v = 0; v < 3; ++v) {
                        for (std::size_t u = 0; u < 2; ++u) {
                            // Fold inverse-transform normalization into this existing write pass.
                            const std::size_t offset = (6 * w + 2 * v + u) * stride + i;
                            const Value value = operations.multiply(result[w][v][u], scale);
                            if constexpr (width == 1) {
                                output[offset] = value;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                            } else {
                                _mm256_storeu_pd(output.data() + offset, value);
#endif
                            }
                        }
                    }
                }
            }
        }

        small_prime_fft_arithmetic arithmetic_;
    };
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
#pragma GCC diagnostic pop
#endif

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_DETAIL_BN254_SMALL_PRIME_FFT_HPP
