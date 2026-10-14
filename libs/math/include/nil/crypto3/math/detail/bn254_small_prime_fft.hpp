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
        static constexpr std::size_t block_size = 32;
        static constexpr std::size_t correction_count = std::bit_width(small_prime_fft_primes.size() - 1);
        // Each CRT term is < M, so their sum is < prime_count*M < 12*2^600 < 2^604.
        static constexpr unsigned crt_bits = 50 * small_prime_fft_primes.size() + correction_count;
        using integer_type = boost::multiprecision::number<boost::multiprecision::cpp_int_modular_backend<crt_bits>>;

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
        template<bool UseFma>
        void reconstruct_impl(std::span<value_type> output, std::span<const double> residues,
                              std::size_t stride) const {
            using base_field_type = algebra::fields::alt_bn128_base_field<254>;
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
                        for (std::size_t i = 0; i < count; ++i) {
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
                            sums[i] += context.cofactor * canonical;
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
                        const base_field_type::value_type value(magnitude);
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
    };

    /**
     * In-place, power-of-two transforms for one prime from the catalogue. The forward transform
     * maps natural-order coefficients to evaluations at omega^k in bit-reversed index order.
     * The inverse consumes that same order, so convolution needs no permutation passes.
     *
     * A reusable table of max(1,N/2) centered twiddles serves both directions. Radix-4 recursion
     * fuses pairs of layers, with radix-2 leaves for odd depths and special cases for unit twiddles.
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
        template<bool UseFma>
        double reduce(double value) const {
            if constexpr (UseFma) {
                return arithmetic_.reduce_fma(value);
            } else {
                return arithmetic_.reduce_integer(value);
            }
        }

        template<bool UseFma>
        double multiply(double left, double right) const {
            if constexpr (UseFma) {
                return arithmetic_.multiply_fma(left, right);
            } else {
                return arithmetic_.multiply_integer(left, right);
            }
        }

        template<bool UseFma, bool UnitTwiddle>
        void forward_radix4(double *values, std::size_t quarter, std::size_t node) const {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
            const double w2 = twiddles_[node];
            const double w = twiddles_[2 * node];
            const double iw = twiddles_[2 * node + 1];
            for (std::size_t i = 0; i < quarter; ++i) {
                const double a0 = reduce<UseFma>(values[i]);
                const double a1 = values[i + quarter];
                const double a2 = UnitTwiddle ? reduce<UseFma>(values[i + 2 * quarter]) :
                                                multiply<UseFma>(values[i + 2 * quarter], w2);
                const double a3 = UnitTwiddle ? reduce<UseFma>(values[i + 3 * quarter]) :
                                                multiply<UseFma>(values[i + 3 * quarter], w2);
                const double even_sum = a0 + a2;
                const double even_difference = a0 - a2;
                const double odd_sum = UnitTwiddle ? reduce<UseFma>(a1 + a3) : multiply<UseFma>(a1 + a3, w);
                const double odd_difference = multiply<UseFma>(a1 - a3, iw);
                // Input magnitudes < 3*q and weights < q/2 imply |a2|,|a3| < q and
                // |a1 +/- a3| < 4*q. Every product is < 2*q^2, so the odd terms are < q.
                // Thus all outputs are < 3*q, with every integer sum exact below 2^53.
                values[i] = even_sum + odd_sum;
                values[i + quarter] = even_sum - odd_sum;
                values[i + 2 * quarter] = even_difference + odd_difference;
                values[i + 3 * quarter] = even_difference - odd_difference;
            }
        }

        template<bool UseFma>
        void forward_recursive(double *values, std::size_t length, std::size_t node) const {
            if (length == 1) {
                return;
            }
            if (length == 2) {
                const double left = reduce<UseFma>(values[0]);
                const double right =
                    node == 0 ? reduce<UseFma>(values[1]) : multiply<UseFma>(values[1], twiddles_[node]);
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

        template<bool UseFma, bool UnitTwiddle>
        void inverse_radix4(double *values, std::size_t quarter, std::size_t node) const {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
            const double w2 = inverse_twiddle(node);
            const double w = inverse_twiddle(2 * node);
            const double iw = inverse_twiddle(2 * node + 1);
            for (std::size_t i = 0; i < quarter; ++i) {
                // Leaf inputs may be forward residues of magnitude < 3*q. Reduce them before
                // adding four entries; non-leaf inputs already satisfy the inverse bound < 2*q.
                const double b0 = quarter == 1 ? reduce<UseFma>(values[i]) : values[i];
                const double b1 = quarter == 1 ? reduce<UseFma>(values[i + quarter]) : values[i + quarter];
                const double b2 = quarter == 1 ? reduce<UseFma>(values[i + 2 * quarter]) : values[i + 2 * quarter];
                const double b3 = quarter == 1 ? reduce<UseFma>(values[i + 3 * quarter]) : values[i + 3 * quarter];
                const double even_sum = b0 + b1;
                const double odd_sum = b2 + b3;
                const double left = UnitTwiddle ? reduce<UseFma>(b0 - b1) : multiply<UseFma>(b0 - b1, w);
                const double right = multiply<UseFma>(b2 - b3, iw);
                // Four inputs < 2*q give sums < 8*q < 2^53. Times a centered weight this is
                // < 4*q^2, so the a2 remainder is < 1.5*q. Left and right are each < q;
                // their sum is < 2*q and their weighted difference is < q. The invariant closes.
                values[i] = reduce<UseFma>(even_sum + odd_sum);
                values[i + quarter] = left + right;
                values[i + 2 * quarter] =
                    UnitTwiddle ? reduce<UseFma>(even_sum - odd_sum) : multiply<UseFma>(even_sum - odd_sum, w2);
                values[i + 3 * quarter] =
                    UnitTwiddle ? reduce<UseFma>(left - right) : multiply<UseFma>(left - right, w2);
            }
        }

        template<bool UseFma>
        void inverse_recursive(double *values, std::size_t length, std::size_t node) const {
            if (length == 1) {
                values[0] = reduce<UseFma>(values[0]);
                return;
            }
            if (length == 2) {
                const double left = reduce<UseFma>(values[0]);
                const double right = reduce<UseFma>(values[1]);
                values[0] = reduce<UseFma>(left + right);
                values[1] =
                    node == 0 ? reduce<UseFma>(left - right) : multiply<UseFma>(left - right, inverse_twiddle(node));
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
                for (double &value : values) {
                    value = multiply<UseFma>(value, size_inverse_);
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

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_DETAIL_BN254_SMALL_PRIME_FFT_HPP
