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

#ifndef CRYPTO3_MATH_DETAIL_BN254_SMALL_PRIME_FFT_CRT_HPP
#define CRYPTO3_MATH_DETAIL_BN254_SMALL_PRIME_FFT_CRT_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/limb_ops.hpp>
#include <nil/crypto3/math/polynomial/backends/detail/small_prime_fft/bn254_small_prime_fft_conversion.hpp>
#include <nil/crypto3/math/polynomial/backends/detail/small_prime_fft/small_prime_fft_arithmetic.hpp>
#include <nil/crypto3/math/polynomial/backends/detail/small_prime_fft/small_prime_fft_params.hpp>
#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>
#include <nil/crypto3/multiprecision/inverse.hpp>

#if defined(CRYPTO3_HAS_AVX512IFMA) && defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_AVAILABLE) && \
    !defined(CRYPTO3_SMALL_PRIME_CRT_DISABLE_IFMA)
#define CRYPTO3_SMALL_PRIME_CRT_IFMA_SELECTED
#endif

namespace nil::crypto3::math::detail {

    // Independent coefficient lanes, using the existing IFMA operation policy. A CRT digit is
    // below 2^50 and each five-digit export weight is below p. With at most twelve primes,
    // every unnormalized column is below 24*2^52 and the complete sum is below 2^54*p.
    template<typename Params, typename Operations>
    struct small_prime_fft_crt_accumulator {
        using radix = algebra::fields::detail::fp12_fast::montgomery_radix52<Params>;
        using vector = typename Operations::value_type;
        std::array<vector, 6> columns;

        small_prime_fft_crt_accumulator() {
            columns.fill(Operations::broadcast(0));
        }

        void add(const typename radix::value_type &weight, vector digits) {
            for (std::size_t limb = 0; limb < weight.size(); ++limb) {
                const auto w = Operations::broadcast(weight[limb]);
                columns[limb] = Operations::multiply_add_low(columns[limb], digits, w);
                columns[limb + 1] = Operations::multiply_add_high(columns[limb + 1], digits, w);
            }
        }

        // The certificate reciprocal has up to 61 bits, so it cannot be an IFMA operand
        // without splitting. Compute floor(d*r/2^50) exactly from four 52-bit split products.
        static vector certificate_term(vector digits, std::uint64_t reciprocal) {
            const auto zero = Operations::broadcast(0);
            const auto low = Operations::broadcast(reciprocal & radix::mask);
            const auto high = Operations::broadcast(reciprocal >> 52);
            const auto a = Operations::multiply_add_low(zero, digits, low);
            const auto b = Operations::multiply_add_high(zero, digits, low);
            const auto c = Operations::multiply_add_low(zero, digits, high);
            const auto d = Operations::multiply_add_high(zero, digits, high);
            return Operations::add(
                Operations::add(Operations::template shift_right<50>(a), Operations::template shift_left<2>(b)),
                Operations::add(Operations::template shift_left<2>(c), Operations::template shift_left<54>(d)));
        }

        struct certificate {
            vector quotients;
            vector accepted;
        };

        static certificate certify(vector lower, std::uint64_t scaled_bound, std::size_t prime_count) {
            constexpr std::uint64_t scale = std::uint64_t(1) << 60;
            const auto zero = Operations::broadcast(0);
            const auto one = Operations::broadcast(1);
            const auto quotient =
                Operations::template shift_right<60>(Operations::add(lower, Operations::broadcast(scale / 2)));
            const auto center = Operations::template shift_left<60>(quotient);
            const auto bound = Operations::broadcast(scaled_bound);
            const auto upper = Operations::add(lower, Operations::broadcast(2 * prime_count));
            const auto lower_ok =
                Operations::select(Operations::less(Operations::add(lower, bound), center), zero, one);
            const auto upper_ok =
                Operations::select(Operations::less(Operations::add(center, bound), upper), zero, one);
            return {quotient, Operations::bit_and(lower_ok, upper_ok)};
        }

        // Add the selected centering correction and divide by S=2^260. The caller prepares
        // weights and corrections with R*S (R=2^256), so the result is already R-Montgomery.
        std::array<vector, 5> reduce(const std::array<vector, 5> &correction) const {
            const auto zero = Operations::broadcast(0);
            const auto mask = Operations::broadcast(radix::mask);
            const auto inverse = Operations::broadcast(radix::p_dash());
            std::array<vector, 10> t;
            t.fill(zero);
            std::copy(columns.begin(), columns.end(), t.begin());
            for (std::size_t i = 0; i < correction.size(); ++i)
                t[i] = Operations::add(t[i], correction[i]);
            for (std::size_t i = 0; i < 9; ++i) {
                t[i + 1] = Operations::add(t[i + 1], Operations::shift52(t[i]));
                t[i] = Operations::bit_and(t[i], mask);
            }
            // Each cancellation adds at most five split products to a column. All columns
            // stay below 2^57; the initial sum plus correction is < 2^54*p < S*p.
            for (std::size_t i = 0; i < 5; ++i) {
                const auto m = Operations::multiply_add_low(zero, t[i], inverse);
                for (std::size_t j = 0; j < 5; ++j) {
                    const auto p = Operations::broadcast(radix::modulus()[j]);
                    t[i + j] = Operations::multiply_add_low(t[i + j], m, p);
                    t[i + j + 1] = Operations::multiply_add_high(t[i + j + 1], m, p);
                }
                t[i + 1] = Operations::add(t[i + 1], Operations::shift52(t[i]));
            }
            std::array<vector, 5> reduced;
            for (std::size_t i = 5; i < 9; ++i) {
                reduced[i - 5] = Operations::bit_and(t[i], mask);
                t[i + 1] = Operations::add(t[i + 1], Operations::shift52(t[i]));
            }
            reduced[4] = t[9];
            // REDC(T)<2*p<S; a single existing vector subtraction makes every lane canonical.
            algebra::fields::detail::fp12_fast::montgomery_radix52_simd<Params, Operations>::subtract_modulus(reduced);
            return reduced;
        }
    };

    /**
     * Reconstruct signed integer coordinates from a prefix of the small-prime catalogue, then
     * reduce them into BN254 Fp12. The supplied bound C must cover the original integer expression,
     * with M > 2*C for the prime product M. Centering must precede reduction into the target field:
     * the CRT representative M-1 denotes -1, not the field element (M-1) mod p.
     *
     * Cofactors, inverses and Montgomery export weights are prepared once. An integer interval
     * certificate permits direct export when it proves both centering and the coefficient bound;
     * uncertain cases use exact wide reconstruction. Both paths use fixed-size scratch and
     * caller-owned output, with no arbitrary-precision setup work in the coefficient loop.
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
        using fast_params = bn254_small_prime_fft_converter::field_type::extension_policy::fast_params;
        using direct_weight_type = std::array<std::uint64_t, fast_params::base_value_limb_count>;
#if defined(CRYPTO3_SMALL_PRIME_CRT_IFMA_SELECTED)
        using radix52 = algebra::fields::detail::fp12_fast::montgomery_radix52<fast_params>;
#endif
        static constexpr bool direct_export_supported = base_field_type::integral_type::backend_type::limb_bits == 64;
        static constexpr unsigned certificate_bits = 60;
        static constexpr std::uint64_t certificate_scale = std::uint64_t(1) << certificate_bits;
        static_assert(std::all_of(small_prime_fft_primes.begin(), small_prime_fft_primes.end(), [](const auto &prime) {
            return prime.modulus > (std::uint64_t(1) << 49) && prime.modulus < (std::uint64_t(1) << 50);
        }));

        struct prime_context {
            small_prime_fft_arithmetic arithmetic;
            integer_type cofactor;
            double inverse;
            std::uint64_t certificate_reciprocal;
            direct_weight_type export_weight;
#if defined(CRYPTO3_SMALL_PRIME_CRT_IFMA_SELECTED)
            radix52::value_type export_weight52;
#endif

            prime_context(std::uint64_t modulus, const integer_type &cofactor, const direct_weight_type &weight) :
                arithmetic(modulus), cofactor(cofactor), export_weight(weight) {
                using wide_type = boost::multiprecision::uint128_modular_t;
                const auto residue = (cofactor % modulus).convert_to<std::uint64_t>();
                inverse = static_cast<double>(boost::multiprecision::inverse_mod(wide_type(residue), wide_type(modulus))
                                                  .convert_to<std::uint64_t>());
                certificate_reciprocal =
                    ((wide_type(1) << (certificate_bits + 50)) / modulus).convert_to<std::uint64_t>();
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
            const boost::multiprecision::cpp_int characteristic(base_field_type::modulus.backend().to_cpp_int());
            const boost::multiprecision::cpp_int r_squared =
                (boost::multiprecision::cpp_int(1) << 512) % characteristic;
            scaled_bound_ = ((coefficient_bound << certificate_bits) / product).convert_to<std::uint64_t>();
            for (std::size_t i = 0; i < params.prime_count; ++i) {
                const auto modulus = small_prime_fft_primes[i].modulus;
                const boost::multiprecision::cpp_int weight = (product / modulus) * r_squared % characteristic;
                direct_weight_type limbs {};
                boost::multiprecision::export_bits(weight, limbs.begin(), 64, false);
                primes_.emplace_back(modulus, multiples_[0] / modulus, limbs);
#if defined(CRYPTO3_SMALL_PRIME_CRT_IFMA_SELECTED)
                // S/R=16. Export after radix-52 REDC must retain the field's R=2^256 scaling.
                const boost::multiprecision::cpp_int vector_weight = (weight * 16) % characteristic;
                limbs.fill(0);
                boost::multiprecision::export_bits(vector_weight, limbs.begin(), 64, false);
                primes_.back().export_weight52 = radix52::pack(limbs.data());
#endif
            }
            // The centering quotient is in [0,prime_count]. Store canonical -k*M*R^2 mod p.
            for (std::size_t k = 0; k <= params.prime_count; ++k) {
                boost::multiprecision::cpp_int correction = (k * product * r_squared) % characteristic;
                if (correction != 0) {
                    correction = characteristic - correction;
                }
                boost::multiprecision::export_bits(correction, export_corrections_[k].begin(), 64, false);
#if defined(CRYPTO3_SMALL_PRIME_CRT_IFMA_SELECTED)
                const boost::multiprecision::cpp_int vector_correction = (correction * 16) % characteristic;
                direct_weight_type limbs {};
                boost::multiprecision::export_bits(vector_correction, limbs.begin(), 64, false);
                const auto packed = radix52::pack(limbs.data());
                for (std::size_t digit = 0; digit < packed.size(); ++digit)
                    export_corrections52_[digit][k] = packed[digit];
#endif
            }
        }

        std::size_t prime_count() const {
            return primes_.size();
        }

        /**
         * Combine an inverse-transform scale with (M/q)^-1, returning a centered residue.
         * Require an integer-valued scale of magnitude < 2^53; invalid prime indices throw.
         */
        double weighted_normalization(std::size_t prime, double scale) const {
            const auto &context = primes_.at(prime);
            const double q = static_cast<double>(context.arithmetic.modulus());
            double result =
                context.arithmetic.multiply_integer(context.arithmetic.reduce_integer(scale), context.inverse);
            if (result > q / 2) {
                result -= q;
            } else if (result < -q / 2) {
                result += q;
            }
            return result;
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
            reconstruct_dispatch<false>(output, residues, stride);
        }

        /**
         * Reconstruct residues already multiplied by (M/q)^-1 modulo each q. All dimensions,
         * residue bounds and failure semantics are those of reconstruct(). The backend can fuse
         * this weight with 1/N in the pointwise pass because the inverse transform is linear.
         */
        void reconstruct_weighted(std::span<value_type> output, std::span<const double> residues,
                                  std::size_t stride) const {
            reconstruct_dispatch<true>(output, residues, stride);
        }

    private:
        template<bool Weighted>
        void reconstruct_dispatch(std::span<value_type> output, std::span<const double> residues,
                                  std::size_t stride) const {
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
#if defined(CRYPTO3_SMALL_PRIME_CRT_IFMA_SELECTED)
                reconstruct_ifma<Weighted>(output, residues, stride);
                return;
#endif
                reconstruct_impl<true, Weighted>(output, residues, stride);
            } else {
                reconstruct_impl<false, Weighted>(output, residues, stride);
            }
        }

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

        base_field_type::value_type export_exact(integer_type magnitude) const {
            // Each canonical CRT term is < M. Reduce their sum to [0,M) without division.
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
                throw std::overflow_error("bn254_small_prime_fft_reconstructor: coefficient exceeds CRT bound");
            }
            const auto value = export_magnitude(magnitude);
            return negative ? -value : value;
        }

        static void add_export_product(fast_params::limb_array &sum, const direct_weight_type &weight,
                                       std::uint64_t digit) {
            using wide_type = algebra::fields::detail::fp12_fast::wide_limb;
            wide_type carry = 0;
            for (std::size_t limb = 0; limb < weight.size(); ++limb) {
                const wide_type value = wide_type(weight[limb]) * digit + sum[limb] + carry;
                sum[limb] = static_cast<std::uint64_t>(value);
                carry = value >> 64;
            }
            // At most twelve digits < 2^50 times weights < p, plus one correction < p.
            // The total is < 2^54*p < 2^308, so the fifth limb cannot overflow.
            sum[weight.size()] += static_cast<std::uint64_t>(carry);
        }

        bool certify_export(std::uint64_t lower, std::size_t &quotient) const {
            // For d < q < 2^50 and a=floor(2^110/q), floor(d*a/2^50) <= 2^60*d/q
            // < floor(d*a/2^50)+2. Summing encloses 2^60*sum(d/q) in [lower,upper].
            // With at most twelve primes, all quantities below fit in uint64_t.
            const std::uint64_t upper = lower + 2 * primes_.size();
            quotient = (lower + certificate_scale / 2) >> certificate_bits;
            const std::uint64_t center = std::uint64_t(quotient) << certificate_bits;
            // Containment in k +/- C/M proves the bound AND the unique centering quotient,
            // since C/M < 1/2. Using floor(2^60*C/M) makes this test conservative.
            return lower + scaled_bound_ >= center && upper <= center + scaled_bound_;
        }

        template<typename Params = fast_params>
        base_field_type::value_type export_direct(typename Params::limb_array &sum, std::size_t quotient) const {
            add_export_product(sum, export_corrections_[quotient], 1);
            // sum = (sum(d_i*M_i)-k*M)*R^2 mod p, and sum < 2^54*p < R*p for R=2^256.
            // One existing REDC therefore yields a canonical Montgomery coordinate directly,
            // without materializing the wide CRT integer, folding, or a second field conversion.
            direct_weight_type reduced;
            algebra::fields::detail::fp12_fast::montgomery_reduce<Params>(reduced.data(), sum);
            base_field_type::value_type result;
            std::copy(reduced.begin(), reduced.end(), result.data.backend().base_data().limbs());
            return result;
        }

        template<bool UseFma, bool Weighted, bool Direct = direct_export_supported>
        void reconstruct_impl(std::span<value_type> output, std::span<const double> residues,
                              std::size_t stride) const {
            for (std::size_t coordinate = 0; coordinate < field_type::arity; ++coordinate) {
                for (std::size_t offset = 0; offset < output.size();) {
                    const std::size_t count = std::min(block_size, output.size() - offset);
                    using accumulator_type = std::conditional_t<Direct, fast_params::limb_array, integer_type>;
                    std::array<accumulator_type, block_size> sums {};
                    std::array<std::uint64_t, block_size> lower {};
                    std::array<std::array<std::uint64_t, block_size>, small_prime_fft_primes.size()> saved_digits;
                    for (std::size_t prime = 0; prime < primes_.size(); ++prime) {
                        const auto &context = primes_[prime];
                        const auto &arithmetic = context.arithmetic;
                        const double modulus = static_cast<double>(arithmetic.modulus());
                        const double *source =
                            residues.data() + (prime * field_type::arity + coordinate) * stride + offset;
                        const auto accumulate = [&](std::size_t i, std::uint64_t digit) {
                            if constexpr (Direct) {
                                using wide_type = algebra::fields::detail::fp12_fast::wide_limb;
                                saved_digits[prime][i] = digit;
                                lower[i] += static_cast<std::uint64_t>(
                                    (wide_type(digit) * context.certificate_reciprocal) >> 50);
                                add_export_product(sums[i], context.export_weight, digit);
                            } else {
                                add_cofactor_product(sums[i], context.cofactor, digit);
                            }
                        };
                        std::size_t i = 0;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                        if constexpr (UseFma) {
                            if (count >= 4) {
                                const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic);
                                const __m256d inverse = _mm256_set1_pd(context.inverse);
                                const __m256d q = _mm256_set1_pd(modulus);
                                for (; count - i >= 4; i += 4) {
                                    __m256d digits = packed.reduce_fma(_mm256_loadu_pd(source + i));
                                    if constexpr (!Weighted) {
                                        digits = packed.multiply_fma(digits, inverse);
                                    }
                                    const __m256d negative = _mm256_cmp_pd(digits, _mm256_setzero_pd(), _CMP_LT_OQ);
                                    digits = _mm256_add_pd(digits, _mm256_and_pd(negative, q));
                                    std::array<double, 4> canonical;
                                    _mm256_storeu_pd(canonical.data(), digits);
                                    for (std::size_t lane = 0; lane < canonical.size(); ++lane) {
                                        accumulate(i + lane, static_cast<std::uint64_t>(canonical[lane]));
                                    }
                                }
                            }
                        }
#endif
                        for (; i < count; ++i) {
                            double digit;
                            if constexpr (UseFma) {
                                digit = arithmetic.reduce_fma(source[i]);
                                if constexpr (!Weighted) {
                                    digit = arithmetic.multiply_fma(digit, context.inverse);
                                }
                            } else {
                                digit = arithmetic.reduce_integer(source[i]);
                                if constexpr (!Weighted) {
                                    digit = arithmetic.multiply_integer(digit, context.inverse);
                                }
                            }
                            // Both operands have magnitude < q, so the product remainder is in (-q,q).
                            // Canonical digits give 0 <= digit*(M/q) < M and bound the complete sum.
                            const auto canonical = static_cast<std::uint64_t>(digit < 0 ? digit + modulus : digit);
                            accumulate(i, canonical);
                        }
                    }
                    for (std::size_t i = 0; i < count; ++i) {
                        auto &destination = output[offset + i].coordinate(coordinate);
                        if constexpr (Direct) {
                            std::size_t quotient;
                            if (certify_export(lower[i], quotient)) {
                                destination = export_direct(sums[i], quotient);
                            } else {
                                integer_type exact_sum = 0;
                                for (std::size_t prime = 0; prime < primes_.size(); ++prime) {
                                    add_cofactor_product(exact_sum, primes_[prime].cofactor, saved_digits[prime][i]);
                                }
                                destination = export_exact(exact_sum);
                            }
                        } else {
                            destination = export_exact(sums[i]);
                        }
                    }
                    offset += count;
                }
            }
        }

#if defined(CRYPTO3_SMALL_PRIME_CRT_IFMA_SELECTED)
        template<bool Weighted>
        void reconstruct_ifma(std::span<value_type> output, std::span<const double> residues,
                              std::size_t stride) const {
            using operations = algebra::fields::detail::fp12_fast::ifma52_operations<8>;
            for (std::size_t coordinate = 0; coordinate < field_type::arity; ++coordinate) {
                for (std::size_t offset = 0; offset < output.size();) {
                    const std::size_t count = std::min(operations::lanes, output.size() - offset);
                    const auto active = static_cast<__mmask8>((1u << count) - 1);
                    using accumulator_type = small_prime_fft_crt_accumulator<fast_params, operations>;
                    accumulator_type accumulator;
                    auto lower = operations::broadcast(0);
                    std::array<std::array<std::uint64_t, operations::lanes>, small_prime_fft_primes.size()>
                        saved_digits;
                    for (std::size_t prime = 0; prime < primes_.size(); ++prime) {
                        const auto &context = primes_[prime];
                        const small_prime_fft_arithmetic::wide_fma_kernel packed(context.arithmetic);
                        const double *source =
                            residues.data() + (prime * field_type::arity + coordinate) * stride + offset;
                        auto digits = packed.reduce_fma(_mm512_maskz_loadu_pd(active, source));
                        if constexpr (!Weighted) {
                            digits = packed.multiply_fma(digits, _mm512_set1_pd(context.inverse));
                        }
                        const auto negative = _mm512_cmp_pd_mask(digits, _mm512_setzero_pd(), _CMP_LT_OQ);
                        digits = _mm512_mask_add_pd(digits, negative, digits,
                                                    _mm512_set1_pd(static_cast<double>(context.arithmetic.modulus())));
                        const operations::value_type canonical {_mm512_cvttpd_epu64(digits)};
                        operations::store(saved_digits[prime].data(), canonical);
                        accumulator.add(context.export_weight52, canonical);
                        lower = operations::add(
                            lower, accumulator_type::certificate_term(canonical, context.certificate_reciprocal));
                    }
                    const auto certificate = accumulator_type::certify(lower, scaled_bound_, primes_.size());
                    const auto accepted = operations::less(operations::broadcast(0), certificate.accepted) & active;
                    std::array<std::array<std::uint64_t, operations::lanes>, 5> reduced;
                    if (accepted != 0) {
                        std::array<operations::value_type, 5> corrections;
                        for (std::size_t digit = 0; digit < corrections.size(); ++digit) {
                            // Every quotient is in [0,prime_count], including uncertified lanes.
                            corrections[digit] = {_mm512_i64gather_epi64(certificate.quotients.data,
                                                                         export_corrections52_[digit].data(), 8)};
                        }
                        const auto values = accumulator.reduce(corrections);
                        for (std::size_t digit = 0; digit < values.size(); ++digit)
                            operations::store(reduced[digit].data(), values[digit]);
                    }
                    for (std::size_t lane = 0; lane < count; ++lane) {
                        auto &destination = output[offset + lane].coordinate(coordinate);
                        if ((accepted & (1u << lane)) != 0) {
                            // Numeric repacking only: the R*S weights already account for S/R.
                            base_field_type::value_type value;
                            auto *words = value.data.backend().base_data().limbs();
                            words[0] = reduced[0][lane] | (reduced[1][lane] << 52);
                            words[1] = (reduced[1][lane] >> 12) | (reduced[2][lane] << 40);
                            words[2] = (reduced[2][lane] >> 24) | (reduced[3][lane] << 28);
                            words[3] = (reduced[3][lane] >> 36) | (reduced[4][lane] << 16);
                            destination = value;
                        } else {
                            integer_type exact_sum = 0;
                            for (std::size_t prime = 0; prime < primes_.size(); ++prime) {
                                add_cofactor_product(exact_sum, primes_[prime].cofactor, saved_digits[prime][lane]);
                            }
                            destination = export_exact(exact_sum);
                        }
                    }
                    offset += count;
                }
            }
        }
#endif

        std::vector<prime_context> primes_;
        std::array<integer_type, correction_count> multiples_;
        integer_type half_modulus_;
        integer_type bound_;
        std::size_t accumulator_limbs_;
        export_type export_radix_;
        std::uint64_t scaled_bound_;
        std::array<direct_weight_type, small_prime_fft_primes.size() + 1> export_corrections_ {};
#if defined(CRYPTO3_SMALL_PRIME_CRT_IFMA_SELECTED)
        std::array<std::array<std::uint64_t, small_prime_fft_primes.size() + 1>, 5> export_corrections52_ {};
#endif
    };

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_DETAIL_BN254_SMALL_PRIME_FFT_CRT_HPP
