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
#include <type_traits>
#include <vector>

#include <boost/config.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/limb_ops.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/math/detail/small_prime_fft_arithmetic.hpp>
#include <nil/crypto3/math/detail/small_prime_fft_params.hpp>
#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>
#include <nil/crypto3/multiprecision/inverse.hpp>

#ifndef CRYPTO3_SMALL_PRIME_FFT_RADIX
#define CRYPTO3_SMALL_PRIME_FFT_RADIX 4
#endif
static_assert(CRYPTO3_SMALL_PRIME_FFT_RADIX == 4 || CRYPTO3_SMALL_PRIME_FFT_RADIX == 8 ||
                  CRYPTO3_SMALL_PRIME_FFT_RADIX == 16,
              "small-prime FFT stage radix must be 4, 8 or 16");

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
        template<bool UseFma, typename Params = field_type::extension_policy::fast_params>
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
                        using base_backend_type = typename Params::base_field_type::integral_type::backend_type;
                        if constexpr (base_backend_type::limb_bits == 64) {
                            namespace fp12_fast = algebra::fields::detail::fp12_fast;
                            // The stored value is a*R mod p, with R = 2^256. Zero-extend it
                            // into the low half of montgomery; one reduction recovers a in [0,p).
                            typename Params::limb_array montgomery {};
                            const auto stored = fp12_fast::load_limbs<base_backend_type>(
                                coefficients[offset + i].coordinate(coordinate).data.backend().base_data());
                            std::copy(stored.begin(), stored.end(), montgomery.begin());
                            std::array<std::uint64_t, Params::base_value_limb_count> canonical;
                            fp12_fast::montgomery_reduce<Params>(canonical.data(), montgomery);
                            // pack only extracts five 52-bit chunks; it does not rescale them.
                            // Packing the canonical value preserves the lift of field -1 to p-1.
                            words = fp12_fast::montgomery_radix52<Params>::pack(canonical.data());
                        } else {
                            const auto canonical = coefficients[offset + i].coordinate(coordinate).to_integral();
                            boost::multiprecision::export_bits(canonical, words.begin(), chunk_bits, false);
                        }
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

    /**
     * In-place, power-of-two transforms for one prime from the catalogue. The forward transform
     * maps natural-order coefficients to evaluations at omega^k in bit-reversed index order.
     * Both directions use this order, so convolution needs no permutation passes.
     *
     * A reusable table of max(1,N/2) centered twiddles serves both directions. Radix-4 recursion
     * fuses pairs of layers, with radix-2 leaves for odd depths and special cases for unit twiddles.
     * AVX2/FMA builds fuse 8-, 16- and 32-point leaves; AVX-512 adds paired 32-/64-point
     * leaves and eight-lane stages. CRYPTO3_SMALL_PRIME_FFT_RADIX selects 4 (default), 8
     * or 16 for large SIMD stages; the scalar/integer fallback keeps the radix-4 schedule.
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

        /**
         * Unscaled inverse for pointwise products whose integer representatives are in (-q,q).
         * The caller has already folded 1/N into its pointwise scale. This stricter input bound
         * permits omitting the first leaf reductions; the output still lies in (-2*q,2*q).
         * Dimension errors throw before writes, as in inverse_fft().
         */
        void inverse_fft_pointwise(std::span<double> values) const {
            check_buffer<1>(values);
            if (arithmetic_.fma_supported()) {
                inverse_impl<true, true>(values, false);
            } else {
                inverse_impl<false, true>(values, false);
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

        template<unsigned Bound = 3>
        void check_buffer(std::span<const double> values) const {
            if (values.size() != size_ || twiddles_.empty()) {
                throw std::invalid_argument("small_prime_fft_plan: invalid buffer dimensions or moved-from plan");
            }
#ifndef NDEBUG
            const double limit = Bound * static_cast<double>(modulus());
            for (double value : values) {
                assert(std::isfinite(value) && std::trunc(value) == value && std::abs(value) < limit);
            }
#endif
        }

        double centered(double value) const {
            const double q = static_cast<double>(modulus());
            return value > q / 2 ? value - q : (value < -q / 2 ? value + q : value);
        }

        BOOST_FORCEINLINE double inverse_twiddle(std::size_t index) const {
            if (index == 0) {
                return 1.0;
            }
            // Within [b,2*b), reflected exponents add to half the root order, hence w^-1=-w_reflected.
            const std::size_t base = std::bit_floor(index);
            return -twiddles_[base + (2 * base - 1 - index)];
        }

        // Keep arithmetic selection outside the butterfly loops, with compile-time dispatch inside them.
        template<bool UseFma, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static Value reduce(const Arithmetic &arithmetic, Value value) {
            if constexpr (UseFma) {
                return arithmetic.reduce_fma(value);
            } else {
                return arithmetic.reduce_integer(value);
            }
        }

        template<bool UseFma, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static Value multiply(const Arithmetic &arithmetic, Value left, Value right) {
            if constexpr (UseFma) {
                return arithmetic.multiply_fma(left, right);
            } else {
                return arithmetic.multiply_integer(left, right);
            }
        }

        template<bool UseFma, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void multiply_pair(const Arithmetic &arithmetic, Value &a, Value &b, Value wa,
                                                    Value wb) {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma && requires { arithmetic.multiply_fma_pair(a, b, wa, wb); }) {
                arithmetic.multiply_fma_pair(a, b, wa, wb);
            } else
#endif
            {
                a = multiply<UseFma>(arithmetic, a, wa);
                b = multiply<UseFma>(arithmetic, b, wb);
            }
        }

        template<bool UseFma, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void forward_butterfly(const Arithmetic &arithmetic, Value &v0, Value &v1, Value &v2,
                                                        Value &v3, Value w2, Value w, Value iw) {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
            const Value a0 = reduce<UseFma>(arithmetic, v0);
            const Value a1 = v1;
            Value a2 = v2;
            Value a3 = v3;
            if constexpr (UnitTwiddle) {
                a2 = reduce<UseFma>(arithmetic, a2);
                a3 = reduce<UseFma>(arithmetic, a3);
            } else {
                multiply_pair<UseFma>(arithmetic, a2, a3, w2, w2);
            }
            const Value even_sum = a0 + a2;
            const Value even_difference = a0 - a2;
            Value odd_sum = a1 + a3;
            Value odd_difference = a1 - a3;
            if constexpr (UnitTwiddle) {
                odd_sum = reduce<UseFma>(arithmetic, odd_sum);
                odd_difference = multiply<UseFma>(arithmetic, odd_difference, iw);
            } else {
                multiply_pair<UseFma>(arithmetic, odd_sum, odd_difference, w, iw);
            }
            // Input magnitudes < 3*q and weights < q/2 imply |a2|,|a3| < q and
            // |a1 +/- a3| < 4*q. Every product is < 2*q^2, so the odd terms are < q.
            // Thus all outputs are < 3*q, with every integer sum exact below 2^53.
            v0 = even_sum + odd_sum;
            v1 = even_sum - odd_sum;
            v2 = even_difference + odd_difference;
            v3 = even_difference - odd_difference;
        }

#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
        BOOST_FORCEINLINE static void transpose4(__m256d &v0, __m256d &v1, __m256d &v2, __m256d &v3) {
            const __m256d a = _mm256_unpacklo_pd(v0, v1);
            const __m256d b = _mm256_unpackhi_pd(v0, v1);
            const __m256d c = _mm256_unpacklo_pd(v2, v3);
            const __m256d d = _mm256_unpackhi_pd(v2, v3);
            v0 = _mm256_permute2f128_pd(a, c, 0x20);
            v1 = _mm256_permute2f128_pd(b, d, 0x20);
            v2 = _mm256_permute2f128_pd(a, c, 0x31);
            v3 = _mm256_permute2f128_pd(b, d, 0x31);
        }

        // Aligned table indices, with unaligned loads: the owning vector needs no special allocator.
        template<bool Inverse>
        BOOST_FORCEINLINE __m256d twiddle4(std::size_t index) const {
            if constexpr (!Inverse) {
                return _mm256_loadu_pd(twiddles_.data() + index);
            } else {
                if (index == 0) {
                    return _mm256_setr_pd(1.0, inverse_twiddle(1), inverse_twiddle(2), inverse_twiddle(3));
                }
                // Four consecutive entries stay in one [b,2*b) segment. Reflect, reverse, negate.
                const std::size_t base = std::bit_floor(index);
                const __m256d reflected = _mm256_loadu_pd(twiddles_.data() + 3 * base - index - 4);
                return _mm256_sub_pd(_mm256_setzero_pd(), _mm256_permute4x64_pd(reflected, 0x1b));
            }
        }

        template<bool Inverse>
        BOOST_FORCEINLINE void leaf_twiddles(std::size_t node, __m256d &w2, __m256d &w, __m256d &iw) const {
            w2 = twiddle4<Inverse>(4 * node);
            const __m256d low = twiddle4<Inverse>(8 * node);
            const __m256d high = twiddle4<Inverse>(8 * node + 4);
            w = _mm256_permute4x64_pd(_mm256_unpacklo_pd(low, high), 0xd8);
            iw = _mm256_permute4x64_pd(_mm256_unpackhi_pd(low, high), 0xd8);
        }

        template<bool UnitTwiddle>
        BOOST_FORCEINLINE void forward16(const small_prime_fft_arithmetic::packed_fma_kernel &packed, __m256d &v0,
                                         __m256d &v1, __m256d &v2, __m256d &v3, std::size_t node) const {
            forward_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, _mm256_set1_pd(twiddles_[node]),
                                                 _mm256_set1_pd(twiddles_[2 * node]),
                                                 _mm256_set1_pd(twiddles_[2 * node + 1]));
            // The four child transforms now run across SIMD lanes, without scalar leaves or reloads.
            transpose4(v0, v1, v2, v3);
            __m256d w2, w, iw;
            leaf_twiddles<false>(node, w2, w, iw);
            forward_butterfly<true, false>(packed, v0, v1, v2, v3, w2, w, iw);
            transpose4(v0, v1, v2, v3);
        }

        template<bool UnitTwiddle, bool BoundedInput>
        BOOST_FORCEINLINE void inverse16(const small_prime_fft_arithmetic::packed_fma_kernel &packed, __m256d &v0,
                                         __m256d &v1, __m256d &v2, __m256d &v3, std::size_t node) const {
            transpose4(v0, v1, v2, v3);
            __m256d w2, w, iw;
            leaf_twiddles<true>(node, w2, w, iw);
            inverse_butterfly<true, false>(packed, v0, v1, v2, v3, w2, w, iw, !BoundedInput);
            transpose4(v0, v1, v2, v3);
            inverse_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, _mm256_set1_pd(inverse_twiddle(node)),
                                                 _mm256_set1_pd(inverse_twiddle(2 * node)),
                                                 _mm256_set1_pd(inverse_twiddle(2 * node + 1)), false);
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput = false>
        void base16(double *values, std::size_t node) const {
            const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
            __m256d v0 = _mm256_loadu_pd(values);
            __m256d v1 = _mm256_loadu_pd(values + 4);
            __m256d v2 = _mm256_loadu_pd(values + 8);
            __m256d v3 = _mm256_loadu_pd(values + 12);
            if constexpr (Inverse) {
                inverse16<UnitTwiddle, BoundedInput>(packed, v0, v1, v2, v3, node);
            } else {
                forward16<UnitTwiddle>(packed, v0, v1, v2, v3, node);
            }
            _mm256_storeu_pd(values, v0);
            _mm256_storeu_pd(values + 4, v1);
            _mm256_storeu_pd(values + 8, v2);
            _mm256_storeu_pd(values + 12, v3);
        }

        template<bool Inverse, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void radix2_packed(const Arithmetic &packed, Value &left, Value &right, Value weight) {
            if constexpr (Inverse) {
                // Both subtransforms already return magnitudes < 2*q. No input reductions are needed.
                const Value sum = left + right;
                const Value difference = left - right;
                left = packed.reduce_fma(sum);
                right = UnitTwiddle ? packed.reduce_fma(difference) : packed.multiply_fma(difference, weight);
            } else {
                const Value a = packed.reduce_fma(left);
                const Value b = UnitTwiddle ? packed.reduce_fma(right) : packed.multiply_fma(right, weight);
                left = a + b;
                right = a - b;
            }
        }

        template<bool Inverse, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void stage4(const Arithmetic &arithmetic, Value &v0, Value &v1, Value &v2, Value &v3,
                                             const double *weights) {
            const Value w2 = Arithmetic::broadcast(weights[0]);
            const Value w = Arithmetic::broadcast(weights[1]);
            const Value iw = Arithmetic::broadcast(weights[2]);
            if constexpr (Inverse) {
                inverse_butterfly<true, UnitTwiddle>(arithmetic, v0, v1, v2, v3, w2, w, iw, false);
            } else {
                forward_butterfly<true, UnitTwiddle>(arithmetic, v0, v1, v2, v3, w2, w, iw);
            }
        }

        template<bool Inverse>
        BOOST_FORCEINLINE void stage_weights(double *weights, std::size_t node) const {
            weights[0] = Inverse ? inverse_twiddle(node) : twiddles_[node];
            weights[1] = Inverse ? inverse_twiddle(2 * node) : twiddles_[2 * node];
            weights[2] = Inverse ? inverse_twiddle(2 * node + 1) : twiddles_[2 * node + 1];
        }

        template<bool Inverse, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE void stage4(const Arithmetic &arithmetic, Value &v0, Value &v1, Value &v2, Value &v3,
                                      std::size_t node) const {
            double weights[3];
            stage_weights<Inverse>(weights, node);
            stage4<Inverse, UnitTwiddle>(arithmetic, v0, v1, v2, v3, weights);
        }

        template<unsigned Radix, bool Inverse, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void fused_outer(const Arithmetic &arithmetic, Value *v, const double *weights) {
            if constexpr (Radix == 8) {
                const Value weight = Arithmetic::broadcast(weights[0]);
                radix2_packed<Inverse, UnitTwiddle>(arithmetic, v[0], v[4], weight);
                radix2_packed<Inverse, UnitTwiddle>(arithmetic, v[1], v[5], weight);
                radix2_packed<Inverse, UnitTwiddle>(arithmetic, v[2], v[6], weight);
                radix2_packed<Inverse, UnitTwiddle>(arithmetic, v[3], v[7], weight);
            } else {
                stage4<Inverse, UnitTwiddle>(arithmetic, v[0], v[4], v[8], v[12], weights);
                stage4<Inverse, UnitTwiddle>(arithmetic, v[1], v[5], v[9], v[13], weights);
                stage4<Inverse, UnitTwiddle>(arithmetic, v[2], v[6], v[10], v[14], weights);
                stage4<Inverse, UnitTwiddle>(arithmetic, v[3], v[7], v[11], v[15], weights);
            }
        }

        // Fuse three/four binary layers without intermediate array passes. Every constituent
        // butterfly retains the forward <3*q and inverse <2*q bounds; no wider lazy range is used.
        // UnitTwiddle is valid only at node zero, hence only the first child inherits it.
        template<unsigned Radix, bool Inverse, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void fused_stage(const Arithmetic &arithmetic, Value *v, const double *weights) {
            static_assert(Radix == 8 || Radix == 16);
            if constexpr (!Inverse) {
                fused_outer<Radix, Inverse, UnitTwiddle>(arithmetic, v, weights);
            }
            constexpr std::size_t outer_weights = Radix == 8 ? 1 : 3;
            stage4<Inverse, UnitTwiddle>(arithmetic, v[0], v[1], v[2], v[3], weights + outer_weights);
            stage4<Inverse, false>(arithmetic, v[4], v[5], v[6], v[7], weights + outer_weights + 3);
            if constexpr (Radix == 16) {
                stage4<Inverse, false>(arithmetic, v[8], v[9], v[10], v[11], weights + outer_weights + 6);
                stage4<Inverse, false>(arithmetic, v[12], v[13], v[14], v[15], weights + outer_weights + 9);
            }
            if constexpr (Inverse) {
                fused_outer<Radix, Inverse, UnitTwiddle>(arithmetic, v, weights);
            }
        }

        template<unsigned Radix, bool Inverse, bool UnitTwiddle>
        void large_radix(double *values, std::size_t stride, std::size_t node) const {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
            using arithmetic_type = small_prime_fft_arithmetic::wide_fma_kernel;
#else
            using arithmetic_type = small_prime_fft_arithmetic::packed_fma_kernel;
#endif
            const arithmetic_type arithmetic(arithmetic_);
            // Dispatch requires stride >=16, a power of two divisible by the SIMD width.
            assert(stride >= 16 && stride % arithmetic_type::lanes == 0);
            // The node is fixed throughout this coefficient loop. Resolve reflections and
            // table lookups once, retaining scalar weights to avoid forcing Radix-1 live vectors.
            double weights[Radix - 1];
            if constexpr (Radix == 8) {
                weights[0] = Inverse ? inverse_twiddle(node) : twiddles_[node];
            } else {
                stage_weights<Inverse>(weights, node);
            }
            constexpr std::size_t outer_weights = Radix == 8 ? 1 : 3;
#pragma GCC unroll 4
            for (std::size_t child = 0; child < Radix / 4; ++child) {
                stage_weights<Inverse>(weights + outer_weights + 3 * child, Radix / 4 * node + child);
            }
            for (std::size_t i = 0; i < stride; i += arithmetic_type::lanes) {
                arithmetic_type::value_type v[Radix];
#pragma GCC unroll 16
                for (std::size_t part = 0; part < Radix; ++part) {
                    v[part] = arithmetic_type::load(values + part * stride + i);
                }
                fused_stage<Radix, Inverse, UnitTwiddle>(arithmetic, v, weights);
#pragma GCC unroll 16
                for (std::size_t part = 0; part < Radix; ++part) {
                    arithmetic_type::store(values + part * stride + i, v[part]);
                }
            }
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput = false>
        void base32(double *values, std::size_t node) const {
            const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
            __m256d v0 = _mm256_loadu_pd(values);
            __m256d v1 = _mm256_loadu_pd(values + 4);
            __m256d v2 = _mm256_loadu_pd(values + 8);
            __m256d v3 = _mm256_loadu_pd(values + 12);
            __m256d v4 = _mm256_loadu_pd(values + 16);
            __m256d v5 = _mm256_loadu_pd(values + 20);
            __m256d v6 = _mm256_loadu_pd(values + 24);
            __m256d v7 = _mm256_loadu_pd(values + 28);
            const __m256d weight = _mm256_set1_pd(Inverse ? inverse_twiddle(node) : twiddles_[node]);
            // One radix-2 layer and two register-resident 16-point transforms cover all five layers.
            if constexpr (Inverse) {
                inverse16<UnitTwiddle, BoundedInput>(packed, v0, v1, v2, v3, 2 * node);
                inverse16<false, BoundedInput>(packed, v4, v5, v6, v7, 2 * node + 1);
            }
            radix2_packed<Inverse, UnitTwiddle>(packed, v0, v4, weight);
            radix2_packed<Inverse, UnitTwiddle>(packed, v1, v5, weight);
            radix2_packed<Inverse, UnitTwiddle>(packed, v2, v6, weight);
            radix2_packed<Inverse, UnitTwiddle>(packed, v3, v7, weight);
            if constexpr (!Inverse) {
                forward16<UnitTwiddle>(packed, v0, v1, v2, v3, 2 * node);
                forward16<false>(packed, v4, v5, v6, v7, 2 * node + 1);
            }
            _mm256_storeu_pd(values, v0);
            _mm256_storeu_pd(values + 4, v1);
            _mm256_storeu_pd(values + 8, v2);
            _mm256_storeu_pd(values + 12, v3);
            _mm256_storeu_pd(values + 16, v4);
            _mm256_storeu_pd(values + 20, v5);
            _mm256_storeu_pd(values + 24, v6);
            _mm256_storeu_pd(values + 28, v7);
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput = false>
        void base8(double *values, std::size_t node) const {
            const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
            const __m256d low = _mm256_loadu_pd(values);
            const __m256d high = _mm256_loadu_pd(values + 4);
            __m256d v0, v1, v2, v3;
            if constexpr (Inverse) {
                __m256d a = _mm256_permute4x64_pd(_mm256_unpacklo_pd(low, high), 0xd8);
                __m256d b = _mm256_permute4x64_pd(_mm256_unpackhi_pd(low, high), 0xd8);
                if constexpr (!BoundedInput) {
                    a = packed.reduce_fma(a);
                    b = packed.reduce_fma(b);
                }
                const __m256d sum = packed.reduce_fma(a + b);
                const __m256d difference = packed.multiply_fma(a - b, twiddle4<true>(4 * node));
                const __m256d even = _mm256_unpacklo_pd(sum, difference);
                const __m256d odd = _mm256_unpackhi_pd(sum, difference);
                v0 = _mm256_permute2f128_pd(even, even, 0x00);
                v1 = _mm256_permute2f128_pd(odd, odd, 0x00);
                v2 = _mm256_permute2f128_pd(even, even, 0x11);
                v3 = _mm256_permute2f128_pd(odd, odd, 0x11);
                inverse_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, _mm256_set1_pd(inverse_twiddle(node)),
                                                     _mm256_set1_pd(inverse_twiddle(2 * node)),
                                                     _mm256_set1_pd(inverse_twiddle(2 * node + 1)), false);
                _mm256_storeu_pd(values, _mm256_permute2f128_pd(v0, v1, 0x20));
                _mm256_storeu_pd(values + 4, _mm256_permute2f128_pd(v2, v3, 0x20));
            } else {
                // Two independent columns occupy the lower lanes; the duplicated upper lanes are discarded.
                v0 = _mm256_permute2f128_pd(low, low, 0x00);
                v1 = _mm256_permute2f128_pd(low, low, 0x11);
                v2 = _mm256_permute2f128_pd(high, high, 0x00);
                v3 = _mm256_permute2f128_pd(high, high, 0x11);
                forward_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, _mm256_set1_pd(twiddles_[node]),
                                                     _mm256_set1_pd(twiddles_[2 * node]),
                                                     _mm256_set1_pd(twiddles_[2 * node + 1]));
                const __m256d a = packed.reduce_fma(
                    _mm256_permute2f128_pd(_mm256_unpacklo_pd(v0, v1), _mm256_unpacklo_pd(v2, v3), 0x20));
                const __m256d b = packed.multiply_fma(
                    _mm256_permute2f128_pd(_mm256_unpackhi_pd(v0, v1), _mm256_unpackhi_pd(v2, v3), 0x20),
                    twiddle4<false>(4 * node));
                const __m256d even = _mm256_unpacklo_pd(a + b, a - b);
                const __m256d odd = _mm256_unpackhi_pd(a + b, a - b);
                _mm256_storeu_pd(values, _mm256_permute2f128_pd(even, odd, 0x20));
                _mm256_storeu_pd(values + 4, _mm256_permute2f128_pd(even, odd, 0x31));
            }
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput = false>
        void basecase(double *values, std::size_t length, std::size_t node) const {
            if (length == 16) {
                base16<Inverse, UnitTwiddle, BoundedInput>(values, node);
            } else if (length == 32) {
                base32<Inverse, UnitTwiddle, BoundedInput>(values, node);
            } else {
                base8<Inverse, UnitTwiddle, BoundedInput>(values, node);
            }
        }

#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
        BOOST_FORCEINLINE static __m512d join4(__m256d low, __m256d high) {
            return _mm512_insertf64x4(_mm512_castpd256_pd512(low), high, 1);
        }

        // Two independent 4x4 transposes, one in each 256-bit half. The stored spectrum
        // remains bit-reversed; these permutations only change the lanes during a leaf.
        BOOST_FORCEINLINE static void transpose4(__m512d &v0, __m512d &v1, __m512d &v2, __m512d &v3) {
            const __m512d a = _mm512_unpacklo_pd(v0, v1);
            const __m512d b = _mm512_unpackhi_pd(v0, v1);
            const __m512d c = _mm512_unpacklo_pd(v2, v3);
            const __m512d d = _mm512_unpackhi_pd(v2, v3);
            const __m512i low = _mm512_setr_epi64(0, 1, 8, 9, 4, 5, 12, 13);
            const __m512i high = _mm512_setr_epi64(2, 3, 10, 11, 6, 7, 14, 15);
            v0 = _mm512_permutex2var_pd(a, low, c);
            v1 = _mm512_permutex2var_pd(b, low, d);
            v2 = _mm512_permutex2var_pd(a, high, c);
            v3 = _mm512_permutex2var_pd(b, high, d);
        }

        template<bool Inverse>
        BOOST_FORCEINLINE __m512d paired_twiddle(std::size_t low, std::size_t high) const {
            return join4(_mm256_set1_pd(Inverse ? inverse_twiddle(low) : twiddles_[low]),
                         _mm256_set1_pd(Inverse ? inverse_twiddle(high) : twiddles_[high]));
        }

        // Keep the complete leaf in one function: outlining this step forces live vectors
        // through memory across the call. Recursive transform boundaries remain uninlined.
        template<bool Inverse>
        BOOST_FORCEINLINE void paired_outer(const small_prime_fft_arithmetic::wide_fma_kernel &arithmetic, __m512d &v0,
                                            __m512d &v1, __m512d &v2, __m512d &v3, std::size_t node) const {
            const __m512d w2 = paired_twiddle<Inverse>(node, node + 1);
            const __m512d w = paired_twiddle<Inverse>(2 * node, 2 * node + 2);
            const __m512d iw = paired_twiddle<Inverse>(2 * node + 1, 2 * node + 3);
            if constexpr (Inverse) {
                inverse_butterfly<true, false>(arithmetic, v0, v1, v2, v3, w2, w, iw, false);
            } else {
                forward_butterfly<true, false>(arithmetic, v0, v1, v2, v3, w2, w, iw);
            }
        }

        template<bool Inverse, bool BoundedInput>
        BOOST_FORCEINLINE void base16_pair(const small_prime_fft_arithmetic::wide_fma_kernel &arithmetic, __m512d &v0,
                                           __m512d &v1, __m512d &v2, __m512d &v3, std::size_t node) const {
            // Group matching rows of two adjacent 16-point transforms in the two vector halves.
            const __m512d a = _mm512_shuffle_f64x2(v0, v2, 0x44);
            const __m512d b = _mm512_shuffle_f64x2(v0, v2, 0xee);
            const __m512d c = _mm512_shuffle_f64x2(v1, v3, 0x44);
            const __m512d d = _mm512_shuffle_f64x2(v1, v3, 0xee);
            v0 = a;
            v1 = b;
            v2 = c;
            v3 = d;
            if constexpr (!Inverse) {
                paired_outer<false>(arithmetic, v0, v1, v2, v3, node);
            }
            transpose4(v0, v1, v2, v3);
            __m256d lo_w2, lo_w, lo_iw, hi_w2, hi_w, hi_iw;
            leaf_twiddles<Inverse>(node, lo_w2, lo_w, lo_iw);
            leaf_twiddles<Inverse>(node + 1, hi_w2, hi_w, hi_iw);
            if constexpr (Inverse) {
                inverse_butterfly<true, false>(arithmetic, v0, v1, v2, v3, join4(lo_w2, hi_w2), join4(lo_w, hi_w),
                                               join4(lo_iw, hi_iw), !BoundedInput);
            } else {
                forward_butterfly<true, false>(arithmetic, v0, v1, v2, v3, join4(lo_w2, hi_w2), join4(lo_w, hi_w),
                                               join4(lo_iw, hi_iw));
            }
            transpose4(v0, v1, v2, v3);
            if constexpr (Inverse) {
                paired_outer<true>(arithmetic, v0, v1, v2, v3, node);
            }
            // Restore contiguous bit-reversed order in the same four input/output vectors.
            const __m512d out0 = _mm512_shuffle_f64x2(v0, v1, 0x44);
            const __m512d out1 = _mm512_shuffle_f64x2(v2, v3, 0x44);
            const __m512d out2 = _mm512_shuffle_f64x2(v0, v1, 0xee);
            const __m512d out3 = _mm512_shuffle_f64x2(v2, v3, 0xee);
            v0 = out0;
            v1 = out1;
            v2 = out2;
            v3 = out3;
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput>
        void base32_wide(double *values, std::size_t node) const {
            const small_prime_fft_arithmetic::wide_fma_kernel arithmetic(arithmetic_);
            __m512d v0 = _mm512_loadu_pd(values);
            __m512d v1 = _mm512_loadu_pd(values + 8);
            __m512d v2 = _mm512_loadu_pd(values + 16);
            __m512d v3 = _mm512_loadu_pd(values + 24);
            const __m512d weight = _mm512_set1_pd(Inverse ? inverse_twiddle(node) : twiddles_[node]);
            if constexpr (Inverse) {
                base16_pair<true, BoundedInput>(arithmetic, v0, v1, v2, v3, 2 * node);
            }
            radix2_packed<Inverse, UnitTwiddle>(arithmetic, v0, v2, weight);
            radix2_packed<Inverse, UnitTwiddle>(arithmetic, v1, v3, weight);
            if constexpr (!Inverse) {
                base16_pair<false, false>(arithmetic, v0, v1, v2, v3, 2 * node);
            }
            _mm512_storeu_pd(values, v0);
            _mm512_storeu_pd(values + 8, v1);
            _mm512_storeu_pd(values + 16, v2);
            _mm512_storeu_pd(values + 24, v3);
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput>
        void base64_wide(double *values, std::size_t node) const {
            const small_prime_fft_arithmetic::wide_fma_kernel arithmetic(arithmetic_);
            __m512d v[8];
#pragma GCC unroll 8
            for (std::size_t i = 0; i < 8; ++i) {
                v[i] = _mm512_loadu_pd(values + 8 * i);
            }
            if constexpr (!Inverse) {
                stage4<false, UnitTwiddle>(arithmetic, v[0], v[2], v[4], v[6], node);
                stage4<false, UnitTwiddle>(arithmetic, v[1], v[3], v[5], v[7], node);
            }
            base16_pair<Inverse, BoundedInput>(arithmetic, v[0], v[1], v[2], v[3], 4 * node);
            base16_pair<Inverse, BoundedInput>(arithmetic, v[4], v[5], v[6], v[7], 4 * node + 2);
            if constexpr (Inverse) {
                stage4<true, UnitTwiddle>(arithmetic, v[0], v[2], v[4], v[6], node);
                stage4<true, UnitTwiddle>(arithmetic, v[1], v[3], v[5], v[7], node);
            }
#pragma GCC unroll 8
            for (std::size_t i = 0; i < 8; ++i) {
                _mm512_storeu_pd(values + 8 * i, v[i]);
            }
        }

        template<bool Inverse, bool BoundedInput>
        void wide_basecase(double *values, std::size_t length, std::size_t node) const {
            if (length == 32) {
                if (node == 0) {
                    base32_wide<Inverse, true, BoundedInput>(values, node);
                } else {
                    base32_wide<Inverse, false, BoundedInput>(values, node);
                }
            } else {
                if (node == 0) {
                    base64_wide<Inverse, true, BoundedInput>(values, node);
                } else {
                    base64_wide<Inverse, false, BoundedInput>(values, node);
                }
            }
        }
#endif
#endif

        template<bool UseFma, bool UnitTwiddle>
        void forward_radix4(double *values, std::size_t quarter, std::size_t node, std::size_t first,
                            std::size_t last) const {
            const double w2 = twiddles_[node];
            const double w = twiddles_[2 * node];
            const double iw = twiddles_[2 * node + 1];
            std::size_t i = first;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
            if constexpr (UseFma) {
                if (quarter >= 8) {
                    const small_prime_fft_arithmetic::wide_fma_kernel packed(arithmetic_);
                    const __m512d packed_w2 = _mm512_set1_pd(w2);
                    const __m512d packed_w = _mm512_set1_pd(w);
                    const __m512d packed_iw = _mm512_set1_pd(iw);
#pragma GCC unroll 2
                    for (; last - i >= 8; i += 8) {
                        __m512d v0 = _mm512_loadu_pd(values + i);
                        __m512d v1 = _mm512_loadu_pd(values + i + quarter);
                        __m512d v2 = _mm512_loadu_pd(values + i + 2 * quarter);
                        __m512d v3 = _mm512_loadu_pd(values + i + 3 * quarter);
                        forward_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, packed_w2, packed_w, packed_iw);
                        _mm512_storeu_pd(values + i, v0);
                        _mm512_storeu_pd(values + i + quarter, v1);
                        _mm512_storeu_pd(values + i + 2 * quarter, v2);
                        _mm512_storeu_pd(values + i + 3 * quarter, v3);
                    }
                }
            }
#endif
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
                if (quarter >= 4) {
                    const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
                    const __m256d packed_w2 = _mm256_set1_pd(w2);
                    const __m256d packed_w = _mm256_set1_pd(w);
                    const __m256d packed_iw = _mm256_set1_pd(iw);
#pragma GCC unroll 2
                    for (; last - i >= 4; i += 4) {
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
            for (; i < last; ++i) {
                forward_butterfly<UseFma, UnitTwiddle>(arithmetic_, values[i], values[i + quarter],
                                                       values[i + 2 * quarter], values[i + 3 * quarter], w2, w, iw);
            }
        }

#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
        template<bool Inverse, bool BoundedInput>
        void four_layer_block(double *values, std::size_t length, std::size_t node) const {
            const std::size_t quarter = length / 4;
            const std::size_t sixteenth = length / 16;
            // Sixteen strips of 64 doubles occupy 8 KiB. Complete both radix-4 passes on
            // each strip group while it is hot, before recursing into its sixteen children.
            // This only changes visitation order; every butterfly retains the same inputs and bounds.
            constexpr std::size_t strip_size = 64;
            const auto pass = [&](double *base, std::size_t stride, std::size_t j, std::size_t first,
                                  std::size_t last) {
                if constexpr (Inverse) {
                    if (j == 0) {
                        inverse_radix4<true, true, BoundedInput>(base, stride, j, first, last);
                    } else {
                        inverse_radix4<true, false, BoundedInput>(base, stride, j, first, last);
                    }
                } else {
                    if (j == 0) {
                        forward_radix4<true, true>(base, stride, j, first, last);
                    } else {
                        forward_radix4<true, false>(base, stride, j, first, last);
                    }
                }
            };
            for (std::size_t first = 0; first < sixteenth; first += strip_size) {
                const std::size_t last = std::min(first + strip_size, sixteenth);
                if constexpr (Inverse) {
                    for (std::size_t child = 0; child < 4; ++child) {
                        pass(values + child * quarter, sixteenth, 4 * node + child, first, last);
                    }
                }
                for (std::size_t part = 0; part < 4; ++part) {
                    pass(values, quarter, node, first + part * sixteenth, last + part * sixteenth);
                }
                if constexpr (!Inverse) {
                    for (std::size_t child = 0; child < 4; ++child) {
                        pass(values + child * quarter, sixteenth, 4 * node + child, first, last);
                    }
                }
            }
        }
#endif

        template<bool UseFma>
        void forward_recursive(double *values, std::size_t length, std::size_t node) const {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
                if (length == 32 || length == 64) {
                    wide_basecase<false, false>(values, length, node);
                    return;
                }
#endif
                if (length >= 8 && length <= 32) {
                    if (node == 0) {
                        basecase<false, true>(values, length, node);
                    } else {
                        basecase<false, false>(values, length, node);
                    }
                    return;
                }
#if CRYPTO3_SMALL_PRIME_FFT_RADIX != 4
                constexpr std::size_t radix = CRYPTO3_SMALL_PRIME_FFT_RADIX;
                if (length >= 16 * radix) {
                    const auto stride = length / radix;
                    if (node == 0) {
                        large_radix<radix, false, true>(values, stride, node);
                    } else {
                        large_radix<radix, false, false>(values, stride, node);
                    }
                    for (std::size_t child = 0; child < radix; ++child) {
                        forward_recursive<true>(values + child * stride, stride, radix * node + child);
                    }
                    return;
                }
#endif
                if (length >= 4096) {
                    four_layer_block<false, false>(values, length, node);
                    for (std::size_t child = 0; child < 16; ++child) {
                        forward_recursive<true>(values + child * (length / 16), length / 16, 16 * node + child);
                    }
                    return;
                }
            }
#endif
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
                forward_radix4<UseFma, true>(values, quarter, node, 0, quarter);
            } else {
                forward_radix4<UseFma, false>(values, quarter, node, 0, quarter);
            }
            if (quarter > 1) {
                for (std::size_t child = 0; child < 4; ++child) {
                    forward_recursive<UseFma>(values + child * quarter, quarter, 4 * node + child);
                }
            }
        }

        template<bool UseFma, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void inverse_butterfly(const Arithmetic &arithmetic, Value &v0, Value &v1, Value &v2,
                                                        Value &v3, Value w2, Value w, Value iw, bool leaf) {
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
            Value left = b0 - b1;
            Value right = b2 - b3;
            if constexpr (UnitTwiddle) {
                left = reduce<UseFma>(arithmetic, left);
                right = multiply<UseFma>(arithmetic, right, iw);
            } else {
                multiply_pair<UseFma>(arithmetic, left, right, w, iw);
            }
            // Four inputs < 2*q give sums < 8*q < 2^53. Times a centered weight this is
            // < 4*q^2, so the a2 remainder is < 1.5*q. Left and right are each < q;
            // their sum is < 2*q and their weighted difference is < q. The invariant closes.
            v0 = reduce<UseFma>(arithmetic, even_sum + odd_sum);
            v1 = left + right;
            v2 = even_sum - odd_sum;
            v3 = left - right;
            if constexpr (UnitTwiddle) {
                v2 = reduce<UseFma>(arithmetic, v2);
                v3 = reduce<UseFma>(arithmetic, v3);
            } else {
                multiply_pair<UseFma>(arithmetic, v2, v3, w2, w2);
            }
        }

        template<bool UseFma, bool UnitTwiddle, bool BoundedInput>
        void inverse_radix4(double *values, std::size_t quarter, std::size_t node, std::size_t first,
                            std::size_t last) const {
            const double w2 = inverse_twiddle(node);
            const double w = inverse_twiddle(2 * node);
            const double iw = inverse_twiddle(2 * node + 1);
            std::size_t i = first;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
            if constexpr (UseFma) {
                if (quarter >= 8) {
                    const small_prime_fft_arithmetic::wide_fma_kernel packed(arithmetic_);
                    const __m512d packed_w2 = _mm512_set1_pd(w2);
                    const __m512d packed_w = _mm512_set1_pd(w);
                    const __m512d packed_iw = _mm512_set1_pd(iw);
#pragma GCC unroll 2
                    for (; last - i >= 8; i += 8) {
                        __m512d v0 = _mm512_loadu_pd(values + i);
                        __m512d v1 = _mm512_loadu_pd(values + i + quarter);
                        __m512d v2 = _mm512_loadu_pd(values + i + 2 * quarter);
                        __m512d v3 = _mm512_loadu_pd(values + i + 3 * quarter);
                        inverse_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, packed_w2, packed_w, packed_iw,
                                                             false);
                        _mm512_storeu_pd(values + i, v0);
                        _mm512_storeu_pd(values + i + quarter, v1);
                        _mm512_storeu_pd(values + i + 2 * quarter, v2);
                        _mm512_storeu_pd(values + i + 3 * quarter, v3);
                    }
                }
            }
#endif
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
                if (quarter >= 4) {
                    const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
                    const __m256d packed_w2 = _mm256_set1_pd(w2);
                    const __m256d packed_w = _mm256_set1_pd(w);
                    const __m256d packed_iw = _mm256_set1_pd(iw);
#pragma GCC unroll 2
                    for (; last - i >= 4; i += 4) {
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
            for (; i < last; ++i) {
                inverse_butterfly<UseFma, UnitTwiddle>(arithmetic_, values[i], values[i + quarter],
                                                       values[i + 2 * quarter], values[i + 3 * quarter], w2, w, iw,
                                                       quarter == 1 && !BoundedInput);
            }
        }

        template<bool UseFma, bool BoundedInput = false>
        void inverse_recursive(double *values, std::size_t length, std::size_t node) const {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
                if (length == 32 || length == 64) {
                    wide_basecase<true, BoundedInput>(values, length, node);
                    return;
                }
#endif
                if (length >= 8 && length <= 32) {
                    if (node == 0) {
                        basecase<true, true, BoundedInput>(values, length, node);
                    } else {
                        basecase<true, false, BoundedInput>(values, length, node);
                    }
                    return;
                }
#if CRYPTO3_SMALL_PRIME_FFT_RADIX != 4
                constexpr std::size_t radix = CRYPTO3_SMALL_PRIME_FFT_RADIX;
                if (length >= 16 * radix) {
                    const auto stride = length / radix;
                    for (std::size_t child = 0; child < radix; ++child) {
                        inverse_recursive<true, BoundedInput>(values + child * stride, stride, radix * node + child);
                    }
                    if (node == 0) {
                        large_radix<radix, true, true>(values, stride, node);
                    } else {
                        large_radix<radix, true, false>(values, stride, node);
                    }
                    return;
                }
#endif
                if (length >= 4096) {
                    for (std::size_t child = 0; child < 16; ++child) {
                        inverse_recursive<true, BoundedInput>(values + child * (length / 16), length / 16,
                                                              16 * node + child);
                    }
                    four_layer_block<true, BoundedInput>(values, length, node);
                    return;
                }
            }
#endif
            if (length == 1) {
                if constexpr (!BoundedInput) {
                    values[0] = reduce<UseFma>(arithmetic_, values[0]);
                }
                return;
            }
            if (length == 2) {
                const double left = BoundedInput ? values[0] : reduce<UseFma>(arithmetic_, values[0]);
                const double right = BoundedInput ? values[1] : reduce<UseFma>(arithmetic_, values[1]);
                values[0] = reduce<UseFma>(arithmetic_, left + right);
                values[1] = node == 0 ? reduce<UseFma>(arithmetic_, left - right) :
                                        multiply<UseFma>(arithmetic_, left - right, inverse_twiddle(node));
                return;
            }
            const std::size_t quarter = length / 4;
            if (quarter > 1) {
                for (std::size_t child = 0; child < 4; ++child) {
                    inverse_recursive<UseFma, BoundedInput>(values + child * quarter, quarter, 4 * node + child);
                }
            }
            if (node == 0) {
                inverse_radix4<UseFma, true, BoundedInput>(values, quarter, node, 0, quarter);
            } else {
                inverse_radix4<UseFma, false, BoundedInput>(values, quarter, node, 0, quarter);
            }
        }

        template<bool UseFma, bool BoundedInput = false>
        void inverse_impl(std::span<double> values, bool normalize) const {
            inverse_recursive<UseFma, BoundedInput>(values.data(), size_, 0);
            if (normalize) {
                const double q = static_cast<double>(modulus());
                std::size_t i = 0;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
                if constexpr (UseFma) {
                    const small_prime_fft_arithmetic::wide_fma_kernel packed(arithmetic_);
                    const __m512d scale = _mm512_set1_pd(size_inverse_);
                    const __m512d packed_q = _mm512_set1_pd(q);
                    for (; values.size() - i >= 8; i += 8) {
                        __m512d value = packed.multiply_fma(_mm512_loadu_pd(values.data() + i), scale);
                        const auto negative = _mm512_cmp_pd_mask(value, _mm512_setzero_pd(), _CMP_LT_OQ);
                        value = _mm512_mask_add_pd(value, negative, value, packed_q);
                        _mm512_storeu_pd(values.data() + i, value);
                    }
                }
#endif
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
     * identity with scalar, four-lane or eight-lane residue kernels. The existing field elements
     * bind their base arithmetic to the target modulus and cannot represent these small-prime intermediates.
     *
     * @pre Input cells are integer-valued doubles of magnitude < 2^53. scale is an integer with
     * |scale| < q/2, normally the FFT plan's centered 1/N. Buffers are disjoint or identical.
     * The thread's rounding mode stays unchanged during the call. Output cells are in (-q,q).
     * No buffers are allocated and no constant-time guarantee is provided.
     */
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
// The packed kernels use vector coordinates in the existing tower arrays.
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

            BOOST_FORCEINLINE void multiply_pair(Value &a, Value &b, Value wa, Value wb) const {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                if constexpr (UseFma && requires { arithmetic.multiply_fma_pair(a, b, wa, wb); }) {
                    arithmetic.multiply_fma_pair(a, b, wa, wb);
                } else
#endif
                {
                    a = multiply(a, wa);
                    b = multiply(b, wb);
                }
            }

            template<bool Reduce = true>
            Value add(Value a, Value b) const {
                if constexpr (Reduce) {
                    return reduce(a + b);
                } else {
                    return a + b;
                }
            }

            // Three reduced coordinates give an exact integer of magnitude < 3*q.
            // Callers can defer this reduction to the next bounded sum or final scaling.
            template<bool Reduce = true>
            Value subtract(Value a, Value b, Value c) const {
                if constexpr (Reduce) {
                    return reduce(a - b - c);
                } else {
                    return a - b - c;
                }
            }

            template<bool Reduce = true, typename T, std::size_t N>
            std::array<T, N> add(const std::array<T, N> &a, const std::array<T, N> &b) const {
                std::array<T, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                    result[i] = add<Reduce>(a[i], b[i]);
                }
                return result;
            }

            template<bool Reduce = true, typename T, std::size_t N>
            std::array<T, N> subtract(const std::array<T, N> &a, const std::array<T, N> &b,
                                      const std::array<T, N> &c) const {
                std::array<T, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                    result[i] = subtract<Reduce>(a[i], b[i], c[i]);
                }
                return result;
            }

            // Quadratic inputs and returned coordinates are in (-q,q). Pair sums are < 2*q;
            // the Karatsuba product is < 4*q^2 and its remainder < 1.5*q. Correcting by the
            // two diagonal products gives a sum < 3.5*q, hence exact in binary64 before reduction.
            BOOST_FORCEINLINE quadratic_type multiply(const quadratic_type &a, const quadratic_type &b) const {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
                // Keep the two independent products in registers and interleave their FMA chains.
                Value diagonal0 = a[0];
                Value diagonal1 = a[1];
                multiply_pair(diagonal0, diagonal1, b[0], b[1]);
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

            quadratic_type add_multiply_xi(const quadratic_type &a, const quadratic_type &b) const {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
                // a + (9+u)*b, with |a| < 3*q and |b| < q. After reducing 8*b,
                // each four-term sum is < 6*q < 2^53, so only its final reduction is needed.
                return {reduce(reduce(8 * b[0]) + b[0] - b[1] + a[0]), reduce(reduce(8 * b[1]) + b[1] + b[0] + a[1])};
            }

            cubic_type multiply_v(const cubic_type &a) const {
                return {multiply_xi(a[2]), a[0], a[1]};
            }

            cubic_type multiply(const cubic_type &a, const cubic_type &b) const {
                const auto diagonal0 = multiply(a[0], b[0]);
                const auto diagonal1 = multiply(a[1], b[1]);
                const auto diagonal2 = multiply(a[2], b[2]);
                // cross12 is reduced before multiplication by 8 inside multiply_xi.
                // The other cross terms stay < 3*q until these final sums; the cubic result is < q.
                // Finish each output pair before forming the next cross term to limit live temporaries.
                cubic_type result;
                {
                    const auto cross12 = subtract(multiply(add(a[1], a[2]), add(b[1], b[2])), diagonal1, diagonal2);
                    result[0] = add_multiply_xi(diagonal0, cross12);
                }
                {
                    const auto cross01 =
                        subtract<false>(multiply(add(a[0], a[1]), add(b[0], b[1])), diagonal0, diagonal1);
                    result[1] = add_multiply_xi(cross01, diagonal2);
                }
                {
                    const auto cross02 =
                        subtract<false>(multiply(add(a[0], a[2]), add(b[0], b[2])), diagonal0, diagonal2);
                    result[2] = add(cross02, diagonal1);
                }
                return result;
            }

            // Both tower operations return coordinates of magnitude < 3*q. The final
            // |scale| < q/2 gives a product < 1.5*q^2, whose FMA remainder is already in (-q,q).
            // Keeping these sums unreduced is also safe when two tower products are added below.
            tower_type multiply(const tower_type &a, const tower_type &b) const {
                const auto diagonal0 = multiply(a[0], b[0]);
                const auto diagonal1 = multiply(a[1], b[1]);
                const auto cross = subtract<false>(multiply(add(a[0], a[1]), add(b[0], b[1])), diagonal0, diagonal1);
                return {add<false>(diagonal0, multiply_v(diagonal1)), cross};
            }

            tower_type square(const tower_type &a) const {
                // The quadratic-extension square needs two cubic products instead of three:
                // c0 = (a0+a1)*(a0+v*a1) - a0*a1 - v*a0*a1; c1 = 2*a0*a1.
                const auto cross = multiply(a[0], a[1]);
                const auto combined = multiply(add(a[0], a[1]), add(a[0], multiply_v(a[1])));
                return {subtract<false>(combined, cross, multiply_v(cross)), add<false>(cross, cross)};
            }

            template<bool Square>
            BOOST_FORCEINLINE void load(tower_type &a, tower_type &b, std::span<const double> left,
                                        std::span<const double> right, std::size_t stride, std::size_t index) const {
                constexpr std::size_t width = sizeof(Value) / sizeof(double);
                for (std::size_t w = 0; w < 2; ++w) {
                    for (std::size_t v = 0; v < 3; ++v) {
                        for (std::size_t u = 0; u < 2; ++u) {
                            const std::size_t offset = (6 * w + 2 * v + u) * stride + index;
                            if constexpr (width == 1) {
                                a[w][v][u] = reduce(left[offset]);
                                if constexpr (!Square) {
                                    b[w][v][u] = reduce(right[offset]);
                                }
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                            } else {
                                a[w][v][u] = reduce(Arithmetic::load(left.data() + offset));
                                if constexpr (!Square) {
                                    b[w][v][u] = reduce(Arithmetic::load(right.data() + offset));
                                }
#endif
                            }
                        }
                    }
                }
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
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
                using wide_kernel = small_prime_fft_arithmetic::wide_fma_kernel;
                offset = stride - stride % wide_kernel::lanes;
                if (offset != 0) {
                    const wide_kernel wide(arithmetic_);
                    apply_range<true, Square, Sum>(wide, output, left, right, extra_left, extra_right, 0, offset,
                                                   wide_kernel::broadcast(scale));
                }
#endif
                using packed_kernel = small_prime_fft_arithmetic::packed_fma_kernel;
                const std::size_t packed_end = stride - stride % packed_kernel::lanes;
                if (offset < packed_end) {
                    const packed_kernel packed(arithmetic_);
                    apply_range<true, Square, Sum>(packed, output, left, right, extra_left, extra_right, offset,
                                                   packed_end, packed_kernel::broadcast(scale));
                    offset = packed_end;
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
                tower_type a, b;
                operations.template load<Square>(a, b, left, right, stride, i);
                tower_type result;
                if constexpr (Square) {
                    result = operations.square(a);
                } else {
                    result = operations.multiply(a, b);
                    if constexpr (Sum) {
                        // Reuse the input slots after the first product. Load both remaining
                        // inputs before storing anything so output may still alias any input.
                        operations.template load<false>(a, b, extra_left, extra_right, stride, i);
                        // Each product coordinate has magnitude < 3*q; their sum is exact below 6*q.
                        // Reduce that sum before the usual inverse-transform normalization.
                        result = operations.add(result, operations.multiply(a, b));
                    }
                }
                for (std::size_t w = 0; w < 2; ++w) {
                    for (std::size_t v = 0; v < 3; ++v) {
                        // Fold inverse-transform normalization into this existing write pass,
                        // pairing the independent scales of the two quadratic coordinates.
                        const std::size_t offset = (6 * w + 2 * v) * stride + i;
                        Value first = result[w][v][0];
                        Value second = result[w][v][1];
                        operations.multiply_pair(first, second, scale, scale);
                        if constexpr (width == 1) {
                            output[offset] = first;
                            output[offset + stride] = second;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                        } else {
                            Arithmetic::store(output.data() + offset, first);
                            Arithmetic::store(output.data() + offset + stride, second);
#endif
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
