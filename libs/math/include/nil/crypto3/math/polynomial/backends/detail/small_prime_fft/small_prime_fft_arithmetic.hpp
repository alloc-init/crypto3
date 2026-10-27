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

#ifndef CRYPTO3_MATH_DETAIL_SMALL_PRIME_FFT_ARITHMETIC_HPP
#define CRYPTO3_MATH_DETAIL_SMALL_PRIME_FFT_ARITHMETIC_HPP

#include <cassert>
#include <cfenv>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>

#if defined(__x86_64__) && defined(__AVX2__) && defined(__FMA__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED
#if defined(__AVX512F__) && defined(__AVX512DQ__)
#define CRYPTO3_SMALL_PRIME_FFT_AVX512_AVAILABLE
#if !defined(CRYPTO3_SMALL_PRIME_FFT_DISABLE_AVX512)
#define CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED
#endif
#endif
#endif

namespace nil::crypto3::math::detail {

    /**
     * Residue arithmetic for an odd 50-bit modulus q. Values are integer-valued binary64
     * numbers, not approximate coefficients. Primality is a caller concern (the prime catalogue
     * is separately certified); these arithmetic operations also work for odd composite moduli.
     *
     * Select the kernel once per operation with fma_supported(). The FMA methods require this
     * predicate to remain true until the operation finishes, including FE_TONEAREST on the calling
     * thread. The integer methods are the exact fallback for other rounding modes/compiler flags.
     * Keeping this dispatch outside each scalar operation avoids a rounding-mode query per butterfly.
     * No method changes the floating-point environment.
     */
    class small_prime_fft_arithmetic {
        using wide_type = boost::multiprecision::uint128_modular_t;

        static_assert(std::numeric_limits<double>::radix == 2 && std::numeric_limits<double>::digits == 53 &&
                          std::numeric_limits<double>::max_exponent == 1024,
                      "small-prime residues require binary64 doubles");

    public:
        explicit small_prime_fft_arithmetic(std::uint64_t modulus) : modulus_(modulus) {
            if (modulus < (std::uint64_t(1) << 49) || modulus >= (std::uint64_t(1) << 50) || modulus % 2 == 0) {
                throw std::invalid_argument("small_prime_fft_arithmetic: expected an odd 50-bit modulus");
            }

            // Since 2^49 < q < 2^50, the nearest binary64 reciprocal is S/2^102, where
            // S = round(2^102/q) fits exactly in 53 bits. An odd q cannot produce a halfway tie.
            // Construct S using integers so that setup is independent of the active rounding mode.
            const wide_type scale = wide_type(1) << 102;
            const wide_type significand = (scale + modulus / 2) / modulus;
            reciprocal_ = std::ldexp(static_cast<double>(significand.convert_to<std::uint64_t>()), -102);

            // For h = RN(a*b), l = fma(a,b,-h), and k = nearbyint(h*reciprocal), let
            // epsilon = reciprocal - 1/q and |a*b| < 2^e. Then h+l is exact and
            // |a*b/q - k| <= |a*b*epsilon| + 2^(e-54)*reciprocal + rounding_error.
            // For |a*b| < 2*q^2 the rounding_error is at most 5/8. With
            // D = |q*S - 2^102|, the remaining numerator over 2^102 is
            // T = 2*q*D + 2^(e-54)*S. Requiring T < 3*2^99 makes the total < 1.
            // For |a*b| < 4*q^2, e increases by one and the rounding_error is 3/4:
            // the condition becomes 2*T < 3*2^100, giving an output strictly inside (-1.5*q,1.5*q).
            // All certificate intermediates fit within 128 bits.
            const wide_type scaled_reciprocal = significand * modulus;
            const wide_type error = scaled_reciprocal > scale ? scaled_reciprocal - scale : scale - scaled_reciprocal;
            const wide_type product_limit = 2 * wide_type(modulus) * modulus;
            const unsigned exponent = boost::multiprecision::msb(product_limit) + 1;
            const wide_type numerator = 2 * wide_type(modulus) * error + (significand << (exponent - 54));
            reciprocal_satisfies_bounds_ = numerator < (wide_type(3) << 99);
        }

        std::uint64_t modulus() const {
            return modulus_;
        }

        double reciprocal() const {
            return reciprocal_;
        }

        /** Check once before calling FMA kernels; use the integer kernels when this returns false. */
        bool fma_supported() const {
#if defined(__FAST_MATH__) || defined(__ASSOCIATIVE_MATH__) || defined(__RECIPROCAL_MATH__) || defined(_M_FP_FAST) || \
    (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0) || FLT_EVAL_METHOD != 0
            return false;
#else
            return std::numeric_limits<double>::is_iec559 && reciprocal_satisfies_bounds_ &&
                   std::fegetround() == FE_TONEAREST;
#endif
        }

        /**
         * Return an integer representative of value modulo q strictly inside (-q,q).
         * @pre value is an integer with |value| < 2^53, and fma_supported() is true.
         */
        double reduce_fma(double value) const {
#if defined(__clang__)
            // Clang does not expose every unsafe floating-point option through predefined macros.
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
            assert(valid_integer(value));
            assert(fma_supported());
            // |value/q| < 16, so quotient rounding error is < 1. The FMA remainder is an
            // exactly representable integer of magnitude < q, even for a large input representative.
            const double quotient = std::nearbyint(value * reciprocal_);
            return std::fma(-quotient, static_cast<double>(modulus_), value);
        }

        /**
         * Return an integer representative of left*right modulo q. The output is strictly inside
         * (-q,q) when |left*right| < 2*q^2, and inside (-1.5*q,1.5*q) when |left*right| < 4*q^2.
         * @pre both operands are integers of magnitude < 2^53, |left*right| < 4*q^2, and
         *      fma_supported() is true. These are integer product bounds, not rounded double tests.
         */
        double multiply_fma(double left, double right) const {
#if defined(__clang__)
            // Preserve the rounded high product even when the caller enables reassociation/contraction.
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
            assert(valid_integer(left) && valid_integer(right));
            assert(product_magnitude(left, right) < 4 * wide_type(modulus_) * modulus_);
            assert(fma_supported());
            const double high = left * right;
            const double low = std::fma(left, right, -high);
            const double quotient = std::nearbyint(high * reciprocal_);
            // Both the compensated subtraction and the final remainder are small enough to be
            // exact integers in binary64. Preserve the split product: reassociation is not allowed.
            return std::fma(-quotient, static_cast<double>(modulus_), high) + low;
        }

#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
        /**
         * Snapshot the two packed constants once before a batch of FMA operations. Keeping them
         * here avoids repeating integer conversion and broadcasts at each arithmetic operation.
         * The source arithmetic must remain unchanged and outlive this kernel. The scalar FMA
         * preconditions still apply.
         *
         * Under the required fixed FE_TONEAREST mode, the assembly computes register outputs
         * from its inputs. It is non-volatile so the compiler can reuse the packed constants
         * across operations; the instructions within each block still stay together.
         */
        class packed_fma_kernel {
        public:
            using value_type = __m256d;
            static constexpr std::size_t lanes = 4;

            static value_type broadcast(double value) {
                return _mm256_set1_pd(value);
            }

            static value_type load(const double *values) {
                return _mm256_loadu_pd(values);
            }

            static void store(double *values, value_type value) {
                _mm256_storeu_pd(values, value);
            }

            explicit packed_fma_kernel(const small_prime_fft_arithmetic &arithmetic) :
                arithmetic_(arithmetic), modulus_(_mm256_set1_pd(static_cast<double>(arithmetic.modulus_))),
                reciprocal_(_mm256_set1_pd(arithmetic.reciprocal_)) {
                assert(arithmetic.fma_supported());
            }

            /** Four independent reductions with the same preconditions and bounds as the scalar overload. */
            __m256d reduce_fma(__m256d value) const {
                assert(arithmetic_.fma_supported());
#ifndef NDEBUG
                alignas(32) double lanes[4];
                _mm256_storeu_pd(lanes, value);
                for (double lane : lanes) {
                    assert(valid_integer(lane));
                }
#endif
                __m256d quotient;
                // Round to nearest, ties to even, suppressing the inexact exception just like nearbyint.
                // All operands are registers; the early-clobber scratch cannot overlap a live input.
                // clang-format off
                asm(
                    "vmulpd %[reciprocal], %[value], %[quotient]\n"
                    "vroundpd $8, %[quotient], %[quotient]\n"
                    "vfnmadd231pd %[modulus], %[quotient], %[value]\n"
                    : [value] "+x"(value), [quotient] "=&x"(quotient)
                    : [modulus] "x"(modulus_), [reciprocal] "x"(reciprocal_)
                );
                // clang-format on
                return value;
            }

            /** Four independent compensated products; the scalar product bounds apply to every lane. */
            __m256d multiply_fma(__m256d left, __m256d right) const {
                assert(arithmetic_.fma_supported());
#ifndef NDEBUG
                alignas(32) double a[4], b[4];
                _mm256_storeu_pd(a, left);
                _mm256_storeu_pd(b, right);
                for (std::size_t i = 0; i < 4; ++i) {
                    assert(valid_integer(a[i]) && valid_integer(b[i]));
                    assert(product_magnitude(a[i], b[i]) < 4 * wide_type(arithmetic_.modulus_) * arithmetic_.modulus_);
                }
#endif
                __m256d high, low, quotient;
                // Keep the rounded high product and its exact FMA residual separate. The quotient
                // estimate and low product are independent once high is available.
                // clang-format off
                asm(
                    "vmulpd %[right], %[left], %[high]\n"
                    "vmulpd %[reciprocal], %[high], %[quotient]\n"
                    "vmovapd %[high], %[low]\n"
                    "vfmsub231pd %[right], %[left], %[low]\n"
                    "vroundpd $8, %[quotient], %[quotient]\n"
                    "vfnmadd231pd %[modulus], %[quotient], %[high]\n"
                    "vaddpd %[low], %[high], %[high]\n"
                    : [high] "=&x"(high), [low] "=&x"(low), [quotient] "=&x"(quotient)
                    : [left] "x"(left), [right] "x"(right),
                      [modulus] "x"(modulus_), [reciprocal] "x"(reciprocal_)
                );
                // clang-format on
                return high;
            }

            /**
             * Two independent compensated products. The distinct left vectors are overwritten;
             * right operands are copied before either write. Interleave the independent product
             * and quotient work while sharing one quotient scratch register. Bounds match multiply_fma.
             */
            void multiply_fma_pair(__m256d &left0, __m256d &left1, __m256d right0, __m256d right1) const {
                assert(&left0 != &left1);
                assert(arithmetic_.fma_supported());
#ifndef NDEBUG
                alignas(32) double a[8], b[8];
                _mm256_storeu_pd(a, left0);
                _mm256_storeu_pd(a + 4, left1);
                _mm256_storeu_pd(b, right0);
                _mm256_storeu_pd(b + 4, right1);
                for (std::size_t i = 0; i < 8; ++i) {
                    assert(valid_integer(a[i]) && valid_integer(b[i]));
                    assert(product_magnitude(a[i], b[i]) < 4 * wide_type(arithmetic_.modulus_) * arithmetic_.modulus_);
                }
#endif
                __m256d high0, high1, quotient;
                // The early-clobber left operands hold their exact low products after vfmsub213pd.
                // Keep right operands distinct from both outputs until both residuals are computed.
                // clang-format off
                asm(
                    "vmulpd %[right0], %[left0], %[high0]\n"
                    "vmulpd %[right1], %[left1], %[high1]\n"
                    "vmulpd %[reciprocal], %[high0], %[quotient]\n"
                    "vfmsub213pd %[high0], %[right0], %[left0]\n"
                    "vfmsub213pd %[high1], %[right1], %[left1]\n"
                    "vroundpd $8, %[quotient], %[quotient]\n"
                    "vfnmadd231pd %[modulus], %[quotient], %[high0]\n"
                    "vmulpd %[reciprocal], %[high1], %[quotient]\n"
                    "vroundpd $8, %[quotient], %[quotient]\n"
                    "vfnmadd231pd %[modulus], %[quotient], %[high1]\n"
                    "vaddpd %[high0], %[left0], %[left0]\n"
                    "vaddpd %[high1], %[left1], %[left1]\n"
                    : [left0] "+&x"(left0), [left1] "+&x"(left1),
                      [high0] "=&x"(high0), [high1] "=&x"(high1), [quotient] "=&x"(quotient)
                    : [right0] "x"(right0), [right1] "x"(right1),
                      [modulus] "x"(modulus_), [reciprocal] "x"(reciprocal_)
                );
                // clang-format on
            }

        private:
            const small_prime_fft_arithmetic &arithmetic_;
            const __m256d modulus_;
            const __m256d reciprocal_;
        };

#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_AVAILABLE)
        // Eight independent lanes use the same reciprocal certificate and exact product bounds.
        // AVX-512 register constraints permit all 32 ZMM registers; no runtime CPU dispatch is
        // needed because the compiler target already requires AVX-512F/DQ.
        class wide_fma_kernel {
        public:
            using value_type = __m512d;
            static constexpr std::size_t lanes = 8;

            static value_type broadcast(double value) {
                return _mm512_set1_pd(value);
            }

            static value_type load(const double *values) {
                return _mm512_loadu_pd(values);
            }

            static void store(double *values, value_type value) {
                _mm512_storeu_pd(values, value);
            }

            explicit wide_fma_kernel(const small_prime_fft_arithmetic &arithmetic) :
                arithmetic_(arithmetic), modulus_(_mm512_set1_pd(static_cast<double>(arithmetic.modulus_))),
                reciprocal_(_mm512_set1_pd(arithmetic.reciprocal_)) {
                assert(arithmetic.fma_supported());
            }

            /** Eight independent reductions with the same preconditions and bounds as the scalar overload. */
            __m512d reduce_fma(__m512d value) const {
                assert(arithmetic_.fma_supported());
#ifndef NDEBUG
                // Keep debug scratch SIMD-aligned even when sanitizers change the stack layout.
                alignas(64) double lanes[8];
                _mm512_storeu_pd(lanes, value);
                for (double lane : lanes) {
                    assert(valid_integer(lane));
                }
#endif
                __m512d quotient;
                // Round to nearest, ties to even, suppressing the inexact exception just like nearbyint.
                // All operands are registers; the early-clobber scratch cannot overlap a live input.
                // clang-format off
                asm(
                    "vmulpd %[reciprocal], %[value], %[quotient]\n"
                    "vrndscalepd $8, %[quotient], %[quotient]\n"
                    "vfnmadd231pd %[modulus], %[quotient], %[value]\n"
                    : [value] "+v"(value), [quotient] "=&v"(quotient)
                    : [modulus] "v"(modulus_), [reciprocal] "v"(reciprocal_)
                );
                // clang-format on
                return value;
            }

            /** Eight independent compensated products; the scalar product bounds apply to every lane. */
            __m512d multiply_fma(__m512d left, __m512d right) const {
                assert(arithmetic_.fma_supported());
#ifndef NDEBUG
                alignas(64) double a[8], b[8];
                _mm512_storeu_pd(a, left);
                _mm512_storeu_pd(b, right);
                for (std::size_t i = 0; i < 8; ++i) {
                    assert(valid_integer(a[i]) && valid_integer(b[i]));
                    assert(product_magnitude(a[i], b[i]) < 4 * wide_type(arithmetic_.modulus_) * arithmetic_.modulus_);
                }
#endif
                __m512d high, low, quotient;
                // Keep the rounded high product and its exact FMA residual separate. The quotient
                // estimate and low product are independent once high is available.
                // clang-format off
                asm(
                    "vmulpd %[right], %[left], %[high]\n"
                    "vmulpd %[reciprocal], %[high], %[quotient]\n"
                    "vmovapd %[high], %[low]\n"
                    "vfmsub231pd %[right], %[left], %[low]\n"
                    "vrndscalepd $8, %[quotient], %[quotient]\n"
                    "vfnmadd231pd %[modulus], %[quotient], %[high]\n"
                    "vaddpd %[low], %[high], %[high]\n"
                    : [high] "=&v"(high), [low] "=&v"(low), [quotient] "=&v"(quotient)
                    : [left] "v"(left), [right] "v"(right),
                      [modulus] "v"(modulus_), [reciprocal] "v"(reciprocal_)
                );
                // clang-format on
                return high;
            }

            /**
             * Two independent compensated products. The distinct left vectors are overwritten;
             * right operands are copied before either write. Interleave the independent product
             * and quotient work with separate scratch registers. Bounds match multiply_fma.
             */
            void multiply_fma_pair(__m512d &left0, __m512d &left1, __m512d right0, __m512d right1) const {
                assert(&left0 != &left1);
                assert(arithmetic_.fma_supported());
#ifndef NDEBUG
                alignas(64) double a[16], b[16];
                _mm512_storeu_pd(a, left0);
                _mm512_storeu_pd(a + 8, left1);
                _mm512_storeu_pd(b, right0);
                _mm512_storeu_pd(b + 8, right1);
                for (std::size_t i = 0; i < 16; ++i) {
                    assert(valid_integer(a[i]) && valid_integer(b[i]));
                    assert(product_magnitude(a[i], b[i]) < 4 * wide_type(arithmetic_.modulus_) * arithmetic_.modulus_);
                }
#endif
                __m512d high0, high1, quotient0, quotient1;
                // The early-clobber left operands hold their exact low products after vfmsub213pd.
                // Keep right operands distinct from both outputs until both residuals are computed.
                // clang-format off
                asm(
                    "vmulpd %[right0], %[left0], %[high0]\n"
                    "vmulpd %[right1], %[left1], %[high1]\n"
                    "vmulpd %[reciprocal], %[high0], %[quotient0]\n"
                    "vmulpd %[reciprocal], %[high1], %[quotient1]\n"
                    "vfmsub213pd %[high0], %[right0], %[left0]\n"
                    "vfmsub213pd %[high1], %[right1], %[left1]\n"
                    "vrndscalepd $8, %[quotient0], %[quotient0]\n"
                    "vrndscalepd $8, %[quotient1], %[quotient1]\n"
                    "vfnmadd231pd %[modulus], %[quotient0], %[high0]\n"
                    "vfnmadd231pd %[modulus], %[quotient1], %[high1]\n"
                    "vaddpd %[high0], %[left0], %[left0]\n"
                    "vaddpd %[high1], %[left1], %[left1]\n"
                    : [left0] "+&v"(left0), [left1] "+&v"(left1),
                      [high0] "=&v"(high0), [high1] "=&v"(high1),
                      [quotient0] "=&v"(quotient0), [quotient1] "=&v"(quotient1)
                    : [right0] "v"(right0), [right1] "v"(right1),
                      [modulus] "v"(modulus_), [reciprocal] "v"(reciprocal_)
                );
                // clang-format on
            }

        private:
            const small_prime_fft_arithmetic &arithmetic_;
            const __m512d modulus_;
            const __m512d reciprocal_;
        };
#endif
#endif

        /** Exact fallback for reduce_fma, with the same input/output bounds and any rounding mode. */
        double reduce_integer(double value) const {
            assert(valid_integer(value));
            return static_cast<double>(static_cast<std::int64_t>(value) % static_cast<std::int64_t>(modulus_));
        }

        /**
         * Exact fallback: both operands must be integers of magnitude < 2^53. Their full product
         * fits in 106 bits, so the existing 128-bit integer type suffices. The result is inside (-q,q).
         */
        double multiply_integer(double left, double right) const {
            assert(valid_integer(left) && valid_integer(right));
            const wide_type residue = product_magnitude(left, right) % modulus_;
            const double result = static_cast<double>(residue.convert_to<std::uint64_t>());
            return (left < 0) != (right < 0) ? -result : result;
        }

    private:
        static bool valid_integer(double value) {
            return std::isfinite(value) && std::abs(value) < 0x1p53 && std::trunc(value) == value;
        }

        // Call only after checking the integer/range preconditions, so both conversions are exact.
        static wide_type product_magnitude(double left, double right) {
            return wide_type(static_cast<std::uint64_t>(std::abs(left))) * static_cast<std::uint64_t>(std::abs(right));
        }

        std::uint64_t modulus_;
        double reciprocal_;
        bool reciprocal_satisfies_bounds_;
    };

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_DETAIL_SMALL_PRIME_FFT_ARITHMETIC_HPP
