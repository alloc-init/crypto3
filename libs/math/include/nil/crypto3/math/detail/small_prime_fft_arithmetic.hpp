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
#include <cstdint>
#include <limits>
#include <stdexcept>

#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>

namespace nil::crypto3::math::detail {

    /**
     * Scalar residue arithmetic for an odd 50-bit modulus q. Values are integer-valued binary64
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
