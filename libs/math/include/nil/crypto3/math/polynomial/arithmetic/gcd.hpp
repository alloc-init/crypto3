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

#ifndef CRYPTO3_MATH_POLYNOMIAL_GCD_HPP
#define CRYPTO3_MATH_POLYNOMIAL_GCD_HPP

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <utility>

#include <nil/crypto3/math/polynomial/arithmetic/half_gcd.hpp>

namespace nil::crypto3::math {

    /**
     * Compute the monic greatest common divisor of two canonical coefficient polynomials. Small inputs use the
     * Euclidean algorithm. When enabled by the configured cutoff, larger inputs use recursive half-GCD reduction.
     * Polynomial products and division steps use the configured arithmetic backend. Output may alias either input.
     * The GCD of two zero polynomials is the zero polynomial.
     *
     * The Euclidean algorithm divides the larger polynomial A by the smaller polynomial B and replaces the pair
     * (A, B) with (B, A mod B). The remainder has lower degree than B, so repeating this step eventually leaves a zero
     * remainder. This replacement preserves all common divisors: a divisor of A and B also divides A - q * B, the
     * remainder, while a divisor of B and the remainder also divides q * B + remainder, which is A. At termination,
     * the common divisors are therefore exactly the divisors of the last nonzero remainder, making its monic form the
     * GCD. Half-GCD batches several of these same quotient steps into one recursively constructed polynomial-matrix
     * transformation.
     */
    template<detail::SupportsDivrem Backend>
        requires std::copy_constructible<typename Backend::polynomial_type> &&
                 std::movable<typename Backend::polynomial_type> &&
                 requires(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input) {
                     make_monic(output, input);
                 }
    void gcd(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &left,
             const typename Backend::polynomial_type &right,
             polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;

        polynomial_type previous_remainder(left);
        polynomial_type current_remainder(right);
        if (is_zero(previous_remainder)) {
            if (is_zero(current_remainder)) {
                output = std::move(current_remainder);
            } else {
                make_monic(output, current_remainder);
            }
            return;
        }
        if (is_zero(current_remainder)) {
            make_monic(output, previous_remainder);
            return;
        }

        if (previous_remainder.size() < current_remainder.size()) {
            std::swap(previous_remainder, current_remainder);
        }

        polynomial_type quotient;
        polynomial_type next_remainder;
        detail::gcd_divrem_step(quotient, next_remainder, previous_remainder, current_remainder, arithmetic_context);
        previous_remainder = std::move(current_remainder);
        current_remainder = std::move(next_remainder);

        const auto &options = arithmetic_context.options();
        while (!is_zero(current_remainder)) {
            if (options.gcd_half_gcd_cutoff == 0 || current_remainder.size() < options.gcd_half_gcd_cutoff) {
                detail::gcd_divrem_step(quotient, next_remainder, previous_remainder, current_remainder,
                                        arithmetic_context);
                previous_remainder = std::move(current_remainder);
                current_remainder = std::move(next_remainder);
                continue;
            }

            polynomial_type reduced_first;
            polynomial_type reduced_second;
            detail::half_gcd_reduce(reduced_first, reduced_second, previous_remainder, current_remainder,
                                    options.half_gcd_basecase_cutoff, arithmetic_context);
            previous_remainder = std::move(reduced_first);
            current_remainder = std::move(reduced_second);
            if (is_zero(current_remainder)) {
                break;
            }

            detail::gcd_divrem_step(quotient, next_remainder, previous_remainder, current_remainder,
                                    arithmetic_context);
            previous_remainder = std::move(current_remainder);
            current_remainder = std::move(next_remainder);
        }

        make_monic(output, previous_remainder);
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_GCD_HPP
