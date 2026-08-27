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

#ifndef CRYPTO3_MATH_POLYNOMIAL_RATIONAL_RECONSTRUCTION_HPP
#define CRYPTO3_MATH_POLYNOMIAL_RATIONAL_RECONSTRUCTION_HPP

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>

#include <nil/crypto3/math/polynomial/arithmetic/half_gcd.hpp>

namespace nil::crypto3::math {

    /**
     * Reconstruct polynomials numerator P and denominator Q such that
     *
     *     P = residue * Q mod modulus,
     *
     * with degree(P) at most maximum_numerator_degree and degree(Q) at most
     * maximum_denominator_degree. The denominator is normalized to monic form, and the numerator is scaled by the
     * same field element. Outputs are replaced only when reconstruction succeeds and may alias either input, but must
     * be distinct from each other.
     *
     * Starting with consecutive Euclidean remainders R0 = modulus and R1 = residue, write
     *
     *     Ri = Ui * modulus + Qi * residue.
     *
     * Here Ui is the coefficient of modulus and Qi is the coefficient of residue. Initially U0 = 1, Q0 = 0 and
     * U1 = 0, Q1 = 1. The Ui terms vanish modulo modulus, so the algorithm stores only Q0 and Q1.
     *
     * A Euclidean quotient q replaces the pairs by
     *
     *     (R0, R1) <- (R1, R0 - q * R1),
     *     (Q0, Q1) <- (Q1, Q0 - q * Q1).
     *
     * Repeat this step while R1 is nonzero and degree(R1) exceeds maximum_numerator_degree. The loop normally stops
     * before the Euclidean sequence reaches zero. At that point R1 is the first remainder within the numerator bound,
     * so R1 is a valid numerator and Q1 is its denominator. Reconstruction fails if Q1 exceeds the denominator bound.
     *
     * The strict condition
     *
     *     maximum_numerator_degree + maximum_denominator_degree < degree(modulus)
     *
     * ensures that at most one rational function can satisfy the requested bounds.
     *
     * @return true on success; false if the Euclidean candidate does not satisfy the denominator bound.
     * @throws std::invalid_argument if the outputs are the same object; an input is empty or noncanonical; modulus is
     *         constant; residue is not reduced modulo modulus; or the degree bounds do not satisfy the strict
     *         uniqueness condition.
     */
    template<detail::SupportsDivrem Backend>
    bool rational_reconstruct(typename Backend::polynomial_type &numerator,
                              typename Backend::polynomial_type &denominator,
                              const typename Backend::polynomial_type &residue,
                              const typename Backend::polynomial_type &modulus,
                              std::size_t maximum_numerator_degree,
                              std::size_t maximum_denominator_degree,
                              polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        if (std::addressof(numerator) == std::addressof(denominator)) {
            throw std::invalid_argument("rational-reconstruction outputs must be distinct objects");
        }
        const auto is_canonical = [](const polynomial_type &polynomial) {
            return !polynomial.empty() &&
                   (polynomial.size() == 1 || polynomial[polynomial.size() - 1] != value_type {});
        };
        if (!is_canonical(residue) || !is_canonical(modulus)) {
            throw std::invalid_argument("rational reconstruction requires canonical nonempty inputs");
        }
        if (modulus.size() == 1) {
            throw std::invalid_argument("rational reconstruction requires a nonconstant modulus");
        }

        const std::size_t modulus_degree = modulus.size() - 1;
        if (residue.size() > modulus_degree) {
            throw std::invalid_argument("rational reconstruction requires a reduced residue");
        }
        // Express the strict sum bound without a potentially overflowing addition.
        if (maximum_numerator_degree >= modulus_degree ||
            maximum_denominator_degree >= modulus_degree - maximum_numerator_degree) {
            throw std::invalid_argument("rational-reconstruction degree bounds do not ensure uniqueness");
        }

        polynomial_type previous_remainder(modulus);
        polynomial_type current_remainder(residue);
        polynomial_type previous_denominator = {value_type {}};
        polynomial_type current_denominator = {value_type::one()};
        polynomial_type quotient;
        polynomial_type next_remainder;
        polynomial_type quotient_times_denominator;
        polynomial_type next_denominator;

        while (!is_zero(current_remainder) && current_remainder.size() - 1 > maximum_numerator_degree) {
            detail::gcd_divrem_step(quotient, next_remainder, previous_remainder, current_remainder,
                                    arithmetic_context);
            arithmetic_context.multiply(quotient_times_denominator, quotient, current_denominator);
            subtraction(next_denominator, previous_denominator, quotient_times_denominator);

            previous_remainder = std::move(current_remainder);
            current_remainder = std::move(next_remainder);
            previous_denominator = std::move(current_denominator);
            current_denominator = std::move(next_denominator);
        }

        if (is_zero(current_denominator) || current_denominator.size() - 1 > maximum_denominator_degree) {
            return false;
        }

        const value_type denominator_leading_coefficient = current_denominator[current_denominator.size() - 1];
        if (denominator_leading_coefficient != value_type::one()) {
            const value_type normalization = denominator_leading_coefficient.inversed();
            scalar_multiplication(current_remainder, current_remainder, normalization);
            scalar_multiplication(current_denominator, current_denominator, normalization);
        }

        numerator = std::move(current_remainder);
        denominator = std::move(current_denominator);
        return true;
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_RATIONAL_RECONSTRUCTION_HPP
