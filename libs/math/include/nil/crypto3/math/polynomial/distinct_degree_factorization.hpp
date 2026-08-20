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

#ifndef CRYPTO3_MATH_DISTINCT_DEGREE_FACTORIZATION_HPP
#define CRYPTO3_MATH_DISTINCT_DEGREE_FACTORIZATION_HPP

#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <utility>

#include <nil/crypto3/math/polynomial/gcd.hpp>
#include <nil/crypto3/math/polynomial/polynomial_factorization.hpp>
#include <nil/crypto3/math/polynomial/polynomial_frobenius.hpp>

namespace nil::crypto3::math {

    /**
     * Split a square-free polynomial into products of irreducible factors of equal degree using the classical
     * distinct-degree algorithm. This implementation is intended as a correctness reference for faster blocked
     * algorithms.
     *
     * Over a finite field with Q elements, X^(Q^d) - X is the product of all monic irreducible polynomials whose
     * degrees divide d. The algorithm visits d = 1, 2, ... and removes the factors found at earlier degrees. Therefore
     *
     *     gcd(remaining, X^(Q^d) - X)
     *
     * contains exactly the remaining irreducible factors of degree d. Frobenius powers are computed modulo the
     * original monic input. This is also valid modulo every factor subsequently removed from that input and avoids
     * rebuilding the Frobenius precomputation after each split.
     *
     * The input is normalized to monic form and its original leading coefficient is preserved in the result. Zero
     * and constant inputs return that coefficient and no factors. Nonconstant input must be square-free; this is
     * checked before the decomposition begins.
     *
     * Factors are emitted in increasing irreducible-factor degree. After each factor is appended to the result,
     * factor_callback may request an early stop. A stopped result includes that factor and has complete set to false.
     *
     * @throws std::invalid_argument if a nonconstant input is not square-free.
     * @pre input is a nonempty coefficient polynomial.
     */
    template<detail::SupportsDivrem Backend, typename FactorCallback>
        requires requires(FactorCallback &callback,
                          const distinct_degree_factor<typename Backend::polynomial_type> &factor) {
            { callback(factor) } -> std::same_as<factorization_control>;
        }
    distinct_degree_factorization_result<typename Backend::polynomial_type>
        distinct_degree_factorization_reference(const typename Backend::polynomial_type &input,
                                                polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                                FactorCallback &&factor_callback) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;
        using result_type = distinct_degree_factorization_result<polynomial_type>;

        result_type result;
        polynomial_type monic_input(input);
        condense(monic_input);
        result.leading_coefficient = monic_input.back();
        if (monic_input.size() == 1) {
            return result;
        }
        make_monic(monic_input, monic_input);

        polynomial_type input_derivative;
        derivative(input_derivative, monic_input);
        polynomial_type repeated_factor;
        gcd(repeated_factor, monic_input, input_derivative, arithmetic_context);
        if (repeated_factor.size() > 1) {
            throw std::invalid_argument("distinct-degree factorization requires a square-free polynomial");
        }

        const polynomial_type x = {value_type {}, value_type::one()};
        polynomial_type frobenius_power(x);
        polynomial_frobenius_context<Backend> frobenius_context(monic_input, arithmetic_context);
        polynomial_type remaining(std::move(monic_input));

        std::size_t irreducible_factor_degree = 1;
        // Every factor of degree below the current degree has already been removed. If remaining has degree less
        // than twice the current degree, it can contain at most one irreducible factor and is emitted after the loop.
        while (irreducible_factor_degree <= (remaining.size() - 1) / 2) {
            frobenius_map(frobenius_power, frobenius_power, frobenius_context, arithmetic_context);

            polynomial_type frobenius_difference;
            subtraction(frobenius_difference, frobenius_power, x);
            polynomial_type factor;
            gcd(factor, remaining, frobenius_difference, arithmetic_context);
            if (factor.size() > 1) {
                result.factors.push_back({std::move(factor), irreducible_factor_degree});
                if (factor_callback(result.factors.back()) == factorization_control::stop_factorization) {
                    result.complete = false;
                    return result;
                }

                polynomial_type quotient;
                detail::factorization_exact_quotient(quotient, remaining, result.factors.back().polynomial,
                                                     arithmetic_context);
                remaining = std::move(quotient);
            }
            ++irreducible_factor_degree;
        }

        if (remaining.size() > 1) {
            const std::size_t remaining_degree = remaining.size() - 1;
            result.factors.push_back({std::move(remaining), remaining_degree});
            if (factor_callback(result.factors.back()) == factorization_control::stop_factorization) {
                result.complete = false;
            }
        }

        return result;
    }

    /** Compute the complete reference distinct-degree factorization without a staged callback. */
    template<detail::SupportsDivrem Backend>
    distinct_degree_factorization_result<typename Backend::polynomial_type> distinct_degree_factorization_reference(
        const typename Backend::polynomial_type &input,
        polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using factor_type = distinct_degree_factor<typename Backend::polynomial_type>;
        return distinct_degree_factorization_reference<Backend>(input, arithmetic_context, [](const factor_type &) {
            return factorization_control::continue_factorization;
        });
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_DISTINCT_DEGREE_FACTORIZATION_HPP
