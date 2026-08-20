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

#ifndef CRYPTO3_MATH_SQUARE_FREE_FACTORIZATION_HPP
#define CRYPTO3_MATH_SQUARE_FREE_FACTORIZATION_HPP

#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <utility>

#include <nil/crypto3/algebra/fields/field_order.hpp>

#include <nil/crypto3/math/polynomial/gcd.hpp>
#include <nil/crypto3/math/polynomial/polynomial_factorization.hpp>

namespace nil::crypto3::math {
    namespace detail {

        template<SupportsDivrem Backend>
        void factorization_exact_quotient(typename Backend::polynomial_type &output,
                                          const typename Backend::polynomial_type &dividend,
                                          const typename Backend::polynomial_type &divisor,
                                          polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
            const std::size_t quotient_coefficient_count =
                dividend.size() >= divisor.size() ? dividend.size() - divisor.size() + 1 : 1;
            polynomial_divisor_context<Backend> divisor_context(divisor, quotient_coefficient_count,
                                                                arithmetic_context);
            exact_division(output, dividend, divisor_context, arithmetic_context);
        }

    }    // namespace detail

    /**
     * Separate a polynomial into monic square-free factors with distinct multiplicities using Yun's decomposition.
     * For a nonconstant input, the returned leading coefficient is the input's original leading coefficient. By
     * convention, the zero polynomial returns scalar zero and no factors; a nonzero constant returns its sole
     * coefficient and no factors.
     *
     * Write a nonconstant monic input as F = product(P_j^m_j), where the P_j are distinct irreducible factors. Each
     * summand of the product-rule expansion of F' differentiates one factor, but every summand remains divisible by
     * P_j^(m_j-1) for every j. Under the characteristic restriction below, G = gcd(F, F') therefore contains exactly
     * m_j-1 copies of each P_j. Consequently W = F/G contains one copy of every distinct factor.
     * In formulas
     * F = \prod_j P_j^m_j(x)
     * F'= \{\prod_{k}P_k^{m_k-1}\}\, \sum_j m_j P'_j(x) \prod_{i\neq j} P_i(x)
     * G = gcd(F,F') = \prod_{k}P_k^{m_k-1}
     * W = F/G = \prod_k P_k(x)
     *
     * The algorithm then loops over factor multiplicities i = 1, 2, 3, ... and stops when W = 1. At the start of the
     * iteration for multiplicity i, W contains one copy of each factor whose original multiplicity is at least i, and G
     * contains its remaining copies. Set Y = gcd(W, G): Y retains the factors whose multiplicity is greater than i, so
     * Z = W/Y contains exactly the factors of multiplicity i. Emit Z when it is nonconstant, replace W with Y and G
     * with G/Y, and continue with the next multiplicity.
     * In formulas, at multiplicity i,
     * W_i = \prod_{m_j >= i} P_j,
     * G_i = \prod_{m_j >= i} P_j^(m_j-i),
     * Y_i = gcd(W_i, G_i) = \prod_{m_j > i} P_j,
     * Z_i = W_i/Y_i = \prod_{m_j = i} P_j.
     * Thus Y_i equals W_i only when no factor has multiplicity exactly i; in that case Z_i = 1 and nothing is emitted.
     *
     * For example, for the non-square-free input F = A * B^2 * C^3, initially W = A * B * C and G = B * C^2. The
     * first iteration emits A and leaves W = B * C, G = C; the second emits B and leaves W = C, G = 1; the third emits
     * C.
     *
     * This initial implementation requires the coefficient-field characteristic to be greater than degree(F). Under
     * that restriction, differentiation cannot erase a nonconstant p-th-power component, so the recurrence accounts
     * for every factor without extracting polynomial p-th roots. Later support for small characteristic can remove
     * this restriction.
     *
     * Factors are emitted in increasing multiplicity. After each nonconstant factor is appended to the result,
     * factor_callback may request an early stop. A stopped result includes that factor and has complete set to false.
     *
     * @throws std::invalid_argument if the coefficient-field characteristic is not greater than the input degree.
     */
    template<detail::SupportsDivrem Backend, typename FactorCallback>
        requires requires(FactorCallback &callback,
                          const polynomial_factor<typename Backend::polynomial_type> &factor) {
            { callback(factor) } -> std::same_as<factorization_control>;
        }
    polynomial_factorization_result<typename Backend::polynomial_type>
        square_free_factorization(const typename Backend::polynomial_type &input,
                                  polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                  FactorCallback &&factor_callback) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;
        using field_type = typename value_type::field_type;
        using result_type = polynomial_factorization_result<polynomial_type>;

        result_type result;
        polynomial_type monic_input(input);
        condense(monic_input);
        result.leading_coefficient = monic_input.back();

        if (monic_input.size() == 1) {
            return result;
        }

        const std::size_t input_degree = monic_input.size() - 1;
        if (algebra::fields::field_characteristic<field_type>() <= input_degree) {
            throw std::invalid_argument(
                "square-free factorization requires characteristic greater than the polynomial degree");
        }
        make_monic(monic_input, monic_input);

        polynomial_type input_derivative;
        derivative(input_derivative, monic_input);

        polynomial_type repeated_part;
        gcd(repeated_part, monic_input, input_derivative, arithmetic_context);

        polynomial_type square_free_part;
        detail::factorization_exact_quotient(square_free_part, monic_input, repeated_part, arithmetic_context);

        std::size_t multiplicity = 1;
        while (square_free_part.size() > 1) {
            polynomial_type next_square_free_part;
            gcd(next_square_free_part, square_free_part, repeated_part, arithmetic_context);

            polynomial_type factor;
            detail::factorization_exact_quotient(factor, square_free_part, next_square_free_part, arithmetic_context);
            if (factor.size() > 1) {
                result.factors.push_back({std::move(factor), multiplicity});
                if (factor_callback(result.factors.back()) == factorization_control::stop_factorization) {
                    result.complete = false;
                    return result;
                }
            }

            detail::factorization_exact_quotient(repeated_part, repeated_part, next_square_free_part,
                                                 arithmetic_context);
            square_free_part = std::move(next_square_free_part);
            ++multiplicity;
        }

        return result;
    }

    /** Compute the complete square-free factorization without a staged callback. */
    template<detail::SupportsDivrem Backend>
    polynomial_factorization_result<typename Backend::polynomial_type>
        square_free_factorization(const typename Backend::polynomial_type &input,
                                  polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using factor_type = polynomial_factor<typename Backend::polynomial_type>;
        return square_free_factorization<Backend>(input, arithmetic_context, [](const factor_type &) {
            return factorization_control::continue_factorization;
        });
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_SQUARE_FREE_FACTORIZATION_HPP
