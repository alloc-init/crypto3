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

#ifndef CRYPTO3_MATH_POLYNOMIAL_FACTORIZATION_HPP
#define CRYPTO3_MATH_POLYNOMIAL_FACTORIZATION_HPP

#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <nil/crypto3/math/polynomial/concepts.hpp>
#include <nil/crypto3/math/polynomial/gcd.hpp>
#include <nil/crypto3/math/polynomial/polynomial_division.hpp>

namespace nil::crypto3::math {

    /** The action requested by a callback after a factorization stage produces a factor. */
    enum class factorization_control { continue_factorization, stop_factorization };

    /** A polynomial factor together with its positive multiplicity in the input polynomial. */
    template<CoefficientPolynomial Polynomial>
    struct polynomial_factor {
        Polynomial polynomial;
        std::size_t multiplicity = 1;

        bool operator==(const polynomial_factor &) const = default;
    };

    /**
     * A collected polynomial factorization. A complete result satisfies
     *
     *     input = leading_coefficient * product(factor.polynomial ^ factor.multiplicity).
     *
     * Nonconstant factors are stored in monic canonical form. When a staged callback requests an early stop, factors
     * contains the produced prefix, including the factor that caused the stop, and complete is false.
     */
    template<CoefficientPolynomial Polynomial>
    struct polynomial_factorization_result {
        using polynomial_type = Polynomial;
        using value_type = typename polynomial_type::value_type;
        using factor_type = polynomial_factor<polynomial_type>;

        value_type leading_coefficient {};
        std::vector<factor_type> factors;
        bool complete = true;

        bool operator==(const polynomial_factorization_result &) const = default;
    };

    /**
     * One group produced by distinct-degree factorization. polynomial is the product of all irreducible input factors
     * whose degree equals irreducible_factor_degree; it need not itself be irreducible.
     */
    template<CoefficientPolynomial Polynomial>
    struct distinct_degree_factor {
        Polynomial polynomial;
        std::size_t irreducible_factor_degree = 1;

        bool operator==(const distinct_degree_factor &) const = default;
    };

    /**
     * A collected distinct-degree factorization. A complete result satisfies
     *
     *     input = leading_coefficient * product(factor.polynomial).
     *
     * Each stored polynomial is monic and square-free, and contains all input factors having the associated
     * irreducible_factor_degree. A stopped result contains the produced prefix, including the factor that caused the
     * stop, and has complete set to false.
     */
    template<CoefficientPolynomial Polynomial>
    struct distinct_degree_factorization_result {
        using polynomial_type = Polynomial;
        using value_type = typename polynomial_type::value_type;
        using factor_type = distinct_degree_factor<polynomial_type>;

        value_type leading_coefficient {};
        std::vector<factor_type> factors;
        bool complete = true;

        bool operator==(const distinct_degree_factorization_result &) const = default;
    };

    namespace detail {

        /** A callback that controls staged factorization after receiving one polynomial factor. */
        template<typename FactorCallback, typename Polynomial>
        concept PolynomialFactorCallback =
            CoefficientPolynomial<Polynomial> &&
            requires(FactorCallback &callback, const polynomial_factor<Polynomial> &factor) {
                { callback(factor) } -> std::same_as<factorization_control>;
            };

        /**
         * Normalize a factorization input and verify that every irreducible factor has multiplicity one. The
         * original leading coefficient is returned separately. Constant inputs return false; nonconstant
         * square-free inputs return true and monic_input receives their monic form.
         *
         * @throws std::invalid_argument if a nonconstant input is not square-free.
         * @pre input is a nonempty coefficient polynomial.
         */
        template<SupportsDivrem Backend>
        bool prepare_square_free_factorization_input(
            typename Backend::polynomial_type &monic_input,
            typename Backend::polynomial_type::value_type &leading_coefficient,
            const typename Backend::polynomial_type &input,
            polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
            monic_input = input;
            condense(monic_input);
            leading_coefficient = monic_input.back();
            if (monic_input.size() == 1) {
                return false;
            }
            make_monic(monic_input, monic_input);

            typename Backend::polynomial_type input_derivative;
            derivative(input_derivative, monic_input);
            typename Backend::polynomial_type repeated_factor;
            gcd(repeated_factor, monic_input, input_derivative, arithmetic_context);
            if (repeated_factor.size() > 1) {
                throw std::invalid_argument("factorization requires a square-free polynomial");
            }
            return true;
        }

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

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_FACTORIZATION_HPP
