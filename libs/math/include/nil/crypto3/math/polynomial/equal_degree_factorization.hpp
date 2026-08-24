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

#ifndef CRYPTO3_MATH_EQUAL_DEGREE_FACTORIZATION_HPP
#define CRYPTO3_MATH_EQUAL_DEGREE_FACTORIZATION_HPP

#include <cstddef>
#include <stdexcept>

#include <nil/crypto3/math/polynomial/polynomial_factorization.hpp>

namespace nil::crypto3::math {

    namespace detail {

        /**
         * Prepare one group produced by distinct-degree factorization for equal-degree splitting. The input is
         * normalized to monic form and its original leading coefficient is returned separately. Constant inputs
         * return false because they contain no factors.
         *
         * Equal-degree factorization requires a square-free input whose irreducible factors all have
         * irreducible_factor_degree. This function checks square-freeness and the necessary divisibility of the total
         * degree. The caller must supply a group produced by distinct-degree factorization. Verifying the degree of
         * every irreducible factor here would require repeating that factorization, so this function does not perform
         * that check.
         *
         * For example, suppose distinct-degree factorization produces the group G = Q1 * Q2 * Q3, where Q1, Q2,
         * and Q3 are distinct irreducible quadratic polynomials. Equal-degree factorization of G with
         * irreducible_factor_degree = 2 separates that group into Q1, Q2, and Q3. Splitting stops when a resulting
         * piece has degree two: under the equal-degree precondition, that piece is one of Q1, Q2, or Q3 and is already
         * irreducible.
         *
         * @throws std::invalid_argument if irreducible_factor_degree is zero, the input is not square-free, or its
         * total degree is not divisible by irreducible_factor_degree.
         * @pre input is a nonempty coefficient polynomial.
         */
        template<SupportsDivrem Backend>
        bool prepare_equal_degree_factorization_input(
            typename Backend::polynomial_type &monic_input,
            typename Backend::polynomial_type::value_type &leading_coefficient,
            const typename Backend::polynomial_type &input, std::size_t irreducible_factor_degree,
            polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
            if (irreducible_factor_degree == 0) {
                throw std::invalid_argument("equal-degree factorization requires a positive factor degree");
            }
            if (!prepare_square_free_factorization_input<Backend>(monic_input, leading_coefficient, input,
                                                                  arithmetic_context)) {
                return false;
            }
            if ((monic_input.size() - 1) % irreducible_factor_degree != 0) {
                throw std::invalid_argument(
                    "equal-degree factorization requires the factor degree to divide the polynomial degree");
            }
            return true;
        }

    }    // namespace detail

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_EQUAL_DEGREE_FACTORIZATION_HPP
