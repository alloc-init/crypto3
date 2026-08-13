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

#ifndef CRYPTO3_MATH_POLYNOMIAL_CONCEPTS_HPP
#define CRYPTO3_MATH_POLYNOMIAL_CONCEPTS_HPP

#include <concepts>
#include <type_traits>

namespace nil::crypto3::math {

    /** Identifies a polynomial whose indexed values are coefficients. */
    struct coefficient_representation { };

    /** Identifies a polynomial whose indexed values are evaluations. */
    struct evaluation_representation { };

    namespace detail {
        template<typename T>
        using unqualified_polynomial_t = std::remove_cvref_t<T>;

        template<typename T, typename Representation>
        concept PolynomialInRepresentation = requires(const unqualified_polynomial_t<T>& polynomial,
                                                      typename unqualified_polynomial_t<T>::size_type index) {
            typename unqualified_polynomial_t<T>::value_type;
            typename unqualified_polynomial_t<T>::representation_type;
            { polynomial.size() } -> std::convertible_to<typename unqualified_polynomial_t<T>::size_type>;
            { polynomial.degree() } -> std::convertible_to<typename unqualified_polynomial_t<T>::size_type>;
            { polynomial[index] } -> std::convertible_to<typename unqualified_polynomial_t<T>::value_type>;
        } && std::same_as<typename unqualified_polynomial_t<T>::representation_type, Representation>;
    }    // namespace detail

    /** A readable polynomial represented by coefficients in ascending degree order. */
    template<typename T>
    concept CoefficientPolynomial = detail::PolynomialInRepresentation<T, coefficient_representation>;

    /** A readable polynomial represented by evaluations over a domain. */
    template<typename T>
    concept EvaluationPolynomial = detail::PolynomialInRepresentation<T, evaluation_representation>;

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_CONCEPTS_HPP
