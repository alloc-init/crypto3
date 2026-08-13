//---------------------------------------------------------------------------//
// Copyright (c) 2020-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
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

#ifndef CRYPTO3_MATH_EVALUATE_HPP
#define CRYPTO3_MATH_EVALUATE_HPP

#include <concepts>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <vector>

#include <boost/math/tools/polynomial.hpp>

#include <nil/crypto3/algebra/type_traits.hpp>

namespace nil {
    namespace crypto3 {
        namespace math {
            /*!
             * @brief
             * Naive evaluation of a *single* polynomial, used for testing purposes.
             *
             * The inputs are:
             * - an integer m
             * - a vector coeff representing monomial P of size m
             * - a field element element t
             * The output is the polynomial P(x) evaluated at x = t.
             *
             * @pre m > 0.
             * @pre [first, last) contains exactly m coefficients.
             */
            template<typename FieldValueType, std::contiguous_iterator ContiguousIterator>
                requires std::same_as<std::iter_value_t<ContiguousIterator>, FieldValueType>
            inline FieldValueType evaluate_polynomial(ContiguousIterator first, ContiguousIterator last,
                                                      const FieldValueType &t, std::size_t m) {
                BOOST_ASSERT(std::size_t(std::distance(first, last)) == m);

                return boost::math::tools::evaluate_polynomial(&*first, t, m);
            }

            template<typename FieldValueType, typename ContiguousContainer>
                requires std::ranges::contiguous_range<const ContiguousContainer> &&
                         std::same_as<std::ranges::range_value_t<const ContiguousContainer>, FieldValueType>
            inline FieldValueType evaluate_polynomial(const ContiguousContainer &coeff, const FieldValueType &t,
                                                      std::size_t m) {
                return evaluate_polynomial(std::ranges::begin(coeff), std::ranges::end(coeff), t, m);
            }

            /*!
             * @brief
             * Naive evaluation of a *single* Lagrange polynomial, used for testing purposes.
             *
             * The inputs are:
             * - an integer m
             * - a domain S = (a_{0},...,a_{m-1}) of size m
             * - a field element element t
             * - an index idx in {0,...,m-1}
             * The output is the polynomial L_{idx,S}(z) evaluated at z = t.
             *
             * @pre The points in [first, last) are pairwise distinct.
             */
            template<typename FieldValueType, std::random_access_iterator InputIterator>
                requires std::same_as<std::iter_value_t<InputIterator>, FieldValueType> &&
                         algebra::is_field_element<FieldValueType>::value
            inline FieldValueType evaluate_lagrange_polynomial(InputIterator first, InputIterator last,
                                                               const FieldValueType &t, std::size_t m,
                                                               std::size_t idx) {
                typedef typename std::iterator_traits<InputIterator>::value_type value_type;

                if (m != std::size_t(std::distance(first, last))) {
                    throw std::invalid_argument("expected m == domain.size()");
                }
                if (idx >= m) {
                    throw std::invalid_argument("expected idx < m");
                }

                value_type num = value_type::one();
                value_type denom = value_type::one();

                for (std::size_t k = 0; k < m; ++k) {
                    if (k == idx) {
                        continue;
                    }

                    num *= t - *(first + k);
                    denom *= *(first + idx) - *(first + k);
                }

                return num * denom.inversed();
            }

            template<typename FieldValueType, typename Range>
                requires std::ranges::random_access_range<const Range> &&
                         std::same_as<std::ranges::range_value_t<const Range>, FieldValueType> &&
                         algebra::is_field_element<FieldValueType>::value
            inline FieldValueType evaluate_lagrange_polynomial(const Range &domain, const FieldValueType &t,
                                                               std::size_t m, std::size_t idx) {
                return evaluate_lagrange_polynomial(std::ranges::begin(domain), std::ranges::end(domain), t, m, idx);
            }
        }    // namespace math
    }    // namespace crypto3
}    // namespace nil

#endif    // ALGEBRA_FFT_NAIVE_EVALUATE_HPP
