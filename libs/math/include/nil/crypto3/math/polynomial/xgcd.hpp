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

#ifndef CRYPTO3_MATH_XGCD_HPP
#define CRYPTO3_MATH_XGCD_HPP

#include <algorithm>
#include <cassert>
#include <concepts>
#include <ranges>
#include <vector>

#include <nil/crypto3/algebra/type_traits.hpp>
#include <nil/crypto3/math/polynomial/basic_operations.hpp>

namespace nil {
    namespace crypto3 {
        namespace math {

            /*!
             * @brief Perform the standard Extended Euclidean Division algorithm.
             * Input: Polynomial A, Polynomial B.
             * Output: Polynomial G, Polynomial U, Polynomial V, such that G = (A * U) + (B * V).
             * This implementation uses quadratic polynomial division and the legacy radix-2 FFT multiplication path.
             *
             * @pre A and B are nonempty canonical coefficient ranges.
             * @pre G, U, and V are distinct output objects.
             */
            template<detail::PolynomialCoefficientRange Range1, detail::PolynomialCoefficientRange Range2,
                     detail::MutablePolynomialCoefficientRange Range3, detail::MutablePolynomialCoefficientRange Range4,
                     detail::MutablePolynomialCoefficientRange Range5>
                requires std::same_as<std::ranges::range_value_t<const Range1>,
                                      std::ranges::range_value_t<const Range2>> &&
                         std::same_as<std::ranges::range_value_t<const Range1>, std::ranges::range_value_t<Range3>> &&
                         std::same_as<std::ranges::range_value_t<const Range1>, std::ranges::range_value_t<Range4>> &&
                         std::same_as<std::ranges::range_value_t<const Range1>, std::ranges::range_value_t<Range5>> &&
                         algebra::is_field_element<std::ranges::range_value_t<const Range1>>::value &&
                         requires(const Range1 &input, Range3 &g, Range4 &u, Range5 &v,
                                  const std::vector<std::ranges::range_value_t<const Range1>> &coefficients) {
                             g = input;
                             g = coefficients;
                             u = coefficients;
                             v = coefficients;
                         }
            void extended_euclidean(const Range1 &a, const Range2 &b, Range3 &g, Range4 &u, Range5 &v) {

                using value_type = std::ranges::range_value_t<const Range1>;

                if (is_zero(b)) {
                    g = a;
                    u = std::vector<value_type>(1, value_type::one());
                    v = std::vector<value_type>(1, value_type::zero());
                    return;
                }

                std::vector<value_type> previous_a_coefficient(1, value_type::one());
                std::vector<value_type> a_coefficient(1, value_type::zero());
                std::vector<value_type> previous_remainder(a);
                std::vector<value_type> remainder(b);

                std::vector<value_type> quotient(1, value_type::zero());
                std::vector<value_type> next_remainder(1, value_type::zero());
                std::vector<value_type> product(1, value_type::zero());
                std::vector<value_type> next_a_coefficient(1, value_type::zero());

                while (!is_zero(remainder)) {
                    division(quotient, next_remainder, previous_remainder, remainder);
                    multiplication(product, a_coefficient, quotient);
                    subtraction(next_a_coefficient, previous_a_coefficient, product);

                    previous_a_coefficient = a_coefficient;
                    previous_remainder = remainder;
                    a_coefficient = next_a_coefficient;
                    remainder = next_remainder;
                }

                // Recover the coefficient of b from G = a * U + b * V once the Euclidean loop has found G and U.
                multiplication(product, a, previous_a_coefficient);
                subtraction(product, previous_remainder, product);
                std::vector<value_type> b_coefficient(1, value_type::zero());
                division(b_coefficient, next_remainder, product, b);
                assert(is_zero(next_remainder));

                const value_type inverse_leading_coefficient = previous_remainder.back().inversed();
                const auto scale_to_monic = [&inverse_leading_coefficient](const value_type &coefficient) {
                    return coefficient * inverse_leading_coefficient;
                };
                std::transform(previous_remainder.begin(), previous_remainder.end(), previous_remainder.begin(),
                               scale_to_monic);
                std::transform(previous_a_coefficient.begin(), previous_a_coefficient.end(),
                               previous_a_coefficient.begin(), scale_to_monic);
                std::transform(b_coefficient.begin(), b_coefficient.end(), b_coefficient.begin(), scale_to_monic);

                g = previous_remainder;
                u = previous_a_coefficient;
                v = b_coefficient;
            }
        }    // namespace math
    }    // namespace crypto3
}    // namespace nil

#endif    // ALGEBRA_FFT_XGCD_HPP
