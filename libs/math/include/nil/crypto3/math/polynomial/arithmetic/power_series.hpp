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

#ifndef CRYPTO3_MATH_POLYNOMIAL_POWER_SERIES_HPP
#define CRYPTO3_MATH_POLYNOMIAL_POWER_SERIES_HPP

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <utility>

#include <nil/crypto3/math/polynomial/operations/basic_operations.hpp>

namespace nil::crypto3::math {

    /**
     * Compute the inverse of input modulo X^coefficient_count using Newton iteration. The constant coefficient of
     * input must be nonzero. The output is canonical and may alias input; a coefficient count of zero produces [0].
     *
     * Let f be input and let g_m satisfy f * g_m = 1 mod X^m. Write
     *
     *     f * g_m = 1 + X^m * e_m.
     *
     * Starting with g_1 = f[0]^-1, an iteration extending the precision from m to n, where m < n <= 2m, retains the
     * first m coefficients and computes only the new high block:
     *
     *     g_n = g_m - X^m * (g_m * e_m mod X^(n-m)).
     *
     * Substitution gives f * g_n = 1 mod X^n. Computing only this high block avoids repeating the already-correct low
     * coefficients in the second multiplication. The final iteration is truncated when coefficient_count is not a
     * power of two.
     *
     * @throws std::invalid_argument if coefficient_count is nonzero and input has zero constant coefficient.
     */
    template<polynomial_arithmetic::PolynomialBackend Backend>
        requires detail::MutableNormalizableCoefficientPolynomial<typename Backend::polynomial_type> &&
                 std::default_initializable<typename Backend::polynomial_type> &&
                 std::movable<typename Backend::polynomial_type> &&
                 requires(const typename Backend::polynomial_type::value_type &left,
                          const typename Backend::polynomial_type::value_type &right) {
                     {
                         Backend::polynomial_type::value_type::one()
                     } -> std::convertible_to<typename Backend::polynomial_type::value_type>;
                     { left.inversed() } -> std::convertible_to<typename Backend::polynomial_type::value_type>;
                     { left - right } -> std::convertible_to<typename Backend::polynomial_type::value_type>;
                 }
    void inverse_series(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input,
                        std::size_t coefficient_count, polynomial_arithmetic::polynomial_context<Backend> &context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        if (coefficient_count == 0) {
            output.resize(1);
            output[0] = value_type {};
            return;
        }
        if (input.size() == 0 || input[0] == value_type {}) {
            throw std::invalid_argument("a power series with zero constant coefficient is not invertible");
        }

        polynomial_type approximation;
        approximation.resize(1);
        approximation[0] = input[0].inversed();

        polynomial_type product;
        polynomial_type high_error;
        polynomial_type new_coefficients;

        std::size_t precision = 1;
        while (precision < coefficient_count) {
            const std::size_t next_precision = precision + std::min(precision, coefficient_count - precision);
            const std::size_t added_coefficient_count = next_precision - precision;

            // The first precision coefficients of input * approximation are already [1, 0, ..., 0]. The following
            // block is e_m from input * approximation = 1 + X^precision * e_m.
            multiply_low(product, input, approximation, next_precision, context);
            product.resize(next_precision, value_type {});
            high_error.resize(added_coefficient_count);
            for (std::size_t i = 0; i < added_coefficient_count; ++i) {
                high_error[i] = product[precision + i];
            }
            condense(high_error);

            // -approximation * high_error gives exactly the new coefficients to append; the known low block is kept.
            multiply_low(new_coefficients, approximation, high_error, added_coefficient_count, context);
            new_coefficients.resize(added_coefficient_count, value_type {});
            approximation.resize(next_precision, value_type {});
            for (std::size_t i = 0; i < added_coefficient_count; ++i) {
                approximation[precision + i] = value_type {} - new_coefficients[i];
            }
            condense(approximation);
            precision = next_precision;
        }

        output = std::move(approximation);
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_POWER_SERIES_HPP
