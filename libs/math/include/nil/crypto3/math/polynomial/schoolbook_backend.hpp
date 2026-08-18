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

#ifndef CRYPTO3_MATH_SCHOOLBOOK_BACKEND_HPP
#define CRYPTO3_MATH_SCHOOLBOOK_BACKEND_HPP

#include <algorithm>
#include <cstddef>
#include <utility>

#include <nil/crypto3/math/polynomial/basic_operations.hpp>
#include <nil/crypto3/math/polynomial/polynomial.hpp>
#include <nil/crypto3/math/polynomial/polynomial_backend.hpp>

namespace nil::crypto3::math::polynomial_arithmetic {

    /**
     * Reference backend using quadratic coefficient multiplication. It provides
     * a simple correctness oracle for higher-level algorithms and avoids the
     * setup cost of faster backends for small products.
     *
     * Every operation builds a private result before assigning output, so output
     * may own storage referenced by either input view without changing the
     * computation.
     */
    template<typename ValueType>
    class schoolbook_backend {
    public:
        using value_type = ValueType;
        using polynomial_type = math::polynomial<value_type>;

        void multiply(polynomial_type &output, coefficient_view<value_type> left,
                      coefficient_view<value_type> right) const {
            const std::size_t coefficient_count = left.empty() || right.empty() ? 0 : left.size() + right.size() - 1;
            multiply_ranges_low(output, left, right, coefficient_count);
        }

        void square(polynomial_type &output, const polynomial_type &input) const {
            const coefficient_view<value_type> input_view(input);
            multiply(output, input_view, input_view);
        }

        void multiply_low(polynomial_type &output, const polynomial_type &left, const polynomial_type &right,
                          std::size_t coefficient_count) const {
            multiply_ranges_low(output, coefficient_view<value_type>(left), coefficient_view<value_type>(right),
                                coefficient_count);
        }

    private:
        static void multiply_ranges_low(polynomial_type &output, coefficient_view<value_type> left,
                                        coefficient_view<value_type> right, std::size_t coefficient_count) {
            if (coefficient_count == 0 || left.empty() || right.empty()) {
                output.assign(1, value_type::zero());
                return;
            }

            const std::size_t product_size = left.size() + right.size() - 1;
            const std::size_t result_size = std::min(coefficient_count, product_size);

            polynomial_type result(result_size, value_type::zero());
            const std::size_t left_count = std::min(left.size(), result_size);
            for (std::size_t left_index = 0; left_index < left_count; ++left_index) {
                if (left[left_index].is_zero()) {
                    continue;
                }
                const std::size_t right_count = std::min(right.size(), result_size - left_index);
                for (std::size_t right_index = 0; right_index < right_count; ++right_index) {
                    result[left_index + right_index] += left[left_index] * right[right_index];
                }
            }

            math::condense(result);
            output = std::move(result);
        }
    };

}    // namespace nil::crypto3::math::polynomial_arithmetic

#endif    // CRYPTO3_MATH_SCHOOLBOOK_BACKEND_HPP
