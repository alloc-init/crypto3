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

#ifndef CRYPTO3_MATH_MIXED_RADIX_BACKEND_HPP
#define CRYPTO3_MATH_MIXED_RADIX_BACKEND_HPP

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/math/algorithms/mixed_radix_fft.hpp>
#include <nil/crypto3/math/polynomial/basic_operations.hpp>
#include <nil/crypto3/math/polynomial/polynomial.hpp>
#include <nil/crypto3/math/polynomial/polynomial_backend.hpp>

namespace nil::crypto3::math::polynomial_arithmetic {

    /**
     * Mixed-radix multiplication backend with one explicitly sized transform plan.
     * The backend reuses its workspace across operations and is therefore intended
     * for sequential use. Concurrent operations require separate backend instances.
     *
     * The configured transform must contain the complete product of the coefficients
     * used by an operation. multiply_low discards input coefficients that cannot
     * affect the requested low coefficients before computing that product.
     */
    template<typename RootFieldType, typename ValueType = typename RootFieldType::value_type>
    class mixed_radix_backend {
    public:
        using value_type = ValueType;
        using polynomial_type = math::polynomial<value_type>;

        explicit mixed_radix_backend(std::size_t transform_size) : plan_(transform_size) {
        }

        void multiply(polynomial_type &output, const polynomial_type &left, const polynomial_type &right) {
            if (math::is_zero(left) || math::is_zero(right)) {
                set_zero(output);
                return;
            }

            const std::size_t result_size = left.size() + right.size() - 1;
            multiply_prefixes(output, left, left.size(), right, right.size(), result_size);
        }

        void square(polynomial_type &output, const polynomial_type &input) {
            if (math::is_zero(input)) {
                set_zero(output);
                return;
            }

            const std::size_t result_size = 2 * input.size() - 1;
            validate_result_size(result_size);

            polynomial_type transformed = input;
            plan_.fft(transformed.get_storage(), workspace_);
            for (value_type &value : transformed) {
                value = value * value;
            }

            finish_transform(output, transformed, result_size);
        }

        void multiply_low(polynomial_type &output, const polynomial_type &left, const polynomial_type &right,
                          std::size_t coefficient_count) {
            if (coefficient_count == 0) {
                set_zero(output);
                return;
            }

            const std::size_t left_size = std::min(left.size(), coefficient_count);
            const std::size_t right_size = std::min(right.size(), coefficient_count);
            if (math::is_zero(left.begin(), left.begin() + left_size) ||
                math::is_zero(right.begin(), right.begin() + right_size)) {
                set_zero(output);
                return;
            }

            const std::size_t prefix_product_size = left_size + right_size - 1;
            const std::size_t result_size = std::min(coefficient_count, prefix_product_size);
            multiply_prefixes(output, left, left_size, right, right_size, result_size);
        }

    private:
        static void set_zero(polynomial_type &output) {
            output.assign(1, value_type::zero());
        }

        void validate_result_size(std::size_t result_size) const {
            if (result_size > plan_.size()) {
                throw std::invalid_argument("mixed_radix_backend: transform is too small for the product");
            }
        }

        void multiply_prefixes(polynomial_type &output, const polynomial_type &left, std::size_t left_size,
                               const polynomial_type &right, std::size_t right_size, std::size_t result_size) {
            validate_result_size(left_size + right_size - 1);

            polynomial_type transformed_left(left.begin(), left.begin() + left_size);
            polynomial_type transformed_right(right.begin(), right.begin() + right_size);
            plan_.fft(transformed_left.get_storage(), workspace_);
            plan_.fft(transformed_right.get_storage(), workspace_);

            for (std::size_t i = 0; i < plan_.size(); ++i) {
                transformed_left[i] = transformed_left[i] * transformed_right[i];
            }

            finish_transform(output, transformed_left, result_size);
        }

        void finish_transform(polynomial_type &output, polynomial_type &transformed, std::size_t result_size) {
            plan_.inverse_fft(transformed.get_storage(), workspace_);
            math::truncate(transformed, result_size);
            output = std::move(transformed);
        }

        mixed_radix_fft_plan<RootFieldType> plan_;
        std::vector<value_type> workspace_;
    };

}    // namespace nil::crypto3::math::polynomial_arithmetic

#endif    // CRYPTO3_MATH_MIXED_RADIX_BACKEND_HPP
