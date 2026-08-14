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

#include <nil/crypto3/math/algorithms/mixed_radix_fft.hpp>
#include <nil/crypto3/math/polynomial/basic_operations.hpp>
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
        using coefficient_vector = polynomial_arithmetic::coefficient_vector<value_type>;

        explicit mixed_radix_backend(std::size_t transform_size) : plan_(transform_size) {
        }

        void multiply(coefficient_vector &output, const coefficient_vector &left, const coefficient_vector &right) {
            if (math::is_zero(left) || math::is_zero(right)) {
                set_zero(output);
                return;
            }

            const std::size_t result_size = left.size() + right.size() - 1;
            multiply_prefixes(output, left, left.size(), right, right.size(), result_size);
        }

        void square(coefficient_vector &output, const coefficient_vector &input) {
            if (math::is_zero(input)) {
                set_zero(output);
                return;
            }

            const std::size_t result_size = 2 * input.size() - 1;
            validate_result_size(result_size);

            coefficient_vector transformed = input;
            plan_.fft(transformed, workspace_);
            for (value_type &value : transformed) {
                value = value * value;
            }

            finish_transform(output, transformed, result_size);
        }

        void multiply_low(coefficient_vector &output, const coefficient_vector &left, const coefficient_vector &right,
                          std::size_t coefficient_count) {
            if (coefficient_count == 0) {
                set_zero(output);
                return;
            }
            if (math::is_zero(left) || math::is_zero(right)) {
                set_zero(output);
                return;
            }

            const std::size_t left_size = std::min(left.size(), coefficient_count);
            const std::size_t right_size = std::min(right.size(), coefficient_count);
            const std::size_t prefix_product_size = left_size + right_size - 1;
            const std::size_t result_size = std::min(coefficient_count, prefix_product_size);
            multiply_prefixes(output, left, left_size, right, right_size, result_size);
        }

    private:
        static void set_zero(coefficient_vector &output) {
            output.assign(1, value_type::zero());
        }

        void validate_result_size(std::size_t result_size) const {
            if (result_size > plan_.size()) {
                throw std::invalid_argument("mixed_radix_backend: transform is too small for the product");
            }
        }

        void multiply_prefixes(coefficient_vector &output, const coefficient_vector &left, std::size_t left_size,
                               const coefficient_vector &right, std::size_t right_size, std::size_t result_size) {
            validate_result_size(left_size + right_size - 1);

            coefficient_vector transformed_left(left.begin(), left.begin() + left_size);
            coefficient_vector transformed_right(right.begin(), right.begin() + right_size);
            plan_.fft(transformed_left, workspace_);
            plan_.fft(transformed_right, workspace_);

            for (std::size_t i = 0; i < plan_.size(); ++i) {
                transformed_left[i] = transformed_left[i] * transformed_right[i];
            }

            finish_transform(output, transformed_left, result_size);
        }

        void finish_transform(coefficient_vector &output, coefficient_vector &transformed, std::size_t result_size) {
            plan_.inverse_fft(transformed, workspace_);
            transformed.resize(result_size);
            math::condense(transformed);
            output = std::move(transformed);
        }

        mixed_radix_fft_plan<RootFieldType> plan_;
        coefficient_vector workspace_;
    };

}    // namespace nil::crypto3::math::polynomial_arithmetic

#endif    // CRYPTO3_MATH_MIXED_RADIX_BACKEND_HPP
