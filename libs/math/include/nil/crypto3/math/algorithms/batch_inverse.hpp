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

#ifndef CRYPTO3_MATH_BATCH_INVERSE_HPP
#define CRYPTO3_MATH_BATCH_INVERSE_HPP

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace nil {
    namespace crypto3 {
        namespace math {

            /**
             * Compute the multiplicative inverse of every nonzero value using Montgomery's batch-inversion trick.
             *
             * For n > 0 this performs one field inversion and exactly 3 * (n - 1) field multiplications. The result
             * vector doubles as prefix-product scratch space, so no additional scratch allocation is required.
             *
             * @throws std::invalid_argument if any input value is zero.
             */
            template<typename ValueType, typename Allocator>
            [[nodiscard]] std::vector<ValueType, Allocator>
                batch_inverse_nonzero(const std::vector<ValueType, Allocator> &values) {
                std::vector<ValueType, Allocator> result(values.get_allocator());
                if (values.empty()) {
                    return result;
                }

                result.reserve(values.size());

                if (values.front().is_zero()) {
                    throw std::invalid_argument("batch_inverse_nonzero: input values must be nonzero");
                }
                result.emplace_back(values.front());

                for (std::size_t i = 1; i < values.size(); ++i) {
                    if (values[i].is_zero()) {
                        throw std::invalid_argument("batch_inverse_nonzero: input values must be nonzero");
                    }
                    result.emplace_back(result.back() * values[i]);
                }

                ValueType inverse = result.back().inversed();
                for (std::size_t i = values.size() - 1; i > 0; --i) {
                    result[i] = result[i - 1] * inverse;
                    inverse = inverse * values[i];
                }
                result[0] = inverse;

                return result;
            }

        }    // namespace math
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MATH_BATCH_INVERSE_HPP
