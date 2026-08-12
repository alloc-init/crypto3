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

#ifndef CRYPTO3_MATH_POLYNOMIAL_BACKEND_HPP
#define CRYPTO3_MATH_POLYNOMIAL_BACKEND_HPP

#include <concepts>
#include <cstddef>
#include <utility>
#include <vector>

namespace nil::crypto3::math::polynomial_arithmetic {

    /**
     * Polynomial coefficients are stored in ascending degree order. A canonical
     * polynomial is nonempty and has no trailing zero coefficients, except that
     * the zero polynomial is represented by the single coefficient [0].
     */
    template<typename ValueType>
    using coefficient_vector = std::vector<ValueType>;

    /**
     * Interface for interchangeable polynomial multiplication implementations.
     * Higher-level polynomial algorithms use these operations without depending
     * on how products are computed.
     *
     * Inputs use canonical vector form and every operation produces canonical
     * output. The output may alias either input. multiply_low computes the product
     * modulo X^coefficient_count; when coefficient_count is zero, it stores [0].
     *
     * Operations are invoked on a mutable backend so implementations may update
     * reusable caches or scratch storage; stateless implementations may still
     * declare their operations const.
     */
    template<typename Backend>
    concept PolynomialBackend =
        requires { typename Backend::value_type; } &&
        requires(Backend &backend, coefficient_vector<typename Backend::value_type> &output,
                 const coefficient_vector<typename Backend::value_type> &left,
                 const coefficient_vector<typename Backend::value_type> &right, std::size_t coefficient_count) {
            { backend.multiply(output, left, right) } -> std::same_as<void>;
            { backend.square(output, left) } -> std::same_as<void>;
            { backend.multiply_low(output, left, right, coefficient_count) } -> std::same_as<void>;
        };

    /**
     * Owns the multiplication implementation used by a sequence of higher-level
     * polynomial operations. Keeping one backend alive lets those operations
     * reuse implementation-specific configuration, precomputed state, and scratch
     * storage instead of rebuilding them for every product.
     *
     * The context is the stable entry point used by higher-level algorithms, so
     * they do not handle backend-specific resources or representations. A context
     * can be reused sequentially, but callers must use separate contexts for
     * concurrent operations.
     */
    template<PolynomialBackend Backend>
    class polynomial_context {
    public:
        using backend_type = Backend;
        using value_type = typename backend_type::value_type;
        using coefficient_vector = polynomial_arithmetic::coefficient_vector<value_type>;

        polynomial_context() = default;

        explicit polynomial_context(backend_type backend) : backend_(std::move(backend)) {
        }

        void multiply(coefficient_vector &output, const coefficient_vector &left, const coefficient_vector &right) {
            backend_.multiply(output, left, right);
        }

        void square(coefficient_vector &output, const coefficient_vector &input) {
            backend_.square(output, input);
        }

        void multiply_low(coefficient_vector &output, const coefficient_vector &left, const coefficient_vector &right,
                          std::size_t coefficient_count) {
            backend_.multiply_low(output, left, right, coefficient_count);
        }

    private:
        backend_type backend_;
    };

}    // namespace nil::crypto3::math::polynomial_arithmetic

#endif    // CRYPTO3_MATH_POLYNOMIAL_BACKEND_HPP
