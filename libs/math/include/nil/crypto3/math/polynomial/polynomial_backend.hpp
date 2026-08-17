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

#include <nil/crypto3/math/polynomial/concepts.hpp>

namespace nil::crypto3::math::polynomial_arithmetic {

    /**
     * Algorithm-selection parameters shared by higher-level polynomial operations.
     * A zero cutoff disables the corresponding basecase-division criterion.
     */
    struct polynomial_context_options {
        std::size_t basecase_divisor_coefficient_cutoff = 10;
        std::size_t basecase_quotient_coefficient_cutoff = 2;
    };

    /**
     * Interface for interchangeable polynomial multiplication implementations.
     * Higher-level polynomial algorithms use these operations without depending
     * on how products are computed.
     *
     * The associated polynomial type must use the coefficient representation.
     * Inputs use canonical form and every operation produces canonical
     * output. The output may alias either input. multiply_low computes the product
     * modulo X^coefficient_count; when coefficient_count is zero, it stores [0].
     *
     * Operations are invoked on a mutable backend so implementations may update
     * reusable caches or scratch storage; stateless implementations may still
     * declare their operations const.
     */
    template<typename Backend>
    concept PolynomialBackend =
        requires {
            typename Backend::polynomial_type;
            requires CoefficientPolynomial<typename Backend::polynomial_type>;
        } && requires(Backend &backend, typename Backend::polynomial_type &output,
                      const typename Backend::polynomial_type &left, const typename Backend::polynomial_type &right,
                      std::size_t coefficient_count) {
            { backend.multiply(output, left, right) } -> std::same_as<void>;
            { backend.square(output, left) } -> std::same_as<void>;
            { backend.multiply_low(output, left, right, coefficient_count) } -> std::same_as<void>;
        };

    /**
     * Owns the multiplication implementation and algorithm-selection parameters
     * used by a sequence of higher-level polynomial operations. Keeping one backend
     * alive lets those operations reuse implementation-specific configuration,
     * precomputed state, and scratch storage instead of rebuilding them for every
     * product.
     *
     * Higher-level algorithms use the context without managing backend-specific
     * plans, configuration, or scratch storage. A context can be reused sequentially,
     * but callers must use separate contexts for concurrent operations.
     */
    template<PolynomialBackend Backend>
    class polynomial_context {
    public:
        using backend_type = Backend;
        using polynomial_type = typename backend_type::polynomial_type;
        using value_type = typename polynomial_type::value_type;
        using options_type = polynomial_context_options;

        polynomial_context() = default;

        explicit polynomial_context(backend_type backend, options_type options = {}) :
            backend_(std::move(backend)), options_(options) {
        }

        void multiply(polynomial_type &output, const polynomial_type &left, const polynomial_type &right) {
            backend_.multiply(output, left, right);
        }

        void square(polynomial_type &output, const polynomial_type &input) {
            backend_.square(output, input);
        }

        void multiply_low(polynomial_type &output, const polynomial_type &left, const polynomial_type &right,
                          std::size_t coefficient_count) {
            backend_.multiply_low(output, left, right, coefficient_count);
        }

        const options_type &options() const {
            return options_;
        }

    private:
        backend_type backend_;
        options_type options_;
    };

}    // namespace nil::crypto3::math::polynomial_arithmetic

#endif    // CRYPTO3_MATH_POLYNOMIAL_BACKEND_HPP
