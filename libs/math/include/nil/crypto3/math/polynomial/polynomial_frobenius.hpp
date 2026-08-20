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

#ifndef CRYPTO3_MATH_POLYNOMIAL_FROBENIUS_HPP
#define CRYPTO3_MATH_POLYNOMIAL_FROBENIUS_HPP

#include <algorithm>
#include <cstddef>
#include <utility>

#include <nil/crypto3/algebra/fields/field_order.hpp>

#include <nil/crypto3/math/polynomial/polynomial_composition.hpp>
#include <nil/crypto3/math/polynomial/polynomial_exponentiation.hpp>

namespace nil::crypto3::math {

    /**
     * Immutable precomputation for Frobenius maps in the quotient ring K[X]/(B), where K is the finite coefficient
     * field and B is a nonzero polynomial. If Q is the number of elements of K, the context stores
     *
     *     X^Q mod B.
     *
     * For A with coefficients in K, raising A to the Q-th power fixes every coefficient, so
     *
     *     A(X)^Q mod B = A(X^Q mod B) mod B.
     *
     * Consequently, constructing the context performs one modular exponentiation and precomputes the Brent-Kung
     * powers of X^Q mod B. Each later Frobenius map reuses those powers for modular composition. B need not be
     * irreducible.
     *
     * @throws std::invalid_argument if B is the zero polynomial or the configured modular-composition cached-power
     *         limit is zero.
     */
    template<detail::SupportsDivrem Backend>
    class polynomial_frobenius_context {
    public:
        using backend_type = Backend;
        using polynomial_type = typename backend_type::polynomial_type;
        using value_type = typename polynomial_type::value_type;
        using field_type = typename value_type::field_type;

        polynomial_frobenius_context(const polynomial_type &divisor,
                                     polynomial_arithmetic::polynomial_context<backend_type> &arithmetic_context) :
            divisor_context_(divisor, required_inverse_precision(divisor), arithmetic_context),
            x_to_field_order_(compute_x_to_field_order(divisor_context_, arithmetic_context)),
            composition_precomputation_(x_to_field_order_, std::max<std::size_t>(1, divisor_context_.degree()),
                                        divisor_context_, arithmetic_context) {
        }

        const polynomial_divisor_context<backend_type> &divisor_context() const {
            return divisor_context_;
        }

        const polynomial_type &x_to_field_order() const {
            return x_to_field_order_;
        }

        const polynomial_composition_precomputation<backend_type> &composition_precomputation() const {
            return composition_precomputation_;
        }

    private:
        static polynomial_type
            compute_x_to_field_order(const polynomial_divisor_context<backend_type> &divisor_context,
                                     polynomial_arithmetic::polynomial_context<backend_type> &arithmetic_context) {
            const polynomial_type x = {value_type {}, value_type::one()};
            polynomial_type result;
            powmod(result, x, algebra::fields::field_order<field_type>(), divisor_context, arithmetic_context);
            return result;
        }

        static std::size_t required_inverse_precision(const polynomial_type &divisor) {
            // Products of representatives of degree below d have quotients with at most d - 1 coefficients. Keep
            // one coefficient for constant and linear divisors because the divisor context requires positive
            // precision.
            return std::max<std::size_t>(1, divisor.size() > 1 ? divisor.size() - 2 : 0);
        }

        polynomial_divisor_context<backend_type> divisor_context_;
        polynomial_type x_to_field_order_;
        polynomial_composition_precomputation<backend_type> composition_precomputation_;
    };

    /**
     * Apply the Q-power Frobenius map to input in K[X]/(B), using the cached value X^Q mod B. The canonical result
     * may alias input.
     */
    template<detail::SupportsDivrem Backend>
    void frobenius_map(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input,
                       const polynomial_frobenius_context<Backend> &frobenius_context,
                       polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        if (input.size() <= frobenius_context.composition_precomputation().maximum_outer_coefficient_count()) {
            compose_mod(output, input, frobenius_context.composition_precomputation(),
                        frobenius_context.divisor_context(), arithmetic_context);
        } else {
            // Quotient-ring representatives use the reusable precomputation. Retain support for larger, unreduced
            // inputs by constructing a one-off composition precomputation sized for that input.
            compose_mod(output, input, frobenius_context.x_to_field_order(), frobenius_context.divisor_context(),
                        arithmetic_context);
        }
    }

    /**
     * Apply the Q-power Frobenius map iteration_count times. A zero iteration count returns the canonical residue of
     * input. The result may alias input.
     *
     * @throws std::invalid_argument if reducing input at iteration count zero requires more inverse coefficients than
     *         the context stores.
     */
    template<detail::SupportsDivrem Backend>
    void frobenius_map(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input,
                       std::size_t iteration_count, const polynomial_frobenius_context<Backend> &frobenius_context,
                       polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;

        if (iteration_count == 0) {
            remainder(output, input, frobenius_context.divisor_context(), arithmetic_context);
            return;
        }

        polynomial_type result;
        frobenius_map(result, input, frobenius_context, arithmetic_context);
        for (std::size_t iteration = 1; iteration < iteration_count; ++iteration) {
            frobenius_map(result, result, frobenius_context, arithmetic_context);
        }
        output = std::move(result);
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_FROBENIUS_HPP
