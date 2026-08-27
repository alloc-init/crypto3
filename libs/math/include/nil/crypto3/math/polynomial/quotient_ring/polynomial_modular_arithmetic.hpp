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

#ifndef CRYPTO3_MATH_POLYNOMIAL_MODULAR_ARITHMETIC_HPP
#define CRYPTO3_MATH_POLYNOMIAL_MODULAR_ARITHMETIC_HPP

#include <nil/crypto3/math/polynomial/arithmetic/polynomial_division.hpp>

namespace nil::crypto3::math {

    /**
     * Compute output = (left * right) mod B, where B is the nonzero polynomial stored in divisor_context. B need not
     * be monic or irreducible. The result is canonical, has degree less than degree(B) when B is nonconstant, and may
     * alias either input.
     *
     * If the canonical product has p coefficients and d = degree(B), reduction requires p - d quotient coefficients
     * when p > d. The context's precomputed inverse must have at least that precision. In particular, for nonconstant
     * B, inputs already reduced modulo B require at most d - 1 inverse coefficients.
     *
     * @throws std::invalid_argument if the precomputed inverse has insufficient precision.
     */
    template<detail::SupportsDivrem Backend>
    void mulmod(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &left,
                const typename Backend::polynomial_type &right,
                const polynomial_divisor_context<Backend> &divisor_context,
                polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        // Every polynomial is zero modulo a nonzero constant, so no product needs to be computed.
        if (divisor_context.degree() == 0) {
            output.resize(1);
            output[0] = value_type {};
            return;
        }

        polynomial_type product;
        multiplication(product, left, right, arithmetic_context);
        remainder(output, product, divisor_context, arithmetic_context);
    }

    /**
     * Compute output = input^2 mod B, where B is the nonzero polynomial stored in divisor_context. This is the
     * squaring counterpart of mulmod and uses the backend's dedicated square operation. The result is canonical, has
     * degree less than degree(B) when B is nonconstant, and may alias input.
     *
     * @throws std::invalid_argument if the precomputed inverse has insufficient precision.
     */
    template<detail::SupportsDivrem Backend>
    void squaremod(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input,
                   const polynomial_divisor_context<Backend> &divisor_context,
                   polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        if (divisor_context.degree() == 0) {
            output.resize(1);
            output[0] = value_type {};
            return;
        }

        polynomial_type square;
        arithmetic_context.square(square, input);
        remainder(output, square, divisor_context, arithmetic_context);
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_MODULAR_ARITHMETIC_HPP
