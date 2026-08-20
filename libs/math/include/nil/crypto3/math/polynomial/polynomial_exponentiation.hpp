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

#ifndef CRYPTO3_MATH_POLYNOMIAL_EXPONENTIATION_HPP
#define CRYPTO3_MATH_POLYNOMIAL_EXPONENTIATION_HPP

#include <concepts>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <nil/crypto3/math/polynomial/polynomial_modular_arithmetic.hpp>

namespace nil::crypto3::math {

    namespace detail {
        /** Integer-like exponent that can be consumed one binary digit at a time. */
        template<typename Exponent>
        concept IntegerExponent =
            std::copy_constructible<Exponent> && requires(Exponent value, const Exponent constant) {
                { constant == 0 } -> std::convertible_to<bool>;
                { constant < 0 } -> std::convertible_to<bool>;
                { (constant % 2) != 0 } -> std::convertible_to<bool>;
                { value >>= 1 } -> std::same_as<Exponent &>;
            };
    }    // namespace detail

    /**
     * Compute output = base^exponent mod B by binary exponentiation, where B is the nonzero polynomial stored in
     * divisor_context. B need not be monic or irreducible. The base is reduced before exponentiation, and every
     * intermediate square and product is reduced before the next operation. Output is canonical and may alias base.
     *
     * If d = degree(B), products of reduced operands require at most d - 1 inverse coefficients. Reducing an
     * unreduced base may require more, so divisor_context must also have sufficient precision for that initial
     * reduction. Exponent zero returns the quotient-ring identity; modulo a nonzero constant this is the zero
     * polynomial because the quotient ring is the zero ring.
     *
     * @throws std::invalid_argument if exponent is negative or the precomputed inverse has insufficient precision.
     */
    template<detail::SupportsDivrem Backend, detail::IntegerExponent Exponent>
    void powmod(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &base,
                const Exponent &exponent, const polynomial_divisor_context<Backend> &divisor_context,
                polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        if constexpr (std::signed_integral<Exponent> || !std::integral<Exponent>) {
            if (exponent < 0) {
                throw std::invalid_argument("polynomial exponent must be nonnegative");
            }
        }

        if (divisor_context.degree() == 0) {
            output.assign(1, value_type {});
            return;
        }
        if (exponent == 0) {
            output.assign(1, value_type::one());
            return;
        }

        polynomial_type power;
        remainder(power, base, divisor_context, arithmetic_context);

        Exponent remaining(exponent);
        polynomial_type result;
        bool found_set_bit = false;
        while (remaining != 0) {
            if ((remaining % 2) != 0) {
                if (found_set_bit) {
                    mulmod(result, result, power, divisor_context, arithmetic_context);
                } else {
                    result = power;
                    found_set_bit = true;
                }
            }

            remaining >>= 1;
            if (remaining != 0) {
                squaremod(power, power, divisor_context, arithmetic_context);
            }
        }

        output = std::move(result);
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_EXPONENTIATION_HPP
