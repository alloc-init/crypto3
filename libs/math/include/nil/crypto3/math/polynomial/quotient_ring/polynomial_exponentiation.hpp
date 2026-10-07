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
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <nil/crypto3/math/polynomial/operations/shift.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_modular_arithmetic.hpp>

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

    /**
     * Compute output = X^exponent mod H, where H is the nonzero polynomial stored in divisor_context. H need not
     * be monic or irreducible. Left-to-right binary exponentiation uses L - 1 modular squarings and popcount(exponent)
     * - 1 coefficient shifts/reductions for a positive L-bit exponent and nonconstant H. Each multiply-by-X step
     * takes O(degree(H)) field operations, without general polynomial multiplication or division.
     *
     * Squarings use squaremod and its existing inverse-precision requirements: products of reduced operands need
     * at most degree(H) - 1 inverse coefficients. Exponent zero returns the quotient-ring identity, including [0]
     * modulo a nonzero constant. Output is canonical and is replaced only after the computation succeeds.
     *
     * @throws std::invalid_argument if exponent is negative or a squaring requires more inverse precision.
     */
    template<detail::SupportsDivrem Backend, detail::IntegerExponent Exponent>
    void powmod_x(typename Backend::polynomial_type &output, const Exponent &exponent,
                  const polynomial_divisor_context<Backend> &divisor_context,
                  polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        if constexpr (std::signed_integral<Exponent> || !std::integral<Exponent>) {
            if (exponent < 0) {
                throw std::invalid_argument("polynomial exponent must be nonnegative");
            }
        }

        const std::size_t divisor_degree = divisor_context.degree();
        if (divisor_degree == 0) {
            output.assign(1, value_type {});
            return;
        }
        if (exponent == 0) {
            output.assign(1, value_type::one());
            return;
        }

        // IntegerExponent only requires right shifts and parity. Record bits without narrowing the exponent or
        // requiring a backend-specific most-significant-bit operation, then consume them from the other end.
        std::vector<bool> bits;
        for (Exponent remaining(exponent); remaining != 0; remaining >>= 1) {
            bits.push_back((remaining % 2) != 0);
        }

        const polynomial_type &modulus = divisor_context.divisor();
        // The constant coefficient of rev(H)^-1 is already 1/lead(H); reuse it without another field inversion.
        const value_type &inverse_lead = divisor_context.reversed_divisor_inverse()[0];
        polynomial_type result;
        if (divisor_degree == 1) {
            result.assign(1, -modulus[0] * inverse_lead);
        } else {
            result.resize(2, value_type {});
            result[1] = value_type::one();
        }
        bits.pop_back();    // The leading set bit initializes result to X mod H.

        while (!bits.empty()) {
            // If result = X^k mod H for the processed prefix k, the next bit b = bits.back() gives
            // X^(2*k+b): square result, then multiply by X only when b is one.
            squaremod(result, result, divisor_context, arithmetic_context);
            if (bits.back()) {
                if (result.size() == divisor_degree) {
                    // With d = divisor_degree, X*result has degree d and leading coefficient result[d-1].
                    // Subtract c*H once, with c = result[d-1]/lead(H). Store its lower d coefficients directly,
                    // descending to preserve the old result[i-1]; the canceled X^d coefficient needs no storage.
                    const value_type c = result[divisor_degree - 1] * inverse_lead;
                    for (std::size_t i = divisor_degree - 1; i > 0; --i) {
                        result[i] = result[i - 1] - c * modulus[i];
                    }
                    result[0] = -c * modulus[0];
                    condense(result);
                } else {
                    // The shift stays below degree(H); this helper also preserves canonical zero.
                    shift_left(result, result, 1);
                }
            }
            bits.pop_back();
        }

        output = std::move(result);
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_EXPONENTIATION_HPP
