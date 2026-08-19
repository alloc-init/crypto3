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

#ifndef CRYPTO3_MATH_POLYNOMIAL_COMPOSITION_HPP
#define CRYPTO3_MATH_POLYNOMIAL_COMPOSITION_HPP

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/math/detail/integer_sqrt.hpp>
#include <nil/crypto3/math/polynomial/polynomial_modular_arithmetic.hpp>

namespace nil::crypto3::math {

    /**
     * Compute output = outer(inner(X)) mod B by Horner's rule, where B is the nonzero polynomial stored in
     * divisor_context. This Horner algorithm performs one modular multiplication for each coefficient of outer except
     * its leading coefficient and is intended as a correctness oracle for faster modular composition algorithms. The
     * result is canonical, has degree less than degree(B) when B is nonconstant, and may alias either input.
     *
     * The inner polynomial is reduced before composition, and every Horner product is reduced before the next
     * coefficient is added. For d = degree(B), the products of reduced polynomials require at most d - 1 inverse
     * coefficients. Reducing an unreduced inner polynomial may require greater precomputed precision. Modulo a
     * nonzero constant, every composition is the zero polynomial.
     *
     * @throws std::invalid_argument if the precomputed inverse has insufficient precision.
     * @pre outer and inner are nonempty canonical coefficient polynomials.
     */
    template<detail::SupportsDivrem Backend>
    void compose_mod_reference(typename Backend::polynomial_type &output,
                               const typename Backend::polynomial_type &outer,
                               const typename Backend::polynomial_type &inner,
                               const polynomial_divisor_context<Backend> &divisor_context,
                               polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        if (divisor_context.degree() == 0) {
            output.assign(1, value_type {});
            return;
        }
        if (outer.size() == 1) {
            output = outer;
            return;
        }

        // Reference an already reduced inner polynomial directly to avoid copying it; use local storage only when
        // reduction is required.
        polynomial_type reduced_inner_storage;
        const polynomial_type *reduced_inner = &inner;
        if (inner.size() >= divisor_context.divisor().size()) {
            remainder(reduced_inner_storage, inner, divisor_context, arithmetic_context);
            reduced_inner = &reduced_inner_storage;
        }

        polynomial_type result {outer.back()};
        for (std::size_t i = outer.size() - 1; i-- > 0;) {
            mulmod(result, result, *reduced_inner, divisor_context, arithmetic_context);
            result[0] = result[0] + outer[i];
        }

        output = std::move(result);
    }

    /**
     * Compute output = outer(inner(X)) mod B using blocked Brent-Kung composition, where B is the nonzero polynomial
     * stored in divisor_context. The result is canonical, has degree less than degree(B) when B is nonconstant, and
     * may alias either input.
     *
     * For a block size k, split outer into polynomials F_j of degree less than k:
     *
     *     outer(Y) = F_0(Y) + F_1(Y) * Y^k + F_2(Y) * Y^(2k) + ...
     *
     * Cache the baby steps 1, inner, ..., inner^(k-1) modulo B. Each F_j(inner) is then a linear combination of these
     * cached powers. Set the giant step to inner^k modulo B and combine the block values by Horner's rule:
     *
     *     (...(F_last(inner) * giant_step + F_previous(inner)) * giant_step + ...) + F_0(inner).
     *
     * If outer has L coefficients, this uses about k + L/k modular multiplications, minimized by
     * k = ceil(sqrt(L)). Forming the block values also requires O(L * degree(B)) coefficient operations; this
     * quadratic work can dominate when L and degree(B) are both large. The cached-power limit may select a smaller k,
     * storing at most k * degree(B) coefficients at the cost of more giant-step multiplications.
     *
     * The inner polynomial is reduced before its powers are constructed. Reducing an unreduced inner polynomial may
     * require more precomputed inverse coefficients than products of reduced operands.
     *
     * @throws std::invalid_argument if the cached-power limit is zero or the precomputed inverse has insufficient
     *         precision.
     * @pre outer and inner are nonempty canonical coefficient polynomials.
     */
    template<detail::SupportsDivrem Backend>
    void compose_mod(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &outer,
                     const typename Backend::polynomial_type &inner,
                     const polynomial_divisor_context<Backend> &divisor_context,
                     polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        const std::size_t cached_power_limit = arithmetic_context.options().modular_composition_cached_power_limit;
        if (cached_power_limit == 0) {
            throw std::invalid_argument("the modular-composition cached-power limit must be positive");
        }

        if (divisor_context.degree() == 0) {
            output.assign(1, value_type {});
            return;
        }

        if (outer.size() == 1) {
            output = outer;
            return;
        }

        // Reference an already reduced inner polynomial directly to avoid copying it; use local storage only when
        // reduction is required.
        polynomial_type reduced_inner_storage;
        const polynomial_type *reduced_inner = &inner;
        if (inner.size() >= divisor_context.divisor().size()) {
            remainder(reduced_inner_storage, inner, divisor_context, arithmetic_context);
            reduced_inner = &reduced_inner_storage;
        }

        const std::size_t block_size = std::min(detail::ceil_sqrt(outer.size()), cached_power_limit);
        // The zeroth power contributes only to a block's constant coefficient, and the first power can reference the
        // reduced inner polynomial directly. Store only the computed powers from inner^2 onward.
        std::vector<polynomial_type> computed_baby_steps;
        computed_baby_steps.reserve(block_size > 2 ? block_size - 2 : 0);
        auto baby_step = [&](std::size_t exponent) -> const polynomial_type & {
            return exponent == 1 ? *reduced_inner : computed_baby_steps[exponent - 2];
        };
        // Build each power from two already cached powers with approximately equal exponents. Even exponents use the
        // dedicated squaring path; odd exponents multiply the floor and ceiling half-powers.
        for (std::size_t exponent = 2; exponent < block_size; ++exponent) {
            polynomial_type power;
            if (exponent % 2 == 0) {
                squaremod(power, baby_step(exponent / 2), divisor_context, arithmetic_context);
            } else {
                mulmod(power, baby_step(exponent / 2), baby_step(exponent / 2 + 1), divisor_context,
                       arithmetic_context);
            }
            computed_baby_steps.emplace_back(std::move(power));
        }

        polynomial_type giant_step_storage;
        const polynomial_type *giant_step = reduced_inner;
        if (block_size > 1) {
            if (block_size % 2 == 0) {
                squaremod(giant_step_storage, baby_step(block_size / 2), divisor_context, arithmetic_context);
            } else {
                mulmod(giant_step_storage, baby_step(block_size / 2), baby_step(block_size / 2 + 1), divisor_context,
                       arithmetic_context);
            }
            giant_step = &giant_step_storage;
        }

        // Evaluate F_block(inner) as a linear combination of the baby steps and add it directly to the Horner
        // accumulator. Streaming one block at a time avoids storing the full matrix of evaluated block polynomials.
        auto add_block = [&](polynomial_type &result, std::size_t block_index) {
            const std::size_t first_coefficient = block_index * block_size;
            const std::size_t block_coefficient_count = std::min(block_size, outer.size() - first_coefficient);
            for (std::size_t exponent = 0; exponent < block_coefficient_count; ++exponent) {
                const value_type &scalar = outer[first_coefficient + exponent];
                if (scalar == value_type {}) {
                    continue;
                }
                if (exponent == 0) {
                    result[0] = result[0] + scalar;
                    continue;
                }

                const polynomial_type &power = baby_step(exponent);
                if (result.size() < power.size()) {
                    result.resize(power.size(), value_type {});
                }
                for (std::size_t coefficient = 0; coefficient < power.size(); ++coefficient) {
                    result[coefficient] = result[coefficient] + power[coefficient] * scalar;
                }
            }
            condense(result);
        };

        const std::size_t block_count = 1 + (outer.size() - 1) / block_size;
        polynomial_type result;
        add_block(result, block_count - 1);
        for (std::size_t block = block_count - 1; block-- > 0;) {
            mulmod(result, result, *giant_step, divisor_context, arithmetic_context);
            add_block(result, block);
        }

        output = std::move(result);
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_COMPOSITION_HPP
