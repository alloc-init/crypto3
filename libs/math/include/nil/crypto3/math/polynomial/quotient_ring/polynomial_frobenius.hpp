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
#include <optional>
#include <stdexcept>
#include <utility>

#include <nil/crypto3/algebra/fields/field_order.hpp>

#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_composition.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_exponentiation.hpp>

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
            polynomial_type result;
            powmod_x(result, algebra::fields::field_order<field_type>(), divisor_context, arithmetic_context);
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

    namespace detail {

        // Shared chain for coefficient-fixed Frobenius products and the characteristic-power square-root seed.
        // In the characteristic case the supplied image is X^p, and coefficients must support Frobenius_map.
        template<bool CharacteristicFrobenius, SupportsDivrem Backend>
        void frobenius_product_impl(typename Backend::polynomial_type &output,
                                    const typename Backend::polynomial_type &input, std::size_t term_count,
                                    const typename Backend::polynomial_type &x_to_frobenius_power,
                                    const polynomial_divisor_context<Backend> &divisor_context,
                                    polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
            using polynomial_type = typename Backend::polynomial_type;
            using value_type = typename polynomial_type::value_type;
            using cache_type = polynomial_composition_precomputation<Backend>;
            const std::size_t degree = divisor_context.degree();
            const auto canonical = [](const polynomial_type &value) {
                return !value.empty() && (value.size() == 1 || value.back() != value_type::zero());
            };
            if (!canonical(input) || !canonical(x_to_frobenius_power) ||
                x_to_frobenius_power.size() > std::max<std::size_t>(1, degree)) {
                throw std::invalid_argument(
                    "Frobenius products require canonical inputs and a reduced Frobenius image");
            }
            if (degree == 0 || term_count == 0) {
                output.assign(1, degree == 0 ? value_type::zero() : value_type::one());
                return;
            }
            polynomial_type base;
            remainder(base, input, divisor_context, arithmetic_context);
            if (term_count == 1) {
                output = std::move(base);
                return;
            }
            if constexpr (!CharacteristicFrobenius) {
                // Constants are fixed by Q-Frobenius, but generally not by p-Frobenius in an extension field.
                if (base.size() == 1) {
                    powmod(output, base, term_count, divisor_context, arithmetic_context);
                    return;
                }
            }

            const auto apply = [&](polynomial_type &result, const polynomial_type &operand, std::size_t power,
                                   const cache_type &cache) {
                if constexpr (CharacteristicFrobenius) {
                    // Absolute Frobenius is semilinear: f^(p^power) = sigma^power(f)(X^(p^power)).
                    // Copy before twisting so result may alias operand, including the cached-image input.
                    polynomial_type twisted = operand;
                    for (auto &coefficient : twisted) {
                        coefficient = coefficient.Frobenius_map(power);
                    }
                    compose_mod(result, twisted, cache, divisor_context, arithmetic_context);
                } else {
                    compose_mod(result, operand, cache, divisor_context, arithmetic_context);
                }
            };

            const cache_type first_cache(x_to_frobenius_power, degree, divisor_context, arithmetic_context);
            const auto chain = [&](auto &&self, std::size_t count,
                                   bool need_image) -> std::pair<polynomial_type, polynomial_type> {
                if (count == 1) {
                    return {base, x_to_frobenius_power};
                }
                auto [product, image] = self(self, count / 2, true);
                std::optional<cache_type> doubling_cache;
                if (count / 2 > 1) {
                    doubling_cache.emplace(image, degree, divisor_context, arithmetic_context);
                }
                const auto &cache = doubling_cache ? *doubling_cache : first_cache;
                // With P_k = product and H_k = image, P_2k = P_k * F^k(P_k) and H_2k = F^k(H_k).
                polynomial_type shifted;
                apply(shifted, product, count / 2, cache);
                mulmod(product, product, shifted, divisor_context, arithmetic_context);
                if (need_image) {
                    apply(image, image, count / 2, cache);
                }
                if (count % 2 != 0) {
                    // P_(2k+1) = base * F(P_2k); use the same initial cache at every odd step.
                    apply(shifted, product, 1, first_cache);
                    mulmod(product, base, shifted, divisor_context, arithmetic_context);
                    if (need_image) {
                        apply(image, image, 1, first_cache);
                    }
                }
                return {std::move(product), std::move(image)};
            };
            // The final image is not used by the caller, so avoid its compositions.
            auto result = chain(chain, term_count, false);
            output = std::move(result.first);
        }

    }    // namespace detail

    /**
     * Compute input * F(input) * ... * F^(term_count - 1)(input) modulo B, where F is a power of the coefficient
     * field's Frobenius map. The supplied x_to_frobenius_power is the canonical reduced image F(X); it can represent
     * X^(Q^s) for any s >= 0, with Q the coefficient-field order. Its algebraic identity is a caller precondition,
     * not checked by repeating the exponentiation. B need not be irreducible.
     *
     * Binary doubling uses O(log(term_count)) compositions and products. Each doubling shares one composition
     * precomputation between the partial product and the next Frobenius image. Odd steps reuse the initial image's
     * precomputation. Only the initial cache and one doubling cache are live at once.
     *
     * An empty product is the quotient-ring identity, and a singleton is the reduced input. The result is canonical
     * and may alias either input polynomial. Reducing an unreduced input may require more inverse coefficients than
     * products of reduced operands; otherwise degree(B) - 1 inverse coefficients suffice.
     *
     * @throws std::invalid_argument if an input is empty or noncanonical, the supplied image is not reduced, a
     *         needed composition has a zero cached-power limit, or the divisor inverse has insufficient precision.
     */
    template<detail::SupportsDivrem Backend>
    void frobenius_product(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input,
                           std::size_t term_count, const typename Backend::polynomial_type &x_to_frobenius_power,
                           const polynomial_divisor_context<Backend> &divisor_context,
                           polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        detail::frobenius_product_impl<false>(output, input, term_count, x_to_frobenius_power, divisor_context,
                                              arithmetic_context);
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_FROBENIUS_HPP
