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

#ifndef CRYPTO3_MATH_POLYNOMIAL_X_NORM_RECONSTRUCTION_HPP
#define CRYPTO3_MATH_POLYNOMIAL_X_NORM_RECONSTRUCTION_HPP

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/fields/field_algorithms.hpp>

#include <nil/crypto3/math/polynomial/operations/basic_operations.hpp>
#include <nil/crypto3/math/polynomial/operations/shift.hpp>
#include <nil/crypto3/math/polynomial/factorization/complete_factorization.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_square_root.hpp>
#include <nil/crypto3/math/polynomial/reconstruction/polynomial_rational_reconstruction.hpp>

namespace nil::crypto3::math {

    /** Coefficients of P + Q * sqrt(X), represented by the polynomials P and Q. */
    template<CoefficientPolynomial Polynomial>
    struct polynomial_x_norm_representation {
        Polynomial p;
        Polynomial q;
    };

    namespace detail {

        template<CoefficientPolynomial Polynomial>
        bool is_canonical_polynomial_x_norm_representation(
            const polynomial_x_norm_representation<Polynomial> &representation) {
            using value_type = typename Polynomial::value_type;
            const auto is_canonical = [](const Polynomial &polynomial) {
                return !polynomial.empty() &&
                       (polynomial.size() == 1 || polynomial[polynomial.size() - 1] != value_type {});
            };
            return is_canonical(representation.p) && is_canonical(representation.q);
        }

    }    // namespace detail

    /**
     * Evaluate the polynomial norm of P + Q * sqrt(X):
     *
     *     (P + Q * sqrt(X)) * (P - Q * sqrt(X)) = P^2 - X * Q^2.
     *
     * Both squares use the caller-owned arithmetic context. Multiplication by X is a coefficient shift and therefore
     * does not require a polynomial product.
     *
     * @throws std::invalid_argument if P or Q is empty or noncanonical.
     */
    template<polynomial_arithmetic::PolynomialBackend Backend>
    typename Backend::polynomial_type evaluate_polynomial_x_norm(
        const polynomial_x_norm_representation<typename Backend::polynomial_type> &representation,
        polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        if (!detail::is_canonical_polynomial_x_norm_representation(representation)) {
            throw std::invalid_argument("polynomial X-norm evaluation requires canonical nonempty inputs");
        }

        polynomial_type p_squared;
        polynomial_type q_squared;
        polynomial_type x_q_squared;
        polynomial_type norm;
        arithmetic_context.square(p_squared, representation.p);
        arithmetic_context.square(q_squared, representation.q);
        shift_left(x_q_squared, q_squared, 1);
        subtraction(norm, p_squared, x_q_squared);
        return norm;
    }

    /**
     * Multiply two polynomial X-norm representations using
     *
     *     P = P1 * P2 + X * Q1 * Q2,
     *     Q = P1 * Q2 + Q1 * P2.
     *
     * The returned representation has norm equal to the product of the input norms. Every polynomial product uses the
     * caller-owned arithmetic context; multiplication by X is a coefficient shift.
     *
     * @throws std::invalid_argument if an input polynomial is empty or noncanonical.
     */
    template<polynomial_arithmetic::PolynomialBackend Backend>
    polynomial_x_norm_representation<typename Backend::polynomial_type> multiply_polynomial_x_norm_representations(
        const polynomial_x_norm_representation<typename Backend::polynomial_type> &left,
        const polynomial_x_norm_representation<typename Backend::polynomial_type> &right,
        polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using representation_type = polynomial_x_norm_representation<polynomial_type>;

        if (!detail::is_canonical_polynomial_x_norm_representation(left) ||
            !detail::is_canonical_polynomial_x_norm_representation(right)) {
            throw std::invalid_argument("polynomial X-norm multiplication requires canonical nonempty inputs");
        }

        polynomial_type p_product;
        polynomial_type q_product;
        polynomial_type shifted_q_product;
        polynomial_type result_p;
        arithmetic_context.multiply(p_product, left.p, right.p);
        arithmetic_context.multiply(q_product, left.q, right.q);
        shift_left(shifted_q_product, q_product, 1);
        addition(result_p, p_product, shifted_q_product);

        polynomial_type left_p_right_q;
        polynomial_type left_q_right_p;
        polynomial_type result_q;
        arithmetic_context.multiply(left_p_right_q, left.p, right.q);
        arithmetic_context.multiply(left_q_right_p, left.q, right.p);
        addition(result_q, left_p_right_q, left_q_right_p);

        return representation_type {std::move(result_p), std::move(result_q)};
    }

    /**
     * Recover P and Q satisfying
     *
     *     P^2 - X * Q^2 = g,
     *
     * with degree(P) at most floor(degree(g) / 2) and degree(Q) at most
     * floor((degree(g) - 1) / 2). The coefficient generator supplies field elements used to find a quadratic
     * nonresidue in K[X]/(g); it remains owned by the caller.
     *
     * Irreducibility of g is a caller precondition and is not tested. The coefficient generator must eventually produce
     * coefficients forming a nonsquare canonical representative of degree below degree(g).
     *
     * @return a normalized representation whose evaluated norm is exactly g; no value if X is nonsquare modulo g,
     *         bounded rational reconstruction fails, or the resulting nonzero scalar multiple of g cannot be
     * normalized.
     * @throws std::invalid_argument if g is empty, noncanonical, zero, or constant, or another documented precondition
     *         of a composed polynomial operation is violated.
     * @throws std::logic_error if a composed operation reports success but its resulting identities are inconsistent.
     */
    template<detail::SupportsDivrem Backend, typename Generator>
        requires algebra::FieldValue<typename Backend::polynomial_type::value_type> &&
                 std::constructible_from<typename Backend::polynomial_type, std::size_t> &&
                 requires(typename Backend::polynomial_type &polynomial, Generator &generator,
                          const typename Backend::polynomial_type::value_type &value) {
                     polynomial[0] = generator();
                     { value.is_square() } -> std::convertible_to<bool>;
                 }
    std::optional<polynomial_x_norm_representation<typename Backend::polynomial_type>>
        recover_irreducible_polynomial_x_norm_representation(
            const typename Backend::polynomial_type &g,
            polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
            Generator &coefficient_generator) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;
        using representation_type = polynomial_x_norm_representation<polynomial_type>;

        if (g.empty() || (g.size() > 1 && g[g.size() - 1] == value_type {})) {
            throw std::invalid_argument("polynomial X-norm recovery requires a canonical nonempty polynomial");
        }
        if (g.size() == 1) {
            throw std::invalid_argument("polynomial X-norm recovery requires a nonconstant polynomial");
        }

        const std::size_t degree = g.size() - 1;
        const std::size_t inverse_precision = std::max<std::size_t>(1, degree - 1);
        polynomial_divisor_context<Backend> divisor_context(g, inverse_precision, arithmetic_context);

        // Work with the canonical representative of X in K[X]/(g). Reduction is necessary when g is linear.
        const polynomial_type x = {value_type::zero(), value_type::one()};
        polynomial_type x_mod_g;
        remainder(x_mod_g, x, divisor_context, arithmetic_context);
        if (!is_square_mod(x_mod_g, divisor_context, arithmetic_context)) {
            return std::nullopt;
        }

        // Adapt caller-generated field coefficients into canonical representatives of K[X]/(g). No random source or
        // polynomial backend is constructed here.
        auto quotient_representative_generator = [&]() {
            polynomial_type representative(degree);
            for (std::size_t index = 0; index < degree; ++index) {
                representative[index] = coefficient_generator();
            }
            condense(representative);
            return representative;
        };
        polynomial_square_root_context<Backend> square_root_context(divisor_context, arithmetic_context,
                                                                    quotient_representative_generator);

        // Recover R with R^2 = X mod g. A failure here contradicts the successful square test above.
        polynomial_type root;
        if (!square_root_mod(root, x_mod_g, square_root_context, arithmetic_context)) {
            throw std::logic_error("polynomial X-norm recovery failed after X was reported square modulo g");
        }

        // Reconstruct P = R * Q mod g with the unique degree bounds required by the norm equation.
        polynomial_type p;
        polynomial_type q;
        if (!rational_reconstruct(p, q, root, g, degree / 2, (degree - 1) / 2, arithmetic_context)) {
            return std::nullopt;
        }

        // The two modular identities imply P^2 - X * Q^2 = lambda * g. The degree bounds ensure that lambda is a
        // scalar. A nonzero remainder or nonconstant quotient would contradict those successful operations.
        representation_type representation {std::move(p), std::move(q)};
        const polynomial_type norm = evaluate_polynomial_x_norm(representation, arithmetic_context);
        polynomial_type scalar_quotient;
        polynomial_type scalar_remainder;
        divrem(scalar_quotient, scalar_remainder, norm, divisor_context, arithmetic_context);
        if (!is_zero(scalar_remainder) || scalar_quotient.size() != 1) {
            throw std::logic_error("polynomial X-norm reconstruction produced an inconsistent scalar multiple");
        }

        const value_type lambda = scalar_quotient[0];
        if (lambda.is_zero() || !lambda.is_square()) {
            return std::nullopt;
        }

        // Multiplying P and Q by sqrt(lambda^-1) changes their norm from lambda * g to exactly g.
        const value_type normalization = algebra::fields::sqrt_known_square(lambda.inversed());
        scalar_multiplication(representation.p, representation.p, normalization);
        scalar_multiplication(representation.q, representation.q, normalization);

        if (evaluate_polynomial_x_norm(representation, arithmetic_context) != g) {
            throw std::logic_error("normalized polynomial X-norm representation failed exact verification");
        }
        return std::move(representation);
    }

    namespace detail {

        /** Raise a canonical polynomial to a nonnegative integer power using the supplied arithmetic context. */
        template<polynomial_arithmetic::PolynomialBackend Backend>
        typename Backend::polynomial_type
            polynomial_x_norm_power(const typename Backend::polynomial_type &base, std::size_t exponent,
                                    polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
            using polynomial_type = typename Backend::polynomial_type;
            using value_type = typename polynomial_type::value_type;

            polynomial_type result = {value_type::one()};
            polynomial_type current_power(base);
            while (exponent != 0) {
                if ((exponent & 1) != 0) {
                    polynomial_type product;
                    arithmetic_context.multiply(product, result, current_power);
                    result = std::move(product);
                }
                exponent >>= 1;
                if (exponent != 0) {
                    polynomial_type square;
                    arithmetic_context.square(square, current_power);
                    current_power = std::move(square);
                }
            }
            return result;
        }

        /**
         * Construct an X-norm representation of one irreducible factor raised to its multiplicity. Complete
         * factorization expresses the input as a leading scalar times a product of powers g^e, so each such power
         * needs a representation before the factor representations can be combined.
         *
         * If e = 2r, then g^e is already a square. The pair (g^r, 0) represents it because
         *
         *     (g^r)^2 - X * 0^2 = g^(2r).
         *
         * If e = 2r + 1, recover (P_g, Q_g) for the one unpaired copy of g, where
         *
         *     P_g^2 - X * Q_g^2 = g.
         *
         * Scaling both components by g^r gives a representation of the complete factor power:
         *
         *     (g^r * P_g)^2 - X * (g^r * Q_g)^2 = g^(2r) * g = g^e.
         *
         * For example, if H contains g^3 * h^2 and (P_g, Q_g) represents g, then (g * P_g, g * Q_g)
         * represents g^3, while (h, 0) represents h^2. Combining those two representations produces the
         * representation of g^3 * h^2. Thus only the odd-multiplicity factor g requires irreducible recovery.
         */
        template<SupportsDivrem Backend, typename Generator>
        std::optional<polynomial_x_norm_representation<typename Backend::polynomial_type>>
            recover_polynomial_x_norm_factor_power(
                const polynomial_factor<typename Backend::polynomial_type> &factor,
                polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                Generator &coefficient_generator) {
            using polynomial_type = typename Backend::polynomial_type;
            using value_type = typename polynomial_type::value_type;
            using representation_type = polynomial_x_norm_representation<polynomial_type>;

            if (factor.multiplicity == 0) {
                throw std::logic_error("complete factorization produced a factor with zero multiplicity");
            }

            const std::size_t half_multiplicity = factor.multiplicity / 2;
            polynomial_type half_power =
                polynomial_x_norm_power<Backend>(factor.polynomial, half_multiplicity, arithmetic_context);
            if ((factor.multiplicity & 1) == 0) {
                return representation_type {std::move(half_power), polynomial_type {value_type::zero()}};
            }

            auto odd_representation = recover_irreducible_polynomial_x_norm_representation<Backend>(
                factor.polynomial, arithmetic_context, coefficient_generator);
            if (!odd_representation) {
                return std::nullopt;
            }
            if (half_multiplicity == 0) {
                return odd_representation;
            }

            polynomial_type lifted_p;
            polynomial_type lifted_q;
            arithmetic_context.multiply(lifted_p, half_power, odd_representation->p);
            arithmetic_context.multiply(lifted_q, half_power, odd_representation->q);
            return representation_type {std::move(lifted_p), std::move(lifted_q)};
        }

        /**
         * Combine X-norm representations in balanced levels so no left-deep product chain is formed. For example,
         * five factor representations are combined as
         *
         *     [A, B, C, D, E]
         *     [A * B, C * D, E]
         *     [(A * B) * (C * D), E]
         *     [(A * B) * (C * D) * E].
         *
         * Each product uses the X-norm product identity. If a level has an odd number of representations, its final
         * representation is carried unchanged to the next level. An empty input represents the empty product and
         * therefore returns the multiplicative identity (1, 0).
         */
        template<polynomial_arithmetic::PolynomialBackend Backend>
        polynomial_x_norm_representation<typename Backend::polynomial_type>
            combine_polynomial_x_norm_representations_balanced(
                std::vector<polynomial_x_norm_representation<typename Backend::polynomial_type>>
                    representations,
                polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
            using polynomial_type = typename Backend::polynomial_type;
            using value_type = typename polynomial_type::value_type;
            using representation_type = polynomial_x_norm_representation<polynomial_type>;

            if (representations.empty()) {
                return representation_type {polynomial_type {value_type::one()}, polynomial_type {value_type::zero()}};
            }

            while (representations.size() > 1) {
                std::vector<representation_type> next_level;
                next_level.reserve((representations.size() + 1) / 2);
                std::size_t index = 0;
                for (; index + 1 < representations.size(); index += 2) {
                    next_level.push_back(multiply_polynomial_x_norm_representations<Backend>(
                        representations[index], representations[index + 1], arithmetic_context));
                }
                if (index < representations.size()) {
                    next_level.push_back(std::move(representations[index]));
                }
                representations = std::move(next_level);
            }
            return std::move(representations.front());
        }

    }    // namespace detail

    /**
     * Recover P and Q satisfying
     *
     *     P^2 - X * Q^2 = h
     *
     * for a canonical polynomial h. Zero and square constants are handled directly. A nonconstant input is factored
     * into monic irreducible factors. Even factor multiplicities are represented as polynomial squares; odd
     * multiplicities use recover_irreducible_polynomial_x_norm_representation. Factor representations are combined in
     * a balanced product tree, then scaled by the square root of the factorization's leading coefficient.
     *
     * The coefficient generator remains caller-owned and is shared by complete factorization and irreducible-factor
     * recovery. It must satisfy the documented requirements of both operations.
     *
     * @return a representation whose evaluated norm is exactly h; no value if a necessary coefficient square test,
     *         odd-factor recovery, or leading-scalar normalization fails.
     * @throws std::invalid_argument if h is empty or noncanonical, or a composed factorization or recovery contract is
     *         violated.
     * @throws std::logic_error if completed internal operations produce an inconsistent identity.
     */
    template<detail::SupportsDivrem Backend, typename Generator>
        requires algebra::FieldValue<typename Backend::polynomial_type::value_type> &&
                 std::constructible_from<typename Backend::polynomial_type, std::size_t> &&
                 requires(typename Backend::polynomial_type &polynomial, Generator &generator,
                          const typename Backend::polynomial_type::value_type &value) {
                     polynomial[0] = generator();
                     { value.is_square() } -> std::convertible_to<bool>;
                 }
    std::optional<polynomial_x_norm_representation<typename Backend::polynomial_type>>
        recover_polynomial_x_norm_representation(const typename Backend::polynomial_type &h,
                                                 polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                                 Generator &coefficient_generator) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;
        using representation_type = polynomial_x_norm_representation<polynomial_type>;

        if (h.empty() || (h.size() > 1 && h[h.size() - 1] == value_type::zero())) {
            throw std::invalid_argument("polynomial X-norm recovery requires a canonical nonempty polynomial");
        }

        const auto verify_exact = [&](representation_type representation) -> std::optional<representation_type> {
            if (evaluate_polynomial_x_norm<Backend>(representation, arithmetic_context) != h) {
                throw std::logic_error("polynomial X-norm recovery failed exact verification");
            }
            return std::move(representation);
        };

        if (is_zero(h)) {
            return verify_exact(
                representation_type {polynomial_type {value_type::zero()}, polynomial_type {value_type::zero()}});
        }
        if (h.size() == 1) {
            if (!h[0].is_square()) {
                return std::nullopt;
            }
            return verify_exact(representation_type {polynomial_type {algebra::fields::sqrt_known_square(h[0])},
                                                     polynomial_type {value_type::zero()}});
        }

        // square filters
        const std::size_t degree = h.size() - 1;
        if (!h[0].is_square()) {
            return std::nullopt;
        }
        value_type signed_leading_coefficient = h[h.size() - 1];
        if ((degree & 1) != 0) {
            signed_leading_coefficient = value_type::zero() - signed_leading_coefficient;
        }
        if (!signed_leading_coefficient.is_square()) {
            return std::nullopt;
        }

        std::vector<representation_type> factor_representations;
        bool factor_recovery_failed = false;
        const auto factorization = complete_factorization<Backend>(
            h, arithmetic_context, coefficient_generator, [&](const polynomial_factor<polynomial_type> &factor) {
                auto representation = detail::recover_polynomial_x_norm_factor_power<Backend>(
                    factor, arithmetic_context, coefficient_generator);
                if (!representation) {
                    factor_recovery_failed = true;
                    return factorization_control::stop_factorization;
                }
                factor_representations.push_back(std::move(*representation));
                return factorization_control::continue_factorization;
            });

        if (factor_recovery_failed) {
            return std::nullopt;
        }
        if (!factorization.complete || factor_representations.empty()) {
            throw std::logic_error("complete factorization did not produce all nonconstant factors");
        }

        representation_type result = detail::combine_polynomial_x_norm_representations_balanced<Backend>(
            std::move(factor_representations), arithmetic_context);
        const value_type leading_coefficient = factorization.leading_coefficient;
        if (leading_coefficient.is_zero()) {
            throw std::logic_error("complete factorization produced a zero leading coefficient");
        }
        if (!leading_coefficient.is_square()) {
            return std::nullopt;
        }
        const value_type scalar = algebra::fields::sqrt_known_square(leading_coefficient);
        scalar_multiplication(result.p, result.p, scalar);
        scalar_multiplication(result.q, result.q, scalar);
        return verify_exact(std::move(result));
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_X_NORM_RECONSTRUCTION_HPP
