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

#include <nil/crypto3/algebra/fields/field_algorithms.hpp>

#include <nil/crypto3/math/polynomial/operations/basic_operations.hpp>
#include <nil/crypto3/math/polynomial/operations/shift.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_square_root.hpp>
#include <nil/crypto3/math/polynomial/reconstruction/polynomial_rational_reconstruction.hpp>

namespace nil::crypto3::math {

    /** Coefficients of P + Q * sqrt(X), represented by the polynomials P and Q. */
    template<CoefficientPolynomial Polynomial>
    struct polynomial_x_norm_representation {
        Polynomial p;
        Polynomial q;
    };

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
        using value_type = typename polynomial_type::value_type;

        const auto is_canonical = [](const polynomial_type &polynomial) {
            return !polynomial.empty() &&
                   (polynomial.size() == 1 || polynomial[polynomial.size() - 1] != value_type {});
        };
        if (!is_canonical(representation.p) || !is_canonical(representation.q)) {
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

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_X_NORM_RECONSTRUCTION_HPP
