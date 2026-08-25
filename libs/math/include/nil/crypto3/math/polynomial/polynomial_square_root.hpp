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

#ifndef CRYPTO3_MATH_POLYNOMIAL_SQUARE_ROOT_HPP
#define CRYPTO3_MATH_POLYNOMIAL_SQUARE_ROOT_HPP

#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <utility>

#include <boost/multiprecision/cpp_int.hpp>

#include <nil/crypto3/algebra/fields/field_order.hpp>
#include <nil/crypto3/algebra/type_traits.hpp>

#include <nil/crypto3/math/polynomial/polynomial_exponentiation.hpp>

namespace nil::crypto3::math {

    namespace detail {

        template<SupportsDivrem Backend>
            requires algebra::FieldValue<typename Backend::polynomial_type::value_type>
        bool is_square_mod_impl(const typename Backend::polynomial_type &input,
                                const polynomial_divisor_context<Backend> &divisor_context,
                                polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                const algebra::fields::multiplicative_group_decomposition *order_decomposition);

    }    // namespace detail

    /**
     * Immutable cache of multiplicative-order parameters for square roots in the quotient field K[X]/(B). If K has Q
     * elements and d = degree(B), the shared field-order utility decomposes
     *
     *     Q^d - 1 = odd_order * 2^two_adicity.
     *
     * Given a quadratic nonresidue z in that quotient field, the context also stores z^odd_order. Tonelli-Shanks reuses
     * the order decomposition and this cached nonresidue power for every square root modulo B. Irreducibility is a
     * caller precondition and is not tested. The referenced divisor context must outlive this context.
     *
     * Tonelli-Shanks needs a nonsquare z. Whether a polynomial is a square depends on B, so the second constructor asks
     * a caller-owned generator for candidates until one is nonsquare modulo B. Each candidate polynomial represents
     * one element of K[X]/(B). The caller owns the generator and controls its random source and seed. The first
     * constructor can be used when the caller already knows a suitable z.
     *
     * @throws std::invalid_argument if B is constant, K has characteristic two, the explicitly supplied nonresidue is
     *         square, or a supplied or generated representative is empty, noncanonical, or not reduced modulo B.
     * @pre B is irreducible. An explicitly supplied nonresidue is canonical. A generator returns canonical reduced
     *      representatives and eventually returns a nonsquare.
     */
    template<detail::SupportsDivrem Backend>
        requires algebra::FieldValue<typename Backend::polynomial_type::value_type>
    class polynomial_square_root_context {
    public:
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;
        using field_type = typename value_type::field_type;

        polynomial_square_root_context(const polynomial_type &quadratic_non_residue,
                                       const polynomial_divisor_context<Backend> &divisor_context,
                                       polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) :
            divisor_context_(&divisor_context) {
            initialize_order_decomposition();
            require_canonical_reduced(quadratic_non_residue);
            if (detail::is_square_mod_impl(quadratic_non_residue, divisor_context, arithmetic_context,
                                           &order_decomposition_)) {
                throw std::invalid_argument("a polynomial square-root context requires a quadratic nonresidue");
            }
            cache_non_residue_power(quadratic_non_residue, arithmetic_context);
        }

        template<typename Generator>
            requires requires(Generator &generator) {
                { generator() } -> std::convertible_to<polynomial_type>;
            }
        polynomial_square_root_context(const polynomial_divisor_context<Backend> &divisor_context,
                                       polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                       Generator &generator) : divisor_context_(&divisor_context) {
            initialize_order_decomposition();

            polynomial_type candidate;
            do {
                candidate = generator();
                require_canonical_reduced(candidate);
            } while (detail::is_square_mod_impl(candidate, divisor_context, arithmetic_context, &order_decomposition_));
            cache_non_residue_power(candidate, arithmetic_context);
        }

        const boost::multiprecision::cpp_int &odd_order() const {
            return order_decomposition_.odd_order;
        }

        std::size_t two_adicity() const {
            return order_decomposition_.two_adicity;
        }

        const polynomial_type &non_residue_to_odd_order() const {
            return non_residue_to_odd_order_;
        }

        const polynomial_divisor_context<Backend> &divisor_context() const {
            return *divisor_context_;
        }

    private:
        void initialize_order_decomposition() {
            const std::size_t extension_degree = divisor_context().degree();
            if (extension_degree == 0) {
                throw std::invalid_argument("polynomial square roots require a nonconstant divisor");
            }
            if (algebra::fields::field_characteristic<field_type>() == 2) {
                throw std::invalid_argument("polynomial square roots in characteristic two are not implemented");
            }
            order_decomposition_ =
                algebra::fields::extension_field_multiplicative_group_decomposition<field_type>(extension_degree);
        }

        void require_canonical_reduced(const polynomial_type &candidate) const {
            if (candidate.empty() || candidate.size() > divisor_context().degree() ||
                (candidate.size() > 1 && candidate.back() == value_type::zero())) {
                throw std::invalid_argument(
                    "a polynomial square-root nonresidue candidate must be a nonempty canonical reduced "
                    "representative");
            }
        }

        void cache_non_residue_power(const polynomial_type &quadratic_non_residue,
                                     polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
            powmod(non_residue_to_odd_order_, quadratic_non_residue, odd_order(), divisor_context(),
                   arithmetic_context);
        }

        const polynomial_divisor_context<Backend> *divisor_context_;
        algebra::fields::multiplicative_group_decomposition order_decomposition_;
        polynomial_type non_residue_to_odd_order_;
    };

    /**
     * Return whether input is a square in K[X]/(B), where K is a finite field and B is the irreducible polynomial
     * stored in divisor_context. If Q is the order of K and d = degree(B), the quotient is a field of order Q^d.
     * For a nonzero residue a in odd characteristic, Euler's criterion gives
     *
     *     a^((Q^d - 1) / 2) = 1
     *
     * exactly when a is a square. Zero is a square. In characteristic two, squaring is an automorphism of every
     * finite field, so every residue is a square.
     *
     * The canonical indeterminate X uses the closed-form field norm
     *
     *     Norm(X) = (-1)^d B(0) / leading_coefficient(B).
     *
     * An element of a finite extension of an odd-order field is a square exactly when its norm is a square in the
     * coefficient field. When the coefficient value type provides is_square(), this avoids quotient-ring
     * exponentiation for X. Other inputs use Euler's criterion as a reference implementation.
     *
     * The irreducibility precondition is essential: for a reducible B the quotient has zero divisors and Euler's
     * criterion does not characterize its squares. The function does not repeat an irreducibility test.
     *
     * @throws std::invalid_argument if B is constant or input is not a reduced quotient-field representative.
     * @pre B is irreducible and input is a nonempty canonical coefficient polynomial.
     */
    namespace detail {

        template<SupportsDivrem Backend>
            requires algebra::FieldValue<typename Backend::polynomial_type::value_type>
        bool is_square_mod_impl(const typename Backend::polynomial_type &input,
                                const polynomial_divisor_context<Backend> &divisor_context,
                                polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                const algebra::fields::multiplicative_group_decomposition *order_decomposition) {
            using polynomial_type = typename Backend::polynomial_type;
            using value_type = typename polynomial_type::value_type;
            using field_type = typename value_type::field_type;

            const std::size_t extension_degree = divisor_context.degree();
            if (extension_degree == 0) {
                throw std::invalid_argument("square testing modulo a polynomial requires a nonconstant divisor");
            }
            if (input.size() > extension_degree) {
                throw std::invalid_argument("square testing requires a reduced quotient-field representative");
            }
            if (input.size() == 1 && input[0] == value_type::zero()) {
                return true;
            }
            if (algebra::fields::field_characteristic<field_type>() == 2) {
                return true;
            }

            // Use the norm shortcut when the coefficient field provides its own square test. Other field-value types
            // continue to the generic Euler test below.
            if constexpr (requires(const value_type &value) {
                              { value.is_square() } -> std::convertible_to<bool>;
                          }) {
                if (input.size() == 2 && input[0] == value_type::zero() && input[1] == value_type::one()) {
                    const polynomial_type &divisor = divisor_context.divisor();
                    // The constant coefficient of the inverse reversed divisor is the cached inverse of the leading
                    // coefficient of the divisor.
                    value_type norm = divisor[0] * divisor_context.reversed_divisor_inverse()[0];
                    if (extension_degree % 2 != 0) {
                        norm = value_type::zero() - norm;
                    }
                    return norm.is_square();
                }
            }

            boost::multiprecision::cpp_int exponent;
            if (order_decomposition != nullptr) {
                // In odd characteristic the two-adicity is positive, so half the group order is
                // odd_order * 2^(two_adicity - 1).
                exponent = order_decomposition->odd_order;
                exponent <<= order_decomposition->two_adicity - 1;
            } else {
                exponent = (algebra::fields::extension_field_order<field_type>(extension_degree) - 1) >> 1;
            }
            polynomial_type quadratic_character;
            powmod(quadratic_character, input, exponent, divisor_context, arithmetic_context);
            return quadratic_character.size() == 1 && quadratic_character[0] == value_type::one();
        }

    }    // namespace detail

    template<detail::SupportsDivrem Backend>
        requires algebra::FieldValue<typename Backend::polynomial_type::value_type>
    bool is_square_mod(const typename Backend::polynomial_type &input,
                       const polynomial_divisor_context<Backend> &divisor_context,
                       polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        return detail::is_square_mod_impl(input, divisor_context, arithmetic_context, nullptr);
    }

    /**
     * Compute a square root of input in K[X]/(B) with Tonelli-Shanks, using the multiplicative-group decomposition and
     * quadratic nonresidue cached in square_root_context. The irreducible polynomial B is stored in the divisor context
     * referenced by square_root_context. The result is canonical and may alias input.
     *
     * The function returns false when input is not a square. In that case output is set to the zero polynomial. Zero
     * is its own square root.
     *
     * @throws std::invalid_argument if input is not a reduced quotient-field representative.
     * @pre B is irreducible and input is a nonempty canonical coefficient polynomial.
     */
    template<detail::SupportsDivrem Backend>
        requires algebra::FieldValue<typename Backend::polynomial_type::value_type>
    bool square_root_mod(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input,
                         const polynomial_square_root_context<Backend> &square_root_context,
                         polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;
        const polynomial_divisor_context<Backend> &divisor_context = square_root_context.divisor_context();

        if (input.size() > divisor_context.degree()) {
            throw std::invalid_argument("a polynomial square root requires a reduced quotient-field representative");
        }
        if (input.size() == 1 && input[0] == value_type::zero()) {
            output.assign(1, value_type::zero());
            return true;
        }

        polynomial_type root;
        polynomial_type input_to_odd_order;
        polynomial_type non_residue_power = square_root_context.non_residue_to_odd_order();

        // Write the quotient field's multiplicative-group order as odd_order * 2^S. For input a and the cached
        // nonresidue z, Tonelli-Shanks starts with
        //
        //     input_to_odd_order = a^odd_order,
        //     root               = a^((odd_order + 1) / 2),
        //     non_residue_power  = z^odd_order.
        //
        // The first two values satisfy root^2 = a * input_to_odd_order.
        powmod(input_to_odd_order, input, square_root_context.odd_order(), divisor_context, arithmetic_context);
        boost::multiprecision::cpp_int root_exponent = square_root_context.odd_order() + 1;
        root_exponent >>= 1;
        powmod(root, input, root_exponent, divisor_context, arithmetic_context);

        std::size_t remaining_two_adicity = square_root_context.two_adicity();
        const polynomial_type one = {value_type::one()};
        while (input_to_odd_order != one) {
            // The odd-order exponent places this value in the power-of-two subgroup. Find the smallest i for which
            //
            //     input_to_odd_order^(2^i) = 1.
            //
            // A square must reach one for some i < M, where M = remaining_two_adicity and its current order divides
            // 2^M.
            polynomial_type power = input_to_odd_order;
            std::size_t first_one_power = 0;
            do {
                squaremod(power, power, divisor_context, arithmetic_context);
                ++first_one_power;
            } while (first_one_power < remaining_two_adicity && power != one);

            if (first_one_power == remaining_two_adicity) {
                output.assign(1, value_type::zero());
                return false;
            }

            // Set correction = non_residue_power^(2^(M - i - 1)), then update
            //
            //     non_residue_power = correction^2,
            //     input_to_odd_order = input_to_odd_order * correction^2,
            //     root               = root * correction,
            //     M                  = i.
            //
            // This preserves root^2 = input * input_to_odd_order while reducing input_to_odd_order's two-power order
            // from 2^M to 2^i.
            // The previous nonresidue power is no longer needed, so move its storage into the correction.
            polynomial_type correction = std::move(non_residue_power);
            for (std::size_t i = 0; i < remaining_two_adicity - first_one_power - 1; ++i) {
                squaremod(correction, correction, divisor_context, arithmetic_context);
            }

            squaremod(non_residue_power, correction, divisor_context, arithmetic_context);
            mulmod(input_to_odd_order, input_to_odd_order, non_residue_power, divisor_context, arithmetic_context);
            mulmod(root, root, correction, divisor_context, arithmetic_context);
            remaining_two_adicity = first_one_power;
        }

        output = std::move(root);
        return true;
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_SQUARE_ROOT_HPP
