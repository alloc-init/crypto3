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

#define BOOST_TEST_MODULE polynomial_square_root_test

#include <cstddef>
#include <stdexcept>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>

#include <nil/crypto3/math/polynomial/polynomial_square_root.hpp>
#include <nil/crypto3/math/polynomial/schoolbook_backend.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace fields = nil::crypto3::algebra::fields;
    namespace polynomial_arithmetic = math::polynomial_arithmetic;

    using field_type = fields::babybear;
    using value_type = field_type::value_type;
    using backend_type = polynomial_arithmetic::schoolbook_backend<value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    value_type first_quadratic_non_residue() {
        value_type candidate(2);
        while (candidate.is_square()) {
            candidate = candidate + value_type::one();
        }
        return candidate;
    }

    value_type first_cubic_non_residue() {
        const field_type::integral_type cubic_character_exponent = (field_type::modulus - 1u) / 3u;
        value_type candidate(2);
        while (candidate.pow(cubic_character_exponent) == value_type::one()) {
            candidate = candidate + value_type::one();
        }
        return candidate;
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(polynomial_square_root_test_suite)

BOOST_AUTO_TEST_CASE(square_root_context_validates_and_caches_tonelli_shanks_parameters) {
    const value_type base_non_residue = first_quadratic_non_residue();
    const polynomial_type irreducible_divisor = {-base_non_residue, value_type::zero(), value_type::one()};
    // BabyBear has order 1 modulo 4, so the square root X of a base-field nonsquare remains a nonsquare in this
    // quadratic extension.
    const polynomial_type quotient_non_residue = {value_type::zero(), value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(irreducible_divisor, 1, arithmetic_context);
    const math::polynomial_square_root_context<backend_type> square_root_context(quotient_non_residue, divisor_context,
                                                                                 arithmetic_context);

    const boost::multiprecision::cpp_int quotient_order = fields::extension_field_order<field_type>(2);
    BOOST_CHECK((square_root_context.odd_order() & 1) == 1);
    BOOST_CHECK((square_root_context.odd_order() << square_root_context.two_adicity()) == quotient_order - 1);

    polynomial_type expected_non_residue_to_odd_order;
    math::powmod(expected_non_residue_to_odd_order, quotient_non_residue, square_root_context.odd_order(),
                 divisor_context, arithmetic_context);
    BOOST_CHECK(square_root_context.non_residue_to_odd_order() == expected_non_residue_to_odd_order);

    BOOST_CHECK_THROW((math::polynomial_square_root_context<backend_type>(polynomial_type {value_type::one()},
                                                                          divisor_context, arithmetic_context)),
                      std::invalid_argument);
    BOOST_CHECK_THROW((math::polynomial_square_root_context<backend_type>(
                          polynomial_type {value_type::zero(), value_type::zero(), value_type::one()}, divisor_context,
                          arithmetic_context)),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(euler_criterion_classifies_squares_in_a_quadratic_quotient_field) {
    const value_type non_residue = first_quadratic_non_residue();
    const polynomial_type irreducible_divisor = {-non_residue, value_type::zero(), value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(irreducible_divisor, 1, arithmetic_context);

    const polynomial_type value = {value_type(3), value_type(5)};
    polynomial_type square;
    math::squaremod(square, value, divisor_context, arithmetic_context);

    BOOST_CHECK(math::is_square_mod(square, divisor_context, arithmetic_context));
    BOOST_CHECK(math::is_square_mod(polynomial_type {value_type::zero()}, divisor_context, arithmetic_context));
    BOOST_CHECK(math::is_square_mod(polynomial_type {value_type::one()}, divisor_context, arithmetic_context));

    // BabyBear has order 1 modulo 4. In a quadratic extension of such a field, the square root X of a base-field
    // nonsquare remains a nonsquare in the extension.
    const polynomial_type outer_generator = {value_type::zero(), value_type::one()};
    BOOST_CHECK(!math::is_square_mod(outer_generator, divisor_context, arithmetic_context));
}

BOOST_AUTO_TEST_CASE(indeterminate_square_test_uses_its_norm_for_monic_and_nonmonic_divisors) {
    const value_type non_residue = first_quadratic_non_residue();
    const polynomial_type indeterminate = {value_type::zero(), value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;

    const polynomial_type monic_quadratic = {-non_residue, value_type::zero(), value_type::one()};
    math::polynomial_divisor_context<backend_type> monic_context(monic_quadratic, 1, arithmetic_context);
    BOOST_CHECK(!math::is_square_mod(indeterminate, monic_context, arithmetic_context));

    const value_type scale(7);
    const polynomial_type nonmonic_quadratic = {-scale * non_residue, value_type::zero(), scale};
    math::polynomial_divisor_context<backend_type> nonmonic_context(nonmonic_quadratic, 1, arithmetic_context);
    BOOST_CHECK(!math::is_square_mod(indeterminate, nonmonic_context, arithmetic_context));

    // BabyBear contains the cube roots of unity, so X^3 - c is irreducible when c is not a cube. For this odd-degree
    // divisor, Norm(X) = c, which also exercises the sign in the constant-term formula.
    const value_type cubic_non_residue = first_cubic_non_residue();
    const polynomial_type monic_cubic = {-cubic_non_residue, value_type::zero(), value_type::zero(), value_type::one()};
    math::polynomial_divisor_context<backend_type> cubic_context(monic_cubic, 2, arithmetic_context);
    BOOST_CHECK_EQUAL(math::is_square_mod(indeterminate, cubic_context, arithmetic_context),
                      cubic_non_residue.is_square());
}

BOOST_AUTO_TEST_CASE(square_testing_rejects_inputs_outside_the_quotient_field_contract) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type linear_divisor = {value_type::one(), value_type::one()};
    math::polynomial_divisor_context<backend_type> linear_context(linear_divisor, 1, arithmetic_context);
    BOOST_CHECK_THROW(
        math::is_square_mod(polynomial_type {value_type::one(), value_type::one()}, linear_context, arithmetic_context),
        std::invalid_argument);

    const polynomial_type constant_divisor = {value_type::one()};
    math::polynomial_divisor_context<backend_type> constant_context(constant_divisor, 1, arithmetic_context);
    BOOST_CHECK_THROW((math::polynomial_square_root_context<backend_type>(polynomial_type {value_type::one()},
                                                                          constant_context, arithmetic_context)),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::is_square_mod(polynomial_type {value_type::one()}, constant_context, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(tonelli_shanks_recovers_square_roots_in_a_quotient_field) {
    const value_type base_non_residue = first_quadratic_non_residue();
    const polynomial_type irreducible_divisor = {-base_non_residue, value_type::zero(), value_type::one()};
    const polynomial_type quotient_non_residue = {value_type::zero(), value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(irreducible_divisor, 1, arithmetic_context);
    const math::polynomial_square_root_context<backend_type> square_root_context(quotient_non_residue, divisor_context,
                                                                                 arithmetic_context);

    const polynomial_type value = {value_type(3), value_type(5)};
    polynomial_type square;
    math::squaremod(square, value, divisor_context, arithmetic_context);

    polynomial_type root;
    BOOST_REQUIRE(math::square_root_mod(root, square, square_root_context, arithmetic_context));
    polynomial_type recovered_square;
    math::squaremod(recovered_square, root, divisor_context, arithmetic_context);
    BOOST_CHECK(recovered_square == square);

    polynomial_type aliased_root = square;
    BOOST_REQUIRE(math::square_root_mod(aliased_root, aliased_root, square_root_context, arithmetic_context));
    math::squaremod(recovered_square, aliased_root, divisor_context, arithmetic_context);
    BOOST_CHECK(recovered_square == square);

    const polynomial_type zero = {value_type::zero()};
    BOOST_REQUIRE(math::square_root_mod(root, zero, square_root_context, arithmetic_context));
    BOOST_CHECK(root == zero);

    BOOST_CHECK(!math::square_root_mod(root, quotient_non_residue, square_root_context, arithmetic_context));
    BOOST_CHECK(root == zero);

    BOOST_CHECK_THROW(math::square_root_mod(root, irreducible_divisor, square_root_context, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
