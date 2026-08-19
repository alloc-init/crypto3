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

#define BOOST_TEST_MODULE polynomial_frobenius_test

#include <cstddef>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/field_order.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/polynomial_frobenius.hpp>
#include <nil/crypto3/math/polynomial/schoolbook_backend.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq_value_type = fq_field_type::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;

    fq12_value_type fq12_value(std::size_t first_coordinate) {
        fq12_value_type value = fq12_value_type::zero();
        for (std::size_t i = 0; i < fq12_field_type::arity; ++i) {
            value.coordinate(i) = fq_value_type(first_coordinate + i);
        }
        return value;
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(polynomial_frobenius_test_suite)

BOOST_AUTO_TEST_CASE(precomputation_matches_direct_exponentiation) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type x = {fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    polynomial_type expected;
    math::powmod(expected, x, fields::field_order<fq_field_type>(), frobenius_context.divisor_context(),
                 arithmetic_context);
    BOOST_CHECK(frobenius_context.x_to_field_order() == expected);
}

BOOST_AUTO_TEST_CASE(frobenius_matches_direct_exponentiation_and_may_alias_the_input) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type input = {fq_value_type(5), fq_value_type(7), fq_value_type(11)};
    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    polynomial_type expected;
    math::powmod(expected, input, fields::field_order<fq_field_type>(), frobenius_context.divisor_context(),
                 arithmetic_context);

    polynomial_type result;
    math::frobenius_map(result, input, frobenius_context, arithmetic_context);
    BOOST_CHECK(result == expected);

    polynomial_type aliased = input;
    math::frobenius_map(aliased, aliased, frobenius_context, arithmetic_context);
    BOOST_CHECK(aliased == expected);

    const polynomial_type unreduced = {fq_value_type(13), fq_value_type(17), fq_value_type(19), fq_value_type(23),
                                       fq_value_type(29)};
    math::powmod(expected, unreduced, fields::field_order<fq_field_type>(), frobenius_context.divisor_context(),
                 arithmetic_context);

    math::frobenius_map(result, unreduced, frobenius_context, arithmetic_context);
    BOOST_CHECK(result == expected);

    aliased = unreduced;
    math::frobenius_map(aliased, aliased, frobenius_context, arithmetic_context);
    BOOST_CHECK(aliased == expected);
}

BOOST_AUTO_TEST_CASE(iterated_frobenius_matches_a_power_of_the_field_order) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type input = {fq_value_type(5), fq_value_type(7), fq_value_type(11)};
    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    const boost::multiprecision::cpp_int order = fields::field_order<fq_field_type>();
    const boost::multiprecision::cpp_int order_squared = order * order;
    polynomial_type expected;
    math::powmod(expected, input, order_squared, frobenius_context.divisor_context(), arithmetic_context);

    polynomial_type result;
    math::frobenius_map(result, input, 2, frobenius_context, arithmetic_context);
    BOOST_CHECK(result == expected);

    polynomial_type aliased = input;
    math::frobenius_map(aliased, aliased, 2, frobenius_context, arithmetic_context);
    BOOST_CHECK(aliased == expected);

    const polynomial_type unreduced = {fq_value_type(5), fq_value_type(7), fq_value_type(11), fq_value_type(13),
                                       fq_value_type(17)};
    polynomial_type reduced;
    math::remainder(reduced, unreduced, frobenius_context.divisor_context(), arithmetic_context);

    polynomial_type zero_iteration_alias = unreduced;
    math::frobenius_map(zero_iteration_alias, zero_iteration_alias, 0, frobenius_context, arithmetic_context);
    BOOST_CHECK(zero_iteration_alias == reduced);
}

BOOST_AUTO_TEST_CASE(frobenius_uses_the_extension_coefficient_field_order) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq12_value_type::one(), fq12_value_type::zero(), fq12_value_type::one()};
    const polynomial_type input = {fq12_value(5), fq12_value(17)};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    // X^2 = -1 and Q = q^12 is 1 modulo 4, so X^Q = X. The Q-power map also fixes every Fq12 coefficient.
    const polynomial_type x = {fq12_value_type::zero(), fq12_value_type::one()};
    BOOST_CHECK(frobenius_context.x_to_field_order() == x);

    polynomial_type result;
    math::frobenius_map(result, input, frobenius_context, arithmetic_context);
    BOOST_CHECK(result == input);
}

BOOST_AUTO_TEST_CASE(frobenius_modulo_a_nonzero_constant_is_zero) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq_value_type(7)};
    const polynomial_type input = {fq_value_type(2), fq_value_type(3)};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    polynomial_type result;
    math::frobenius_map(result, input, frobenius_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));
}

BOOST_AUTO_TEST_CASE(frobenius_context_rejects_invalid_precomputation_inputs) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type zero_divisor = {fq_value_type::zero()};
    BOOST_CHECK_THROW(math::polynomial_frobenius_context<backend_type>(zero_divisor, arithmetic_context),
                      std::invalid_argument);

    polynomial_arithmetic::polynomial_context_options invalid_options;
    invalid_options.modular_composition_cached_power_limit = 0;
    polynomial_arithmetic::polynomial_context<backend_type> invalid_arithmetic_context(backend_type {},
                                                                                       invalid_options);
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    BOOST_CHECK_THROW(math::polynomial_frobenius_context<backend_type>(divisor, invalid_arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
