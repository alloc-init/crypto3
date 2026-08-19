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

#define BOOST_TEST_MODULE polynomial_composition_test

#include <cstddef>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/polynomial_composition.hpp>
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

BOOST_AUTO_TEST_SUITE(polynomial_composition_test_suite)

BOOST_AUTO_TEST_CASE(reference_composition_matches_a_known_result_and_may_alias_inputs) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    // Modulo X^2 + 1, (1 + 2Y + 3Y^2) evaluated at Y = X is -2 + 2X.
    const polynomial_type outer = {fq_value_type(1), fq_value_type(2), fq_value_type(3)};
    const polynomial_type inner = {fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type expected = {-fq_value_type(2), fq_value_type(2)};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);

    polynomial_type result;
    math::compose_mod_reference(result, outer, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(result == expected);

    polynomial_type outer_alias = outer;
    math::compose_mod_reference(outer_alias, outer_alias, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(outer_alias == expected);

    polynomial_type inner_alias = inner;
    math::compose_mod_reference(inner_alias, outer, inner_alias, divisor_context, arithmetic_context);
    BOOST_CHECK(inner_alias == expected);
}

BOOST_AUTO_TEST_CASE(reference_composition_reduces_the_inner_polynomial) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type outer = {fq_value_type(2), fq_value_type(3), fq_value_type(4)};
    const polynomial_type inner = {fq_value_type(5), fq_value_type(6), fq_value_type(7), fq_value_type(8)};
    const polynomial_type divisor = {fq_value_type(1), fq_value_type(1), fq_value_type(1)};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 2, arithmetic_context);

    polynomial_type reduced_inner;
    math::remainder(reduced_inner, inner, divisor_context, arithmetic_context);

    polynomial_type from_unreduced;
    polynomial_type from_reduced;
    math::compose_mod_reference(from_unreduced, outer, inner, divisor_context, arithmetic_context);
    math::compose_mod_reference(from_reduced, outer, reduced_inner, divisor_context, arithmetic_context);
    BOOST_CHECK(from_unreduced == from_reduced);

    polynomial_type inner_alias = inner;
    math::compose_mod_reference(inner_alias, outer, inner_alias, divisor_context, arithmetic_context);
    BOOST_CHECK(inner_alias == from_unreduced);
}

BOOST_AUTO_TEST_CASE(reference_composition_handles_zero_constant_and_constant_modulus_cases) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type(1)};
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);
    const polynomial_type inner = {fq_value_type(4), fq_value_type(5)};

    polynomial_type result;
    math::compose_mod_reference(result, polynomial_type {fq_value_type::zero()}, inner, divisor_context,
                                arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));

    math::compose_mod_reference(result, polynomial_type {fq_value_type(7)}, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type(7)}));

    const polynomial_type constant_divisor = {fq_value_type(11)};
    math::polynomial_divisor_context<backend_type> constant_context(constant_divisor, 1, arithmetic_context);
    math::compose_mod_reference(result, polynomial_type {fq_value_type(1), fq_value_type(2)}, inner, constant_context,
                                arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));
}

BOOST_AUTO_TEST_CASE(reference_composition_rejects_insufficient_inverse_precision_without_changing_output) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type {}, options);

    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);
    const polynomial_type outer = {fq_value_type(2), fq_value_type(3)};
    const polynomial_type inner = {fq_value_type(4), fq_value_type(5), fq_value_type(6), fq_value_type(7)};

    polynomial_type output = {fq_value_type(17)};
    BOOST_CHECK_THROW(math::compose_mod_reference(output, outer, inner, divisor_context, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK(output == polynomial_type({fq_value_type(17)}));
}

BOOST_AUTO_TEST_CASE(reference_composition_matches_between_schoolbook_and_mixed_radix_backends) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    const polynomial_type outer = {fq12_value(1), fq12_value(13), fq12_value(25), fq12_value(37)};
    const polynomial_type inner = {fq12_value(49), fq12_value(61), fq12_value(73)};
    const polynomial_type divisor = {fq12_value(85), fq12_value(97), fq12_value(109), fq12_value(121)};

    polynomial_arithmetic::polynomial_context<schoolbook_backend> schoolbook_context;
    math::polynomial_divisor_context<schoolbook_backend> schoolbook_divisor(divisor, 2, schoolbook_context);
    polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, schoolbook_divisor, schoolbook_context);

    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<mixed_radix_backend> mixed_radix_context(mixed_radix_backend(18),
                                                                                       options);
    math::polynomial_divisor_context<mixed_radix_backend> mixed_radix_divisor(divisor, 2, mixed_radix_context);
    polynomial_type result;
    math::compose_mod_reference(result, outer, inner, mixed_radix_divisor, mixed_radix_context);
    BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_SUITE_END()
