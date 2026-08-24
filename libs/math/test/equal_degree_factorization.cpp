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

#define BOOST_TEST_MODULE equal_degree_factorization_test

#include <array>
#include <cstddef>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/equal_degree_factorization.hpp>
#include <nil/crypto3/math/polynomial/schoolbook_backend.hpp>
#include <nil/crypto3/random/algebraic_engine.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using field_type = fields::alt_bn128_base_field<254>;
    using value_type = field_type::value_type;
    using extension_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using extension_value_type = extension_field_type::value_type;
    using backend_type = polynomial_arithmetic::schoolbook_backend<value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
}    // namespace

BOOST_AUTO_TEST_SUITE(equal_degree_factorization_test_suite)

BOOST_AUTO_TEST_CASE(random_polynomial_sampling_is_seeded) {
    nil::crypto3::random::algebraic_engine<field_type> first_generator(17);
    nil::crypto3::random::algebraic_engine<field_type> second_generator(17);

    const polynomial_type first = math::detail::sample_random_polynomial<polynomial_type>(8, first_generator);
    const polynomial_type second = math::detail::sample_random_polynomial<polynomial_type>(8, second_generator);

    BOOST_CHECK(first == second);
    BOOST_CHECK_LE(first.size(), 8);
    BOOST_CHECK_THROW(math::detail::sample_random_polynomial<polynomial_type>(0, first_generator),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(random_polynomial_sampling_returns_canonical_polynomials) {
    const std::array coefficients = {value_type(2), value_type(3), value_type::zero(), value_type::zero()};
    std::size_t next_coefficient = 0;
    auto generator = [&] { return coefficients[next_coefficient++]; };

    const polynomial_type sampled = math::detail::sample_random_polynomial<polynomial_type>(4, generator);

    BOOST_CHECK(sampled == polynomial_type({value_type(2), value_type(3)}));
    BOOST_CHECK_EQUAL(next_coefficient, coefficients.size());
}

BOOST_AUTO_TEST_CASE(random_polynomial_sampling_supports_extension_field_coefficients) {
    using extension_polynomial_type = math::polynomial<extension_value_type>;

    nil::crypto3::random::algebraic_engine<extension_field_type> first_generator(29);
    nil::crypto3::random::algebraic_engine<extension_field_type> second_generator(29);
    const extension_polynomial_type first =
        math::detail::sample_random_polynomial<extension_polynomial_type>(5, first_generator);
    const extension_polynomial_type second =
        math::detail::sample_random_polynomial<extension_polynomial_type>(5, second_generator);

    BOOST_CHECK(first == second);
    BOOST_CHECK_LE(first.size(), 5);
}

BOOST_AUTO_TEST_CASE(one_cantor_zassenhaus_trial_finds_a_proper_factor) {
    backend_type backend;
    const polynomial_type first_factor = {-value_type::one(), value_type::one()};
    const polynomial_type second_factor = {value_type::one(), value_type::one()};
    polynomial_type group;
    backend.multiply(group, first_factor, second_factor);

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(group, group.size() - 1, arithmetic_context);
    math::detail::cantor_zassenhaus_context<backend_type> split_context(1);

    const std::array coefficients = {value_type::zero(), value_type::one()};
    std::size_t next_coefficient = 0;
    auto generator = [&] { return coefficients[next_coefficient++]; };
    polynomial_type factor;
    const bool split = math::detail::try_cantor_zassenhaus_split<backend_type>(factor, split_context, divisor_context,
                                                                               arithmetic_context, generator);

    BOOST_CHECK(split);
    BOOST_CHECK(factor == first_factor || factor == second_factor);
    BOOST_CHECK_EQUAL(next_coefficient, coefficients.size());
}

BOOST_AUTO_TEST_CASE(one_cantor_zassenhaus_trial_reports_an_unsuccessful_sample) {
    backend_type backend;
    const polynomial_type first_factor = {-value_type::one(), value_type::one()};
    const polynomial_type second_factor = {value_type::one(), value_type::one()};
    polynomial_type group;
    backend.multiply(group, first_factor, second_factor);

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(group, group.size() - 1, arithmetic_context);
    math::detail::cantor_zassenhaus_context<backend_type> split_context(1);

    const std::array coefficients = {value_type::one(), value_type::zero()};
    std::size_t next_coefficient = 0;
    auto generator = [&] { return coefficients[next_coefficient++]; };
    polynomial_type factor;
    const bool split = math::detail::try_cantor_zassenhaus_split<backend_type>(factor, split_context, divisor_context,
                                                                               arithmetic_context, generator);

    BOOST_CHECK(!split);
    BOOST_CHECK(factor == polynomial_type({value_type::zero()}));
}

BOOST_AUTO_TEST_CASE(cantor_zassenhaus_trial_rejects_invalid_inputs_before_sampling) {
    BOOST_CHECK_THROW(static_cast<void>(math::detail::cantor_zassenhaus_context<backend_type>(0)),
                      std::invalid_argument);

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    polynomial_type factor;
    std::size_t generated_coefficient_count = 0;
    auto generator = [&] {
        ++generated_coefficient_count;
        return value_type::one();
    };

    const polynomial_type degree_two_group = {-value_type::one(), value_type::zero(), value_type::one()};
    math::polynomial_divisor_context<backend_type> degree_two_divisor(degree_two_group, degree_two_group.size() - 1,
                                                                      arithmetic_context);
    math::detail::cantor_zassenhaus_context<backend_type> degree_two_split_context(2);
    BOOST_CHECK_THROW(math::detail::try_cantor_zassenhaus_split<backend_type>(
                          factor, degree_two_split_context, degree_two_divisor, arithmetic_context, generator),
                      std::invalid_argument);

    const polynomial_type degree_three_group = {value_type::one(), value_type::zero(), value_type::zero(),
                                                value_type::one()};
    math::polynomial_divisor_context<backend_type> degree_three_divisor(
        degree_three_group, degree_three_group.size() - 1, arithmetic_context);
    BOOST_CHECK_THROW(math::detail::try_cantor_zassenhaus_split<backend_type>(
                          factor, degree_two_split_context, degree_three_divisor, arithmetic_context, generator),
                      std::invalid_argument);

    const polynomial_type nonmonic_group = {-value_type(2), value_type::zero(), value_type(2)};
    math::polynomial_divisor_context<backend_type> nonmonic_divisor(nonmonic_group, nonmonic_group.size() - 1,
                                                                    arithmetic_context);
    math::detail::cantor_zassenhaus_context<backend_type> degree_one_split_context(1);
    BOOST_CHECK_THROW(math::detail::try_cantor_zassenhaus_split<backend_type>(
                          factor, degree_one_split_context, nonmonic_divisor, arithmetic_context, generator),
                      std::invalid_argument);

    BOOST_CHECK_EQUAL(generated_coefficient_count, 0);
}

BOOST_AUTO_TEST_CASE(preparation_normalizes_a_square_free_equal_degree_input) {
    backend_type backend;
    const polynomial_type first_factor = {-value_type(3), value_type::zero(), value_type::one()};
    // Three is a nonsquare in BN254 Fq. Multiplying it by the nonzero square four preserves that property, so both
    // quadratics are irreducible.
    const polynomial_type second_factor = {-value_type(12), value_type::zero(), value_type::one()};
    polynomial_type input;
    backend.multiply(input, first_factor, second_factor);

    const value_type leading_coefficient(11);
    for (value_type &coefficient : input) {
        coefficient *= leading_coefficient;
    }

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    polynomial_type monic_input;
    value_type recovered_leading_coefficient;
    const bool has_factors = math::detail::prepare_equal_degree_factorization_input<backend_type>(
        monic_input, recovered_leading_coefficient, input, 2, arithmetic_context);

    polynomial_type expected;
    backend.multiply(expected, first_factor, second_factor);
    BOOST_CHECK(has_factors);
    BOOST_CHECK(recovered_leading_coefficient == leading_coefficient);
    BOOST_CHECK(monic_input == expected);
}

BOOST_AUTO_TEST_CASE(preparation_rejects_invalid_equal_degree_inputs) {
    backend_type backend;
    const polynomial_type linear_factor = {value_type(1), value_type::one()};
    const polynomial_type quadratic_factor = {-value_type(3), value_type::zero(), value_type::one()};
    polynomial_type mixed_degree_input;
    backend.multiply(mixed_degree_input, linear_factor, quadratic_factor);

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    polynomial_type monic_input;
    value_type leading_coefficient;
    BOOST_CHECK_THROW(math::detail::prepare_equal_degree_factorization_input<backend_type>(
                          monic_input, leading_coefficient, mixed_degree_input, 0, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::prepare_equal_degree_factorization_input<backend_type>(
                          monic_input, leading_coefficient, mixed_degree_input, 2, arithmetic_context),
                      std::invalid_argument);

    polynomial_type repeated_input;
    backend.multiply(repeated_input, linear_factor, linear_factor);
    BOOST_CHECK_THROW(math::detail::prepare_equal_degree_factorization_input<backend_type>(
                          monic_input, leading_coefficient, repeated_input, 1, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(preparation_reports_that_constants_have_no_factors) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    polynomial_type monic_input;
    value_type leading_coefficient;
    const bool has_factors = math::detail::prepare_equal_degree_factorization_input<backend_type>(
        monic_input, leading_coefficient, polynomial_type {value_type(13)}, 1, arithmetic_context);

    BOOST_CHECK(!has_factors);
    BOOST_CHECK(leading_coefficient == value_type(13));
    BOOST_CHECK(monic_input == polynomial_type({value_type(13)}));
}

BOOST_AUTO_TEST_SUITE_END()
