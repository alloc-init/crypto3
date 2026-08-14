//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// MIT License
//---------------------------------------------------------------------------//

#define BOOST_TEST_MODULE linear_combination_test

#include <boost/test/unit_test.hpp>

#include <sstream>
#include <vector>

#include <nil/crypto3/algebra/fields/arithmetic_params/bls12.hpp>
#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/fp2.hpp>
#include <nil/crypto3/algebra/fields/fp6_3over2.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/math/linear_combination.hpp>
#include <nil/crypto3/math/linear_variable.hpp>

using field_type = nil::crypto3::algebra::fields::bls12_fr<381>;
using value_type = field_type::value_type;
using variable_type = nil::crypto3::math::linear_variable<field_type>;
using term_type = nil::crypto3::math::linear_term<variable_type>;
using combination_type = nil::crypto3::math::linear_combination<variable_type>;
using explicit_constant_combination_type =
    nil::crypto3::math::linear_combination<variable_type, nil::crypto3::math::assignment_layout::explicit_constant>;

using base_field_type = nil::crypto3::algebra::fields::alt_bn128_base_field<254>;
using fp2_field_type = nil::crypto3::algebra::fields::fp2<base_field_type>;
using fp6_field_type = nil::crypto3::algebra::fields::fp6_3over2<base_field_type>;
using fp12_field_type = nil::crypto3::algebra::fields::fp12_2over3over2<base_field_type>;
using base_value_type = base_field_type::value_type;
using fp2_value_type = fp2_field_type::value_type;
using fp6_value_type = fp6_field_type::value_type;
using fp12_value_type = fp12_field_type::value_type;

BOOST_AUTO_TEST_CASE(construction_and_normalization) {
    combination_type from_index(variable_type(3));
    BOOST_REQUIRE_EQUAL(from_index.terms.size(), 1);
    BOOST_CHECK_EQUAL(from_index.terms[0].index, 3);
    BOOST_CHECK_EQUAL(from_index.terms[0].coeff, value_type(1));

    const variable_type x1(1);
    const variable_type x3(3);
    combination_type normalized(std::vector<term_type> {x3 * value_type(2), x1 * value_type(4), x3 * value_type(5)});

    BOOST_REQUIRE_EQUAL(normalized.terms.size(), 2);
    BOOST_CHECK_EQUAL(normalized.terms[0].index, 1);
    BOOST_CHECK_EQUAL(normalized.terms[0].coeff, value_type(4));
    BOOST_CHECK_EQUAL(normalized.terms[1].index, 3);
    BOOST_CHECK_EQUAL(normalized.terms[1].coeff, value_type(7));
}

BOOST_AUTO_TEST_CASE(evaluation_layouts) {
    combination_type combination(value_type(2));
    combination.add_term(variable_type(2), value_type(3));

    const std::vector<value_type> assignment {value_type(5), value_type(7)};
    BOOST_CHECK_EQUAL(combination.evaluate(assignment), value_type(23));

    explicit_constant_combination_type explicit_combination(value_type(2));
    explicit_combination.add_term(variable_type(2), value_type(3));
    const std::vector<value_type> witness {value_type(10), value_type(20), value_type(30)};
    BOOST_CHECK_EQUAL(explicit_combination.evaluate(witness), value_type(110));
}

BOOST_AUTO_TEST_CASE(evaluation_over_extension_field) {
    using base_variable_type = nil::crypto3::math::linear_variable<base_field_type>;
    using base_combination_type =
        nil::crypto3::math::linear_combination<base_variable_type,
                                               nil::crypto3::math::assignment_layout::explicit_constant>;

    base_combination_type combination(base_value_type(2));
    combination.add_term(base_variable_type(2), base_value_type(3));

    const fp6_value_type extension_component(fp2_value_type::zero(), fp2_value_type::one(), fp2_value_type::zero());
    const std::vector<fp12_value_type> witness {
        fp12_value_type(fp6_value_type::one(), extension_component),
        fp12_value_type::one(),
        fp12_value_type(extension_component, fp6_value_type::one()),
    };

    const auto expected = base_value_type(2) * witness[0] + base_value_type(3) * witness[2];
    BOOST_CHECK_EQUAL(combination.evaluate(witness), expected);
}

BOOST_AUTO_TEST_CASE(arithmetic_and_scaled_append) {
    combination_type left(variable_type(1) * value_type(2));
    combination_type right(variable_type(2) * value_type(3));

    const auto sum = left + right;
    BOOST_REQUIRE_EQUAL(sum.terms.size(), 2);
    BOOST_CHECK_EQUAL(sum.terms[0].coeff, value_type(2));
    BOOST_CHECK_EQUAL(sum.terms[1].coeff, value_type(3));

    combination_type scaled;
    nil::crypto3::math::add_scaled(scaled, value_type(4), sum);
    BOOST_REQUIRE_EQUAL(scaled.terms.size(), 2);
    BOOST_CHECK_EQUAL(scaled.terms[0].coeff, value_type(8));
    BOOST_CHECK_EQUAL(scaled.terms[1].coeff, value_type(12));

    const auto negated = -variable_type(1);
    BOOST_CHECK_EQUAL(negated.index, 1);
    BOOST_CHECK_EQUAL(negated.coeff, value_type(-1));
}

BOOST_AUTO_TEST_CASE(validity_and_stream_output) {
    combination_type empty;
    BOOST_CHECK(empty.is_valid(2));

    combination_type sorted(std::vector<term_type> {variable_type(0), variable_type(1)});
    BOOST_CHECK(sorted.is_valid(2));

    combination_type unsorted;
    unsorted.add_term(variable_type(2));
    unsorted.add_term(variable_type(1));
    BOOST_CHECK(!unsorted.is_valid(3));

    std::ostringstream out;
    out << sorted;
    std::ostringstream expected;
    expected << value_type(1) << " * v0 + " << value_type(1) << " * v1";
    BOOST_CHECK_EQUAL(out.str(), expected.str());
}
