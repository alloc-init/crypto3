//---------------------------------------------------------------------------//
// Copyright (c) 2026 Alloc Init
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

#define BOOST_TEST_MODULE r1cs_test

#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <utility>

#include <nil/crypto3/algebra/fields/alt_bn128/scalar_field.hpp>
#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/r1cs.hpp>

namespace {
    using field_type = nil::crypto3::algebra::fields::alt_bn128_scalar_field<254>;
    using value_type = field_type::value_type;
    using system_type = nil::crypto3::zk::snark::r1cs_constraint_system<field_type>;
    using constraint_type = nil::crypto3::zk::snark::r1cs_constraint<field_type>;
    using variable_type = nil::crypto3::math::linear_variable<field_type>;
    using combination_type = nil::crypto3::math::linear_combination<variable_type>;
}    // namespace

BOOST_AUTO_TEST_SUITE(r1cs_test_suite)

BOOST_AUTO_TEST_CASE(accepts_constant_only_system) {
    system_type cs;
    BOOST_CHECK(cs.is_valid());

    cs.add_constraint({variable_type(0), variable_type(0), variable_type(0)});
    BOOST_CHECK(cs.is_valid());
    BOOST_CHECK(cs.is_satisfied({}, {}));
}

BOOST_AUTO_TEST_CASE(accepts_last_variable_in_each_combination) {
    system_type cs;
    cs.primary_input_size = 1;
    cs.auxiliary_input_size = 2;

    const std::array<constraint_type, 3> constraints = {{
        {variable_type(3), variable_type(0), value_type(3)},
        {variable_type(0), variable_type(3), value_type(3)},
        {value_type(3), variable_type(0), variable_type(3)},
    }};
    for (std::size_t i = 0; i < constraints.size(); ++i) {
        BOOST_TEST_CONTEXT("combination " << i) {
            cs.constraints = {constraints[i]};
            BOOST_CHECK(cs.is_valid());
            BOOST_CHECK(cs.is_satisfied({value_type(1)}, {value_type(2), value_type(3)}));
        }
    }
}

BOOST_AUTO_TEST_CASE(accepts_public_only_and_private_only_assignments) {
    system_type cs;
    cs.add_constraint({variable_type(1), variable_type(0), variable_type(1)});

    cs.primary_input_size = 1;
    BOOST_CHECK(cs.is_valid());
    BOOST_CHECK(cs.is_satisfied({value_type(7)}, {}));

    cs.primary_input_size = 0;
    cs.auxiliary_input_size = 1;
    BOOST_CHECK(cs.is_valid());
    BOOST_CHECK(cs.is_satisfied({}, {value_type(7)}));
}

BOOST_AUTO_TEST_CASE(rejects_out_of_range_terms_in_each_combination) {
    system_type cs;
    cs.primary_input_size = 1;
    cs.auxiliary_input_size = 2;

    const std::array<combination_type constraint_type::*, 3> combinations = {&constraint_type::a, &constraint_type::b,
                                                                             &constraint_type::c};
    for (std::size_t i = 0; i < combinations.size(); ++i) {
        BOOST_TEST_CONTEXT("combination " << i) {
            cs.constraints = {constraint_type()};
            auto &combination = cs.constraints.front().*combinations[i];
            combination.add_term(variable_type(4));
            BOOST_CHECK(!cs.is_valid());

            // An invalid index is still invalid when its coefficient is zero.
            combination.terms.front().coeff = value_type::zero();
            BOOST_CHECK(!cs.is_valid());
        }
    }
}

BOOST_AUTO_TEST_CASE(preserves_linear_combination_ordering_checks) {
    system_type cs;
    cs.primary_input_size = 1;
    cs.auxiliary_input_size = 1;
    cs.constraints = {constraint_type()};
    auto &a = cs.constraints.front().a;

    a.add_term(variable_type(2));
    a.add_term(variable_type(1));
    BOOST_CHECK(!cs.is_valid());

    a.terms.clear();
    a.add_term(variable_type(1), value_type::one());
    a.add_term(variable_type(1), -value_type::one());
    BOOST_CHECK(!cs.is_valid());
}

BOOST_AUTO_TEST_CASE(rejects_dimension_overflow) {
    const auto max_size = std::numeric_limits<std::size_t>::max();
    const std::array<std::pair<std::size_t, std::size_t>, 7> dimensions = {{
        {max_size, 1},
        {1, max_size},
        {max_size, max_size},
        {max_size, 0},
        {0, max_size},
        {1, max_size - 1},
        {max_size - 1, 1},
    }};
    for (const auto &[primary_size, auxiliary_size] : dimensions) {
        BOOST_TEST_CONTEXT("primary size " << primary_size << ", auxiliary size " << auxiliary_size) {
            system_type cs;
            cs.primary_input_size = primary_size;
            cs.auxiliary_input_size = auxiliary_size;
            // Reject overflow even without rows, without allocating an assignment.
            BOOST_CHECK(!cs.is_valid());
        }
    }
}

BOOST_AUTO_TEST_CASE(accepts_largest_representable_index_bound) {
    const auto max_size = std::numeric_limits<std::size_t>::max();
    system_type cs;
    cs.primary_input_size = 1;
    cs.auxiliary_input_size = max_size - 2;
    cs.add_constraint({variable_type(max_size - 1), {}, {}});
    BOOST_CHECK(cs.is_valid());

    cs.constraints.front().a = variable_type(max_size);
    BOOST_CHECK(!cs.is_valid());
}

BOOST_AUTO_TEST_SUITE_END()
