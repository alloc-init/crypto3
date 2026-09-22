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

#define BOOST_TEST_MODULE modified_sqap_test

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/modified_sqap.hpp>

namespace {
    using field_type = nil::crypto3::algebra::curves::alt_bn128_254::scalar_field_type;
    using value_type = field_type::value_type;
    using system_type = nil::crypto3::zk::snark::modified_sqap_constraint_system<field_type>;
    using constraint_type = system_type::constraint_type;
    using variable_type = constraint_type::variable_type;
    using combination_type = constraint_type::linear_combination_type;

    combination_type raw_combination(std::initializer_list<std::pair<std::size_t, int>> terms) {
        combination_type result;
        // Preserve the supplied order and duplicates to exercise normalization.
        for (const auto &[index, coefficient] : terms) {
            result.add_term(variable_type(index), value_type(coefficient));
        }
        return result;
    }

    constraint_type binding_row() {
        return {{}, combination_type(-variable_type(0))};
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(modified_sqap_test_suite)

BOOST_AUTO_TEST_CASE(valid_structure_and_index_boundaries) {
    const system_type minimum {1, {binding_row()}};
    BOOST_CHECK(minimum.is_valid());
    BOOST_CHECK(minimum.normalized() == minimum);

    system_type system {3, {binding_row(), {raw_combination({{2, 1}}), raw_combination({{0, 2}})}, {}}};
    // Structural validity does not require a witness or exclude unsatisfied rows.
    BOOST_CHECK(system.is_valid());
    const auto normalized = system.normalized();
    BOOST_CHECK_EQUAL(normalized.num_variables(), 3);
    BOOST_CHECK_EQUAL(normalized.num_constraints(), 3);
    BOOST_CHECK(normalized.constraints.back().a.terms.empty());
    BOOST_CHECK(normalized.constraints.back().c.terms.empty());

    // Declared witness entries need not all occur in a constraint.
    system.witness_size = 4;
    BOOST_CHECK(system.is_valid());
    BOOST_CHECK_EQUAL(system.normalized().num_variables(), 4);
}

BOOST_AUTO_TEST_CASE(binding_row_is_checked_after_normalization) {
    const constraint_type binding {raw_combination({{2, 0}, {1, 5}, {0, 2}, {1, -5}, {0, -2}}),
                                   raw_combination({{1, 4}, {0, 2}, {1, -4}, {0, -3}, {2, 0}})};
    const system_type system {3, {binding}};
    BOOST_CHECK(system.is_valid());
    const auto normalized = system.normalized();
    BOOST_REQUIRE_EQUAL(normalized.constraints.size(), 1);
    BOOST_CHECK(normalized.constraints.front().a.terms.empty());
    BOOST_REQUIRE_EQUAL(normalized.constraints.front().c.terms.size(), 1);
    BOOST_CHECK_EQUAL(normalized.constraints.front().c.terms.front().index, 0);
    BOOST_CHECK(normalized.constraints.front().c.terms.front().coeff == -value_type::one());
    BOOST_CHECK(system.constraints.front().a.terms == binding.a.terms);
    BOOST_CHECK(system.constraints.front().c.terms == binding.c.terms);
}

BOOST_AUTO_TEST_CASE(normalization_preserves_evaluations_order_and_source) {
    const system_type source {
        3,
        {binding_row(),
         {raw_combination({{2, 3}, {0, 2}, {1, 0}, {2, -1}, {0, -2}, {1, 4}}),
          raw_combination({{2, 7}, {0, -3}, {2, -7}, {1, 0}})},
         {raw_combination({{1, 9}, {0, 1}, {1, -9}}), raw_combination({{2, -2}, {1, 6}, {2, 5}})}}};
    const auto original_rows = source.constraints;
    auto normalized = source.normalized();
    BOOST_REQUIRE_EQUAL(normalized.constraints.size(), original_rows.size());
    BOOST_CHECK_EQUAL(normalized.witness_size, source.witness_size);
    BOOST_CHECK(normalized.is_valid());
    BOOST_CHECK(normalized.constraints[1].a.terms == raw_combination({{1, 4}, {2, 2}}).terms);
    BOOST_CHECK(normalized.constraints[1].c.terms == raw_combination({{0, -3}}).terms);
    BOOST_CHECK(normalized.constraints[2].a.terms == raw_combination({{0, 1}}).terms);
    BOOST_CHECK(normalized.constraints[2].c.terms == raw_combination({{1, 6}, {2, 3}}).terms);

    const std::vector<std::vector<value_type>> assignments = {
        {value_type(3), value_type(5), value_type(7)},
        {value_type(11), -value_type(2), value_type::zero()},
        {value_type::zero(), value_type::one(), -value_type::one()}};
    for (std::size_t i = 0; i < source.constraints.size(); ++i) {
        BOOST_TEST_CONTEXT("row " << i) {
            BOOST_CHECK(source.constraints[i].a.terms == original_rows[i].a.terms);
            BOOST_CHECK(source.constraints[i].c.terms == original_rows[i].c.terms);
            for (const auto &assignment : assignments) {
                BOOST_CHECK(normalized.constraints[i].a.evaluate(assignment) ==
                            source.constraints[i].a.evaluate(assignment));
                BOOST_CHECK(normalized.constraints[i].c.evaluate(assignment) ==
                            source.constraints[i].c.evaluate(assignment));
            }
        }
    }

    const auto again = normalized.normalized();
    BOOST_CHECK_EQUAL(again.witness_size, normalized.witness_size);
    BOOST_REQUIRE_EQUAL(again.constraints.size(), normalized.constraints.size());
    for (std::size_t i = 0; i < again.constraints.size(); ++i) {
        // Compare the stored sequences directly: linear-combination equality sorts its operands.
        BOOST_CHECK(again.constraints[i].a.terms == normalized.constraints[i].a.terms);
        BOOST_CHECK(again.constraints[i].c.terms == normalized.constraints[i].c.terms);
    }
    normalized.constraints[1].a.terms.front().coeff += value_type::one();
    BOOST_CHECK(source.constraints[1].a.terms == original_rows[1].a.terms);
}

BOOST_AUTO_TEST_CASE(rejects_missing_dimensions_or_incorrect_binding_row) {
    const std::vector<system_type> invalid_systems = {{},
                                                      {0, {binding_row()}},
                                                      {1, {}},
                                                      {2, {{}}},
                                                      {2, {{raw_combination({{0, 1}}), raw_combination({{0, -1}})}}},
                                                      {2, {{{}, raw_combination({{0, 1}})}}},
                                                      {2, {{{}, raw_combination({{0, -2}})}}},
                                                      {2, {{{}, raw_combination({{1, -1}})}}},
                                                      {2, {{{}, raw_combination({{0, -1}, {1, 1}})}}},
                                                      {2, {{}, binding_row()}}};
    for (std::size_t i = 0; i < invalid_systems.size(); ++i) {
        BOOST_TEST_CONTEXT("invalid structure " << i) {
            BOOST_CHECK(!invalid_systems[i].is_valid());
            BOOST_CHECK_THROW(invalid_systems[i].normalized(), std::invalid_argument);
        }
    }
}

BOOST_AUTO_TEST_CASE(rejects_out_of_range_indices_before_terms_can_disappear) {
    const std::vector<combination_type> invalid_combinations = {
        raw_combination({{3, 1}}), raw_combination({{3, 0}}), raw_combination({{3, 5}, {3, -5}}),
        raw_combination({{std::numeric_limits<std::size_t>::max(), 0}})};
    for (std::size_t i = 0; i < invalid_combinations.size(); ++i) {
        for (bool in_a : {true, false}) {
            BOOST_TEST_CONTEXT("invalid combination " << i << ", in A: " << in_a) {
                constraint_type row;
                (in_a ? row.a : row.c) = invalid_combinations[i];
                const system_type system {3, {binding_row(), row}};
                BOOST_CHECK(!system.is_valid());
                BOOST_CHECK_THROW(system.normalized(), std::invalid_argument);
                BOOST_CHECK(system.constraints[1].a.terms == row.a.terms);
                BOOST_CHECK(system.constraints[1].c.terms == row.c.terms);

                auto binding = binding_row();
                auto &terms = (in_a ? binding.a : binding.c).terms;
                terms.insert(terms.end(), invalid_combinations[i].terms.begin(), invalid_combinations[i].terms.end());
                const system_type bad_binding {3, {binding}};
                BOOST_CHECK(!bad_binding.is_valid());
                BOOST_CHECK_THROW(bad_binding.normalized(), std::invalid_argument);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
