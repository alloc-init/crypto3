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

#define BOOST_TEST_MODULE modified_sap_test

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/operations/lagrange_interpolation.hpp>
#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/modified_sap.hpp>
#include <nil/crypto3/zk/snark/reductions/modified_sap_to_polynomials.hpp>

namespace {
    using field_type = nil::crypto3::algebra::curves::alt_bn128_254::scalar_field_type;
    using value_type = field_type::value_type;
    using system_type = nil::crypto3::zk::snark::modified_sap_constraint_system<field_type>;
    using constraint_type = system_type::constraint_type;
    using variable_type = constraint_type::variable_type;
    using combination_type = constraint_type::linear_combination_type;
    using reduction_type = nil::crypto3::zk::snark::reductions::modified_sap_to_polynomials<field_type>;
    using polynomial_type = reduction_type::polynomial_type;
    using reference_backend_type = nil::crypto3::math::polynomial_arithmetic::schoolbook_backend<value_type>;

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

    system_type basis_evaluation_system() {
        // Slot 3 is unused. Raw terms include duplicates, cancellations and explicit zeros.
        return {4,
                {{raw_combination({{2, 0}, {1, 2}, {1, -2}}), raw_combination({{0, 2}, {2, 5}, {0, -3}, {2, -5}})},
                 {raw_combination({{2, 2}, {0, 3}, {1, 1}, {0, -1}, {2, -1}}),
                  raw_combination({{2, 4}, {0, -1}, {2, -3}, {1, 2}})},
                 {raw_combination({{1, 2}, {0, -1}}), raw_combination({{0, 3}, {2, -1}, {0, -1}})},
                 {},
                 {raw_combination({{2, 3}, {1, -2}, {2, -3}, {0, 0}}), raw_combination({{1, 4}, {0, -1}, {1, -4}})}}};
    }

    system_type witness_interpolation_system() {
        // Satisfied by u = 3, w[1] = 4, w[2] = 13; slot 3 is unused.
        // Logical row evaluations are A = (0, 4, 2, 1, -4), C = (-3, 13, 1, -2, 13).
        return {4,
                {{raw_combination({{3, 0}, {1, 5}, {0, 2}, {1, -5}, {0, -2}}),
                  raw_combination({{1, 4}, {0, 2}, {1, -4}, {0, -3}, {2, 0}})},
                 {raw_combination({{2, 0}, {1, 3}, {0, 5}, {1, -2}, {0, -5}}),
                  raw_combination({{0, 7}, {2, 4}, {0, -7}, {2, -3}})},
                 {raw_combination({{1, -2}, {0, 2}, {1, 1}}), raw_combination({{2, 2}, {1, -3}, {2, -1}, {0, 0}})},
                 {raw_combination({{1, 1}, {0, -1}}), raw_combination({{1, -2}, {0, 2}})},
                 {raw_combination({{1, -1}}), raw_combination({{2, 1}})}}};
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(modified_sap_test_suite)

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
            const std::vector<value_type> witness(invalid_systems[i].witness_size, value_type::one());
            BOOST_CHECK(!invalid_systems[i].is_satisfied(value_type::one(), witness));
            BOOST_CHECK_THROW(reduction_type::get_domain(invalid_systems[i]), std::invalid_argument);
            BOOST_CHECK_THROW(reduction_type::get_padded_constraints(invalid_systems[i]), std::invalid_argument);
            BOOST_CHECK_THROW(reduction_type::instance_map_with_evaluation(invalid_systems[i], value_type::zero()),
                              std::invalid_argument);
            BOOST_CHECK_THROW(reduction_type::witness_map(invalid_systems[i], value_type::one(), witness),
                              std::invalid_argument);
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
                const std::vector<value_type> witness(3, value_type::one());
                BOOST_CHECK(!system.is_satisfied(value_type::one(), witness));
                BOOST_CHECK_THROW(reduction_type::get_domain(system), std::invalid_argument);
                BOOST_CHECK_THROW(reduction_type::get_padded_constraints(system), std::invalid_argument);
                BOOST_CHECK_THROW(reduction_type::instance_map_with_evaluation(system, value_type::zero()),
                                  std::invalid_argument);
                BOOST_CHECK_THROW(reduction_type::witness_map(system, value_type::one(), witness),
                                  std::invalid_argument);
                BOOST_CHECK(system.constraints[1].a.terms == row.a.terms);
                BOOST_CHECK(system.constraints[1].c.terms == row.c.terms);

                auto binding = binding_row();
                auto &terms = (in_a ? binding.a : binding.c).terms;
                terms.insert(terms.end(), invalid_combinations[i].terms.begin(), invalid_combinations[i].terms.end());
                const system_type bad_binding {3, {binding}};
                BOOST_CHECK(!bad_binding.is_valid());
                BOOST_CHECK_THROW(bad_binding.normalized(), std::invalid_argument);
                BOOST_CHECK(!bad_binding.is_satisfied(value_type::one(), witness));
                BOOST_CHECK_THROW(reduction_type::get_domain(bad_binding), std::invalid_argument);
                BOOST_CHECK_THROW(reduction_type::get_padded_constraints(bad_binding), std::invalid_argument);
                BOOST_CHECK_THROW(reduction_type::instance_map_with_evaluation(bad_binding, value_type::zero()),
                                  std::invalid_argument);
                BOOST_CHECK_THROW(reduction_type::witness_map(bad_binding, value_type::one(), witness),
                                  std::invalid_argument);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(satisfaction_checks_every_row_with_explicit_indices) {
    const system_type system {3,
                              {binding_row(),
                               {raw_combination({{1, 1}}), raw_combination({{2, 1}})},
                               {raw_combination({{0, 2}, {1, -1}}), raw_combination({{1, -3}, {2, 1}})}}};
    // 4^2 - 13 = 3, and (2*3 - 4)^2 - (13 - 3*4) = 3.
    const value_type u(3);
    const std::vector<value_type> witness = {u, value_type(4), value_type(13)};
    BOOST_CHECK(system.is_satisfied(u, witness));

    auto changed = witness;
    changed[1] = value_type(5);
    BOOST_CHECK(!system.is_satisfied(u, changed));
    BOOST_CHECK_THROW(reduction_type::witness_map(system, u, changed), std::invalid_argument);
    changed = witness;
    changed[2] += value_type::one();
    BOOST_CHECK(!system.is_satisfied(u, changed));
    BOOST_CHECK_THROW(reduction_type::witness_map(system, u, changed), std::invalid_argument);

    // Negating w[1] preserves the first squaring equation but violates the last row.
    changed = witness;
    changed[1] = -changed[1];
    BOOST_CHECK(!system.is_satisfied(u, changed));
    BOOST_CHECK_THROW(reduction_type::witness_map(system, u, changed), std::invalid_argument);

    auto with_empty_row = system;
    with_empty_row.constraints.push_back({});
    BOOST_CHECK(with_empty_row.is_valid());
    BOOST_CHECK(!with_empty_row.is_satisfied(u, witness));
    BOOST_CHECK_THROW(reduction_type::witness_map(with_empty_row, u, witness), std::invalid_argument);

    // The two-row circuit also accepts field subtraction across zero.
    auto single_equation = system;
    single_equation.constraints.pop_back();
    BOOST_CHECK(single_equation.is_satisfied(u, {u, value_type::zero(), -u}));
}

BOOST_AUTO_TEST_CASE(satisfaction_rejects_wrong_witness_length_or_public_input) {
    const system_type system {3, {binding_row()}};
    const value_type u(3);
    const std::vector<value_type> witness = {u, value_type(4), value_type(13)};
    BOOST_REQUIRE(system.is_satisfied(u, witness));

    const std::vector<std::vector<value_type>> wrong_lengths = {
        {}, {u}, {u, value_type(4)}, {u, value_type(4), value_type(13), value_type::zero()}};
    for (const auto &wrong_length : wrong_lengths) {
        BOOST_TEST_CONTEXT("witness length " << wrong_length.size()) {
            BOOST_CHECK(!system.is_satisfied(u, wrong_length));
            BOOST_CHECK_THROW(reduction_type::witness_map(system, u, wrong_length), std::invalid_argument);
        }
    }

    BOOST_CHECK(!system.is_satisfied(value_type(5), witness));
    BOOST_CHECK_THROW(reduction_type::witness_map(system, value_type(5), witness), std::invalid_argument);
    auto changed = witness;
    changed[0] = value_type(5);
    BOOST_CHECK(!system.is_satisfied(u, changed));
    BOOST_CHECK_THROW(reduction_type::witness_map(system, u, changed), std::invalid_argument);
    BOOST_CHECK(system.is_satisfied(value_type(5), changed));
}

BOOST_AUTO_TEST_CASE(satisfaction_uses_canonical_public_input_parity) {
    const system_type system {1, {binding_row()}};
    // The binding equation holds for every u, isolating the canonical-parity requirement.
    const std::vector<std::pair<value_type, bool>> inputs = {
        {value_type::zero(), false},
        {value_type::one(), true},
        {value_type(2), false},
        {value_type(3), true},
        {-value_type::one(), false},                    // r - 1 is even.
        {-value_type(2), true},                         // r - 2 is odd.
        {value_type(field_type::modulus + 1), true},    // Reduces to 1, despite the even integer input.
        {value_type(field_type::modulus + 2), false}    // Reduces to 2, despite the odd integer input.
    };
    for (std::size_t i = 0; i < inputs.size(); ++i) {
        BOOST_TEST_CONTEXT("public input case " << i) {
            const auto &[u, expected] = inputs[i];
            BOOST_CHECK_EQUAL(system.is_satisfied(u, {u}), expected);
            if (expected) {
                BOOST_CHECK_NO_THROW(reduction_type::witness_map(system, u, {u}));
            } else {
                BOOST_CHECK_THROW(reduction_type::witness_map(system, u, {u}), std::invalid_argument);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(satisfaction_agrees_before_and_after_normalization) {
    const system_type source {
        3,
        {{raw_combination({{2, 0}, {1, 5}, {0, 2}, {1, -5}, {0, -2}}),
          raw_combination({{1, 4}, {0, 2}, {1, -4}, {0, -3}, {2, 0}})},
         {raw_combination({{2, 0}, {1, 3}, {0, 5}, {1, -2}, {0, -5}}),
          raw_combination({{0, 7}, {2, 4}, {0, -7}, {2, -3}})},
         {raw_combination({{1, -2}, {0, 2}, {1, 1}}), raw_combination({{2, 2}, {1, -3}, {2, -1}, {0, 0}})}}};
    const auto original_rows = source.constraints;
    const auto normalized = source.normalized();
    const value_type u(3);
    const std::vector<value_type> witness = {u, value_type(4), value_type(13)};
    BOOST_CHECK(source.is_satisfied(u, witness));
    BOOST_CHECK(normalized.is_satisfied(u, witness));

    auto changed = witness;
    changed[2] += value_type::one();
    BOOST_CHECK(!source.is_satisfied(u, changed));
    BOOST_CHECK(!normalized.is_satisfied(u, changed));
    for (std::size_t i = 0; i < source.constraints.size(); ++i) {
        BOOST_CHECK(source.constraints[i].a.terms == original_rows[i].a.terms);
        BOOST_CHECK(source.constraints[i].c.terms == original_rows[i].c.terms);
    }
}

BOOST_AUTO_TEST_CASE(radix2_domain_selects_minimum_size_and_natural_order) {
    const std::vector<std::pair<std::size_t, std::size_t>> sizes = {{1, 2}, {2, 2}, {3, 4}, {4, 4},
                                                                    {5, 8}, {7, 8}, {8, 8}, {9, 16}};
    for (const auto &[logical_rows, expected_size] : sizes) {
        BOOST_TEST_CONTEXT("logical rows " << logical_rows) {
            const system_type system {1, std::vector<constraint_type>(logical_rows, binding_row())};
            BOOST_CHECK_EQUAL(reduction_type::get_domain_size(logical_rows), expected_size);
            const auto domain = reduction_type::get_domain(system);
            BOOST_REQUIRE_EQUAL(domain->size(), expected_size);
            const auto omega = nil::crypto3::math::unity_root<field_type>(expected_size);
            BOOST_CHECK(domain->get_unity_root() == omega);
            BOOST_CHECK(omega.pow(expected_size / 2) != value_type::one());
            auto point = value_type::one();
            for (std::size_t j = 0; j < expected_size; ++j) {
                BOOST_CHECK(domain->get_domain_element(j) == point);
                point *= omega;
            }
            BOOST_CHECK(point == value_type::one());
        }
    }
}

BOOST_AUTO_TEST_CASE(radix2_domain_rejects_unsupported_counts_without_allocation) {
    constexpr std::size_t max_field_domain = std::size_t(1)
                                             << nil::crypto3::algebra::fields::arithmetic_params<field_type>::s;
    BOOST_CHECK_EQUAL(reduction_type::get_domain_size(max_field_domain - 1), max_field_domain);
    BOOST_CHECK_EQUAL(reduction_type::get_domain_size(max_field_domain), max_field_domain);
    const std::size_t highest_power_of_two = std::size_t(1) << (std::numeric_limits<std::size_t>::digits - 1);
    const std::vector<std::size_t> invalid_counts = {0, max_field_domain + 1, highest_power_of_two,
                                                     highest_power_of_two + 1, std::numeric_limits<std::size_t>::max()};
    for (const auto count : invalid_counts) {
        BOOST_TEST_CONTEXT("logical rows " << count) {
            BOOST_CHECK_THROW(reduction_type::get_domain_size(count), std::invalid_argument);
        }
    }
}

BOOST_AUTO_TEST_CASE(radix2_domain_respects_field_capacity) {
    // BN254's base field supports radix-two size 2 only, allowing a small unsupported-size fixture.
    using base_field_type = nil::crypto3::algebra::curves::alt_bn128_254::base_field_type;
    using base_value_type = base_field_type::value_type;
    using base_reduction_type = nil::crypto3::zk::snark::reductions::modified_sap_to_polynomials<base_field_type>;
    using base_system_type = base_reduction_type::constraint_system_type;
    using base_constraint_type = base_system_type::constraint_type;
    const base_constraint_type binding {{}, base_constraint_type::linear_combination_type(-base_value_type::one())};
    base_system_type system {1, {binding}};
    BOOST_CHECK_EQUAL(base_reduction_type::get_domain(system)->size(), 2);
    system.constraints.resize(3, binding);
    BOOST_REQUIRE(system.is_valid());
    BOOST_CHECK_THROW(base_reduction_type::get_domain(system), std::invalid_argument);
    BOOST_CHECK_THROW(base_reduction_type::get_padded_constraints(system), std::invalid_argument);
    BOOST_CHECK_THROW(base_reduction_type::instance_map_with_evaluation(system, base_value_type::zero()),
                      std::invalid_argument);
    BOOST_CHECK_THROW(base_reduction_type::witness_map(system, base_value_type::one(), {base_value_type::one()}),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(padding_copies_binding_row_and_preserves_logical_rows) {
    const constraint_type raw_binding {raw_combination({{2, 0}, {1, 5}, {1, -5}}),
                                       raw_combination({{1, 4}, {0, 2}, {1, -4}, {0, -3}})};
    for (std::size_t logical_rows : {1, 2, 3, 4, 5}) {
        BOOST_TEST_CONTEXT("logical rows " << logical_rows) {
            system_type source {3, std::vector<constraint_type>(logical_rows, binding_row())};
            source.constraints[0] = raw_binding;
            if (logical_rows > 1) {
                source.constraints[1] = {raw_combination({{1, 3}, {0, 0}, {1, -2}}), raw_combination({{2, 1}})};
            }
            if (logical_rows > 2) {
                source.constraints[2] = {raw_combination({{1, -1}, {0, 2}}), raw_combination({{2, 1}, {1, -3}})};
            }
            const auto original_rows = source.constraints;
            const auto canonical = source.normalized();
            const auto padded = reduction_type::get_padded_constraints(source);
            BOOST_REQUIRE_EQUAL(padded.size(), reduction_type::get_domain_size(logical_rows));
            BOOST_CHECK_EQUAL(source.witness_size, 3);
            BOOST_REQUIRE_EQUAL(source.constraints.size(), logical_rows);
            for (std::size_t j = 0; j < logical_rows; ++j) {
                BOOST_CHECK(padded[j].a.terms == canonical.constraints[j].a.terms);
                BOOST_CHECK(padded[j].c.terms == canonical.constraints[j].c.terms);
                BOOST_CHECK(source.constraints[j].a.terms == original_rows[j].a.terms);
                BOOST_CHECK(source.constraints[j].c.terms == original_rows[j].c.terms);
            }
            const auto binding = binding_row();
            for (std::size_t j = logical_rows; j < padded.size(); ++j) {
                BOOST_CHECK(padded[j].a.terms.empty());
                BOOST_CHECK(padded[j].c.terms == binding.c.terms);
            }
            // Padding must remain satisfied for the nonzero, odd public input u = 3.
            const system_type padded_system {source.witness_size, padded};
            BOOST_CHECK(padded_system.is_satisfied(value_type(3), {value_type(3), value_type(4), value_type(13)}));
        }
    }
}

BOOST_AUTO_TEST_CASE(vanishing_polynomial_matches_domain_and_root_product) {
    for (std::size_t logical_rows : {1, 3, 5}) {
        BOOST_TEST_CONTEXT("logical rows " << logical_rows) {
            const system_type system {1, std::vector<constraint_type>(logical_rows, binding_row())};
            const auto domain = reduction_type::get_domain(system);
            const auto z = domain->get_vanishing_polynomial();
            BOOST_REQUIRE_EQUAL(z.size(), domain->size() + 1);
            BOOST_CHECK_EQUAL(z.degree(), domain->size());
            BOOST_CHECK(z[0] == -value_type::one());
            BOOST_CHECK(z[domain->size()] == value_type::one());
            for (std::size_t i = 1; i < domain->size(); ++i) {
                BOOST_CHECK(z[i].is_zero());
            }
            for (std::size_t j = 0; j < domain->size(); ++j) {
                BOOST_CHECK(z.evaluate(domain->get_domain_element(j)).is_zero());
            }
            for (const value_type &t : {value_type::zero(), value_type::one(), value_type(2), -value_type::one()}) {
                auto product = value_type::one();
                for (std::size_t j = 0; j < domain->size(); ++j) {
                    product *= t - domain->get_domain_element(j);
                }
                BOOST_CHECK(z.evaluate(t) == product);
                BOOST_CHECK(domain->compute_vanishing_polynomial(t) == product);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(basis_evaluation_at_domain_points_recovers_padded_rows) {
    const auto binding = binding_row();
    for (std::size_t logical_rows : {1, 2, 3, 4, 5}) {
        BOOST_TEST_CONTEXT("logical rows " << logical_rows) {
            auto system = basis_evaluation_system();
            system.constraints.resize(logical_rows);
            const auto domain = reduction_type::get_domain(system);
            for (std::size_t j = 0; j < domain->size(); ++j) {
                BOOST_TEST_CONTEXT("domain row " << j) {
                    const auto result =
                        reduction_type::instance_map_with_evaluation(system, domain->get_domain_element(j));
                    BOOST_REQUIRE_EQUAL(result.At.size(), system.num_variables());
                    BOOST_REQUIRE_EQUAL(result.Ct.size(), system.num_variables());
                    BOOST_CHECK(result.Zt.is_zero());
                    const auto &row = j < logical_rows ? system.constraints[j] : binding;
                    for (std::size_t i = 0; i < system.num_variables(); ++i) {
                        // A unit assignment extracts the coefficient at witness slot i.
                        std::vector<value_type> unit(system.num_variables(), value_type::zero());
                        unit[i] = value_type::one();
                        BOOST_CHECK(result.At[i] == row.a.evaluate(unit));
                        BOOST_CHECK(result.Ct[i] == row.c.evaluate(unit));
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(basis_evaluation_matches_coefficient_interpolation) {
    const auto binding = binding_row();
    const std::vector<value_type> points = {value_type::zero(), value_type(2), value_type(42), -value_type(2)};
    for (std::size_t logical_rows : {1, 2, 3, 4, 5}) {
        BOOST_TEST_CONTEXT("logical rows " << logical_rows) {
            auto system = basis_evaluation_system();
            system.constraints.resize(logical_rows);
            const auto domain = reduction_type::get_domain(system);
            std::vector<reduction_type::instance_evaluation> results;
            for (const auto &t : points) {
                BOOST_REQUIRE(!domain->compute_vanishing_polynomial(t).is_zero());
                results.push_back(reduction_type::instance_map_with_evaluation(system, t));
                BOOST_REQUIRE_EQUAL(results.back().At.size(), system.num_variables());
                BOOST_REQUIRE_EQUAL(results.back().Ct.size(), system.num_variables());
                BOOST_CHECK(results.back().Zt == t.pow(domain->size()) - value_type::one());
            }
            for (std::size_t i = 0; i < system.num_variables(); ++i) {
                std::vector<value_type> unit(system.num_variables(), value_type::zero());
                unit[i] = value_type::one();
                std::vector<std::pair<value_type, value_type>> a_points, c_points;
                for (std::size_t j = 0; j < domain->size(); ++j) {
                    const auto &row = j < logical_rows ? system.constraints[j] : binding;
                    const auto x = domain->get_domain_element(j);
                    a_points.emplace_back(x, row.a.evaluate(unit));
                    c_points.emplace_back(x, row.c.evaluate(unit));
                }
                // This oracle constructs coefficient polynomials using explicit products of linear factors.
                const auto a = nil::crypto3::math::lagrange_interpolation(a_points);
                const auto c = nil::crypto3::math::lagrange_interpolation(c_points);
                for (std::size_t k = 0; k < points.size(); ++k) {
                    BOOST_CHECK(results[k].At[i] == a.evaluate(points[k]));
                    BOOST_CHECK(results[k].Ct[i] == c.evaluate(points[k]));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(basis_evaluation_agrees_after_normalization_and_preserves_source) {
    const auto source = basis_evaluation_system();
    const auto original_rows = source.constraints;
    const auto canonical = source.normalized();
    const auto domain = reduction_type::get_domain(source);
    for (const value_type &t : {value_type::zero(), value_type(2), domain->get_domain_element(domain->size() - 1)}) {
        const auto raw_result = reduction_type::instance_map_with_evaluation(source, t);
        const auto canonical_result = reduction_type::instance_map_with_evaluation(canonical, t);
        BOOST_CHECK(raw_result.At == canonical_result.At);
        BOOST_CHECK(raw_result.Ct == canonical_result.Ct);
        BOOST_CHECK(raw_result.Zt == canonical_result.Zt);
    }
    BOOST_REQUIRE_EQUAL(source.constraints.size(), original_rows.size());
    for (std::size_t j = 0; j < source.constraints.size(); ++j) {
        BOOST_CHECK(source.constraints[j].a.terms == original_rows[j].a.terms);
        BOOST_CHECK(source.constraints[j].c.terms == original_rows[j].c.terms);
    }
}

BOOST_AUTO_TEST_CASE(basis_evaluation_of_binding_only_circuit_is_constant) {
    const system_type system {1, {binding_row()}};
    for (const value_type &t : {value_type::zero(), value_type::one(), -value_type::one(), value_type(2)}) {
        const auto result = reduction_type::instance_map_with_evaluation(system, t);
        BOOST_REQUIRE_EQUAL(result.At.size(), 1);
        BOOST_REQUIRE_EQUAL(result.Ct.size(), 1);
        BOOST_CHECK(result.At[0].is_zero());
        BOOST_CHECK(result.Ct[0] == -value_type::one());
        BOOST_CHECK(result.Zt == t.squared() - value_type::one());
    }
}

BOOST_AUTO_TEST_CASE(basis_evaluation_rejects_unrepresentable_witness_size) {
    const std::vector<std::size_t> invalid_sizes = {std::vector<value_type>().max_size() + 1,
                                                    std::numeric_limits<std::size_t>::max()};
    for (const auto n : invalid_sizes) {
        BOOST_TEST_CONTEXT("witness size " << n) {
            const system_type system {n, {binding_row()}};
            BOOST_REQUIRE(system.is_valid());
            BOOST_CHECK_THROW(reduction_type::instance_map_with_evaluation(system, value_type::zero()),
                              std::invalid_argument);
            BOOST_CHECK_THROW(reduction_type::witness_map(system, value_type::one(), {value_type::one()}),
                              std::invalid_argument);
        }
    }
}

BOOST_AUTO_TEST_CASE(witness_interpolation_matches_reference_and_padded_rows) {
    const value_type u(3);
    const std::vector<value_type> witness = {u, value_type(4), value_type(13), value_type(42)};
    for (std::size_t logical_rows : {1, 2, 3, 4, 5, 7, 8, 9}) {
        BOOST_TEST_CONTEXT("logical rows " << logical_rows) {
            auto system = witness_interpolation_system();
            const auto repeated_row = system.constraints[1];
            system.constraints.resize(logical_rows, repeated_row);
            const auto domain = reduction_type::get_domain(system);
            const auto padded = reduction_type::get_padded_constraints(system);
            const auto result = reduction_type::witness_map(system, u, witness);
            BOOST_REQUIRE_EQUAL(result.A.size(), domain->size());
            BOOST_REQUIRE_EQUAL(result.C.size(), domain->size());

            std::vector<std::pair<value_type, value_type>> a_points, c_points;
            for (std::size_t j = 0; j < domain->size(); ++j) {
                const auto x = domain->get_domain_element(j);
                const auto a = padded[j].a.evaluate(witness);
                const auto c = padded[j].c.evaluate(witness);
                BOOST_CHECK(result.A.evaluate(x) == a);
                BOOST_CHECK(result.C.evaluate(x) == c);
                BOOST_CHECK(a.squared() - c == u);
                a_points.emplace_back(x, a);
                c_points.emplace_back(x, c);
            }
            // Independent coefficient interpolation uses products of linear factors, without the inverse FFT.
            auto expected_a = nil::crypto3::math::lagrange_interpolation(a_points);
            auto expected_c = nil::crypto3::math::lagrange_interpolation(c_points);
            expected_a.resize(domain->size(), value_type::zero());
            expected_c.resize(domain->size(), value_type::zero());
            BOOST_CHECK(result.A == expected_a);
            BOOST_CHECK(result.C == expected_c);

            // Compute the quotient and full polynomial identity using coefficient arithmetic, without coset FFTs.
            const reference_backend_type reference;
            polynomial_type numerator, expected_h, remainder, zh;
            reference.square(numerator, expected_a);
            numerator -= expected_c;
            numerator[0] -= u;
            numerator.condense();
            const auto z = domain->get_vanishing_polynomial();
            nil::crypto3::math::division(expected_h, remainder, numerator, z);
            BOOST_REQUIRE(remainder.is_zero());
            BOOST_REQUIRE_LE(expected_h.size(), domain->size() - 1);
            expected_h.resize(domain->size() - 1, value_type::zero());
            BOOST_REQUIRE_EQUAL(result.H.size(), domain->size() - 1);
            BOOST_CHECK(result.H == expected_h);
            reference.multiply(zh, z, result.H);
            BOOST_CHECK(numerator == zh);
        }
    }
}

BOOST_AUTO_TEST_CASE(witness_interpolation_matches_weighted_basis_evaluations) {
    const auto system = witness_interpolation_system();
    const value_type u(3);
    const std::vector<value_type> witness = {u, value_type(4), value_type(13), value_type(42)};
    const auto domain = reduction_type::get_domain(system);
    const auto result = reduction_type::witness_map(system, u, witness);
    for (const value_type &t :
         {value_type::zero(), value_type(2), -value_type(2), domain->get_domain_element(domain->size() - 1)}) {
        const auto basis = reduction_type::instance_map_with_evaluation(system, t);
        auto a = value_type::zero();
        auto c = value_type::zero();
        for (std::size_t i = 0; i < witness.size(); ++i) {
            a += witness[i] * basis.At[i];
            c += witness[i] * basis.Ct[i];
        }
        BOOST_CHECK(result.A.evaluate(t) == a);
        BOOST_CHECK(result.C.evaluate(t) == c);
        BOOST_CHECK(a.squared() - c - u == basis.Zt * result.H.evaluate(t));
    }
}

BOOST_AUTO_TEST_CASE(witness_interpolation_agrees_after_normalization_and_preserves_inputs) {
    const auto source = witness_interpolation_system();
    const auto original_rows = source.constraints;
    const value_type u(3);
    const std::vector<value_type> witness = {u, value_type(4), value_type(13), value_type(42)};
    const auto original_witness = witness;
    const auto raw_result = reduction_type::witness_map(source, u, witness);
    const auto canonical_result = reduction_type::witness_map(source.normalized(), u, witness);
    BOOST_CHECK(raw_result.A == canonical_result.A);
    BOOST_CHECK(raw_result.C == canonical_result.C);
    BOOST_CHECK(raw_result.H == canonical_result.H);

    auto changed = witness;
    changed[3] += value_type::one();
    const auto unused_slot_result = reduction_type::witness_map(source, u, changed);
    BOOST_CHECK(raw_result.A == unused_slot_result.A);
    BOOST_CHECK(raw_result.C == unused_slot_result.C);
    BOOST_CHECK(raw_result.H == unused_slot_result.H);
    BOOST_CHECK(witness == original_witness);
    BOOST_CHECK_EQUAL(source.witness_size, witness.size());
    BOOST_REQUIRE_EQUAL(source.constraints.size(), original_rows.size());
    for (std::size_t j = 0; j < source.constraints.size(); ++j) {
        BOOST_CHECK(source.constraints[j].a.terms == original_rows[j].a.terms);
        BOOST_CHECK(source.constraints[j].c.terms == original_rows[j].c.terms);
    }
}

BOOST_AUTO_TEST_CASE(witness_interpolation_of_binding_only_circuit_is_constant) {
    for (std::size_t logical_rows : {1, 2, 3, 5}) {
        BOOST_TEST_CONTEXT("logical rows " << logical_rows) {
            const system_type system {1, std::vector<constraint_type>(logical_rows, binding_row())};
            for (const value_type &u : {value_type::one(), value_type(3), -value_type(2)}) {
                const auto result = reduction_type::witness_map(system, u, {u});
                const auto m = reduction_type::get_domain_size(logical_rows);
                BOOST_REQUIRE_EQUAL(result.A.size(), m);
                BOOST_REQUIRE_EQUAL(result.C.size(), m);
                BOOST_REQUIRE_EQUAL(result.H.size(), m - 1);
                BOOST_CHECK(result.H.is_zero());
                for (std::size_t i = 0; i < m; ++i) {
                    BOOST_CHECK(result.A[i].is_zero());
                    BOOST_CHECK(result.C[i] == (i == 0 ? -u : value_type::zero()));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(quotient_on_minimum_domain_is_constant) {
    const system_type system {3, {binding_row(), {raw_combination({{1, 1}}), raw_combination({{2, 1}})}}};
    for (const value_type &u : {value_type::one(), value_type(3), -value_type(2)}) {
        const std::vector<value_type> witness = {u, value_type(4), value_type(16) - u};
        const auto result = reduction_type::witness_map(system, u, witness);
        // A = 2 - 2X, C = (8 - u) - 8X, so A^2 - C - u = 4(X^2 - 1).
        const polynomial_type expected_a {value_type(2), -value_type(2)};
        const polynomial_type expected_c {value_type(8) - u, -value_type(8)};
        const polynomial_type expected_h {value_type(4)};
        BOOST_CHECK(result.A == expected_a);
        BOOST_CHECK(result.C == expected_c);
        BOOST_CHECK(result.H == expected_h);
    }
}

BOOST_AUTO_TEST_CASE(quotient_reaches_degree_bound) {
    const auto u = value_type::one();
    for (std::size_t m : {2, 4, 8, 16}) {
        BOOST_TEST_CONTEXT("domain size " << m) {
            system_type system {1, std::vector<constraint_type>(m, binding_row())};
            const auto domain = reduction_type::get_domain(system);
            // Prescribe A = X^(m-1) - 1. Since u = w[0] = 1, row coefficients are row values.
            for (std::size_t j = 1; j < m; ++j) {
                const auto a = domain->get_domain_element(j).pow(m - 1) - value_type::one();
                system.constraints[j] = {combination_type(a), combination_type(a.squared() - u)};
            }
            const auto result = reduction_type::witness_map(system, u, {u});
            polynomial_type expected_a(m, value_type::zero());
            expected_a[0] = -value_type::one();
            expected_a[m - 1] = value_type::one();
            polynomial_type expected_c(m, value_type::zero());
            expected_c[m - 2] = value_type::one();
            expected_c[m - 1] = -value_type(2);
            // A^2 - C - 1 = X^(2m-2) - X^(m-2) = (X^m - 1) * X^(m-2).
            polynomial_type expected_h(m - 1, value_type::zero());
            expected_h[m - 2] = value_type::one();
            BOOST_CHECK(result.A == expected_a);
            BOOST_CHECK(result.C == expected_c);
            BOOST_CHECK(result.H == expected_h);
        }
    }
}

BOOST_AUTO_TEST_CASE(quotient_rejects_nonzero_remainder) {
    const auto system = witness_interpolation_system();
    const auto domain = reduction_type::get_domain(system);
    const auto padded = reduction_type::get_padded_constraints(system);
    const value_type u(3);
    // This sign change preserves the first squaring equation but violates later rows.
    const std::vector<value_type> witness = {u, -value_type(4), value_type(13), value_type(42)};
    std::vector<std::pair<value_type, value_type>> a_points, c_points;
    for (std::size_t j = 0; j < domain->size(); ++j) {
        const auto x = domain->get_domain_element(j);
        a_points.emplace_back(x, padded[j].a.evaluate(witness));
        c_points.emplace_back(x, padded[j].c.evaluate(witness));
    }
    const auto a = nil::crypto3::math::lagrange_interpolation(a_points);
    const auto c = nil::crypto3::math::lagrange_interpolation(c_points);
    polynomial_type numerator, quotient, remainder;
    reference_backend_type().square(numerator, a);
    numerator -= c;
    numerator[0] -= u;
    numerator.condense();
    nil::crypto3::math::division(quotient, remainder, numerator, domain->get_vanishing_polynomial());
    BOOST_REQUIRE(!remainder.is_zero());
    BOOST_CHECK_THROW(reduction_type::witness_map(system, u, witness), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
