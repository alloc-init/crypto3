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
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/modified_sqap.hpp>
#include <nil/crypto3/zk/snark/reductions/modified_sqap_to_polynomials.hpp>

namespace {
    using field_type = nil::crypto3::algebra::curves::alt_bn128_254::scalar_field_type;
    using value_type = field_type::value_type;
    using system_type = nil::crypto3::zk::snark::modified_sqap_constraint_system<field_type>;
    using constraint_type = system_type::constraint_type;
    using variable_type = constraint_type::variable_type;
    using combination_type = constraint_type::linear_combination_type;
    using reduction_type = nil::crypto3::zk::snark::reductions::modified_sqap_to_polynomials<field_type>;

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
            const std::vector<value_type> witness(invalid_systems[i].witness_size, value_type::one());
            BOOST_CHECK(!invalid_systems[i].is_satisfied(value_type::one(), witness));
            BOOST_CHECK_THROW(reduction_type::get_domain(invalid_systems[i]), std::invalid_argument);
            BOOST_CHECK_THROW(reduction_type::get_padded_constraints(invalid_systems[i]), std::invalid_argument);
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
    changed = witness;
    changed[2] += value_type::one();
    BOOST_CHECK(!system.is_satisfied(u, changed));

    // Negating w[1] preserves the first squaring equation but violates the last row.
    changed = witness;
    changed[1] = -changed[1];
    BOOST_CHECK(!system.is_satisfied(u, changed));

    auto with_empty_row = system;
    with_empty_row.constraints.push_back({});
    BOOST_CHECK(with_empty_row.is_valid());
    BOOST_CHECK(!with_empty_row.is_satisfied(u, witness));

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
        }
    }

    BOOST_CHECK(!system.is_satisfied(value_type(5), witness));
    auto changed = witness;
    changed[0] = value_type(5);
    BOOST_CHECK(!system.is_satisfied(u, changed));
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
    using base_reduction_type = nil::crypto3::zk::snark::reductions::modified_sqap_to_polynomials<base_field_type>;
    using base_system_type = base_reduction_type::constraint_system_type;
    using base_constraint_type = base_system_type::constraint_type;
    const base_constraint_type binding {{}, base_constraint_type::linear_combination_type(-base_value_type::one())};
    base_system_type system {1, {binding}};
    BOOST_CHECK_EQUAL(base_reduction_type::get_domain(system)->size(), 2);
    system.constraints.resize(3, binding);
    BOOST_REQUIRE(system.is_valid());
    BOOST_CHECK_THROW(base_reduction_type::get_domain(system), std::invalid_argument);
    BOOST_CHECK_THROW(base_reduction_type::get_padded_constraints(system), std::invalid_argument);
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

BOOST_AUTO_TEST_SUITE_END()
