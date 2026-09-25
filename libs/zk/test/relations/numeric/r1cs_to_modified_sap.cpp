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

#define BOOST_TEST_MODULE r1cs_to_modified_sap_test

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/alt_bn128/scalar_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/zk/snark/reductions/r1cs_to_modified_sap.hpp>

namespace {
    namespace fields = nil::crypto3::algebra::fields;
    namespace snark = nil::crypto3::zk::snark;
    using field_type = fields::alt_bn128_scalar_field<254>;
    using value_type = field_type::value_type;
    using recovery_type = snark::reductions::detail::modified_sap_constant_recovery<field_type>;
    using reduction_type = snark::reductions::r1cs_to_modified_sap<field_type>;
    using source_system_type = snark::r1cs_constraint_system<field_type>;
    using source_constraint_type = snark::r1cs_constraint<field_type>;
    using primary_input_type = snark::r1cs_primary_input<field_type>;
    using auxiliary_input_type = snark::r1cs_auxiliary_input<field_type>;
    using variable_type = recovery_type::variable_type;

    source_system_type source_example() {
        source_system_type cs;
        cs.primary_input_size = 1;
        cs.auxiliary_input_size = 2;
        // Source indices are 0 = one, 1 = u, 2 = x, 3 = y; enforce (x + 1)*y = u.
        source_constraint_type row {variable_type(0), variable_type(3), variable_type(1)};
        row.a.add_term(variable_type(2));
        cs.add_constraint(row);
        return cs;
    }

    recovery_type::combination_type target_combination(std::initializer_list<std::pair<std::size_t, int>> terms) {
        recovery_type::combination_type result;
        for (const auto &[index, coefficient] : terms) {
            result.add_term(variable_type(index), value_type(coefficient));
        }
        return result;
    }

    // Use Crypto3's prime-field arithmetic with small moduli for exhaustive checks.
    template<unsigned Modulus>
    class test_prime_field : public fields::field<8> {
    public:
        using integral_type = fields::field<8>::integral_type;
        constexpr static integral_type modulus = Modulus;
        constexpr static integral_type group_order_minus_one_half = (modulus - 1u) / 2;
        constexpr static modular_params_type modulus_params = modulus.backend();
        using modular_type = boost::multiprecision::number<boost::multiprecision::backends::modular_adaptor<
            modular_backend, boost::multiprecision::backends::modular_params_ct<modular_backend, modulus_params>>>;
        using value_type = fields::detail::element_fp<fields::params<test_prime_field>>;
    };

    template<typename FieldType>
    snark::modified_sap_constraint_system<FieldType> recovery_system(std::size_t constant_index) {
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<FieldType>;
        using variable = typename recovery::variable_type;
        snark::modified_sap_constraint_system<FieldType> cs;
        cs.witness_size = constant_index + recovery::witness_size;
        cs.constraints.push_back({{}, -variable(0)});
        recovery::append_constraints(cs.constraints, constant_index);
        return cs;
    }

    template<typename FieldType>
    bool rows_satisfied(const snark::modified_sap_constraint_system<FieldType> &cs,
                        const typename FieldType::value_type &u,
                        const std::vector<typename FieldType::value_type> &witness, std::size_t begin,
                        std::size_t end) {
        for (std::size_t row = begin; row < end; ++row) {
            if (cs.constraints[row].a.evaluate(witness).squared() - cs.constraints[row].c.evaluate(witness) != u) {
                return false;
            }
        }
        return true;
    }

    template<unsigned Modulus>
    void test_all_recovery_assignments() {
        using field = test_prime_field<Modulus>;
        using value = typename field::value_type;
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<field>;
        const auto cs = recovery_system<field>(1);
        BOOST_REQUIRE(cs.is_valid());
        for (unsigned u = 0; u < Modulus; ++u) {
            const auto honest_auxiliary = (value::one() - value(u)).squared();
            if ((u & 1) != 0) {
                std::vector<value> generated = {value(u)};
                recovery::append_witness(generated);
                const std::vector<value> expected = {value(u), value::one(), honest_auxiliary};
                BOOST_CHECK(generated == expected);
            }
            for (unsigned e = 0; e < Modulus; ++e) {
                for (unsigned t = 0; t < Modulus; ++t) {
                    BOOST_TEST_CONTEXT("modulus " << Modulus << ", u " << u << ", e " << e << ", t " << t) {
                        const std::vector<value> witness = {value(u), value(e), value(t)};
                        // At u = 0, the rows permit any e with t = e^2; the public-input rule excludes this case.
                        const bool valid_rows =
                            u == 0 ? value(t) == value(e).squared() : e == 1 && value(t) == honest_auxiliary;
                        BOOST_CHECK_EQUAL(rows_satisfied(cs, value(u), witness, 0, cs.constraints.size()), valid_rows);
                        BOOST_CHECK_EQUAL(cs.is_satisfied(value(u), witness), valid_rows && (u & 1) != 0);
                    }
                }
            }
        }
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(r1cs_to_modified_sap_test_suite)

BOOST_AUTO_TEST_CASE(recovery_dimensions_and_source_prefix) {
    BOOST_CHECK_EQUAL(recovery_type::witness_size, 2);
    BOOST_CHECK_EQUAL(recovery_type::constraint_count, 2);
    for (const std::size_t constant_index : {1, 3}) {
        const auto cs = recovery_system<field_type>(constant_index);
        BOOST_CHECK_EQUAL(cs.num_variables(), constant_index + 2);
        BOOST_CHECK_EQUAL(cs.num_constraints(), 3);
        BOOST_REQUIRE(cs.is_valid());
        const auto normalized = cs.normalized();
        const std::vector<recovery_type::constraint_type> expected_rows = {
            {{}, target_combination({{0, -1}})},
            {target_combination({{0, 1}, {constant_index, 1}}), target_combination({{0, 3}, {constant_index + 1, 1}})},
            {target_combination({{0, -1}, {constant_index, 1}}),
             target_combination({{0, -1}, {constant_index + 1, 1}})},
        };
        BOOST_CHECK(normalized.constraints == expected_rows);
        for (const auto &u : {value_type(1), value_type(3), value_type(field_type::modulus - 2)}) {
            std::vector<value_type> witness(constant_index, value_type(13));
            witness[0] = u;
            const auto prefix = witness;
            recovery_type::append_witness(witness);
            auto expected = prefix;
            expected.push_back(value_type::one());
            expected.push_back((value_type::one() - u).squared());
            BOOST_CHECK(witness == expected);
            BOOST_CHECK_EQUAL_COLLECTIONS(witness.begin(), witness.begin() + constant_index, prefix.begin(),
                                          prefix.end());
            BOOST_CHECK_EQUAL(witness.size(), cs.num_variables());
            BOOST_CHECK(cs.is_satisfied(u, witness));
            BOOST_CHECK(normalized.is_satisfied(u, witness));
        }
    }
}

BOOST_AUTO_TEST_CASE(rejects_forged_constant_even_with_adapted_auxiliary) {
    const auto cs = recovery_system<field_type>(1);
    for (const auto &u : {value_type(1), value_type(3), value_type(field_type::modulus - 2)}) {
        for (const auto &e : {value_type(0), value_type(2), -value_type::one()}) {
            // Satisfy the second row with a forged e and its matching auxiliary: the first row must reject.
            std::vector<value_type> witness = {u, e, (e - u).squared()};
            BOOST_CHECK(rows_satisfied(cs, u, witness, 2, 3));
            BOOST_CHECK(!rows_satisfied(cs, u, witness, 1, 2));
            BOOST_CHECK(!cs.is_satisfied(u, witness));
            // Satisfying the first row instead must leave the second row unsatisfied.
            witness[2] = (e + u).squared() - value_type(4) * u;
            BOOST_CHECK(rows_satisfied(cs, u, witness, 1, 2));
            BOOST_CHECK(!rows_satisfied(cs, u, witness, 2, 3));
            BOOST_CHECK(!cs.is_satisfied(u, witness));
        }
    }
}

BOOST_AUTO_TEST_CASE(rejects_changed_recovery_auxiliary_and_public_input) {
    const auto cs = recovery_system<field_type>(1);
    const value_type u(field_type::modulus - 2);
    std::vector<value_type> witness = {u};
    recovery_type::append_witness(witness);
    BOOST_REQUIRE(cs.is_satisfied(u, witness));

    auto changed = witness;
    changed[2] += value_type::one();
    BOOST_CHECK(!cs.is_satisfied(u, changed));
    changed = witness;
    changed[0] += value_type::one();
    BOOST_CHECK(!rows_satisfied(cs, u, changed, 0, 1));
    BOOST_CHECK(!cs.is_satisfied(u, changed));
    BOOST_CHECK(!cs.is_satisfied(value_type(1), witness));
}

BOOST_AUTO_TEST_CASE(zero_input_does_not_recover_one_and_is_rejected) {
    const auto cs = recovery_system<field_type>(1);
    for (const auto &e : {value_type(0), value_type(1), value_type(2)}) {
        const std::vector<value_type> witness = {value_type::zero(), e, e.squared()};
        BOOST_CHECK(rows_satisfied(cs, value_type::zero(), witness, 0, cs.constraints.size()));
        BOOST_CHECK(!cs.is_satisfied(value_type::zero(), witness));
    }
}

BOOST_AUTO_TEST_CASE(invalid_inputs_fail_before_mutation) {
    for (const auto &u : {value_type(0), value_type(2), value_type(field_type::modulus - 1)}) {
        std::vector<value_type> witness = {u, value_type(13)};
        const auto original = witness;
        BOOST_CHECK_THROW(recovery_type::append_witness(witness), std::invalid_argument);
        BOOST_CHECK(witness == original);
    }
    std::vector<value_type> empty;
    BOOST_CHECK_THROW(recovery_type::append_witness(empty), std::invalid_argument);
    BOOST_CHECK(empty.empty());

    auto cs = recovery_system<field_type>(1);
    const auto original = cs.constraints;
    const auto max_size = std::numeric_limits<std::size_t>::max();
    for (const auto constant_index : {std::size_t(0), max_size, max_size - recovery_type::witness_size + 1,
                                      std::vector<value_type>().max_size() - recovery_type::witness_size + 1}) {
        BOOST_CHECK_THROW(recovery_type::append_constraints(cs.constraints, constant_index), std::invalid_argument);
        BOOST_CHECK(cs.constraints == original);
    }
}

BOOST_AUTO_TEST_CASE(exhaustive_small_field_recovery_assignments) {
    test_all_recovery_assignments<3>();
    test_all_recovery_assignments<5>();
    test_all_recovery_assignments<7>();
    test_all_recovery_assignments<11>();
    test_all_recovery_assignments<13>();
    test_all_recovery_assignments<17>();
    test_all_recovery_assignments<31>();
}

BOOST_AUTO_TEST_CASE(instance_map_matches_fixed_rows_and_supports_multiple_public_inputs) {
    const auto source = source_example();
    const auto original = source;
    BOOST_REQUIRE(source.is_valid());
    const auto converted = reduction_type::instance_map(source);
    BOOST_CHECK(source == original);
    BOOST_REQUIRE(converted.is_valid());
    BOOST_CHECK_EQUAL(converted.num_variables(), 6);
    BOOST_CHECK_EQUAL(converted.num_constraints(), 5);

    // N = 3. The recovered one is w[3], and the source product's auxiliary is w[5].
    auto expected = recovery_system<field_type>(3).normalized();
    ++expected.witness_size;
    expected.constraints.push_back(
        {target_combination({{1, 1}, {2, 1}, {3, 1}}), target_combination({{0, 3}, {5, 1}})});
    expected.constraints.push_back(
        {target_combination({{1, 1}, {2, -1}, {3, 1}}), target_combination({{0, -1}, {5, 1}})});
    BOOST_CHECK(converted == expected);
    BOOST_CHECK(converted.normalized() == converted);

    // The same circuit accepts several odd public inputs; x = u - 1 and y = 1.
    for (const auto &u : {value_type(3), value_type(5), value_type(field_type::modulus - 2)}) {
        const value_type x = u - value_type::one();
        BOOST_REQUIRE(source.is_satisfied({u}, {x, value_type::one()}));
        std::vector<value_type> witness = {u, x, value_type::one()};
        recovery_type::append_witness(witness);
        witness.push_back(x.squared());    // ((x + 1) - y)^2 = x^2.
        BOOST_CHECK(converted.is_satisfied(u, witness));
        witness.back() += value_type::one();
        BOOST_CHECK(!converted.is_satisfied(u, witness));
    }
}

BOOST_AUTO_TEST_CASE(instance_map_normalizes_raw_terms_without_changing_the_source) {
    auto source = source_example();
    auto &row = source.constraints.front();
    row.a.terms.clear();
    row.a.add_term(variable_type(2), value_type(3));
    row.a.add_term(variable_type(0), value_type(2));
    row.a.add_term(variable_type(2), -value_type(2));
    row.a.add_term(variable_type(0), -value_type::one());
    row.a.add_term(variable_type(1), value_type::zero());
    row.b.terms.clear();
    row.b.add_term(variable_type(3), value_type(2));
    row.b.add_term(variable_type(2), value_type(7));
    row.b.add_term(variable_type(3), -value_type::one());
    row.b.add_term(variable_type(2), -value_type(7));
    row.b.add_term(variable_type(0), value_type::zero());
    row.c.terms.clear();
    row.c.add_term(variable_type(1), value_type(4));
    row.c.add_term(variable_type(0), value_type(9));
    row.c.add_term(variable_type(1), -value_type(3));
    row.c.add_term(variable_type(0), -value_type(9));

    BOOST_CHECK(!source.is_valid());    // Ordinary validation requires sorted, unique indices.
    const auto original = source;
    BOOST_CHECK(reduction_type::instance_map(source) == reduction_type::instance_map(source_example()));
    BOOST_CHECK(source == original);
}

BOOST_AUTO_TEST_CASE(instance_map_empty_sources_keep_binding_recovery_and_unused_variables) {
    for (const std::size_t auxiliary_size : {0, 3}) {
        source_system_type source;
        source.primary_input_size = 1;
        source.auxiliary_input_size = auxiliary_size;
        const auto converted = reduction_type::instance_map(source);
        BOOST_CHECK(converted == recovery_system<field_type>(source.num_variables()).normalized());
        BOOST_CHECK_EQUAL(converted.num_variables(), source.num_variables() + 2);
        BOOST_CHECK_EQUAL(converted.num_constraints(), 3);
        std::vector<value_type> witness(source.num_variables(), value_type(13));
        witness[0] = value_type(3);
        recovery_type::append_witness(witness);
        BOOST_CHECK(converted.is_satisfied(value_type(3), witness));
    }
}

BOOST_AUTO_TEST_CASE(instance_map_preserves_empty_and_constant_products_in_source_order) {
    source_system_type source;
    source.primary_input_size = 1;
    source.auxiliary_input_size = 2;
    source.add_constraint({{}, variable_type(2), {}});                                // 0*x = 0.
    source.add_constraint({variable_type(0), variable_type(0), variable_type(0)});    // 1*1 = 1.
    source.add_constraint(source_constraint_type());                                  // 0*0 = 0.
    const auto converted = reduction_type::instance_map(source);
    BOOST_REQUIRE(converted.is_valid());
    BOOST_CHECK_EQUAL(converted.num_variables(), 8);
    BOOST_CHECK_EQUAL(converted.num_constraints(), 9);

    const std::vector<recovery_type::constraint_type> expected_rows = {
        {target_combination({{1, 1}}), target_combination({{0, -1}, {5, 1}})},
        {target_combination({{1, -1}}), target_combination({{0, -1}, {5, 1}})},
        {target_combination({{3, 2}}), target_combination({{0, -1}, {3, 4}, {6, 1}})},
        {{}, target_combination({{0, -1}, {6, 1}})},
        {{}, target_combination({{0, -1}, {7, 1}})},
        {{}, target_combination({{0, -1}, {7, 1}})},
    };
    const std::vector<recovery_type::constraint_type> source_rows(converted.constraints.begin() + 3,
                                                                  converted.constraints.end());
    BOOST_CHECK(source_rows == expected_rows);
    std::vector<value_type> witness = {value_type(3), value_type(4), value_type(7)};
    recovery_type::append_witness(witness);
    witness.insert(witness.end(), {value_type(16), value_type::zero(), value_type::zero()});
    BOOST_CHECK(converted.is_satisfied(value_type(3), witness));

    // Structurally valid, unsatisfiable rows are converted too: changing the last row to 0*0 = 1 must reject.
    source.constraints.back().c.add_term(variable_type(0));
    const auto impossible = reduction_type::instance_map(source);
    BOOST_CHECK(impossible.is_valid());
    BOOST_CHECK(!impossible.is_satisfied(value_type(3), witness));
}

BOOST_AUTO_TEST_CASE(maps_reject_invalid_public_input_counts) {
    for (const auto public_size : {std::size_t(0), std::size_t(2), std::numeric_limits<std::size_t>::max()}) {
        auto source = source_example();
        source.primary_input_size = public_size;
        const auto original = source;
        BOOST_CHECK_THROW(reduction_type::instance_map(source), std::invalid_argument);
        BOOST_CHECK_THROW(reduction_type::witness_map(source, {value_type(3)}, {value_type(2), value_type(1)}),
                          std::invalid_argument);
        BOOST_CHECK(source == original);
    }
}

BOOST_AUTO_TEST_CASE(maps_check_every_raw_source_index_before_normalization_or_evaluation) {
    using source_combination_type = nil::crypto3::math::linear_combination<variable_type>;
    const std::array<source_combination_type source_constraint_type::*, 3> combinations = {
        &source_constraint_type::a, &source_constraint_type::b, &source_constraint_type::c};
    // With N = 3, even index 4 is invalid. Index 5 could otherwise alias the newly introduced auxiliary.
    for (const auto index : {std::size_t(4), std::size_t(5), std::numeric_limits<std::size_t>::max()}) {
        for (std::size_t combination = 0; combination < combinations.size(); ++combination) {
            for (unsigned form = 0; form < 3; ++form) {
                BOOST_TEST_CONTEXT("index=" << index << ", combination=" << combination << ", form=" << form) {
                    auto source = source_example();
                    auto &terms = source.constraints.front().*combinations[combination];
                    terms.add_term(variable_type(index), form == 0 ? value_type::zero() : value_type::one());
                    if (form == 2) {
                        terms.add_term(variable_type(index), -value_type::one());
                    }
                    const auto original = source;
                    BOOST_CHECK_THROW(reduction_type::instance_map(source), std::invalid_argument);
                    BOOST_CHECK_THROW(
                        reduction_type::witness_map(source, {value_type(3)}, {value_type(2), value_type(1)}),
                        std::invalid_argument);
                    BOOST_CHECK(source == original);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(maps_reject_witness_size_overflow_without_allocation) {
    const auto max_size = std::numeric_limits<std::size_t>::max();
    const auto max_witness_size = std::vector<value_type>().max_size();
    for (const auto auxiliary_size : {max_size, max_size - 1, max_witness_size - recovery_type::witness_size,
                                      max_witness_size - recovery_type::witness_size - 1}) {
        auto source = source_example();
        source.auxiliary_input_size = auxiliary_size;
        const auto original = source;
        BOOST_CHECK_THROW(reduction_type::instance_map(source), std::invalid_argument);
        BOOST_CHECK_THROW(reduction_type::witness_map(source, {value_type(3)}, {value_type(2), value_type(1)}),
                          std::invalid_argument);
        BOOST_CHECK(source == original);
    }
}

BOOST_AUTO_TEST_CASE(maps_reject_unsupported_domain) {
    // The BN254 base field supports only a size-two radix-two domain, smaller than constant recovery needs.
    using base_field = fields::alt_bn128_base_field<254>;
    snark::r1cs_constraint_system<base_field> source;
    source.primary_input_size = 1;
    BOOST_CHECK_THROW(snark::reductions::r1cs_to_modified_sap<base_field>::instance_map(source), std::invalid_argument);
    BOOST_CHECK_THROW(
        snark::reductions::r1cs_to_modified_sap<base_field>::witness_map(source, {base_field::value_type(3)}, {}),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(witness_map_matches_reference_and_preserves_inputs) {
    const auto source = source_example();
    const auto original = source;
    const auto converted = reduction_type::instance_map(source);
    for (const auto &u : {value_type(1), value_type(3), value_type(5), value_type(field_type::modulus - 2)}) {
        const value_type x = u - value_type::one();
        const primary_input_type primary = {u};
        const auxiliary_input_type auxiliary = {x, value_type::one()};
        const auto original_primary = primary;
        const auto original_auxiliary = auxiliary;
        const auto witness = reduction_type::witness_map(source, primary, auxiliary);

        // Independently compute recovery and the source auxiliary ((x + 1) - 1)^2.
        const std::vector<value_type> expected = {
            u, x, value_type::one(), value_type::one(), (value_type::one() - u).squared(), x.squared()};
        BOOST_CHECK(witness == expected);
        BOOST_CHECK_EQUAL(witness.size(), converted.num_variables());
        BOOST_CHECK(converted.is_satisfied(u, witness));
        BOOST_CHECK(source == original);
        BOOST_CHECK(primary == original_primary);
        BOOST_CHECK(auxiliary == original_auxiliary);
    }
}

BOOST_AUTO_TEST_CASE(witness_map_preserves_product_order_and_unused_variables) {
    auto source = source_example();
    source.auxiliary_input_size = 3;                      // The last original auxiliary is unused.
    source.add_constraint({{}, variable_type(2), {}});    // 0*x = 0.
    source.add_constraint({variable_type(0), variable_type(0), variable_type(0)});    // 1*1 = 1.
    source.add_constraint(source_constraint_type());                                  // 0*0 = 0.
    source_constraint_type signed_product {variable_type(0), variable_type(3), -variable_type(1)};
    signed_product.a.add_term(variable_type(2), -value_type(2));    // (1 - 2*x)*y = -u.
    source.add_constraint(signed_product);

    const primary_input_type primary = {value_type(3)};
    const auxiliary_input_type auxiliary = {value_type(2), value_type(1), value_type(7)};
    BOOST_REQUIRE(source.is_satisfied(primary, auxiliary));
    const auto converted = reduction_type::instance_map(source);
    const auto witness = reduction_type::witness_map(source, primary, auxiliary);
    BOOST_REQUIRE_EQUAL(witness.size(), 11);
    const std::vector<value_type> prefix = {value_type(3), value_type(2), value_type(1), value_type(7)};
    BOOST_CHECK_EQUAL_COLLECTIONS(witness.begin(), witness.begin() + 4, prefix.begin(), prefix.end());
    BOOST_CHECK(witness[4].is_one());    // Recovery starts after all original variables, including unused ones.
    const std::vector<value_type> products = {value_type(4), value_type(4), value_type(0), value_type(0),
                                              value_type(16)};
    BOOST_CHECK_EQUAL_COLLECTIONS(witness.end() - 5, witness.end(), products.begin(), products.end());
    BOOST_CHECK(converted.is_satisfied(primary[0], witness));
}

BOOST_AUTO_TEST_CASE(witness_map_accepts_unsorted_and_repeated_terms) {
    auto source = source_example();
    auto &row = source.constraints.front();
    for (auto *combination : {&row.a, &row.b, &row.c}) {
        std::reverse(combination->terms.begin(), combination->terms.end());
        combination->add_term(variable_type(0), value_type(9));
        combination->add_term(variable_type(0), -value_type(9));
        combination->add_term(variable_type(2), value_type::zero());
    }
    BOOST_CHECK(!source.is_valid());    // Raw terms need not satisfy the ordinary ordering requirement.
    const auto original = source;
    const primary_input_type primary = {value_type(3)};
    const auxiliary_input_type auxiliary = {value_type(2), value_type(1)};
    const auto witness = reduction_type::witness_map(source, primary, auxiliary);
    BOOST_CHECK(witness == reduction_type::witness_map(source_example(), primary, auxiliary));
    BOOST_CHECK(reduction_type::instance_map(source).is_satisfied(primary[0], witness));
    BOOST_CHECK(source == original);
}

BOOST_AUTO_TEST_CASE(witness_map_handles_sources_without_constraints) {
    for (const std::size_t auxiliary_size : {0, 3}) {
        source_system_type source;
        source.primary_input_size = 1;
        source.auxiliary_input_size = auxiliary_size;
        const primary_input_type primary = {value_type(3)};
        const auxiliary_input_type auxiliary(auxiliary_size, value_type(13));
        const auto witness = reduction_type::witness_map(source, primary, auxiliary);
        BOOST_REQUIRE_EQUAL(witness.size(), source.num_variables() + recovery_type::witness_size);
        BOOST_CHECK(witness[0] == primary[0]);
        BOOST_CHECK_EQUAL_COLLECTIONS(witness.begin() + 1, witness.begin() + source.num_variables(), auxiliary.begin(),
                                      auxiliary.end());
        BOOST_CHECK(witness[source.num_variables()].is_one());
        BOOST_CHECK(reduction_type::instance_map(source).is_satisfied(primary[0], witness));
    }
}

BOOST_AUTO_TEST_CASE(witness_map_rejects_input_length_mismatches) {
    const auto source = source_example();
    const auto original = source;
    for (const std::size_t primary_size : {0, 1, 2}) {
        for (const std::size_t auxiliary_size : {0, 1, 2, 3}) {
            if (primary_size == 1 && auxiliary_size == 2) {
                continue;
            }
            BOOST_TEST_CONTEXT("primary size " << primary_size << ", auxiliary size " << auxiliary_size) {
                const primary_input_type primary(primary_size, value_type(3));
                auxiliary_input_type auxiliary(auxiliary_size, value_type(1));
                if (!auxiliary.empty()) {
                    auxiliary[0] = value_type(2);
                }
                const auto original_primary = primary;
                const auto original_auxiliary = auxiliary;
                BOOST_CHECK_THROW(reduction_type::witness_map(source, primary, auxiliary), std::invalid_argument);
                BOOST_CHECK(source == original);
                BOOST_CHECK(primary == original_primary);
                BOOST_CHECK(auxiliary == original_auxiliary);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(witness_map_rejects_even_canonical_public_inputs) {
    const auto source = source_example();
    const auto original = source;
    for (const auto &u : {value_type(0), value_type(2), value_type(field_type::modulus - 1)}) {
        const primary_input_type primary = {u};
        const auxiliary_input_type auxiliary = {u - value_type::one(), value_type::one()};
        const auto original_primary = primary;
        const auto original_auxiliary = auxiliary;
        BOOST_REQUIRE(source.is_satisfied(primary, auxiliary));
        BOOST_CHECK_THROW(reduction_type::witness_map(source, primary, auxiliary), std::invalid_argument);
        BOOST_CHECK(source == original);
        BOOST_CHECK(primary == original_primary);
        BOOST_CHECK(auxiliary == original_auxiliary);
    }
}

BOOST_AUTO_TEST_CASE(witness_map_agrees_with_source_on_small_assignments) {
    const auto source = source_example();
    const auto converted = reduction_type::instance_map(source);
    for (const unsigned u : {1, 3, 5, 7}) {
        for (unsigned x = 0; x < 5; ++x) {
            for (unsigned y = 0; y < 5; ++y) {
                BOOST_TEST_CONTEXT("u=" << u << ", x=" << x << ", y=" << y) {
                    const primary_input_type primary = {value_type(u)};
                    const auxiliary_input_type auxiliary = {value_type(x), value_type(y)};
                    const bool satisfied = (x + 1) * y == u;
                    BOOST_CHECK_EQUAL(source.is_satisfied(primary, auxiliary), satisfied);
                    if (satisfied) {
                        const auto witness = reduction_type::witness_map(source, primary, auxiliary);
                        BOOST_CHECK(converted.is_satisfied(primary[0], witness));
                    } else {
                        BOOST_CHECK_THROW(reduction_type::witness_map(source, primary, auxiliary),
                                          std::invalid_argument);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(witness_map_checks_later_source_rows) {
    auto source = source_example();
    source.add_constraint({variable_type(0), variable_type(0), {}});    // Unsatisfiable: 1*1 = 0.
    const auto original = source;
    const primary_input_type primary = {value_type(3)};
    const auxiliary_input_type auxiliary = {value_type(2), value_type(1)};
    BOOST_REQUIRE(!source.is_satisfied(primary, auxiliary));
    BOOST_CHECK_THROW(reduction_type::witness_map(source, primary, auxiliary), std::invalid_argument);
    BOOST_CHECK(source == original);
}

BOOST_AUTO_TEST_SUITE_END()
