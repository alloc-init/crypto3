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
    snark::modified_sap_constraint_system<FieldType> recovery_system(std::size_t first_bit) {
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<FieldType>;
        using variable = typename recovery::variable_type;
        snark::modified_sap_constraint_system<FieldType> cs;
        cs.witness_size = first_bit + recovery::witness_size;
        cs.constraints.push_back({{}, -variable(0)});
        recovery::append_constraints(cs.constraints, first_bit);
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

    template<typename FieldType>
    void complete_product_auxiliaries(const snark::modified_sap_constraint_system<FieldType> &cs,
                                      std::vector<typename FieldType::value_type> &witness,
                                      std::size_t first_bit) {
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<FieldType>;
        std::size_t auxiliary = first_bit + recovery::bit_count + recovery::marker_count;
        // Every second product row forces t = A(w)^2 when w[0] = u. Satisfy it even for forged bits/markers.
        for (std::size_t row = recovery::bit_count + 3; row < cs.constraints.size(); row += 2) {
            witness[auxiliary++] = cs.constraints[row].a.evaluate(witness).squared();
        }
    }

    template<typename FieldType>
    void complete_comparison_witness(const snark::modified_sap_constraint_system<FieldType> &cs,
                                     std::vector<typename FieldType::value_type> &witness,
                                     std::size_t first_bit) {
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<FieldType>;
        auto marker = witness[first_bit + recovery::bit_count - 1];
        std::size_t next_marker = first_bit + recovery::bit_count;
        // Reference markers are prefix products over the modulus's internal one bits.
        for (std::size_t bit = recovery::bit_count - 1; bit > 1;) {
            --bit;
            if (boost::multiprecision::bit_test(FieldType::modulus, bit)) {
                marker *= witness[first_bit + bit];
                witness[next_marker++] = marker;
            }
        }
        complete_product_auxiliaries(cs, witness, first_bit);
    }

    template<typename FieldType>
    std::vector<typename FieldType::value_type>
        witness_from_integer(const snark::modified_sap_constraint_system<FieldType> &cs,
                             const typename FieldType::integral_type &integer, std::size_t first_bit) {
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<FieldType>;
        using value = typename FieldType::value_type;
        std::vector<value> witness(cs.witness_size, value::zero());
        witness[0] = value(integer);
        for (std::size_t bit = 0; bit < recovery::bit_count; ++bit) {
            witness[first_bit + bit] = value(boost::multiprecision::bit_test(integer, bit));
        }
        complete_comparison_witness(cs, witness, first_bit);
        return witness;
    }

    template<unsigned Modulus>
    void test_all_bit_patterns() {
        using field = test_prime_field<Modulus>;
        using value = typename field::value_type;
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<field>;
        const auto cs = recovery_system<field>(1);
        BOOST_REQUIRE(cs.is_valid());
        for (unsigned integer = 0; integer < (1u << recovery::bit_count); ++integer) {
            BOOST_TEST_CONTEXT("modulus " << Modulus << ", integer " << integer) {
                const auto witness = witness_from_integer(cs, typename field::integral_type(integer), 1);
                BOOST_CHECK_EQUAL(rows_satisfied(cs, witness[0], witness, 0, cs.constraints.size()), integer < Modulus);
                BOOST_CHECK_EQUAL(cs.is_satisfied(witness[0], witness), integer < Modulus && (integer & 1) != 0);
                if (integer < Modulus && (integer & 1) != 0) {
                    std::vector<value> generated = {value(integer)};
                    recovery::append_witness(generated);
                    BOOST_CHECK(generated == witness);
                }
            }
        }
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(r1cs_to_modified_sap_test_suite)

BOOST_AUTO_TEST_CASE(recovery_dimensions_and_source_prefix) {
    BOOST_CHECK_EQUAL(recovery_type::bit_count, 254);
    BOOST_CHECK_EQUAL(recovery_type::marker_count, 99);
    BOOST_CHECK_EQUAL(recovery_type::product_count, 153);
    BOOST_CHECK_EQUAL(recovery_type::witness_size, 506);
    BOOST_CHECK_EQUAL(recovery_type::constraint_count, 561);
    for (const std::size_t first_bit : {1, 3}) {
        const auto cs = recovery_system<field_type>(first_bit);
        BOOST_CHECK_EQUAL(cs.num_variables(), first_bit + 506);
        BOOST_CHECK_EQUAL(cs.num_constraints(), 562);
        BOOST_REQUIRE(cs.is_valid());
        const auto normalized = cs.normalized();
        for (const auto &u : {value_type(1), value_type(3), value_type(field_type::modulus - 2)}) {
            std::vector<value_type> witness(first_bit, value_type(13));
            witness[0] = u;
            const auto prefix = witness;
            recovery_type::append_witness(witness);
            auto expected = witness_from_integer(cs, u.to_integral(), first_bit);
            std::copy(prefix.begin(), prefix.end(), expected.begin());
            BOOST_CHECK(witness == expected);
            BOOST_CHECK_EQUAL_COLLECTIONS(witness.begin(), witness.begin() + first_bit, prefix.begin(), prefix.end());
            BOOST_CHECK_EQUAL(witness.size(), cs.num_variables());
            BOOST_CHECK(witness[first_bit].is_one());
            BOOST_CHECK(cs.is_satisfied(u, witness));
            BOOST_CHECK(normalized.is_satisfied(u, witness));
        }
    }
}

BOOST_AUTO_TEST_CASE(rejects_noncanonical_representatives_even_when_reconstruction_passes) {
    const auto cs = recovery_system<field_type>(1);
    for (const unsigned offset : {0, 1, 3}) {
        BOOST_TEST_CONTEXT("modulus plus " << offset) {
            const field_type::integral_type integer = field_type::modulus + offset;
            const auto witness = witness_from_integer(cs, integer, 1);
            BOOST_CHECK(witness[0] == value_type(offset));
            // Binding, Booleanity and field reconstruction all hold; comparison must reject the lift.
            BOOST_CHECK(rows_satisfied(cs, witness[0], witness, 0, recovery_type::bit_count + 2));
            BOOST_CHECK(!rows_satisfied(cs, witness[0], witness, recovery_type::bit_count + 2, cs.constraints.size()));
            BOOST_CHECK(!cs.is_satisfied(witness[0], witness));
            if (offset != 0) {
                BOOST_CHECK(witness[1].is_zero());
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(rejects_nonboolean_bits_even_when_reconstruction_and_comparison_pass) {
    const auto cs = recovery_system<field_type>(1);
    std::vector<value_type> witness(cs.num_variables(), value_type::zero());
    witness[0] = value_type(3);
    witness[1] = value_type(3);    // Represents 3 through a non-Boolean low bit.
    complete_comparison_witness(cs, witness, 1);
    BOOST_CHECK(rows_satisfied(cs, witness[0], witness, recovery_type::bit_count + 1, cs.constraints.size()));
    BOOST_CHECK(!rows_satisfied(cs, witness[0], witness, 1, recovery_type::bit_count + 1));
    BOOST_CHECK(!cs.is_satisfied(witness[0], witness));
}

BOOST_AUTO_TEST_CASE(rejects_changed_bits_markers_and_product_auxiliaries) {
    const auto cs = recovery_system<field_type>(1);
    const value_type u(field_type::modulus - 2);
    std::vector<value_type> witness = {u};
    recovery_type::append_witness(witness);
    BOOST_REQUIRE(cs.is_satisfied(u, witness));

    for (const std::size_t bit : {0, 1, 127, 253}) {
        auto changed = witness;
        changed[1 + bit] = value_type::one() - changed[1 + bit];
        complete_comparison_witness(cs, changed, 1);
        BOOST_CHECK(!cs.is_satisfied(u, changed));
    }
    for (std::size_t marker = 0; marker < recovery_type::marker_count; ++marker) {
        BOOST_TEST_CONTEXT("marker " << marker) {
            auto changed = witness;
            changed[1 + recovery_type::bit_count + marker] += value_type::one();
            complete_product_auxiliaries(cs, changed, 1);
            BOOST_CHECK(!cs.is_satisfied(u, changed));
        }
    }
    for (std::size_t auxiliary = 0; auxiliary < recovery_type::product_count; ++auxiliary) {
        BOOST_TEST_CONTEXT("product auxiliary " << auxiliary) {
            auto changed = witness;
            changed[1 + recovery_type::bit_count + recovery_type::marker_count + auxiliary] += value_type::one();
            BOOST_CHECK(!cs.is_satisfied(u, changed));
        }
    }
    auto changed = witness;
    changed[0] += value_type::one();
    BOOST_CHECK(!rows_satisfied(cs, u, changed, 0, 1));
    BOOST_CHECK(!cs.is_satisfied(u, changed));
    BOOST_CHECK(!cs.is_satisfied(value_type(1), witness));
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
    for (const auto first_bit : {std::size_t(0), max_size, max_size - recovery_type::witness_size + 1,
                                 std::vector<value_type>().max_size() - recovery_type::witness_size + 1}) {
        BOOST_CHECK_THROW(recovery_type::append_constraints(cs.constraints, first_bit), std::invalid_argument);
        BOOST_CHECK(cs.constraints == original);
    }
}

BOOST_AUTO_TEST_CASE(exhaustive_small_field_bit_patterns) {
    test_all_bit_patterns<3>();
    test_all_bit_patterns<5>();
    test_all_bit_patterns<7>();
    test_all_bit_patterns<11>();
    test_all_bit_patterns<13>();
    test_all_bit_patterns<17>();
    test_all_bit_patterns<31>();
}

BOOST_AUTO_TEST_CASE(exhaustive_bits_and_marker_over_field_seven) {
    using field = test_prime_field<7>;
    using value = field::value_type;
    const auto cs = recovery_system<field>(1);
    // Exhaust field values for all three bits and the single marker; product auxiliaries are forced by the minus rows.
    for (unsigned b0 = 0; b0 < 7; ++b0) {
        for (unsigned b1 = 0; b1 < 7; ++b1) {
            for (unsigned b2 = 0; b2 < 7; ++b2) {
                for (unsigned marker = 0; marker < 7; ++marker) {
                    BOOST_TEST_CONTEXT("bits " << b0 << ',' << b1 << ',' << b2 << ", marker " << marker) {
                        std::vector<value> witness(cs.num_variables(), value::zero());
                        witness[0] = value(b0 + 2 * b1 + 4 * b2);
                        witness[1] = value(b0);
                        witness[2] = value(b1);
                        witness[3] = value(b2);
                        witness[4] = value(marker);
                        complete_product_auxiliaries(cs, witness, 1);
                        const bool valid =
                            b0 <= 1 && b1 <= 1 && b2 <= 1 && b0 + 2 * b1 + 4 * b2 < 7 && marker == b2 * b1;
                        BOOST_CHECK_EQUAL(rows_satisfied(cs, witness[0], witness, 0, cs.constraints.size()), valid);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(instance_map_matches_fixed_rows_and_supports_multiple_public_inputs) {
    const auto source = source_example();
    const auto original = source;
    BOOST_REQUIRE(source.is_valid());
    const auto converted = reduction_type::instance_map(source);
    BOOST_CHECK(source == original);
    BOOST_REQUIRE(converted.is_valid());
    BOOST_CHECK_EQUAL(converted.num_variables(), 510);
    BOOST_CHECK_EQUAL(converted.num_constraints(), 564);

    // N = 3. The recovered one is w[3], and the source product's auxiliary is w[509].
    auto expected = recovery_system<field_type>(3).normalized();
    ++expected.witness_size;
    expected.constraints.push_back(
        {target_combination({{1, 1}, {2, 1}, {3, 1}}), target_combination({{0, 3}, {509, 1}})});
    expected.constraints.push_back(
        {target_combination({{1, 1}, {2, -1}, {3, 1}}), target_combination({{0, -1}, {509, 1}})});
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
        BOOST_CHECK_EQUAL(converted.num_variables(), source.num_variables() + 506);
        BOOST_CHECK_EQUAL(converted.num_constraints(), 562);
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
    BOOST_CHECK_EQUAL(converted.num_variables(), 512);
    BOOST_CHECK_EQUAL(converted.num_constraints(), 568);

    const std::vector<recovery_type::constraint_type> expected_rows = {
        {target_combination({{1, 1}}), target_combination({{0, -1}, {509, 1}})},
        {target_combination({{1, -1}}), target_combination({{0, -1}, {509, 1}})},
        {target_combination({{3, 2}}), target_combination({{0, -1}, {3, 4}, {510, 1}})},
        {{}, target_combination({{0, -1}, {510, 1}})},
        {{}, target_combination({{0, -1}, {511, 1}})},
        {{}, target_combination({{0, -1}, {511, 1}})},
    };
    const std::vector<recovery_type::constraint_type> source_rows(converted.constraints.begin() + 562,
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

BOOST_AUTO_TEST_CASE(instance_map_rejects_invalid_public_input_counts) {
    for (const auto public_size : {std::size_t(0), std::size_t(2), std::numeric_limits<std::size_t>::max()}) {
        auto source = source_example();
        source.primary_input_size = public_size;
        const auto original = source;
        BOOST_CHECK_THROW(reduction_type::instance_map(source), std::invalid_argument);
        BOOST_CHECK(source == original);
    }
}

BOOST_AUTO_TEST_CASE(instance_map_checks_every_raw_source_index_before_normalization) {
    using source_combination_type = nil::crypto3::math::linear_combination<variable_type>;
    const std::array<source_combination_type source_constraint_type::*, 3> combinations = {
        &source_constraint_type::a, &source_constraint_type::b, &source_constraint_type::c};
    // With N = 3, even index 4 is invalid. Index 509 could otherwise alias the newly introduced auxiliary.
    for (const auto index : {std::size_t(4), std::size_t(509), std::numeric_limits<std::size_t>::max()}) {
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
                    BOOST_CHECK(source == original);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(instance_map_rejects_witness_size_overflow_without_allocation) {
    const auto max_size = std::numeric_limits<std::size_t>::max();
    const auto max_witness_size = std::vector<value_type>().max_size();
    for (const auto auxiliary_size : {max_size, max_size - 1, max_witness_size - recovery_type::witness_size,
                                      max_witness_size - recovery_type::witness_size - 1}) {
        auto source = source_example();
        source.auxiliary_input_size = auxiliary_size;
        const auto original = source;
        BOOST_CHECK_THROW(reduction_type::instance_map(source), std::invalid_argument);
        BOOST_CHECK(source == original);
    }
}

BOOST_AUTO_TEST_CASE(instance_map_rejects_unsupported_domain) {
    // The BN254 base field supports only a size-two radix-two domain, smaller than constant recovery needs.
    using base_field = fields::alt_bn128_base_field<254>;
    snark::r1cs_constraint_system<base_field> source;
    source.primary_input_size = 1;
    BOOST_CHECK_THROW(snark::reductions::r1cs_to_modified_sap<base_field>::instance_map(source), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
