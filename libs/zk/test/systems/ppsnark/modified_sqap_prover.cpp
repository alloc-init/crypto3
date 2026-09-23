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

#define BOOST_TEST_MODULE modified_sqap_prover_test

#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>

#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/prover.hpp>

#include <nil/crypto3/random/chacha_urbg.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/generator.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/policy.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/transcript.hpp>

namespace {

    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using scalar_value_type = curve_type::scalar_field_type::value_type;
    using base_value_type = curve_type::base_field_type::value_type;
    using g1_value_type = curve_type::g1_type<>::value_type;
    using pairing_policy_type = nil::crypto3::zk::snark::modified_sqap_bn254_exact_pairing_policy;
    using transcript_policy_type = nil::crypto3::zk::snark::modified_sqap_bn254_poseidon_transcript_policy;
    using policy_type =
        nil::crypto3::zk::snark::modified_sqap_policy<curve_type, pairing_policy_type, transcript_policy_type>;
    using prover_type = nil::crypto3::zk::snark::modified_sqap_prover<policy_type>;
    using generator_type = nil::crypto3::zk::snark::modified_sqap_generator<policy_type>;
    using proving_key_type = policy_type::proving_key_type;
    using system_type = policy_type::constraint_system_type;
    using constraint_type = system_type::constraint_type;
    using variable_type = constraint_type::variable_type;
    using combination_type = constraint_type::linear_combination_type;

    combination_type raw_combination(std::initializer_list<std::pair<std::size_t, int>> terms) {
        combination_type result;
        for (const auto &[index, coefficient] : terms) {
            result.add_term(variable_type(index), scalar_value_type(coefficient));
        }
        return result;
    }

    constraint_type binding_row() {
        return {{}, combination_type(-variable_type(0))};
    }

    system_type test_system() {
        return {3,
                {binding_row(),
                 {raw_combination({{1, 2}, {1, -1}}), raw_combination({{2, 1}})},
                 {raw_combination({{2, 1}}), raw_combination({{2, -2}})}}};
    }

    proving_key_type generate_proving_key(const system_type &source) {
        // Fixed seed for reproducible tests with the production setup and transcript.
        const std::array<std::uint8_t, 32> seed = {};
        nil::crypto3::random::chacha_urbg<> random_source(seed);
        return generator_type::process(source, random_source).first;
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(modified_sqap_prover_test_suite)

BOOST_AUTO_TEST_CASE(generated_keys_validate_without_modification) {
    for (const std::size_t logical_rows : {1, 2, 3, 4, 5}) {
        BOOST_TEST_CONTEXT("logical rows: " << logical_rows) {
            auto source = test_system();
            source.constraints.resize(logical_rows, binding_row());
            const auto key = generate_proving_key(source);
            const auto original = key;

            BOOST_CHECK_NO_THROW(prover_type::validate(key));
            BOOST_CHECK(key == original);
        }
    }
}

BOOST_AUTO_TEST_CASE(minimum_domain_accepts_only_an_empty_h_opening_query) {
    auto key = generate_proving_key({1, {binding_row()}});
    BOOST_REQUIRE_EQUAL(key.verification_key.domain_size, 2);
    BOOST_REQUIRE(key.T_H.empty());
    BOOST_CHECK_NO_THROW(prover_type::validate(key));

    key.T_H.push_back(g1_value_type::zero());
    BOOST_CHECK_THROW(prover_type::validate(key), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(rejects_short_and_long_queries) {
    const auto original = generate_proving_key(test_system());
    const std::array queries = {
        std::pair {"W", &proving_key_type::W},     std::pair {"H_query", &proving_key_type::H_query},
        std::pair {"T_A", &proving_key_type::T_A}, std::pair {"T_C", &proving_key_type::T_C},
        std::pair {"T_H", &proving_key_type::T_H}, std::pair {"T_Z", &proving_key_type::T_Z}};

    for (const auto &[name, member] : queries) {
        BOOST_TEST_CONTEXT("query: " << name) {
            BOOST_REQUIRE(!(original.*member).empty());
            auto shorter = original;
            (shorter.*member).pop_back();
            BOOST_CHECK_THROW(prover_type::validate(shorter), std::invalid_argument);

            auto longer = original;
            (longer.*member).push_back(g1_value_type::zero());
            BOOST_CHECK_THROW(prover_type::validate(longer), std::invalid_argument);
        }
    }
}

BOOST_AUTO_TEST_CASE(rejects_inconsistent_verification_key_dimensions) {
    const auto original = generate_proving_key(test_system());
    for (const auto n : {std::size_t(0), std::size_t(2), std::size_t(4), std::numeric_limits<std::size_t>::max()}) {
        BOOST_TEST_CONTEXT("num_variables: " << n) {
            auto key = original;
            key.verification_key.num_variables = n;
            BOOST_CHECK_THROW(prover_type::validate(key), std::invalid_argument);
        }
    }
    for (const auto m : {std::size_t(0), std::size_t(1), std::size_t(2), std::size_t(3), std::size_t(5), std::size_t(8),
                         std::numeric_limits<std::size_t>::max()}) {
        BOOST_TEST_CONTEXT("domain_size: " << m) {
            auto key = original;
            key.verification_key.domain_size = m;
            BOOST_CHECK_THROW(prover_type::validate(key), std::invalid_argument);
        }
    }
}

BOOST_AUTO_TEST_CASE(rejects_invalid_constraint_systems) {
    BOOST_CHECK_THROW(prover_type::validate(proving_key_type {}), std::invalid_argument);
    const auto original = generate_proving_key(test_system());

    auto zero_variables = original;
    zero_variables.constraint_system.witness_size = 0;
    BOOST_CHECK_THROW(prover_type::validate(zero_variables), std::invalid_argument);

    auto empty_rows = original;
    empty_rows.constraint_system.constraints.clear();
    BOOST_CHECK_THROW(prover_type::validate(empty_rows), std::invalid_argument);

    auto missing_binding = original;
    missing_binding.constraint_system.constraints.front() = original.constraint_system.constraints[1];
    BOOST_CHECK_THROW(prover_type::validate(missing_binding), std::invalid_argument);

    auto out_of_range = original;
    out_of_range.constraint_system.constraints[1].a.add_term(
        variable_type(out_of_range.constraint_system.num_variables()), scalar_value_type::zero());
    BOOST_CHECK_THROW(prover_type::validate(out_of_range), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(rejects_changed_circuit_dimensions) {
    const auto original = generate_proving_key(test_system());

    auto more_variables = original;
    ++more_variables.constraint_system.witness_size;
    BOOST_CHECK_THROW(prover_type::validate(more_variables), std::invalid_argument);

    auto more_rows = original;
    more_rows.constraint_system.constraints.resize(5, binding_row());
    BOOST_CHECK_THROW(prover_type::validate(more_rows), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(rejects_changed_digest_or_same_shape_circuit) {
    const auto original = generate_proving_key(test_system());

    auto changed_digest = original;
    changed_digest.verification_key.circuit_digest += base_value_type::one();
    BOOST_CHECK_THROW(prover_type::validate(changed_digest), std::invalid_argument);

    auto changed_coefficient = original;
    changed_coefficient.constraint_system.constraints[1].a.terms.front().coeff += scalar_value_type::one();
    BOOST_REQUIRE(changed_coefficient.constraint_system.is_valid());
    BOOST_CHECK_THROW(prover_type::validate(changed_coefficient), std::invalid_argument);

    auto changed_index = original;
    changed_index.constraint_system.constraints[1].a.terms.front().index = 2;
    BOOST_REQUIRE(changed_index.constraint_system.is_valid());
    BOOST_CHECK_THROW(prover_type::validate(changed_index), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(digest_binds_logical_rows_and_witness_dimension) {
    const auto original = generate_proving_key(test_system());

    // A fourth logical row preserves the padded domain size and every query length.
    auto more_rows = original;
    more_rows.constraint_system.constraints.push_back(binding_row());
    BOOST_REQUIRE(more_rows.constraint_system.is_valid());
    BOOST_CHECK_THROW(prover_type::validate(more_rows), std::invalid_argument);

    // Matching the new witness dimension in the metadata and W leaves the digest stale.
    auto more_variables = original;
    ++more_variables.constraint_system.witness_size;
    ++more_variables.verification_key.num_variables;
    more_variables.W.push_back(g1_value_type::zero());
    BOOST_REQUIRE(more_variables.constraint_system.is_valid());
    BOOST_CHECK_THROW(prover_type::validate(more_variables), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
