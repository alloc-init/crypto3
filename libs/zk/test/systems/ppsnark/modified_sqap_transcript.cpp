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

#define BOOST_TEST_MODULE modified_sqap_transcript_test

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <utility>

#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/transcript.hpp>

namespace {

    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using scalar_field_type = curve_type::scalar_field_type;
    using scalar_value_type = scalar_field_type::value_type;
    using transcript_policy_type = nil::crypto3::zk::snark::modified_sqap_bn254_poseidon_transcript_policy;
    using system_type = transcript_policy_type::constraint_system_type;
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
        return {raw_combination({{2, 0}}), raw_combination({{0, -2}, {1, 0}, {0, 1}})};
    }

    system_type test_system() {
        return {3,
                {binding_row(),
                 {raw_combination({{2, 2}, {1, 3}, {2, 0}, {1, -2}}),
                  raw_combination({{2, 4}, {0, 4}, {0, -1}})},
                 {raw_combination({{1, 5}}), raw_combination({{2, 6}})}}};
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(modified_sqap_transcript_test_suite)

BOOST_AUTO_TEST_CASE(circuit_digest_matches_fixed_vector) {
    const auto digest = transcript_policy_type::circuit_digest(test_system());
    const transcript_policy_type::digest_type expected(
        0x29c79eb7645dd1169cd6c49be2583b99df91b4415a9b86dfeeb04749ec696e4f_cppui_modular254);
    BOOST_CHECK_EQUAL(digest, expected);
}

BOOST_AUTO_TEST_CASE(circuit_digest_normalizes_without_mutating_the_source) {
    const auto source = test_system();
    const auto original = source;

    BOOST_CHECK_EQUAL(transcript_policy_type::circuit_digest(source),
                      transcript_policy_type::circuit_digest(source.normalized()));
    BOOST_CHECK(source == original);
}

BOOST_AUTO_TEST_CASE(circuit_digest_binds_structure_and_domain) {
    const auto source = test_system().normalized();
    const auto digest = transcript_policy_type::circuit_digest(source);

    auto coefficient_change = source;
    coefficient_change.constraints[1].a.terms[0].coeff += scalar_value_type::one();
    BOOST_CHECK_NE(transcript_policy_type::circuit_digest(coefficient_change), digest);

    auto index_change = source;
    index_change.constraints[2].a.terms[0].index = 2;
    BOOST_CHECK_NE(transcript_policy_type::circuit_digest(index_change), digest);

    auto row_order_change = source;
    std::swap(row_order_change.constraints[1], row_order_change.constraints[2]);
    BOOST_CHECK_NE(transcript_policy_type::circuit_digest(row_order_change), digest);

    auto witness_dimension_change = source;
    witness_dimension_change.witness_size = 4;
    BOOST_CHECK_NE(transcript_policy_type::circuit_digest(witness_dimension_change), digest);

    auto domain_change = source;
    domain_change.constraints.push_back(binding_row());
    domain_change.constraints.push_back(binding_row());
    BOOST_CHECK_NE(transcript_policy_type::circuit_digest(domain_change), digest);
}

BOOST_AUTO_TEST_CASE(circuit_digest_rejects_an_invalid_system) {
    BOOST_CHECK_THROW(transcript_policy_type::circuit_digest(system_type()), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
