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

#define BOOST_TEST_MODULE modified_sap_transcript_test

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <utility>

#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/transcript.hpp>

#include <nil/crypto3/algebra/algorithms/pair.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/policy.hpp>

namespace {

    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using scalar_field_type = curve_type::scalar_field_type;
    using scalar_value_type = scalar_field_type::value_type;
    using base_value_type = curve_type::base_field_type::value_type;
    using g1_value_type = curve_type::g1_type<>::value_type;
    using g2_value_type = curve_type::g2_type<>::value_type;
    using gt_value_type = curve_type::gt_type::value_type;
    using pairing_policy_type = nil::crypto3::zk::snark::modified_sap_bn254_exact_pairing_policy;
    using transcript_policy_type = nil::crypto3::zk::snark::modified_sap_bn254_poseidon_transcript_policy;
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

    gt_value_type gt_generator() {
        const auto result = nil::crypto3::algebra::pair_reduced<curve_type, pairing_policy_type>(
            g1_value_type::one(), g2_value_type::one());
        if (!result) {
            throw std::runtime_error("modified SAP test pairing failed");
        }
        return *result;
    }

    transcript_policy_type::verification_key_type test_verification_key() {
        const auto gT = gt_generator();
        transcript_policy_type::verification_key_type result;
        result.g2_one = g2_value_type::one();
        result.tau_g2 = scalar_value_type(2) * g2_value_type::one();
        result.gamma_inverse_g2 = scalar_value_type(3) * g2_value_type::one();
        result.alpha_z_vanishing_gt = gT.pow(5);
        result.alpha_gt = {gT.pow(7), gT.pow(11), gT.pow(13), gT.pow(17)};
        result.num_variables = 3;
        result.domain_size = 4;
        result.circuit_digest = transcript_policy_type::circuit_digest(test_system());
        return result;
    }

    g1_value_type alternate_coordinates(const g1_value_type &point) {
        const auto affine = point.to_affine();
        const base_value_type scale(2);
        const auto scale_squared = scale.squared();
        return {affine.X * scale_squared, affine.Y * scale_squared * scale, scale};
    }

    g2_value_type alternate_coordinates(const g2_value_type &point) {
        const auto affine = point.to_affine();
        const typename g2_value_type::field_type::value_type scale(base_value_type(2));
        const auto scale_squared = scale.squared();
        return {affine.X * scale_squared, affine.Y * scale_squared * scale, scale};
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(modified_sap_transcript_test_suite)

BOOST_AUTO_TEST_CASE(circuit_digest_matches_fixed_vector) {
    const auto digest = transcript_policy_type::circuit_digest(test_system());
    const transcript_policy_type::digest_type expected(
        0x304c9ffefccc2d862eacb58cd202a7aa62596aed78f08c394c4e2f35b01935c1_cppui_modular254);
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

BOOST_AUTO_TEST_CASE(proof_challenge_matches_fixed_vector) {
    const auto challenge = transcript_policy_type::proof_challenge(
        test_verification_key(), scalar_value_type(19) * g1_value_type::one(), scalar_value_type(23));
    const scalar_value_type expected(
        0x06b57277e5dacbb0f3e4959c180b6e8b28e8bbe1f9d0765e1256265265595941_cppui_modular254);
    BOOST_CHECK_EQUAL(challenge, expected);
}

BOOST_AUTO_TEST_CASE(proof_challenge_is_coordinate_independent) {
    const auto verification_key = test_verification_key();
    const auto P = scalar_value_type(19) * g1_value_type::one();
    const auto expected = transcript_policy_type::proof_challenge(verification_key, P, scalar_value_type(23));

    const auto alternate_P = alternate_coordinates(P);
    BOOST_REQUIRE(alternate_P == P);
    BOOST_CHECK_EQUAL(transcript_policy_type::proof_challenge(verification_key, alternate_P, scalar_value_type(23)),
                      expected);

    auto alternate_key = verification_key;
    alternate_key.g2_one = alternate_coordinates(verification_key.g2_one);
    alternate_key.tau_g2 = alternate_coordinates(verification_key.tau_g2);
    alternate_key.gamma_inverse_g2 = alternate_coordinates(verification_key.gamma_inverse_g2);
    BOOST_REQUIRE(alternate_key.g2_one == verification_key.g2_one);
    BOOST_REQUIRE(alternate_key.tau_g2 == verification_key.tau_g2);
    BOOST_REQUIRE(alternate_key.gamma_inverse_g2 == verification_key.gamma_inverse_g2);
    BOOST_CHECK_EQUAL(transcript_policy_type::proof_challenge(alternate_key, P, scalar_value_type(23)), expected);
}

BOOST_AUTO_TEST_CASE(proof_challenge_handles_identity_points) {
    transcript_policy_type::verification_key_type identity_key;
    identity_key.g2_one = g2_value_type::zero();
    identity_key.tau_g2 = g2_value_type::zero();
    identity_key.gamma_inverse_g2 = g2_value_type::zero();
    identity_key.num_variables = 1;
    identity_key.domain_size = 2;

    const auto challenge = transcript_policy_type::proof_challenge(
        identity_key, g1_value_type::zero(), scalar_value_type::one());
    BOOST_CHECK_EQUAL(challenge,
                      transcript_policy_type::proof_challenge(
                          identity_key, g1_value_type::zero(), scalar_value_type::one()));
}

BOOST_AUTO_TEST_CASE(proof_challenge_binds_every_input) {
    const auto verification_key = test_verification_key();
    const auto P = scalar_value_type(19) * g1_value_type::one();
    const scalar_value_type u(23);
    const auto expected = transcript_policy_type::proof_challenge(verification_key, P, u);
    const auto gT = gt_generator();

    const auto check_key_change = [&](const auto &changed_key) {
        BOOST_CHECK_NE(transcript_policy_type::proof_challenge(changed_key, P, u), expected);
    };

    auto changed_key = verification_key;
    changed_key.circuit_digest += base_value_type::one();
    check_key_change(changed_key);

    changed_key = verification_key;
    ++changed_key.num_variables;
    check_key_change(changed_key);

    changed_key = verification_key;
    changed_key.domain_size *= 2;
    check_key_change(changed_key);

    changed_key = verification_key;
    changed_key.g2_one += g2_value_type::one();
    check_key_change(changed_key);

    changed_key = verification_key;
    changed_key.tau_g2 += g2_value_type::one();
    check_key_change(changed_key);

    changed_key = verification_key;
    changed_key.gamma_inverse_g2 += g2_value_type::one();
    check_key_change(changed_key);

    changed_key = verification_key;
    changed_key.alpha_z_vanishing_gt *= gT;
    check_key_change(changed_key);

    for (std::size_t i = 0; i < verification_key.alpha_gt.size(); ++i) {
        changed_key = verification_key;
        changed_key.alpha_gt[i] *= gT;
        check_key_change(changed_key);
    }

    BOOST_CHECK_NE(transcript_policy_type::proof_challenge(verification_key, P + g1_value_type::one(), u),
                   expected);
    BOOST_CHECK_NE(transcript_policy_type::proof_challenge(
                       verification_key, P, u + scalar_value_type::one()),
                   expected);
}

BOOST_AUTO_TEST_SUITE_END()
