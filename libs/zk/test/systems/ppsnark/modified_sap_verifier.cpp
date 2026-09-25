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

#define BOOST_TEST_MODULE modified_sap_verifier_test

#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <stdexcept>

#include <nil/crypto3/algebra/algorithms/pair.hpp>
#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/curves/detail/scalar_mul.hpp>

#include <nil/crypto3/zk/algorithms/verify.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/policy.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/transcript.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/verifier.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap_snark.hpp>

namespace {

    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using native_pairing_policy_type = nil::crypto3::algebra::pairing::pairing_policy<curve_type>;
    using pairing_policy_type = nil::crypto3::zk::snark::modified_sap_bn254_exact_pairing_policy;
    using transcript_policy_type = nil::crypto3::zk::snark::modified_sap_bn254_poseidon_transcript_policy;
    using policy_type =
        nil::crypto3::zk::snark::modified_sap_policy<curve_type, pairing_policy_type, transcript_policy_type>;
    using verifier_type = nil::crypto3::zk::snark::modified_sap_verifier<policy_type>;
    using scheme_type = nil::crypto3::zk::snark::modified_sap_snark<policy_type>;
    using native_policy_type =
        nil::crypto3::zk::snark::modified_sap_policy<curve_type, native_pairing_policy_type, transcript_policy_type>;
    using native_verifier_type = nil::crypto3::zk::snark::modified_sap_verifier<native_policy_type>;
    using scalar_field_type = curve_type::scalar_field_type;
    using scalar_value_type = scalar_field_type::value_type;
    using base_value_type = curve_type::base_field_type::value_type;
    using g1_value_type = curve_type::g1_type<>::value_type;
    using g2_value_type = curve_type::g2_type<>::value_type;
    using g2_field_value_type = g2_value_type::field_type::value_type;
    using gt_value_type = curve_type::gt_type::value_type;

    gt_value_type gt_generator() {
        const auto result = nil::crypto3::algebra::pair_reduced<curve_type, pairing_policy_type>(g1_value_type::one(),
                                                                                                 g2_value_type::one());
        if (!result) {
            throw std::runtime_error("modified SAP test pairing failed");
        }
        return *result;
    }

    policy_type::verification_key_type valid_verification_key() {
        const auto gT = gt_generator();
        policy_type::verification_key_type result;
        result.g2_one = g2_value_type::one();
        result.tau_g2 = scalar_value_type(2) * g2_value_type::one();
        result.gamma_inverse_g2 = scalar_value_type(3) * g2_value_type::one();
        result.alpha_z_vanishing_gt = gT.pow(5);
        result.alpha_gt = {gT.pow(7), gT.pow(11), gT.pow(13), gT.pow(17)};
        result.num_variables = 3;
        result.domain_size = 4;
        result.circuit_digest = base_value_type(19);
        return result;
    }

    policy_type::proof_type valid_shape_proof() {
        policy_type::proof_type result;
        result.P = scalar_value_type(23) * g1_value_type::one();
        result.Q = scalar_value_type(29) * g1_value_type::one();
        return result;
    }

    struct verifier_vector {
        policy_type::verification_key_type verification_key;
        scalar_value_type public_input;
        policy_type::proof_type proof;
    };

    verifier_vector independent_verifier_vector() {
        const scalar_value_type tau(2);
        const scalar_value_type gamma_inverse(3);
        const std::array<scalar_value_type, 4> alpha = {scalar_value_type(5), scalar_value_type(7),
                                                        scalar_value_type(11), scalar_value_type(13)};
        const scalar_value_type alpha_z_vanishing = alpha[3] * (tau.pow(4) - scalar_value_type::one());
        const auto gT = gt_generator();

        verifier_vector result;
        result.verification_key.g2_one = g2_value_type::one();
        result.verification_key.tau_g2 = tau * g2_value_type::one();
        result.verification_key.gamma_inverse_g2 = gamma_inverse * g2_value_type::one();
        result.verification_key.alpha_z_vanishing_gt = gT.pow(alpha_z_vanishing.to_integral());
        for (std::size_t i = 0; i < alpha.size(); ++i) {
            result.verification_key.alpha_gt[i] = gT.pow(alpha[i].to_integral());
        }
        result.verification_key.num_variables = 3;
        result.verification_key.domain_size = 4;
        result.verification_key.circuit_digest = base_value_type(31);

        result.public_input = scalar_value_type(3);
        result.proof.P = scalar_value_type(17) * g1_value_type::one();
        result.proof.v_A = scalar_value_type(5);
        result.proof.v_C = scalar_value_type(7);
        result.proof.v_H = scalar_value_type(3);
        result.proof.v_Z = scalar_value_type(5);
        // Fixed q solving the target-group exponent equation for this transcript vector.
        result.proof.Q =
            scalar_value_type(0x0a91ceaa34d1bc72e8a7707e083fe3f3aea4127294778619ed36a148ef36374e_cppui_modular254) *
            g1_value_type::one();
        return result;
    }

    g2_value_type wrong_subgroup_g2() {
        for (std::size_t i = 0; i < 1024; ++i) {
            const g2_field_value_type x(base_value_type(i), base_value_type::zero());
            const auto rhs = x * x.squared() + g2_value_type::params_type::b;
            if (!rhs.is_square()) {
                continue;
            }
            const g2_value_type point(x, rhs.sqrt(), g2_field_value_type::one());
            if (point.is_well_formed() && !nil::crypto3::algebra::curves::detail::subgroup_check(point)) {
                return point;
            }
        }
        throw std::runtime_error("modified SAP test could not construct a non-subgroup G2 point");
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(modified_sap_verifier_test_suite)

BOOST_AUTO_TEST_CASE(accepts_valid_shapes_and_identity_proof_points) {
    const auto verification_key = valid_verification_key();
    const scalar_value_type u(3);

    BOOST_CHECK(verifier_type::validate(verification_key, u, valid_shape_proof()));
    BOOST_CHECK(verifier_type::validate(verification_key, u, policy_type::proof_type()));

    auto identity_gt_key = verification_key;
    identity_gt_key.alpha_z_vanishing_gt = gt_value_type::one();
    identity_gt_key.alpha_gt.fill(gt_value_type::one());
    BOOST_CHECK(verifier_type::validate(identity_gt_key, u, valid_shape_proof()));
}

BOOST_AUTO_TEST_CASE(rejects_even_public_input_and_invalid_metadata) {
    const auto proof = valid_shape_proof();
    const auto verification_key = valid_verification_key();
    BOOST_CHECK(!verifier_type::validate(verification_key, scalar_value_type(2), proof));

    auto changed_key = verification_key;
    changed_key.num_variables = 0;
    BOOST_CHECK(!verifier_type::validate(changed_key, scalar_value_type(3), proof));

    for (const std::size_t domain_size : {0, 1, 3}) {
        changed_key = verification_key;
        changed_key.domain_size = domain_size;
        BOOST_CHECK(!verifier_type::validate(changed_key, scalar_value_type(3), proof));
    }

    changed_key = verification_key;
    constexpr std::size_t max_domain_log = nil::crypto3::algebra::fields::arithmetic_params<scalar_field_type>::s;
    changed_key.domain_size = std::size_t(1) << max_domain_log;
    BOOST_CHECK(verifier_type::validate(changed_key, scalar_value_type(3), proof));

    changed_key.domain_size = std::size_t(1) << (max_domain_log + 1);
    BOOST_CHECK(!verifier_type::validate(changed_key, scalar_value_type(3), proof));
}

BOOST_AUTO_TEST_CASE(rejects_malformed_and_degenerate_g2_parameters) {
    const auto proof = valid_shape_proof();
    const auto verification_key = valid_verification_key();
    const scalar_value_type u(3);

    auto changed_key = verification_key;
    changed_key.g2_one = scalar_value_type(2) * g2_value_type::one();
    BOOST_CHECK(!verifier_type::validate(changed_key, u, proof));

    changed_key = verification_key;
    changed_key.tau_g2 = g2_value_type::zero();
    BOOST_CHECK(!verifier_type::validate(changed_key, u, proof));

    changed_key = verification_key;
    changed_key.gamma_inverse_g2 = g2_value_type::zero();
    BOOST_CHECK(!verifier_type::validate(changed_key, u, proof));

    changed_key = verification_key;
    changed_key.tau_g2.X += g2_field_value_type::one();
    BOOST_REQUIRE(!changed_key.tau_g2.is_well_formed());
    BOOST_CHECK(!verifier_type::validate(changed_key, u, proof));

    changed_key = verification_key;
    changed_key.tau_g2 = wrong_subgroup_g2();
    BOOST_REQUIRE(changed_key.tau_g2.is_well_formed());
    BOOST_CHECK(!verifier_type::validate(changed_key, u, proof));
}

BOOST_AUTO_TEST_CASE(rejects_invalid_gt_parameters) {
    const auto proof = valid_shape_proof();
    const auto verification_key = valid_verification_key();
    const scalar_value_type u(3);

    auto check_rejected = [&](const auto &changed_key) {
        BOOST_CHECK(!verifier_type::validate(changed_key, u, proof));
    };

    auto changed_key = verification_key;
    changed_key.alpha_z_vanishing_gt = gt_value_type::zero();
    check_rejected(changed_key);

    for (std::size_t i = 0; i < verification_key.alpha_gt.size(); ++i) {
        changed_key = verification_key;
        changed_key.alpha_gt[i] = gt_value_type::zero();
        check_rejected(changed_key);
    }

    auto outside_subgroup = gt_value_type::one();
    outside_subgroup.data[0].data[0].data[0] += base_value_type::one();
    BOOST_REQUIRE(!outside_subgroup.is_zero());
    BOOST_REQUIRE_NE(outside_subgroup.pow(scalar_field_type::modulus), gt_value_type::one());

    changed_key = verification_key;
    changed_key.alpha_z_vanishing_gt = outside_subgroup;
    check_rejected(changed_key);
    for (std::size_t i = 0; i < verification_key.alpha_gt.size(); ++i) {
        changed_key = verification_key;
        changed_key.alpha_gt[i] = outside_subgroup;
        check_rejected(changed_key);
    }
}

BOOST_AUTO_TEST_CASE(rejects_malformed_proof_points) {
    const auto verification_key = valid_verification_key();
    const scalar_value_type u(3);

    auto proof = valid_shape_proof();
    proof.P.X += base_value_type::one();
    BOOST_REQUIRE(!proof.P.is_well_formed());
    BOOST_CHECK(!verifier_type::validate(verification_key, u, proof));

    proof = valid_shape_proof();
    proof.Q.X += base_value_type::one();
    BOOST_REQUIRE(!proof.Q.is_well_formed());
    BOOST_CHECK(!verifier_type::validate(verification_key, u, proof));
}

BOOST_AUTO_TEST_CASE(verifies_an_independently_constructed_equation_vector) {
    const auto vector = independent_verifier_vector();
    const scalar_value_type expected_challenge(
        0x23c9f909c9627ebe6e047c9d238fe7b429aa3ee6131b4ebab795dab7a975ada1_cppui_modular254);
    const scalar_value_type q(0x0a91ceaa34d1bc72e8a7707e083fe3f3aea4127294778619ed36a148ef36374e_cppui_modular254);
    const scalar_value_type expected_left =
        scalar_value_type(17) * scalar_value_type(3) + scalar_value_type(13) * scalar_value_type(15);
    const scalar_value_type expected_evaluations =
        scalar_value_type(5) * scalar_value_type(5) + scalar_value_type(7) * scalar_value_type(7) +
        scalar_value_type(11) * scalar_value_type(3) + scalar_value_type(13) * scalar_value_type(5);

    BOOST_REQUIRE_EQUAL(vector.proof.v_A.squared() - vector.proof.v_C,
                        vector.proof.v_H * vector.proof.v_Z + vector.public_input);
    BOOST_REQUIRE_EQUAL(
        transcript_policy_type::proof_challenge(vector.verification_key, vector.proof.P, vector.public_input),
        expected_challenge);
    BOOST_REQUIRE_EQUAL(expected_left, q * (scalar_value_type(2) - expected_challenge) + expected_evaluations);
    BOOST_REQUIRE_EQUAL(vector.proof.Q, q * g1_value_type::one());
    BOOST_CHECK(verifier_type::process(vector.verification_key, vector.public_input, vector.proof));
    BOOST_CHECK(scheme_type::verify(vector.verification_key, vector.public_input, vector.proof));
    BOOST_CHECK(nil::crypto3::zk::verify<scheme_type>(vector.verification_key, vector.public_input, vector.proof));
    BOOST_CHECK(!native_verifier_type::process(vector.verification_key, vector.public_input, vector.proof));
}

BOOST_AUTO_TEST_CASE(verifies_identity_commitments_without_affine_inversion) {
    auto vector = independent_verifier_vector();
    vector.verification_key.alpha_z_vanishing_gt = gt_value_type::one();
    vector.verification_key.alpha_gt.fill(gt_value_type::one());
    vector.proof = policy_type::proof_type();
    vector.proof.v_A = scalar_value_type(2);
    vector.proof.v_C = scalar_value_type(1);

    BOOST_CHECK(verifier_type::process(vector.verification_key, vector.public_input, vector.proof));
}

BOOST_AUTO_TEST_CASE(rejects_failed_arithmetic_and_pairing_equations) {
    const auto vector = independent_verifier_vector();

    auto changed_proof = vector.proof;
    changed_proof.v_A += scalar_value_type::one();
    BOOST_CHECK(!verifier_type::process(vector.verification_key, vector.public_input, changed_proof));

    changed_proof = vector.proof;
    changed_proof.Q += g1_value_type::one();
    BOOST_CHECK(!verifier_type::process(vector.verification_key, vector.public_input, changed_proof));

    changed_proof = vector.proof;
    changed_proof.v_H = scalar_value_type(4);
    changed_proof.v_C = scalar_value_type(2);
    BOOST_REQUIRE_EQUAL(changed_proof.v_A.squared() - changed_proof.v_C,
                        changed_proof.v_H * changed_proof.v_Z + vector.public_input);
    BOOST_CHECK(!verifier_type::process(vector.verification_key, vector.public_input, changed_proof));
}

BOOST_AUTO_TEST_CASE(rejects_changed_public_input_commitment_and_parameters) {
    const auto vector = independent_verifier_vector();

    auto changed_proof = vector.proof;
    changed_proof.v_C = scalar_value_type(5);
    const scalar_value_type changed_input(5);
    BOOST_REQUIRE_EQUAL(changed_proof.v_A.squared() - changed_proof.v_C,
                        changed_proof.v_H * changed_proof.v_Z + changed_input);
    BOOST_CHECK(!verifier_type::process(vector.verification_key, changed_input, changed_proof));

    changed_proof = vector.proof;
    changed_proof.P += g1_value_type::one();
    BOOST_CHECK(!verifier_type::process(vector.verification_key, vector.public_input, changed_proof));

    auto changed_key = vector.verification_key;
    changed_key.circuit_digest += base_value_type::one();
    BOOST_CHECK(!verifier_type::process(changed_key, vector.public_input, vector.proof));

    changed_key = vector.verification_key;
    changed_key.tau_g2 += g2_value_type::one();
    BOOST_CHECK(!verifier_type::process(changed_key, vector.public_input, vector.proof));

    changed_key = vector.verification_key;
    changed_key.alpha_gt[0] *= gt_generator();
    BOOST_CHECK(!verifier_type::process(changed_key, vector.public_input, vector.proof));
}

BOOST_AUTO_TEST_SUITE_END()
