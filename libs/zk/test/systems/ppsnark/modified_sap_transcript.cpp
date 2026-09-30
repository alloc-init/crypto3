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

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/transcript.hpp>

#include <nil/crypto3/algebra/algorithms/pair.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_policy.hpp>
#include <nil/crypto3/hash/poseidon.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/policy.hpp>

namespace {

    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using scalar_field_type = curve_type::scalar_field_type;
    using scalar_value_type = scalar_field_type::value_type;
    using base_value_type = curve_type::base_field_type::value_type;
    using base_integral_type = curve_type::base_field_type::integral_type;
    using g1_value_type = curve_type::g1_type<>::value_type;
    using g2_value_type = curve_type::g2_type<>::value_type;
    using gt_value_type = curve_type::gt_type::value_type;
    using pairing_policy_type = nil::crypto3::zk::snark::modified_sap_bn254_exact_pairing_policy;
    using transcript_policy_type = nil::crypto3::zk::snark::modified_sap_bn254_poseidon_transcript_policy;
    using truncated_128_policy = nil::crypto3::zk::snark::modified_sap_bn254_poseidon_truncated_transcript_policy<128>;
    using truncated_192_policy = nil::crypto3::zk::snark::modified_sap_bn254_poseidon_truncated_transcript_policy<192>;
    using transcript_policies = boost::mpl::list<transcript_policy_type, truncated_128_policy, truncated_192_policy>;
    using boundary_policies =
        boost::mpl::list<nil::crypto3::zk::snark::modified_sap_bn254_poseidon_truncated_transcript_policy<1>,
                         nil::crypto3::zk::snark::modified_sap_bn254_poseidon_truncated_transcript_policy<253>>;

    static_assert(truncated_128_policy::challenge_bits == 128);
    static_assert(truncated_192_policy::challenge_bits == 192);
    using poseidon_policy_type = nil::crypto3::hashes::detail::poseidon1_policy<curve_type::base_field_type, 128, 2>;
    using poseidon_type = nil::crypto3::hashes::poseidon1<poseidon_policy_type>;
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
                 {raw_combination({{2, 2}, {1, 3}, {2, 0}, {1, -2}}), raw_combination({{2, 4}, {0, 4}, {0, -1}})},
                 {raw_combination({{1, 5}}), raw_combination({{2, 6}})}}};
    }

    gt_value_type gt_generator() {
        const auto result = nil::crypto3::algebra::pair_reduced<curve_type, pairing_policy_type>(g1_value_type::one(),
                                                                                                 g2_value_type::one());
        if (!result) {
            throw std::runtime_error("modified SAP test pairing failed");
        }
        return *result;
    }

    template<typename TranscriptPolicy = transcript_policy_type>
    typename TranscriptPolicy::verification_key_type test_verification_key() {
        const auto gT = gt_generator();
        typename TranscriptPolicy::verification_key_type result;
        result.g2_one = g2_value_type::one();
        result.tau_g2 = scalar_value_type(2) * g2_value_type::one();
        result.gamma_inverse_g2 = scalar_value_type(3) * g2_value_type::one();
        result.alpha_z_vanishing_gt = gT.pow(5);
        result.alpha_gt = {gT.pow(7), gT.pow(11), gT.pow(13), gT.pow(17)};
        result.num_variables = 3;
        result.domain_size = 4;
        result.circuit_digest = TranscriptPolicy::circuit_digest(test_system());
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

    // Apply the two suffix permutations directly, independently of the full-message sponge path.
    // Copy the prepared state so it can be reused for subsequent commitments and public inputs.
    base_value_type continue_proof_prefix(poseidon_policy_type::state_type state,
                                          const g1_value_type &P,
                                          const scalar_value_type &u) {
        base_value_type y = base_value_type::zero();
        if (P.is_zero()) {
            state[0] = base_value_type::zero();
            state[1] = base_value_type::zero();
        } else {
            const auto affine = P.to_affine();
            state[0] = base_value_type::one();
            state[1] = affine.X;
            y = affine.Y;
        }
        poseidon_type::permutation_type::permute(state);

        state[0] = y;
        state[1] = base_value_type(u.to_integral());
        state[2] += base_value_type::one();    // Existing final-full-block padding in the capacity cell.
        poseidon_type::permutation_type::permute(state);
        return state[0];
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(modified_sap_transcript_test_suite)

BOOST_AUTO_TEST_CASE(circuit_digest_matches_fixed_vector) {
    const auto digest = transcript_policy_type::circuit_digest(test_system());
    const transcript_policy_type::digest_type expected(
        0x304c9ffefccc2d862eacb58cd202a7aa62596aed78f08c394c4e2f35b01935c1_cppui_modular254);
    BOOST_CHECK_EQUAL(digest, expected);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(circuit_digest_normalizes_without_mutating_the_source,
                              TranscriptPolicy,
                              transcript_policies) {
    using transcript_policy_type = TranscriptPolicy;
    const auto source = test_system();
    const auto original = source;

    BOOST_CHECK_EQUAL(transcript_policy_type::circuit_digest(source),
                      transcript_policy_type::circuit_digest(source.normalized()));
    BOOST_CHECK(source == original);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(circuit_digest_binds_structure_and_domain, TranscriptPolicy, transcript_policies) {
    using transcript_policy_type = TranscriptPolicy;
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

BOOST_AUTO_TEST_CASE_TEMPLATE(circuit_digest_rejects_an_invalid_system, TranscriptPolicy, transcript_policies) {
    using transcript_policy_type = TranscriptPolicy;
    BOOST_CHECK_THROW(transcript_policy_type::circuit_digest(system_type()), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(proof_challenge_matches_fixed_vector) {
    const auto challenge = transcript_policy_type::proof_challenge(
        test_verification_key(), scalar_value_type(19) * g1_value_type::one(), scalar_value_type(23));
    const scalar_value_type expected(
        0x06b57277e5dacbb0f3e4959c180b6e8b28e8bbe1f9d0765e1256265265595941_cppui_modular254);
    BOOST_CHECK_EQUAL(challenge, expected);
}

BOOST_AUTO_TEST_CASE(prepared_prefix_matches_fixed_vector) {
    const auto prefix = transcript_policy_type::prepare_proof_prefix(test_verification_key());
    // Derived independently from the documented 86-element encoding using the dense Poseidon1 permutation,
    // without padding. Pin all three cells: continuation overwrites the rate cells before reading them.
    const transcript_policy_type::state_type expected = {
        base_value_type(0x2116e8d872a53ebf26737b7c35555afb17724b900c10359f92277eb799198c30_cppui_modular254),
        base_value_type(0x1cd453dc71eecee560f98da448f773187e0929a3e8ec258412c751a8d217e193_cppui_modular254),
        base_value_type(0x072ec2ca101c92aae4ad233485d5d9c4e475312b241de2c8f33ec3bf46b5d0e6_cppui_modular254)};
    for (std::size_t i = 0; i < prefix.size(); ++i) {
        BOOST_CHECK_EQUAL(prefix[i], expected[i]);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(prepared_prefix_matches_complete_message_digest, TranscriptPolicy, transcript_policies) {
    using transcript_policy_type = TranscriptPolicy;
    const auto verification_key = test_verification_key<transcript_policy_type>();
    const auto prefix = transcript_policy_type::prepare_proof_prefix(verification_key);
    const auto original_prefix = prefix;

    for (const auto &P : {g1_value_type::zero(), g1_value_type::one(), scalar_value_type(19) * g1_value_type::one()}) {
        for (const auto &u :
             {scalar_value_type::zero(), scalar_value_type::one(), scalar_value_type(23), -scalar_value_type::one()}) {
            // proof_digest uses ordinary full-message hashing, independently of this continuation.
            const auto digest = transcript_policy_type::proof_digest(verification_key, P, u);
            BOOST_CHECK_EQUAL(continue_proof_prefix(prefix, P, u), digest);
            const auto challenge = transcript_policy_type::proof_challenge(verification_key, P, u);
            if constexpr (std::is_same_v<TranscriptPolicy,
                                         nil::crypto3::zk::snark::modified_sap_bn254_poseidon_transcript_policy>) {
                BOOST_CHECK_EQUAL(challenge, scalar_value_type(digest.to_integral()));
            } else {
                const base_integral_type bound = base_integral_type(1) << transcript_policy_type::challenge_bits;
                BOOST_CHECK_EQUAL(challenge.to_integral(), digest.to_integral() % bound);
            }
            BOOST_CHECK(prefix == original_prefix);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(proof_challenge_is_coordinate_independent, TranscriptPolicy, transcript_policies) {
    using transcript_policy_type = TranscriptPolicy;
    const auto verification_key = test_verification_key<transcript_policy_type>();
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
    BOOST_CHECK(transcript_policy_type::prepare_proof_prefix(alternate_key) ==
                transcript_policy_type::prepare_proof_prefix(verification_key));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(proof_challenge_handles_identity_points, TranscriptPolicy, transcript_policies) {
    using transcript_policy_type = TranscriptPolicy;
    typename transcript_policy_type::verification_key_type identity_key;
    identity_key.g2_one = g2_value_type::zero();
    identity_key.tau_g2 = g2_value_type::zero();
    identity_key.gamma_inverse_g2 = g2_value_type::zero();
    identity_key.num_variables = 1;
    identity_key.domain_size = 2;

    const auto prefix = transcript_policy_type::prepare_proof_prefix(identity_key);
    const auto digest =
        transcript_policy_type::proof_digest(identity_key, g1_value_type::zero(), scalar_value_type::one());
    BOOST_CHECK_EQUAL(continue_proof_prefix(prefix, g1_value_type::zero(), scalar_value_type::one()), digest);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(proof_challenge_binds_every_input, TranscriptPolicy, transcript_policies) {
    using transcript_policy_type = TranscriptPolicy;
    const auto verification_key = test_verification_key<transcript_policy_type>();
    const auto P = scalar_value_type(19) * g1_value_type::one();
    const scalar_value_type u(23);
    const auto expected = transcript_policy_type::proof_challenge(verification_key, P, u);
    const auto prefix = transcript_policy_type::prepare_proof_prefix(verification_key);
    const auto gT = gt_generator();

    const auto check_key_change = [&](const auto &changed_key) {
        BOOST_CHECK_NE(transcript_policy_type::proof_challenge(changed_key, P, u), expected);
        BOOST_CHECK(transcript_policy_type::prepare_proof_prefix(changed_key) != prefix);
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

    BOOST_CHECK_NE(transcript_policy_type::proof_challenge(verification_key, P + g1_value_type::one(), u), expected);
    BOOST_CHECK_NE(transcript_policy_type::proof_challenge(verification_key, P, u + scalar_value_type::one()),
                   expected);
}

BOOST_AUTO_TEST_CASE(truncation_profile_and_width_bind_both_digests) {
    const auto system = test_system();
    const auto original = transcript_policy_type::circuit_digest(system);
    const auto bits128 = truncated_128_policy::circuit_digest(system);
    const auto bits192 = truncated_192_policy::circuit_digest(system);
    BOOST_CHECK_NE(original, bits128);
    BOOST_CHECK_NE(original, bits192);
    BOOST_CHECK_NE(bits128, bits192);

    // Hold every VK field fixed, including circuit_digest: the proof transcript must bind the
    // conversion and width directly, as well as through the circuit digest stored by setup.
    const auto key = test_verification_key();
    const auto P = scalar_value_type(19) * g1_value_type::one();
    const scalar_value_type u(23);
    const auto original_digest = transcript_policy_type::proof_digest(key, P, u);
    const auto digest128 = truncated_128_policy::proof_digest(key, P, u);
    const auto digest192 = truncated_192_policy::proof_digest(key, P, u);
    BOOST_CHECK_NE(original_digest, digest128);
    BOOST_CHECK_NE(original_digest, digest192);
    BOOST_CHECK_NE(digest128, digest192);
}

BOOST_AUTO_TEST_CASE(truncated_transcripts_match_fixed_vectors) {
    // Independently derived from the normalized circuit and VK encodings with the dense Poseidon1
    // schedule. The reference also reproduces the original profile's unchanged fixed vectors.
    const auto check = []<typename TranscriptPolicy>(const base_value_type &expected_circuit,
                                                     const poseidon_policy_type::state_type &expected_prefix,
                                                     const base_value_type &expected_digest,
                                                     const scalar_value_type &expected_challenge) {
        BOOST_TEST_CONTEXT("challenge bits " << TranscriptPolicy::challenge_bits) {
            const auto key = test_verification_key<TranscriptPolicy>();
            BOOST_CHECK_EQUAL(key.circuit_digest, expected_circuit);
            const auto prefix = TranscriptPolicy::prepare_proof_prefix(key);
            for (std::size_t i = 0; i < prefix.size(); ++i) {
                BOOST_CHECK_EQUAL(prefix[i], expected_prefix[i]);
            }
            const auto P = scalar_value_type(19) * g1_value_type::one();
            const scalar_value_type u(23);
            BOOST_CHECK_EQUAL(TranscriptPolicy::proof_digest(key, P, u), expected_digest);
            BOOST_CHECK_EQUAL(TranscriptPolicy::proof_challenge(key, P, u), expected_challenge);
        }
    };
    check.template operator()<truncated_128_policy>(
        base_value_type(0x037d8fb7bad20135a91cd03e865df8b78a7a6dffc0ddd90e7c10544264551aed_cppui_modular254),
        {base_value_type(0x23f659e5b2e16cc6207244689cd307da933a50ed8fb8061cb411c1e6a792de2_cppui_modular254),
         base_value_type(0x13bb63bdb8ba6057e4947f449cb6281efd8431dfea8478370e245645c2e1fb1f_cppui_modular254),
         base_value_type(0x178a4657899a6fd19b37234c7b30cfcfc0efb816ba32b7a1d61bb9dfdda622a1_cppui_modular254)},
        base_value_type(0x2f65370367edf0565b37a6930518c290443c2397f4b1631b9610ee4f1f015483_cppui_modular254),
        scalar_value_type(0x443c2397f4b1631b9610ee4f1f015483_cppui_modular254));
    check.template operator()<truncated_192_policy>(
        base_value_type(0x2a2e9a1308d3d4321f1ecfca712f47e69670238e67a8ea4c2e4ab628626ebf8f_cppui_modular254),
        {base_value_type(0x1394c0d81987a073ad5752b65f1a6f76eff4f2bb97a8ff1df996a714201a6724_cppui_modular254),
         base_value_type(0x1e90ca47b1ab99dc5c569eb1edde1eecf0c43bda430761393842cfb9d0486baf_cppui_modular254),
         base_value_type(0x15b3b079f37f44b3b77d678dc7e73302e21b4a3e767493b7eeb525e75cc53892_cppui_modular254)},
        base_value_type(0x2a8b320b844036ae4a193e34ac353c7df7a20dd64dcc274582a94984692345e9_cppui_modular254),
        scalar_value_type(0x4a193e34ac353c7df7a20dd64dcc274582a94984692345e9_cppui_modular254));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(truncated_challenge_accepts_boundary_widths, TranscriptPolicy, boundary_policies) {
    const auto key = test_verification_key<TranscriptPolicy>();
    const auto prefix = TranscriptPolicy::prepare_proof_prefix(key);
    const base_integral_type bound = base_integral_type(1) << TranscriptPolicy::challenge_bits;
    bool removed_high_bits = false;
    for (const auto &P : {g1_value_type::zero(), g1_value_type::one(), scalar_value_type(19) * g1_value_type::one()}) {
        for (const auto &u :
             {scalar_value_type::zero(), scalar_value_type::one(), scalar_value_type(23), -scalar_value_type::one()}) {
            const auto digest = TranscriptPolicy::proof_digest(key, P, u);
            const auto challenge = TranscriptPolicy::proof_challenge(key, P, u);
            BOOST_CHECK_EQUAL(continue_proof_prefix(prefix, P, u), digest);
            BOOST_CHECK_LT(challenge.to_integral(), bound);
            BOOST_CHECK_EQUAL(challenge.to_integral(), digest.to_integral() % bound);
            removed_high_bits |= digest.to_integral() >= bound;
        }
    }
    // Ensure the fixtures actually exercise truncation, including at the 253-bit upper boundary.
    BOOST_CHECK(removed_high_bits);
}

BOOST_AUTO_TEST_SUITE_END()
