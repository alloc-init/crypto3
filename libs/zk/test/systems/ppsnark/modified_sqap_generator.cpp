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

#define BOOST_TEST_MODULE modified_sqap_generator_test

#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/algorithms/pair.hpp>
#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/pairing/alt_bn128.hpp>
#include <nil/crypto3/math/linear_combination.hpp>

#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/generator.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/policy.hpp>

namespace {

    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using scalar_field_type = curve_type::scalar_field_type;
    using scalar_value_type = scalar_field_type::value_type;
    using base_value_type = curve_type::base_field_type::value_type;
    using native_pairing_policy_type = nil::crypto3::algebra::pairing::pairing_policy<curve_type>;
    using exact_pairing_policy_type =
        nil::crypto3::zk::snark::modified_sqap_bn254_exact_pairing_policy;

    struct transcript_policy;

    using exact_policy_type = nil::crypto3::zk::snark::modified_sqap_policy<
        curve_type, exact_pairing_policy_type, transcript_policy>;
    using native_policy_type = nil::crypto3::zk::snark::modified_sqap_policy<
        curve_type, native_pairing_policy_type, transcript_policy>;
    using generator_type =
        nil::crypto3::zk::snark::detail::modified_sqap_deterministic_generator<exact_policy_type>;
    using native_generator_type =
        nil::crypto3::zk::snark::detail::modified_sqap_deterministic_generator<native_policy_type>;
    using system_type = exact_policy_type::constraint_system_type;
    using constraint_type = system_type::constraint_type;
    using variable_type = constraint_type::variable_type;
    using combination_type = constraint_type::linear_combination_type;
    using g1_value_type = curve_type::g1_type<>::value_type;
    using g2_value_type = curve_type::g2_type<>::value_type;
    using gt_value_type = curve_type::gt_type::value_type;

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

    system_type two_row_system() {
        return {3,
                {binding_row(),
                 {raw_combination({{1, 3}, {1, -2}, {2, 2}}),
                  raw_combination({{0, 4}, {0, -1}, {2, 4}})}}};
    }

    generator_type::trapdoor_type test_trapdoor() {
        return {scalar_value_type(2),
                scalar_value_type(3),
                {scalar_value_type(5), scalar_value_type(7), scalar_value_type(11), scalar_value_type(13)}};
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(modified_sqap_generator_test_suite)

BOOST_AUTO_TEST_CASE(deterministic_setup_matches_small_independent_example) {
    const auto source = two_row_system();
    const auto original_rows = source.constraints;
    const auto trapdoor = test_trapdoor();
    const base_value_type digest(17);
    const auto keys = generator_type::process(source, digest, trapdoor);
    const auto &pk = keys.first;
    const auto &vk = keys.second;

    BOOST_REQUIRE_EQUAL(pk.W.size(), 3);
    BOOST_REQUIRE_EQUAL(pk.H_query.size(), 1);
    BOOST_REQUIRE_EQUAL(pk.T_A.size(), 1);
    BOOST_REQUIRE_EQUAL(pk.T_C.size(), 1);
    BOOST_CHECK(pk.T_H.empty());
    BOOST_REQUIRE_EQUAL(pk.T_Z.size(), 2);

    const scalar_value_type two_inverse = scalar_value_type(2).inversed();
    const scalar_value_type lagrange_0 = (trapdoor.tau + scalar_value_type::one()) * two_inverse;
    const scalar_value_type lagrange_1 = (scalar_value_type::one() - trapdoor.tau) * two_inverse;
    const std::array<scalar_value_type, 3> at = {
        scalar_value_type::zero(), lagrange_1, scalar_value_type(2) * lagrange_1};
    const std::array<scalar_value_type, 3> ct = {
        -lagrange_0 + scalar_value_type(3) * lagrange_1,
        scalar_value_type::zero(),
        scalar_value_type(4) * lagrange_1};
    for (std::size_t i = 0; i < pk.W.size(); ++i) {
        const auto expected_scalar = trapdoor.gamma * (trapdoor.alpha[0] * at[i] + trapdoor.alpha[1] * ct[i]);
        BOOST_CHECK_EQUAL(pk.W[i], expected_scalar * g1_value_type::one());
    }

    BOOST_CHECK_EQUAL(pk.H_query[0], (trapdoor.gamma * trapdoor.alpha[2]) * g1_value_type::one());
    BOOST_CHECK_EQUAL(pk.T_A[0], trapdoor.alpha[0] * g1_value_type::one());
    BOOST_CHECK_EQUAL(pk.T_C[0], trapdoor.alpha[1] * g1_value_type::one());
    BOOST_CHECK_EQUAL(pk.T_Z[0], trapdoor.alpha[3] * g1_value_type::one());
    BOOST_CHECK_EQUAL(pk.T_Z[1], (trapdoor.alpha[3] * trapdoor.tau) * g1_value_type::one());

    const auto gt_generator = nil::crypto3::algebra::pair_reduced<curve_type, exact_pairing_policy_type>(
        g1_value_type::one(), g2_value_type::one());
    BOOST_REQUIRE(gt_generator);
    const scalar_value_type z_at_tau = trapdoor.tau.squared() - scalar_value_type::one();
    BOOST_CHECK_EQUAL(vk.g2_one, g2_value_type::one());
    BOOST_CHECK_EQUAL(vk.tau_g2, trapdoor.tau * g2_value_type::one());
    BOOST_CHECK_EQUAL(vk.gamma_inverse_g2, trapdoor.gamma.inversed() * g2_value_type::one());
    BOOST_CHECK_EQUAL(vk.alpha_z_vanishing_gt,
                      gt_generator->pow((trapdoor.alpha[3] * z_at_tau).to_integral()));
    for (std::size_t i = 0; i < vk.alpha_gt.size(); ++i) {
        BOOST_CHECK_EQUAL(vk.alpha_gt[i], gt_generator->pow(trapdoor.alpha[i].to_integral()));
    }

    BOOST_CHECK_EQUAL(vk.num_variables, 3);
    BOOST_CHECK_EQUAL(vk.domain_size, 2);
    BOOST_CHECK_EQUAL(vk.circuit_digest, digest);
    BOOST_CHECK(pk.verification_key == vk);
    BOOST_CHECK(pk.constraint_system == source.normalized());
    BOOST_CHECK(source.constraints == original_rows);
}

BOOST_AUTO_TEST_CASE(query_lengths_follow_the_padded_domain) {
    const system_type source {1, {binding_row(), binding_row(), binding_row()}};
    const auto keys = generator_type::process(source, base_value_type(19), test_trapdoor());
    const auto &pk = keys.first;

    BOOST_CHECK_EQUAL(keys.second.domain_size, 4);
    BOOST_CHECK_EQUAL(pk.W.size(), 1);
    BOOST_CHECK_EQUAL(pk.H_query.size(), 3);
    BOOST_CHECK_EQUAL(pk.T_A.size(), 3);
    BOOST_CHECK_EQUAL(pk.T_C.size(), 3);
    BOOST_CHECK_EQUAL(pk.T_H.size(), 2);
    BOOST_CHECK_EQUAL(pk.T_Z.size(), 4);
    BOOST_CHECK_EQUAL(pk.constraint_system.num_constraints(), 3);
}

BOOST_AUTO_TEST_CASE(pairing_convention_is_selected_by_the_scheme_policy) {
    const auto source = two_row_system();
    const auto trapdoor = test_trapdoor();
    const base_value_type digest(23);
    const auto exact_keys = generator_type::process(source, digest, trapdoor);
    const auto native_keys = native_generator_type::process(source, digest, trapdoor);

    BOOST_CHECK(exact_keys.first.W == native_keys.first.W);
    BOOST_CHECK(exact_keys.first.H_query == native_keys.first.H_query);
    BOOST_CHECK(exact_keys.first.T_A == native_keys.first.T_A);
    BOOST_CHECK(exact_keys.first.T_C == native_keys.first.T_C);
    BOOST_CHECK(exact_keys.first.T_H == native_keys.first.T_H);
    BOOST_CHECK(exact_keys.first.T_Z == native_keys.first.T_Z);
    BOOST_CHECK_NE(exact_keys.second.alpha_gt[0], native_keys.second.alpha_gt[0]);
    BOOST_CHECK_NE(exact_keys.second.alpha_z_vanishing_gt, native_keys.second.alpha_z_vanishing_gt);
}

BOOST_AUTO_TEST_CASE(rejects_invalid_circuits_and_zero_gamma) {
    auto zero_gamma = test_trapdoor();
    zero_gamma.gamma = scalar_value_type::zero();
    BOOST_CHECK_THROW(generator_type::process(two_row_system(), base_value_type::zero(), zero_gamma),
                      std::invalid_argument);

    const system_type invalid_system;
    BOOST_CHECK_THROW(generator_type::process(invalid_system, base_value_type::zero(), test_trapdoor()),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
