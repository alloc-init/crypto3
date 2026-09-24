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

#define BOOST_TEST_MODULE modified_sap_pairing_policy_test

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/algorithms/pair.hpp>
#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/pairing/alt_bn128.hpp>
#include <nil/crypto3/algebra/pairing/detail/alt_bn128/params.hpp>

#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/policy.hpp>

namespace {

    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using native_pairing_policy_type = nil::crypto3::algebra::pairing::pairing_policy<curve_type>;
    using exact_pairing_policy_type = nil::crypto3::zk::snark::modified_sap_bn254_exact_pairing_policy;
    using pairing_params_type = nil::crypto3::algebra::pairing::detail::pairing_params<curve_type>;
    using g1_value_type = typename curve_type::g1_type<>::value_type;
    using g2_value_type = typename curve_type::g2_type<>::value_type;
    using gt_value_type = typename curve_type::gt_type::value_type;

}    // namespace

BOOST_AUTO_TEST_SUITE(modified_sap_pairing_policy_test_suite)

BOOST_AUTO_TEST_CASE(exact_result_matches_direct_final_exponentiation) {
    const auto miller_result = nil::crypto3::algebra::pair<curve_type>(g1_value_type::one(), g2_value_type::one());
    const auto expected = miller_result.pow(pairing_params_type::final_exponent);
    const auto exact =
        nil::crypto3::algebra::final_exponentiation<curve_type, exact_pairing_policy_type>(miller_result);
    const auto native =
        nil::crypto3::algebra::final_exponentiation<curve_type, native_pairing_policy_type>(miller_result);

    BOOST_REQUIRE(exact);
    BOOST_REQUIRE(native);
    BOOST_CHECK_EQUAL(*exact, expected);
    BOOST_CHECK_NE(*exact, *native);
}

BOOST_AUTO_TEST_CASE(zero_identity_and_pairing_identity_follow_existing_contracts) {
    const auto zero_result =
        nil::crypto3::algebra::final_exponentiation<curve_type, exact_pairing_policy_type>(gt_value_type::zero());
    const auto one_result =
        nil::crypto3::algebra::final_exponentiation<curve_type, exact_pairing_policy_type>(gt_value_type::one());
    const auto identity_pairing = nil::crypto3::algebra::pair_reduced<curve_type, exact_pairing_policy_type>(
        g1_value_type::zero(), g2_value_type::one());

    BOOST_CHECK(!zero_result);
    BOOST_REQUIRE(one_result);
    BOOST_REQUIRE(identity_pairing);
    BOOST_CHECK_EQUAL(*one_result, gt_value_type::one());
    BOOST_CHECK_EQUAL(*identity_pairing, gt_value_type::one());
}

BOOST_AUTO_TEST_SUITE_END()
