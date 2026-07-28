//---------------------------------------------------------------------------//
// Copyright (c) 2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2021 Ilias Khairullin <ilias@nil.foundation>
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

#define BOOST_TEST_MODULE zk_transcript_test

#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <nil/crypto3/algebra/curves/bls12.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/bls12.hpp>
#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/curves/mnt4.hpp>
#include <nil/crypto3/algebra/curves/mnt6.hpp>
#include <nil/crypto3/algebra/curves/pallas.hpp>

#include <nil/crypto3/hash/block_to_field_elements_wrapper.hpp>
#include <nil/crypto3/hash/poseidon.hpp>
#include <nil/crypto3/hash/keccak.hpp>
#include <nil/crypto3/hash/sha2.hpp>
#include <nil/crypto3/hash/h2f.hpp>

#include <nil/crypto3/zk/transcript/fiat_shamir.hpp>

#include <nil/crypto3/marshalling/algebra/processing/curve_element.hpp>

using namespace nil::crypto3;
using namespace nil::crypto3::zk;

BOOST_AUTO_TEST_SUITE(zk_transcript_test_suite)

BOOST_AUTO_TEST_CASE(zk_transcript_manual_test) {
    using field_type = algebra::curves::alt_bn128_254::scalar_field_type;
    std::vector<std::uint8_t> init_blob {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    transcript::fiat_shamir_heuristic_sequential<hashes::keccak_1600<256>> tr(init_blob);
    auto ch1 = tr.challenge<field_type>();
    auto ch2 = tr.challenge<field_type>();
    auto ch_n = tr.challenges<field_type, 3>();

    BOOST_CHECK_EQUAL(ch1, field_type::value_type(0xe858ba005424eabd6d97de7e930779def59a85c1a9ff7e8a5d001cdb07f6e4_cppui_modular254));
    BOOST_CHECK_EQUAL(ch2, field_type::value_type(0xf61f38f58a55b3bbee0480fc5ec3cf8df81603579f4f7134f764bfd3ca5938b_cppui_modular254));

    BOOST_CHECK_EQUAL(ch_n[0], field_type::value_type(0x4f6b97a9bc99d6996fab5e03d1cd0b418a9b3c97ed64cca070e15777e7cc99a_cppui_modular254));
    BOOST_CHECK_EQUAL(ch_n[1], field_type::value_type(0x2414ddf7ecff246500beb2c01b0c5912a400bc3cdca6d7f24bd2bd4987b21e04_cppui_modular254));
    BOOST_CHECK_EQUAL(ch_n[2], field_type::value_type(0x10bfe2f4a414eec551dda5fd9899e9b46e327648b4fa564ed0517b6a99396aec_cppui_modular254));
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(zk_poseidon_transcript_test_suite)

// Keep the standard Poseidon transcript behavior stable across refactoring.
BOOST_AUTO_TEST_CASE(zk_poseidon_transcript_init_test) {
    using field_type = algebra::curves::alt_bn128_254::scalar_field_type;
    using poseidon_type =
        hashes::poseidon<nil::crypto3::hashes::detail::poseidon1_policy<field_type, 128, 2>>;

    std::vector<std::uint8_t> init_blob {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    transcript::fiat_shamir_heuristic_sequential<poseidon_type> tr{
        hashes::block_to_field_elements_wrapper<field_type, std::vector<std::uint8_t>>(init_blob)
    };
    auto ch1 = tr.challenge<field_type>();
    auto ch2 = tr.challenge<field_type>();
    int ch_int = tr.int_challenge<int>();

    BOOST_CHECK_EQUAL(
        ch1,
        field_type::value_type(
            0x13442afffe340d61e9965271e0b4ff5760c90543d2cedb6f55c50c6cfec28ff9_cppui_modular254));
    BOOST_CHECK_EQUAL(
        ch2,
        field_type::value_type(
            0x2089a04399956cc5c98374dfb958b7af3e5b6e0cf1ce3549372bf544951c7259_cppui_modular254));
    BOOST_CHECK_EQUAL(ch_int, 0xfef0);

    init_blob = {};
    tr = transcript::fiat_shamir_heuristic_sequential<poseidon_type>(init_blob);
    ch1 = tr.challenge<field_type>();
    ch2 = tr.challenge<field_type>();
    ch_int = tr.int_challenge<int>();

    BOOST_CHECK_EQUAL(
        ch1,
        field_type::value_type(
            0x0ee069e6aa796ef0e46cbd51d10468393d443a00f5affe72898d9ab62e335e16_cppui_modular254));
    BOOST_CHECK_EQUAL(
        ch2,
        field_type::value_type(
            0x0b6bcbfc77af6d252b1b1f46f0549f88337fc46446c17c75e30b47d5a681e88a_cppui_modular254));
    BOOST_CHECK_EQUAL(ch_int, 0x740);
}

BOOST_AUTO_TEST_CASE(zk_poseidon_transcript_no_init_test) {
    using field_type = algebra::curves::alt_bn128_254::scalar_field_type;
    using poseidon_type =
        hashes::poseidon<nil::crypto3::hashes::detail::poseidon1_policy<field_type, 128, 2>>;

    transcript::fiat_shamir_heuristic_sequential<poseidon_type> tr;
    auto ch1 = tr.challenge<field_type>();
    auto ch2 = tr.challenge<field_type>();
    int ch_int = tr.int_challenge<int>();

    BOOST_CHECK_EQUAL(
        ch1,
        field_type::value_type(
            0x0ee069e6aa796ef0e46cbd51d10468393d443a00f5affe72898d9ab62e335e16_cppui_modular254));
    BOOST_CHECK_EQUAL(
        ch2,
        field_type::value_type(
            0x0b6bcbfc77af6d252b1b1f46f0549f88337fc46446c17c75e30b47d5a681e88a_cppui_modular254));
    BOOST_CHECK_EQUAL(ch_int, 0x740);
}

BOOST_AUTO_TEST_SUITE_END()

/* TODO: Write more elaborate tests for transcript of curve elements */
BOOST_AUTO_TEST_SUITE(transcript_test_curves)

template<typename curve_type, typename hash_type>
void test_transcript(typename curve_type::base_field_type::value_type const& expected_value)
{
    using field_type = typename curve_type::base_field_type;
    using g1_type = typename curve_type::template g1_type<>;

    transcript::fiat_shamir_heuristic_sequential<hash_type> transcript;

    transcript(g1_type::value_type::one());
    auto challenge = transcript.template challenge<field_type>();

    BOOST_CHECK_EQUAL(challenge, expected_value);
}

BOOST_AUTO_TEST_CASE(mnt4_keccak) {
    test_transcript<algebra::curves::mnt4_298, hashes::keccak_1600<256>>
        (0x2b4e9c317f18745b6b89cbb97728923a2d797e261f8320d90f204192c7aabd2b397a0cc155c_cppui_modular298);
}

BOOST_AUTO_TEST_CASE(mnt6_keccak) {
    test_transcript<algebra::curves::mnt6_298, hashes::keccak_1600<256>>
        (0x25a45c6b7d107961d135e640abfb1840cefd9c8ea318f7f33cd327cd55dabdd18c125d6c6b_cppui_modular298);
}

// bc: no bls12_381 keccak h2f suite
// BOOST_AUTO_TEST_CASE(bls12_keccak) {
//     test_transcript<algebra::curves::bls12_381, hashes::keccak_1600<256>>
//         (0x7cc24317960f68f067e0a1cfe610fe3db024d52b064ff2115ea0f594602f0784_cppui_modular381);
//     /* TODO: no marshalling for bls12-377 curve
//     test_transcript<algebra::curves::bls12_377, hashes::keccak_1600<256>>
//         (0x0_cppui_modular377);
//     */
// }

BOOST_AUTO_TEST_SUITE_END()
