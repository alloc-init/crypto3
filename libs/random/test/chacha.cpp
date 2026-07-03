//---------------------------------------------------------------------------//
// Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
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

#define BOOST_TEST_MODULE chacha_rng_test

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/mac/algorithm/compute.hpp>
#include <nil/crypto3/random/chacha.hpp>
#include <nil/crypto3/stream/algorithm/encrypt.hpp>

using namespace nil::crypto3;

namespace {
    typedef random::chacha<> rng_type;
    typedef rng_type::mac_type mac_type;
    typedef rng_type::cipher_key_type cipher_key_type;
    typedef rng_type::nonce_type nonce_type;

    struct test_rng : rng_type {
        using rng_type::set_output_offset_for_testing;
    };

    struct rng_state {
        cipher_key_type hk;
        cipher_key_type ck;
    };

    template<typename InputRange>
    rng_state reseed_reference(const cipher_key_type &hk, const InputRange &seed_material) {
        rng_state state = {};
        state.hk = hk;

        mac::mac_key<mac_type> hmac_key(state.hk);
        typename mac_type::digest_type digest = static_cast<typename mac_type::digest_type>(
            compute<mac_type>(seed_material.begin(), seed_material.end(), hmac_key));
        std::copy(digest.begin(), digest.end(), state.ck.begin());

        nonce_type nonce = {0};
        std::array<std::uint8_t, std::tuple_size<cipher_key_type>::value> zeros = {0};
        stream::chacha20_cipher hk_cipher(state.ck, nonce, 0);
        hk_cipher.process(zeros.begin(), zeros.end(), state.hk.begin());

        return state;
    }

    std::vector<std::uint8_t> chacha20_output_reference(const cipher_key_type &key, std::size_t size) {
        nonce_type nonce = {0};
        std::vector<std::uint8_t> zeros(size, 0);
        std::vector<std::uint8_t> output(size, 0);
        encrypt(zeros.begin(), zeros.end(), key, nonce, output.begin(), 1);
        return output;
    }

    template<typename InputRange>
    std::vector<std::uint8_t> output_after_seed_reference(const InputRange &seed_material, std::size_t size) {
        cipher_key_type hk = {0};
        const rng_state state = reseed_reference(hk, seed_material);
        return chacha20_output_reference(state.ck, size);
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(chacha_rng_tests)

BOOST_AUTO_TEST_CASE(default_seed_is_deterministic) {
    rng_type rng;
    const std::vector<std::uint8_t> empty_seed;

    const rng_type::result_type output = rng();
    const std::vector<std::uint8_t> expected = output_after_seed_reference(empty_seed, rng_type::result_size);

    BOOST_TEST(output.size() == rng_type::result_size);
    BOOST_TEST(std::equal(expected.begin(), expected.end(), output.begin()));
    BOOST_TEST(rng.entropy() == 0.0);
    BOOST_TEST(rng_type::min() == 0u);
}

BOOST_AUTO_TEST_CASE(seed_uses_hmac_output_as_chacha20_key) {
    const std::vector<std::uint8_t> seed = {0x00, 0x01, 0x02, 0x03, 0x20, 0x21, 0x22, 0x23};
    rng_type rng(seed);

    std::vector<std::uint8_t> output(96, 0);
    rng.generate(output.begin(), output.end());

    const std::vector<std::uint8_t> expected = output_after_seed_reference(seed, output.size());
    BOOST_TEST(std::equal(expected.begin(), expected.end(), output.begin()));
}

BOOST_AUTO_TEST_CASE(generate_continues_across_calls_and_chacha20_blocks) {
    const std::vector<std::uint8_t> seed = {0xa0, 0xa1, 0xa2, 0xa3, 0xa4};

    rng_type single(seed);
    std::vector<std::uint8_t> single_output(150, 0);
    single.generate(single_output.begin(), single_output.end());

    rng_type split(seed);
    std::vector<std::uint8_t> split_output(150, 0);
    split.generate(split_output.begin(), split_output.begin() + 17);
    split.generate(split_output.begin() + 17, split_output.begin() + 80);
    split.generate(split_output.begin() + 80, split_output.end());

    const std::vector<std::uint8_t> expected = output_after_seed_reference(seed, single_output.size());

    BOOST_TEST(std::equal(expected.begin(), expected.end(), single_output.begin()));
    BOOST_TEST(std::equal(single_output.begin(), single_output.end(), split_output.begin()));
}

BOOST_AUTO_TEST_CASE(zero_length_generate_preserves_state) {
    const std::vector<std::uint8_t> seed = {0xb0, 0xb1, 0xb2, 0xb3};

    rng_type uninterrupted(seed);
    std::vector<std::uint8_t> uninterrupted_output(96, 0);
    uninterrupted.generate(uninterrupted_output.begin(), uninterrupted_output.end());

    rng_type with_empty_generate(seed);
    std::vector<std::uint8_t> output(96, 0);
    with_empty_generate.generate(output.begin(), output.begin() + 31);
    with_empty_generate.generate(output.begin() + 31, output.begin() + 31);
    with_empty_generate.generate(output.begin() + 31, output.end());

    BOOST_TEST(std::equal(uninterrupted_output.begin(), uninterrupted_output.end(), output.begin()));
}

BOOST_AUTO_TEST_CASE(reseed_after_partial_output_restarts_from_new_seed) {
    const std::vector<std::uint8_t> first_seed = {0xc0, 0xc1, 0xc2, 0xc3};
    const std::vector<std::uint8_t> second_seed = {0xd0, 0xd1, 0xd2, 0xd3};

    rng_type rng(first_seed);
    std::array<std::uint8_t, 17> discarded = {0};
    rng.generate(discarded.begin(), discarded.end());
    rng.seed(second_seed);

    std::vector<std::uint8_t> output(96, 0);
    rng.generate(output.begin(), output.end());

    cipher_key_type zero_hk = {0};
    const rng_state first_state = reseed_reference(zero_hk, first_seed);
    const rng_state second_state = reseed_reference(first_state.hk, second_seed);
    const std::vector<std::uint8_t> expected = chacha20_output_reference(second_state.ck, output.size());

    BOOST_TEST(std::equal(expected.begin(), expected.end(), output.begin()));
}

BOOST_AUTO_TEST_CASE(repeated_operator_call_continues_stream) {
    const std::vector<std::uint8_t> seed = {0xe0, 0xe1, 0xe2, 0xe3};

    rng_type rng(seed);
    const rng_type::result_type first = rng();
    const rng_type::result_type second = rng();

    const std::vector<std::uint8_t> expected = output_after_seed_reference(seed, first.size() + second.size());

    BOOST_REQUIRE(first.size() == rng_type::result_size);
    BOOST_REQUIRE(second.size() == rng_type::result_size);
    BOOST_TEST(std::equal(expected.begin(), expected.begin() + first.size(), first.begin()));
    BOOST_TEST(std::equal(expected.begin() + first.size(), expected.end(), second.begin()));
}

BOOST_AUTO_TEST_CASE(reseed_uses_refreshed_hmac_key) {
    const std::vector<std::uint8_t> first_seed = {0x10, 0x11, 0x12, 0x13};
    const std::vector<std::uint8_t> second_seed = {0x70, 0x71, 0x72, 0x73};

    rng_type rng(first_seed);
    rng.seed(second_seed);

    std::vector<std::uint8_t> output(64, 0);
    rng.generate(output.begin(), output.end());

    cipher_key_type zero_hk = {0};
    const rng_state first_state = reseed_reference(zero_hk, first_seed);
    const rng_state second_state = reseed_reference(first_state.hk, second_seed);
    const std::vector<std::uint8_t> expected = chacha20_output_reference(second_state.ck, output.size());

    rng_type fresh(second_seed);
    std::vector<std::uint8_t> fresh_output(64, 0);
    fresh.generate(fresh_output.begin(), fresh_output.end());

    BOOST_TEST(std::equal(expected.begin(), expected.end(), output.begin()));
    BOOST_TEST(!std::equal(fresh_output.begin(), fresh_output.end(), output.begin()));
}

BOOST_AUTO_TEST_CASE(output_limit_allows_final_byte_and_rejects_wrap) {
    test_rng rng;
    std::array<std::uint8_t, 1> output = {0};

    rng.set_output_offset_for_testing(rng_type::max_output_bytes_per_seed - 1);
    rng.generate(output.begin(), output.end());

    BOOST_CHECK_THROW(rng.generate(output.begin(), output.end()), std::out_of_range);

    rng.seed();
    rng.set_output_offset_for_testing(rng_type::max_output_bytes_per_seed);
    BOOST_CHECK_THROW(rng.generate(output.begin(), output.end()), std::out_of_range);
}

BOOST_AUTO_TEST_SUITE_END()
