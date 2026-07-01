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

#define BOOST_TEST_MODULE chacha_cipher_test

#include <algorithm>
#include <array>
#include <cstdint>

#include <boost/test/unit_test.hpp>
#include <boost/endian/conversion.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include <nil/crypto3/stream/algorithm/decrypt.hpp>
#include <nil/crypto3/stream/algorithm/encrypt.hpp>
#include <nil/crypto3/stream/algorithm/seek.hpp>
#include <nil/crypto3/stream/chacha.hpp>
#include <nil/crypto3/stream/detail/chacha/chacha_impl.hpp>

using namespace nil::crypto3::stream;

namespace {
    using block4_type = std::array<std::uint8_t, 256>;
    using schedule_type = std::array<std::uint32_t, 16>;

    bool equal_block(const block4_type &actual, std::size_t block_index,
                     const std::array<std::uint8_t, 64> &expected) {
        return std::equal(expected.begin(), expected.end(), actual.begin() + block_index * 64);
    }

    std::uint32_t rotl(std::uint32_t x, std::size_t n) {
        return static_cast<std::uint32_t>((x << n) | (x >> (32 - n)));
    }

    void quarter_round(std::uint32_t &a, std::uint32_t &b, std::uint32_t &c, std::uint32_t &d) {
        a += b;
        d ^= a;
        d = rotl(d, 16);
        c += d;
        b ^= c;
        b = rotl(b, 12);
        a += b;
        d ^= a;
        d = rotl(d, 8);
        c += d;
        b ^= c;
        b = rotl(b, 7);
    }

    std::array<std::uint8_t, 64> reference_block(const schedule_type &input) {
        schedule_type x = input;

        for (std::size_t r = 0; r != 10; ++r) {
            quarter_round(x[0], x[4], x[8], x[12]);
            quarter_round(x[1], x[5], x[9], x[13]);
            quarter_round(x[2], x[6], x[10], x[14]);
            quarter_round(x[3], x[7], x[11], x[15]);

            quarter_round(x[0], x[5], x[10], x[15]);
            quarter_round(x[1], x[6], x[11], x[12]);
            quarter_round(x[2], x[7], x[8], x[13]);
            quarter_round(x[3], x[4], x[9], x[14]);
        }

        std::array<std::uint8_t, 64> out = {0};
        for (std::size_t i = 0; i != x.size(); ++i) {
            boost::endian::store_little_u32(out.data() + 4 * i, x[i] + input[i]);
        }
        return out;
    }

    schedule_type rfc8439_block_state() {
        return {{0x61707865, 0x3320646e, 0x79622d32, 0x6b206574,
                 0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c,
                 0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c,
                 0x00000001, 0x09000000, 0x4a000000, 0x00000000}};
    }

    std::array<std::uint8_t, 32> rfc8439_key() {
        return {{0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
                 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
                 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f}};
    }

    std::array<std::uint8_t, 12> rfc8439_iv() {
        return {{0x00, 0x00, 0x00, 0x09, 0x00, 0x00, 0x00, 0x4a,
                 0x00, 0x00, 0x00, 0x00}};
    }

    schedule_type original_chacha_state() {
        return {{0x61707865, 0x3320646e, 0x79622d32, 0x6b206574,
                 0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c,
                 0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c,
                 0x00000000, 0x00000000, 0x4a000000, 0x00000000}};
    }

    const std::array<std::uint8_t, 64> rfc8439_expected_block = {
        {0x10, 0xf1, 0xe7, 0xe4, 0xd1, 0x3b, 0x59, 0x15, 0x50, 0x0f, 0xdd, 0x1f, 0xa3, 0x20, 0x71, 0xc4,
         0xc7, 0xd1, 0xf4, 0xc7, 0x33, 0xc0, 0x68, 0x03, 0x04, 0x22, 0xaa, 0x9a, 0xc3, 0xd4, 0x6c, 0x4e,
         0xd2, 0x82, 0x64, 0x46, 0x07, 0x9f, 0xaa, 0x09, 0x14, 0xc2, 0xd7, 0x05, 0xd9, 0x8b, 0x02, 0xa2,
         0xb5, 0x12, 0x9c, 0xd1, 0xde, 0x16, 0x4e, 0xb9, 0xcb, 0xd0, 0x83, 0xe8, 0xa2, 0x50, 0x3c, 0x4e}};
}    // namespace

BOOST_AUTO_TEST_SUITE(chacha_test_suite)

BOOST_AUTO_TEST_CASE(chacha_policy_accepts_supported_round_counts) {
    BOOST_TEST((detail::chacha_policy<8, 96, 256>::rounds == 8u));
    BOOST_TEST((detail::chacha_policy<12, 96, 256>::rounds == 12u));
    BOOST_TEST((detail::chacha_policy<20, 96, 256>::rounds == 20u));
}

BOOST_AUTO_TEST_CASE(ietf_chacha20_x4_matches_rfc8439_first_block) {
    using impl_type = detail::chacha_impl<20, 96, 256>;

    impl_type::key_schedule_type state = rfc8439_block_state();
    block4_type out = {0};

    impl_type::chacha_x4_ietf(out, state);

    BOOST_TEST(equal_block(out, 0, rfc8439_expected_block));
    BOOST_TEST(state[12] == 5u);
    BOOST_TEST(state[13] == 0x09000000u);
    BOOST_TEST(state[14] == 0x4a000000u);
    BOOST_TEST(state[15] == 0x00000000u);
}

BOOST_AUTO_TEST_CASE(ietf_chacha20_x4_does_not_carry_into_nonce_word) {
    using impl_type = detail::chacha_impl<20, 96, 256>;

    impl_type::key_schedule_type state = rfc8439_block_state();
    state[12] = 0xffffffff;
    state[13] = 0x11223344;
    block4_type out = {0};

    impl_type::chacha_x4_ietf(out, state);

    BOOST_TEST(state[12] == 3u);
    BOOST_TEST(state[13] == 0x11223344u);
}

BOOST_AUTO_TEST_CASE(original_chacha_x4_carries_counter_into_second_counter_word) {
    using impl_type = detail::chacha_impl<20, 64, 256>;

    impl_type::key_schedule_type state = original_chacha_state();
    state[12] = 0xffffffff;
    state[13] = 7;
    block4_type out = {0};

    impl_type::chacha_x4_original(out, state);

    BOOST_TEST(state[12] == 3u);
    BOOST_TEST(state[13] == 8u);
}

BOOST_AUTO_TEST_CASE(original_chacha_x4_generates_consecutive_64_bit_counter_blocks) {
    using impl_type = detail::chacha_impl<20, 64, 256>;

    impl_type::key_schedule_type state = original_chacha_state();
    state[12] = 0xfffffffe;
    state[13] = 3;
    block4_type out = {0};

    impl_type::chacha_x4_original(out, state);

    schedule_type expected_state = original_chacha_state();
    expected_state[12] = 0xfffffffe;
    expected_state[13] = 3;

    for (std::size_t i = 0; i != 4; ++i) {
        BOOST_TEST(equal_block(out, i, reference_block(expected_state)));
        ++expected_state[12];
        if (expected_state[12] == 0) {
            ++expected_state[13];
        }
    }
}

BOOST_AUTO_TEST_CASE(chacha_functions_generates_single_block) {
    using functions_type = detail::chacha_functions<20, 96, 256>;

    functions_type::key_schedule_type state = rfc8439_block_state();
    functions_type::block_type out = {0};

    functions_type::generate_block(out, state);

    BOOST_TEST(std::equal(rfc8439_expected_block.begin(), rfc8439_expected_block.end(), out.begin()));
    BOOST_TEST(state[12] == 2u);
    BOOST_TEST(state[13] == 0x09000000u);
    BOOST_TEST(state[14] == 0x4a000000u);
    BOOST_TEST(state[15] == 0x00000000u);
}

BOOST_AUTO_TEST_CASE(public_chacha_constructor_uses_schedule_iv_block_argument) {
    using cipher_type = chacha<96, 256, 20>;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();

    cipher_type cipher(block, schedule, key, iv);
    schedule_type expected_state = rfc8439_block_state();
    expected_state[12] = 0;
    const std::array<std::uint8_t, 64> expected_block = reference_block(expected_state);

    (void)cipher;
    BOOST_TEST(std::equal(expected_block.begin(), expected_block.end(), block.begin()));
    BOOST_TEST(schedule[12] == 1u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);
}

BOOST_AUTO_TEST_CASE(public_chacha_facade_matches_rfc8439_block_vector) {
    using cipher_type = chacha<96, 256, 20>;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();
    cipher_type cipher(block, schedule, key, iv);

    cipher.seek(block, schedule, cipher_type::block_size);

    cipher_type::block_type plaintext = {0};
    cipher_type::block_type ciphertext = {0};
    std::uint8_t *in = plaintext.data();
    std::uint8_t *out = ciphertext.data();

    cipher.process(in, out, schedule, block);

    BOOST_TEST(std::equal(rfc8439_expected_block.begin(), rfc8439_expected_block.end(), ciphertext.begin()));
    BOOST_TEST(schedule[12] == 3u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);
}

BOOST_AUTO_TEST_CASE(public_chacha_seek_ietf_sets_32_bit_counter_and_preserves_nonce) {
    using cipher_type = chacha<96, 256, 20>;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();
    cipher_type cipher(block, schedule, key, iv);

    schedule = rfc8439_block_state();
    cipher.seek(block, schedule, 7 * cipher_type::block_size);

    schedule_type expected_state = rfc8439_block_state();
    expected_state[12] = 7;
    const std::array<std::uint8_t, 64> expected_block = reference_block(expected_state);

    BOOST_TEST(std::equal(expected_block.begin(), expected_block.end(), block.begin()));
    BOOST_TEST(schedule[12] == 8u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);
}

BOOST_AUTO_TEST_CASE(public_chacha_seek_ietf_accepts_maximum_counter_block) {
    using cipher_type = chacha<96, 256, 20>;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();
    cipher_type cipher(block, schedule, key, iv);

    const std::uint64_t counter = 0xffffffffULL;
    schedule = rfc8439_block_state();
    cipher.seek(block, schedule, counter * cipher_type::block_size);

    schedule_type expected_state = rfc8439_block_state();
    expected_state[12] = 0xffffffff;
    const std::array<std::uint8_t, 64> expected_block = reference_block(expected_state);

    BOOST_TEST(std::equal(expected_block.begin(), expected_block.end(), block.begin()));
    BOOST_TEST(schedule[12] == 0u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);
}

BOOST_AUTO_TEST_CASE(public_chacha_seek_original_sets_64_bit_counter) {
    using cipher_type = chacha<64, 256, 20>;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = {0};
    cipher_type::iv_type iv = {0};
    cipher_type cipher(block, schedule, key, iv);

    const std::uint64_t counter = (std::uint64_t(1) << 32) - 1;
    schedule = original_chacha_state();
    cipher.seek(block, schedule, counter * cipher_type::block_size);

    schedule_type expected_state = original_chacha_state();
    expected_state[12] = 0xffffffff;
    const std::array<std::uint8_t, 64> expected_block = reference_block(expected_state);

    BOOST_TEST(std::equal(expected_block.begin(), expected_block.end(), block.begin()));
    BOOST_TEST(schedule[12] == 0u);
    BOOST_TEST(schedule[13] == 1u);
    BOOST_TEST(schedule[14] == original_chacha_state()[14]);
    BOOST_TEST(schedule[15] == original_chacha_state()[15]);
}

BOOST_AUTO_TEST_CASE(public_chacha_seek_original_sets_high_counter_word) {
    using cipher_type = chacha<64, 256, 20>;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = {0};
    cipher_type::iv_type iv = {0};
    cipher_type cipher(block, schedule, key, iv);

    const std::uint64_t counter = std::uint64_t(1) << 32;
    schedule = original_chacha_state();
    cipher.seek(block, schedule, counter * cipher_type::block_size);

    schedule_type expected_state = original_chacha_state();
    expected_state[12] = 0;
    expected_state[13] = 1;
    const std::array<std::uint8_t, 64> expected_block = reference_block(expected_state);

    BOOST_TEST(std::equal(expected_block.begin(), expected_block.end(), block.begin()));
    BOOST_TEST(schedule[12] == 1u);
    BOOST_TEST(schedule[13] == 1u);
    BOOST_TEST(schedule[14] == original_chacha_state()[14]);
    BOOST_TEST(schedule[15] == original_chacha_state()[15]);
}

BOOST_AUTO_TEST_CASE(public_chacha_process_advances_across_consecutive_blocks) {
    using cipher_type = chacha<96, 256, 20>;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();
    cipher_type cipher(block, schedule, key, iv);

    cipher_type::block_type plaintext = {0};
    cipher_type::block_type first_ciphertext = {0};
    cipher_type::block_type second_ciphertext = {0};

    std::uint8_t *first_in = plaintext.data();
    std::uint8_t *first_out = first_ciphertext.data();
    cipher.process(first_in, first_out, schedule, block);

    std::uint8_t *second_in = plaintext.data();
    std::uint8_t *second_out = second_ciphertext.data();
    cipher.process(second_in, second_out, schedule, block);

    schedule_type expected_first_state = rfc8439_block_state();
    expected_first_state[12] = 0;
    const std::array<std::uint8_t, 64> expected_first_block = reference_block(expected_first_state);

    BOOST_TEST(std::equal(expected_first_block.begin(), expected_first_block.end(), first_ciphertext.begin()));
    BOOST_TEST(std::equal(rfc8439_expected_block.begin(), rfc8439_expected_block.end(), second_ciphertext.begin()));
    BOOST_TEST(schedule[12] == 3u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);
}

BOOST_DATA_TEST_CASE(chacha_single_range_encrypt, boost::unit_test::data::xrange(7), index) {
    using cipher_type = chacha<96, 256, 20>;

    cipher_type::block_type encrypt_block = {0};
    cipher_type::block_type decrypt_block = {0};
    cipher_type::key_schedule_type encrypt_schedule = {0};
    cipher_type::key_schedule_type decrypt_schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();

    cipher_type encrypt_cipher(encrypt_block, encrypt_schedule, key, iv);
    cipher_type decrypt_cipher(decrypt_block, decrypt_schedule, key, iv);

    cipher_type::block_type plaintext = {0};
    cipher_type::block_type ciphertext = {0};
    cipher_type::block_type decrypted = {0};

    for (std::size_t i = 0; i != plaintext.size(); ++i) {
        plaintext[i] = static_cast<std::uint8_t>((i * 3 + index) & 0xff);
    }

    std::uint8_t *encrypt_in = plaintext.data();
    std::uint8_t *encrypt_out = ciphertext.data();
    encrypt_cipher.process(encrypt_in, encrypt_out, encrypt_schedule, encrypt_block);

    std::uint8_t *decrypt_in = ciphertext.data();
    std::uint8_t *decrypt_out = decrypted.data();
    decrypt_cipher.process(decrypt_in, decrypt_out, decrypt_schedule, decrypt_block);

    BOOST_TEST(!std::equal(plaintext.begin(), plaintext.end(), ciphertext.begin()));
    BOOST_TEST(std::equal(plaintext.begin(), plaintext.end(), decrypted.begin()));
    BOOST_TEST(encrypt_schedule[12] == 2u);
    BOOST_TEST(decrypt_schedule[12] == 2u);
}

BOOST_AUTO_TEST_SUITE_END()
