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
#include <stdexcept>

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

    std::array<std::uint8_t, 12> rfc8439_encryption_iv() {
        return {{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x4a,
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

    const std::array<std::uint8_t, 114> rfc8439_sunscreen_plaintext = {
        {0x4c, 0x61, 0x64, 0x69, 0x65, 0x73, 0x20, 0x61, 0x6e, 0x64, 0x20, 0x47, 0x65, 0x6e, 0x74, 0x6c,
         0x65, 0x6d, 0x65, 0x6e, 0x20, 0x6f, 0x66, 0x20, 0x74, 0x68, 0x65, 0x20, 0x63, 0x6c, 0x61, 0x73,
         0x73, 0x20, 0x6f, 0x66, 0x20, 0x27, 0x39, 0x39, 0x3a, 0x20, 0x49, 0x66, 0x20, 0x49, 0x20, 0x63,
         0x6f, 0x75, 0x6c, 0x64, 0x20, 0x6f, 0x66, 0x66, 0x65, 0x72, 0x20, 0x79, 0x6f, 0x75, 0x20, 0x6f,
         0x6e, 0x6c, 0x79, 0x20, 0x6f, 0x6e, 0x65, 0x20, 0x74, 0x69, 0x70, 0x20, 0x66, 0x6f, 0x72, 0x20,
         0x74, 0x68, 0x65, 0x20, 0x66, 0x75, 0x74, 0x75, 0x72, 0x65, 0x2c, 0x20, 0x73, 0x75, 0x6e, 0x73,
         0x63, 0x72, 0x65, 0x65, 0x6e, 0x20, 0x77, 0x6f, 0x75, 0x6c, 0x64, 0x20, 0x62, 0x65, 0x20, 0x69,
         0x74, 0x2e}};

    const std::array<std::uint8_t, 114> rfc8439_sunscreen_ciphertext = {
        {0x6e, 0x2e, 0x35, 0x9a, 0x25, 0x68, 0xf9, 0x80, 0x41, 0xba, 0x07, 0x28, 0xdd, 0x0d, 0x69, 0x81,
         0xe9, 0x7e, 0x7a, 0xec, 0x1d, 0x43, 0x60, 0xc2, 0x0a, 0x27, 0xaf, 0xcc, 0xfd, 0x9f, 0xae, 0x0b,
         0xf9, 0x1b, 0x65, 0xc5, 0x52, 0x47, 0x33, 0xab, 0x8f, 0x59, 0x3d, 0xab, 0xcd, 0x62, 0xb3, 0x57,
         0x16, 0x39, 0xd6, 0x24, 0xe6, 0x51, 0x52, 0xab, 0x8f, 0x53, 0x0c, 0x35, 0x9f, 0x08, 0x61, 0xd8,
         0x07, 0xca, 0x0d, 0xbf, 0x50, 0x0d, 0x6a, 0x61, 0x56, 0xa3, 0x8e, 0x08, 0x8a, 0x22, 0xb6, 0x5e,
         0x52, 0xbc, 0x51, 0x4d, 0x16, 0xcc, 0xf8, 0x06, 0x81, 0x8c, 0xe9, 0x1a, 0xb7, 0x79, 0x37, 0x36,
         0x5a, 0xf9, 0x0b, 0xbf, 0x74, 0xa3, 0x5b, 0xe6, 0xb4, 0x0b, 0x8e, 0xed, 0xf2, 0x78, 0x5e, 0x42,
         0x87, 0x4d}};
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
    state[12] = 0xfffffffb;
    state[13] = 0x11223344;
    block4_type out = {0};

    impl_type::chacha_x4_ietf(out, state);

    BOOST_TEST(state[12] == 0xffffffffu);
    BOOST_TEST(state[13] == 0x11223344u);
}

BOOST_AUTO_TEST_CASE(ietf_chacha20_x4_rejects_counter_wrap) {
    using impl_type = detail::chacha_impl<20, 96, 256>;

    impl_type::key_schedule_type state = rfc8439_block_state();
    state[12] = 0xfffffffc;
    state[13] = 0x11223344;
    block4_type out = {0};

    BOOST_CHECK_THROW(impl_type::chacha_x4_ietf(out, state), std::out_of_range);
    BOOST_TEST(state[12] == 0xfffffffcu);
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

BOOST_AUTO_TEST_CASE(chacha_functions_rejects_ietf_counter_wrap) {
    using functions_type = detail::chacha_functions<20, 96, 256>;

    functions_type::key_schedule_type state = rfc8439_block_state();
    state[12] = 0xffffffff;
    functions_type::block_type out = {0};

    BOOST_CHECK_THROW(functions_type::generate_block(out, state), std::out_of_range);
    BOOST_TEST(state[12] == 0xffffffffu);
    BOOST_TEST(state[13] == rfc8439_block_state()[13]);
    BOOST_TEST(state[14] == rfc8439_block_state()[14]);
    BOOST_TEST(state[15] == rfc8439_block_state()[15]);
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
    BOOST_TEST(schedule[12] == 2u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);
}

BOOST_AUTO_TEST_CASE(public_chacha_facade_matches_rfc8439_encryption_vector) {
    using cipher_type = chacha<96, 256, 20>;

    cipher_type::block_type encrypt_block = {0};
    cipher_type::key_schedule_type encrypt_schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_encryption_iv();
    cipher_type encrypt_cipher(encrypt_block, encrypt_schedule, key, iv);

    encrypt_cipher.seek(encrypt_block, encrypt_schedule, cipher_type::block_size);

    std::array<std::uint8_t, rfc8439_sunscreen_plaintext.size()> ciphertext = {0};
    encrypt_cipher.process(rfc8439_sunscreen_plaintext.begin(), rfc8439_sunscreen_plaintext.end(),
                           ciphertext.begin(), encrypt_schedule, encrypt_block);

    BOOST_TEST(std::equal(rfc8439_sunscreen_ciphertext.begin(), rfc8439_sunscreen_ciphertext.end(),
                          ciphertext.begin()));
    BOOST_TEST(encrypt_schedule[12] == 3u);

    cipher_type::block_type decrypt_block = {0};
    cipher_type::key_schedule_type decrypt_schedule = {0};
    cipher_type decrypt_cipher(decrypt_block, decrypt_schedule, key, iv);

    decrypt_cipher.seek(decrypt_block, decrypt_schedule, cipher_type::block_size);

    std::array<std::uint8_t, rfc8439_sunscreen_ciphertext.size()> decrypted = {0};
    decrypt_cipher.process(rfc8439_sunscreen_ciphertext.begin(), rfc8439_sunscreen_ciphertext.end(),
                           decrypted.begin(), decrypt_schedule, decrypt_block);

    BOOST_TEST(std::equal(rfc8439_sunscreen_plaintext.begin(), rfc8439_sunscreen_plaintext.end(), decrypted.begin()));
    BOOST_TEST(decrypt_schedule[12] == 3u);
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

BOOST_AUTO_TEST_CASE(public_chacha_seek_ietf_rejects_counter_overflow) {
    using cipher_type = chacha<96, 256, 20>;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();
    cipher_type cipher(block, schedule, key, iv);

    const std::uint64_t offset = (std::uint64_t(1) << 32) * cipher_type::block_size;

    BOOST_CHECK_THROW(cipher.seek(block, schedule, offset), std::out_of_range);
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
    BOOST_TEST(schedule[12] == 0xffffffffu);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);

    cipher_type::block_type plaintext = {0};
    cipher_type::block_type ciphertext = {0};
    std::uint8_t *in = plaintext.data();
    std::uint8_t *out = ciphertext.data();

    cipher.process(in, out, schedule, block);

    BOOST_TEST(std::equal(expected_block.begin(), expected_block.end(), ciphertext.begin()));
    BOOST_TEST(schedule[12] == 0xffffffffu);

    std::array<std::uint8_t, 1> extra_plaintext = {0};
    std::array<std::uint8_t, 1> extra_ciphertext = {0};

    BOOST_CHECK_THROW(cipher.process(extra_plaintext.begin(), extra_plaintext.end(), extra_ciphertext.begin(), schedule,
                                     block),
                      std::out_of_range);
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
    BOOST_TEST(schedule[12] == 2u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);
}

BOOST_AUTO_TEST_CASE(public_chacha_process_handles_partial_range_and_resume) {
    using cipher_type = chacha<96, 256, 20>;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();
    cipher_type cipher(block, schedule, key, iv);

    std::array<std::uint8_t, 70> plaintext = {0};
    std::array<std::uint8_t, 70> ciphertext = {0};
    std::array<std::uint8_t, 70> expected_ciphertext = {0};

    for (std::size_t i = 0; i != plaintext.size(); ++i) {
        plaintext[i] = static_cast<std::uint8_t>((i * 5 + 1) & 0xff);
    }

    schedule_type first_state = rfc8439_block_state();
    first_state[12] = 0;
    const std::array<std::uint8_t, 64> first_block = reference_block(first_state);
    const std::array<std::uint8_t, 64> second_block = rfc8439_expected_block;

    for (std::size_t i = 0; i != 64; ++i) {
        expected_ciphertext[i] = plaintext[i] ^ first_block[i];
    }
    for (std::size_t i = 64; i != plaintext.size(); ++i) {
        expected_ciphertext[i] = plaintext[i] ^ second_block[i - 64];
    }

    cipher.process(plaintext.begin(), plaintext.end(), ciphertext.begin(), schedule, block);

    BOOST_TEST(std::equal(expected_ciphertext.begin(), expected_ciphertext.end(), ciphertext.begin()));
    BOOST_TEST(schedule[12] == 2u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);

    std::array<std::uint8_t, 58> continuation = {0};
    std::array<std::uint8_t, 58> continuation_ciphertext = {0};
    std::array<std::uint8_t, 58> expected_continuation = {0};

    for (std::size_t i = 0; i != continuation.size(); ++i) {
        continuation[i] = static_cast<std::uint8_t>((i * 7 + 3) & 0xff);
        expected_continuation[i] = continuation[i] ^ second_block[i + 6];
    }

    cipher.process(continuation.begin(), continuation.end(), continuation_ciphertext.begin(), schedule, block);

    BOOST_TEST(std::equal(expected_continuation.begin(), expected_continuation.end(), continuation_ciphertext.begin()));
    BOOST_TEST(schedule[12] == 2u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);
}

BOOST_AUTO_TEST_CASE(public_chacha_seek_supports_byte_offsets_inside_block) {
    using cipher_type = chacha<96, 256, 20>;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();
    cipher_type cipher(block, schedule, key, iv);

    const std::uint64_t offset = cipher_type::block_size + 7;
    cipher.seek(block, schedule, offset);

    std::array<std::uint8_t, 10> plaintext = {0};
    std::array<std::uint8_t, 10> ciphertext = {0};
    std::array<std::uint8_t, 10> expected_ciphertext = {0};

    for (std::size_t i = 0; i != plaintext.size(); ++i) {
        plaintext[i] = static_cast<std::uint8_t>(0xa0 + i);
        expected_ciphertext[i] = plaintext[i] ^ rfc8439_expected_block[i + 7];
    }

    cipher.process(plaintext.begin(), plaintext.end(), ciphertext.begin(), schedule, block);

    BOOST_TEST(std::equal(expected_ciphertext.begin(), expected_ciphertext.end(), ciphertext.begin()));
    BOOST_TEST(schedule[12] == 2u);
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
    BOOST_TEST(encrypt_schedule[12] == 1u);
    BOOST_TEST(decrypt_schedule[12] == 1u);
}

BOOST_AUTO_TEST_SUITE_END()
