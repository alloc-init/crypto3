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
#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

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

    std::uint8_t hex_value(char c) {
        if (c >= '0' && c <= '9') {
            return static_cast<std::uint8_t>(c - '0');
        }
        if (c >= 'a' && c <= 'f') {
            return static_cast<std::uint8_t>(c - 'a' + 10);
        }
        if (c >= 'A' && c <= 'F') {
            return static_cast<std::uint8_t>(c - 'A' + 10);
        }
        throw std::invalid_argument("invalid hex digit");
    }

    std::vector<std::uint8_t> hex_to_bytes(const char *hex) {
        std::vector<std::uint8_t> bytes;
        int high_nibble = -1;

        for (const char *p = hex; *p != '\0'; ++p) {
            const unsigned char c = static_cast<unsigned char>(*p);
            if (std::isspace(c)) {
                continue;
            }

            const std::uint8_t value = hex_value(*p);
            if (high_nibble < 0) {
                high_nibble = value;
            } else {
                bytes.push_back(static_cast<std::uint8_t>((high_nibble << 4) | value));
                high_nibble = -1;
            }
        }

        if (high_nibble >= 0) {
            throw std::invalid_argument("odd number of hex digits");
        }

        return bytes;
    }

    template<std::size_t Size>
    std::array<std::uint8_t, Size> hex_to_array(const char *hex) {
        const std::vector<std::uint8_t> bytes = hex_to_bytes(hex);
        if (bytes.size() != Size) {
            throw std::invalid_argument("unexpected vector size");
        }

        std::array<std::uint8_t, Size> out = {0};
        std::copy(bytes.begin(), bytes.end(), out.begin());
        return out;
    }

    std::vector<std::uint8_t> chacha20_encrypt(const chacha20::key_type &key, const chacha20::iv_type &iv,
                                               std::uint32_t counter,
                                               const std::vector<std::uint8_t> &plaintext) {
        std::vector<std::uint8_t> ciphertext(plaintext.size());
        nil::crypto3::encrypt<chacha20>(plaintext.begin(), plaintext.end(), key, iv, ciphertext.begin(), counter);
        return ciphertext;
    }

    bool equal_block(const block4_type &actual, std::size_t block_index,
                     const std::array<std::uint8_t, 64> &expected) {
        return std::equal(expected.begin(), expected.end(), actual.begin() + block_index * 64);
    }

    std::uint32_t rotl(std::uint32_t x, std::size_t n) {
        return static_cast<std::uint32_t>((x << n) | (x >> (32 - n)));
    }

    void reference_quarter_round(std::uint32_t &a, std::uint32_t &b, std::uint32_t &c, std::uint32_t &d) {
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
            reference_quarter_round(x[0], x[4], x[8], x[12]);
            reference_quarter_round(x[1], x[5], x[9], x[13]);
            reference_quarter_round(x[2], x[6], x[10], x[14]);
            reference_quarter_round(x[3], x[7], x[11], x[15]);

            reference_quarter_round(x[0], x[5], x[10], x[15]);
            reference_quarter_round(x[1], x[6], x[11], x[12]);
            reference_quarter_round(x[2], x[7], x[8], x[13]);
            reference_quarter_round(x[3], x[4], x[9], x[14]);
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

BOOST_AUTO_TEST_CASE(public_chacha_aliases_document_standard_variants) {
    BOOST_TEST((std::is_same<chacha20, chacha<96, 256, 20>>::value));
    BOOST_TEST((std::is_same<original_chacha20, chacha<64, 256, 20>>::value));
    BOOST_TEST((std::is_same<ietf_chacha<128, 12>, chacha<96, 128, 12>>::value));
    BOOST_TEST((std::is_same<original_chacha<128, 8>, chacha<64, 128, 8>>::value));
}

BOOST_AUTO_TEST_CASE(chacha_quarter_round_matches_rfc8439_vector) {
    using impl_type = detail::chacha_impl<20, 96, 256>;

    std::uint32_t a = 0x11111111;
    std::uint32_t b = 0x01020304;
    std::uint32_t c = 0x9b8d6f43;
    std::uint32_t d = 0x01234567;

    impl_type::quarter_round(a, b, c, d);

    BOOST_TEST(a == 0xea2a92f4u);
    BOOST_TEST(b == 0xcb1cf8ceu);
    BOOST_TEST(c == 0x4581472eu);
    BOOST_TEST(d == 0x5881c4bbu);
}

BOOST_AUTO_TEST_CASE(chacha_quarter_round_on_state_matches_rfc8439_vector) {
    using impl_type = detail::chacha_impl<20, 96, 256>;

    schedule_type state = {{0x879531e0, 0xc5ecf37d, 0x516461b1, 0xc9a62f8a,
                            0x44c20ef3, 0x3390af7f, 0xd9fc690b, 0x2a5f714c,
                            0x53372767, 0xb00a5631, 0x974c541a, 0x359e9963,
                            0x5c971061, 0x3d631689, 0x2098d9d6, 0x91dbd320}};
    const schedule_type expected = {{0x879531e0, 0xc5ecf37d, 0xbdb886dc, 0xc9a62f8a,
                                     0x44c20ef3, 0x3390af7f, 0xd9fc690b, 0xcfacafd2,
                                     0xe46bea80, 0xb00a5631, 0x974c541a, 0x359e9963,
                                     0x5c971061, 0xccc07c79, 0x2098d9d6, 0x91dbd320}};

    impl_type::quarter_round(state[2], state[7], state[8], state[13]);

    BOOST_TEST(std::equal(expected.begin(), expected.end(), state.begin()));
}

BOOST_AUTO_TEST_CASE(chacha_unimplemented_simd_path_throws) {
    using simd_impl_type = detail::chacha_unimplemented_simd_impl<20, 96, 256>;
    using avx2_impl_type = detail::chacha_unimplemented_avx2_impl<20, 96, 256>;
    using sse2_impl_type = detail::chacha_unimplemented_sse2_impl<20, 96, 256>;

    simd_impl_type::key_schedule_type state = rfc8439_block_state();
    block4_type out = {0};

    BOOST_CHECK_EXCEPTION(simd_impl_type::chacha_x4(out, state), std::logic_error, [](const std::logic_error &e) {
        return std::string(e.what()) == "ChaCha SIMD implementation is not implemented";
    });
    BOOST_CHECK_EXCEPTION(avx2_impl_type::chacha_x4(out, state), std::logic_error, [](const std::logic_error &e) {
        return std::string(e.what()) == "ChaCha AVX2 implementation is not implemented";
    });
    BOOST_CHECK_EXCEPTION(sse2_impl_type::chacha_x4_ietf(out, state), std::logic_error, [](const std::logic_error &e) {
        return std::string(e.what()) == "ChaCha SSE2 implementation is not implemented";
    });
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

BOOST_AUTO_TEST_CASE(chacha_functions_schedule_iv_uses_ietf_initial_counter) {
    using functions_type = detail::chacha_functions<20, 96, 256>;

    functions_type::key_schedule_type schedule = {0};
    functions_type::block_type block = {0};
    functions_type::key_type key = rfc8439_key();
    functions_type::iv_type iv = rfc8439_iv();

    functions_type::schedule_key(schedule, key);
    functions_type::schedule_iv(block, schedule, iv, 1);

    BOOST_TEST(std::equal(rfc8439_expected_block.begin(), rfc8439_expected_block.end(), block.begin()));
    BOOST_TEST(schedule[12] == 2u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);
}

BOOST_AUTO_TEST_CASE(chacha_functions_schedule_iv_rejects_ietf_initial_counter_overflow) {
    using functions_type = detail::chacha_functions<20, 96, 256>;

    functions_type::key_schedule_type schedule = {0};
    functions_type::block_type block = {0};
    functions_type::key_type key = rfc8439_key();
    functions_type::iv_type iv = rfc8439_iv();

    functions_type::schedule_key(schedule, key);

    BOOST_CHECK_THROW(functions_type::schedule_iv(block, schedule, iv, std::uint64_t(1) << 32), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(chacha_functions_schedule_iv_accepts_maximum_ietf_initial_counter) {
    {
        using functions_type = detail::chacha_functions<20, 96, 128>;

        functions_type::key_schedule_type expected_schedule = {0};
        functions_type::key_schedule_type schedule = {0};
        functions_type::block_type block = {0};
        functions_type::key_type key = {0};
        functions_type::iv_type iv = rfc8439_iv();

        functions_type::schedule_key(expected_schedule, key);
        functions_type::schedule_iv(expected_schedule, iv, 0xffffffffULL);

        functions_type::schedule_key(schedule, key);
        functions_type::schedule_iv(block, schedule, iv, 0xffffffffULL);

        const std::array<std::uint8_t, 64> expected_block = reference_block(expected_schedule);

        BOOST_TEST(std::equal(expected_block.begin(), expected_block.end(), block.begin()));
        BOOST_TEST(std::equal(expected_schedule.begin(), expected_schedule.end(), schedule.begin()));
    }

    {
        using functions_type = detail::chacha_functions<20, 96, 256>;

        functions_type::key_schedule_type expected_schedule = {0};
        functions_type::key_schedule_type schedule = {0};
        functions_type::block_type block = {0};
        functions_type::key_type key = rfc8439_key();
        functions_type::iv_type iv = rfc8439_iv();

        functions_type::schedule_key(expected_schedule, key);
        functions_type::schedule_iv(expected_schedule, iv, 0xffffffffULL);

        functions_type::schedule_key(schedule, key);
        functions_type::schedule_iv(block, schedule, iv, 0xffffffffULL);

        const std::array<std::uint8_t, 64> expected_block = reference_block(expected_schedule);

        BOOST_TEST(std::equal(expected_block.begin(), expected_block.end(), block.begin()));
        BOOST_TEST(std::equal(expected_schedule.begin(), expected_schedule.end(), schedule.begin()));
    }
}

BOOST_AUTO_TEST_CASE(chacha_functions_schedule_iv_splits_original_initial_counter) {
    using functions_type = detail::chacha_functions<20, 64, 256>;

    functions_type::key_schedule_type schedule = original_chacha_state();
    functions_type::block_type block = {0};
    functions_type::iv_type iv = {{0x00, 0x00, 0x00, 0x4a, 0x00, 0x00, 0x00, 0x00}};
    const std::uint64_t counter = std::uint64_t(1) << 32;

    functions_type::schedule_iv(block, schedule, iv, counter);

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

BOOST_AUTO_TEST_CASE(public_chacha_constructor_uses_schedule_iv_block_argument) {
    using cipher_type = chacha20;

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

BOOST_AUTO_TEST_CASE(public_chacha_constructor_uses_initial_counter) {
    using cipher_type = chacha20;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();

    cipher_type cipher(block, schedule, key, iv, 1);

    (void)cipher;
    BOOST_TEST(std::equal(rfc8439_expected_block.begin(), rfc8439_expected_block.end(), block.begin()));
    BOOST_TEST(schedule[12] == 2u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);
}

BOOST_AUTO_TEST_CASE(public_chacha_constructor_accepts_maximum_ietf_initial_counter) {
    using cipher_type = chacha20;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();
    cipher_type cipher(block, schedule, key, iv, 0xffffffffULL);

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

BOOST_AUTO_TEST_CASE(public_chacha_facade_matches_rfc8439_block_vector) {
    using cipher_type = chacha20;

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

    std::uint8_t *returned_out = cipher.process(in, out, schedule, block);

    BOOST_TEST(std::equal(rfc8439_expected_block.begin(), rfc8439_expected_block.end(), ciphertext.begin()));
    BOOST_TEST(returned_out == ciphertext.data() + ciphertext.size());
    BOOST_TEST(schedule[12] == 2u);
    BOOST_TEST(schedule[13] == rfc8439_block_state()[13]);
    BOOST_TEST(schedule[14] == rfc8439_block_state()[14]);
    BOOST_TEST(schedule[15] == rfc8439_block_state()[15]);
}

BOOST_AUTO_TEST_CASE(public_chacha_facade_matches_rfc8439_encryption_vector) {
    using cipher_type = chacha20;

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

BOOST_AUTO_TEST_CASE(public_chacha20_cipher_is_not_copyable) {
    BOOST_TEST(!std::is_copy_constructible<chacha20_cipher>::value);
    BOOST_TEST(!std::is_copy_assignable<chacha20_cipher>::value);
}

BOOST_AUTO_TEST_CASE(public_chacha20_cipher_matches_rfc8439_encryption_vector) {
    chacha20::key_type key = rfc8439_key();
    chacha20::iv_type iv = rfc8439_encryption_iv();

    chacha20_cipher encrypt_cipher(key, iv, 1);
    std::array<std::uint8_t, rfc8439_sunscreen_plaintext.size()> ciphertext = {0};

    const auto encrypt_out = encrypt_cipher.process(rfc8439_sunscreen_plaintext.begin(),
                                                    rfc8439_sunscreen_plaintext.end(), ciphertext.begin());

    BOOST_TEST((encrypt_out == ciphertext.end()));
    BOOST_TEST(std::equal(rfc8439_sunscreen_ciphertext.begin(), rfc8439_sunscreen_ciphertext.end(),
                          ciphertext.begin()));

    chacha20_cipher decrypt_cipher(key, iv, 1);
    std::array<std::uint8_t, rfc8439_sunscreen_ciphertext.size()> decrypted = {0};

    const auto decrypt_out = decrypt_cipher.process(rfc8439_sunscreen_ciphertext.begin(),
                                                    rfc8439_sunscreen_ciphertext.end(), decrypted.begin());

    BOOST_TEST((decrypt_out == decrypted.end()));
    BOOST_TEST(std::equal(rfc8439_sunscreen_plaintext.begin(), rfc8439_sunscreen_plaintext.end(), decrypted.begin()));
}

BOOST_AUTO_TEST_CASE(public_chacha20_cipher_processes_resumable_byte_streams) {
    chacha20::key_type key = rfc8439_key();
    chacha20::iv_type iv = rfc8439_iv();

    std::vector<std::uint8_t> plaintext(129);
    for (std::size_t i = 0; i != plaintext.size(); ++i) {
        plaintext[i] = static_cast<std::uint8_t>((i * 11 + 5) & 0xff);
    }

    std::vector<std::uint8_t> expected(plaintext.size());
    nil::crypto3::encrypt<chacha20>(plaintext.begin(), plaintext.end(), key, iv, expected.begin(), 1);

    chacha20_cipher cipher(key, iv, 1);
    std::vector<std::uint8_t> ciphertext(plaintext.size());

    const std::size_t first_chunk = 70;
    auto out = cipher.process(plaintext.begin(), plaintext.begin() + first_chunk, ciphertext.begin());
    out = cipher.process(plaintext.begin() + first_chunk, plaintext.end(), out);

    BOOST_TEST((out == ciphertext.end()));
    BOOST_TEST(std::equal(expected.begin(), expected.end(), ciphertext.begin()));
}

BOOST_AUTO_TEST_CASE(public_chacha20_cipher_seek_is_relative_to_initial_counter) {
    chacha20::key_type key = rfc8439_key();
    chacha20::iv_type iv = rfc8439_iv();
    chacha20_cipher cipher(key, iv, 1);

    const std::uint64_t offset = 7;
    cipher.seek(offset);

    std::array<std::uint8_t, 10> plaintext = {0};
    std::array<std::uint8_t, 10> ciphertext = {0};
    std::array<std::uint8_t, 10> expected_ciphertext = {0};

    for (std::size_t i = 0; i != plaintext.size(); ++i) {
        plaintext[i] = static_cast<std::uint8_t>(0xa0 + i);
        expected_ciphertext[i] = plaintext[i] ^ rfc8439_expected_block[i + offset];
    }

    cipher.process(plaintext.begin(), plaintext.end(), ciphertext.begin());

    BOOST_TEST(std::equal(expected_ciphertext.begin(), expected_ciphertext.end(), ciphertext.begin()));
}

BOOST_AUTO_TEST_CASE(public_chacha20_cipher_supports_in_place_encryption) {
    chacha20::key_type key = rfc8439_key();
    chacha20::iv_type iv = rfc8439_iv();

    std::vector<std::uint8_t> buffer(130);
    for (std::size_t i = 0; i != buffer.size(); ++i) {
        buffer[i] = static_cast<std::uint8_t>((i * 13 + 9) & 0xff);
    }

    const std::vector<std::uint8_t> plaintext = buffer;
    std::vector<std::uint8_t> expected(buffer.size());
    nil::crypto3::encrypt<chacha20>(plaintext.begin(), plaintext.end(), key, iv, expected.begin(), 0);

    nil::crypto3::encrypt<chacha20>(buffer.begin(), buffer.end(), key, iv, buffer.begin(), 0);

    BOOST_TEST(std::equal(expected.begin(), expected.end(), buffer.begin()));

    nil::crypto3::decrypt<chacha20>(buffer.begin(), buffer.end(), key, iv, buffer.begin(), 0);

    BOOST_TEST(std::equal(plaintext.begin(), plaintext.end(), buffer.begin()));
}

BOOST_AUTO_TEST_CASE(public_chacha20_encrypt_decrypt_helpers_match_rfc8439_vector) {
    chacha20::key_type key = rfc8439_key();
    chacha20::iv_type iv = rfc8439_encryption_iv();

    std::array<std::uint8_t, rfc8439_sunscreen_plaintext.size()> ciphertext = {0};
    std::array<std::uint8_t, rfc8439_sunscreen_ciphertext.size()> decrypted = {0};
    std::array<std::uint8_t, rfc8439_sunscreen_plaintext.size()> ciphertext_without_explicit_cipher = {0};

    nil::crypto3::encrypt<chacha20>(rfc8439_sunscreen_plaintext.begin(), rfc8439_sunscreen_plaintext.end(), key, iv,
                                    ciphertext.begin(), 1);
    nil::crypto3::decrypt<chacha20>(ciphertext.begin(), ciphertext.end(), key, iv, decrypted.begin(), 1);
    nil::crypto3::encrypt(rfc8439_sunscreen_plaintext.begin(), rfc8439_sunscreen_plaintext.end(), key, iv,
                          ciphertext_without_explicit_cipher.begin(), 1);

    BOOST_TEST(std::equal(rfc8439_sunscreen_ciphertext.begin(), rfc8439_sunscreen_ciphertext.end(),
                          ciphertext.begin()));
    BOOST_TEST(std::equal(rfc8439_sunscreen_plaintext.begin(), rfc8439_sunscreen_plaintext.end(), decrypted.begin()));
    BOOST_TEST(std::equal(ciphertext.begin(), ciphertext.end(), ciphertext_without_explicit_cipher.begin()));
}

BOOST_AUTO_TEST_CASE(public_chacha20_cipher_rejects_ietf_counter_wrap) {
    chacha20::key_type key = rfc8439_key();
    chacha20::iv_type iv = rfc8439_iv();
    chacha20_cipher cipher(key, iv, 0xffffffff);

    chacha20::block_type plaintext = {0};
    chacha20::block_type ciphertext = {0};

    cipher.process(plaintext.begin(), plaintext.end(), ciphertext.begin());

    std::array<std::uint8_t, 1> extra_plaintext = {0};
    std::array<std::uint8_t, 1> extra_ciphertext = {0};
    std::array<std::uint8_t, chacha20::block_size + 1> too_long_plaintext = {0};
    std::array<std::uint8_t, chacha20::block_size + 1> too_long_ciphertext = {0};

    BOOST_CHECK_THROW(cipher.process(extra_plaintext.begin(), extra_plaintext.end(), extra_ciphertext.begin()),
                      std::out_of_range);
    BOOST_CHECK_THROW(cipher.seek(chacha20::block_size), std::out_of_range);
    BOOST_CHECK_THROW(nil::crypto3::encrypt<chacha20>(too_long_plaintext.begin(), too_long_plaintext.end(), key, iv,
                                                      too_long_ciphertext.begin(), 0xffffffff),
                      std::out_of_range);
}

BOOST_AUTO_TEST_CASE(public_chacha20_matches_rfc8439_additional_block_vectors) {
    struct block_vector {
        const char *key;
        const char *iv;
        std::uint32_t counter;
        const char *block;
    };

    const block_vector vectors[] = {
        {"00000000000000000000000000000000"
         "00000000000000000000000000000000",
         "000000000000000000000000",
         0,
         "76b8e0ada0f13d90405d6ae55386bd28"
         "bdd219b8a08ded1aa836efcc8b770dc7"
         "da41597c5157488d7724e03fb8d84a37"
         "6a43b8f41518a11cc387b669b2ee6586"},
        {"00000000000000000000000000000000"
         "00000000000000000000000000000000",
         "000000000000000000000000",
         1,
         "9f07e7be5551387a98ba977c732d080d"
         "cb0f29a048e3656912c6533e32ee7aed"
         "29b721769ce64e43d57133b074d839d5"
         "31ed1f28510afb45ace10a1f4b794d6f"},
        {"00000000000000000000000000000000"
         "00000000000000000000000000000001",
         "000000000000000000000000",
         1,
         "3aeb5224ecf849929b9d828db1ced4dd"
         "832025e8018b8160b82284f3c949aa5a"
         "8eca00bbb4a73bdad192b5c42f73f2fd"
         "4e273644c8b36125a64addeb006c13a0"},
        {"00ff0000000000000000000000000000"
         "00000000000000000000000000000000",
         "000000000000000000000000",
         2,
         "72d54dfbf12ec44b362692df94137f32"
         "8fea8da73990265ec1bbbea1ae9af0ca"
         "13b25aa26cb4a648cb9b9d1be65b2c09"
         "24a66c54d545ec1b7374f4872e99f096"},
        {"00000000000000000000000000000000"
         "00000000000000000000000000000000",
         "000000000000000000000002",
         0,
         "c2c64d378cd536374ae204b9ef933fcd"
         "1a8b2288b3dfa49672ab765b54ee27c7"
         "8a970e0e955c14f3a88e741b97c286f7"
         "5f8fc299e8148362fa198a39531bed6d"}};

    const std::vector<std::uint8_t> plaintext(64, 0);

    for (std::size_t i = 0; i != sizeof(vectors) / sizeof(vectors[0]); ++i) {
        const block_vector &vector = vectors[i];
        const chacha20::key_type key = hex_to_array<32>(vector.key);
        const chacha20::iv_type iv = hex_to_array<12>(vector.iv);
        const std::vector<std::uint8_t> expected = hex_to_bytes(vector.block);

        const std::vector<std::uint8_t> ciphertext = chacha20_encrypt(key, iv, vector.counter, plaintext);

        BOOST_TEST_CONTEXT("RFC 8439 Appendix A.1 vector " << (i + 1)) {
            BOOST_REQUIRE(ciphertext.size() == expected.size());
            BOOST_TEST(std::equal(expected.begin(), expected.end(), ciphertext.begin()));
        }
    }
}

BOOST_AUTO_TEST_CASE(public_chacha20_matches_rfc8439_additional_encryption_vectors) {
    struct encryption_vector {
        const char *key;
        const char *iv;
        std::uint32_t counter;
        const char *plaintext;
        const char *ciphertext;
    };

    const encryption_vector vectors[] = {
        {"00000000000000000000000000000000"
         "00000000000000000000000000000000",
         "000000000000000000000000",
         0,
         "00000000000000000000000000000000"
         "00000000000000000000000000000000"
         "00000000000000000000000000000000"
         "00000000000000000000000000000000",
         "76b8e0ada0f13d90405d6ae55386bd28"
         "bdd219b8a08ded1aa836efcc8b770dc7"
         "da41597c5157488d7724e03fb8d84a37"
         "6a43b8f41518a11cc387b669b2ee6586"},
        {"00000000000000000000000000000000"
         "00000000000000000000000000000001",
         "000000000000000000000002",
         1,
         "416e79207375626d697373696f6e2074"
         "6f20746865204945544620696e74656e"
         "6465642062792074686520436f6e7472"
         "696275746f7220666f72207075626c69"
         "636174696f6e20617320616c6c206f72"
         "2070617274206f6620616e2049455446"
         "20496e7465726e65742d447261667420"
         "6f722052464320616e6420616e792073"
         "746174656d656e74206d616465207769"
         "7468696e2074686520636f6e74657874"
         "206f6620616e20494554462061637469"
         "7669747920697320636f6e7369646572"
         "656420616e20224945544620436f6e74"
         "7269627574696f6e222e205375636820"
         "73746174656d656e747320696e636c75"
         "6465206f72616c2073746174656d656e"
         "747320696e20494554462073657373696f6e73"
         "2c2061732077656c6c20617320777269"
         "7474656e20616e6420656c656374726f"
         "6e696320636f6d6d756e69636174696f"
         "6e73206d61646520617420616e792074"
         "696d65206f7220706c6163652c207768"
         "69636820617265206164647265737365"
         "6420746f",
         "a3fbf07df3fa2fde4f376ca23e827370"
         "41605d9f4f4f57bd8cff2c1d4b7955ec"
         "2a97948bd3722915c8f3d337f7d37005"
         "0e9e96d647b7c39f56e031ca5eb6250d"
         "4042e02785ececfa4b4bb5e8ead0440e"
         "20b6e8db09d881a7c6132f420e527950"
         "42bdfa7773d8a9051447b3291ce1411c"
         "680465552aa6c405b7764d5e87bea85a"
         "d00f8449ed8f72d0d662ab052691ca66"
         "424bc86d2df80ea41f43abf937d3259d"
         "c4b2d0dfb48a6c9139ddd7f76966e928"
         "e635553ba76c5c879d7b35d49eb2e62b"
         "0871cdac638939e25e8a1e0ef9d5280f"
         "a8ca328b351c3c765989cbcf3daa8b6c"
         "cc3aaf9f3979c92b3720fc88dc95ed84"
         "a1be059c6499b9fda236e7e818b04b0b"
         "c39c1e876b193bfe5569753f88128cc0"
         "8aaa9b63d1a16f80ef2554d7189c411f"
         "5869ca52c5b83fa36ff216b9c1d30062"
         "bebcfd2dc5bce0911934fda79a86f6e6"
         "98ced759c3ff9b6477338f3da4f9cd85"
         "14ea9982ccafb341b2384dd902f3d1ab"
         "7ac61dd29c6f21ba5b862f3730e37cfd"
         "c4fd806c22f221"},
        {"1c9240a5eb55d38af333888604f6b5f0"
         "473917c1402b80099dca5cbc207075c0",
         "000000000000000000000002",
         42,
         "2754776173206272696c6c69672c2061"
         "6e642074686520736c6974687920746f"
         "7665730a446964206779726520616e64"
         "2067696d626c6520696e207468652077"
         "6162653a0a416c6c206d696d737920"
         "776572652074686520626f726f676f76"
         "65732c0a416e6420746865206d6f6d"
         "65207261746873206f75746772616265"
         "2e",
         "62e6347f95ed87a45ffae7426f27a1df"
         "5fb69110044c0d73118effa95b01e5cf"
         "166d3df2d721caf9b21e5fb14c616871"
         "fd84c54f9d65b283196c7fe4f60553eb"
         "f39c6402c42234e32a356b3e764312a6"
         "1a5532055716ead6962568f87d3f3f77"
         "04c6a8d1bcd1bf4d50d6154b6da731b1"
         "87b58dfd728afa36757a797ac188d1"}};

    for (std::size_t i = 0; i != sizeof(vectors) / sizeof(vectors[0]); ++i) {
        const encryption_vector &vector = vectors[i];
        const chacha20::key_type key = hex_to_array<32>(vector.key);
        const chacha20::iv_type iv = hex_to_array<12>(vector.iv);
        const std::vector<std::uint8_t> plaintext = hex_to_bytes(vector.plaintext);
        const std::vector<std::uint8_t> expected_ciphertext = hex_to_bytes(vector.ciphertext);

        const std::vector<std::uint8_t> ciphertext = chacha20_encrypt(key, iv, vector.counter, plaintext);
        const std::vector<std::uint8_t> decrypted = chacha20_encrypt(key, iv, vector.counter, ciphertext);

        BOOST_TEST_CONTEXT("RFC 8439 Appendix A.2 vector " << (i + 1)) {
            BOOST_REQUIRE(ciphertext.size() == expected_ciphertext.size());
            BOOST_TEST(std::equal(expected_ciphertext.begin(), expected_ciphertext.end(), ciphertext.begin()));
            BOOST_TEST(std::equal(plaintext.begin(), plaintext.end(), decrypted.begin()));
        }
    }
}

BOOST_AUTO_TEST_CASE(public_chacha_seek_ietf_sets_32_bit_counter_and_preserves_nonce) {
    using cipher_type = chacha20;

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
    using cipher_type = chacha20;

    cipher_type::block_type block = {0};
    cipher_type::key_schedule_type schedule = {0};
    cipher_type::key_type key = rfc8439_key();
    cipher_type::iv_type iv = rfc8439_iv();
    cipher_type cipher(block, schedule, key, iv);

    const std::uint64_t offset = (std::uint64_t(1) << 32) * cipher_type::block_size;

    BOOST_CHECK_THROW(cipher.seek(block, schedule, offset), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(public_chacha_seek_ietf_accepts_maximum_counter_block) {
    using cipher_type = chacha20;

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
    using cipher_type = original_chacha20;

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
    using cipher_type = original_chacha20;

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
    using cipher_type = chacha20;

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
    using cipher_type = chacha20;

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
    using cipher_type = chacha20;

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
    using cipher_type = chacha20;

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
