//---------------------------------------------------------------------------//
//
// Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
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

#define BOOST_TEST_MODULE aead_chacha20poly1305_test

#include <nil/crypto3/modes/aead/chacha20poly1305.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nil::crypto3;

namespace {
    std::vector<std::uint8_t> hex_to_bytes(const std::string &hex) {
        std::vector<std::uint8_t> out;
        int high_nibble = -1;

        for (char c : hex) {
            if (std::isspace(static_cast<unsigned char>(c))) {
                continue;
            }

            int value = -1;
            if (c >= '0' && c <= '9') {
                value = c - '0';
            } else if (c >= 'a' && c <= 'f') {
                value = c - 'a' + 10;
            } else if (c >= 'A' && c <= 'F') {
                value = c - 'A' + 10;
            } else {
                throw std::invalid_argument("invalid hex character");
            }

            if (high_nibble < 0) {
                high_nibble = value;
            } else {
                out.push_back(static_cast<std::uint8_t>((high_nibble << 4) | value));
                high_nibble = -1;
            }
        }

        if (high_nibble >= 0) {
            throw std::invalid_argument("odd hex length");
        }
        return out;
    }

    template<std::size_t Size>
    std::array<std::uint8_t, Size> array_from_hex(const std::string &hex) {
        const std::vector<std::uint8_t> bytes = hex_to_bytes(hex);
        if (bytes.size() != Size) {
            throw std::invalid_argument("hex vector has unexpected size");
        }

        std::array<std::uint8_t, Size> out = {0};
        std::copy(bytes.begin(), bytes.end(), out.begin());
        return out;
    }

    std::vector<std::uint8_t> bytes_from_text(const std::string &text) {
        return std::vector<std::uint8_t>(text.begin(), text.end());
    }

    typedef modes::aead::chacha20poly1305 aead_type;

    aead_type::key_type rfc8439_aead_key() {
        return array_from_hex<aead_type::key_size>(
            "808182838485868788898a8b8c8d8e8f"
            "909192939495969798999a9b9c9d9e9f");
    }

    aead_type::nonce_type rfc8439_aead_nonce() {
        return array_from_hex<aead_type::nonce_size>("070000004041424344454647");
    }

    std::vector<std::uint8_t> rfc8439_aead_aad() {
        return hex_to_bytes("50515253c0c1c2c3c4c5c6c7");
    }

    std::vector<std::uint8_t> rfc8439_aead_plaintext() {
        return bytes_from_text(
            "Ladies and Gentlemen of the class of '99: If I could offer you only one tip for the future, sunscreen "
            "would be it.");
    }

    std::vector<std::uint8_t> rfc8439_aead_ciphertext() {
        return hex_to_bytes(
            "d31a8d34648e60db7b86afbc53ef7ec2"
            "a4aded51296e08fea9e2b5a736ee62d6"
            "3dbea45e8ca9671282fafb69da92728b"
            "1a71de0a9e060b2905d6a5b67ecd3b36"
            "92ddbd7f2d778b8c9803aee328091b58"
            "fab324e4fad675945585808b4831d7bc"
            "3ff4def08e4b7a9de576d26586cec64b"
            "6116");
    }

    aead_type::tag_type rfc8439_aead_tag() {
        return array_from_hex<aead_type::tag_size>("1ae10b594f09e26a7e902ecbd0600691");
    }

    aead_type::key_type rfc8439_decryption_key() {
        return array_from_hex<aead_type::key_size>(
            "1c9240a5eb55d38af333888604f6b5f0"
            "473917c1402b80099dca5cbc207075c0");
    }

    aead_type::nonce_type rfc8439_decryption_nonce() {
        return array_from_hex<aead_type::nonce_size>("000000000102030405060708");
    }

    std::vector<std::uint8_t> rfc8439_decryption_aad() {
        return hex_to_bytes("f33388860000000000004e91");
    }

    std::vector<std::uint8_t> rfc8439_decryption_ciphertext() {
        return hex_to_bytes(
            "64a0861575861af460f062c79be643bd"
            "5e805cfd345cf389f108670ac76c8cb2"
            "4c6cfc18755d43eea09ee94e382d26b0"
            "bdb7b73c321b0100d4f03b7f355894cf"
            "332f830e710b97ce98c8a84abd0b9481"
            "14ad176e008d33bd60f982b1ff37c855"
            "9797a06ef4f0ef61c186324e2b350638"
            "3606907b6a7c02b0f9f6157b53c867e4"
            "b9166c767b804d46a59b5216cde7a4e9"
            "9040c5a40433225ee282a1b0a06c523e"
            "af4534d7f83fa1155b0047718cbc546a"
            "0d072b04b3564eea1b422273f548271a"
            "0bb2316053fa76991955ebd63159434e"
            "cebb4e466dae5a1073a6727627097a10"
            "49e617d91d361094fa68f0ff77987130"
            "305beaba2eda04df997b714d6c6f2c29"
            "a6ad5cb4022b02709b");
    }

    std::vector<std::uint8_t> rfc8439_decryption_plaintext() {
        return hex_to_bytes(
            "496e7465726e65742d44726166747320"
            "61726520647261667420646f63756d65"
            "6e74732076616c696420666f72206120"
            "6d6178696d756d206f6620736978206d"
            "6f6e74687320616e64206d6179206265"
            "20757064617465642c207265706c6163"
            "65642c206f72206f62736f6c65746564"
            "206279206f7468657220646f63756d65"
            "6e747320617420616e792074696d652e"
            "20497420697320696e617070726f7072"
            "6961746520746f2075736520496e7465"
            "726e65742d4472616674732061732072"
            "65666572656e6365206d617465726961"
            "6c206f7220746f206369746520746865"
            "6d206f74686572207468616e20617320"
            "2fe2809c776f726b20696e2070726f67"
            "726573732e2fe2809d");
    }

    aead_type::tag_type rfc8439_decryption_tag() {
        return array_from_hex<aead_type::tag_size>("eead9d67890cbb22392336fea1851f38");
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(chacha20poly1305_mode_test_suite)

BOOST_AUTO_TEST_CASE(poly1305_key_gen_matches_rfc8439_aead_vector) {
    const aead_type::poly1305_key_type expected_key = array_from_hex<aead_type::key_size>(
        "7bac2b252db447af09b67a55a4e95584"
        "0ae1d6731075d9eb2a9375783ed553ff");

    const aead_type::poly1305_key_type actual_key =
        aead_type::poly1305_key_gen(rfc8439_aead_key(), rfc8439_aead_nonce());

    BOOST_TEST(std::equal(actual_key.begin(), actual_key.end(), expected_key.begin()));
}

BOOST_AUTO_TEST_CASE(encrypt_matches_rfc8439_aead_vector) {
    const std::vector<std::uint8_t> plaintext = rfc8439_aead_plaintext();
    const std::vector<std::uint8_t> aad = rfc8439_aead_aad();
    const std::vector<std::uint8_t> expected_ciphertext = rfc8439_aead_ciphertext();
    const aead_type::tag_type expected_tag = rfc8439_aead_tag();

    std::vector<std::uint8_t> ciphertext(plaintext.size());
    const aead_type::tag_type tag =
        aead_type::encrypt(plaintext, aad, rfc8439_aead_key(), rfc8439_aead_nonce(), ciphertext.begin());

    BOOST_TEST(std::equal(ciphertext.begin(), ciphertext.end(), expected_ciphertext.begin()));
    BOOST_TEST(std::equal(tag.begin(), tag.end(), expected_tag.begin()));
}

BOOST_AUTO_TEST_CASE(decrypt_matches_rfc8439_aead_vector) {
    const std::vector<std::uint8_t> ciphertext = rfc8439_aead_ciphertext();
    const std::vector<std::uint8_t> aad = rfc8439_aead_aad();
    const std::vector<std::uint8_t> expected_plaintext = rfc8439_aead_plaintext();

    std::vector<std::uint8_t> plaintext(ciphertext.size());
    const bool authenticated = aead_type::decrypt(ciphertext, aad, rfc8439_aead_tag(), rfc8439_aead_key(),
                                                  rfc8439_aead_nonce(), plaintext.begin());

    BOOST_TEST(authenticated);
    BOOST_TEST(std::equal(plaintext.begin(), plaintext.end(), expected_plaintext.begin()));
}

BOOST_AUTO_TEST_CASE(decrypt_matches_rfc8439_appendix_a5_vector) {
    const std::vector<std::uint8_t> ciphertext = rfc8439_decryption_ciphertext();
    const std::vector<std::uint8_t> aad = rfc8439_decryption_aad();
    const std::vector<std::uint8_t> expected_plaintext = rfc8439_decryption_plaintext();

    std::vector<std::uint8_t> plaintext(ciphertext.size());
    const bool authenticated = aead_type::decrypt(ciphertext, aad, rfc8439_decryption_tag(), rfc8439_decryption_key(),
                                                  rfc8439_decryption_nonce(), plaintext.begin());

    BOOST_TEST(authenticated);
    BOOST_TEST(std::equal(plaintext.begin(), plaintext.end(), expected_plaintext.begin()));
}

BOOST_AUTO_TEST_CASE(decrypt_rejects_modified_tag_without_writing_plaintext) {
    const std::vector<std::uint8_t> ciphertext = rfc8439_aead_ciphertext();
    const std::vector<std::uint8_t> aad = rfc8439_aead_aad();
    aead_type::tag_type tag = rfc8439_aead_tag();
    tag[0] ^= 1;

    std::vector<std::uint8_t> plaintext(ciphertext.size(), 0xa5);
    const std::vector<std::uint8_t> original_output = plaintext;
    const bool authenticated =
        aead_type::decrypt(ciphertext, aad, tag, rfc8439_aead_key(), rfc8439_aead_nonce(), plaintext.begin());

    BOOST_TEST(!authenticated);
    BOOST_TEST(std::equal(plaintext.begin(), plaintext.end(), original_output.begin()));
}

BOOST_AUTO_TEST_CASE(encrypt_supports_in_place_operation) {
    std::vector<std::uint8_t> buffer = rfc8439_aead_plaintext();
    const std::vector<std::uint8_t> aad = rfc8439_aead_aad();
    const std::vector<std::uint8_t> expected_ciphertext = rfc8439_aead_ciphertext();
    const std::vector<std::uint8_t> expected_plaintext = rfc8439_aead_plaintext();
    const aead_type::tag_type expected_tag = rfc8439_aead_tag();

    const aead_type::tag_type tag =
        aead_type::encrypt(buffer, aad, rfc8439_aead_key(), rfc8439_aead_nonce(), buffer.begin());

    BOOST_TEST(std::equal(buffer.begin(), buffer.end(), expected_ciphertext.begin()));
    BOOST_TEST(std::equal(tag.begin(), tag.end(), expected_tag.begin()));

    const bool authenticated =
        aead_type::decrypt(buffer, aad, tag, rfc8439_aead_key(), rfc8439_aead_nonce(), buffer.begin());

    BOOST_TEST(authenticated);
    BOOST_TEST(std::equal(buffer.begin(), buffer.end(), expected_plaintext.begin()));
}

BOOST_AUTO_TEST_SUITE_END()
