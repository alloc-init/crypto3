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

#define BOOST_TEST_MODULE cmac_test

#include <nil/crypto3/block/aes.hpp>
#include <nil/crypto3/mac/algorithm/compute.hpp>
#include <nil/crypto3/mac/cmac.hpp>

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

    std::array<std::uint8_t, 16> tag_from_hex(const std::string &hex) {
        const std::vector<std::uint8_t> bytes = hex_to_bytes(hex);
        if (bytes.size() != 16) {
            throw std::invalid_argument("CMAC tag test vector must be 16 bytes");
        }

        std::array<std::uint8_t, 16> tag = {0};
        std::copy(bytes.begin(), bytes.end(), tag.begin());
        return tag;
    }

    template<typename MacType>
    void check_cmac_vector(const std::vector<std::uint8_t> &key_bytes,
                           const std::vector<std::uint8_t> &message,
                           const std::array<std::uint8_t, 16> &expected) {
        mac::mac_key<MacType> key(key_bytes);

        const typename MacType::digest_type tag = compute<MacType>(message, key);
        BOOST_TEST(std::equal(tag.begin(), tag.end(), expected.begin()));

        mac::computation_accumulator_set<mac::computation_policy<MacType>> acc(key);
        for (std::size_t i = 0; i < message.size(); ++i) {
            compute<MacType>(message.begin() + i, message.begin() + i + 1, acc);
        }

        const typename MacType::digest_type streamed_tag =
            accumulators::extract::mac<mac::computation_policy<MacType>>(acc);
        BOOST_TEST(std::equal(streamed_tag.begin(), streamed_tag.end(), expected.begin()));
    }

    const std::string message_64 =
        "6bc1bee22e409f96e93d7e117393172a"
        "ae2d8a571e03ac9c9eb76fac45af8e51"
        "30c81c46a35ce411e5fbc1191a0a52ef"
        "f69f2445df4f9b17ad2b417be66c3710";
}    // namespace

BOOST_AUTO_TEST_SUITE(cmac_test_suite)

BOOST_AUTO_TEST_CASE(aes128_cmac_matches_nist_sp800_38b_examples) {
    typedef mac::cmac<block::aes<128>> cmac_type;

    const std::vector<std::uint8_t> key = hex_to_bytes("2b7e151628aed2a6abf7158809cf4f3c");

    check_cmac_vector<cmac_type>(key, hex_to_bytes(""), tag_from_hex("bb1d6929e95937287fa37d129b756746"));
    check_cmac_vector<cmac_type>(key, hex_to_bytes(message_64.substr(0, 32)),
                                 tag_from_hex("070a16b46b4d4144f79bdd9dd04a287c"));
    check_cmac_vector<cmac_type>(key, hex_to_bytes(message_64.substr(0, 80)),
                                 tag_from_hex("dfa66747de9ae63030ca32611497c827"));
    check_cmac_vector<cmac_type>(key, hex_to_bytes(message_64), tag_from_hex("51f0bebf7e3b9d92fc49741779363cfe"));
}

BOOST_AUTO_TEST_CASE(aes192_cmac_matches_nist_sp800_38b_examples) {
    typedef mac::cmac<block::aes<192>> cmac_type;

    const std::vector<std::uint8_t> key = hex_to_bytes("8e73b0f7da0e6452c810f32b809079e562f8ead2522c6b7b");

    check_cmac_vector<cmac_type>(key, hex_to_bytes(""), tag_from_hex("d17ddf46adaacde531cac483de7a9367"));
    check_cmac_vector<cmac_type>(key, hex_to_bytes(message_64.substr(0, 32)),
                                 tag_from_hex("9e99a7bf31e710900662f65e617c5184"));
    check_cmac_vector<cmac_type>(key, hex_to_bytes(message_64.substr(0, 80)),
                                 tag_from_hex("8a1de5be2eb31aad089a82e6ee908b0e"));
    check_cmac_vector<cmac_type>(key, hex_to_bytes(message_64), tag_from_hex("a1d5df0eed790f794d77589659f39a11"));
}

BOOST_AUTO_TEST_CASE(aes256_cmac_matches_nist_sp800_38b_examples) {
    typedef mac::cmac<block::aes<256>> cmac_type;

    const std::vector<std::uint8_t> key =
        hex_to_bytes("603deb1015ca71be2b73aef0857d77811f352c073b6108d72d9810a30914dff4");

    check_cmac_vector<cmac_type>(key, hex_to_bytes(""), tag_from_hex("028962f61b7bf89efc6b551f4667d983"));
    check_cmac_vector<cmac_type>(key, hex_to_bytes(message_64.substr(0, 32)),
                                 tag_from_hex("28a7023f452e8f82bd4bf28d8c37c35c"));
    check_cmac_vector<cmac_type>(key, hex_to_bytes(message_64.substr(0, 80)),
                                 tag_from_hex("aaf3d8f1de5640c232f5b169b9c911e6"));
    check_cmac_vector<cmac_type>(key, hex_to_bytes(message_64), tag_from_hex("e1992190549f6ed5696a2c056c315410"));
}

BOOST_AUTO_TEST_CASE(cmac_rejects_wrong_key_size) {
    std::vector<std::uint8_t> short_key(15, 0);
    std::vector<std::uint8_t> long_key(17, 0);

    BOOST_CHECK_THROW(mac::mac_key<mac::cmac<block::aes<128>>> key(short_key), std::invalid_argument);
    BOOST_CHECK_THROW(mac::mac_key<mac::cmac<block::aes<128>>> key(long_key), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
