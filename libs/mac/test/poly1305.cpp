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

#define BOOST_TEST_MODULE poly1305_test

#include <nil/crypto3/mac/algorithm/compute.hpp>
#include <nil/crypto3/mac/poly1305.hpp>

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
            throw std::invalid_argument("Poly1305 tag test vector must be 16 bytes");
        }

        std::array<std::uint8_t, 16> tag = {0};
        std::copy(bytes.begin(), bytes.end(), tag.begin());
        return tag;
    }

    std::vector<std::uint8_t> bytes_from_text(const std::string &text) {
        return std::vector<std::uint8_t>(text.begin(), text.end());
    }

    std::vector<std::uint8_t> key_from_r_s(const std::string &r, const std::string &s) {
        std::vector<std::uint8_t> key = hex_to_bytes(r + s);
        if (key.size() != 32) {
            throw std::invalid_argument("Poly1305 key test vector must be 32 bytes");
        }
        return key;
    }

    template<typename MacType>
    void check_poly1305_vector_for(const std::vector<std::uint8_t> &key_bytes,
                                   const std::vector<std::uint8_t> &message,
                                   const std::array<std::uint8_t, 16> &expected) {
        mac::mac_key<MacType> key(key_bytes);

        const typename MacType::digest_type tag = compute<MacType>(message, key);
        BOOST_TEST(std::equal(tag.begin(), tag.end(), expected.begin()));

        mac::computation_accumulator_set<mac::computation_policy<MacType>> acc(key);
        const std::size_t split = message.size() / 2;
        compute<MacType>(message.begin(), message.begin() + split, acc);
        compute<MacType>(message.begin() + split, message.end(), acc);

        const typename MacType::digest_type streamed_tag =
            accumulators::extract::mac<mac::computation_policy<MacType>>(acc);
        BOOST_TEST(std::equal(streamed_tag.begin(), streamed_tag.end(), expected.begin()));
    }

    void check_poly1305_vector(const std::vector<std::uint8_t> &key_bytes,
                               const std::vector<std::uint8_t> &message,
                               const std::array<std::uint8_t, 16> &expected) {
        check_poly1305_vector_for<mac::poly1305>(key_bytes, message, expected);
        check_poly1305_vector_for<mac::poly1305_reference>(key_bytes, message, expected);
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(poly1305_test_suite)

BOOST_AUTO_TEST_CASE(poly1305_matches_rfc8439_main_test_vector) {
    const std::vector<std::uint8_t> key =
        hex_to_bytes("85d6be7857556d337f4452fe42d506a80103808afb0db2fd4abff6af4149f51b");
    const std::vector<std::uint8_t> message = bytes_from_text("Cryptographic Forum Research Group");
    const std::array<std::uint8_t, 16> expected = tag_from_hex("a8061dc1305136c6c22b8baf0c0127a9");

    check_poly1305_vector(key, message, expected);
}

BOOST_AUTO_TEST_CASE(poly1305_matches_rfc8439_appendix_a3_vectors) {
    const std::string ietf_text =
        "Any submission to the IETF intended by the Contributor for publication as all or part of an IETF "
        "Internet-Draft or RFC and any statement made within the context of an IETF activity is considered an "
        "\"IETF Contribution\". Such statements include oral statements in IETF sessions, as well as written and "
        "electronic communications made at any time or place, which are addressed to";

    const std::string jabberwocky =
        "'Twas brillig, and the slithy toves\n"
        "Did gyre and gimble in the wabe:\n"
        "All mimsy were the borogoves,\n"
        "And the mome raths outgrabe.";

    check_poly1305_vector(hex_to_bytes("0000000000000000000000000000000000000000000000000000000000000000"),
                          hex_to_bytes("00000000000000000000000000000000"
                                       "00000000000000000000000000000000"
                                       "00000000000000000000000000000000"
                                       "00000000000000000000000000000000"),
                          tag_from_hex("00000000000000000000000000000000"));

    check_poly1305_vector(key_from_r_s("00000000000000000000000000000000", "36e5f6b5c5e06070f0efca96227a863e"),
                          bytes_from_text(ietf_text), tag_from_hex("36e5f6b5c5e06070f0efca96227a863e"));

    check_poly1305_vector(key_from_r_s("36e5f6b5c5e06070f0efca96227a863e", "00000000000000000000000000000000"),
                          bytes_from_text(ietf_text), tag_from_hex("f3477e7cd95417af89a6b8794c310cf0"));

    check_poly1305_vector(key_from_r_s("1c9240a5eb55d38af333888604f6b5f0", "473917c1402b80099dca5cbc207075c0"),
                          bytes_from_text(jabberwocky), tag_from_hex("4541669a7eaaee61e708dc7cbcc5eb62"));

    check_poly1305_vector(key_from_r_s("02000000000000000000000000000000", "00000000000000000000000000000000"),
                          hex_to_bytes("ffffffffffffffffffffffffffffffff"),
                          tag_from_hex("03000000000000000000000000000000"));

    check_poly1305_vector(key_from_r_s("02000000000000000000000000000000", "ffffffffffffffffffffffffffffffff"),
                          hex_to_bytes("02000000000000000000000000000000"),
                          tag_from_hex("03000000000000000000000000000000"));

    check_poly1305_vector(key_from_r_s("01000000000000000000000000000000", "00000000000000000000000000000000"),
                          hex_to_bytes("ffffffffffffffffffffffffffffffff"
                                       "f0ffffffffffffffffffffffffffffff"
                                       "11000000000000000000000000000000"),
                          tag_from_hex("05000000000000000000000000000000"));

    check_poly1305_vector(key_from_r_s("01000000000000000000000000000000", "00000000000000000000000000000000"),
                          hex_to_bytes("ffffffffffffffffffffffffffffffff"
                                       "fbfefefefefefefefefefefefefefefe"
                                       "01010101010101010101010101010101"),
                          tag_from_hex("00000000000000000000000000000000"));

    check_poly1305_vector(key_from_r_s("02000000000000000000000000000000", "00000000000000000000000000000000"),
                          hex_to_bytes("fdffffffffffffffffffffffffffffff"),
                          tag_from_hex("faffffffffffffffffffffffffffffff"));

    check_poly1305_vector(key_from_r_s("01000000000000000400000000000000", "00000000000000000000000000000000"),
                          hex_to_bytes("e33594d7505e43b90000000000000000"
                                       "3394d7505e4379cd0100000000000000"
                                       "00000000000000000000000000000000"
                                       "01000000000000000000000000000000"),
                          tag_from_hex("14000000000000005500000000000000"));

    check_poly1305_vector(key_from_r_s("01000000000000000400000000000000", "00000000000000000000000000000000"),
                          hex_to_bytes("e33594d7505e43b90000000000000000"
                                       "3394d7505e4379cd0100000000000000"
                                       "00000000000000000000000000000000"),
                          tag_from_hex("13000000000000000000000000000000"));
}

BOOST_AUTO_TEST_CASE(poly1305_rejects_wrong_key_size) {
    std::vector<std::uint8_t> short_key(31, 0);

    BOOST_CHECK_THROW(mac::mac_key<mac::poly1305> key(short_key), std::invalid_argument);
    BOOST_CHECK_THROW(mac::mac_key<mac::poly1305_reference> key(short_key), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
