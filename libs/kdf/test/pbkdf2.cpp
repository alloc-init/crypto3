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

#define BOOST_TEST_MODULE pbkdf2_test

#include <nil/crypto3/hash/sha1.hpp>
#include <nil/crypto3/hash/sha2.hpp>
#include <nil/crypto3/mac/hmac.hpp>
#include <nil/crypto3/kdf/algorithm/derive.hpp>
#include <nil/crypto3/kdf/pbkdf2.hpp>

#include <boost/test/unit_test.hpp>

#include <cctype>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nil::crypto3;

namespace {
    typedef kdf::pbkdf2<mac::hmac<hashes::sha1>> pbkdf2_hmac_sha1;
    typedef kdf::pbkdf2<mac::hmac<hashes::sha2<256>>> pbkdf2_hmac_sha256;

    struct rfc6070_test_vector {
        std::vector<std::uint8_t> password;
        std::vector<std::uint8_t> salt;
        std::size_t iterations;
        std::size_t derived_key_octets;
        const char *expected;
    };

    std::vector<std::uint8_t> bytes(const char *text) {
        return std::vector<std::uint8_t>(text, text + std::char_traits<char>::length(text));
    }

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

    const rfc6070_test_vector rfc6070_vectors[] = {
        {bytes("password"), bytes("salt"), 1, 20, "0c60c80f961f0e71f3a9b524af6012062fe037a6"},
        {bytes("password"), bytes("salt"), 2, 20, "ea6c014dc72d6f8ccd1ed92ace1d41f0d8de8957"},
        {bytes("password"), bytes("salt"), 4096, 20, "4b007901b765489abead49d926f721d065a429c1"},
        {bytes("passwordPASSWORDpassword"), bytes("saltSALTsaltSALTsaltSALTsaltSALTsalt"), 4096, 25,
         "3d2eec4fe41c849b80c8d83662c0e44a8b291a964cf2f07038"},
        {{'p', 'a', 's', 's', '\0', 'w', 'o', 'r', 'd'},
         {'s', 'a', '\0', 'l', 't'},
         4096,
         16,
         "56fa6aa75548099dcc37d7f03425e0c3"},
    };
}    // namespace

BOOST_AUTO_TEST_SUITE(pbkdf2_test_suite)

BOOST_AUTO_TEST_CASE(pbkdf2_hmac_sha1_matches_fast_rfc6070_vectors) {
    for (const rfc6070_test_vector &test_case : rfc6070_vectors) {
        const std::vector<std::uint8_t> expected = hex_to_bytes(test_case.expected);
        const std::vector<std::uint8_t> actual = derive<pbkdf2_hmac_sha1>(
            test_case.password, test_case.salt, test_case.iterations, test_case.derived_key_octets);

        BOOST_TEST(actual == expected);
    }
}

BOOST_AUTO_TEST_CASE(pbkdf2_hmac_sha1_output_iterator_api_matches_vector_api) {
    const rfc6070_test_vector &test_case = rfc6070_vectors[0];

    std::vector<std::uint8_t> actual;
    derive<pbkdf2_hmac_sha1>(test_case.password, test_case.salt, test_case.iterations, test_case.derived_key_octets,
                             std::back_inserter(actual));

    BOOST_TEST(actual == derive<pbkdf2_hmac_sha1>(test_case.password, test_case.salt, test_case.iterations,
                                                  test_case.derived_key_octets));
}

BOOST_AUTO_TEST_CASE(pbkdf2_hmac_sha256_matches_known_vector) {
    const std::vector<std::uint8_t> expected =
        hex_to_bytes("120fb6cffcf8b32c43e7225256c4f837a86548c92ccc35480805987cb70be17b");

    BOOST_TEST(derive<pbkdf2_hmac_sha256>(bytes("password"), bytes("salt"), 1, 32) == expected);
}

BOOST_AUTO_TEST_CASE(pbkdf2_hmac_sha1_rejects_zero_iterations) {
    BOOST_CHECK_THROW(derive<pbkdf2_hmac_sha1>(bytes("password"), bytes("salt"), 0, 20), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(pbkdf2_hmac_sha1_rejects_zero_output_length) {
    BOOST_CHECK_THROW(derive<pbkdf2_hmac_sha1>(bytes("password"), bytes("salt"), 1, 0), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
