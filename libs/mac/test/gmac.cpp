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

#define BOOST_TEST_MODULE gmac_test

#include <nil/crypto3/block/aes.hpp>
#include <nil/crypto3/mac/algorithm/compute.hpp>
#include <nil/crypto3/mac/algorithm/verify.hpp>
#include <nil/crypto3/mac/detail/ghash.hpp>
#include <nil/crypto3/mac/gmac.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nil::crypto3;

namespace {
    struct official_gmac_vector {
        std::size_t key_bits;
        std::size_t tag_bits;
        const char *key;
        const char *iv;
        const char *associated_data;
        const char *tag;
    };

    // Selected from NIST CAVP GCM encrypt response files. These GMAC-shaped cases have PTlen = 0 and cover:
    // AES-128/192/256, every supported tag length, 96-bit and non-96-bit IVs, empty AAD, exact-block AAD,
    // partial-block AAD, and multi-block AAD.
    const official_gmac_vector selected_official_gmac_vectors[] = {
        {128, 128, "77be63708971c4e240d1cb79e8d77feb", "e0e00f19fed7ba0136a797f3", "7a43ec1d9c0a5a78a0b16533a6213cab",
         "209fcc8d3675ed938e9c7166709dd946"},
        {128, 120, "dcee093edc5724d30c8f46dc3f8dccfd", "69", "0b1df0fd5c5ea793d8d5ad9427c316704e77f9bd",
         "29be0d2f94ab5aec58f3d6e47784ee"},
        {128, 112, "8d06bed41844c9522d252098244e4381",
         "b1a1dc553e571d68537ff95badcfe86d8be01a7077e089b46c824156a80dcaa4683a0d31ed28589214611ee81e6ad544e"
         "bc874f8a72099948d3f34fa453f8ea52498a9bb469d41c2e081c011e3da9ed193946ee1440cdf27250ce2431ad1444"
         "e322980b04b53e3195821e4865ec74664dcafba910a4114b5d470bf99a4b349cf",
         "adf79c70baa2113997204234b6bbb246c974df8ac004fcd551f0d9c7bf45782bd5679905b1720800ba01bf0877fd"
         "719b",
         "3c3f648ad9557e5a0db0ba41feda"},
        {128, 104, "7674e09a02011df80cb66abf646c4448", "c9191abd553ca614e09b3784",
         "3911e9620a863cb435a74305be2d104d569a6e4ef796de23e8ec3b79b013244713dd4a2b1d3310de6ca74ecf51c44"
         "a2c3de573c8dfb6cf8c3c2763aa0add823b4426374591869ab5a5c08252d81dad73384d8965eb4a810be4e1",
         "92cf1e6b89d22b456bd8e1d63d"},
        {128, 96, "3275cab187425ec9f607215ba3e2f780", "81", "", "86cc4cc6111d10254244e84d"},
        {128, 64, "5f056b1881bc8c51567ebbc0af7f53b4",
         "5f09d04f7ca13b9e5290e005c225c1880ac2aa95c050a4bb91e864d163eb6d39ec18cb5f073f226efa221ffb1c490"
         "4f8346d1cbed485e669a4ac8c7ff5775b5e9e04dfa5b81e2c10df3eaef06957d848d2f9f289063284c2e6b3"
         "dc0a9e801b64a8c95c5cf1c17939418c38f95f64863008464cec06a22e69b42b99a639c3e30d",
         "169eb077a0ea36dd5c5d7174f3d86a62", "2b5ace0369d477e0"},
        {128, 32, "1ea8bc7bc395877342e911d476ae32e2", "3333c4f484aa2f906a8e7c28",
         "30a5a9270a017047eab5bf10effc5055532b874b", "bd514ad3"},
        {192, 128, "8eec76dd2be43bbb5bc7e02ca274126185d1d10063c7e536", "68", "4aca89cc2265c64c12f6b9c43003054b",
         "61b500419ccd8b13bf1f6f15d8c49446"},
        {192, 120, "d2f2b7f5af5fd06fe8e37b7a07abc764fa061095b2fbc1a0",
         "23ac2450a6a9e7c8168bf9e182f55fa3840ac4787fddaab9103697abdaa54878273664d6c41bbbab4ffe9ba9ca8"
         "eabf05fd63edaff770fb2ee3754a09b624c124eb79c17daf3633587e234994b267842d40450169be7673a191"
         "40cbfd1e41f0bd72a84b4e636e0ea0612c96eb3b1686093ff3aafaa1103c30889add2987f3941",
         "b8ada3d292d8098fee5ccadfa6359c79d3f7c4a0", "5ab15c1431595e6bc21f38089b3cf9"},
        {192, 112, "2399adc61851d6603589191d7826dab8f73e61665f40bab3", "af58ba14f7a0fd41f30b82c1",
         "d6256004d1a53446171f28f05ae533d04c1c5501d3332e9c25100f7c4d30220474c666973cc8965aee2d8edde20"
         "acf2a",
         "ea3e141f08f761576b5fd57f3158"},
        {192, 104, "38c3c03f588664036e9b92ec830807d87ea50ea771d4d4d3", "17",
         "5aada6d4a95836315ddfdeee497f2958652d704eb1f4ca4501db10d68951e5754707a8e41dc31eb3582c8fe13fc"
         "11107b3f48d0424408412d8e57f367b9090c9f43e5c5d216129be246a7e731f36bb680edb25a409e6dfbfb603",
         "b2f65f0a958007cc0f7c46e60f"},
        {192, 96, "bf1107662c9d3d2032360301aaad2949e253850e6d85bdbb",
         "a23cb98c19be0d4768b16cd77fa4bdbffd9ebac5ccd6c920c25c8c061c45374e73110b8f6c6147c8b446d6fd3aa5"
         "d03d946017c466a81c677292118b0b05466b9369c6f434a5db5de1cb1265ab77f0ef8bcf3b36ac695b915354"
         "85b0ee8df815ffb0bc09d317eae863a0e9debc52f34ecd7cbf38ef2f7fae5df8f0ec4de6d063",
         "", "b0f3ca0b16b1603a8873714c"},
        {192, 64, "b8a3d93121300c7009aa1dd6ffaaa22c170f077a51bf0f59", "c8f6dd62b3b1d786e07f7c18",
         "e067341f5ae51d6f710599c77dcfa5c2", "9e23d2a8b2dff3bf"},
        {192, 32, "333a855188acea05cfb1295fd3f696d58fb1e081fa9784ca", "e7", "246fc16d9756d7a8187e4fa4638383c2944920ec",
         "6995feef"},
        {256, 128, "f43e7ecc3634beccb40763b6fd5344957597ea1aa831753a0c3a56fe6f4b7c5a",
         "2d80ad96187ed28295fcaf6f780765a1df7aea9d1203844c03416c9d4867fb06ad5461ffbccf691141d5f37d408"
         "c54c4f973393c77c1edfa1004acaaeb6cdddd97b00b8e04e58f1324090e369d3149f20df143da68c521dc552"
         "ef4edba016133d040bdc331b703225de6927a3735a98750b5d4c7e968d16b55843f9f4f341c42",
         "2acd2a55be81b12763aa9803b5a835b6", "b619670f80a37317eb99d7c2f41de176"},
        {256, 120, "dc5d9151df042881721c840e7b1bc1890feb1ece43df04269e373696f82a5dbc", "9be5b02b8282d9b92286a365",
         "027a4b31fa71c0a5b5ef5b85775275236ccfaaf8", "68c9d4a19b362dd20cabc18c24b128"},
        {256, 112, "95421d106dd56ee646f4c753fa33e78339d1be9115eb0dd506d2cd9c2d3e56dc", "5c",
         "d10e6280de199c4c378f3c38a6af761b30c125cae3645fc6e9c65cddc79a37251df515ea376ccb47f472ca5d7"
         "ed4d274",
         "b5789ccc06e852f3bae44916ddd9"},
        {256, 104, "07e38981b8f8f4950680914c7806509e98e02fb95afdfe8ab4dfba3d4a83d99d",
         "52a727bb1f0de0ef3ffd3b7cea05a9edc42dffdc70482d5632f4f046e599c9afe56bd6e197be05229feda8bcb"
         "549f506867fb47caf26127cdbb0eefd8af1e9c5be9dbc4fce813069b46664f16191df45f4eeb6f674682d0"
         "e54529663b0b3026774b1d2d2b7adbb5309a2c120cc071b38196a860384a3f6ef7e936004bcda63cd",
         "2bfec0a84f70b3752f78283f25a233a793638d5dc47dd3e9bade5cb3b7e8fb7a53ae651476b5a3885e0775804"
         "a80d33c17708c749246eb4513c9dc6a40127526072b48f615e55902b154c78e4363a4dd4894f2a57c03e05565e8",
         "452bb36914488448d7bbe88916"},
        {256, 96, "98ebf7a58db8b8371d9069171190063cc1fdc1927e49a3385f890d41a838619c", "3e6db953bd4e641de644e50a", "",
         "2fb9c3e41fff24ef07437c47"},
        {256, 64, "bbcbe10b4c5b2062913badb8dcc39fd174f1a2792b30947a76b816be879d4990", "67",
         "78fcaff51369e8a284e8970932ec07e0", "adcec7915ed4722b"},
        {256, 32, "c5d17247a46c8733858caae85580087d1dd61579eb3430f31e818514f76e54b3",
         "5f7eb902c3b2a8cc2e908dde2558c94d2f11d7597fda31b98ebfeb876480a81e740ae97be3187dc9f388b35a95"
         "fcd6edc0660ebe258df9732bdde8ca44e2d9d5600b4c1923cc9190d1fceae91e6b4564ea5d7c791b703ace"
         "15edf6c065a1400ad3a6093b745c185db719ad635ce48737350d319a14edcda253d85b71041084ae",
         "02e49c0bfdedbbeb790000b15ebb54d140ef7c02", "1b372f04"},
    };

    // Selected from NIST CAVP GCM decrypt response files. These are GMAC-shaped FAIL cases.
    const official_gmac_vector rejected_official_gmac_vectors[] = {
        {128, 128, "d1f6af919cde85661208bdce0c27cb22", "898c6929b435017bf031c3c5", "7c5faa40e636bbc91107e68010c92b9f",
         "ae45f11777540a2caeb128be8092468a"},
        {192, 96, "9bf81228dd4c5a6242b1c45692461a4de3d8ce585abd6cfb", "da", "4546f07aee0d078dfecf1f31847bdbb2ce05a511",
         "5620f3e98d80cbde78c0cc90"},
        {256, 32, "685253ee063933134994a21d37ff4e02b038773c05001300291b6576a26b4476",
         "53b005fdd71666aa71670b553eb7fb72fff5ed4163d3521c39680e65b6198eb19430120aa2b58a9a92b618f30e"
         "119b0585d5dd6d88031dde0c2e40514b368c5574a43c753d27c71ff1112e9e1e9b9735c4cdb48f0cd737c3"
         "feba3e8e6a0d08aad9530544f4573d4bbe98dd4e5b3e2a0eb7b919e7dcb4b2af8cc6a0e7afeeb911",
         "45d730f4c640d0e62fbae2ffa56fbaf8c0c5b86fc42e9b0d4a3ab20519a3891ce6c4cbe2d2647cc7e1cfa068"
         "c357518ed5d57059347bebd17f4ed6307f57d82ae8aa62b60ef6870905a172091fd3560be0c3c458cc13416e9ffa",
         "01787afe"},
    };

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

    mac::detail::ghash::block_type block_from_hex(const std::string &hex) {
        const std::vector<std::uint8_t> bytes = hex_to_bytes(hex);
        if (bytes.size() != mac::detail::ghash::block_octets) {
            throw std::invalid_argument("GHASH block test vector must be 16 bytes");
        }

        mac::detail::ghash::block_type block = {0};
        std::copy(bytes.begin(), bytes.end(), block.begin());
        return block;
    }

    void check_equal(const mac::detail::ghash::block_type &actual, const mac::detail::ghash::block_type &expected) {
        BOOST_TEST(std::equal(actual.begin(), actual.end(), expected.begin()));
    }

    template<std::size_t KeyBits, std::size_t TagBits>
    std::vector<std::uint8_t> compute_public_gmac_tag_for(const official_gmac_vector &test_case) {
        typedef mac::gmac<block::aes<KeyBits>, TagBits> mac_type;

        const std::vector<std::uint8_t> key_bytes = hex_to_bytes(test_case.key);
        const std::vector<std::uint8_t> iv = hex_to_bytes(test_case.iv);
        const std::vector<std::uint8_t> associated_data = hex_to_bytes(test_case.associated_data);

        mac::mac_key<mac_type> key(key_bytes, iv);
        const typename mac_type::digest_type tag = compute<mac_type>(associated_data, key);
        return std::vector<std::uint8_t>(tag.begin(), tag.end());
    }

    template<std::size_t KeyBits, typename Operation>
    auto dispatch_gmac_tag_size(const official_gmac_vector &test_case, Operation operation)
        -> decltype(operation.template operator()<KeyBits, 128>(test_case)) {
        switch (test_case.tag_bits) {
            case 128:
                return operation.template operator()<KeyBits, 128>(test_case);
            case 120:
                return operation.template operator()<KeyBits, 120>(test_case);
            case 112:
                return operation.template operator()<KeyBits, 112>(test_case);
            case 104:
                return operation.template operator()<KeyBits, 104>(test_case);
            case 96:
                return operation.template operator()<KeyBits, 96>(test_case);
            case 64:
                return operation.template operator()<KeyBits, 64>(test_case);
            case 32:
                return operation.template operator()<KeyBits, 32>(test_case);
            default:
                throw std::invalid_argument("unsupported GMAC tag size");
        }
    }

    template<typename Operation>
    auto dispatch_gmac_key_size(const official_gmac_vector &test_case, Operation operation)
        -> decltype(operation.template operator()<128, 128>(test_case)) {
        switch (test_case.key_bits) {
            case 128:
                return dispatch_gmac_tag_size<128>(test_case, operation);
            case 192:
                return dispatch_gmac_tag_size<192>(test_case, operation);
            case 256:
                return dispatch_gmac_tag_size<256>(test_case, operation);
            default:
                throw std::invalid_argument("unsupported GMAC key size");
        }
    }

    struct compute_gmac_operation {
        template<std::size_t KeyBits, std::size_t TagBits>
        std::vector<std::uint8_t> operator()(const official_gmac_vector &test_case) const {
            return compute_public_gmac_tag_for<KeyBits, TagBits>(test_case);
        }
    };

    template<std::size_t KeyBits, std::size_t TagBits>
    bool verify_public_gmac_tag_for(const official_gmac_vector &test_case) {
        typedef mac::gmac<block::aes<KeyBits>, TagBits> mac_type;

        const std::vector<std::uint8_t> key_bytes = hex_to_bytes(test_case.key);
        const std::vector<std::uint8_t> iv = hex_to_bytes(test_case.iv);
        const std::vector<std::uint8_t> associated_data = hex_to_bytes(test_case.associated_data);
        const std::vector<std::uint8_t> tag = hex_to_bytes(test_case.tag);

        mac::mac_key<mac_type> key(key_bytes, iv);
        return verify<mac_type>(associated_data, key, tag);
    }

    struct verify_gmac_operation {
        template<std::size_t KeyBits, std::size_t TagBits>
        bool operator()(const official_gmac_vector &test_case) const {
            return verify_public_gmac_tag_for<KeyBits, TagBits>(test_case);
        }
    };

    official_gmac_vector aes128_gmac_vector() {
        return selected_official_gmac_vectors[0];
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(gmac_test_suite)

BOOST_AUTO_TEST_CASE(gf128_multiply_treats_zero_as_absorbing) {
    const mac::detail::ghash::block_type zero = mac::detail::ghash::zero_block();
    const mac::detail::ghash::block_type h = block_from_hex("66e94bd4ef8a2c3b884cfa59ca342b2e");

    check_equal(mac::detail::ghash::multiply(zero, h), zero);
    check_equal(mac::detail::ghash::multiply(h, zero), zero);
}

BOOST_AUTO_TEST_CASE(ghash_matches_nist_gcm_test_case_2) {
    const mac::detail::ghash::block_type h = block_from_hex("66e94bd4ef8a2c3b884cfa59ca342b2e");
    const std::vector<std::uint8_t> ciphertext = hex_to_bytes("0388dace60b6a392f328c2b971b2fe78");
    const mac::detail::ghash::block_type expected = block_from_hex("f38cbb1ad69223dcc3457ae5b6b0f885");

    mac::detail::ghash::state ghash(h);
    ghash.update(ciphertext.begin(), ciphertext.end());
    ghash.pad();
    ghash.update_lengths(0, 128);

    check_equal(ghash.digest(), expected);
}

BOOST_AUTO_TEST_CASE(ghash_streaming_updates_match_single_update) {
    const mac::detail::ghash::block_type h = block_from_hex("66e94bd4ef8a2c3b884cfa59ca342b2e");
    const std::vector<std::uint8_t> ciphertext = hex_to_bytes(
        "0388dace60b6a392f328c2b971b2fe78"
        "42831ec2217774244b7221b784d0d49c"
        "e3aa212f2c02a4e035c17e2329aca12e");

    mac::detail::ghash::state single_update(h);
    single_update.update(ciphertext.begin(), ciphertext.end());
    single_update.pad();
    single_update.update_lengths(0, static_cast<std::uint64_t>(ciphertext.size() * 8));

    mac::detail::ghash::state streamed(h);
    streamed.update(ciphertext.begin(), ciphertext.begin() + 5);
    streamed.update(ciphertext.begin() + 5, ciphertext.begin() + 37);
    streamed.update(ciphertext.begin() + 37, ciphertext.end());
    streamed.pad();
    streamed.update_lengths(0, static_cast<std::uint64_t>(ciphertext.size() * 8));

    check_equal(streamed.digest(), single_update.digest());
}

BOOST_AUTO_TEST_CASE(public_gmac_matches_selected_nist_cavp_vectors) {
    for (const official_gmac_vector &test_case : selected_official_gmac_vectors) {
        const std::vector<std::uint8_t> expected = hex_to_bytes(test_case.tag);
        const std::vector<std::uint8_t> actual = dispatch_gmac_key_size(test_case, compute_gmac_operation());

        BOOST_TEST_CONTEXT("key_bits " << test_case.key_bits << " tag_bits " << test_case.tag_bits) {
            BOOST_TEST(actual == expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(public_gmac_verifies_selected_nist_cavp_vectors) {
    for (const official_gmac_vector &test_case : selected_official_gmac_vectors) {
        BOOST_TEST_CONTEXT("key_bits " << test_case.key_bits << " tag_bits " << test_case.tag_bits) {
            BOOST_TEST(dispatch_gmac_key_size(test_case, verify_gmac_operation()));
        }
    }
}

BOOST_AUTO_TEST_CASE(public_gmac_rejects_selected_nist_cavp_fail_vectors) {
    for (const official_gmac_vector &test_case : rejected_official_gmac_vectors) {
        BOOST_TEST_CONTEXT("key_bits " << test_case.key_bits << " tag_bits " << test_case.tag_bits) {
            BOOST_TEST(!dispatch_gmac_key_size(test_case, verify_gmac_operation()));
        }
    }
}

BOOST_AUTO_TEST_CASE(public_gmac_rejects_wrong_key_size) {
    typedef mac::gmac<block::aes<128>> mac_type;

    const official_gmac_vector test_case = aes128_gmac_vector();
    const std::vector<std::uint8_t> short_key(15, 0);
    const std::vector<std::uint8_t> iv = hex_to_bytes(test_case.iv);

    BOOST_CHECK_THROW(mac::mac_key<mac_type> key(short_key, iv), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(public_gmac_rejects_empty_iv) {
    typedef mac::gmac<block::aes<128>> mac_type;

    const official_gmac_vector test_case = aes128_gmac_vector();
    const std::vector<std::uint8_t> key_bytes = hex_to_bytes(test_case.key);
    const std::vector<std::uint8_t> empty_iv;

    BOOST_CHECK_THROW(mac::mac_key<mac_type> key(key_bytes, empty_iv), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(public_gmac_verify_rejects_wrong_tag_length) {
    typedef mac::gmac<block::aes<128>> mac_type;

    const official_gmac_vector test_case = aes128_gmac_vector();
    const std::vector<std::uint8_t> key_bytes = hex_to_bytes(test_case.key);
    const std::vector<std::uint8_t> iv = hex_to_bytes(test_case.iv);
    const std::vector<std::uint8_t> associated_data = hex_to_bytes(test_case.associated_data);
    std::vector<std::uint8_t> short_tag = hex_to_bytes(test_case.tag);
    short_tag.pop_back();

    mac::mac_key<mac_type> key(key_bytes, iv);
    BOOST_TEST(!verify<mac_type>(associated_data, key, short_tag));
}

BOOST_AUTO_TEST_CASE(public_gmac_verify_rejects_modified_tag) {
    typedef mac::gmac<block::aes<128>> mac_type;

    const official_gmac_vector test_case = aes128_gmac_vector();
    const std::vector<std::uint8_t> key_bytes = hex_to_bytes(test_case.key);
    const std::vector<std::uint8_t> iv = hex_to_bytes(test_case.iv);
    const std::vector<std::uint8_t> associated_data = hex_to_bytes(test_case.associated_data);
    std::vector<std::uint8_t> tag = hex_to_bytes(test_case.tag);
    tag[0] ^= 1;

    mac::mac_key<mac_type> key(key_bytes, iv);
    BOOST_TEST(!verify<mac_type>(associated_data, key, tag));
}

BOOST_AUTO_TEST_SUITE_END()
