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

#define BOOST_TEST_MODULE crypto3_marshalling_modified_sap_verification_key_test

#include <boost/algorithm/hex.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <forward_list>
#include <iterator>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <nil/crypto3/marshalling/zk/types/modified_sap/proof.hpp>
#include <nil/crypto3/marshalling/zk/types/modified_sap/verification_key.hpp>
#include <nil/crypto3/random/chacha_urbg.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/policy.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/transcript.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap_snark.hpp>

#include "detail/marshalling.hpp"

namespace {
    namespace types = nil::crypto3::marshalling::types;
    namespace snark = nil::crypto3::zk::snark;
    namespace test_tools = nil::crypto3::marshalling::test_tools;
    using test_tools::encode;
    using test_tools::write_integer;
    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using scalar_field_type = curve_type::scalar_field_type;
    using scalar_type = scalar_field_type::value_type;
    using base_field_type = curve_type::base_field_type;
    using base_type = base_field_type::value_type;
    using integral_type = base_field_type::integral_type;
    using point_type = curve_type::g2_type<>::value_type;
    using gt_type = curve_type::gt_type::value_type;
    using policy_type = snark::modified_sap_policy<curve_type,
                                                   snark::modified_sap_bn254_exact_pairing_policy,
                                                   snark::modified_sap_bn254_poseidon_transcript_policy>;
    using scheme_type = snark::modified_sap_snark<policy_type>;
    using verifier_type = snark::modified_sap_verifier<policy_type>;
    using key_type = policy_type::verification_key_type;
    using endianness = nil::marshalling::option::big_endian;
    using type_base = nil::marshalling::field_type<endianness>;
    using marshalled_type = types::modified_sap_verification_key<type_base, key_type>;
    using gt_codec = types::detail::validated_field_element<type_base, gt_type>;
    using status_type = nil::marshalling::status_type;
    constexpr std::size_t encoded_size = 2144 + 2 * sizeof(std::size_t);
    constexpr std::size_t dimensions_offset = 2112;
    constexpr std::size_t digest_offset = dimensions_offset + 2 * sizeof(std::size_t);

    // A valid parameter-shape fixture; generated setup keys are tested separately below.
    key_type test_key() {
        key_type key;
        key.tau_g2 = point_type::one();
        key.gamma_inverse_g2 = point_type::one();
        key.num_variables = 3;
        key.domain_size = 4;
        key.circuit_digest = base_type(19);
        return key;
    }

    key_type decode(const std::vector<std::uint8_t> &bytes) {
        return types::make_modified_sap_verification_key<policy_type, endianness>(
            test_tools::decode<marshalled_type>(bytes));
    }

    void check_invalid_encoding(const std::vector<std::uint8_t> &bytes) {
        marshalled_type filled;
        auto input = bytes.begin();
        BOOST_CHECK(filled.read(input, bytes.size()) == status_type::invalid_msg_data);
    }

    void check_invalid_key(const key_type &key, const marshalled_type &filled) {
        BOOST_CHECK(!verifier_type::validate_verification_key(key));
        BOOST_CHECK_THROW((types::fill_modified_sap_verification_key<policy_type, endianness>(key)),
                          std::invalid_argument);
        BOOST_CHECK_THROW((types::make_modified_sap_verification_key<policy_type, endianness>(filled)),
                          std::invalid_argument);
        // The bytes can be canonical while the key parameters violate native invariants.
        BOOST_CHECK_THROW(decode(encode(filled)), std::invalid_argument);
    }
}    // namespace

BOOST_AUTO_TEST_CASE(verification_key_wire_format_and_identity_gt) {
    static_assert(marshalled_type::max_length() == encoded_size);
    const auto key = test_key();
    const auto filled = types::fill_modified_sap_verification_key<policy_type, endianness>(key);
    BOOST_CHECK_EQUAL(filled.length(), encoded_size);
    const auto bytes = encode(filled);

    std::vector<std::uint8_t> generator;
    boost::algorithm::unhex(std::string("198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2"
                                        "1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed"),
                            std::back_inserter(generator));
    std::vector<std::uint8_t> expected(encoded_size, 0);
    for (std::size_t i = 0; i < 3; ++i) {
        std::copy(generator.begin(), generator.end(), expected.begin() + 64 * i);
    }
    // Each GT identity has first Fp component one and eleven zero components. No alpha-array count is stored.
    for (std::size_t i = 0; i < 5; ++i) {
        expected[192 + 384 * i + 31] = 1;
    }
    expected[dimensions_offset + sizeof(std::size_t) - 1] = 3;
    expected[digest_offset - 1] = 4;
    expected.back() = 19;
    BOOST_CHECK(bytes == expected);
    BOOST_CHECK(decode(bytes) == key);
}

BOOST_AUTO_TEST_CASE(verification_key_preserves_gt_component_and_alpha_order) {
    const auto generator = nil::crypto3::algebra::pair_reduced<curve_type, policy_type::pairing_policy_type>(
        curve_type::g1_type<>::value_type::one(), point_type::one());
    BOOST_REQUIRE(generator.has_value());
    auto key = test_key();
    key.alpha_z_vanishing_gt = *generator;
    key.alpha_gt = {generator->pow(2), generator->pow(3), generator->pow(5), generator->pow(7)};
    const auto bytes = encode(types::fill_modified_sap_verification_key<policy_type, endianness>(key));
    for (std::size_t i = 0; i < 5; ++i) {
        const auto &value = i == 0 ? key.alpha_z_vanishing_gt : key.alpha_gt[i - 1];
        std::vector<std::uint8_t> expected(384);
        std::size_t offset = 0;
        // Spell out the tower order independently of the extension-field marshaller's flattening helper.
        for (const auto &outer : value.data) {
            for (const auto &middle : outer.data) {
                for (const auto &component : middle.data) {
                    write_integer<endianness>(expected, offset, component.to_integral());
                    offset += 32;
                }
            }
        }
        BOOST_CHECK(std::equal(expected.begin(), expected.end(), bytes.begin() + 192 + 384 * i));
    }
    BOOST_CHECK(decode(bytes) == key);
}

BOOST_AUTO_TEST_CASE(verification_key_accepts_boundary_metadata_and_fp_digest) {
    auto key = test_key();
    key.num_variables = 1;
    constexpr auto max_domain_log = nil::crypto3::algebra::fields::arithmetic_params<scalar_field_type>::s;
    for (std::size_t domain_size : {std::size_t(2), std::size_t(1) << max_domain_log}) {
        key.domain_size = domain_size;
        for (const auto &digest : {base_type::zero(), base_type(scalar_field_type::modulus), -base_type::one()}) {
            key.circuit_digest = digest;    // Fr's modulus is a valid Fp digest; the upper boundary is p - 1.
            BOOST_CHECK(decode(encode(types::fill_modified_sap_verification_key<policy_type, endianness>(key))) == key);
        }
    }
}

BOOST_AUTO_TEST_CASE(verification_key_read_leaves_following_bytes) {
    const auto key = test_key();
    auto bytes = encode(types::fill_modified_sap_verification_key<policy_type, endianness>(key));
    bytes.push_back(0xa5);
    marshalled_type filled;
    const auto *input = bytes.data();
    BOOST_REQUIRE(filled.read(input, bytes.size()) == status_type::success);
    BOOST_CHECK(input == bytes.data() + encoded_size);
    BOOST_CHECK_EQUAL(*input, 0xa5);
    BOOST_CHECK((types::make_modified_sap_verification_key<policy_type, endianness>(filled) == key));

    const std::forward_list<std::uint8_t> sequential(bytes.begin(), bytes.end());
    auto sequential_input = sequential.begin();
    BOOST_REQUIRE(filled.read(sequential_input, bytes.size()) == status_type::success);
    BOOST_REQUIRE(sequential_input != sequential.end());
    BOOST_CHECK_EQUAL(*sequential_input, 0xa5);
    BOOST_CHECK((types::make_modified_sap_verification_key<policy_type, endianness>(filled) == key));
}

BOOST_AUTO_TEST_CASE(verification_key_rejects_noncanonical_components) {
    const auto canonical = encode(types::fill_modified_sap_verification_key<policy_type, endianness>(test_key()));
    // All sixty GT components, then the independent Fp digest.
    for (std::size_t component = 0; component <= 60; ++component) {
        const std::size_t offset = component == 60 ? digest_offset : 192 + 32 * component;
        BOOST_TEST_CONTEXT("component at byte " << offset) {
            for (const integral_type &raw :
                 {integral_type(base_field_type::modulus), integral_type(base_field_type::modulus + 1)}) {
                auto bytes = canonical;
                write_integer<endianness>(bytes, offset, raw);
                check_invalid_encoding(bytes);
            }
            for (std::uint8_t padding : {0x40, 0x80}) {
                auto bytes = canonical;
                bytes[offset] |= padding;
                check_invalid_encoding(bytes);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(verification_key_rejects_malformed_g2_encodings) {
    const auto canonical = encode(types::fill_modified_sap_verification_key<policy_type, endianness>(test_key()));
    for (std::size_t offset : {0, 64, 128}) {
        BOOST_TEST_CONTEXT("G2 at byte " << offset) {
            for (std::size_t component : {0, 32}) {
                auto bytes = canonical;
                write_integer<endianness>(bytes, offset + component, integral_type(base_field_type::modulus));
                check_invalid_encoding(bytes);
            }
            auto bytes = canonical;
            bytes[offset + 32] |= 0x80;    // The second coordinate component has padding bits, not flags.
            check_invalid_encoding(bytes);
            bytes = canonical;
            bytes[offset] |= 0x80;    // Infinity cannot carry the finite point's nonzero coordinates.
            check_invalid_encoding(bytes);
        }
    }
}

BOOST_AUTO_TEST_CASE(verification_key_rejects_invalid_dimensions_and_g2_parameters) {
    const auto original = test_key();
    const auto original_filled = types::fill_modified_sap_verification_key<policy_type, endianness>(original);
    auto key = original;
    auto filled = original_filled;
    key.num_variables = 0;
    std::get<5>(filled.value()).value() = 0;
    check_invalid_key(key, filled);

    constexpr auto max_domain_log = nil::crypto3::algebra::fields::arithmetic_params<scalar_field_type>::s;
    for (std::size_t domain_size :
         {std::size_t(0), std::size_t(1), std::size_t(3), std::size_t(1) << (max_domain_log + 1)}) {
        key = original;
        filled = original_filled;
        key.domain_size = domain_size;
        std::get<6>(filled.value()).value() = domain_size;
        check_invalid_key(key, filled);
    }

    key = original;
    filled = original_filled;
    key.g2_one = scalar_type(2) * point_type::one();
    std::get<0>(filled.value()).value() = key.g2_one;
    check_invalid_key(key, filled);
    const std::array members {&key_type::g2_one, &key_type::tau_g2, &key_type::gamma_inverse_g2};
    for (std::size_t i = 0; i < members.size(); ++i) {
        key = original;
        filled = original_filled;
        key.*members[i] = point_type::zero();
        const std::array fields {&std::get<0>(filled.value()), &std::get<1>(filled.value()),
                                 &std::get<2>(filled.value())};
        fields[i]->value() = point_type::zero();
        check_invalid_key(key, filled);

        key = original;
        filled = original_filled;
        (key.*members[i]).X += point_type::field_type::value_type::one();
        BOOST_REQUIRE(!(key.*members[i]).is_well_formed());
        fields[i]->value() = key.*members[i];
        BOOST_CHECK_THROW((types::fill_modified_sap_verification_key<policy_type, endianness>(key)),
                          std::invalid_argument);
        BOOST_CHECK_THROW((types::make_modified_sap_verification_key<policy_type, endianness>(filled)),
                          std::invalid_argument);
    }
}

BOOST_AUTO_TEST_CASE(verification_key_rejects_zero_and_non_subgroup_gt) {
    auto outside_subgroup = gt_type::one();
    outside_subgroup.data[0].data[0].data[0] = base_type(2);
    BOOST_REQUIRE(!outside_subgroup.is_zero());
    BOOST_REQUIRE(outside_subgroup.pow(scalar_field_type::modulus) != gt_type::one());
    for (std::size_t i = 0; i < 5; ++i) {
        for (const auto &invalid : {gt_type::zero(), outside_subgroup}) {
            BOOST_TEST_CONTEXT("GT member " << i) {
                auto key = test_key();
                auto filled = types::fill_modified_sap_verification_key<policy_type, endianness>(key);
                (i == 0 ? key.alpha_z_vanishing_gt : key.alpha_gt[i - 1]) = invalid;
                auto &field = i == 0 ? std::get<3>(filled.value()) : std::get<4>(filled.value()).value()[i - 1];
                field = gt_codec(invalid);
                check_invalid_key(key, filled);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(verification_key_rejects_invalid_assembled_field_storage) {
    const auto original = types::fill_modified_sap_verification_key<policy_type, endianness>(test_key());
    for (std::size_t count = 0; count < 4; ++count) {
        auto filled = original;
        std::get<4>(filled.value()).value().resize(count);
        BOOST_CHECK(!filled.valid());
        BOOST_CHECK_THROW((types::make_modified_sap_verification_key<policy_type, endianness>(filled)),
                          std::invalid_argument);
    }
    for (std::size_t i = 0; i < 5; ++i) {
        auto filled = original;
        auto &field = i == 0 ? std::get<3>(filled.value()) : std::get<4>(filled.value()).value()[i - 1];
        field = gt_codec();    // A default extension field has no stored components.
        BOOST_CHECK(!filled.valid());
        BOOST_CHECK_THROW((types::make_modified_sap_verification_key<policy_type, endianness>(filled)),
                          std::invalid_argument);

        std::vector<std::uint8_t> raw(384, 0);
        write_integer<endianness>(raw, 11 * 32, integral_type(base_field_type::modulus));
        auto input = raw.begin();
        BOOST_REQUIRE(field.read(input, raw.size()) == status_type::success);
        BOOST_CHECK(!filled.valid());
        BOOST_CHECK_THROW((types::make_modified_sap_verification_key<policy_type, endianness>(filled)),
                          std::invalid_argument);
    }
    auto filled = original;
    std::vector<std::uint8_t> raw(32);
    write_integer<endianness>(raw, 0, integral_type(base_field_type::modulus));
    auto input = raw.begin();
    BOOST_REQUIRE(std::get<7>(filled.value()).read(input, raw.size()) == status_type::success);
    BOOST_CHECK_THROW((types::make_modified_sap_verification_key<policy_type, endianness>(filled)),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(verification_key_rejects_short_buffers) {
    const auto filled = types::fill_modified_sap_verification_key<policy_type, endianness>(test_key());
    const auto bytes = encode(filled);
    for (std::size_t size = 0; size < bytes.size(); ++size) {
        BOOST_TEST_CONTEXT("available bytes " << size) {
            const std::vector<std::uint8_t> truncated(bytes.begin(), bytes.begin() + size);
            auto input = truncated.begin();
            marshalled_type decoded;
            BOOST_CHECK(decoded.read(input, size) == status_type::not_enough_data);
            BOOST_CHECK(input == truncated.begin());
            // Exercise output bounds immediately before each encoded component and dimension boundary.
            if (size % 32 == 31 || size == 0 || size == dimensions_offset + sizeof(std::size_t) - 1 ||
                size == digest_offset - 1 || size == bytes.size() - 1) {
                std::vector<std::uint8_t> short_output(size);
                auto output = short_output.begin();
                BOOST_CHECK(filled.write(output, size) == status_type::buffer_overflow);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(generated_key_and_proof_survive_round_trip) {
    using system_type = scheme_type::constraint_system_type;
    using variable_type = system_type::constraint_type::variable_type;
    using combination_type = system_type::constraint_type::linear_combination_type;
    // Binding row, then x^2 - y = u, satisfied by (u, x, y) = (5, 3, 4).
    const system_type system {3,
                              {{{}, combination_type(-variable_type(0))},
                               {combination_type(variable_type(1)), combination_type(variable_type(2))}}};
    const scalar_type u(5);
    const std::array<std::uint8_t, 32> seed = {};    // Fixed setup seed for reproducible tests only.
    nil::crypto3::random::chacha_urbg<> random_source(seed);
    const auto keys = scheme_type::generate(system, random_source);
    const auto proof = scheme_type::prove(keys.first, u, {u, scalar_type(3), scalar_type(4)});
    BOOST_REQUIRE(scheme_type::verify(keys.second, u, proof));

    const auto restored_key =
        decode(encode(types::fill_modified_sap_verification_key<policy_type, endianness>(keys.second)));
    const auto proof_bytes = encode(types::fill_modified_sap_proof<policy_type::proof_type, endianness>(proof));
    const auto filled_proof =
        test_tools::decode<types::modified_sap_proof<type_base, policy_type::proof_type>>(proof_bytes);
    const auto restored_proof = types::make_modified_sap_proof<policy_type::proof_type, endianness>(filled_proof);
    BOOST_CHECK(restored_key == keys.second);
    BOOST_CHECK(restored_proof == proof);
    BOOST_CHECK(scheme_type::verify(restored_key, u, restored_proof));
    BOOST_CHECK(!scheme_type::verify(restored_key, scalar_type(7), restored_proof));

    // A standalone decoder can validate the digest as Fp, but has no circuit with which to recompute it.
    auto changed_key = restored_key;
    changed_key.circuit_digest += base_type::one();
    const auto changed_restored =
        decode(encode(types::fill_modified_sap_verification_key<policy_type, endianness>(changed_key)));
    BOOST_CHECK(!scheme_type::verify(changed_restored, u, restored_proof));
}
