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

#define BOOST_TEST_MODULE crypto3_marshalling_modified_sap_proving_key_test

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <forward_list>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <nil/crypto3/marshalling/zk/types/modified_sap/proving_key.hpp>
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
    using point_type = curve_type::g1_type<>::value_type;
    using policy_type = snark::modified_sap_policy<curve_type,
                                                   snark::modified_sap_bn254_exact_pairing_policy,
                                                   snark::modified_sap_bn254_poseidon_transcript_policy>;
    using scheme_type = snark::modified_sap_snark<policy_type>;
    using key_type = policy_type::proving_key_type;
    using system_type = key_type::constraint_system_type;
    using combination_type = system_type::linear_combination_type;
    using variable_type = system_type::constraint_type::variable_type;
    using endianness = nil::marshalling::option::big_endian;
    using type_base = nil::marshalling::field_type<endianness>;
    using marshalled_type = types::modified_sap_proving_key<type_base, key_type>;
    using point_codec = types::curve_element<type_base, key_type::g1_type>;
    using coefficient_codec = types::detail::validated_field_element<type_base, scalar_type>;
    using status_type = nil::marshalling::status_type;
    constexpr std::size_t S = sizeof(std::size_t);
    constexpr std::size_t verification_key_size = 2144 + 2 * S;
    constexpr std::array query_members = {&key_type::W,   &key_type::H_query, &key_type::T_A,
                                          &key_type::T_C, &key_type::T_H,     &key_type::T_Z};

    combination_type combination(std::initializer_list<std::pair<std::size_t, int>> terms) {
        combination_type result;
        for (const auto &[index, coefficient] : terms) {
            result.add_term(variable_type(index), scalar_type(coefficient));
        }
        return result;
    }

    system_type test_system(std::size_t row_count = 3) {
        // Witness (u, x, y) = (5, 3, 4): a binding row, followed by copies of x^2 - y = u.
        // A binding-only system needs just the single witness entry u.
        system_type system {row_count == 1 ? std::size_t(1) : std::size_t(3), {{{}, combination({{0, -1}})}}};
        for (std::size_t row = 1; row < row_count; ++row) {
            system.constraints.push_back({combination({{1, 1}}), combination({{2, 1}})});
        }
        return system;
    }

    // A shape-valid fixture with distinct query points; actual generated setup keys are tested separately.
    key_type test_key(std::size_t row_count = 3) {
        key_type key;
        key.constraint_system = test_system(row_count);
        auto &vk = key.verification_key;
        vk.tau_g2 = key_type::verification_key_type::g2_value_type::one();
        vk.gamma_inverse_g2 = key_type::verification_key_type::g2_value_type::one();
        vk.num_variables = key.constraint_system.witness_size;
        vk.domain_size = row_count <= 2 ? 2 : 4;
        vk.circuit_digest = policy_type::transcript_policy_type::circuit_digest(key.constraint_system);
        const std::size_t m = vk.domain_size;
        const std::array lengths = {vk.num_variables, m - 1, m - 1, m - 1, m - 2, m};
        std::size_t scalar = 0;
        for (std::size_t i = 0; i < query_members.size(); ++i) {
            for (std::size_t j = 0; j < lengths[i]; ++j) {
                (key.*query_members[i]).push_back(scalar_type(scalar++) * point_type::one());
            }
        }
        return key;
    }

    auto query_fields(marshalled_type &key) {
        auto &members = key.value();
        return std::array {&std::get<0>(members), &std::get<1>(members), &std::get<2>(members),
                           &std::get<3>(members), &std::get<4>(members), &std::get<5>(members)};
    }

    std::array<std::size_t, 6> query_offsets(const marshalled_type &key) {
        std::array<std::size_t, 6> offsets;
        std::size_t offset = 0;
        std::size_t index = 0;
        nil::marshalling::processing::tuple_for_each_until<6>(key.value(), [&](const auto &query) {
            offsets[index++] = offset;
            offset += query.length();
        });
        return offsets;
    }

    key_type decode(const std::vector<std::uint8_t> &bytes) {
        return types::make_modified_sap_proving_key<policy_type, endianness>(
            test_tools::decode<marshalled_type>(bytes));
    }

    void check_invalid_key(const key_type &key, const marshalled_type &filled) {
        BOOST_CHECK_THROW((types::fill_modified_sap_proving_key<policy_type, endianness>(key)), std::invalid_argument);
        BOOST_CHECK_THROW((types::make_modified_sap_proving_key<policy_type, endianness>(filled)),
                          std::invalid_argument);
        BOOST_CHECK_THROW(decode(encode(filled)), std::invalid_argument);
    }
}    // namespace

BOOST_AUTO_TEST_CASE(proving_key_wire_format_and_round_trip) {
    const auto key = test_key();
    const auto filled = types::fill_modified_sap_proving_key<policy_type, endianness>(key);
    const auto bytes = encode(filled);
    // n = 3, m = 4, q = 3, T = 5: six query counts, eighteen G1 points, system, then verification key.
    const std::size_t system_size = 2 * S + 2 * 3 * S + (S + 32) * 5;
    BOOST_CHECK_EQUAL(bytes.size(), 6 * S + 32 * (3 + 5 * 4 - 5) + system_size + verification_key_size);
    std::vector<std::uint8_t> expected;
    for (const auto member : query_members) {
        const auto &query = key.*member;
        const std::size_t offset = expected.size();
        expected.resize(offset + S, 0);
        for (std::size_t i = 0; i < S; ++i) {
            expected[offset + S - 1 - i] = static_cast<std::uint8_t>(query.size() >> (8 * i));
        }
        for (const auto &point : query) {
            const auto point_bytes = encode(point_codec(point));
            expected.insert(expected.end(), point_bytes.begin(), point_bytes.end());
        }
    }
    const auto system_bytes =
        encode(types::fill_modified_sap_constraint_system<system_type, endianness>(key.constraint_system));
    const auto vk_bytes =
        encode(types::fill_modified_sap_verification_key<policy_type, endianness>(key.verification_key));
    expected.insert(expected.end(), system_bytes.begin(), system_bytes.end());
    expected.insert(expected.end(), vk_bytes.begin(), vk_bytes.end());
    BOOST_CHECK(bytes == expected);
    BOOST_CHECK(decode(bytes) == key);
    BOOST_REQUIRE(key.W.front().is_zero());    // Identity query elements have their usual compressed encoding.
    BOOST_CHECK_EQUAL(bytes[S], 0x80);
}

BOOST_AUTO_TEST_CASE(proving_key_fill_normalizes_an_owned_system) {
    auto key = test_key();
    key.constraint_system.constraints[0].c = combination({{0, -2}, {2, 0}, {0, 1}});
    key.constraint_system.constraints[1].a = combination({{2, 0}, {1, 2}, {1, -1}});
    const auto original = key;
    const auto bytes = encode(types::fill_modified_sap_proving_key<policy_type, endianness>(key));
    BOOST_CHECK(bytes == encode(types::fill_modified_sap_proving_key<policy_type, endianness>(test_key())));
    for (std::size_t i = 0; i < key.constraint_system.constraints.size(); ++i) {
        BOOST_CHECK(key.constraint_system.constraints[i].a.terms == original.constraint_system.constraints[i].a.terms);
        BOOST_CHECK(key.constraint_system.constraints[i].c.terms == original.constraint_system.constraints[i].c.terms);
    }
}

BOOST_AUTO_TEST_CASE(proving_key_rejects_wrong_query_lengths) {
    const auto original = test_key();
    const auto canonical = types::fill_modified_sap_proving_key<policy_type, endianness>(original);
    for (std::size_t i = 0; i < query_members.size(); ++i) {
        for (bool shorter : {true, false}) {
            auto key = original;
            auto filled = canonical;
            auto &query = key.*query_members[i];
            auto &marshalled_query = query_fields(filled)[i]->value();
            if (shorter) {
                query.pop_back();
                marshalled_query.pop_back();
            } else {
                query.push_back(point_type::zero());
                marshalled_query.emplace_back(point_type::zero());
            }
            check_invalid_key(key, filled);
        }
    }
}

BOOST_AUTO_TEST_CASE(proving_key_rejects_inconsistent_metadata_and_digest) {
    const auto original = test_key();
    const auto canonical = types::fill_modified_sap_proving_key<policy_type, endianness>(original);
    for (std::size_t mutation = 0; mutation < 4; ++mutation) {
        auto key = original;
        auto filled = canonical;
        auto &vk = std::get<7>(filled.value()).value();
        if (mutation == 0) {
            ++key.verification_key.num_variables;
            ++std::get<5>(vk).value();
        } else if (mutation == 1) {
            key.verification_key.domain_size = 8;    // Supported, but inconsistent with three logical rows.
            std::get<6>(vk).value() = 8;
        } else if (mutation == 2) {
            key.verification_key.circuit_digest += base_type::one();
            std::get<7>(vk) =
                types::detail::validated_field_element<type_base, base_type>(key.verification_key.circuit_digest);
        } else {
            key.constraint_system.constraints[1].a.terms.front().coeff = scalar_type(2);
            auto &row = std::get<1>(std::get<6>(filled.value()).value()).value()[1];
            std::get<1>(std::get<0>(row.value()).value().front().value()) = coefficient_codec(scalar_type(2));
        }
        check_invalid_key(key, filled);
    }
}

BOOST_AUTO_TEST_CASE(proving_key_make_rejects_noncanonical_nested_system) {
    const auto canonical = types::fill_modified_sap_proving_key<policy_type, endianness>(test_key());
    for (const auto &invalid : {combination({{1, 2}, {1, -1}}), combination({{2, 0}, {1, 1}})}) {
        auto filled = canonical;
        auto &row = std::get<1>(std::get<6>(filled.value()).value()).value()[1];
        std::get<0>(row.value()) =
            types::fill_linear_combination<combination_type, endianness, coefficient_codec>(invalid);
        // Both combinations represent the original x: matching the digest cannot excuse a noncanonical payload.
        BOOST_CHECK_THROW((types::make_modified_sap_proving_key<policy_type, endianness>(filled)),
                          std::invalid_argument);
        BOOST_CHECK_THROW(decode(encode(filled)), std::invalid_argument);
    }
}

BOOST_AUTO_TEST_CASE(proving_key_rejects_invalid_query_points) {
    const auto original = test_key();
    const auto canonical = types::fill_modified_sap_proving_key<policy_type, endianness>(original);
    const auto canonical_bytes = encode(canonical);
    const auto offsets = query_offsets(canonical);
    const point_type malformed(base_type::one(), base_type::zero(), base_type::one());
    BOOST_REQUIRE(!malformed.is_well_formed());
    for (std::size_t i = 0; i < query_members.size(); ++i) {
        auto key = original;
        (key.*query_members[i]).front() = malformed;
        BOOST_CHECK_THROW((types::fill_modified_sap_proving_key<policy_type, endianness>(key)), std::invalid_argument);
        auto filled = canonical;
        query_fields(filled)[i]->value().front().value() = malformed;
        BOOST_CHECK_THROW((types::make_modified_sap_proving_key<policy_type, endianness>(filled)),
                          std::invalid_argument);
        for (bool bad_flags : {false, true}) {
            auto bytes = canonical_bytes;
            if (bad_flags) {
                bytes[offsets[i] + S] = 0xc0;    // Infinity cannot also carry a sign bit.
            } else {
                write_integer<endianness>(
                    bytes, offsets[i] + S, base_field_type::integral_type(base_field_type::modulus));
            }
            marshalled_type decoded;
            auto input = bytes.begin();
            BOOST_CHECK(decoded.read(input, bytes.size()) == status_type::invalid_msg_data);
        }
    }
}

BOOST_AUTO_TEST_CASE(proving_key_reuses_nested_verification_key_validation) {
    const auto original = test_key();
    const auto canonical = types::fill_modified_sap_proving_key<policy_type, endianness>(original);
    auto key = original;
    auto filled = canonical;
    key.verification_key.tau_g2 = key_type::verification_key_type::g2_value_type::zero();
    std::get<1>(std::get<7>(filled.value()).value()).value() = key.verification_key.tau_g2;
    check_invalid_key(key, filled);
    key = original;
    filled = canonical;
    key.verification_key.alpha_z_vanishing_gt = key_type::verification_key_type::gt_value_type::zero();
    std::get<3>(std::get<7>(filled.value()).value()) =
        types::detail::validated_field_element<type_base, key_type::verification_key_type::gt_value_type>(
            key.verification_key.alpha_z_vanishing_gt);
    check_invalid_key(key, filled);

    // Even a directly assembled raw digest must be rejected before conversion reduces it modulo p.
    filled = canonical;
    std::vector<std::uint8_t> raw(32);
    write_integer<endianness>(raw, 0, base_field_type::integral_type(base_field_type::modulus));
    auto input = raw.begin();
    BOOST_REQUIRE(std::get<7>(std::get<7>(filled.value()).value()).read(input, raw.size()) == status_type::success);
    BOOST_CHECK_THROW((types::make_modified_sap_proving_key<policy_type, endianness>(filled)), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(proving_key_reader_checks_nested_encodings) {
    const auto filled = types::fill_modified_sap_proving_key<policy_type, endianness>(test_key());
    const auto canonical = encode(filled);
    const std::size_t vk_offset = canonical.size() - verification_key_size;
    const std::size_t system_offset = vk_offset - std::get<6>(filled.value()).length();
    for (std::size_t offset : {system_offset + 5 * S, vk_offset + 192, canonical.size() - 32}) {
        auto bytes = canonical;
        if (offset == system_offset + 5 * S) {
            write_integer<endianness>(bytes, offset, scalar_field_type::integral_type(scalar_field_type::modulus));
        } else {
            write_integer<endianness>(bytes, offset, base_field_type::integral_type(base_field_type::modulus));
        }
        marshalled_type decoded;
        auto input = bytes.begin();
        BOOST_CHECK(decoded.read(input, bytes.size()) == status_type::invalid_msg_data);
    }
}

BOOST_AUTO_TEST_CASE(proving_key_reader_bounds_counts_and_preserves_nested_bytes) {
    const auto filled = types::fill_modified_sap_proving_key<policy_type, endianness>(test_key());
    const auto canonical = encode(filled);
    for (std::size_t offset : query_offsets(filled)) {
        for (std::size_t count : {std::numeric_limits<std::size_t>::max(), canonical.size()}) {
            auto bytes = canonical;
            write_integer<endianness>(bytes, offset, count);
            marshalled_type decoded;
            auto input = bytes.begin();
            const auto expected =
                count == canonical.size() ? status_type::not_enough_data : status_type::invalid_msg_data;
            BOOST_CHECK(decoded.read(input, bytes.size()) == expected);
        }
    }
    // One query point would leave one byte too few for the later prefixes and nested objects.
    std::vector<std::uint8_t> bytes(8 * S + verification_key_size + 32 - 1, 0);
    write_integer<endianness>(bytes, 0, std::size_t(1));
    marshalled_type decoded;
    auto input = bytes.begin();
    BOOST_CHECK(decoded.read(input, bytes.size()) == status_type::not_enough_data);
    BOOST_CHECK(input == bytes.begin() + S);
    BOOST_CHECK(std::get<0>(decoded.value()).value().empty());

    // Six empty queries, then n = q = 1, empty a and one advertised c term with no term bytes.
    // The system reader must not consume bytes reserved for the following verification key.
    bytes.assign(10 * S + verification_key_size, 0);
    write_integer<endianness>(bytes, 6 * S, std::size_t(1));
    write_integer<endianness>(bytes, 7 * S, std::size_t(1));
    write_integer<endianness>(bytes, 9 * S, std::size_t(1));
    input = bytes.begin();
    BOOST_CHECK(decoded.read(input, bytes.size()) == status_type::not_enough_data);
    BOOST_CHECK(input == bytes.begin() + 10 * S);
}

BOOST_AUTO_TEST_CASE(proving_key_rejects_short_buffers) {
    auto key = test_key(1);
    for (auto member : query_members) {
        std::fill((key.*member).begin(), (key.*member).end(), point_type::zero());
    }
    const auto filled = types::fill_modified_sap_proving_key<policy_type, endianness>(key);
    const auto bytes = encode(filled);
    BOOST_REQUIRE(key.T_H.empty());
    for (std::size_t size = 0; size < bytes.size(); ++size) {
        BOOST_TEST_CONTEXT("available bytes " << size) {
            const std::vector<std::uint8_t> truncated(bytes.begin(), bytes.begin() + size);
            auto input = truncated.begin();
            marshalled_type decoded;
            BOOST_CHECK(decoded.read(input, size) == status_type::not_enough_data);
            std::vector<std::uint8_t> output_bytes(size);
            auto output = output_bytes.begin();
            BOOST_CHECK(filled.write(output, size) == status_type::buffer_overflow);
        }
    }
}

BOOST_AUTO_TEST_CASE(proving_key_read_leaves_following_bytes) {
    const auto key = test_key();
    auto bytes = encode(types::fill_modified_sap_proving_key<policy_type, endianness>(key));
    bytes.push_back(0xa5);
    marshalled_type decoded;
    const auto *input = bytes.data();
    BOOST_REQUIRE(decoded.read(input, bytes.size()) == status_type::success);
    BOOST_CHECK(input == bytes.data() + bytes.size() - 1);
    BOOST_CHECK_EQUAL(*input, 0xa5);
    BOOST_CHECK((types::make_modified_sap_proving_key<policy_type, endianness>(decoded) == key));
    const std::forward_list<std::uint8_t> sequential(bytes.begin(), bytes.end());
    auto sequential_input = sequential.begin();
    BOOST_REQUIRE(decoded.read(sequential_input, bytes.size()) == status_type::success);
    BOOST_REQUIRE(sequential_input != sequential.end());
    BOOST_CHECK_EQUAL(*sequential_input, 0xa5);
    BOOST_CHECK((types::make_modified_sap_proving_key<policy_type, endianness>(decoded) == key));
}

BOOST_AUTO_TEST_CASE(generated_proving_keys_remain_usable_after_round_trip) {
    const std::array<std::uint8_t, 32> seed = {};    // Fixed setup seed for reproducible tests only.
    nil::crypto3::random::chacha_urbg<> random_source(seed);
    const scalar_type u(5);
    for (std::size_t rows : {1, 2, 3}) {
        const auto keys = scheme_type::generate(test_system(rows), random_source);
        const auto filled = types::fill_modified_sap_proving_key<policy_type, endianness>(keys.first);
        const auto restored = decode(encode(filled));
        BOOST_CHECK(restored == keys.first);
        BOOST_CHECK(restored.verification_key == keys.second);
        BOOST_CHECK_EQUAL(restored.constraint_system.constraints.size(), rows);
        BOOST_CHECK_EQUAL(restored.T_H.size(), rows <= 2 ? 0 : 2);
        if (rows <= 2) {
            // The empty T_H query still occupies its size_t zero-count prefix.
            BOOST_CHECK_EQUAL(std::get<4>(filled.value()).length(), S);
        }
        const std::vector<scalar_type> witness =
            rows == 1 ? std::vector<scalar_type> {u} : std::vector<scalar_type> {u, scalar_type(3), scalar_type(4)};
        const auto proof = scheme_type::prove(restored, u, witness);
        BOOST_CHECK(proof == scheme_type::prove(keys.first, u, witness));
        BOOST_CHECK(scheme_type::verify(restored.verification_key, u, proof));
        BOOST_CHECK(!scheme_type::verify(restored.verification_key, scalar_type(7), proof));
    }
}
