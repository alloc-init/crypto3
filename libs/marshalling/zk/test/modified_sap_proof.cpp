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

#define BOOST_TEST_MODULE crypto3_marshalling_modified_sap_proof_test

#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <forward_list>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <nil/crypto3/marshalling/zk/types/modified_sap/proof.hpp>
#include <nil/crypto3/random/chacha_urbg.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/policy.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/transcript.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap_snark.hpp>

namespace {
    namespace types = nil::crypto3::marshalling::types;
    namespace snark = nil::crypto3::zk::snark;
    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using scalar_field_type = curve_type::scalar_field_type;
    using scalar_type = scalar_field_type::value_type;
    using base_field_type = curve_type::base_field_type;
    using base_type = base_field_type::value_type;
    using point_type = curve_type::g1_type<>::value_type;
    using proof_type = snark::modified_sap_proof<curve_type>;
    using endianness = nil::marshalling::option::big_endian;
    using type_base = nil::marshalling::field_type<endianness>;
    using marshalled_type = types::modified_sap_proof<type_base, proof_type>;
    using status_type = nil::marshalling::status_type;

    proof_type test_proof() {
        return {point_type::one(), -point_type::one(), scalar_type(2), scalar_type(3), scalar_type(5), scalar_type(7)};
    }

    std::vector<std::uint8_t> encode(const proof_type &proof) {
        const auto filled = types::fill_modified_sap_proof<proof_type, endianness>(proof);
        std::vector<std::uint8_t> bytes(filled.length());
        auto output = bytes.begin();
        BOOST_REQUIRE(filled.write(output, bytes.size()) == status_type::success);
        BOOST_CHECK(output == bytes.end());
        return bytes;
    }

    proof_type decode(const std::vector<std::uint8_t> &bytes) {
        marshalled_type filled;
        auto input = bytes.begin();
        BOOST_REQUIRE(filled.read(input, bytes.size()) == status_type::success);
        BOOST_CHECK(input == bytes.end());
        return types::make_modified_sap_proof<proof_type, endianness>(filled);
    }

    void check_invalid_encoding(const std::vector<std::uint8_t> &bytes) {
        marshalled_type filled;
        auto input = bytes.begin();
        BOOST_CHECK(filled.read(input, bytes.size()) == status_type::invalid_msg_data);
        // A failed read must not be followed by make_modified_sap_proof.
    }

    template<typename Integral>
    void write_raw_integer(std::vector<std::uint8_t> &bytes, std::size_t offset, const Integral &value) {
        const types::integral<type_base, Integral> raw(value);
        auto output = bytes.begin() + offset;
        BOOST_REQUIRE(raw.write(output, bytes.size() - offset) == status_type::success);
    }
}    // namespace

BOOST_AUTO_TEST_CASE(proof_wire_format_and_round_trip) {
    static_assert(marshalled_type::min_length() == 192);
    static_assert(marshalled_type::max_length() == 192);
    const auto proof = test_proof();
    const auto bytes = encode(proof);

    // G1 generator (1, 2), its negation (sign bit set), then evaluations 2, 3, 5, 7.
    std::vector<std::uint8_t> expected(192, 0);
    expected[31] = 1;
    expected[32] = 0x40;
    expected[63] = 1;
    expected[95] = 2;
    expected[127] = 3;
    expected[159] = 5;
    expected[191] = 7;
    BOOST_CHECK(bytes == expected);
    BOOST_CHECK(decode(bytes) == proof);
}

BOOST_AUTO_TEST_CASE(proof_identity_and_scalar_boundaries) {
    proof_type proof;
    proof.v_C = -scalar_type::one();
    proof.v_H = scalar_type::one();
    const auto bytes = encode(proof);
    BOOST_CHECK_EQUAL(bytes[0], 0x80);
    BOOST_CHECK_EQUAL(bytes[32], 0x80);
    BOOST_CHECK(decode(bytes) == proof);
}

BOOST_AUTO_TEST_CASE(proof_read_leaves_following_bytes) {
    auto bytes = encode(test_proof());
    bytes.push_back(0xa5);
    marshalled_type filled;
    const auto *input = bytes.data();
    BOOST_REQUIRE(filled.read(input, bytes.size()) == status_type::success);
    BOOST_CHECK(input == bytes.data() + 192);
    BOOST_CHECK_EQUAL(*input, 0xa5);
    BOOST_CHECK((types::make_modified_sap_proof<proof_type, endianness>(filled) == test_proof()));

    // The shared checked-value reader also supports sequential input iterators.
    const std::forward_list<std::uint8_t> sequential(bytes.begin(), bytes.end());
    auto sequential_input = sequential.begin();
    BOOST_REQUIRE(filled.read(sequential_input, bytes.size()) == status_type::success);
    BOOST_REQUIRE(sequential_input != sequential.end());
    BOOST_CHECK_EQUAL(*sequential_input, 0xa5);
    BOOST_CHECK((types::make_modified_sap_proof<proof_type, endianness>(filled) == test_proof()));
}

BOOST_AUTO_TEST_CASE(proof_rejects_noncanonical_scalars) {
    using integral_type = scalar_field_type::integral_type;
    const auto canonical = encode(test_proof());
    for (std::size_t offset : {64, 96, 128, 160}) {
        BOOST_TEST_CONTEXT("scalar at byte " << offset) {
            for (const integral_type &raw :
                 {integral_type(scalar_field_type::modulus), integral_type(scalar_field_type::modulus + 1)}) {
                auto bytes = canonical;
                write_raw_integer(bytes, offset, raw);
                check_invalid_encoding(bytes);
            }
            // Fixed-width integer decoding could discard either of these unused high bits.
            for (std::uint8_t padding : {0x40, 0x80}) {
                auto bytes = canonical;
                bytes[offset] |= padding;
                check_invalid_encoding(bytes);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(proof_rejects_malformed_points) {
    const auto canonical = encode(test_proof());
    base_type nonresidue_x = base_type::zero();
    bool found = false;
    for (std::size_t attempt = 0; attempt < 128; ++attempt) {
        if (!(nonresidue_x.pow(3) + base_type(3)).is_square()) {
            found = true;
            break;
        }
        nonresidue_x += base_type::one();
    }
    BOOST_REQUIRE(found);
    for (std::size_t offset : {0, 32}) {
        BOOST_TEST_CONTEXT("point at byte " << offset) {
            auto bytes = canonical;
            write_raw_integer(bytes, offset, base_field_type::integral_type(base_field_type::modulus));
            check_invalid_encoding(bytes);
            bytes = canonical;
            write_raw_integer(bytes, offset, nonresidue_x.to_integral());
            check_invalid_encoding(bytes);
            bytes = canonical;
            bytes[offset] |= 0x80;    // Both finite fixture points have a nonzero infinity payload.
            check_invalid_encoding(bytes);
        }
    }
}

BOOST_AUTO_TEST_CASE(proof_conversions_reject_invalid_assembled_values) {
    const point_type malformed(base_type::one(), base_type::zero(), base_type::one());
    BOOST_REQUIRE(!malformed.is_well_formed());
    for (auto member : {&proof_type::P, &proof_type::Q}) {
        auto proof = test_proof();
        proof.*member = malformed;
        BOOST_CHECK_THROW((types::fill_modified_sap_proof<proof_type, endianness>(proof)), std::invalid_argument);
    }

    auto filled = types::fill_modified_sap_proof<proof_type, endianness>(test_proof());
    nil::marshalling::processing::tuple_for_each_until<2>(filled.value(), [&](auto &member) {
        const auto original = member;
        member.value() = malformed;
        BOOST_CHECK(!filled.valid());
        BOOST_CHECK_THROW((types::make_modified_sap_proof<proof_type, endianness>(filled)), std::invalid_argument);
        std::vector<std::uint8_t> bytes(filled.length());
        auto output = bytes.begin();
        BOOST_CHECK(filled.write(output, bytes.size()) == status_type::invalid_msg_data);
        member = original;
    });

    std::vector<std::uint8_t> noncanonical(32);
    write_raw_integer(noncanonical, 0, scalar_field_type::integral_type(scalar_field_type::modulus));
    nil::marshalling::processing::tuple_for_each_from<2>(filled.value(), [&](auto &member) {
        const auto original = member;
        auto input = noncanonical.begin();
        // Assemble raw storage through the underlying field reader, bypassing the proof's checked reader.
        BOOST_REQUIRE(member.read(input, noncanonical.size()) == status_type::success);
        BOOST_CHECK(!filled.valid());
        BOOST_CHECK_THROW((types::make_modified_sap_proof<proof_type, endianness>(filled)), std::invalid_argument);
        member = original;
    });
}

BOOST_AUTO_TEST_CASE(proof_rejects_short_buffers) {
    const auto bytes = encode(test_proof());
    const auto filled = types::fill_modified_sap_proof<proof_type, endianness>(test_proof());
    for (std::size_t size = 0; size < bytes.size(); ++size) {
        BOOST_TEST_CONTEXT("available bytes " << size) {
            const std::vector<std::uint8_t> truncated(bytes.begin(), bytes.begin() + size);
            auto input = truncated.begin();
            marshalled_type decoded;
            BOOST_CHECK(decoded.read(input, size) == status_type::not_enough_data);
            BOOST_CHECK(input == truncated.begin());

            std::vector<std::uint8_t> short_output(size);
            auto output = short_output.begin();
            BOOST_CHECK(filled.write(output, size) == status_type::buffer_overflow);
        }
    }
}

BOOST_AUTO_TEST_CASE(generated_proof_survives_round_trip) {
    using policy_type = snark::modified_sap_policy<curve_type,
                                                   snark::modified_sap_bn254_exact_pairing_policy,
                                                   snark::modified_sap_bn254_poseidon_transcript_policy>;
    using scheme_type = snark::modified_sap_snark<policy_type>;
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

    const auto restored = decode(encode(proof));
    BOOST_CHECK(restored == proof);
    BOOST_CHECK(scheme_type::verify(keys.second, u, restored));
    BOOST_CHECK(!scheme_type::verify(keys.second, scalar_type(7), restored));

    // Canonical encoding alone does not establish the proof equations.
    auto altered = proof;
    altered.v_C += scalar_type::one();
    BOOST_CHECK(!scheme_type::verify(keys.second, u, decode(encode(altered))));
}
