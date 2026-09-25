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

#define BOOST_TEST_MODULE crypto3_marshalling_r1cs_test

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/marshalling/zk/types/r1cs_gg_ppzksnark/r1cs.hpp>
#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/modified_sap.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace types = nil::crypto3::marshalling::types;
    namespace snark = nil::crypto3::zk::snark;
    using field_type = nil::crypto3::algebra::curves::alt_bn128_254::scalar_field_type;
    using value_type = field_type::value_type;
    using variable_type = math::linear_variable<field_type>;
    using implicit_combination = math::linear_combination<variable_type>;
    using explicit_combination = snark::modified_sap_constraint<field_type>::linear_combination_type;
    using endiannesses =
        boost::mpl::list<nil::marshalling::option::big_endian, nil::marshalling::option::little_endian>;
    using status_type = nil::marshalling::status_type;

    template<typename Marshalled>
    std::vector<std::uint8_t> encode(const Marshalled &value) {
        std::vector<std::uint8_t> bytes(value.length());
        auto output = bytes.begin();
        BOOST_REQUIRE(value.write(output, bytes.size()) == status_type::success);
        BOOST_CHECK(output == bytes.end());
        return bytes;
    }

    template<typename Marshalled>
    Marshalled decode(const std::vector<std::uint8_t> &bytes) {
        Marshalled value;
        auto input = bytes.begin();
        BOOST_REQUIRE(value.read(input, bytes.size()) == status_type::success);
        BOOST_CHECK(input == bytes.end());
        return value;
    }

    template<typename Combination, typename Endianness>
    void check_term_preservation() {
        using type_base = nil::marshalling::field_type<Endianness>;
        using marshalled_type = types::linear_combination<type_base, Combination>;
        const Combination empty;
        const auto empty_bytes = encode(types::fill_linear_combination<Combination, Endianness>(empty));
        BOOST_CHECK(empty_bytes == std::vector<std::uint8_t>(sizeof(std::size_t), 0));
        const auto restored_empty =
            types::make_linear_combination<Combination, Endianness>(decode<marshalled_type>(empty_bytes));
        BOOST_CHECK(restored_empty.terms.empty());

        // The generic codec preserves raw terms; system-level validation can then inspect their original form.
        Combination raw;
        raw.add_term(variable_type(3), value_type(2));
        raw.add_term(variable_type(0), value_type(5));
        raw.add_term(variable_type(3), -value_type(2));
        raw.add_term(variable_type(1), value_type::zero());
        const auto original_terms = raw.terms;
        const auto bytes = encode(types::fill_linear_combination<Combination, Endianness>(raw));
        const auto restored = types::make_linear_combination<Combination, Endianness>(decode<marshalled_type>(bytes));
        // Native combination equality sorts terms, so compare the vectors to check order as well.
        BOOST_CHECK(restored.terms == original_terms);
        BOOST_CHECK(raw.terms == original_terms);
    }
}    // namespace

BOOST_AUTO_TEST_CASE_TEMPLATE(linear_combinations_preserve_layout_and_wire_format, Endianness, endiannesses) {
    using type_base = nil::marshalling::field_type<Endianness>;
    using implicit_marshalled_type = types::linear_combination<type_base, implicit_combination>;
    using explicit_marshalled_type = types::linear_combination<type_base, explicit_combination>;
    static_assert(std::is_same_v<implicit_marshalled_type, explicit_marshalled_type>);

    implicit_combination implicit;
    implicit.add_term(variable_type(0), value_type(2));
    implicit.add_term(variable_type(1), value_type(3));
    explicit_combination explicit_value;
    explicit_value.terms = implicit.terms;
    const auto implicit_bytes = encode(types::fill_linear_combination<implicit_combination, Endianness>(implicit));
    const auto explicit_bytes =
        encode(types::fill_linear_combination<explicit_combination, Endianness>(explicit_value));

    // Fixed wire fixture: size_t count 2, followed by (size_t index, 32-byte scalar) terms (0, 2) and (1, 3).
    std::vector<std::uint8_t> expected;
    const auto append_small_integer = [&expected](std::uint8_t value, std::size_t width) {
        const std::size_t offset = expected.size();
        expected.resize(offset + width, 0);
        if constexpr (std::is_same_v<Endianness, nil::marshalling::option::big_endian>) {
            expected[offset + width - 1] = value;
        } else {
            expected[offset] = value;
        }
    };
    append_small_integer(2, sizeof(std::size_t));
    append_small_integer(0, sizeof(std::size_t));
    append_small_integer(2, 32);
    append_small_integer(1, sizeof(std::size_t));
    append_small_integer(3, 32);
    BOOST_CHECK(implicit_bytes == expected);
    BOOST_CHECK(explicit_bytes == expected);

    const auto restored_implicit = types::make_linear_combination<implicit_combination, Endianness>(
        decode<implicit_marshalled_type>(implicit_bytes));
    const auto restored_explicit = types::make_linear_combination<explicit_combination, Endianness>(
        decode<explicit_marshalled_type>(explicit_bytes));
    BOOST_CHECK(restored_implicit.terms == implicit.terms);
    BOOST_CHECK(restored_explicit.terms == explicit_value.terms);

    const std::vector<value_type> assignment = {value_type(5), value_type(7)};
    // Implicit index zero contributes 2 * 1; explicit index zero contributes 2 * assignment[0].
    BOOST_CHECK(restored_implicit.evaluate(assignment) == value_type(17));    // 2 * 1 + 3 * 5
    BOOST_CHECK(restored_explicit.evaluate(assignment) == value_type(31));    // 2 * 5 + 3 * 7
}

BOOST_AUTO_TEST_CASE_TEMPLATE(linear_combinations_preserve_empty_and_raw_terms, Endianness, endiannesses) {
    check_term_preservation<implicit_combination, Endianness>();
    check_term_preservation<explicit_combination, Endianness>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ordinary_r1cs_round_trip_preserves_implicit_constant, Endianness, endiannesses) {
    using system_type = snark::r1cs_constraint_system<field_type>;
    using type_base = nil::marshalling::field_type<Endianness>;
    using marshalled_type = types::r1cs_constraint_system<type_base, system_type>;

    // (u + 1) * x = y: source indices are 0 = one, 1 = u, 2 = x, 3 = y.
    system_type system;
    system.primary_input_size = 1;
    system.auxiliary_input_size = 2;
    snark::r1cs_constraint<field_type> constraint;
    constraint.a.add_term(variable_type(0), value_type::one());
    constraint.a.add_term(variable_type(1), value_type::one());
    constraint.b.add_term(variable_type(2), value_type::one());
    constraint.c.add_term(variable_type(3), value_type::one());
    system.add_constraint(constraint);
    BOOST_REQUIRE(system.is_valid());
    BOOST_REQUIRE(system.is_satisfied({value_type(3)}, {value_type(2), value_type(8)}));

    const auto bytes = encode(types::fill_r1cs_constraint_system<system_type, Endianness>(system));
    // Three system counts, three term counts, four indices, and four scalar coefficients.
    BOOST_CHECK_EQUAL(bytes.size(), 10 * sizeof(std::size_t) + 4 * 32);
    const auto restored = types::make_r1cs_constraint_system<system_type, Endianness>(decode<marshalled_type>(bytes));
    BOOST_CHECK(restored == system);
    BOOST_REQUIRE(restored.is_valid());
    BOOST_CHECK(restored.is_satisfied({value_type(3)}, {value_type(2), value_type(8)}));
    BOOST_CHECK(!restored.is_satisfied({value_type(5)}, {value_type(2), value_type(8)}));
}
