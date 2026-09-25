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

#define BOOST_TEST_MODULE crypto3_marshalling_modified_sap_constraint_system_test

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <cstdint>
#include <forward_list>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/marshalling/zk/types/modified_sap/constraint_system.hpp>

#include "detail/marshalling.hpp"

namespace {
    namespace types = nil::crypto3::marshalling::types;
    namespace test_tools = nil::crypto3::marshalling::test_tools;
    using test_tools::encode;
    using test_tools::write_integer;
    using field_type = nil::crypto3::algebra::curves::alt_bn128_254::scalar_field_type;
    using value_type = field_type::value_type;
    using system_type = nil::crypto3::zk::snark::modified_sap_constraint_system<field_type>;
    using row_type = system_type::constraint_type;
    using combination_type = system_type::linear_combination_type;
    using variable_type = row_type::variable_type;
    using endiannesses =
        boost::mpl::list<nil::marshalling::option::big_endian, nil::marshalling::option::little_endian>;
    using status_type = nil::marshalling::status_type;
    constexpr std::size_t S = sizeof(std::size_t);
    constexpr std::size_t F = 32;

    template<typename Endianness>
    using codec = types::modified_sap_constraint_system<nil::marshalling::field_type<Endianness>, system_type>;

    template<typename Endianness>
    using coefficient_codec =
        types::detail::validated_field_element<nil::marshalling::field_type<Endianness>, value_type>;

    combination_type combination(std::initializer_list<std::pair<std::size_t, int>> terms) {
        combination_type result;
        for (const auto &[index, coefficient] : terms) {
            result.add_term(variable_type(index), value_type(coefficient));
        }
        return result;
    }

    system_type test_system() {
        // Three witness entries (u, x, y), a binding row, and two copies of x^2 - y = u.
        // The logical row count stays three; no fourth row is serialized as domain padding.
        return {3,
                {{{}, combination({{0, -1}})},
                 {combination({{1, 1}}), combination({{2, 1}})},
                 {combination({{1, 1}}), combination({{2, 1}})}}};
    }

    template<typename Endianness>
    system_type decode(const std::vector<std::uint8_t> &bytes) {
        return types::make_modified_sap_constraint_system<system_type, Endianness>(
            test_tools::decode<codec<Endianness>>(bytes));
    }

    template<typename Endianness>
    void check_invalid_system(const codec<Endianness> &filled) {
        BOOST_CHECK_THROW((types::make_modified_sap_constraint_system<system_type, Endianness>(filled)),
                          std::invalid_argument);
        BOOST_CHECK_THROW(decode<Endianness>(encode(filled)), std::invalid_argument);
    }
}    // namespace

BOOST_AUTO_TEST_CASE_TEMPLATE(constraint_system_wire_format_and_satisfaction, Endianness, endiannesses) {
    const auto system = test_system();
    const auto filled = types::fill_modified_sap_constraint_system<system_type, Endianness>(system);
    const auto bytes = encode(filled);
    // n = 3, q = 3, T = 5: two system counts, two list prefixes per row, and five (index, coefficient) terms.
    BOOST_CHECK_EQUAL(bytes.size(), 2 * S + 2 * 3 * S + (S + F) * 5);
    std::vector<std::uint8_t> expected;
    const auto append = [&expected](auto integer, std::size_t width) {
        const std::size_t offset = expected.size();
        expected.resize(offset + width, 0);
        for (std::size_t i = 0; i < width; ++i) {
            const std::size_t position =
                std::is_same_v<Endianness, nil::marshalling::option::big_endian> ? width - 1 - i : i;
            expected[offset + position] = static_cast<std::uint8_t>(integer & 0xff);
            integer >>= 8;
        }
    };
    append(3, S);    // Witness dimension.
    append(3, S);    // Logical row count.
    append(0, S);    // Binding row: empty a, c = -w[0].
    append(1, S);
    append(0, S);
    append(field_type::integral_type(field_type::modulus - 1), F);
    for (std::size_t row = 0; row < 2; ++row) {
        append(1, S);
        append(1, S);
        append(1, F);
        append(1, S);
        append(2, S);
        append(1, F);
    }
    BOOST_CHECK(bytes == expected);
    const auto restored = decode<Endianness>(bytes);
    BOOST_CHECK(restored == system);
    BOOST_CHECK_EQUAL(restored.constraints.size(), 3);
    BOOST_CHECK(restored.is_satisfied(value_type(5), {value_type(5), value_type(3), value_type(4)}));
    BOOST_CHECK(!restored.is_satisfied(value_type(5), {value_type(5), value_type(3), value_type(5)}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fill_normalizes_an_owned_copy, Endianness, endiannesses) {
    auto source = test_system();
    source.constraints[0] = {combination({{0, 2}, {0, -2}}), combination({{0, -2}, {2, 0}, {0, 1}})};
    source.constraints[1] = {combination({{2, 0}, {1, 2}, {1, -1}}), combination({{2, 2}, {0, 3}, {2, -1}, {0, -3}})};
    source.constraints[2] = source.constraints[1];
    const auto original = source;
    const auto bytes = encode(types::fill_modified_sap_constraint_system<system_type, Endianness>(source));
    BOOST_CHECK(bytes == encode(types::fill_modified_sap_constraint_system<system_type, Endianness>(test_system())));
    const auto restored = decode<Endianness>(bytes);
    // Native combination equality sorts copies; compare term vectors to detect reordering of the caller's data.
    for (std::size_t i = 0; i < source.constraints.size(); ++i) {
        BOOST_CHECK(source.constraints[i].a.terms == original.constraints[i].a.terms);
        BOOST_CHECK(source.constraints[i].c.terms == original.constraints[i].c.terms);
        BOOST_CHECK(restored.constraints[i].a.terms == test_system().constraints[i].a.terms);
        BOOST_CHECK(restored.constraints[i].c.terms == test_system().constraints[i].c.terms);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fill_checks_raw_indices_before_normalization, Endianness, endiannesses) {
    for (std::size_t row : {0, 1}) {
        for (const auto &extra : {combination({{3, 0}}), combination({{3, 1}, {3, -1}})}) {
            auto source = test_system();
            source.constraints[row].a.terms.insert(source.constraints[row].a.terms.end(), extra.terms.begin(),
                                                   extra.terms.end());
            BOOST_CHECK_THROW((types::fill_modified_sap_constraint_system<system_type, Endianness>(source)),
                              std::invalid_argument);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(make_rejects_noncanonical_terms_and_bad_indices, Endianness, endiannesses) {
    const auto original = types::fill_modified_sap_constraint_system<system_type, Endianness>(test_system());
    for (const auto &invalid : {combination({{2, 1}, {1, 1}}), combination({{1, 1}, {1, 1}}),
                                combination({{1, 1}, {1, -1}}), combination({{1, 0}}), combination({{3, 1}})}) {
        auto filled = original;
        std::get<0>(std::get<1>(filled.value()).value()[1].value()) =
            types::fill_linear_combination<combination_type, Endianness, coefficient_codec<Endianness>>(invalid);
        check_invalid_system<Endianness>(filled);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(make_rejects_missing_dimensions_and_bad_binding, Endianness, endiannesses) {
    const auto original = types::fill_modified_sap_constraint_system<system_type, Endianness>(test_system());
    auto filled = original;
    std::get<0>(filled.value()).value() = 0;
    check_invalid_system<Endianness>(filled);
    filled = original;
    std::get<1>(filled.value()).value().clear();
    check_invalid_system<Endianness>(filled);
    for (const auto &invalid_c : {combination({}), combination({{0, -2}}), combination({{1, -1}})}) {
        filled = original;
        std::get<1>(std::get<1>(filled.value()).value()[0].value()) =
            types::fill_linear_combination<combination_type, Endianness, coefficient_codec<Endianness>>(invalid_c);
        check_invalid_system<Endianness>(filled);
    }
    filled = original;
    std::get<0>(std::get<1>(filled.value()).value()[0].value()) =
        types::fill_linear_combination<combination_type, Endianness, coefficient_codec<Endianness>>(
            combination({{0, 1}}));
    check_invalid_system<Endianness>(filled);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(reader_rejects_noncanonical_coefficients, Endianness, endiannesses) {
    const auto canonical = encode(types::fill_modified_sap_constraint_system<system_type, Endianness>(test_system()));
    for (std::size_t offset : {5 * S, 7 * S + F, 9 * S + 2 * F, 11 * S + 3 * F, 13 * S + 4 * F}) {
        for (const auto &raw :
             {field_type::integral_type(field_type::modulus), field_type::integral_type(field_type::modulus + 1)}) {
            auto bytes = canonical;
            write_integer<Endianness>(bytes, offset, raw);
            codec<Endianness> filled;
            auto input = bytes.begin();
            BOOST_CHECK(filled.read(input, bytes.size()) == status_type::invalid_msg_data);
        }
        for (std::uint8_t padding : {0x40, 0x80}) {
            auto bytes = canonical;
            const std::size_t high_byte =
                std::is_same_v<Endianness, nil::marshalling::option::big_endian> ? offset : offset + F - 1;
            bytes[high_byte] |= padding;
            codec<Endianness> filled;
            auto input = bytes.begin();
            BOOST_CHECK(filled.read(input, bytes.size()) == status_type::invalid_msg_data);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(make_checks_assembled_raw_coefficients, Endianness, endiannesses) {
    auto filled = types::fill_modified_sap_constraint_system<system_type, Endianness>(test_system());
    auto &a = std::get<0>(std::get<1>(filled.value()).value()[1].value());
    auto &coefficient = std::get<1>(a.value()[0].value());
    std::vector<std::uint8_t> bytes(F);
    write_integer<Endianness>(bytes, 0, field_type::integral_type(field_type::modulus + 1));
    auto input = bytes.begin();
    // Bypass the system reader to assemble an invalid raw field that would otherwise reduce to one.
    BOOST_REQUIRE(coefficient.read(input, bytes.size()) == status_type::success);
    BOOST_CHECK(!filled.valid());
    BOOST_CHECK_THROW((types::make_modified_sap_constraint_system<system_type, Endianness>(filled)),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(reader_bounds_counts_and_reserves_remaining_prefixes, Endianness, endiannesses) {
    const auto canonical = encode(types::fill_modified_sap_constraint_system<system_type, Endianness>(test_system()));
    for (std::size_t offset : {S, 2 * S, 3 * S, 5 * S + F, 7 * S + 2 * F, 9 * S + 3 * F, 11 * S + 4 * F}) {
        auto bytes = canonical;
        write_integer<Endianness>(bytes, offset, std::numeric_limits<std::size_t>::max());
        codec<Endianness> filled;
        auto input = bytes.begin();
        BOOST_CHECK(filled.read(input, bytes.size()) == status_type::invalid_msg_data);
        bytes = canonical;
        write_integer<Endianness>(bytes, offset, canonical.size());    // Fits a container, but not this byte budget.
        input = bytes.begin();
        BOOST_CHECK(filled.read(input, bytes.size()) == status_type::not_enough_data);
    }
    for (std::size_t row_count : {1, 2}) {
        // Enough bytes for one term, but none for the remaining combination/row prefixes.
        std::vector<std::uint8_t> bytes(3 * S + S + F, 0);
        write_integer<Endianness>(bytes, 0, std::size_t(1));
        write_integer<Endianness>(bytes, S, row_count);
        write_integer<Endianness>(bytes, 2 * S, std::size_t(1));
        codec<Endianness> filled;
        auto input = bytes.begin();
        BOOST_CHECK(filled.read(input, bytes.size()) == status_type::not_enough_data);
        BOOST_CHECK(input == bytes.begin() + 3 * S);    // Reject before reading or allocating the advertised term.
        BOOST_CHECK(std::get<1>(filled.value()).value().empty());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(constraint_system_rejects_short_buffers, Endianness, endiannesses) {
    const auto filled = types::fill_modified_sap_constraint_system<system_type, Endianness>(test_system());
    const auto bytes = encode(filled);
    for (std::size_t size = 0; size < bytes.size(); ++size) {
        BOOST_TEST_CONTEXT("available bytes " << size) {
            const std::vector<std::uint8_t> truncated(bytes.begin(), bytes.begin() + size);
            auto input = truncated.begin();
            codec<Endianness> decoded;
            BOOST_CHECK(decoded.read(input, size) == status_type::not_enough_data);
            std::vector<std::uint8_t> short_output(size);
            auto output = short_output.begin();
            BOOST_CHECK(filled.write(output, size) == status_type::buffer_overflow);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(constraint_system_read_leaves_following_bytes, Endianness, endiannesses) {
    const auto source = test_system();
    auto bytes = encode(types::fill_modified_sap_constraint_system<system_type, Endianness>(source));
    bytes.push_back(0xa5);
    codec<Endianness> filled;
    const auto *input = bytes.data();
    BOOST_REQUIRE(filled.read(input, bytes.size()) == status_type::success);
    BOOST_CHECK(input == bytes.data() + bytes.size() - 1);
    BOOST_CHECK_EQUAL(*input, 0xa5);
    BOOST_CHECK((types::make_modified_sap_constraint_system<system_type, Endianness>(filled) == source));
    const std::forward_list<std::uint8_t> sequential(bytes.begin(), bytes.end());
    auto sequential_input = sequential.begin();
    BOOST_REQUIRE(filled.read(sequential_input, bytes.size()) == status_type::success);
    BOOST_REQUIRE(sequential_input != sequential.end());
    BOOST_CHECK_EQUAL(*sequential_input, 0xa5);
    BOOST_CHECK((types::make_modified_sap_constraint_system<system_type, Endianness>(filled) == source));
}
