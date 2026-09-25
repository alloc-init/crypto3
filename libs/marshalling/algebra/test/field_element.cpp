//---------------------------------------------------------------------------//
// Copyright (c) 2018-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
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

#define BOOST_TEST_MODULE crypto3_marshalling_field_element_test

#include <nil/crypto3/marshalling/algebra/types/detail/checked_field_element.hpp>

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <list>
#include <vector>

#include <nil/marshalling/status_type.hpp>
#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/endianness.hpp>

#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>
#include <boost/multiprecision/number.hpp>

#include <nil/crypto3/algebra/random_element.hpp>
#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/curves/bls12.hpp>
#include <nil/crypto3/algebra/curves/detail/marshalling.hpp>

#include <nil/marshalling/algorithms/pack.hpp>
#include <nil/crypto3/marshalling/algebra/types/field_element.hpp>

template<typename TIter>
void print_byteblob(TIter iter_begin, TIter iter_end) {
    for (TIter it = iter_begin; it != iter_end; it++) {
        std::cout << std::hex << int(*it) << std::endl;
    }
}

template<typename T, typename Endianness>
void test_field_element(T val) {

    using namespace nil::crypto3::marshalling;

    std::size_t units_bits = 8;
    using unit_type = unsigned char;
    using field_element_type = types::field_element<nil::marshalling::field_type<Endianness>, T>;

    static_assert(nil::crypto3::algebra::FieldValue<T>);
    static_assert(nil::marshalling::MarshallingFieldElement<field_element_type>);
    static_assert(nil::marshalling::is_compatible<T>::value);

    using inferenced_type = typename nil::marshalling::is_compatible<T>::template type<Endianness>;

    static_assert(std::is_same<inferenced_type, field_element_type>::value);

    nil::marshalling::status_type status;
    std::vector<unit_type> cv = nil::marshalling::pack<Endianness>(val, status);

    BOOST_CHECK(status == nil::marshalling::status_type::success);

    T test_val = nil::marshalling::pack<Endianness>(cv, status);

    BOOST_CHECK(val == test_val);
    BOOST_CHECK(status == nil::marshalling::status_type::success);
}

template<typename FieldType, typename Endianness>
void test_field_element() {
    std::cout << std::hex;
    std::cerr << std::hex;
    for (unsigned i = 0; i < 128; ++i) {
        if (!(i % 16) && i) {
            std::cout << std::dec << i << " tested" << std::endl;
        }
        typename FieldType::value_type val = nil::crypto3::algebra::random_element<FieldType>();
        test_field_element<typename FieldType::value_type, Endianness>(val);
    }
}

template<typename FieldType, typename Endianness>
void test_checked_field_element() {
    namespace types = nil::crypto3::marshalling::types;
    using status_type = nil::marshalling::status_type;
    using type_base = nil::marshalling::field_type<Endianness>;
    using value_type = typename FieldType::value_type;
    using integral_type = typename FieldType::integral_type;
    using field_element_type = types::detail::validated_field_element<type_base, value_type>;
    using encoding_validator = types::detail::canonical_field_element_encoding_validator<type_base, value_type>;
    using component_type = types::integral<type_base, integral_type>;

    const auto check_read = [](const std::vector<std::uint8_t> &bytes, status_type expected_status,
                               const value_type *expected_value = nullptr) {
        const auto check_input = [&](auto input, auto end) {
            field_element_type decoded;
            std::size_t remaining = bytes.size();
            BOOST_REQUIRE(types::detail::read_marshaled_value(decoded, input, remaining, encoding_validator()) ==
                          expected_status);
            BOOST_CHECK(input == end);
            BOOST_CHECK_EQUAL(remaining, 0);
            if (expected_value != nullptr) {
                BOOST_CHECK(decoded.value() == *expected_value);
            }
        };
        check_input(bytes.begin(), bytes.end());
        const std::list<std::uint8_t> sequential_bytes(bytes.begin(), bytes.end());
        check_input(sequential_bytes.begin(), sequential_bytes.end());
    };

    for (const value_type &value : {value_type::zero(), value_type::one(), -value_type::one()}) {
        const field_element_type encoded(value);
        std::vector<std::uint8_t> expected(encoded.length());
        auto direct_output = expected.begin();
        BOOST_REQUIRE(encoded.write(direct_output, expected.size()) == status_type::success);

        std::vector<std::uint8_t> bytes;
        auto output = std::back_inserter(bytes);
        BOOST_REQUIRE(types::detail::write_marshaled_value(encoded, output) == status_type::success);
        BOOST_CHECK(bytes == expected);
        check_read(bytes, status_type::success, &value);
    }

    // Write raw integers so field construction cannot reduce malformed inputs before they reach the reader.
    for (std::size_t component = 0; component < FieldType::arity; ++component) {
        BOOST_TEST_CONTEXT("component " << component) {
            for (int offset : {-1, 0, 1}) {
                std::vector<std::uint8_t> bytes(field_element_type::max_length(), 0);
                integral_type raw_value = FieldType::modulus;
                if (offset < 0) {
                    --raw_value;
                } else if (offset > 0) {
                    ++raw_value;
                }
                const component_type raw_component(raw_value);
                auto output = bytes.begin() + component * component_type::max_length();
                BOOST_REQUIRE(raw_component.write(output, raw_component.length()) == status_type::success);
                check_read(bytes, offset < 0 ? status_type::success : status_type::invalid_msg_data);
            }

            // BN254 components use 254 bits; neither of the two unused high bits may be set.
            for (std::uint8_t padding_bit : {0x40, 0x80}) {
                std::vector<std::uint8_t> bytes(field_element_type::max_length(), 0);
                const std::size_t high_byte_offset =
                    component * component_type::max_length() +
                    (std::is_same_v<Endianness, nil::marshalling::option::little_endian> ?
                         component_type::max_length() - 1 :
                         0);
                bytes[high_byte_offset] = padding_bit;
                check_read(bytes, status_type::invalid_msg_data);
            }
        }
    }

    const std::vector<std::uint8_t> bytes(field_element_type::max_length(), 0);
    const std::list<std::uint8_t> sequential_bytes(bytes.begin(), bytes.end());
    const auto check_truncated_input = [&](auto input, std::size_t available) {
        const auto initial_input = input;
        field_element_type decoded(value_type::one());
        std::size_t remaining = available;
        BOOST_CHECK(types::detail::read_marshaled_value(decoded, input, remaining, encoding_validator()) ==
                    status_type::not_enough_data);
        BOOST_CHECK(input == initial_input);
        BOOST_CHECK_EQUAL(remaining, available);
        BOOST_CHECK(decoded.value() == value_type::one());
    };
    for (std::size_t available : {std::size_t(0), bytes.size() - 1}) {
        check_truncated_input(bytes.begin(), available);
        check_truncated_input(sequential_bytes.begin(), available);
    }
}

using bn254 = nil::crypto3::algebra::curves::alt_bn128_254;
using checked_field_types =
    boost::mpl::list<bn254::scalar_field_type, bn254::base_field_type, bn254::g2_type<>::field_type, bn254::gt_type>;

BOOST_AUTO_TEST_SUITE(field_element_test_suite)

BOOST_AUTO_TEST_CASE_TEMPLATE(checked_field_element_bn254_be, FieldType, checked_field_types) {
    test_checked_field_element<FieldType, nil::marshalling::option::big_endian>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(checked_field_element_bn254_le, FieldType, checked_field_types) {
    test_checked_field_element<FieldType, nil::marshalling::option::little_endian>();
}

BOOST_AUTO_TEST_CASE(field_element_bls12_381_g1_field_be) {
    std::cout << "BLS12-381 g1 group field big-endian test started" << std::endl;
    test_field_element<nil::crypto3::algebra::curves::bls12<381>::g1_type<>::field_type,
                       nil::marshalling::option::big_endian>();
    std::cout << "BLS12-381 g1 group field big-endian test finished" << std::endl;
}

BOOST_AUTO_TEST_CASE(field_element_bls12_381_g1_field_le) {
    std::cout << "BLS12-381 g1 group field little-endian test started" << std::endl;
    test_field_element<nil::crypto3::algebra::curves::bls12<381>::g1_type<>::field_type,
                       nil::marshalling::option::little_endian>();
    std::cout << "BLS12-381 g1 group field little-endian test finished" << std::endl;
}

BOOST_AUTO_TEST_CASE(field_element_bls12_381_g2_field_be) {
    std::cout << "BLS12-381 g2 group field big-endian test started" << std::endl;
    test_field_element<nil::crypto3::algebra::curves::bls12<381>::g2_type<>::field_type,
                       nil::marshalling::option::big_endian>();
    std::cout << "BLS12-381 g2 group field big-endian test finished" << std::endl;
}

BOOST_AUTO_TEST_CASE(field_element_bls12_381_g2_field_le) {
    std::cout << "BLS12-381 g2 group field little-endian test started" << std::endl;
    test_field_element<nil::crypto3::algebra::curves::bls12<381>::g2_type<>::field_type,
                       nil::marshalling::option::little_endian>();
    std::cout << "BLS12-381 g2 group field little-endian test finished" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
