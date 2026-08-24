//---------------------------------------------------------------------------//
// Copyright (c) 2026
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

#define BOOST_TEST_MODULE field_order_test

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/field_order.hpp>
#include <nil/crypto3/algebra/fields/fp2.hpp>
#include <nil/crypto3/algebra/fields/fp6_3over2.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

namespace {
    namespace fields = nil::crypto3::algebra::fields;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq2_field_type = fields::fp2<fq_field_type>;
    using fq6_field_type = fields::fp6_3over2<fq_field_type>;
    using fq12_field_type = fields::fp12_2over3over2<fq_field_type>;
}    // namespace

BOOST_AUTO_TEST_SUITE(field_order_test_suite)

BOOST_AUTO_TEST_CASE(prime_field_order_equals_its_characteristic) {
    const boost::multiprecision::cpp_int expected(
        "21888242871839275222246405745257275088696311157297823662689037894645226208583");

    BOOST_CHECK(fields::field_characteristic<fq_field_type>() == expected);
    BOOST_CHECK(fields::field_characteristic<fq12_field_type>() == expected);
    BOOST_CHECK(fields::field_order<fq_field_type>() == expected);
}

BOOST_AUTO_TEST_CASE(extension_field_order_is_characteristic_to_the_extension_degree) {
    const boost::multiprecision::cpp_int characteristic = fields::field_order<fq_field_type>();
    boost::multiprecision::cpp_int expected_fq12_order = 1;
    for (std::size_t i = 0; i < fq12_field_type::arity; ++i) {
        expected_fq12_order *= characteristic;
    }

    BOOST_CHECK(fields::field_order<fq2_field_type>() == characteristic * characteristic);
    BOOST_CHECK(fields::field_order<fq6_field_type::value_type>() == fields::field_order<fq6_field_type>());
    BOOST_CHECK(fields::field_order<fq12_field_type>() == expected_fq12_order);
    BOOST_CHECK(fields::field_characteristic<fq12_field_type::value_type>() == characteristic);
}

BOOST_AUTO_TEST_SUITE_END()
