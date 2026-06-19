//---------------------------------------------------------------------------//
// Copyright (c) 2024 Mikhail Komarov <nemo@allocin.it>
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

#define BOOST_TEST_MODULE ct_random_engine_test

#include <string>
#include <tuple>
#include <unordered_map>
#include <sstream>
#include <cstdlib>
#include <ctime>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/curves/bls12.hpp>
#include <nil/crypto3/algebra/curves/mnt4.hpp>

#include <nil/crypto3/random/ct_random_device.hpp>

using namespace nil::crypto3;

template<typename FieldParams>
void print_field_element(std::ostream &os, const typename algebra::fields::detail::element_fp<FieldParams> &e) {
    os << std::hex << e.data;
}

template<typename FieldParams>
void print_field_element(std::ostream &os, const typename algebra::fields::detail::element_fp2<FieldParams> &e) {
    os << "[" << e.data[0].data << ", " << e.data[1].data << "]";
}

template<typename FpCurveGroupElement>
void print_fp_curve_group_element(std::ostream &os, const FpCurveGroupElement &e) {
    os << "( " << e.X.data << " : " << e.Y.data << " : " << e.Z.data << " )";
}

template<typename Fp2CurveGroupElement>
void print_fp2_curve_group_element(std::ostream &os, const Fp2CurveGroupElement &e) {
    os << "(" << e.X.data[0].data << " , " << e.X.data[1].data << ") : (" << e.Y.data[0].data << " , "
       << e.Y.data[1].data << ") : (" << e.Z.data[0].data << " , " << e.Z.data[1].data << ")";
}

template<typename Fp3CurveGroupElement>
void print_fp3_curve_group_element(std::ostream &os, const Fp3CurveGroupElement &e) {
    os << "(" << e.X.data[0].data << " , " << e.X.data[1].data << " , " << e.X.data[2].data << ") : ("
       << e.Y.data[0].data << " , " << e.Y.data[1].data << " , " << e.Y.data[2].data << ") : (" << e.Z.data[0].data
       << " , " << e.Z.data[1].data << " , " << e.Z.data[2].data << ")";
}

namespace boost {
    namespace test_tools {
        namespace tt_detail {
            template<typename FieldParams>
            struct print_log_value<typename algebra::fields::detail::element_fp<FieldParams>> {
                void operator()(std::ostream &os, typename algebra::fields::detail::element_fp<FieldParams> const &e) {
                    print_field_element(os, e);
                }
            };

            template<typename FieldParams>
            struct print_log_value<typename algebra::fields::detail::element_fp2<FieldParams>> {
                void operator()(std::ostream &os, typename algebra::fields::detail::element_fp2<FieldParams> const &e) {
                    print_field_element(os, e);
                }
            };

            template<typename FieldParams>
            struct print_log_value<typename algebra::fields::detail::element_fp6_3over2<FieldParams>> {
                void operator()(std::ostream &os,
                                typename algebra::fields::detail::element_fp6_3over2<FieldParams> const &e) {
                    print_field_element(os, e);
                }
            };

            template<typename FieldParams>
            struct print_log_value<typename algebra::fields::detail::element_fp12_2over3over2<FieldParams>> {
                void operator()(std::ostream &os,
                                typename algebra::fields::detail::element_fp12_2over3over2<FieldParams> const &e) {
                    print_field_element(os, e);
                }
            };

            template<>
            struct print_log_value<typename algebra::curves::mnt4<298>::g1_type<>::value_type> {
                void operator()(std::ostream &os, typename algebra::curves::mnt4<298>::g1_type<>::value_type const &e) {
                    print_fp_curve_group_element(os, e);
                }
            };

            template<>
            struct print_log_value<typename algebra::curves::mnt4<298>::g2_type<>::value_type> {
                void operator()(std::ostream &os, typename algebra::curves::mnt4<298>::g2_type<>::value_type const &e) {
                    print_fp2_curve_group_element(os, e);
                }
            };

            template<>
            struct print_log_value<typename algebra::curves::bls12<381>::g1_type<>::value_type> {
                void operator()(std::ostream &os,
                                typename algebra::curves::bls12<381>::g1_type<>::value_type const &e) {
                    print_fp_curve_group_element(os, e);
                }
            };

            template<>
            struct print_log_value<typename algebra::curves::bls12<381>::g2_type<>::value_type> {
                void operator()(std::ostream &os,
                                typename algebra::curves::bls12<381>::g2_type<>::value_type const &e) {
                    print_fp2_curve_group_element(os, e);
                }
            };

            template<template<typename, typename> class P, typename K, typename V>
            struct print_log_value<P<K, V>> {
                void operator()(std::ostream &, P<K, V> const &) {
                }
            };
        }    // namespace tt_detail
    }    // namespace test_tools
}    // namespace boost

BOOST_AUTO_TEST_SUITE(ct_random_device)

BOOST_AUTO_TEST_CASE(ct_random_device_test) {
    constexpr uint32_t random_value =
        nil::crypto3::random::ct_lcg<std::size_t, 1664525, 1013904223, 4294967296>::value(12345);

    BOOST_CHECK(random_value < 4294967296);
}

BOOST_AUTO_TEST_SUITE_END()