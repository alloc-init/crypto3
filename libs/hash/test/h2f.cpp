//---------------------------------------------------------------------------//
// Copyright (c) 2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2021 Ilias Khairullin <ilias@nil.foundation>
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

#define BOOST_TEST_MODULE hash_h2f_test

#include <iostream>
#include <cstdint>
#include <vector>
#include <string>
#include <type_traits>
#include <tuple>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <nil/crypto3/algebra/curves/bls12.hpp>

#include <nil/crypto3/hash/h2f.hpp>
#include <nil/crypto3/hash/detail/h2f/h2f_functions.hpp>
#include <nil/crypto3/hash/algorithm/hash.hpp>
#include <nil/crypto3/hash/type_traits.hpp>

using namespace nil::crypto3;
using namespace nil::crypto3::algebra;

template<typename FieldParams>
void print_field_element(std::ostream &os, const typename fields::detail::element_fp<FieldParams> &e) {
    std::cout << e.data << std::endl;
}

template<typename FieldParams>
void print_field_element(std::ostream &os, const typename fields::detail::element_fp2<FieldParams> &e) {
    std::cout << e.data[0].data << ", " << e.data[1].data << std::endl;
}

namespace boost {
    namespace test_tools {
        namespace tt_detail {
            template<typename FieldParams>
            struct print_log_value<typename fields::detail::element_fp<FieldParams>> {
                void operator()(std::ostream &os, typename fields::detail::element_fp<FieldParams> const &e) {
                    print_field_element(os, e);
                }
            };

            template<typename FieldParams>
            struct print_log_value<typename fields::detail::element_fp2<FieldParams>> {
                void operator()(std::ostream &os, typename fields::detail::element_fp2<FieldParams> const &e) {
                    print_field_element(os, e);
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

template<typename Expander,
        typename DstType,
        typename MsgType,
        typename ResultType,
        typename = typename std::enable_if<std::is_same<std::uint8_t, typename DstType::value_type>::value &&
                                           std::is_same<std::uint8_t, typename MsgType::value_type>::value &&
                                           std::is_same<std::uint8_t, typename ResultType::value_type>::value>::type>
void check_expand_message(const DstType &dst, const MsgType &msg, const ResultType &result) {
    auto result_compare = [&result](auto my_result) {
        if (result.size() != my_result.size()) {
            return false;
        }
        bool ret = true;
        for (std::size_t i = 0; i < result.size(); i++) {
            ret &= result[i] == my_result[i];
        }
        return ret;
    };

    typename Expander::accumulator_type acc;
    Expander::init_accumulator(acc);
    Expander::update(acc, msg);
    typename Expander::result_type uniform_bytes = Expander::process(acc, dst);
    BOOST_CHECK(result_compare(uniform_bytes));
}

template<typename HashType>
typename std::enable_if<hashes::is_h2f<HashType>::value>::type
check_hash_to_field_ro(const std::string &msg_str, const typename HashType::digest_type &result) {

    std::vector<std::uint8_t> msg(msg_str.begin(), msg_str.end());
    typename HashType::digest_type u = hash<HashType>(msg);
    for (std::size_t i = 0; i < HashType::count; i++) {
        BOOST_CHECK_EQUAL(u[i], result[i]);
    }
}

BOOST_AUTO_TEST_SUITE(hash_h2f_manual_tests)

    BOOST_AUTO_TEST_CASE(expand_message_xmd_sha256_test) {
        // https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-10#appendix-K.1
        using hash_type = hashes::sha2<256>;

        std::string DST_str("QUUX-V01-CS02-with-expander");
        std::vector<std::uint8_t> DST(DST_str.begin(), DST_str.end());

        // {BytesLength, msg, uniform_bytes}
        using samples_type = std::vector<std::tuple<std::size_t, std::vector<std::uint8_t>, std::vector<std::uint8_t>>>;
        samples_type samples{
                {0x20, {},                 {0xf6, 0x59, 0x81, 0x9a, 0x64, 0x73, 0xc1, 0x83, 0x5b, 0x25, 0xea, 0x59, 0xe3, 0xd3, 0x89, 0x14,
                                                                                                                                            0xc9, 0x8b, 0x37, 0x4f, 0x9,  0x70, 0xb7, 0xe4, 0xc9, 0x21, 0x81, 0xdf, 0x92, 0x8f, 0xca, 0x88}},
                {0x20, {0x61, 0x62, 0x63}, {0x1c, 0x38, 0xf7, 0xc2, 0x11, 0xef, 0x23, 0x33, 0x67, 0xb2, 0x42,
                                                                                                              0xd,  0x4,  0x79, 0x8f, 0xa4, 0x69, 0x80, 0x80, 0xa8, 0x90, 0x10,
                                                                                                                                                                                0x21, 0xa7, 0x95, 0xa1, 0x15, 0x17, 0x75, 0xfe, 0x4d, 0xa7}},
                {0x20,
                       {0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39},
                                           {0x8f, 0x7e, 0x7b, 0x66, 0x79, 0x1f, 0xd,  0xa0, 0xdb, 0xb5, 0xec, 0x7c, 0x22, 0xec, 0x63, 0x7f,
                                                                                                                                            0x79, 0x75, 0x8c, 0xa,  0x48, 0x17, 0xb,  0xfb, 0x7c, 0x46, 0x11, 0xbd, 0x30, 0x4e, 0xce, 0x89}},
                {0x20,
                       {0x71, 0x31, 0x32, 0x38, 0x5f, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                              0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                    0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                          0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71},
                                           {0x72, 0xd5, 0xaa, 0x5e, 0xc8, 0x10, 0x37, 0xd,  0x1f, 0x0,  0x13, 0xc0, 0xdf, 0x2f, 0x1d, 0x65,
                                                                                                                                            0x69, 0x94, 0x94, 0xee, 0x2a, 0x39, 0xf7, 0x2e, 0x17, 0x16, 0xb1, 0xb9, 0x64, 0xe1, 0xc6, 0x42}},
                {0x20,
                       {0x61, 0x35, 0x31, 0x32, 0x5f, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                    0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                            0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61},
                                           {0x3b, 0x8e, 0x70, 0x4f, 0xc4, 0x83, 0x36, 0xac, 0xa4, 0xc2, 0xa1, 0x21, 0x95, 0xb7, 0x20, 0x88,
                                                                                                                                            0x2f, 0x21, 0x62, 0xa4, 0xb7, 0xb1, 0x3a, 0x9c, 0x35, 0xd,  0xb4, 0x6f, 0x42, 0x9b, 0x77, 0x1b}},
                {0x80, {},                 {0x8b, 0xcf, 0xfd, 0x1a, 0x3c, 0xae, 0x24, 0xcf, 0x9c, 0xd7, 0xab, 0x85, 0x62, 0x8f, 0xd1, 0x11,
                                                                                                                                            0xbb, 0x17, 0xe3, 0x73, 0x9d, 0x3b, 0x53, 0xf8, 0x95, 0x80, 0xd2, 0x17, 0xaa, 0x79, 0x52, 0x6f,
                                                   0x17, 0x8,  0x35, 0x4a, 0x76, 0xa4, 0x2,  0xd3, 0x56, 0x9d, 0x6a, 0x9d, 0x19, 0xef, 0x3d, 0xe4,
                                                   0xd0, 0xb9, 0x91, 0xe4, 0xf5, 0x4b, 0x9f, 0x20, 0xdc, 0xde, 0x9b, 0x95, 0xa6, 0x68, 0x24, 0xcb,
                                                   0xdf, 0x6c, 0x1a, 0x96, 0x3a, 0x19, 0x13, 0xd4, 0x3f, 0xd7, 0xac, 0x44, 0x3a, 0x2,  0xfc, 0x5d,
                                                   0x9d, 0x8d, 0x77, 0xe2, 0x7,  0x1b, 0x86, 0xab, 0x11, 0x4a, 0x9f, 0x34, 0x15, 0x9,  0x54, 0xa7,
                                                   0x53, 0x1d, 0xa5, 0x68, 0xa1, 0xea, 0x8c, 0x76, 0x8,  0x61, 0xc0, 0xcd, 0xe2, 0x0,  0x5a, 0xfc,
                                                   0x2c, 0x11, 0x40, 0x42, 0xee, 0x7b, 0x58, 0x48, 0xf5, 0x30, 0x3f, 0x6,  0x11, 0xcf, 0x29, 0x7f}},
                {0x80,
                       {0x61, 0x62, 0x63},
                                           {0xfe, 0x99, 0x4e, 0xc5, 0x1b, 0xda, 0xa8, 0x21, 0x59, 0x80, 0x47, 0xb3, 0x12, 0x1c, 0x14, 0x9b,
                                                                                                                                            0x36, 0x4b, 0x17, 0x86, 0x6,  0xd5, 0xe7, 0x2b, 0xfb, 0xb7, 0x13, 0x93, 0x3a, 0xcc, 0x29, 0xc1,
                                                   0x86, 0xf3, 0x16, 0xba, 0xec, 0xf7, 0xea, 0x22, 0x21, 0x2f, 0x24, 0x96, 0xef, 0x3f, 0x78, 0x5a,
                                                   0x27, 0xe8, 0x4a, 0x40, 0xd8, 0xb2, 0x99, 0xce, 0xc5, 0x60, 0x32, 0x76, 0x3e, 0xce, 0xef, 0xf4,
                                                   0xc6, 0x1b, 0xd1, 0xfe, 0x65, 0xed, 0x81, 0xde, 0xca, 0xff, 0xf4, 0xa3, 0x1d, 0x1,  0x98, 0x61,
                                                   0x9c, 0xa,  0xa0, 0xc6, 0xc5, 0x1f, 0xca, 0x15, 0x52, 0x7,  0x89, 0x92, 0x5e, 0x81, 0x3d, 0xcf,
                                                   0xd3, 0x18, 0xb5, 0x42, 0xf8, 0x79, 0x94, 0x41, 0x27, 0x1f, 0x4d, 0xb9, 0xee, 0x3b, 0x80, 0x92,
                                                   0xa7, 0xa2, 0xe8, 0xd5, 0xb7, 0x5b, 0x73, 0xe2, 0x8f, 0xb1, 0xab, 0x6b, 0x45, 0x73, 0xc1, 0x92}},
                {0x80,
                       {0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39},
                                           {0xc9, 0xec, 0x79, 0x41, 0x81, 0x1b, 0x1e, 0x19, 0xce, 0x98, 0xe2, 0x1d, 0xb2, 0x8d, 0x22, 0x25,
                                                                                                                                            0x93, 0x54, 0xd4, 0xd0, 0x64, 0x3e, 0x30, 0x11, 0x75, 0xe2, 0xf4, 0x74, 0xe0, 0x30, 0xd3, 0x26,
                                                   0x94, 0xe9, 0xdd, 0x55, 0x20, 0xdd, 0xe9, 0x3f, 0x36, 0x0,  0xd8, 0xed, 0xad, 0x94, 0xe5, 0xc3,
                                                   0x64, 0x90, 0x30, 0x88, 0xa7, 0x22, 0x8c, 0xc9, 0xef, 0xf6, 0x85, 0xd7, 0xea, 0xac, 0x50, 0xd5,
                                                   0xa5, 0xa8, 0x22, 0x9d, 0x8,  0x3b, 0x51, 0xde, 0x4c, 0xcc, 0x37, 0x33, 0x91, 0x7f, 0x4b, 0x95,
                                                   0x35, 0xa8, 0x19, 0xb4, 0x45, 0x81, 0x48, 0x90, 0xb7, 0x2,  0x9b, 0x5d, 0xe8, 0x5,  0xbf, 0x62,
                                                   0xb3, 0x3a, 0x4d, 0xc7, 0xe2, 0x4a, 0xcd, 0xf2, 0xc9, 0x24, 0xe9, 0xfe, 0x50, 0xd5, 0x5a, 0x6b,
                                                   0x83, 0x2c, 0x8c, 0x84, 0xc7, 0xf8, 0x24, 0x74, 0xb3, 0x4e, 0x48, 0xc6, 0xd4, 0x38, 0x67, 0xbe}},
                {0x80,
                       {0x71, 0x31, 0x32, 0x38, 0x5f, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                              0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                    0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                          0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71},
                                           {0x48, 0xe2, 0x56, 0xdd, 0xba, 0x72, 0x20, 0x53, 0xba, 0x46, 0x2b, 0x2b, 0x93, 0x35, 0x1f, 0xc9,
                                                                                                                                            0x66, 0x2,  0x6e, 0x6d, 0x6d, 0xb4, 0x93, 0x18, 0x97, 0x98, 0x18, 0x1c, 0x5f, 0x3f, 0xee, 0xa3,
                                                   0x77, 0xb5, 0xa6, 0xf1, 0xd8, 0x36, 0x8d, 0x74, 0x53, 0xfa, 0xef, 0x71, 0x5f, 0x9a, 0xec, 0xb0,
                                                   0x78, 0xcd, 0x40, 0x2c, 0xbd, 0x54, 0x8c, 0xe,  0x17, 0x9c, 0x4e, 0xd1, 0xe4, 0xc7, 0xe5, 0xb0,
                                                   0x48, 0xe0, 0xa3, 0x9d, 0x31, 0x81, 0x7b, 0x5b, 0x24, 0xf5, 0xd,  0xb5, 0x8b, 0xb3, 0x72, 0xf,
                                                   0xe9, 0x6b, 0xa5, 0x3d, 0xb9, 0x47, 0x84, 0x21, 0x20, 0xa0, 0x68, 0x81, 0x6a, 0xc0, 0x5c, 0x15,
                                                   0x9b, 0xb5, 0x26, 0x6c, 0x63, 0x65, 0x8b, 0x4f, 0x0,  0xc,  0xbf, 0x87, 0xb1, 0x20, 0x9a, 0x22,
                                                   0x5d, 0xef, 0x8e, 0xf1, 0xdc, 0xa9, 0x17, 0xbc, 0xda, 0x79, 0xa1, 0xe4, 0x2a, 0xcd, 0x80, 0x69}},
                {0x80,
                       {0x61, 0x35, 0x31, 0x32, 0x5f, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                    0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                            0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61},
                                           {0x39, 0x69, 0x62, 0xdb, 0x47, 0xf7, 0x49, 0xec, 0x3b, 0x50, 0x42, 0xce, 0x24, 0x52, 0xb6, 0x19,
                                                                                                                                            0x60, 0x7f, 0x27, 0xfd, 0x39, 0x39, 0xec, 0xe2, 0x74, 0x6a, 0x76, 0x14, 0xfb, 0x83, 0xa1, 0xd0,
                                                   0x97, 0xf5, 0x54, 0xdf, 0x39, 0x27, 0xb0, 0x84, 0xe5, 0x5d, 0xe9, 0x2c, 0x78, 0x71, 0x43, 0xd,
                                                   0x6b, 0x95, 0xc2, 0xa1, 0x38, 0x96, 0xd8, 0xa3, 0x3b, 0xc4, 0x85, 0x87, 0xb1, 0xf6, 0x6d, 0x21,
                                                   0xb1, 0x28, 0xa1, 0xa8, 0x24, 0xd,  0x5b, 0xc,  0x26, 0xdf, 0xe7, 0x95, 0xa1, 0xa8, 0x42, 0xa0,
                                                   0x80, 0x7b, 0xb1, 0x48, 0xb7, 0x7c, 0x2e, 0xf8, 0x2e, 0xd4, 0xb6, 0xc9, 0xf7, 0xfc, 0xb7, 0x32,
                                                   0xe7, 0xf9, 0x44, 0x66, 0xc8, 0xb5, 0x1e, 0x52, 0xbf, 0x37, 0x8f, 0xba, 0x4,  0x4a, 0x31, 0xf5,
                                                   0xcb, 0x44, 0x58, 0x3a, 0x89, 0x2f, 0x59, 0x69, 0xdc, 0xd7, 0x3b, 0x3f, 0xa1, 0x28, 0x81, 0x6e}}};

        for (const auto &s: samples) {

            if (std::get<0>(s) == 0x20) {
                using expand_message = typename hashes::detail::expand_message_xmd<128, 32, hash_type>;
                check_expand_message<expand_message>(DST, std::get<1>(s), std::get<2>(s));
            }

            if (std::get<0>(s) == 0x80) {
                using expand_message = hashes::detail::expand_message_xmd<128, 128, hash_type>;
                check_expand_message<expand_message>(DST, std::get<1>(s), std::get<2>(s));
            }
        }
    }

    BOOST_AUTO_TEST_CASE(expand_message_xof_shake128_test) {
        // https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-10#appendix-K.3

        std::string DST_str("QUUX-V01-CS02-with-expander");
        std::vector<std::uint8_t> DST(DST_str.begin(), DST_str.end());

        // {BytesLength, msg, uniform_bytes}
        using samples_type = std::vector<std::tuple<std::size_t, std::vector<std::uint8_t>, std::vector<std::uint8_t>>>;
        samples_type samples{
                {0x20, {},                 {0xec, 0xa3, 0xfe, 0x8f, 0x7f, 0x5f, 0x1d, 0x52, 0xd7, 0xed, 0x36, 0x91, 0xc3, 0x21, 0xad, 0xc7,
                                                                                                                                            0xd2, 0xa0, 0xfe, 0xf1, 0xf8, 0x43, 0xd2, 0x21, 0xf7, 0x00, 0x25, 0x30, 0x07, 0x007, 0x46, 0xde}},
                {0x20, {0x61, 0x62, 0x63}, {0xc7, 0x9b, 0x8e, 0xa0, 0xaf, 0x10, 0xfd, 0x88, 0x71, 0xed, 0xa9,
                                                                                                              0x83, 0x34, 0xea, 0x9d, 0x54, 0xe9, 0xe5, 0x28, 0x2b, 0xe9, 0x75,
                                                                                                                                                                                0x21, 0x67, 0x8f, 0x98, 0x77, 0x18, 0xb1, 0x87,  0xbc, 0x08}},
                {0x20,
                       {0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39},
                                           {0xfb, 0x6f, 0x4a, 0xf2, 0xa8, 0x3f, 0x62, 0x76, 0xe9, 0xd4, 0x17, 0x84, 0xf1, 0xe2, 0x9d, 0xa5,
                                                                                                                                            0xe2, 0x75, 0x66, 0x16, 0x7c, 0x33, 0xe5, 0xcf, 0x26, 0x82, 0xc3, 0x00, 0x96, 0x87,  0x8b, 0x73}},
                {0x20,
                       {0x71, 0x31, 0x32, 0x38, 0x5f, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                              0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                    0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                          0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71},
                                           {0x12, 0x5d, 0x05, 0x85, 0x0d, 0xb9, 0x15, 0xe0, 0x68, 0x3d, 0x17, 0xd0, 0x44, 0xd8, 0x74, 0x77,
                                                                                                                                            0xe6, 0xe7, 0xb3, 0xf7, 0x0a, 0x45, 0x0d, 0xd0, 0x97, 0x76, 0x1e, 0x18, 0xd1, 0xd1,  0xdc, 0xdf}},
                {0x20,
                       {0x61, 0x35, 0x31, 0x32, 0x5f, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                    0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                            0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61},
                                           {0xbe, 0xaf, 0xd0, 0x26, 0xcb, 0x94, 0x2c, 0x86, 0xf6, 0xa2, 0xb3, 0x1b, 0xb8, 0xe6, 0xbf, 0x71,
                                                                                                                                            0x73, 0xfb, 0x1b, 0x0c, 0xaf, 0x3c, 0x21, 0xea, 0x4b, 0x3b, 0x9d, 0x05, 0xd9, 0x04,  0xfd, 0x23}},
                {0x80, {},                 {0x15, 0x73, 0x3b, 0x3f, 0xb2, 0x2f, 0xac, 0x0e, 0x09, 0x02, 0xc2, 0x20, 0xae, 0xea, 0x48, 0xe5,
                                                                                                                                            0xe4, 0x7d, 0x39, 0xf3, 0x6c, 0x2c, 0xc0, 0x3e, 0xac, 0x34, 0x36, 0x7c, 0x48, 0xf2,  0xa3, 0xeb,
                                                   0xbc, 0xb3, 0xba, 0xa8, 0xa0, 0xcf, 0x17, 0xab, 0x12, 0xff, 0xf4, 0xde, 0xfc, 0x7c, 0xe2, 0x2a,
                                                   0xed, 0x47, 0x18, 0x8b, 0x6c, 0x16, 0x3e, 0x82, 0x87, 0x41, 0x47, 0x3b, 0xd8, 0x9c, 0xc6, 0x46,
                                                   0xa0, 0x82, 0xcb, 0x68, 0xb8, 0xe8, 0x35, 0xb1, 0x37, 0x4e, 0xa9, 0xa6, 0x31, 0x5d, 0x61, 0xdb,
                                                   0x00, 0x43, 0xf4, 0xab, 0xf5, 0x06, 0xc2, 0x63, 0x86, 0xe8, 0x46, 0x68, 0xe0, 0x77, 0xc8, 0x5e,
                                                   0xbd, 0x9d, 0x63, 0x2f, 0x43, 0x90, 0x55, 0x9b, 0x97, 0x9e, 0x70, 0xe9, 0xe7, 0xaf, 0xfb, 0xd0,
                                                   0xac, 0x2a, 0x21, 0x2c, 0x03, 0xb6, 0x98, 0xef, 0xbb, 0xe9, 0x40, 0xf2, 0xd1, 0x64, 0x73, 0x2b}},
                {0x80,
                       {0x61, 0x62, 0x63},
                                           {0x4c, 0xca, 0xfb, 0x6d, 0x95, 0xb9, 0x15, 0x37, 0x79, 0x8d, 0x1f, 0xbb, 0x25, 0xb9, 0xfb, 0xe1,
                                                                                                                                            0xa5, 0xbb, 0xe1, 0x68, 0x3f, 0x43, 0xa4, 0xf6, 0xf0, 0x3e, 0xf5, 0x40, 0xb8, 0x11,  0x23, 0x53,
                                                   0x17, 0xbf, 0xc0, 0xae, 0xfb, 0x21, 0x7f, 0xac, 0xa0, 0x55, 0xe1, 0xb8, 0xf3, 0x2d, 0xfd, 0xe9,
                                                   0xeb, 0x10, 0x2c, 0xdc, 0x02, 0x6e, 0xd2, 0x7c, 0xaa, 0x71, 0x53, 0x0e, 0x36, 0x1b, 0x3a, 0xdb,
                                                   0xb9, 0x2c, 0xcf, 0x68, 0xda, 0x35, 0xae, 0xd8, 0xb9, 0xdc, 0x7e, 0x4e, 0x6b, 0x5d, 0xb0, 0x66,
                                                   0x6c, 0x60, 0x7a, 0x31, 0xdf, 0x05, 0x51, 0x3d, 0xda, 0xf4, 0xc8, 0xee, 0x23, 0xb0, 0xee, 0x7f,
                                                   0x39, 0x5a, 0x6e, 0x8b, 0xe3, 0x2e, 0xb1, 0x3c, 0xa9, 0x7d, 0xa2, 0x89, 0xf2, 0x64, 0x36, 0x16,
                                                   0xac, 0x30, 0xfe, 0x91, 0x04, 0xbb, 0x0d, 0x3a, 0x67, 0xa0, 0xa5, 0x25, 0x83, 0x7c, 0x2d, 0xc6}},
                {0x80,
                       {0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39},
                                           {0xc8, 0xee, 0x0e, 0x12, 0x73, 0x6e, 0xfb, 0xc9, 0xb4, 0x77, 0x81, 0xdb, 0x9d, 0x1e, 0x5d, 0xb9,
                                                                                                                                            0xc8, 0x53, 0x68, 0x43, 0x44, 0xa6, 0x77, 0x6e, 0xb3, 0x62, 0xd7, 0x5b, 0x35, 0x4f,  0x4b, 0x74,
                                                   0xcf, 0x60, 0xba, 0x13, 0x73, 0xdc, 0x2e, 0x22, 0xc6, 0x8e, 0xfb, 0x76, 0xa0, 0x22, 0xed, 0x53,
                                                   0x91, 0xf6, 0x7c, 0x77, 0x99, 0x08, 0x02, 0x01, 0x8c, 0x8c, 0xdc, 0x7a, 0xf6, 0xd0, 0x0c, 0x86,
                                                   0xb6, 0x6a, 0x3b, 0x3c, 0xca, 0xd3, 0xf1, 0x8d, 0x90, 0xf4, 0x43, 0x7a, 0x16, 0x51, 0x86, 0xf6,
                                                   0x60, 0x1c, 0xf0, 0xbb, 0x28, 0x1e, 0xa5, 0xd8, 0x0d, 0x1d, 0xe2, 0x0f, 0xe2, 0x2b, 0xb2, 0xe2,
                                                   0xd8, 0xac, 0xab, 0x0c, 0x04, 0x3e, 0x76, 0xe3, 0xa0, 0xf3, 0x4e, 0x0a, 0x1e, 0x66, 0xc9, 0xad,
                                                   0xe4, 0xfe, 0xf9, 0xef, 0x3b, 0x43, 0x11, 0x30, 0xad, 0x6f, 0x23, 0x2b, 0xab, 0xe9, 0xfe, 0x68}},
                {0x80,
                       {0x71, 0x31, 0x32, 0x38, 0x5f, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                              0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                    0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                          0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71, 0x71},
                                           {0x3e, 0xeb, 0xe6, 0x72, 0x1b, 0x2e, 0xc7, 0x46, 0x62, 0x98, 0x56, 0xdc, 0x2d, 0xd3, 0xf0, 0x3a,
                                                                                                                                            0x83, 0x0d, 0xab, 0xfe, 0xfd, 0x7e, 0x2d, 0x1e, 0x72, 0xaa, 0xf2, 0x12, 0x7d, 0x6a,  0xd1, 0x7c,
                                                   0x98, 0x8b, 0x57, 0x62, 0xf3, 0x2e, 0x6e, 0xdf, 0x61, 0x97, 0x23, 0x78, 0xa4, 0x10, 0x6d, 0xc4,
                                                   0xb6, 0x3f, 0xa1, 0x08, 0xad, 0x03, 0xb7, 0x93, 0xee, 0xdf, 0x45, 0x88, 0xf3, 0x4c, 0x4d, 0xf2,
                                                   0xa9, 0x5b, 0x30, 0x99, 0x5a, 0x46, 0x4c, 0xb3, 0xee, 0x31, 0xd6, 0xdc, 0xa3, 0x0a, 0xdb, 0xfc,
                                                   0x90, 0xff, 0xdf, 0x54, 0x14, 0xd7, 0x89, 0x30, 0x82, 0xc5, 0x5b, 0x26, 0x9d, 0x9e, 0xc9, 0xcd,
                                                   0x6d, 0x2a, 0x71, 0x5b, 0x9c, 0x4f, 0xad, 0x4e, 0xb7, 0x0e, 0xd5, 0x6f, 0x87, 0x8b, 0x55, 0xa1,
                                                   0x7b, 0x59, 0x94, 0xef, 0x0d, 0xe5, 0xb3, 0x38, 0x67, 0x5a, 0xad, 0x35, 0x35, 0x41, 0x95, 0xcd}},
                {0x80,
                       {0x61, 0x35, 0x31, 0x32, 0x5f, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                    0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                            0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61,
                               0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61, 0x61},
                                           {0x85, 0x8c, 0xb4, 0xa6, 0xa5, 0x66, 0x8a, 0x97, 0xd0, 0xf7, 0x03, 0x9b, 0x5d, 0x6d, 0x57, 0x4d,
                                                                                                                                            0xde, 0x18, 0xdd, 0x23, 0x23, 0xcf, 0x6b, 0x20, 0x39, 0x45, 0xc6, 0x6d, 0xf8, 0x64,  0x77, 0xd1,
                                                   0xf7, 0x47, 0xb4, 0x64, 0x01, 0x90, 0x3b, 0x3f, 0xa6, 0x6d, 0x12, 0x76, 0x10, 0x8e, 0xa7, 0x18,
                                                   0x7b, 0x44, 0x11, 0xb7, 0x49, 0x9a, 0xcf, 0x46, 0x00, 0x08, 0x0c, 0xe3, 0x4f, 0xf6, 0xd2, 0x15,
                                                   0x55, 0xc2, 0xaf, 0x16, 0xf0, 0x91, 0xad, 0xf8, 0xb2, 0x85, 0xc8, 0x43, 0x9f, 0x2e, 0x47, 0xfa,
                                                   0x05, 0x53, 0xc3, 0xa6, 0xef, 0x5a, 0x42, 0x27, 0xa1, 0x3f, 0x34, 0x40, 0x62, 0x41, 0xb7, 0xd7,
                                                   0xfd, 0x88, 0x53, 0xa0, 0x80, 0xba, 0xd2, 0x5e, 0xc4, 0x80, 0x4c, 0xdf, 0xe4, 0xfd, 0xa5, 0x00,
                                                   0xe1, 0xc8, 0x72, 0xe7, 0x1b, 0x8c, 0x61, 0xa8, 0xe1, 0x60, 0x69, 0x18, 0x94, 0xb9, 0x60, 0x58}}};

        for (const auto &s: samples) {

            if (std::get<0>(s) == 0x20) {
                using expand_message = typename hashes::detail::expand_message_xof<128, 32, hashes::shake<128, 256>>;
                check_expand_message<expand_message>(DST, std::get<1>(s), std::get<2>(s));
            }

            if (std::get<0>(s) == 0x80) {
                using expand_message = hashes::detail::expand_message_xof<128, 128, hashes::shake<128, 1024>>;
                check_expand_message<expand_message>(DST, std::get<1>(s), std::get<2>(s));
            }
        }
    }

    BOOST_AUTO_TEST_CASE(hash_to_field_bls12_381_g1_h2c_sha256_test) {
        // https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-10#appendix-J.9.1
        using curve_type = curves::bls12_381;
        using field_type = typename curve_type::base_field_type;
        using field_value_type = typename field_type::value_type;
        using integral_type = typename field_type::integral_type;
        using hash_type = hashes::h2f<field_type>;

        using samples_type = std::vector<std::tuple<std::string, std::array<field_value_type, 2>>>;
        samples_type samples = {
                {"",
                        {field_value_type(
                                integral_type(
                                        "1790030616568561980207134218344899338736900885118493183248255875682123737756800"
                                        "213955590674957414534085508415116879")),
                                field_value_type(
                                        integral_type(
                                                "2474702583317621523708233292803940741700450584532633563728739973751669085848991"
                                                "00434893060702108665825589810322121"))}},
                {"abc",
                        {field_value_type(
                                integral_type(
                                        "2088728490498894818688784437928579501848367107744050576780266498473771518428420"
                                        "173373487118890161663886009635645777")),
                                field_value_type(
                                        integral_type(
                                                "3213892493831086209316960640873433141017158792584421675273329354360198845384332"
                                                "7878077294514665889481436558332217"))}},
                {"abcdef0123456789",
                        {field_value_type(
                                integral_type(
                                        "9505970308164648217789710156734861296414103440788614747505275085378045493860586"
                                        "12983484048401731236595379325781716")),
                                field_value_type(
                                        integral_type(
                                                "1979385000937648348925653198641340374887185657649818450486460034420643425685140"
                                                "133042050299078521896600910613745210"))}},
                {"q128_"
                 "qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq"
                 "qqqqqqqqqqqqqqqqqqq",
                        {field_value_type(
                                integral_type(
                                        "1565983848840546529547071507571383550794102107851138573768250148104411885455485"
                                        "95465313883035731540725116276838022")),
                                field_value_type(
                                        integral_type(
                                                "1709027689043323463259398100486189187238532958310276339146988040422594808842792"
                                                "053521671901476006506290292962489454"))}},
                {"a512_"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                        {field_value_type(
                                integral_type(
                                        "1625704516324785166230868561544190006281306318060308039760768255839116494270087"
                                        "378351796462565313509233883467016390")),
                                field_value_type(
                                        integral_type(
                                                "8973476190440398924261230730510508241136153370908604317306021021786458550458325"
                                                "65883684732229117125155988066429111"))}}
                // {"", {field_value_type(integral_type("")), field_value_type(integral_type(""))}}
        };

        for (auto &s: samples) {
            check_hash_to_field_ro<hash_type>(std::get<0>(s), std::get<1>(s));
        }
    }

    BOOST_AUTO_TEST_CASE(hash_to_field_bls12_381_g2_h2c_sha256_test) {
        // https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-10#appendix-J.10.1
        using curve_type = curves::bls12_381;
        using group_type = typename curve_type::g2_type<>;
        using field_type = typename group_type::field_type;
        using field_value_type = typename field_type::value_type;
        using integral_type = typename field_type::integral_type;
        using hash_type = hashes::h2f<field_type>;

        using samples_type = std::vector<std::tuple<std::string, std::array<field_value_type, 2>>>;
        samples_type samples = {
                {"",
                        {field_value_type(
                                integral_type(
                                        "5938684483100054485611722523870295163604099457864574398759743150316400213898356"
                                        "49561235021338510064922970633805048"),
                                integral_type(
                                        "8673753094890675127974598608873659518770540387638184480573261903027016498888499"
                                        "97836339069389536967202878289851290")),
                                field_value_type(
                                        integral_type(
                                                "4578897045199488434740260225626419694433157155954591591128744980829534319713238"
                                                "09145630315884223143822925947137684"),
                                        integral_type(
                                                "3132697209754082586339430915081913810572071485832539443682634025529375380328136"
                                                "128542015469873094481703191673087029"))}},
                {"abc",
                        {field_value_type(
                                integral_type(
                                        "3381151350286428005095780827831774583653641216459357823974407145557165174365389"
                                        "989442078766443621078367363453769585"),
                                integral_type(
                                        "2741746953704442638534180707453397316404679193551841082537168795196953970699630"
                                        "34977795744692362177212201505728989")),
                                field_value_type(
                                        integral_type(
                                                "3761918608077574755256083960277010506684793456226386707192711779006489497410866"
                                                "269311252402421709839991039401264868"),
                                        integral_type(
                                                "1342131492846344403298252211066711749849099599627623100864413228392326132610002"
                                                "371925674088601653350525231531947366"))}},
                {"abcdef0123456789",
                        {field_value_type(
                                integral_type(
                                        "4736756665618245326244300857865191860224326611904114213007749037224882541543738"
                                        "95989233527517731907580580706354657"),
                                integral_type(
                                        "9520540557415691916362510867127307131683791692159959526594533787977337613244945"
                                        "87640793580119096894387397115436943")),
                                field_value_type(
                                        integral_type(
                                                "3574336717567028224405133950386477048284620456829914449302272757384276784667241"
                                                "972055005113408837488328262928878231"),
                                        integral_type(
                                                "2365602345707797244937763470382803726723577073883311775921418854730692345417958"
                                                "26215789679703490403053611203549557"))}},
                {"q128_"
                 "qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq"
                 "qqqqqqqqqqqqqqqqqqq",
                        {field_value_type(
                                integral_type(
                                        "3608131929677217503005188861991391449980483387988256142839334881042292007389752"
                                        "19634991234987258997592645316502099"),
                                integral_type(
                                        "5009904385311786096049606538586133897401985944525938049981086917265658825017777"
                                        "15476408413735192405455364595747963")),
                                field_value_type(
                                        integral_type(
                                                "1414201600433038156752401103621159164529164806638579329495300394501933973057103"
                                                "319123042671630779248244072674138005"),
                                        integral_type(
                                                "2580989994757912640015815541704972436791025324967858519264081257257405036397177"
                                                "981572950833626047365407639272235247"))}},
                {"a512_"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                        {field_value_type(
                                integral_type(
                                        "3854656460966118202185202795415969034444473478700041637108073179423449626727403"
                                        "291696647010132509133525205314259253"),
                                integral_type(
                                        "2873494353363126311409085895530381085174075451844000378947122252646711114869905"
                                        "923958066527312260237781765269081913")),
                                field_value_type(
                                        integral_type(
                                                "2218682278840147973132952196327912255143646871258838127959845658885016361690895"
                                                "544274403462155614933990666846598837"),
                                        integral_type(
                                                "2692054640040186323570630735219910885988179020142391687801930252786130591827501"
                                                "100656577702624784849458500251540952"))}},

                // {"",
                //  {field_value_type(integral_type(""),
                //                    integral_type("")),
                //   field_value_type(integral_type(""),
                //                    integral_type(""))}},
        };

        for (auto &s: samples) {
            check_hash_to_field_ro<hash_type>(std::get<0>(s), std::get<1>(s));
        }
    }

    BOOST_AUTO_TEST_CASE(hash_to_field_bls12_381_g1_h2f_shake128_test) {
        // field elements were generated using the SageMath code from https://github.com/cfrg/draft-irtf-cfrg-hash-to-curve
        using curve_type = curves::bls12_381;
        using field_type = typename curve_type::base_field_type;
        using field_value_type = typename field_type::value_type;
        using integral_type = typename field_type::integral_type;
        using hash_type = hashes::h2f<field_type,
                hashes::shake<128, 1024>,
                hashes::h2f_default_params<field_type,
                        hashes::shake<128, 1024>,
                        128,
                        hashes::uniformity_count_t::uniform_count,
                        hashes::expand_msg_variant_t::rfc_xof>>;

        using samples_type = std::vector<std::tuple<std::string, std::array<field_value_type, 2>>>;
        samples_type samples = {
                {"",
                        {field_value_type(integral_type(
                                "34026336386083886614678315776001410990195491156146619751680086308817537707641"
                                "2511339581237807169328025849613607738")),
                                field_value_type(integral_type(
                                        "35487283142688730820577587057340419437479971382706358249913182139983177321994"
                                        "22366369249146063877280813873779598614"))}},
                {"abc",
                        {field_value_type(integral_type(
                                "10299742258189792071225127476569579055000139249317173341006693762374554218361"
                                "99735720081420774878644714920465448063")),
                                field_value_type(integral_type(
                                        "26449079556749456547763240104029704617197602953243811123341487416058796720985"
                                        "30793364734851520804076509786066300688"))}},
                {"abcdef0123456789",
                        {field_value_type(integral_type(
                                "25719599515144616098862759599781083245320818219725699327118811068540459040324"
                                "01497739611428750762925771009739672760")),
                                field_value_type(integral_type(
                                        "48152350957069938418503459090703763601658613290098485815305620880423150293107"
                                        "3834859207979148213406515909166856065"))}},
                {"q128_"
                 "qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq"
                 "qqqqqqqqqqqqqqqqqqq",
                        {field_value_type(integral_type(
                                "16746159338985009492917568934031554970366258991142848241908114029620654980095"
                                "63818288046948236719654148630540204822")),
                                field_value_type(integral_type(
                                        "25816575372444607942054642448052352758295864608466091303597503924099696329006"
                                        "72887772936290792382770151227422762659"))}},
                {"a512_"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                        {field_value_type(integral_type(
                                "99872323850197704338735859571501164993467800869061716944285932457586706162378"
                                "994013296878816389726622832817306892")),
                                field_value_type(integral_type(
                                        "75080915218496132457991278757568526706694821515781542981734814458320346931034"
                                        "5473734242537722102188109771665077578"))}}
                // {"", {field_value_type(integral_type("")), field_value_type(integral_type(""))}}
        };

        for (auto &s: samples) {
            check_hash_to_field_ro<hash_type>(std::get<0>(s), std::get<1>(s));
        }
    }

    BOOST_AUTO_TEST_CASE(hash_to_field_bls12_381_g2_h2f_shake128_test) {
        // https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-10#appendix-J.10.1
        using curve_type = curves::bls12_381;
        using group_type = typename curve_type::g2_type<>;
        using field_type = typename group_type::field_type;
        using field_value_type = typename field_type::value_type;
        using integral_type = typename field_type::integral_type;
        using hash_type = hashes::h2f<field_type,
                hashes::shake<128, 2048>,
                hashes::h2f_default_params<field_type,
                        hashes::shake<128, 2048>,
                        128,
                        hashes::uniformity_count_t::uniform_count,
                        hashes::expand_msg_variant_t::rfc_xof>>;

        using samples_type = std::vector<std::tuple<std::string, std::array<field_value_type, 2>>>;
        samples_type samples = {
                {"",
                        {field_value_type(integral_type(
                                                  "21761119417140609572957506711724462644916292942987174604015733000612083054138"
                                                  "77868370459391860290442637242685798287"),
                                          integral_type(
                                                  "21727803743924291497680999138943700735898742815755951136165147607237967756386"
                                                  "41300727609117779478067420156555937182")),
                                field_value_type(integral_type(
                                                         "19802774580858068514278700330311615707558747455624782896791901258745532214031"
                                                         "68804685581811930976162756736377819200"),
                                                 integral_type(
                                                         "30577482064538241631055790494561748316708062368356777820792422645161630521893"
                                                         "90202713360090520912664585635089365842"))}},
                {"abc",
                        {field_value_type(integral_type(
                                                  "28095779702340744474188222080117336804030388680055502193888265667402470600038"
                                                  "53667807528408681210194383505927321734"),
                                          integral_type(
                                                  "13480905020538778032639559738025205303070065891728990555656586377754061695581"
                                                  "74504663705980198627002645830947565535")),
                                field_value_type(integral_type(
                                                         "36128837908863750412786113353103990083658782138245371922196435006135331142272"
                                                         "55727791939831017201995232349524111099"),
                                                 integral_type(
                                                         "37031718815542562878057377216607493299341531037255524917868423138783142182529"
                                                         "72527052232848741612107413733184042704"))}},
                {"abcdef0123456789",
                        {field_value_type(integral_type(
                                                  "34179051443373825021627221887887003686762026645123594023800727335353607515140"
                                                  "90249139105898909846372034808498893681"),
                                          integral_type(
                                                  "18602926449762697029401068753653039886784602628583881790458942559844506726399"
                                                  "97323424814684516295366461881861836534")),
                                field_value_type(integral_type(
                                                         "13725831502352924801704573350984719213058089778288294729248198070849411653182"
                                                         "66758515896635824313576096391677734981"),
                                                 integral_type(
                                                         "36908193745246653651064804888056323355497509997526047701009643674942412752742"
                                                         "19038923167463498905572382098560535117"))}},
                {"q128_"
                 "qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq"
                 "qqqqqqqqqqqqqqqqqqq",
                        {field_value_type(integral_type(
                                                  "39599762187870594643873274175209553764202031577019663999380556114139666401464"
                                                  "89507356572246555322990727985468040245"),
                                          integral_type(
                                                  "25676933779138698592679819807400247892045745812751340998780556381675609475573"
                                                  "09341944273444974680463976608552489852")),
                                field_value_type(integral_type(
                                                         "61454319429737261771070011869457935134248576585039662467468431549997609359625"
                                                         "6444120666417881167834118723680453590"),
                                                 integral_type(
                                                         "23955509127144125773005085188669236528208316300227858678645065700992358659384"
                                                         "54674609909635896529114812122026503129"))}},
                {"a512_"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                        {field_value_type(integral_type(
                                                  "42031133646911283970924437529711587195926963466037233984948270009038220365544"
                                                  "5945666063808629865323003385081140252"),
                                          integral_type(
                                                  "32382326666448732822661723450194270503118452492444836120698898827195542329051"
                                                  "10784475616069209850193778954228519664")),
                                field_value_type(integral_type(
                                                         "12731197979403258068394380733350947413090025310423363781063761800501631636886"
                                                         "27388411388878626108888409981407889320"),
                                                 integral_type(
                                                         "15029123207785810142032112101030614655730726369761948395001369689467321863650"
                                                         "57656474980051133773877281944958080905"))}},

                // {"",
                //  {field_value_type(integral_type(""),
                //                    integral_type("")),
                //   field_value_type(integral_type(""),
                //                    integral_type(""))}},
        };

        for (auto &s: samples) {
            check_hash_to_field_ro<hash_type>(std::get<0>(s), std::get<1>(s));
        }
    }

    BOOST_AUTO_TEST_CASE(hash_to_field_bls12_381_g1_h2f_shake256_test) {
        // field elements were generated using the SageMath code from https://github.com/cfrg/draft-irtf-cfrg-hash-to-curve
        using curve_type = curves::bls12_381;
        using field_type = typename curve_type::base_field_type;
        using field_value_type = typename field_type::value_type;
        using integral_type = typename field_type::integral_type;
        using hash_type = hashes::h2f<field_type,
                hashes::shake<256, 1024>,
                hashes::h2f_default_params<field_type,
                        hashes::shake<256, 1024>,
                        128,
                        hashes::uniformity_count_t::uniform_count,
                        hashes::expand_msg_variant_t::rfc_xof>>;

        using samples_type = std::vector<std::tuple<std::string, std::array<field_value_type, 2>>>;
        samples_type samples = {
                {"",
                        {field_value_type(integral_type(
                                "35649413618876074642626806365137625720526779547687130805164279468893543348816"
                                "81149428744233836573075940920773556605")),
                                field_value_type(integral_type(
                                        "32191598794272535593202043296699690894078346767469611555064692257309991804181"
                                        "2868917504156030442189555841682260203"))}},
                {"abc",
                        {field_value_type(integral_type(
                                "36297070445601380522183773418551096131819097989283494533273211638846992997294"
                                "86359633941000674339086857170774012891")),
                                field_value_type(integral_type(
                                        "21336734814454478909919797987655643716116374700180214090544223636910007807464"
                                        "09314360354255496017833770999585612365"))}},
                {"abcdef0123456789",
                        {field_value_type(integral_type(
                                "89562627513935962103233763683462050306922022114977595123138177164512223976408"
                                "2253394057881856751461508824060407068")),
                                field_value_type(integral_type(
                                        "26390405268558144935650411401535229593735188494011206900329278546785106492536"
                                        "24259827299638645593197965314405990223"))}},
                {"q128_"
                 "qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq"
                 "qqqqqqqqqqqqqqqqqqq",
                        {field_value_type(integral_type(
                                "22916397665345561448194985375368991003209854878197341043956812081252109413817"
                                "62551072041089523763472685246541885147")),
                                field_value_type(integral_type(
                                        "29915744450377534368569430780604829283215920795851547158261172894680433208614"
                                        "31194687978643395782359715681390717481"))}},
                {"a512_"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                        {field_value_type(integral_type(
                                "24547068103260723032692377284190124714469178784622226928569646544131198218365"
                                "12158594441511421199579041183882320701")),
                                field_value_type(integral_type(
                                        "32777440595410240006428031495965340366546892470064030857364419905126509912456"
                                        "78958004541570391981799856710821139287"))}}
                // {"", {field_value_type(integral_type("")), field_value_type(integral_type(""))}}
        };

        for (auto &s: samples) {
            check_hash_to_field_ro<hash_type>(std::get<0>(s), std::get<1>(s));
        }
    }

    BOOST_AUTO_TEST_CASE(hash_to_field_bls12_381_g2_h2f_shake256_test) {
        // https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-10#appendix-J.10.1
        using curve_type = curves::bls12_381;
        using group_type = typename curve_type::g2_type<>;
        using field_type = typename group_type::field_type;
        using field_value_type = typename field_type::value_type;
        using integral_type = typename field_type::integral_type;
        using hash_type = hashes::h2f<field_type,
                hashes::shake<256, 2048>,
                hashes::h2f_default_params<field_type,
                        hashes::shake<256, 2048>,
                        128,
                        hashes::uniformity_count_t::uniform_count,
                        hashes::expand_msg_variant_t::rfc_xof>>;

        using samples_type = std::vector<std::tuple<std::string, std::array<field_value_type, 2>>>;
        samples_type samples = {
                {"",
                        {field_value_type(integral_type(
                                                  "21963274071619108128097264722317071301606486779929326898554083672936128235770"
                                                  "4338072415018530247540214055249982132"),
                                          integral_type(
                                                  "10590351498334984848373840030203446950121731278087575746096550910348395558555"
                                                  "1845368665743236877541939312930239681")),
                                field_value_type(integral_type(
                                                         "28091471905464214407695546469896700624773788481198476844345274031808508333074"
                                                         "97632818619741493854405328865108683275"),
                                                 integral_type(
                                                         "30726366767257068092443502510012136173512294464519767218243017912162500672713"
                                                         "34346694333455334523250618559647124233"))}},
                {"abc",
                        {field_value_type(integral_type(
                                                  "23225806639179019890338382797325715438756600031652555722973808596165865668617"
                                                  "34165851712176216313509198818511530980"),
                                          integral_type(
                                                  "32557071023082220791916058169457485364134976942295371950651552550860380978496"
                                                  "48899273319477344766487385190825894003")),
                                field_value_type(integral_type(
                                                         "32194294443752251899867889497549719798956652948195916403399424005829214370036"
                                                         "91862946903896256965180933287123015934"),
                                                 integral_type(
                                                         "28319630043143611222959994456903801307817254469122958785870712904900231118040"
                                                         "53145378372949657374193108898073660733"))}},
                {"abcdef0123456789",
                        {field_value_type(integral_type(
                                                  "17902175988373130734012933586467293569612298137867080532295054921824759111317"
                                                  "92567425530299937389989900138225324599"),
                                          integral_type(
                                                  "12969006743881285697018500257524475024956462836065688657073044500573635075580"
                                                  "27645375674431249396430197967082812936")),
                                field_value_type(integral_type(
                                                         "34406529846799115946535633828620145267368064157958635700227172206263842971745"
                                                         "02640487083363023891622316284584216690"),
                                                 integral_type(
                                                         "28213835204899604742609273129426581197007138205118023703504592620948962976872"
                                                         "75979972257250936391927795870833532960"))}},
                {"q128_"
                 "qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq"
                 "qqqqqqqqqqqqqqqqqqq",
                        {field_value_type(integral_type(
                                                  "10092932943035553727564168499381428496298566886087918016632369681584459692382"
                                                  "01675258277772727721484256251128244060"),
                                          integral_type(
                                                  "17168004723799863694173561785055094879931106899957934052491965588500285242621"
                                                  "24538058198857628689491715908814880409")),
                                field_value_type(integral_type(
                                                         "25719222627124083185067835874983100017932189620376996637741387888990217611488"
                                                         "60467345466576023662564107715706227472"),
                                                 integral_type(
                                                         "36331609886282561421947333390942915867716324305897456844342657503519500457553"
                                                         "83718695873025565472006812136924244469"))}},
                {"a512_"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                        {field_value_type(integral_type(
                                                  "35495382464173722991366714917826997437591972701057816539555858479132737998763"
                                                  "6201901377834680498342583094405742933"),
                                          integral_type(
                                                  "16036052280914287202478767463603457793835943149293145136061129905834330994782"
                                                  "51445849500798421048684965287365764494")),
                                field_value_type(integral_type(
                                                         "19494580623629045697224611727282567482926073920822674841457409389744912206228"
                                                         "80992714498950980540123997447376922328"),
                                                 integral_type(
                                                         "98171940250611754694458305588405551184270341958151133249793446650353577775558"
                                                         "2925716311702110312401790485459582370"))}},

                // {"",
                //  {field_value_type(integral_type(""),
                //                    integral_type("")),
                //   field_value_type(integral_type(""),
                //                    integral_type(""))}},
        };

        for (auto &s: samples) {
            check_hash_to_field_ro<hash_type>(std::get<0>(s), std::get<1>(s));
        }
    }

BOOST_AUTO_TEST_SUITE_END()
