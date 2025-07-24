//---------------------------------------------------------------------------//
// Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020 Nikita Kaskov <nbering@nil.foundation>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#define BOOST_TEST_MODULE des_cipher_test

#include <iostream>
#include <cstdint>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include <nil/crypto3/block/algorithm/encrypt.hpp>
#include <nil/crypto3/block/algorithm/decrypt.hpp>

#include <nil/crypto3/block/des.hpp>

#include <nil/crypto3/block/cipher_value.hpp>

using namespace nil::crypto3;
using namespace nil::crypto3::block;
using namespace nil::crypto3::detail;

BOOST_AUTO_TEST_SUITE(des_stream_processor_filedriven_test_suite)

BOOST_AUTO_TEST_CASE(des_1) {

    std::vector<char> input = {'\x05', '\x9b', '\x5e', '\x08', '\x51', '\xcf', '\x14', '\x3a'};
    std::vector<char> key = {'\x01', '\x13', '\xb9', '\x70', '\xfd', '\x34', '\xf2', '\xce'};
    
    auto res = encrypt<block::des>(input, key);
    std::string out = res;
    
    BOOST_CHECK_EQUAL(out, "86a560f10ec6d85b");
}
/*
BOOST_DATA_TEST_CASE(des_ecb, string_data("ecb_fixed_key"), triples) {

    byte_string const p(triples.first);

    BOOST_FOREACH(boost::property_tree::ptree::value_type pair, triples.second) {
        byte_string const K(pair.first);

        std::string out = encrypt<block::des>(p, K);

        BOOST_CHECK_EQUAL(out, pair.second.data());
    }
}*/

BOOST_AUTO_TEST_SUITE_END()
