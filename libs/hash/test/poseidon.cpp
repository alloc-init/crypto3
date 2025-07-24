//---------------------------------------------------------------------------//
// Copyright (c) 2020 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020 Ilias Khairullin <ilias@nil.foundation>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#define BOOST_TEST_MODULE poseidon_test

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <boost/property_tree/ptree.hpp>

#include <nil/crypto3/hash/algorithm/hash.hpp>
#include <nil/crypto3/hash/block_to_field_elements_wrapper.hpp>
#include <nil/crypto3/hash/detail/poseidon/poseidon_policy.hpp>
#include <nil/crypto3/hash/detail/poseidon/poseidon_permutation.hpp>
#include <nil/crypto3/hash/hash_state.hpp>
#include <nil/crypto3/hash/poseidon.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/scalar_field.hpp>
#include <nil/crypto3/algebra/fields/bls12/scalar_field.hpp>
#include <nil/crypto3/algebra/fields/pallas/base_field.hpp>
#include <nil/crypto3/algebra/fields/pallas/scalar_field.hpp>

using namespace nil::crypto3;
using namespace nil::crypto3::accumulators;
using namespace nil::crypto3::algebra;
using namespace nil::crypto3::hashes::detail;


namespace boost {
    namespace test_tools {
        namespace tt_detail {

            // Functions required by boost, to be able to print the compared values, when assertion fails.
            // TODO(martun): it would be better to implement operator<< for each field element type.
            template<typename FieldParams>
            struct print_log_value<typename fields::detail::element_fp<FieldParams>> {
                void operator()(std::ostream &os, typename fields::detail::element_fp<FieldParams> const &e) {
                    os << e.data << std::endl;
                }
            };

            template<typename FieldParams, size_t array_size>
            struct print_log_value<std::array<typename fields::detail::element_fp<FieldParams>, array_size>> {
                void operator()(std::ostream &os,
                                std::array<typename fields::detail::element_fp<FieldParams>, array_size> const &arr) {
                    for (auto &e: arr) {
                        os << e.data << std::endl;
                    }
                }
            };

            //template<template<typename, typename> class P, typename K, typename V>
            //struct print_log_value<P<K, V>> {
            //    void operator()(std::ostream &, P<K, V> const &) {
            //    }
            //};

        } // namespace tt_detail
    }     // namespace test_tools
}         // namespace boost

template<typename field_type>
void test_pasta_poseidon(std::vector<typename field_type::value_type> input,
                         typename field_type::value_type expected_result) {
    using policy = pasta_poseidon_policy<field_type>;
    using hash_t = hashes::original_poseidon<policy>;

    typename policy::digest_type d = hash<hash_t>(input);
    BOOST_CHECK_EQUAL(d, expected_result);

    accumulator_set<hash_t> acc;

    for (auto &val : input) {
        acc(val);
    }
    typename hash_t::digest_type d_acc = extract::hash<hash_t>(acc);
    BOOST_CHECK_EQUAL(d_acc, expected_result);
}

template<typename FieldType, size_t Rate>
void test_poseidon_permutation(typename poseidon_policy<FieldType, 128, Rate>::state_type input,
                               typename poseidon_policy<FieldType, 128, Rate>::state_type expected_result) {
    using policy = poseidon_policy<FieldType, 128, Rate>;

    // This permutes in place.
    poseidon_permutation<policy>::permute(input);
    BOOST_CHECK_EQUAL(input, expected_result);
}

BOOST_AUTO_TEST_SUITE(poseidon_tests)

BOOST_AUTO_TEST_CASE(poseidon_without_padding_test) {
    using field_type = fields::pallas_scalar_field;
    using policy = pasta_poseidon_policy<field_type>;
    using hash_type = hashes::poseidon<policy>;

    std::vector<uint8_t> hash_input1 {0x00, 0x01, 0x02};
    std::vector<uint8_t> hash_input2 {0x01, 0x02};

    /* Explicit non-padding behavior: input bytes converted into field elements
     * without padding, both inputs produce same element:
     * 0x0000000000000000000000000000000000000000000000000000000000000102
     */
    typename policy::digest_type result1 = nil::crypto3::hash<hash_type>(
        nil::crypto3::hashes::conditional_block_to_field_elements_wrapper<typename hash_type::word_type,
                                                                          decltype(hash_input1),
                                                                          true>(hash_input1));

    typename policy::digest_type result2 = nil::crypto3::hash<hash_type>(
        nil::crypto3::hashes::conditional_block_to_field_elements_wrapper<typename hash_type::word_type,
                                                                          decltype(hash_input2),
                                                                          true>(hash_input2));

    /* Results should be equal */
    BOOST_CHECK_EQUAL(result1, result2);
}

// Poseidon permutation test vectors are taken from:
//   https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/208b5a164c6a252b137997694d90931b2bb851c5/code/test_vectors.txt
BOOST_AUTO_TEST_CASE(poseidon_permutation_254_2) {
    test_poseidon_permutation<fields::alt_bn128_scalar_field<254>, 2>(
        {0x0000000000000000000000000000000000000000000000000000000000000000_cppui_modular254,
         0x0000000000000000000000000000000000000000000000000000000000000001_cppui_modular254,
         0x0000000000000000000000000000000000000000000000000000000000000002_cppui_modular254},
        {0x115cc0f5e7d690413df64c6b9662e9cf2a3617f2743245519e19607a4417189a_cppui_modular254,
         0x0fca49b798923ab0239de1c9e7a4a9a2210312b6a2f616d18b5a87f9b628ae29_cppui_modular254,
         0x0e7ae82e40091e63cbd4f16a6d16310b3729d4b6e138fcf54110e2867045a30c_cppui_modular254});
}

BOOST_AUTO_TEST_CASE(poseidon_permutation_254_4) {
    test_poseidon_permutation<fields::alt_bn128_scalar_field<254>, 4>(
        {0x0000000000000000000000000000000000000000000000000000000000000000_cppui_modular254,
         0x0000000000000000000000000000000000000000000000000000000000000001_cppui_modular254,
         0x0000000000000000000000000000000000000000000000000000000000000002_cppui_modular254,
         0x0000000000000000000000000000000000000000000000000000000000000003_cppui_modular254,
         0x0000000000000000000000000000000000000000000000000000000000000004_cppui_modular254},
        {0x299c867db6c1fdd79dcefa40e4510b9837e60ebb1ce0663dbaa525df65250465_cppui_modular254,
         0x1148aaef609aa338b27dafd89bb98862d8bb2b429aceac47d86206154ffe053d_cppui_modular254,
         0x24febb87fed7462e23f6665ff9a0111f4044c38ee1672c1ac6b0637d34f24907_cppui_modular254,
         0x0eb08f6d809668a981c186beaf6110060707059576406b248e5d9cf6e78b3d3e_cppui_modular254,
         0x07748bc6877c9b82c8b98666ee9d0626ec7f5be4205f79ee8528ef1c4a376fc7_cppui_modular254});
}

BOOST_AUTO_TEST_CASE(poseidon_permutation_255_3) {
    test_poseidon_permutation<fields::bls12_scalar_field<381>, 2>(
        {0x0000000000000000000000000000000000000000000000000000000000000000_cppui_modular255,
         0x0000000000000000000000000000000000000000000000000000000000000001_cppui_modular255,
         0x0000000000000000000000000000000000000000000000000000000000000002_cppui_modular255},
        {0x28ce19420fc246a05553ad1e8c98f5c9d67166be2c18e9e4cb4b4e317dd2a78a_cppui_modular255,
         0x51f3e312c95343a896cfd8945ea82ba956c1118ce9b9859b6ea56637b4b1ddc4_cppui_modular255,
         0x3b2b69139b235626a0bfb56c9527ae66a7bf486ad8c11c14d1da0c69bbe0f79a_cppui_modular255});
}

BOOST_AUTO_TEST_CASE(poseidon_permutation_255_4) {
    test_poseidon_permutation<fields::bls12_scalar_field<381>, 4>(
        {0x0000000000000000000000000000000000000000000000000000000000000000_cppui_modular255,
         0x0000000000000000000000000000000000000000000000000000000000000001_cppui_modular255,
         0x0000000000000000000000000000000000000000000000000000000000000002_cppui_modular255,
         0x0000000000000000000000000000000000000000000000000000000000000003_cppui_modular255,
         0x0000000000000000000000000000000000000000000000000000000000000004_cppui_modular255},
        {0x2a918b9c9f9bd7bb509331c81e297b5707f6fc7393dcee1b13901a0b22202e18_cppui_modular255,
         0x65ebf8671739eeb11fb217f2d5c5bf4a0c3f210e3f3cd3b08b5db75675d797f7_cppui_modular255,
         0x2cc176fc26bc70737a696a9dfd1b636ce360ee76926d182390cdb7459cf585ce_cppui_modular255,
         0x4dc4e29d283afd2a491fe6aef122b9a968e74eff05341f3cc23fda1781dcb566_cppui_modular255,
         0x03ff622da276830b9451b88b85e6184fd6ae15c8ab3ee25a5667be8592cce3b1_cppui_modular255});
}

BOOST_AUTO_TEST_CASE(poseidon_accumulator_255_4) {
    using policy = poseidon_policy<fields::bls12_scalar_field<381>, 128, /*Rate=*/4>;
    using hash_type = hashes::poseidon<policy>;
    accumulator_set<hash_type> acc;

    policy::word_type val = 0u;

    acc(val);

    hash_type::digest_type s = extract::hash<hash_type>(acc);

    BOOST_CHECK_EQUAL(s, 0x20CDA7B88718C51A894AE697F804FACD408616B1A7811A55023EA0E6060AA61C_cppui_modular255);
}

BOOST_AUTO_TEST_CASE(poseidon_stream_255_4) {
    // Since we don't have any test vectors for such a custom structure, just make sure
    // it produces something consistent
    using field_type = fields::bls12_scalar_field<381>;
    using policy = poseidon_policy<field_type, 128, /*Rate=*/4>;
    using hash_type = hashes::poseidon<policy>;

    std::vector<typename field_type::value_type> input = {
        0x0_cppui_modular255, 0x0_cppui_modular255, 0x0_cppui_modular255, 0x0_cppui_modular255, 0x0_cppui_modular255};

    typename policy::digest_type d = hash<hash_type>(input);
    BOOST_CHECK_EQUAL(d, 0x44753e7f86d80790e762345ff8cb156be18eb0318f8846641193f815fbd64038_cppui_modular255);

    input = {0x2a918b9c9f9bd7bb509331c81e297b5707f6fc7393dcee1b13901a0b22202e18_cppui_modular255,
             0x65ebf8671739eeb11fb217f2d5c5bf4a0c3f210e3f3cd3b08b5db75675d797f7_cppui_modular255,
             0x2cc176fc26bc70737a696a9dfd1b636ce360ee76926d182390cdb7459cf585ce_cppui_modular255,
             0x4dc4e29d283afd2a491fe6aef122b9a968e74eff05341f3cc23fda1781dcb566_cppui_modular255,
             0x03ff622da276830b9451b88b85e6184fd6ae15c8ab3ee25a5667be8592cce3b1_cppui_modular255};

    d = hash<hash_type>(input);
    BOOST_CHECK_EQUAL(d, 0x44bff12d3a4713b18bd79c17eaabf8e69e29ce45ca48d7afb702baa1c37f3695_cppui_modular255);
}

BOOST_AUTO_TEST_CASE(poseidon_wrapped_255_4) {
    // Make sure nil_block_poseidon converts non-field input to field elements as we expect
    using field_type = fields::bls12_scalar_field<381>;
    using policy = poseidon_policy<field_type, 128, /*Rate=*/4>;
    using hash_type = hashes::poseidon<policy>;

    std::vector<std::uint8_t> uint8_input = {
        0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF,
        0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD,
        0xEF,    // 256 bits up to this place, the last value should be moved to
                 // the next field element.

    };

    std::vector<typename field_type::value_type> field_input = {
        0x000123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCD_cppui_modular255,
        0x00000000000000000000000000000000000000000000000000000000000000EF_cppui_modular255,
    };

    typename policy::digest_type d_uint8 = hash<hash_type>(
        hashes::conditional_block_to_field_elements_wrapper<hash_type::word_type, decltype(uint8_input)>(uint8_input));
    typename policy::digest_type d_field = hash<hash_type>(field_input);
    BOOST_CHECK_EQUAL(d_uint8, d_field);
}

// This test can be useful for constants generation in the future.
// BOOST_AUTO_TEST_CASE(poseidon_generate_pallas_constants) {
//
//    typedef poseidon_policy<nil::crypto3::algebra::fields::pallas_base_field, 128, 2> PolicyType;
//    typedef poseidon_constants_generator<PolicyType> generator_type;
//    generator_type generator;
//    auto constants = generator.generate_constants();
//}

BOOST_AUTO_TEST_SUITE_END()
