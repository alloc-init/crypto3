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
#include <vector>
#include <type_traits>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/hash/algorithm/hash.hpp>
#include <nil/crypto3/hash/block_to_field_elements_wrapper.hpp>
#include <nil/crypto3/hash/detail/poseidon_common/poseidon_sponge.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_optimized_permutation.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_permutation.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_policy.hpp>
#include <nil/crypto3/hash/detail/poseidon2/poseidon2_permutation.hpp>
#include <nil/crypto3/hash/detail/poseidon2/poseidon2_policy.hpp>
#include <nil/crypto3/hash/poseidon.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/scalar_field.hpp>
#include <nil/crypto3/algebra/fields/bls12/scalar_field.hpp>

using namespace nil::crypto3;
using namespace nil::crypto3::algebra;
using namespace nil::crypto3::hashes::detail;

template<typename PolicyType,
         typename PermutationType,
         poseidon_sponge_absorb_mode AbsorbMode,
         poseidon_sponge_padding_mode PaddingMode>
class poseidon_sponge_test_hash {
public:
    using policy_type = PolicyType;
    using permutation_type = PermutationType;

    using word_type = typename policy_type::word_type;
    constexpr static const std::size_t word_bits = policy_type::word_bits;

    constexpr static const std::size_t block_words = policy_type::block_words;
    using block_type = typename policy_type::block_type;

    constexpr static const std::size_t digest_bits = policy_type::digest_bits;
    using digest_type = typename policy_type::digest_type;

    struct construction {
        struct params_type { };

        using type = poseidon_sponge_construction<policy_type, permutation_type, AbsorbMode, PaddingMode>;
    };

    constexpr static nil::crypto3::hashes::detail::stream_processor_type stream_processor =
        nil::crypto3::hashes::detail::stream_processor_type::raw;
    using accumulator_tag = accumulators::tag::algebraic_hash<
        poseidon_sponge_test_hash<policy_type, permutation_type, AbsorbMode, PaddingMode>>;
};

template<typename PolicyType, typename PermutationType>
using poseidon_pad10_test_hash = poseidon_sponge_test_hash<PolicyType,
                                                           PermutationType,
                                                           poseidon_sponge_absorb_mode::overwrite,
                                                           poseidon_sponge_padding_mode::pad10>;

template<typename PolicyType, typename PermutationType>
using poseidon_padding_free_test_hash = poseidon_sponge_test_hash<PolicyType,
                                                                  PermutationType,
                                                                  poseidon_sponge_absorb_mode::overwrite,
                                                                  poseidon_sponge_padding_mode::padding_free>;

template<typename PolicyType, typename PermutationType>
using poseidon_add_pad10_test_hash = poseidon_sponge_test_hash<PolicyType,
                                                               PermutationType,
                                                               poseidon_sponge_absorb_mode::add,
                                                               poseidon_sponge_padding_mode::pad10>;

template<typename PolicyType, typename PermutationType>
using poseidon_add_padding_free_test_hash = poseidon_sponge_test_hash<PolicyType,
                                                                      PermutationType,
                                                                      poseidon_sponge_absorb_mode::add,
                                                                      poseidon_sponge_padding_mode::padding_free>;

namespace boost {
    namespace test_tools {
        namespace tt_detail {

            // Functions required by boost, to be able to print the compared values, when assertion fails.
            // TODO(martun): it would be better to implement operator<< for each field element type.
            template<typename FieldParams>
            struct print_log_value<typename fields::detail::element_fp<FieldParams>> {
                void operator()(std::ostream &os, typename fields::detail::element_fp<FieldParams> const &e) {
                    os << e << std::endl;
                }
            };

            template<typename FieldParams, size_t array_size>
            struct print_log_value<std::array<typename fields::detail::element_fp<FieldParams>, array_size>> {
                void operator()(std::ostream &os,
                                std::array<typename fields::detail::element_fp<FieldParams>, array_size> const &arr) {
                    for (auto &e : arr) {
                        os << e << std::endl;
                    }
                }
            };

            // template<template<typename, typename> class P, typename K, typename V>
            // struct print_log_value<P<K, V>> {
            //     void operator()(std::ostream &, P<K, V> const &) {
            //     }
            // };

        }    // namespace tt_detail
    }    // namespace test_tools
}    // namespace boost

template<typename FieldType, size_t Rate>
void test_poseidon1_permutation(typename poseidon1_policy<FieldType, 128, Rate>::state_type input,
                                typename poseidon1_policy<FieldType, 128, Rate>::state_type expected_result) {
    using poseidon1_policy_t = poseidon1_policy<FieldType, 128, Rate>;
    BOOST_STATIC_ASSERT_MSG(poseidon1_policy_type<poseidon1_policy_t>,
                            "poseidon1_policy must satisfy the Poseidon1 policy concept");
    BOOST_STATIC_ASSERT_MSG(!poseidon2_policy_type<poseidon1_policy_t>,
                            "poseidon1_policy must not satisfy the Poseidon2 policy concept");

    typename poseidon1_policy_t::state_type poseidon1_input = input;
    typename poseidon1_policy_t::state_type optimized_poseidon1_input = input;

    poseidon1_permutation<poseidon1_policy_t>::permute(poseidon1_input);
    BOOST_CHECK_EQUAL(poseidon1_input, expected_result);

    poseidon1_optimized_permutation<poseidon1_policy_t>::permute(optimized_poseidon1_input);
    BOOST_CHECK_EQUAL(optimized_poseidon1_input, expected_result);
    BOOST_CHECK_EQUAL(optimized_poseidon1_input, poseidon1_input);
}

template<typename HashType>
typename HashType::digest_type test_hash_field_elements(std::vector<typename HashType::word_type> input) {
    return hash<HashType>(input);
}

BOOST_AUTO_TEST_SUITE(poseidon_tests)

BOOST_AUTO_TEST_CASE(poseidon1_pad10_sponge_handles_empty_and_variable_length_inputs) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using sponge_type = poseidon_pad10_sponge_construction<policy, permutation>;
    using hash_type = poseidon_pad10_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    BOOST_STATIC_ASSERT_MSG(hashes::is_sponge_construction<sponge_type>::value,
                            "Poseidon Pad10 sponge must satisfy the crypto3 sponge construction trait");

    const typename hash_type::digest_type empty_digest = test_hash_field_elements<hash_type>({});
    const typename hash_type::digest_type zero_digest = test_hash_field_elements<hash_type>({word_type(0u)});
    const typename hash_type::digest_type full_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u)});
    const typename hash_type::digest_type full_plus_zero_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u), word_type(0u)});

    BOOST_CHECK_NE(empty_digest, zero_digest);
    BOOST_CHECK_NE(full_digest, full_plus_zero_digest);
}

BOOST_AUTO_TEST_CASE(poseidon1_pad10_sponge_matches_accumulator_for_partial_final_block) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using sponge_type = poseidon_pad10_sponge_construction<policy, permutation>;
    using hash_type = poseidon_pad10_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    typename policy::block_type first_block = {word_type(0u), word_type(1u)};
    typename policy::block_type final_block = {word_type(2u), word_type(0u)};

    sponge_type sponge;
    sponge.absorb(first_block);
    sponge.absorb_with_padding(final_block, 1);

    const typename hash_type::digest_type sponge_digest = sponge.digest();
    const typename hash_type::digest_type accumulator_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u), word_type(2u)});

    BOOST_CHECK_EQUAL(sponge_digest, accumulator_digest);
}

BOOST_AUTO_TEST_CASE(poseidon1_pad10_sponge_matches_accumulator_for_full_final_block) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using sponge_type = poseidon_pad10_sponge_construction<policy, permutation>;
    using hash_type = poseidon_pad10_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    typename policy::block_type final_block = {word_type(0u), word_type(1u)};

    sponge_type sponge;
    sponge.absorb_with_padding(final_block, policy::block_words);

    const typename hash_type::digest_type sponge_digest = sponge.digest();
    const typename hash_type::digest_type accumulator_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u)});

    BOOST_CHECK_EQUAL(sponge_digest, accumulator_digest);
}

BOOST_AUTO_TEST_CASE(poseidon1_pad10_sponge_two_full_blocks_pads_only_the_final_block) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using hash_type = poseidon_pad10_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    typename policy::state_type expected_state = policy::iv_generator::generate();
    expected_state[0] = word_type(0u);
    expected_state[1] = word_type(1u);
    permutation::permute(expected_state);
    expected_state[0] = word_type(2u);
    expected_state[1] = word_type(3u);
    expected_state[policy::block_words] += word_type(1u);
    permutation::permute(expected_state);

    const typename hash_type::digest_type accumulator_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u), word_type(2u), word_type(3u)});

    BOOST_CHECK_EQUAL(accumulator_digest, expected_state[0]);
}

BOOST_AUTO_TEST_CASE(poseidon1_pad10_sponge_supports_bit_based_multi_word_digest) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2, /*Capacity=*/1, /*DigestBits=*/256>;
    using permutation = poseidon1_permutation<policy>;
    using sponge_type = poseidon_pad10_sponge_construction<policy, permutation>;
    using hash_type = poseidon_pad10_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    BOOST_STATIC_ASSERT_MSG(policy::digest_bits == 256, "Poseidon1 policy should expose requested digest bits");
    BOOST_STATIC_ASSERT_MSG(policy::digest_words == 2,
                            "Poseidon1 policy should round requested digest bits up to field elements");
    BOOST_STATIC_ASSERT_MSG((std::is_same<typename policy::digest_type, std::array<word_type, 2>>::value),
                            "Multi-word Poseidon1 digest should be an array of field elements");

    typename policy::state_type expected_state = policy::iv_generator::generate();
    expected_state[0] = word_type(0u);
    expected_state[1] = word_type(1u);
    expected_state[policy::block_words] += word_type(1u);
    permutation::permute(expected_state);

    typename policy::digest_type expected_digest = {expected_state[0], expected_state[1]};

    sponge_type sponge;
    typename policy::block_type final_block = {word_type(0u), word_type(1u)};
    sponge.absorb_with_padding(final_block, policy::block_words);

    const typename hash_type::digest_type sponge_digest = sponge.digest();
    const typename hash_type::digest_type accumulator_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u)});

    BOOST_CHECK_EQUAL(sponge_digest, expected_digest);
    BOOST_CHECK_EQUAL(accumulator_digest, expected_digest);
}

BOOST_AUTO_TEST_CASE(poseidon1_pad10_add_mode_partial_final_block_adds_into_existing_rate) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using sponge_type = poseidon_sponge_construction<policy,
                                                     permutation,
                                                     poseidon_sponge_absorb_mode::add,
                                                     poseidon_sponge_padding_mode::pad10>;
    using hash_type = poseidon_add_pad10_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    BOOST_STATIC_ASSERT_MSG(hashes::is_sponge_construction<sponge_type>::value,
                            "Poseidon add-mode Pad10 sponge must satisfy the crypto3 sponge construction trait");

    typename policy::state_type expected_state = policy::iv_generator::generate();
    expected_state[0] += word_type(0u);
    expected_state[1] += word_type(1u);
    permutation::permute(expected_state);
    expected_state[0] += word_type(2u);
    expected_state[1] += word_type(1u);
    permutation::permute(expected_state);

    const typename hash_type::digest_type accumulator_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u), word_type(2u)});

    BOOST_CHECK_EQUAL(accumulator_digest, expected_state[0]);
}

BOOST_AUTO_TEST_CASE(poseidon1_pad10_add_mode_full_final_block_adds_into_existing_rate) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using hash_type = poseidon_add_pad10_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    typename policy::state_type expected_state = policy::iv_generator::generate();
    expected_state[0] += word_type(0u);
    expected_state[1] += word_type(1u);
    permutation::permute(expected_state);
    expected_state[0] += word_type(2u);
    expected_state[1] += word_type(3u);
    expected_state[policy::block_words] += word_type(1u);
    permutation::permute(expected_state);

    const typename hash_type::digest_type accumulator_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u), word_type(2u), word_type(3u)});

    BOOST_CHECK_EQUAL(accumulator_digest, expected_state[0]);
}

BOOST_AUTO_TEST_CASE(poseidon1_pad10_sponge_dense_and_optimized_permutations_match) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using dense_permutation = poseidon1_permutation<policy>;
    using optimized_permutation = poseidon1_optimized_permutation<policy>;
    using dense_hash_type = poseidon_pad10_test_hash<policy, dense_permutation>;
    using optimized_hash_type = poseidon_pad10_test_hash<policy, optimized_permutation>;
    using word_type = typename policy::word_type;

    std::vector<word_type> input = {word_type(0u), word_type(1u), word_type(2u), word_type(3u), word_type(4u)};

    const typename dense_hash_type::digest_type dense_digest = test_hash_field_elements<dense_hash_type>(input);
    const typename optimized_hash_type::digest_type optimized_digest =
        test_hash_field_elements<optimized_hash_type>(input);

    BOOST_CHECK_EQUAL(dense_digest, optimized_digest);
}

BOOST_AUTO_TEST_CASE(poseidon1_public_wrapper_defaults_to_optimized_pad10) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using optimized_permutation = poseidon1_optimized_permutation<policy>;
    using optimized_hash_type = poseidon_pad10_test_hash<policy, optimized_permutation>;
    using poseidon_hash_type = hashes::poseidon<policy>;
    using public_hash_type = hashes::poseidon1<policy>;
    using dense_public_hash_type = hashes::poseidon1_dense<policy>;
    using word_type = typename policy::word_type;

    BOOST_STATIC_ASSERT_MSG(hashes::is_poseidon<poseidon_hash_type>::value,
                            "poseidon public wrapper should be recognized as Poseidon");
    BOOST_STATIC_ASSERT_MSG(hashes::is_poseidon<public_hash_type>::value,
                            "poseidon1 public wrapper should be recognized as Poseidon");
    BOOST_STATIC_ASSERT_MSG(hashes::is_poseidon<dense_public_hash_type>::value,
                            "poseidon1 dense public wrapper should be recognized as Poseidon");

    std::vector<word_type> input = {word_type(0u), word_type(1u), word_type(2u), word_type(3u), word_type(4u)};

    const typename poseidon_hash_type::digest_type poseidon_digest =
        test_hash_field_elements<poseidon_hash_type>(input);
    const typename public_hash_type::digest_type public_digest = test_hash_field_elements<public_hash_type>(input);
    const typename dense_public_hash_type::digest_type dense_public_digest =
        test_hash_field_elements<dense_public_hash_type>(input);
    const typename optimized_hash_type::digest_type optimized_digest =
        test_hash_field_elements<optimized_hash_type>(input);

    BOOST_CHECK_EQUAL(poseidon_digest, optimized_digest);
    BOOST_CHECK_EQUAL(poseidon_digest, public_digest);
    BOOST_CHECK_EQUAL(public_digest, optimized_digest);
    BOOST_CHECK_EQUAL(public_digest, dense_public_digest);
}

BOOST_AUTO_TEST_CASE(poseidon1_padding_free_overwrite_handles_empty_input_without_permutation) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using sponge_type = poseidon_padding_free_sponge_construction<policy, permutation>;
    using hash_type = poseidon_padding_free_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    BOOST_STATIC_ASSERT_MSG(hashes::is_sponge_construction<sponge_type>::value,
                            "Poseidon padding-free sponge must satisfy the crypto3 sponge construction trait");

    const typename hash_type::digest_type empty_digest = test_hash_field_elements<hash_type>({});

    BOOST_CHECK_EQUAL(empty_digest, word_type(0u));
}

BOOST_AUTO_TEST_CASE(poseidon1_padding_free_overwrite_matches_accumulator_for_full_final_block) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using sponge_type = poseidon_padding_free_sponge_construction<policy, permutation>;
    using hash_type = poseidon_padding_free_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    typename policy::block_type final_block = {word_type(0u), word_type(1u)};

    sponge_type sponge;
    sponge.absorb_with_padding(final_block, policy::block_words);

    const typename hash_type::digest_type sponge_digest = sponge.digest();
    const typename hash_type::digest_type accumulator_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u)});

    BOOST_CHECK_EQUAL(sponge_digest, accumulator_digest);
}

BOOST_AUTO_TEST_CASE(poseidon1_padding_free_overwrite_two_full_blocks_are_permuted_without_padding) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using hash_type = poseidon_padding_free_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    typename policy::state_type expected_state = policy::iv_generator::generate();
    expected_state[0] = word_type(0u);
    expected_state[1] = word_type(1u);
    permutation::permute(expected_state);
    expected_state[0] = word_type(2u);
    expected_state[1] = word_type(3u);
    permutation::permute(expected_state);

    const typename hash_type::digest_type accumulator_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u), word_type(2u), word_type(3u)});

    BOOST_CHECK_EQUAL(accumulator_digest, expected_state[0]);
}

BOOST_AUTO_TEST_CASE(poseidon1_padding_free_overwrite_partial_final_block_leaves_suffix_unchanged) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using hash_type = poseidon_padding_free_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    typename policy::state_type expected_state = policy::iv_generator::generate();
    expected_state[0] = word_type(0u);
    expected_state[1] = word_type(1u);
    permutation::permute(expected_state);
    expected_state[0] = word_type(2u);
    permutation::permute(expected_state);

    const typename hash_type::digest_type accumulator_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u), word_type(2u)});

    BOOST_CHECK_EQUAL(accumulator_digest, expected_state[0]);
}

BOOST_AUTO_TEST_CASE(poseidon1_padding_free_add_mode_partial_final_block_adds_into_existing_rate) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using sponge_type = poseidon_sponge_construction<policy,
                                                     permutation,
                                                     poseidon_sponge_absorb_mode::add,
                                                     poseidon_sponge_padding_mode::padding_free>;
    using hash_type = poseidon_add_padding_free_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    BOOST_STATIC_ASSERT_MSG(hashes::is_sponge_construction<sponge_type>::value,
                            "Poseidon add-mode padding-free sponge must satisfy the crypto3 sponge construction trait");

    typename policy::state_type expected_state = policy::iv_generator::generate();
    expected_state[0] += word_type(0u);
    expected_state[1] += word_type(1u);
    permutation::permute(expected_state);
    expected_state[0] += word_type(2u);
    permutation::permute(expected_state);

    const typename hash_type::digest_type accumulator_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u), word_type(2u)});

    BOOST_CHECK_EQUAL(accumulator_digest, expected_state[0]);
}

BOOST_AUTO_TEST_CASE(poseidon1_padding_free_add_mode_full_final_block_adds_into_existing_rate) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using hash_type = poseidon_add_padding_free_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    typename policy::state_type expected_state = policy::iv_generator::generate();
    expected_state[0] += word_type(0u);
    expected_state[1] += word_type(1u);
    permutation::permute(expected_state);
    expected_state[0] += word_type(2u);
    expected_state[1] += word_type(3u);
    permutation::permute(expected_state);

    const typename hash_type::digest_type accumulator_digest =
        test_hash_field_elements<hash_type>({word_type(0u), word_type(1u), word_type(2u), word_type(3u)});

    BOOST_CHECK_EQUAL(accumulator_digest, expected_state[0]);
}

BOOST_AUTO_TEST_CASE(poseidon1_padding_free_overwrite_differs_from_pad10_on_full_final_block) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon1_permutation<policy>;
    using pad10_hash_type = poseidon_pad10_test_hash<policy, permutation>;
    using padding_free_hash_type = poseidon_padding_free_test_hash<policy, permutation>;
    using word_type = typename policy::word_type;

    const typename pad10_hash_type::digest_type pad10_digest =
        test_hash_field_elements<pad10_hash_type>({word_type(0u), word_type(1u)});
    const typename padding_free_hash_type::digest_type padding_free_digest =
        test_hash_field_elements<padding_free_hash_type>({word_type(0u), word_type(1u)});

    BOOST_CHECK_NE(pad10_digest, padding_free_digest);
}

BOOST_AUTO_TEST_CASE(poseidon1_padding_free_overwrite_dense_and_optimized_permutations_match) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using dense_permutation = poseidon1_permutation<policy>;
    using optimized_permutation = poseidon1_optimized_permutation<policy>;
    using dense_hash_type = poseidon_padding_free_test_hash<policy, dense_permutation>;
    using optimized_hash_type = poseidon_padding_free_test_hash<policy, optimized_permutation>;
    using word_type = typename policy::word_type;

    std::vector<word_type> input = {word_type(0u), word_type(1u), word_type(2u), word_type(3u), word_type(4u)};

    const typename dense_hash_type::digest_type dense_digest = test_hash_field_elements<dense_hash_type>(input);
    const typename optimized_hash_type::digest_type optimized_digest =
        test_hash_field_elements<optimized_hash_type>(input);

    BOOST_CHECK_EQUAL(dense_digest, optimized_digest);
}

BOOST_AUTO_TEST_CASE(poseidon1_padding_free_public_wrapper_matches_optimized_padding_free) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon1_policy<field_type, 128, /*Rate=*/2>;
    using optimized_permutation = poseidon1_optimized_permutation<policy>;
    using optimized_hash_type = poseidon_padding_free_test_hash<policy, optimized_permutation>;
    using public_hash_type = hashes::poseidon1_padding_free<policy>;
    using word_type = typename policy::word_type;

    BOOST_STATIC_ASSERT_MSG(hashes::is_poseidon<public_hash_type>::value,
                            "poseidon1 padding-free public wrapper should be recognized as Poseidon");

    std::vector<word_type> input = {word_type(0u), word_type(1u), word_type(2u), word_type(3u), word_type(4u)};

    const typename public_hash_type::digest_type public_digest = test_hash_field_elements<public_hash_type>(input);
    const typename optimized_hash_type::digest_type optimized_digest =
        test_hash_field_elements<optimized_hash_type>(input);

    BOOST_CHECK_EQUAL(public_digest, optimized_digest);
}

BOOST_AUTO_TEST_CASE(poseidon2_external_linear_layer_supports_multiple_widths) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using word_type = typename field_type::value_type;

    using width2_policy =
        base_poseidon2_policy<field_type, 128, /*Rate=*/1, /*Capacity=*/1, 5, 8, 56, field_type::value_bits>;
    using width3_policy =
        base_poseidon2_policy<field_type, 128, /*Rate=*/2, /*Capacity=*/1, 5, 8, 56, field_type::value_bits>;
    using width4_policy =
        base_poseidon2_policy<field_type, 128, /*Rate=*/3, /*Capacity=*/1, 5, 8, 56, field_type::value_bits>;

    BOOST_STATIC_ASSERT_MSG(poseidon2_policy_type<width2_policy>, "Poseidon2 round helpers should support width 2");
    BOOST_STATIC_ASSERT_MSG(poseidon2_policy_type<width3_policy>, "Poseidon2 round helpers should support width 3");
    BOOST_STATIC_ASSERT_MSG(poseidon2_policy_type<width4_policy>, "Poseidon2 round helpers should support width 4");

    typename width2_policy::state_type width2_state = {word_type(0u), word_type(1u)};
    poseidon2_round_functions<width2_policy>::external_linear_layer(width2_state);
    BOOST_CHECK_EQUAL(width2_state, (typename width2_policy::state_type {word_type(1u), word_type(2u)}));

    typename width3_policy::state_type width3_state = {word_type(0u), word_type(1u), word_type(2u)};
    poseidon2_round_functions<width3_policy>::external_linear_layer(width3_state);
    BOOST_CHECK_EQUAL(width3_state, (typename width3_policy::state_type {word_type(3u), word_type(4u), word_type(5u)}));

    typename width4_policy::state_type width4_state = {word_type(0u), word_type(1u), word_type(2u), word_type(3u)};
    poseidon2_round_functions<width4_policy>::external_linear_layer(width4_state);
    BOOST_CHECK_EQUAL(
        width4_state,
        (typename width4_policy::state_type {word_type(36u), word_type(22u), word_type(68u), word_type(54u)}));
}

BOOST_AUTO_TEST_CASE(poseidon2_bn254_width3_permutation_matches_reference_vector) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon2_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon2_permutation<policy>;

    BOOST_STATIC_ASSERT_MSG(poseidon2_policy_type<policy>, "poseidon2_policy must satisfy the Poseidon2 concept");
    BOOST_STATIC_ASSERT_MSG(!poseidon1_policy_type<policy>,
                            "poseidon2_policy must not satisfy the Poseidon1 policy concept");
    BOOST_STATIC_ASSERT_MSG(policy::state_words == 3, "The checked-in BN254 Poseidon2 policy uses width 3");
    BOOST_STATIC_ASSERT_MSG(policy::sbox_power == 5, "BN254 Poseidon2 uses x^5");
    BOOST_STATIC_ASSERT_MSG(policy::full_rounds == 8, "BN254 Poseidon2 uses RF=8");
    BOOST_STATIC_ASSERT_MSG(policy::part_rounds == 56, "BN254 Poseidon2 width 3 uses RP=56");

    typename policy::state_type state = {
        0x0000000000000000000000000000000000000000000000000000000000000000_cppui_modular254,
        0x0000000000000000000000000000000000000000000000000000000000000001_cppui_modular254,
        0x0000000000000000000000000000000000000000000000000000000000000002_cppui_modular254};
    const typename policy::state_type expected_state = {
        0x0bb61d24daca55eebcb1929a82650f328134334da98ea4f847f760054f4a3033_cppui_modular254,
        0x303b6f7c86d043bfcbcc80214f26a30277a15d3f74ca654992defe7ff8d03570_cppui_modular254,
        0x1ed25194542b12eef8617361c3ba7c52e660b145994427cc86296242cf766ec8_cppui_modular254};

    permutation::permute(state);

    BOOST_CHECK_EQUAL(state, expected_state);
}

BOOST_AUTO_TEST_CASE(poseidon2_public_wrapper_defaults_to_pad10) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon2_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon2_permutation<policy>;
    using sponge_hash_type = poseidon_pad10_test_hash<policy, permutation>;
    using public_hash_type = hashes::poseidon2<policy>;
    using word_type = typename policy::word_type;

    BOOST_STATIC_ASSERT_MSG(hashes::is_poseidon<public_hash_type>::value,
                            "poseidon2 public wrapper should be recognized as Poseidon");

    std::vector<word_type> input = {word_type(0u), word_type(1u), word_type(2u), word_type(3u), word_type(4u)};

    const typename public_hash_type::digest_type public_digest = test_hash_field_elements<public_hash_type>(input);
    const typename sponge_hash_type::digest_type sponge_digest = test_hash_field_elements<sponge_hash_type>(input);

    BOOST_CHECK_EQUAL(public_digest, sponge_digest);
}

BOOST_AUTO_TEST_CASE(poseidon2_padding_free_public_wrapper_matches_padding_free_sponge) {
    using field_type = fields::alt_bn128_scalar_field<254>;
    using policy = poseidon2_policy<field_type, 128, /*Rate=*/2>;
    using permutation = poseidon2_permutation<policy>;
    using sponge_hash_type = poseidon_padding_free_test_hash<policy, permutation>;
    using public_hash_type = hashes::poseidon2_padding_free<policy>;
    using word_type = typename policy::word_type;

    BOOST_STATIC_ASSERT_MSG(hashes::is_poseidon<public_hash_type>::value,
                            "poseidon2 padding-free wrapper should be recognized as Poseidon");

    std::vector<word_type> input = {word_type(0u), word_type(1u), word_type(2u), word_type(3u)};

    const typename public_hash_type::digest_type public_digest = test_hash_field_elements<public_hash_type>(input);
    const typename sponge_hash_type::digest_type sponge_digest = test_hash_field_elements<sponge_hash_type>(input);

    BOOST_CHECK_EQUAL(public_digest, sponge_digest);
}

// Poseidon1 permutation test vectors are taken from:
//   https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/208b5a164c6a252b137997694d90931b2bb851c5/code/test_vectors.txt
BOOST_AUTO_TEST_CASE(poseidon1_permutation_254_2) {
    test_poseidon1_permutation<fields::alt_bn128_scalar_field<254>, 2>(
        {0x0000000000000000000000000000000000000000000000000000000000000000_cppui_modular254,
         0x0000000000000000000000000000000000000000000000000000000000000001_cppui_modular254,
         0x0000000000000000000000000000000000000000000000000000000000000002_cppui_modular254},
        {0x115cc0f5e7d690413df64c6b9662e9cf2a3617f2743245519e19607a4417189a_cppui_modular254,
         0x0fca49b798923ab0239de1c9e7a4a9a2210312b6a2f616d18b5a87f9b628ae29_cppui_modular254,
         0x0e7ae82e40091e63cbd4f16a6d16310b3729d4b6e138fcf54110e2867045a30c_cppui_modular254});
}

BOOST_AUTO_TEST_CASE(poseidon1_permutation_254_4) {
    test_poseidon1_permutation<fields::alt_bn128_scalar_field<254>, 4>(
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

BOOST_AUTO_TEST_CASE(poseidon1_permutation_255_3) {
    test_poseidon1_permutation<fields::bls12_scalar_field<381>, 2>(
        {0x0000000000000000000000000000000000000000000000000000000000000000_cppui_modular255,
         0x0000000000000000000000000000000000000000000000000000000000000001_cppui_modular255,
         0x0000000000000000000000000000000000000000000000000000000000000002_cppui_modular255},
        {0x28ce19420fc246a05553ad1e8c98f5c9d67166be2c18e9e4cb4b4e317dd2a78a_cppui_modular255,
         0x51f3e312c95343a896cfd8945ea82ba956c1118ce9b9859b6ea56637b4b1ddc4_cppui_modular255,
         0x3b2b69139b235626a0bfb56c9527ae66a7bf486ad8c11c14d1da0c69bbe0f79a_cppui_modular255});
}

BOOST_AUTO_TEST_CASE(poseidon1_permutation_255_4) {
    test_poseidon1_permutation<fields::bls12_scalar_field<381>, 4>(
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

BOOST_AUTO_TEST_SUITE_END()
