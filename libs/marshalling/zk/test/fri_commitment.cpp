//---------------------------------------------------------------------------//
// Copyright (c) 2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2021 Nikita Kaskov <nbering@nil.foundation>
// Copyright (c) 2021 Ilias Khairullin <ilias@nil.foundation>
// Copyright (c) 2022-2023 Elena Tatuzova <e.tatuzova@nil.foundation>
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

#define BOOST_TEST_MODULE crypto3_marshalling_fri_commitment_test

#include <boost/test/unit_test.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <iostream>
#include <iomanip>
#include <random>
#include <regex>

#include <nil/marshalling/status_type.hpp>
#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/endianness.hpp>

#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>
#include <boost/multiprecision/number.hpp>

#include <nil/crypto3/algebra/random_element.hpp>
#include <nil/crypto3/algebra/curves/bls12.hpp>
#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/curves/pallas.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/pallas.hpp>
#include <nil/crypto3/algebra/curves/mnt4.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/mnt4.hpp>

#include <nil/crypto3/math/polynomial/polynomial.hpp>
#include <nil/crypto3/math/polynomial/lagrange_interpolation.hpp>
#include <nil/crypto3/math/algorithms/unity_root.hpp>
#include <nil/crypto3/math/domains/evaluation_domain.hpp>
#include <nil/crypto3/math/algorithms/make_evaluation_domain.hpp>
#include <nil/crypto3/math/algorithms/calculate_domain_set.hpp>

#include <nil/crypto3/hash/type_traits.hpp>
#include <nil/crypto3/hash/sha2.hpp>
#include <nil/crypto3/hash/keccak.hpp>

#include <nil/crypto3/random/algebraic_random_device.hpp>

#include <nil/crypto3/zk/commitments/detail/polynomial/basic_fri.hpp>
#include <nil/crypto3/zk/commitments/polynomial/fri.hpp>
#include <nil/crypto3/zk/test_tools/random_test_initializer.hpp>
#include <nil/crypto3/marshalling/zk/types/commitments/fri.hpp>

using namespace nil::crypto3;

/*
template<typename TIter>
void print_byteblob(std::ostream &os, TIter iter_begin, TIter iter_end) {
    for (TIter it = iter_begin; it != iter_end; it++) {
        os << std::hex << int(*it) << std::endl;
    }
}

template<typename FieldParams>
void print_field_element(std::ostream &os,
                         const typename nil::crypto3::algebra::fields::detail::element_fp<FieldParams> &e) {
    os << e.data << std::endl;
}

template<typename FieldParams>
void print_field_element(std::ostream &os,
                         const typename nil::crypto3::algebra::fields::detail::element_fp2<FieldParams> &e) {
    os << "[" << e.data[0].data << "," << e.data[1].data << "]" << std::endl;
}
*/
template<typename ValueType, std::size_t N>
typename std::enable_if<std::is_unsigned<ValueType>::value, std::vector<std::array<ValueType, N>>>::type
    generate_random_data(std::size_t leaf_number, boost::random::mt11213b &rnd) {
    std::vector<std::array<ValueType, N>> v;
    for (std::size_t i = 0; i < leaf_number; ++i) {
        std::array<ValueType, N> leaf;
        std::generate(std::begin(leaf), std::end(leaf),
                      [&]() { return rnd() % (std::numeric_limits<ValueType>::max() + 1); });
        v.emplace_back(leaf);
    }
    return v;
}

std::vector<std::vector<std::uint8_t>> generate_random_data_for_merkle_tree(size_t leafs_number, size_t leaf_bytes,
                                                                            boost::random::mt11213b &rnd) {
    std::vector<std::vector<std::uint8_t>> rdata(leafs_number, std::vector<std::uint8_t>(leaf_bytes));

    for (std::size_t i = 0; i < leafs_number; ++i) {
        std::vector<uint8_t> leaf(leaf_bytes);
        for (size_t i = 0; i < leaf_bytes; i++) {
            leaf[i] = rnd() % (std::numeric_limits<std::uint8_t>::max() + 1);
        }
        rdata.emplace_back(leaf);
    }
    return rdata;
}

template<typename FRI>
typename FRI::merkle_proof_type generate_random_merkle_proof(std::size_t tree_depth, boost::random::mt11213b &rnd) {
    std::size_t leafs_number = 1 << tree_depth;
    std::size_t leaf_size = 32;

    auto rdata1 = generate_random_data_for_merkle_tree(leafs_number, leaf_size, rnd);
    auto tree1 =
        containers::make_merkle_tree<typename FRI::merkle_tree_hash_type, FRI::m>(rdata1.begin(), rdata1.end());
    std::size_t idx1 = rnd() % leafs_number;
    typename FRI::merkle_proof_type mp1(tree1, idx1);
    return mp1;
}

inline std::vector<std::size_t> generate_random_step_list(const std::size_t r, const std::size_t max_step,
                                                          boost::random::mt11213b &rnd) {
    using dist_type = std::uniform_int_distribution<int>;

    std::vector<std::size_t> step_list;
    std::size_t steps_sum = 0;
    while (steps_sum != r) {
        if (r - steps_sum <= max_step) {
            while (r - steps_sum != 1) {
                step_list.emplace_back(r - steps_sum - 1);
                steps_sum += step_list.back();
            }
            step_list.emplace_back(1);
            steps_sum += step_list.back();
        } else {
            step_list.emplace_back(dist_type(1, max_step)(rnd));
            steps_sum += step_list.back();
        }
    }

    return step_list;
}

template<typename FRI>
typename FRI::polynomial_values_type
    generate_random_polynomial_values(size_t step,
                                      nil::crypto3::random::algebraic_engine<typename FRI::field_type> &alg_rnd) {
    typename FRI::polynomial_values_type values;

    std::size_t coset_size = 1 << (step - 1);
    values.resize(coset_size);
    for (size_t i = 0; i < coset_size; i++) {
        for (size_t j = 0; j < FRI::m; j++) {
            values[i][j] = alg_rnd();
            values[i][j] = alg_rnd();
        }
    }
    return values;
}

template<typename FieldType>
math::polynomial<typename FieldType::value_type>
    generate_random_polynomial(size_t degree, nil::crypto3::random::algebraic_engine<FieldType> &d) {
    math::polynomial<typename FieldType::value_type> poly;
    poly.resize(degree);

    for (std::size_t i = 0; i < degree; ++i) {
        poly[i] = d();
    }
    return poly;
}

template<typename FRI>
typename FRI::round_proof_type
    generate_random_fri_round_proof(std::size_t r_i,
                                    nil::crypto3::random::algebraic_engine<typename FRI::field_type> &alg_rnd,
                                    boost::random::mt11213b &rnd) {
    typename FRI::round_proof_type res;
    res.p = generate_random_merkle_proof<FRI>(3, rnd);
    res.y = generate_random_polynomial_values<FRI>(r_i, alg_rnd);

    return res;
}

template<typename FRI>
typename FRI::initial_proof_type
    generate_random_fri_initial_proof(std::size_t polynomial_number,
                                      std::size_t r0,
                                      nil::crypto3::random::algebraic_engine<typename FRI::field_type> &alg_rnd,
                                      boost::random::mt11213b &rnd) {
    typename FRI::initial_proof_type res;

    std::size_t coset_size = 1 << r0;
    res.p = generate_random_merkle_proof<FRI>(3, rnd);
    res.values.resize(polynomial_number);
    for (std::size_t i = 0; i < polynomial_number; i++) {
        res.values[i].resize(coset_size / FRI::m);
        for (std::size_t j = 0; j < coset_size / FRI::m; j++) {
            res.values[i][j][0] = alg_rnd();
            res.values[i][j][1] = alg_rnd();
        }
    }

    return res;
}

template<typename FRI>
typename FRI::query_proof_type
    generate_random_fri_query_proof(std::size_t max_batch_size,
                                    std::vector<std::size_t>
                                        step_list,
                                    nil::crypto3::marshalling::types::batch_info_type batch_info,
                                    nil::crypto3::random::algebraic_engine<typename FRI::field_type> &alg_rnd,
                                    boost::random::mt11213b &rnd) {
    typename FRI::query_proof_type res;

    for (const auto &it : batch_info) {
        res.initial_proof[it.first] = generate_random_fri_initial_proof<FRI>(it.second, step_list[0], alg_rnd, rnd);
    }
    res.round_proofs.resize(step_list.size());
    for (std::size_t i = 1; i < step_list.size(); i++) {
        res.round_proofs[i - 1] = generate_random_fri_round_proof<FRI>(step_list[i], alg_rnd, rnd);
    }
    res.round_proofs[step_list.size() - 1] = generate_random_fri_round_proof<FRI>(1, alg_rnd, rnd);
    return res;
}

template<typename FRI>
typename FRI::proof_type
    generate_random_fri_proof(std::size_t d,    // final polynomial degree
                              std::size_t max_batch_size,
                              std::vector<std::size_t>
                                  step_list,
                              std::size_t lambda,
                              bool use_grinding,
                              nil::crypto3::marshalling::types::batch_info_type batch_info,
                              nil::crypto3::random::algebraic_engine<typename FRI::field_type> &alg_rnd,
                              boost::random::mt11213b &rnd) {
    typename FRI::proof_type res;
    res.query_proofs.resize(lambda);
    for (std::size_t k = 0; k < lambda; k++) {
        res.query_proofs[k] = generate_random_fri_query_proof<FRI>(max_batch_size, step_list, batch_info, alg_rnd, rnd);
    }
    res.fri_roots.resize(step_list.size());
    for (std::size_t k = 0; k < step_list.size(); k++) {
        res.fri_roots[k] = nil::crypto3::hash<typename FRI::merkle_tree_hash_type>(
            generate_random_data<std::uint8_t, 32>(1, rnd).at(0));
    }
    if (use_grinding) {
        res.proof_of_work = rnd();
    }
    res.final_polynomial = generate_random_polynomial<typename FRI::field_type>(d, alg_rnd);
    return res;
}

template<typename Endianness, typename FRI>
void test_fri_proof(typename FRI::proof_type &proof,
                    typename nil::crypto3::marshalling::types::batch_info_type batch_info,
                    const typename FRI::params_type &params) {
    using TTypeBase = nil::marshalling::field_type<Endianness>;

    auto filled_proof = nil::crypto3::marshalling::types::fill_fri_proof<Endianness, FRI>(proof, batch_info, params);
    auto _proof = nil::crypto3::marshalling::types::make_fri_proof<Endianness, FRI>(filled_proof, batch_info);
    BOOST_CHECK(proof.fri_roots == _proof.fri_roots);
    BOOST_CHECK(proof.final_polynomial == _proof.final_polynomial);
    BOOST_CHECK(proof.query_proofs[0].initial_proof == _proof.query_proofs[0].initial_proof);
    BOOST_CHECK(proof.query_proofs[0].round_proofs.size() == _proof.query_proofs[0].round_proofs.size());
    for (std::size_t i = 0; i < proof.query_proofs[0].round_proofs.size(); i++) {
        if (proof.query_proofs[0].round_proofs[i] != _proof.query_proofs[0].round_proofs[i]) {
            if (proof.query_proofs[0].round_proofs[i].p != _proof.query_proofs[0].round_proofs[i].p)
                std::cout << "round proof " << i << "merkle proof is not equal" << std::endl;
            if (proof.query_proofs[0].round_proofs[i].y != _proof.query_proofs[0].round_proofs[i].y)
                std::cout << "round proof " << i << "poly values are not equal" << std::endl;
        }
    }
    BOOST_CHECK(proof.query_proofs[0] == _proof.query_proofs[0]);
    BOOST_CHECK(proof.proof_of_work == _proof.proof_of_work);
    BOOST_CHECK(proof == _proof);

    std::vector<std::uint8_t> cv;
    cv.resize(filled_proof.length(), 0x00);
    auto write_iter = cv.begin();
    auto status = filled_proof.write(write_iter, cv.size());
    BOOST_CHECK(status == nil::marshalling::status_type::success);

    typename nil::crypto3::marshalling::types::fri_proof<TTypeBase, FRI>::type test_val_read;
    auto read_iter = cv.begin();
    status = test_val_read.read(read_iter, cv.size());
    BOOST_CHECK(status == nil::marshalling::status_type::success);
    typename FRI::proof_type constructed_val_read =
        nil::crypto3::marshalling::types::make_fri_proof<Endianness, FRI>(test_val_read, batch_info);
    BOOST_CHECK(proof == constructed_val_read);
}

BOOST_FIXTURE_TEST_SUITE(marshalling_fri_proof_elements,
                         zk::test_tools::random_test_initializer<algebra::curves::bls12<381>::scalar_field_type>)
static constexpr std::size_t lambda = 40;
static constexpr std::size_t m = 2;

using curve_type = nil::crypto3::algebra::curves::bls12<381>;
using field_type = typename curve_type::scalar_field_type;
using value_type = typename field_type::value_type;
using hash_type = nil::crypto3::hashes::keccak_1600<256>;

using Endianness = nil::marshalling::option::big_endian;
using TTypeBase = nil::marshalling::field_type<Endianness>;
using FRI = typename nil::crypto3::zk::commitments::detail::basic_batched_fri<field_type, hash_type, hash_type, m>;

BOOST_AUTO_TEST_CASE(polynomial_test) {
    using polynomial_type = math::polynomial<typename field_type::value_type>;
    polynomial_type f = {{1u, 3u, 4u, 1u, 5u, 6u, 7u, 2u, 8u, 7u, 5u, 6u, 1u, 2u, 1u, 1u}};
    auto filled_polynomial = nil::crypto3::marshalling::types::fill_fri_math_polynomial<Endianness, polynomial_type>(f);

    auto _f =
        nil::crypto3::marshalling::types::make_fri_math_polynomial<Endianness, polynomial_type>(filled_polynomial);
    BOOST_CHECK(f == _f);

    f = generate_random_polynomial<field_type>(2048, alg_random_engines.template get_alg_engine<field_type>());
    filled_polynomial = nil::crypto3::marshalling::types::fill_fri_math_polynomial<Endianness, polynomial_type>(f);

    _f = nil::crypto3::marshalling::types::make_fri_math_polynomial<Endianness, polynomial_type>(filled_polynomial);
    BOOST_CHECK(f == _f);
}

BOOST_AUTO_TEST_CASE(merkle_proof_vector_test) {
    std::vector<typename FRI::merkle_proof_type> mp;
    mp.push_back(generate_random_merkle_proof<FRI>(5, generic_random_engine));
    mp.push_back(generate_random_merkle_proof<FRI>(5, generic_random_engine));
    mp.push_back(generate_random_merkle_proof<FRI>(5, generic_random_engine));
    mp.push_back(generate_random_merkle_proof<FRI>(5, generic_random_engine));
    mp.push_back(generate_random_merkle_proof<FRI>(5, generic_random_engine));
    mp.push_back(generate_random_merkle_proof<FRI>(5, generic_random_engine));
    mp.push_back(generate_random_merkle_proof<FRI>(5, generic_random_engine));

    auto filled = nil::crypto3::marshalling::types::fill_merkle_proof_vector<Endianness, FRI>(mp);
    auto _f = nil::crypto3::marshalling::types::make_merkle_proof_vector<Endianness, FRI>(filled);
    BOOST_CHECK(mp == _f);

    using TTypeBase = nil::marshalling::field_type<Endianness>;
    std::vector<std::uint8_t> cv;
    cv.resize(filled.length(), 0x00);
    auto write_iter = cv.begin();
    auto status = filled.write(write_iter, cv.size());
    BOOST_CHECK(status == nil::marshalling::status_type::success);

    nil::crypto3::marshalling::types::merkle_proof_vector_type<TTypeBase, FRI> test_val_read;
    auto read_iter = cv.begin();
    test_val_read.read(read_iter, cv.size());
    BOOST_CHECK(status == nil::marshalling::status_type::success);
    auto constructed_val_read =
        nil::crypto3::marshalling::types::make_merkle_proof_vector<Endianness, FRI>(test_val_read);
    BOOST_CHECK(mp == constructed_val_read);
}

BOOST_AUTO_TEST_CASE(fri_proof_test) {
    nil::crypto3::marshalling::types::batch_info_type batch_info;
    batch_info[0] = 1;
    batch_info[1] = 5;
    batch_info[3] = 6;
    batch_info[4] = 3;

    typename FRI::params_type fri_params(1, 11, lambda, 4);

    auto proof =
        generate_random_fri_proof<FRI>(2, 5, fri_params.step_list, lambda, false, batch_info,
                                       alg_random_engines.template get_alg_engine<field_type>(), generic_random_engine);
    test_fri_proof<Endianness, FRI>(proof, batch_info, fri_params);
}

BOOST_AUTO_TEST_CASE(fri_grinding_proof_test) {
    nil::crypto3::marshalling::types::batch_info_type batch_info;
    batch_info[0] = 1;
    batch_info[1] = 5;
    batch_info[3] = 6;
    batch_info[4] = 3;

    typename FRI::params_type fri_params(1, 11, lambda, 4, true);

    auto proof =
        generate_random_fri_proof<FRI>(2, 5, fri_params.step_list, lambda, true, batch_info,
                                       alg_random_engines.template get_alg_engine<field_type>(), generic_random_engine);
    test_fri_proof<Endianness, FRI>(proof, batch_info, fri_params);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(marshalling_real_fri_proofs,
                         zk::test_tools::random_test_initializer<algebra::curves::pallas::base_field_type>)
using Endianness = nil::marshalling::option::big_endian;

BOOST_AUTO_TEST_CASE(marshalling_fri_basic_test) {
    // setup
    using curve_type = algebra::curves::pallas;
    using field_type = typename curve_type::base_field_type;

    typedef hashes::sha2<256> merkle_hash_type;
    typedef hashes::sha2<256> transcript_hash_type;

    constexpr static const std::size_t d = 16;

    constexpr static const std::size_t r = boost::static_log2<d>::value;
    constexpr static const std::size_t m = 2;
    constexpr static const std::size_t lambda = 40;

    typedef zk::commitments::fri<field_type, merkle_hash_type, transcript_hash_type, m> fri_type;

    static_assert(zk::is_commitment<fri_type>::value);
    static_assert(!zk::is_commitment<merkle_hash_type>::value);

    typedef typename fri_type::proof_type proof_type;
    typedef typename fri_type::params_type params_type;

    std::size_t extended_log = boost::static_log2<d>::value;
    std::vector<std::shared_ptr<math::evaluation_domain<field_type>>> D =
        math::calculate_domain_set<field_type>(extended_log, r);

    params_type params(d - 1,    // max_degree
                       D,
                       generate_random_step_list(r, 3, generic_random_engine),
                       2,    // expand_factor
                       lambda);

    BOOST_CHECK(D[1]->m == D[0]->m / 2);
    BOOST_CHECK(D[1]->get_domain_element(1) == D[0]->get_domain_element(1).squared());

    // commit
    math::polynomial<typename field_type::value_type> f = {
        {1u, 3u, 4u, 1u, 5u, 6u, 7u, 2u, 8u, 7u, 5u, 6u, 1u, 2u, 1u, 1u}};
    std::array<std::vector<math::polynomial<typename field_type::value_type>>, 1> fs;
    fs[0].resize(1);
    fs[0][0] = f;
    typename fri_type::merkle_tree_type tree =
        zk::algorithms::precommit<fri_type>(fs[0], params.D[0], params.step_list[0]);
    auto root = zk::algorithms::commit<fri_type>(tree);

    // eval
    std::vector<std::uint8_t> init_blob {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    zk::transcript::fiat_shamir_heuristic_sequential<transcript_hash_type> transcript(init_blob);

    proof_type proof = zk::algorithms::proof_eval<fri_type>(f, tree, params, transcript);
    nil::crypto3::marshalling::types::batch_info_type batch_info;
    batch_info[0] = 1;
    test_fri_proof<Endianness, fri_type>(proof, batch_info, params);

    // verify
    // zk::transcript::fiat_shamir_heuristic_sequential<transcript_hash_type> transcript_verifier(init_blob);

    // BOOST_CHECK(zk::algorithms::verify_eval<fri_type>(proof, root, params, transcript_verifier));
}
BOOST_AUTO_TEST_SUITE_END()
