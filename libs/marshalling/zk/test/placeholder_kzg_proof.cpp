//---------------------------------------------------------------------------//
// Copyright (c) 2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2021 Nikita Kaskov <nbering@nil.foundation>
// Copyright (c) 2021 Ilias Khairullin <ilias@nil.foundation>
// Copyright (c) 2022-2023 Elena Tatuzova <e.tatuzova@nil.foundation>
// Copyright (c) 2024 Vasiliy Olekhov <vasiliy.olekhov@nil.foundation>
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

#define BOOST_TEST_MODULE crypto3_marshalling_placeholder_kzg_proof_test

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <nil/marshalling/endianness.hpp>
#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/status_type.hpp>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/curves/bls12.hpp>
#include <nil/crypto3/algebra/curves/mnt4.hpp>
#include <nil/crypto3/algebra/curves/mnt6.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/bls12.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/mnt4.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/mnt6.hpp>
#include <nil/crypto3/algebra/pairing/alt_bn128.hpp>
#include <nil/crypto3/algebra/pairing/bls12.hpp>
#include <nil/crypto3/algebra/pairing/mnt4.hpp>
#include <nil/crypto3/algebra/pairing/mnt6.hpp>

#include <nil/crypto3/hash/keccak.hpp>

#include <nil/crypto3/zk/snark/arithmetization/plonk/table_description.hpp>
#include <nil/crypto3/zk/snark/systems/plonk/placeholder/params.hpp>
#include <nil/crypto3/zk/snark/systems/plonk/placeholder/preprocessor.hpp>
#include <nil/crypto3/zk/snark/systems/plonk/placeholder/prover.hpp>
#include <nil/crypto3/zk/snark/systems/plonk/placeholder/verifier.hpp>

#include <nil/crypto3/marshalling/zk/types/commitments/kzg.hpp>
#include <nil/crypto3/marshalling/zk/types/placeholder/proof.hpp>

#include "./detail/circuits.hpp"
#include "./detail/random_test_initializer.hpp"

using namespace nil;
using namespace nil::crypto3;
using namespace nil::crypto3::zk;
using namespace nil::crypto3::zk::snark;

template<typename curve_type, typename merkle_hash_type, typename transcript_hash_type>
struct placeholder_kzg_v2_proof_test_runner
    : public test_tools::random_test_initializer<typename curve_type::scalar_field_type> {
    using field_type = typename curve_type::scalar_field_type;

    using transcript_type = typename transcript::fiat_shamir_heuristic_sequential<transcript_hash_type>;

    using circuit_params = placeholder_circuit_params<field_type>;

    using kzg_type = commitments::batched_kzg<curve_type, transcript_hash_type>;
    using kzg_scheme_type = typename commitments::kzg_commitment_scheme_v2<kzg_type>;
    using kzg_placeholder_params_type = nil::crypto3::zk::snark::placeholder_params<circuit_params, kzg_scheme_type>;

    using policy_type = zk::snark::detail::placeholder_policy<field_type, kzg_placeholder_params_type>;

    using circuit_type = circuit_description<field_type, placeholder_circuit_params<field_type>>;

    placeholder_kzg_v2_proof_test_runner(circuit_type const& circuit) : circuit(circuit) {
    }

    using PlaceholderParams = kzg_placeholder_params_type;
    using ProofType = placeholder_proof<typename PlaceholderParams::field_type, PlaceholderParams>;
    using CommitmentParamsType = typename kzg_type::params_type;

    void test_placeholder_proof(const ProofType& proof, const CommitmentParamsType& params,
                                std::string output_file = "") {
        using namespace nil::crypto3::marshalling;
        using Endianness = nil::marshalling::option::big_endian;
        using TTypeBase = nil::marshalling::field_type<Endianness>;

        using proof_marshalling_type = nil::crypto3::marshalling::types::placeholder_proof<TTypeBase, ProofType>;

        auto filled_placeholder_proof = types::fill_placeholder_proof<Endianness, ProofType>(proof, params);
        ProofType _proof = types::make_placeholder_proof<Endianness, ProofType>(filled_placeholder_proof);
        BOOST_CHECK(_proof == proof);

        std::vector<std::uint8_t> cv;
        cv.resize(filled_placeholder_proof.length(), 0x00);
        auto write_iter = cv.begin();
        auto status = filled_placeholder_proof.write(write_iter, cv.size());
        BOOST_CHECK(status == nil::marshalling::status_type::success);

        proof_marshalling_type test_val_read;
        auto read_iter = cv.begin();
        status = test_val_read.read(read_iter, cv.size());
        BOOST_CHECK(status == nil::marshalling::status_type::success);
        auto constructed_val_read = types::make_placeholder_proof<Endianness, ProofType>(test_val_read);
        BOOST_CHECK(proof == constructed_val_read);
    }

    bool run_test() {
        std::size_t table_rows_log = std::log2(circuit.table_rows);

        typename policy_type::constraint_system_type constraint_system(circuit.gates, circuit.copy_constraints,
                                                                       circuit.lookup_gates);
        typename policy_type::variable_assignment_type assignments = circuit.table;

        // KZG commitment scheme
        typename kzg_type::field_type::value_type alpha(7u);
        auto kzg_params = kzg_scheme_type::create_params(1 << table_rows_log, alpha);
        kzg_scheme_type kzg_scheme(kzg_params);

        plonk_table_description<field_type> desc = circuit.table.get_description();
        desc.usable_rows_amount = circuit.usable_rows;

        typename placeholder_public_preprocessor<field_type, kzg_placeholder_params_type>::preprocessed_data_type
            kzg_preprocessed_public_data =
                placeholder_public_preprocessor<field_type, kzg_placeholder_params_type>::process(
                    constraint_system, assignments.public_table(), desc, kzg_scheme);

        typename placeholder_private_preprocessor<field_type, kzg_placeholder_params_type>::preprocessed_data_type
            kzg_preprocessed_private_data =
                placeholder_private_preprocessor<field_type, kzg_placeholder_params_type>::process(
                    constraint_system, assignments.private_table(), desc);

        auto kzg_proof = placeholder_prover<field_type, kzg_placeholder_params_type>::process(
            kzg_preprocessed_public_data, std::move(kzg_preprocessed_private_data), desc, constraint_system,
            kzg_scheme);

        using common_data_type = typename placeholder_public_preprocessor<
            field_type, kzg_placeholder_params_type>::preprocessed_data_type::common_data_type;
        using Endianness = nil::marshalling::option::big_endian;

        test_placeholder_proof(kzg_proof, kzg_params);

        kzg_scheme = kzg_scheme_type(kzg_params);
        bool verifier_res = placeholder_verifier<field_type, kzg_placeholder_params_type>::process(
            *kzg_preprocessed_public_data.common_data, kzg_proof, desc, constraint_system, kzg_scheme);
        BOOST_CHECK(verifier_res);
        return true;
    }

    circuit_type circuit;
};

BOOST_AUTO_TEST_SUITE(placeholder_kzg_v2_proof)

using keccak_256 = hashes::keccak_1600<256>;

using TestRunners =
    boost::mpl::list<placeholder_kzg_v2_proof_test_runner<curves::alt_bn128_254, keccak_256, keccak_256>>;

template<typename TestRunner>
using proof_marshalling_type = nil::crypto3::marshalling::types::placeholder_proof<
    nil::marshalling::field_type<nil::marshalling::option::big_endian>, typename TestRunner::ProofType>;

static_assert(std::is_default_constructible_v<
              proof_marshalling_type<placeholder_kzg_v2_proof_test_runner<curves::bls12_381, keccak_256, keccak_256>>>);
static_assert(std::is_default_constructible_v<
              proof_marshalling_type<placeholder_kzg_v2_proof_test_runner<curves::mnt4_298, keccak_256, keccak_256>>>);
static_assert(std::is_default_constructible_v<
              proof_marshalling_type<placeholder_kzg_v2_proof_test_runner<curves::mnt6_298, keccak_256, keccak_256>>>);

using canonical_kzg_test_runner = placeholder_kzg_v2_proof_test_runner<curves::alt_bn128_254, keccak_256, keccak_256>;

BOOST_AUTO_TEST_CASE_TEMPLATE(circuit_1, TestRunner, TestRunners) {
    using field_type = typename TestRunner::field_type;
    test_tools::random_test_initializer<field_type> random_test_initializer;
    auto circuit =
        circuit_test_1<field_type>(random_test_initializer.alg_random_engines.template get_alg_engine<field_type>(),
                                   random_test_initializer.generic_random_engine);

    TestRunner test_runner(circuit);

    BOOST_CHECK(test_runner.run_test());
}

BOOST_AUTO_TEST_CASE(circuit_2) {
    using field_type = typename canonical_kzg_test_runner::field_type;
    test_tools::random_test_initializer<field_type> random_test_initializer;
    auto pi0 = random_test_initializer.alg_random_engines.template get_alg_engine<field_type>()();
    auto circuit =
        circuit_test_t<field_type>(pi0,
                                   random_test_initializer.alg_random_engines.template get_alg_engine<field_type>(),
                                   random_test_initializer.generic_random_engine);

    canonical_kzg_test_runner test_runner(circuit);

    BOOST_CHECK(test_runner.run_test());
}

BOOST_AUTO_TEST_CASE(circuit_5) {
    using field_type = typename canonical_kzg_test_runner::field_type;
    test_tools::random_test_initializer<field_type> random_test_initializer;
    auto circuit =
        circuit_test_5<field_type>(random_test_initializer.alg_random_engines.template get_alg_engine<field_type>(),
                                   random_test_initializer.generic_random_engine);

    canonical_kzg_test_runner test_runner(circuit);

    BOOST_CHECK(test_runner.run_test());
}

BOOST_AUTO_TEST_CASE(circuit_fib) {
    using field_type = typename canonical_kzg_test_runner::field_type;
    test_tools::random_test_initializer<field_type> random_test_initializer;
    auto circuit = circuit_test_fib<field_type, 32>(
        random_test_initializer.alg_random_engines.template get_alg_engine<field_type>());

    canonical_kzg_test_runner test_runner(circuit);

    BOOST_CHECK(test_runner.run_test());
}

BOOST_AUTO_TEST_SUITE_END()
