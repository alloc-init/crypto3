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

#include <array>
#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/options.hpp>
#include <nil/marshalling/status_type.hpp>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/marshalling/zk/types/modified_sap/constraint_system.hpp>
#include <nil/crypto3/marshalling/zk/types/modified_sap/proof.hpp>
#include <nil/crypto3/marshalling/zk/types/modified_sap/proving_key.hpp>
#include <nil/crypto3/marshalling/zk/types/modified_sap/verification_key.hpp>
#include <nil/crypto3/math/linear_variable.hpp>
#include <nil/crypto3/random/chacha_urbg.hpp>
#include <nil/crypto3/zk/snark/reductions/r1cs_to_modified_sap.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/policy.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/transcript.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap_snark.hpp>

int main() {
    namespace snark = nil::crypto3::zk::snark;
    namespace types = nil::crypto3::marshalling::types;
    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using field_type = curve_type::scalar_field_type;
    using value_type = field_type::value_type;
    using variable_type = nil::crypto3::math::linear_variable<field_type>;
    using policy_type = snark::modified_sap_policy<curve_type,
                                                   snark::modified_sap_bn254_exact_pairing_policy,
                                                   snark::modified_sap_bn254_poseidon_transcript_policy>;
    using scheme_type = snark::modified_sap_snark<policy_type>;
    using frontend = snark::reductions::r1cs_to_modified_sap<field_type>;
    using endianness = nil::marshalling::option::big_endian;
    using type_base = nil::marshalling::field_type<endianness>;
    using status_type = nil::marshalling::status_type;

    try {
        snark::r1cs_constraint_system<field_type> source;
        source.primary_input_size = 1;
        source.auxiliary_input_size = 2;
        // Source indices: 0 = one, 1 = public u, 2 = private x, 3 = private y.
        snark::r1cs_constraint<field_type> row {variable_type(0), variable_type(3), variable_type(1)};
        row.a.add_term(variable_type(2));    // (1 + x)*y = u.
        source.add_constraint(row);

        const auto system = frontend::instance_map(source);

        // Encode the modified SAP constraint system into bytes.
        auto encoded_system =
            types::fill_modified_sap_constraint_system<scheme_type::constraint_system_type, endianness>(system);
        std::vector<std::uint8_t> system_bytes(encoded_system.length());
        auto system_output = system_bytes.begin();
        if (encoded_system.write(system_output, system_bytes.size()) != status_type::success) {
            throw std::runtime_error("Failed to serialize the constraint system");
        }

        // Decode the bytes and reconstruct the constraint system.
        types::modified_sap_constraint_system<type_base, scheme_type::constraint_system_type> decoded_system;
        auto system_input = system_bytes.cbegin();
        // Each buffer holds one complete object. Check both read status and full frame consumption before make_*.
        if (decoded_system.read(system_input, system_bytes.size()) != status_type::success ||
            system_input != system_bytes.cend()) {
            throw std::runtime_error("Failed to deserialize the constraint system");
        }
        const auto restored_system =
            types::make_modified_sap_constraint_system<scheme_type::constraint_system_type, endianness>(decoded_system);

        // Fixed seed for a reproducible demonstration only: these setup keys are insecure.
        // A real setup must seed ChaCha from cryptographic entropy and keep the seed secret.
        const std::array<std::uint8_t, 32> demonstration_seed = {1};
        nil::crypto3::random::chacha_urbg<> random_source(demonstration_seed);
        const auto keys = scheme_type::generate(restored_system, random_source);

        // Encode the proving key into bytes.
        auto encoded_pk = types::fill_modified_sap_proving_key<policy_type, endianness>(keys.first);
        std::vector<std::uint8_t> pk_bytes(encoded_pk.length());
        auto pk_output = pk_bytes.begin();
        if (encoded_pk.write(pk_output, pk_bytes.size()) != status_type::success) {
            throw std::runtime_error("Failed to serialize the proving key");
        }

        // Decode the bytes and reconstruct the proving key.
        types::modified_sap_proving_key<type_base, scheme_type::proving_key_type> decoded_pk;
        auto pk_input = pk_bytes.cbegin();
        if (decoded_pk.read(pk_input, pk_bytes.size()) != status_type::success || pk_input != pk_bytes.cend()) {
            throw std::runtime_error("Failed to deserialize the proving key");
        }
        const auto restored_pk = types::make_modified_sap_proving_key<policy_type, endianness>(decoded_pk);

        // Encode the verification key into bytes.
        auto encoded_vk = types::fill_modified_sap_verification_key<policy_type, endianness>(keys.second);
        std::vector<std::uint8_t> vk_bytes(encoded_vk.length());
        auto vk_output = vk_bytes.begin();
        if (encoded_vk.write(vk_output, vk_bytes.size()) != status_type::success) {
            throw std::runtime_error("Failed to serialize the verification key");
        }

        // Decode the bytes and reconstruct the verification key.
        types::modified_sap_verification_key<type_base, scheme_type::verification_key_type> decoded_vk;
        auto vk_input = vk_bytes.cbegin();
        if (decoded_vk.read(vk_input, vk_bytes.size()) != status_type::success || vk_input != vk_bytes.cend()) {
            throw std::runtime_error("Failed to deserialize the verification key");
        }
        const auto restored_vk = types::make_modified_sap_verification_key<policy_type, endianness>(decoded_vk);

        std::cout << "Circuit: " << restored_system.witness_size << " witness entries, "
                  << restored_system.num_constraints() << " logical rows, domain size " << restored_vk.domain_size
                  << '\n'
                  << "Constraint system: " << system_bytes.size() << " bytes\n"
                  << "Proving key: " << pk_bytes.size() << " bytes\n"
                  << "Verification key: " << vk_bytes.size() << " bytes\n";

        // Prove (x + 1)*y = u for u = 3, x = 2, y = 1 using the restored proving key.
        const value_type u(3);
        const value_type x(2);
        const value_type y(1);
        const auto witness = frontend::witness_map(source, {u}, {x, y});
        const auto proof = scheme_type::prove(restored_pk, u, witness);

        // Encode the generated proof into bytes.
        auto encoded_proof = types::fill_modified_sap_proof<scheme_type::proof_type, endianness>(proof);
        std::vector<std::uint8_t> proof_bytes(encoded_proof.length());
        auto proof_output = proof_bytes.begin();
        if (encoded_proof.write(proof_output, proof_bytes.size()) != status_type::success) {
            throw std::runtime_error("Failed to serialize the proof");
        }

        // Decode the bytes and reconstruct the proof for verification.
        types::modified_sap_proof<type_base, scheme_type::proof_type> decoded_proof;
        auto proof_input = proof_bytes.cbegin();
        if (decoded_proof.read(proof_input, proof_bytes.size()) != status_type::success ||
            proof_input != proof_bytes.cend()) {
            throw std::runtime_error("Failed to deserialize the proof");
        }
        const auto restored_proof = types::make_modified_sap_proof<scheme_type::proof_type, endianness>(decoded_proof);

        // Verify with the restored verification key and the caller's expected u, which is outside the proof bytes.
        const bool valid = scheme_type::verify(restored_vk, u, restored_proof);
        std::cout << "Proof: " << proof_bytes.size() << " bytes\n"
                  << "Verification: " << (valid ? "passed" : "failed") << '\n';
        return valid ? 0 : 1;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
