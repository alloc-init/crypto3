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

#ifndef CRYPTO3_ZK_MODIFIED_SAP_POLICY_HPP
#define CRYPTO3_ZK_MODIFIED_SAP_POLICY_HPP

#include <utility>
#include <vector>

#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/proof.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/proving_key.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/verification_key.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {

                /**
                 * Compile-time types for the modified SAP SNARK.
                 * PairingPolicy must use the exact final exponent (p^12 - 1) / r.
                 * TranscriptPolicy specifies the protocol transcript.
                 */
                template<typename CurveType, typename PairingPolicy, typename TranscriptPolicy>
                struct modified_sap_policy {
                    using curve_type = CurveType;
                    using base_field_type = typename curve_type::base_field_type;
                    using scalar_field_type = typename curve_type::scalar_field_type;
                    using g1_type = typename curve_type::template g1_type<>;
                    using g2_type = typename curve_type::template g2_type<>;
                    using gt_type = typename curve_type::gt_type;

                    using pairing_policy_type = PairingPolicy;
                    using transcript_policy_type = TranscriptPolicy;

                    using primary_input_type = typename scalar_field_type::value_type;
                    // The full explicitly indexed witness, including u at index zero.
                    using auxiliary_input_type = std::vector<primary_input_type>;
                    using constraint_system_type = modified_sap_constraint_system<scalar_field_type>;
                    using proof_type = modified_sap_proof<curve_type>;
                    using proving_key_type = modified_sap_proving_key<curve_type, constraint_system_type>;
                    using verification_key_type = modified_sap_verification_key<curve_type>;
                    using keypair_type = std::pair<proving_key_type, verification_key_type>;
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SAP_POLICY_HPP
