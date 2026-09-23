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

#ifndef CRYPTO3_ZK_MODIFIED_SQAP_SNARK_HPP
#define CRYPTO3_ZK_MODIFIED_SQAP_SNARK_HPP

#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/generator.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {

                /**
                 * Preprocessing SNARK for A(X)^2 - C(X) = Z(X)*H(X) + u.
                 * Policy selects the curve, pairing convention and transcript.
                 */
                template<typename Policy>
                class modified_sqap_snark {
                public:
                    using policy_type = Policy;
                    using constraint_system_type = typename policy_type::constraint_system_type;
                    using primary_input_type = typename policy_type::primary_input_type;
                    using auxiliary_input_type = typename policy_type::auxiliary_input_type;
                    using proving_key_type = typename policy_type::proving_key_type;
                    using verification_key_type = typename policy_type::verification_key_type;
                    using keypair_type = typename policy_type::keypair_type;
                    using proof_type = typename policy_type::proof_type;

                    // The randomness source is consumed by reference and is not retained.
                    template<typename RandomSource>
                    static keypair_type generate(const constraint_system_type &constraint_system,
                                                 RandomSource &random_source) {
                        return modified_sqap_generator<policy_type>::process(constraint_system, random_source);
                    }

                    // auxiliary_input is the full witness, including primary_input at index zero.
                    static proof_type prove(const proving_key_type &pk,
                                            const primary_input_type &primary_input,
                                            const auxiliary_input_type &auxiliary_input);

                    // primary_input is the caller's expected public input.
                    static bool verify(const verification_key_type &vk,
                                       const primary_input_type &primary_input,
                                       const proof_type &proof);
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SQAP_SNARK_HPP
