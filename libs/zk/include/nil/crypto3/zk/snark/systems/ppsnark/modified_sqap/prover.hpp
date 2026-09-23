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

#ifndef CRYPTO3_ZK_MODIFIED_SQAP_PROVER_HPP
#define CRYPTO3_ZK_MODIFIED_SQAP_PROVER_HPP

#include <cstddef>
#include <stdexcept>

#include <nil/crypto3/zk/snark/reductions/modified_sqap_to_polynomials.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {

                /**
                 * Prover for the modified SQAP SNARK.
                 */
                template<typename Policy>
                class modified_sqap_prover {
                    using policy_type = Policy;
                    using scalar_field_type = typename policy_type::scalar_field_type;
                    using transcript_policy_type = typename policy_type::transcript_policy_type;
                    using reduction_type = reductions::modified_sqap_to_polynomials<scalar_field_type>;

                public:
                    using proving_key_type = typename policy_type::proving_key_type;

                    /**
                     * Check circuit structure, metadata, query lengths and circuit-digest consistency.
                     * Inconsistent inputs throw std::invalid_argument; the key is not modified.
                     */
                    static void validate(const proving_key_type &proving_key) {
                        const auto &constraint_system = proving_key.constraint_system;
                        if (!constraint_system.is_valid()) {
                            throw std::invalid_argument("modified_sqap: invalid proving-key constraint system");
                        }

                        const std::size_t n = constraint_system.num_variables();
                        const std::size_t m = reduction_type::get_domain_size(constraint_system.num_constraints());
                        const auto &verification_key = proving_key.verification_key;
                        if (verification_key.num_variables != n || verification_key.domain_size != m) {
                            throw std::invalid_argument("modified_sqap: inconsistent proving-key dimensions");
                        }

                        if (proving_key.W.size() != n || proving_key.H_query.size() != m - 1 ||
                            proving_key.T_A.size() != m - 1 || proving_key.T_C.size() != m - 1 ||
                            proving_key.T_H.size() != m - 2 || proving_key.T_Z.size() != m) {
                            throw std::invalid_argument("modified_sqap: inconsistent proving-key query lengths");
                        }

                        if (verification_key.circuit_digest !=
                            transcript_policy_type::circuit_digest(constraint_system)) {
                            throw std::invalid_argument("modified_sqap: proving-key circuit digest mismatch");
                        }
                    }
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SQAP_PROVER_HPP
