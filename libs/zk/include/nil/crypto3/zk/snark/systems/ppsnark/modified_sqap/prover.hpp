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

#include <array>
#include <cstddef>
#include <stdexcept>

#include <nil/crypto3/algebra/multiexp/multiexp.hpp>
#include <nil/crypto3/math/polynomial/operations/basic_operations.hpp>
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
                    using scalar_value_type = typename scalar_field_type::value_type;
                    using g1_value_type = typename policy_type::g1_type::value_type;
                    using transcript_policy_type = typename policy_type::transcript_policy_type;
                    using reduction_type = reductions::modified_sqap_to_polynomials<scalar_field_type>;
                    using polynomial_type = typename reduction_type::polynomial_type;

                public:
                    using primary_input_type = typename policy_type::primary_input_type;
                    using auxiliary_input_type = typename policy_type::auxiliary_input_type;
                    using proving_key_type = typename policy_type::proving_key_type;

                    /**
                     * Owned witness polynomials, vanishing polynomial and mixed-basis commitment.
                     */
                    struct commitment_result {
                        typename reduction_type::witness_polynomials polynomials;
                        polynomial_type Z;
                        g1_value_type P;
                    };

                    /**
                     * Evaluation U(z) and the canonical coefficient polynomial (U(X) - U(z)) / (X - z).
                     * The zero quotient has one zero coefficient.
                     */
                    struct polynomial_opening {
                        scalar_value_type evaluation;
                        polynomial_type quotient;
                    };

                    /**
                     * Evaluations in A, C, H, Z order and their mixed opening commitment.
                     */
                    struct opening_result {
                        std::array<scalar_value_type, 4> evaluations;
                        g1_value_type Q;
                    };

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

                    /**
                     * Validate the key and witness, derive A, C, H, Z, and compute
                     * P = sum_i w[i] * W[i] + sum_i h[i] * H_query[i].
                     * Invalid inputs throw std::invalid_argument; inputs are not modified.
                     */
                    static commitment_result commit(const proving_key_type &proving_key,
                                                    const primary_input_type &primary_input,
                                                    const auxiliary_input_type &auxiliary_input) {
                        validate(proving_key);
                        commitment_result result {
                            reduction_type::witness_map(proving_key.constraint_system, primary_input, auxiliary_input),
                            reduction_type::get_domain(proving_key.constraint_system)->get_vanishing_polynomial(),
                            g1_value_type::zero()};

                        // Validated ranges match in length and are nonempty: n >= 1 and m - 1 >= 1.
                        result.P = algebra::multiexp<algebra::policies::multiexp_method_BDLO12>(
                            proving_key.W.begin(), proving_key.W.end(), auxiliary_input.begin(), auxiliary_input.end(),
                            1);
                        result.P += algebra::multiexp<algebra::policies::multiexp_method_BDLO12>(
                            proving_key.H_query.begin(), proving_key.H_query.end(), result.polynomials.H.begin(),
                            result.polynomials.H.end(), 1);
                        return result;
                    }

                    /**
                     * Evaluate and construct an opening quotient for any challenge, including zero and domain points.
                     * Normalize a temporary copy; the source polynomial's coefficient storage is preserved.
                     * A nonzero division remainder signals an internal inconsistency and throws std::logic_error.
                     */
                    static polynomial_opening open_polynomial(const polynomial_type &polynomial,
                                                              const scalar_value_type &challenge) {
                        auto numerator = polynomial;
                        // Division requires canonical, nonempty coefficients; condense also maps empty storage to [0].
                        numerator.condense();
                        polynomial_opening result {numerator.evaluate(challenge), {}};
                        numerator[0] -= result.evaluation;

                        // Coefficients are in ascending powers: {-z, 1} represents X - z.
                        const polynomial_type divisor {-challenge, scalar_value_type::one()};
                        polynomial_type remainder;
                        math::division(result.quotient, remainder, numerator, divisor);
                        if (!remainder.is_zero()) {
                            throw std::logic_error("modified_sqap: opening quotient has a nonzero remainder");
                        }
                        return result;
                    }

                    /**
                     * Derive the transcript challenge and commit to the four opening quotients.
                     * Requires the validated key and matching public input and result from commit().
                     * A nonzero quotient exceeding its query length throws std::invalid_argument.
                     */
                    static opening_result open(const proving_key_type &proving_key,
                                               const primary_input_type &primary_input,
                                               const commitment_result &committed) {
                        const auto challenge = transcript_policy_type::proof_challenge(proving_key.verification_key,
                                                                                       committed.P, primary_input);
                        const std::array polynomials = {&committed.polynomials.A, &committed.polynomials.C,
                                                        &committed.polynomials.H, &committed.Z};
                        const std::array queries = {&proving_key.T_A, &proving_key.T_C, &proving_key.T_H,
                                                    &proving_key.T_Z};

                        opening_result result {{}, g1_value_type::zero()};
                        for (std::size_t i = 0; i < polynomials.size(); ++i) {
                            const auto opening = open_polynomial(*polynomials[i], challenge);
                            result.evaluations[i] = opening.evaluation;
                            // The canonical zero quotient is [0]; it needs no bases, even when T_H is empty.
                            if (opening.quotient.is_zero()) {
                                continue;
                            }
                            const auto &query = *queries[i];
                            if (opening.quotient.size() > query.size()) {
                                throw std::invalid_argument("modified_sqap: opening quotient exceeds query length");
                            }
                            result.Q += algebra::multiexp<algebra::policies::multiexp_method_BDLO12>(
                                query.begin(), query.begin() + opening.quotient.size(), opening.quotient.begin(),
                                opening.quotient.end(), 1);
                        }
                        return result;
                    }
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SQAP_PROVER_HPP
