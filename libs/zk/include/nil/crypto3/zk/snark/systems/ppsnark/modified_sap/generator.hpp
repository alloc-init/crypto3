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

#ifndef CRYPTO3_ZK_MODIFIED_SAP_GENERATOR_HPP
#define CRYPTO3_ZK_MODIFIED_SAP_GENERATOR_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/random/uniform_int_distribution.hpp>

#include <nil/crypto3/algebra/algorithms/pair.hpp>
#include <nil/crypto3/algebra/multiexp/multiexp.hpp>
#include <nil/crypto3/zk/snark/reductions/modified_sap_to_polynomials.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/policy.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {
                namespace detail {

                    template<typename ScalarValueType>
                    struct modified_sap_setup_trapdoor {
                        ScalarValueType tau;
                        ScalarValueType gamma;
                        // A, C, H, Z order.
                        std::array<ScalarValueType, 4> alpha;
                    };

                    /**
                     * Deterministic core of circuit-specific modified SAP setup.
                     * Random sampling and circuit-digest construction are caller responsibilities.
                     */
                    template<typename Policy>
                    class modified_sap_deterministic_generator {
                        using policy_type = Policy;
                        using curve_type = typename policy_type::curve_type;
                        using pairing_policy_type = typename policy_type::pairing_policy_type;
                        using scalar_field_type = typename policy_type::scalar_field_type;
                        using scalar_value_type = typename scalar_field_type::value_type;
                        using base_value_type = typename policy_type::base_field_type::value_type;
                        using g1_type = typename policy_type::g1_type;
                        using g1_value_type = typename g1_type::value_type;
                        using g2_value_type = typename policy_type::g2_type::value_type;
                        using reduction_type = reductions::modified_sap_to_polynomials<scalar_field_type>;

                        static constexpr std::size_t alpha_a_index = 0;
                        static constexpr std::size_t alpha_c_index = 1;
                        static constexpr std::size_t alpha_h_index = 2;
                        static constexpr std::size_t alpha_z_index = 3;

                        static std::vector<g1_value_type> encode_g1(const std::vector<scalar_value_type> &scalars,
                                                                    const algebra::window_table<g1_type> &window_table,
                                                                    std::size_t window_size) {
                            return algebra::batch_exp<g1_type, scalar_field_type>(
                                scalar_field_type::value_bits, window_size, window_table, scalars);
                        }

                    public:
                        using constraint_system_type = typename policy_type::constraint_system_type;
                        using keypair_type = typename policy_type::keypair_type;
                        using trapdoor_type = modified_sap_setup_trapdoor<scalar_value_type>;

                        static keypair_type process(const constraint_system_type &constraint_system,
                                                    const base_value_type &circuit_digest,
                                                    const trapdoor_type &trapdoor) {
                            const auto canonical_constraint_system = constraint_system.normalized();
                            if (trapdoor.gamma.is_zero()) {
                                throw std::invalid_argument("modified_sap: setup gamma must be nonzero");
                            }
                            if (trapdoor.tau.is_zero()) {
                                throw std::invalid_argument("modified_sap: setup tau must be nonzero");
                            }

                            const auto evaluation =
                                reduction_type::instance_map_with_evaluation(canonical_constraint_system, trapdoor.tau);
                            if (evaluation.Zt.is_zero()) {
                                throw std::invalid_argument("modified_sap: setup tau must lie outside the domain");
                            }
                            const std::size_t n = canonical_constraint_system.num_variables();
                            const std::size_t m =
                                reduction_type::get_domain_size(canonical_constraint_system.num_constraints());

                            std::vector<scalar_value_type> tau_powers(m, scalar_value_type::one());
                            for (std::size_t i = 1; i < m; ++i) {
                                tau_powers[i] = tau_powers[i - 1] * trapdoor.tau;
                            }

                            std::vector<scalar_value_type> w_scalars(n);
                            for (std::size_t i = 0; i < n; ++i) {
                                w_scalars[i] = trapdoor.gamma * (trapdoor.alpha[alpha_a_index] * evaluation.At[i] +
                                                                 trapdoor.alpha[alpha_c_index] * evaluation.Ct[i]);
                            }

                            const auto scaled_powers = [&tau_powers](const scalar_value_type &scale,
                                                                     std::size_t count) {
                                std::vector<scalar_value_type> result(count);
                                for (std::size_t i = 0; i < count; ++i) {
                                    result[i] = scale * tau_powers[i];
                                }
                                return result;
                            };

                            const auto h_scalars = scaled_powers(trapdoor.gamma * trapdoor.alpha[alpha_h_index], m - 1);
                            const auto t_a_scalars = scaled_powers(trapdoor.alpha[alpha_a_index], m - 1);
                            const auto t_c_scalars = scaled_powers(trapdoor.alpha[alpha_c_index], m - 1);
                            const auto t_h_scalars = scaled_powers(trapdoor.alpha[alpha_h_index], m - 2);
                            const auto t_z_scalars = scaled_powers(trapdoor.alpha[alpha_z_index], m);

                            const g1_value_type g1_generator = g1_value_type::one();
                            const std::size_t window_size = algebra::get_exp_window_size<g1_type>(std::max(n, m));
                            const auto window_table = algebra::get_window_table<g1_type>(
                                scalar_field_type::value_bits, window_size, g1_generator);

                            typename policy_type::proving_key_type proving_key;
                            proving_key.W = encode_g1(w_scalars, window_table, window_size);
                            proving_key.H_query = encode_g1(h_scalars, window_table, window_size);
                            proving_key.T_A = encode_g1(t_a_scalars, window_table, window_size);
                            proving_key.T_C = encode_g1(t_c_scalars, window_table, window_size);
                            proving_key.T_H = encode_g1(t_h_scalars, window_table, window_size);
                            proving_key.T_Z = encode_g1(t_z_scalars, window_table, window_size);
                            proving_key.constraint_system = canonical_constraint_system;

                            const g2_value_type g2_generator = g2_value_type::one();
                            const auto gt_generator =
                                algebra::pair_reduced<curve_type, pairing_policy_type>(g1_generator, g2_generator);
                            if (!gt_generator) {
                                throw std::logic_error("modified_sap: pairing generators produced no GT value");
                            }

                            typename policy_type::verification_key_type verification_key;
                            verification_key.g2_one = g2_generator;
                            verification_key.tau_g2 = trapdoor.tau * g2_generator;
                            verification_key.gamma_inverse_g2 = trapdoor.gamma.inversed() * g2_generator;
                            verification_key.alpha_z_vanishing_gt =
                                gt_generator->pow((trapdoor.alpha[alpha_z_index] * evaluation.Zt).to_integral());
                            for (std::size_t i = 0; i < verification_key.alpha_gt.size(); ++i) {
                                verification_key.alpha_gt[i] = gt_generator->pow(trapdoor.alpha[i].to_integral());
                            }
                            verification_key.num_variables = n;
                            verification_key.domain_size = m;
                            verification_key.circuit_digest = circuit_digest;

                            proving_key.verification_key = verification_key;
                            return {std::move(proving_key), std::move(verification_key)};
                        }
                    };

                }    // namespace detail

                /**
                 * Randomized circuit-specific setup for the modified SAP SNARK.
                 */
                template<typename Policy>
                class modified_sap_generator {
                    using policy_type = Policy;
                    using scalar_field_type = typename policy_type::scalar_field_type;
                    using scalar_value_type = typename scalar_field_type::value_type;
                    using scalar_integral_type = typename scalar_field_type::integral_type;
                    using transcript_policy_type = typename policy_type::transcript_policy_type;
                    using reduction_type = reductions::modified_sap_to_polynomials<scalar_field_type>;
                    using deterministic_generator_type = detail::modified_sap_deterministic_generator<policy_type>;

                    template<typename RandomSource>
                    static scalar_value_type sample_scalar(RandomSource &random_source) {
                        boost::random::uniform_int_distribution<scalar_integral_type> distribution(
                            scalar_integral_type(0), scalar_field_type::modulus - 1);
                        return scalar_value_type(distribution(random_source));
                    }

                public:
                    using constraint_system_type = typename policy_type::constraint_system_type;
                    using keypair_type = typename policy_type::keypair_type;

                    template<typename RandomSource>
                    static keypair_type process(const constraint_system_type &constraint_system,
                                                RandomSource &random_source) {
                        const auto canonical_constraint_system = constraint_system.normalized();
                        // Validate the padded domain before consuming caller-owned randomness.
                        const auto domain = reduction_type::get_domain(canonical_constraint_system);
                        const auto circuit_digest = transcript_policy_type::circuit_digest(canonical_constraint_system);

                        typename deterministic_generator_type::trapdoor_type trapdoor;
                        do {
                            trapdoor.tau = sample_scalar(random_source);
                        } while (trapdoor.tau.is_zero() ||
                                 domain->compute_vanishing_polynomial(trapdoor.tau).is_zero());
                        do {
                            trapdoor.gamma = sample_scalar(random_source);
                        } while (trapdoor.gamma.is_zero());
                        for (auto &alpha : trapdoor.alpha) {
                            alpha = sample_scalar(random_source);
                        }

                        return deterministic_generator_type::process(
                            canonical_constraint_system, circuit_digest, trapdoor);
                    }
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SAP_GENERATOR_HPP
