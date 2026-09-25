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

#ifndef CRYPTO3_ZK_MODIFIED_SAP_VERIFIER_HPP
#define CRYPTO3_ZK_MODIFIED_SAP_VERIFIER_HPP

#include <bit>
#include <cstddef>

#include <nil/crypto3/algebra/algorithms/pair.hpp>
#include <nil/crypto3/algebra/curves/detail/scalar_mul.hpp>
#include <nil/crypto3/algebra/fields/params.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {

                /**
                 * Native verifier for the modified SAP SNARK.
                 */
                template<typename Policy>
                class modified_sap_verifier {
                    using policy_type = Policy;
                    using curve_type = typename policy_type::curve_type;
                    using pairing_policy_type = typename policy_type::pairing_policy_type;
                    using transcript_policy_type = typename policy_type::transcript_policy_type;
                    using scalar_field_type = typename policy_type::scalar_field_type;
                    using gt_value_type = typename policy_type::gt_type::value_type;

                    template<typename GroupValueType>
                    static bool is_valid_group_element(const GroupValueType &point) {
                        return point.is_well_formed() && algebra::curves::detail::subgroup_check(point);
                    }

                    static bool is_valid_gt_element(const gt_value_type &value) {
                        return !value.is_zero() && value.pow(scalar_field_type::modulus) == gt_value_type::one();
                    }

                public:
                    using primary_input_type = typename policy_type::primary_input_type;
                    using verification_key_type = typename policy_type::verification_key_type;
                    using proof_type = typename policy_type::proof_type;

                    static bool validate(const verification_key_type &verification_key,
                                         const primary_input_type &primary_input,
                                         const proof_type &proof) {
                        if ((primary_input.to_integral() & 1) == 0) {
                            return false;
                        }

                        constexpr std::size_t max_domain_log = algebra::fields::arithmetic_params<scalar_field_type>::s;
                        if (verification_key.num_variables == 0 || verification_key.domain_size < 2 ||
                            !std::has_single_bit(verification_key.domain_size) ||
                            std::bit_width(verification_key.domain_size) - 1 > max_domain_log) {
                            return false;
                        }

                        if (!is_valid_group_element(verification_key.g2_one) ||
                            !is_valid_group_element(verification_key.tau_g2) ||
                            !is_valid_group_element(verification_key.gamma_inverse_g2) ||
                            verification_key.g2_one != policy_type::g2_type::value_type::one() ||
                            verification_key.tau_g2.is_zero() || verification_key.gamma_inverse_g2.is_zero()) {
                            return false;
                        }

                        if (!is_valid_gt_element(verification_key.alpha_z_vanishing_gt)) {
                            return false;
                        }
                        for (const auto &alpha : verification_key.alpha_gt) {
                            if (!is_valid_gt_element(alpha)) {
                                return false;
                            }
                        }

                        return is_valid_group_element(proof.P) && is_valid_group_element(proof.Q);
                    }

                    static bool process(const verification_key_type &verification_key,
                                        const primary_input_type &primary_input,
                                        const proof_type &proof) {
                        if (!validate(verification_key, primary_input, proof)) {
                            return false;
                        }
                        if (proof.v_A.squared() - proof.v_C != proof.v_H * proof.v_Z + primary_input) {
                            return false;
                        }

                        const auto challenge =
                            transcript_policy_type::proof_challenge(verification_key, proof.P, primary_input);
                        const auto pairing_P = algebra::pair_reduced<curve_type, pairing_policy_type>(
                            proof.P, verification_key.gamma_inverse_g2);
                        const auto pairing_Q_tau =
                            algebra::pair_reduced<curve_type, pairing_policy_type>(proof.Q, verification_key.tau_g2);
                        const auto pairing_Q_one =
                            algebra::pair_reduced<curve_type, pairing_policy_type>(proof.Q, verification_key.g2_one);
                        if (!pairing_P || !pairing_Q_tau || !pairing_Q_one) {
                            return false;
                        }

                        const gt_value_type left = *pairing_P * verification_key.alpha_z_vanishing_gt;
                        gt_value_type right = *pairing_Q_tau * pairing_Q_one->pow((-challenge).to_integral());
                        right *= verification_key.alpha_gt[0].pow(proof.v_A.to_integral());
                        right *= verification_key.alpha_gt[1].pow(proof.v_C.to_integral());
                        right *= verification_key.alpha_gt[2].pow(proof.v_H.to_integral());
                        right *= verification_key.alpha_gt[3].pow(proof.v_Z.to_integral());
                        return left == right;
                    }
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SAP_VERIFIER_HPP
