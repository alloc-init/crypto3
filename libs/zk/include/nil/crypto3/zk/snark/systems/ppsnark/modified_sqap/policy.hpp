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

#ifndef CRYPTO3_ZK_MODIFIED_SQAP_POLICY_HPP
#define CRYPTO3_ZK_MODIFIED_SQAP_POLICY_HPP

#include <optional>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/curves/params/multiexp/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/pairing/alt_bn128.hpp>
#include <nil/crypto3/algebra/pairing/detail/alt_bn128/params.hpp>

#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/proof.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/proving_key.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/verification_key.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {

                namespace detail {

                    /**
                     * Convert Crypto3's optimized BN254 final exponentiation to
                     * the exact exponent (p^12 - 1) / r.
                     */
                    class modified_sqap_bn254_exact_final_exponentiation {
                        using curve_type = algebra::curves::alt_bn128_254;
                        using native_pairing_policy_type = algebra::pairing::pairing_policy<curve_type>;
                        using pairing_params_type = algebra::pairing::detail::pairing_params<curve_type>;
                        using scalar_value_type = typename curve_type::scalar_field_type::value_type;
                        using scalar_integral_type = typename curve_type::scalar_field_type::integral_type;
                        using gt_value_type = typename curve_type::gt_type::value_type;

                        static const scalar_integral_type &normalization_exponent() {
                            static const scalar_integral_type exponent = []() {
                                const scalar_value_type t(pairing_params_type::final_exponent_z);
                                const scalar_value_type c =
                                    scalar_value_type(2) * t *
                                    (scalar_value_type(6) * t.squared() + scalar_value_type(3) * t +
                                     scalar_value_type::one());
                                return c.inversed().to_integral();
                            }();
                            return exponent;
                        }

                    public:
                        static std::optional<gt_value_type> process(const gt_value_type &element) {
                            auto optimized_result = native_pairing_policy_type::final_exponentiation::process(element);
                            if (!optimized_result) {
                                return std::nullopt;
                            }
                            return optimized_result->pow(normalization_exponent());
                        }
                    };

                }    // namespace detail

                /**
                 * Crypto3's BN254 pairing with the exact final exponent.
                 * Selecting algebra::pairing::pairing_policy<curve_type> instead
                 * restores Crypto3's native optimized convention.
                 */
                struct modified_sqap_bn254_exact_pairing_policy :
                    algebra::pairing::pairing_policy<algebra::curves::alt_bn128_254> {
                    using final_exponentiation = detail::modified_sqap_bn254_exact_final_exponentiation;
                };

                /**
                 * Compile-time types for the modified SQAP SNARK.
                 * PairingPolicy fixes the pairing convention for all scheme operations.
                 * TranscriptPolicy specifies the protocol transcript.
                 */
                template<typename CurveType, typename PairingPolicy, typename TranscriptPolicy>
                struct modified_sqap_policy {
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
                    using constraint_system_type = modified_sqap_constraint_system<scalar_field_type>;
                    using proof_type = modified_sqap_proof<curve_type>;
                    using proving_key_type = modified_sqap_proving_key<curve_type, constraint_system_type>;
                    using verification_key_type = modified_sqap_verification_key<curve_type>;
                    using keypair_type = std::pair<proving_key_type, verification_key_type>;
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SQAP_POLICY_HPP
