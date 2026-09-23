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

#ifndef CRYPTO3_ZK_MODIFIED_SQAP_TRANSCRIPT_HPP
#define CRYPTO3_ZK_MODIFIED_SQAP_TRANSCRIPT_HPP

#include <cstddef>
#include <string_view>
#include <vector>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/hash/algorithm/hash.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_policy.hpp>
#include <nil/crypto3/hash/poseidon.hpp>
#include <nil/crypto3/zk/snark/reductions/modified_sqap_to_polynomials.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/verification_key.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {

                /**
                 * Standalone BN254 Poseidon1 transcript profile for the modified SQAP SNARK.
                 */
                class modified_sqap_bn254_poseidon_transcript_policy {
                    using curve_type = algebra::curves::alt_bn128_254;
                    using base_field_type = typename curve_type::base_field_type;
                    using base_value_type = typename base_field_type::value_type;
                    using base_integral_type = typename base_field_type::integral_type;
                    using scalar_field_type = typename curve_type::scalar_field_type;
                    using reduction_type = reductions::modified_sqap_to_polynomials<scalar_field_type>;
                    using poseidon_policy_type = hashes::detail::poseidon1_policy<base_field_type, 128, 2>;
                    using poseidon_type = hashes::poseidon1<poseidon_policy_type>;

                    static base_value_type tag(std::string_view text) {
                        base_integral_type value = 0;
                        for (const unsigned char byte : text) {
                            value <<= 8;
                            value += byte;
                        }
                        return base_value_type(value);
                    }

                    static base_value_type encode_size(std::size_t value) {
                        return base_value_type(value);
                    }

                    static base_value_type encode_scalar(const typename scalar_field_type::value_type &value) {
                        return base_value_type(value.to_integral());
                    }

                    using g1_value_type = typename curve_type::template g1_type<>::value_type;
                    using g2_value_type = typename curve_type::template g2_type<>::value_type;
                    using gt_value_type = typename curve_type::gt_type::value_type;

                    static void append_g1(std::vector<base_value_type> &input,
                                          const g1_value_type &point) {
                        if (point.is_zero()) {
                            input.insert(input.end(), 3, base_value_type::zero());
                            return;
                        }
                        const auto affine = point.to_affine();
                        input.push_back(base_value_type::one());
                        input.push_back(affine.X);
                        input.push_back(affine.Y);
                    }

                    static void append_g2(std::vector<base_value_type> &input,
                                          const g2_value_type &point) {
                        if (point.is_zero()) {
                            input.insert(input.end(), 5, base_value_type::zero());
                            return;
                        }
                        const auto affine = point.to_affine();
                        input.push_back(base_value_type::one());
                        input.push_back(affine.X.data[0]);
                        input.push_back(affine.X.data[1]);
                        input.push_back(affine.Y.data[0]);
                        input.push_back(affine.Y.data[1]);
                    }

                    static void append_gt(std::vector<base_value_type> &input,
                                          const gt_value_type &value) {
                        for (std::size_t outer = 0; outer < 2; ++outer) {
                            for (std::size_t middle = 0; middle < 3; ++middle) {
                                for (std::size_t inner = 0; inner < 2; ++inner) {
                                    input.push_back(value.data[outer].data[middle].data[inner]);
                                }
                            }
                        }
                    }

                public:
                    using constraint_system_type = typename reduction_type::constraint_system_type;
                    using digest_type = base_value_type;
                    using verification_key_type = modified_sqap_verification_key<curve_type>;
                    using commitment_type = g1_value_type;
                    using challenge_type = typename scalar_field_type::value_type;

                    static digest_type circuit_digest(const constraint_system_type &constraint_system) {
                        const auto canonical = constraint_system.normalized();
                        const std::size_t domain_size = reduction_type::get_domain_size(canonical.num_constraints());

                        std::vector<digest_type> input;
                        input.push_back(tag("modified-sqap-snark"));
                        input.push_back(tag("v1"));
                        input.push_back(tag("bn254"));
                        input.push_back(tag("poseidon1-fp-128-r2-c1"));
                        input.push_back(tag("circuit-digest"));
                        input.push_back(encode_size(canonical.num_variables()));
                        input.push_back(encode_size(canonical.num_constraints()));
                        input.push_back(encode_size(domain_size));

                        for (std::size_t row = 0; row < canonical.constraints.size(); ++row) {
                            const auto &constraint = canonical.constraints[row];
                            input.push_back(tag("row"));
                            input.push_back(encode_size(row));
                            input.push_back(tag("A"));
                            input.push_back(encode_size(constraint.a.terms.size()));
                            for (const auto &term : constraint.a) {
                                input.push_back(encode_size(term.index));
                                input.push_back(encode_scalar(term.coeff));
                            }
                            input.push_back(tag("C"));
                            input.push_back(encode_size(constraint.c.terms.size()));
                            for (const auto &term : constraint.c) {
                                input.push_back(encode_size(term.index));
                                input.push_back(encode_scalar(term.coeff));
                            }
                        }

                        return hash<poseidon_type>(input);
                    }

                    static challenge_type proof_challenge(const verification_key_type &verification_key,
                                                          const commitment_type &P,
                                                          const challenge_type &u) {
                        std::vector<digest_type> input;
                        input.push_back(tag("modified-sqap-snark"));
                        input.push_back(tag("v1"));
                        input.push_back(tag("bn254"));
                        input.push_back(tag("poseidon1-fp-128-r2-c1"));
                        input.push_back(tag("pairing-exact-e"));
                        input.push_back(tag("proof-challenge"));
                        input.push_back(verification_key.circuit_digest);
                        input.push_back(encode_size(verification_key.num_variables));
                        input.push_back(encode_size(verification_key.domain_size));

                        input.push_back(encode_size(3));
                        append_g2(input, verification_key.g2_one);
                        append_g2(input, verification_key.tau_g2);
                        append_g2(input, verification_key.gamma_inverse_g2);

                        input.push_back(encode_size(5));
                        append_gt(input, verification_key.alpha_z_vanishing_gt);
                        for (const auto &alpha : verification_key.alpha_gt) {
                            append_gt(input, alpha);
                        }

                        append_g1(input, P);
                        input.push_back(encode_scalar(u));

                        const digest_type digest = hash<poseidon_type>(input);
                        // Constructing an Fr element reduces the canonical Fp representative modulo r.
                        return challenge_type(digest.to_integral());
                    }
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SQAP_TRANSCRIPT_HPP
