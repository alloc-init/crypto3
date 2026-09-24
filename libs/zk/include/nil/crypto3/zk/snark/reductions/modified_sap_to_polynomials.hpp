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

#ifndef CRYPTO3_ZK_MODIFIED_SAP_TO_POLYNOMIALS_HPP
#define CRYPTO3_ZK_MODIFIED_SAP_TO_POLYNOMIALS_HPP

#include <algorithm>
#include <bit>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/fields/params.hpp>
#include <nil/crypto3/math/domains/basic_radix2_domain.hpp>
#include <nil/crypto3/math/polynomial/polynomial.hpp>
#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/modified_sap.hpp>
#include <nil/crypto3/zk/snark/reductions/detail/sap_quotient.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {
                namespace reductions {

                    /**
                     * Reduction of modified squaring constraints to polynomials on a radix-two domain.
                     */
                    template<typename FieldType>
                    class modified_sap_to_polynomials {
                    public:
                        using field_type = FieldType;
                        using field_value_type = typename field_type::value_type;
                        using constraint_system_type = modified_sap_constraint_system<field_type>;
                        using constraint_type = typename constraint_system_type::constraint_type;
                        using domain_type = math::basic_radix2_domain<field_type>;
                        using polynomial_type = math::polynomial<field_value_type>;

                        /**
                         * Circuit basis evaluations At[i] = A^i(t), Ct[i] = C^i(t), and Zt = Z(t).
                         * Both vectors have one entry per witness slot (the source constraint system's witness_size),
                         * including slot zero.
                         */
                        struct instance_evaluation {
                            std::vector<field_value_type> At;
                            std::vector<field_value_type> Ct;
                            field_value_type Zt;
                        };

                        /**
                         * Witness polynomials with coefficients in ascending powers.
                         * A and C retain the domain size in coefficients; H retains one fewer, including trailing
                         * zeros.
                         */
                        struct witness_polynomials {
                            polynomial_type A;
                            polynomial_type C;
                            polynomial_type H;
                        };

                        /**
                         * Smallest supported power of two >= max(2, logical_rows), without allocation.
                         * Reject zero rows, unsupported roots and unrepresentable storage sizes.
                         */
                        static std::size_t get_domain_size(std::size_t logical_rows) {
                            if (logical_rows == 0) {
                                throw std::invalid_argument("modified_sap: expected at least one logical row");
                            }

                            const std::size_t log_m = std::max<std::size_t>(1, std::bit_width(logical_rows - 1));
                            if (log_m >= std::numeric_limits<std::size_t>::digits) {
                                throw std::invalid_argument("modified_sap: domain size overflow");
                            }
                            // Existing unity-root and radix-two FFT helpers use shifts of 1u.
                            if (log_m > algebra::fields::arithmetic_params<field_type>::s ||
                                log_m >= std::numeric_limits<unsigned int>::digits) {
                                throw std::invalid_argument("modified_sap: unsupported radix-two domain size");
                            }

                            const std::size_t m = std::size_t(1) << log_m;
                            // Z needs m + 1 coefficients; padding needs m constraint rows.
                            if (m >= std::vector<field_value_type>().max_size() ||
                                m > std::vector<constraint_type>().max_size()) {
                                throw std::invalid_argument("modified_sap: domain storage size overflow");
                            }
                            return m;
                        }

                        /**
                         * Construct the domain in natural order (1, omega, ..., omega^(m-1)).
                         * Its existing vanishing-polynomial methods provide Z(X) = X^m - 1.
                         */
                        static std::shared_ptr<domain_type> get_domain(const constraint_system_type &cs) {
                            const std::size_t m = get_domain_size(cs.num_constraints());
                            if (!cs.is_valid()) {
                                throw std::invalid_argument("modified_sap: invalid constraint system");
                            }
                            return std::make_shared<domain_type>(m);
                        }

                        /**
                         * Return canonical rows in domain order, appending copies of the binding row.
                         * The caller's logical constraint system is not modified.
                         */
                        static std::vector<constraint_type> get_padded_constraints(const constraint_system_type &cs) {
                            const std::size_t m = get_domain_size(cs.num_constraints());
                            auto normalized = cs.normalized();
                            const auto binding = normalized.constraints.front();
                            normalized.constraints.resize(m, binding);
                            return std::move(normalized.constraints);
                        }

                        /**
                         * Evaluate the circuit basis polynomials using sparse Lagrange accumulation.
                         * The point t may be any field value, including zero and domain elements.
                         * Invalid circuit structure or unsupported dimensions throw std::invalid_argument.
                         */
                        static instance_evaluation instance_map_with_evaluation(const constraint_system_type &cs,
                                                                                const field_value_type &t) {
                            const std::size_t n = cs.num_variables();
                            if (n > std::vector<field_value_type>().max_size()) {
                                throw std::invalid_argument("modified_sap: basis evaluation size overflow");
                            }
                            const auto domain = get_domain(cs);
                            instance_evaluation result {std::vector<field_value_type>(n, field_value_type::zero()),
                                                        std::vector<field_value_type>(n, field_value_type::zero()),
                                                        field_value_type::zero()};
                            // TODO: Optimize the radix-two Lagrange evaluator with math::batch_inverse_nonzero()
                            // and cached weights, preserving this domain and Z(X) = X^m - 1. Benchmark the change.
                            const auto lagrange = domain->evaluate_all_lagrange_polynomials(t, result.Zt);

                            for (std::size_t j = 0; j < cs.num_constraints(); ++j) {
                                for (const auto &term : cs.constraints[j].a) {
                                    result.At[term.index] += lagrange[j] * term.coeff;
                                }
                                for (const auto &term : cs.constraints[j].c) {
                                    result.Ct[term.index] += lagrange[j] * term.coeff;
                                }
                            }
                            // Each padding row has a = 0 and c = -w[0].
                            for (std::size_t j = cs.num_constraints(); j < domain->size(); ++j) {
                                result.Ct[0] -= lagrange[j];
                            }
                            return result;
                        }

                        /**
                         * Interpolate A and C from a satisfying witness using the domain's inverse FFT.
                         * Compute H = (A^2 - C - u) / (X^m - 1) on a disjoint coset, where m is the domain size.
                         * The canonical integer representative of u must be odd, and witness[0] must equal u.
                         * Invalid circuit structure, dimensions, public input or witness throw std::invalid_argument.
                         * An unsupported coset or a quotient exceeding degree m - 2 also throws std::invalid_argument.
                         * The caller's constraint system and witness are not modified.
                         */
                        static witness_polynomials witness_map(const constraint_system_type &cs,
                                                               const field_value_type &u,
                                                               const std::vector<field_value_type> &witness) {
                            if (!cs.is_satisfied(u, witness)) {
                                throw std::invalid_argument("modified_sap: invalid or unsatisfied witness");
                            }
                            const auto domain = get_domain(cs);
                            const std::size_t m = domain->size();
                            // divide_by_z_on_coset uses this same multiplicative generator.
                            const field_value_type coset(
                                algebra::fields::arithmetic_params<field_type>::multiplicative_generator);
                            if (coset.is_zero() || domain->compute_vanishing_polynomial(coset).is_zero()) {
                                throw std::invalid_argument("modified_sap: unsupported quotient coset");
                            }
                            // Padding copies the binding row, whose evaluations are a = 0 and c = -u.
                            witness_polynomials result {
                                polynomial_type(m, field_value_type::zero()), polynomial_type(m, -u), {}};
                            for (std::size_t j = 0; j < cs.num_constraints(); ++j) {
                                result.A[j] = cs.constraints[j].a.evaluate(witness);
                                result.C[j] = cs.constraints[j].c.evaluate(witness);
                            }
                            domain->inverse_fft(result.A.get_storage());
                            domain->inverse_fft(result.C.get_storage());

                            // The row checks ensure exact divisibility, so H has degree at most m - 2.
                            // The shared quotient routine consumes coefficients of A and C + u.
                            result.H = result.A;
                            auto c_for_quotient = result.C;
                            c_for_quotient[0] += u;
                            detail::compute_sap_quotient(*domain, result.H.get_storage(), c_for_quotient.get_storage());
                            if (!result.H[m - 1].is_zero()) {
                                throw std::invalid_argument("modified_sap: quotient exceeds degree bound");
                            }
                            result.H.resize(m - 1);
                            return result;
                        }
                    };

                }    // namespace reductions
            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SAP_TO_POLYNOMIALS_HPP
