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

#ifndef CRYPTO3_ZK_R1CS_TO_MODIFIED_SAP_HPP
#define CRYPTO3_ZK_R1CS_TO_MODIFIED_SAP_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/modified_sap.hpp>
#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/r1cs.hpp>
#include <nil/crypto3/zk/snark/reductions/detail/r1cs_to_sap.hpp>
#include <nil/crypto3/zk/snark/reductions/modified_sap_to_polynomials.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {
                namespace reductions {
                    namespace detail {

                        /**
                         * Recover constant one as e by enforcing e*w[0] = w[0].
                         * The caller supplies the binding row w[0] = u and enforces canonical(u) odd.
                         * Thus u != 0, so the two square-conversion rows force e = 1.
                         * Appended entries are e and t = (e - w[0])^2.
                         */
                        template<typename FieldType>
                        class modified_sap_constant_recovery {
                        public:
                            using field_type = FieldType;
                            using value_type = typename field_type::value_type;
                            using constraint_type = modified_sap_constraint<field_type>;
                            using variable_type = typename constraint_type::variable_type;
                            using combination_type = typename constraint_type::linear_combination_type;

                            static_assert(field_type::modulus > 2 && (field_type::modulus & 1) != 0,
                                          "modified_sap: constant recovery requires an odd prime field");

                            constexpr static std::size_t witness_size = 2;
                            constexpr static std::size_t constraint_count = 2;

                            /**
                             * Append recovery rows with e at constant_index >= 1 and t at constant_index + 1.
                             * The caller reserves witness_size fresh entries and normalizes the complete system.
                             * Invalid offsets or unrepresentable sizes throw std::invalid_argument before mutation.
                             */
                            static void append_constraints(std::vector<constraint_type> &constraints,
                                                           std::size_t constant_index) {
                                check_constant_index(constant_index);
                                if (constraint_count > constraints.max_size() ||
                                    constraints.size() > constraints.max_size() - constraint_count) {
                                    throw std::invalid_argument("modified_sap: constant-recovery row size overflow");
                                }
                                constraints.reserve(constraints.size() + constraint_count);
                                const combination_type constant {variable_type(constant_index)};
                                const combination_type input {variable_type(0)};
                                std::array<constraint_type, 2> rows;
                                r1cs_to_sap_constraint(
                                    constant, input, input, constant_index + 1,
                                    [&](std::size_t row, std::size_t index, const auto &value) {
                                        rows[row].a.add_term(variable_type(index), value);
                                    },
                                    [&](std::size_t row, std::size_t index, const auto &value) {
                                        rows[row].c.add_term(variable_type(index), value);
                                    });
                                // Subtracting the row equations gives 4*w[0]*(e - 1) = 0.
                                for (auto &row : rows) {
                                    row.c.add_term(variable_type(0), -value_type::one());
                                    constraints.push_back(std::move(row));
                                }
                            }

                            /**
                             * Append recovery values to an assignment beginning with u, preserving its prefix.
                             * Empty input, even canonical(u), or size overflow throws std::invalid_argument
                             * before mutation. Append e = 1 and t = (1 - u)^2.
                             */
                            static void append_witness(std::vector<value_type> &assignment) {
                                const auto constant_index = assignment.size();
                                check_constant_index(constant_index);
                                if ((assignment.front().to_integral() & 1) == 0) {
                                    throw std::invalid_argument("modified_sap: constant recovery requires odd u");
                                }
                                assignment.reserve(constant_index + witness_size);
                                assignment.push_back(value_type::one());
                                const combination_type constant {variable_type(constant_index)};
                                const combination_type input {variable_type(0)};
                                assignment.push_back(r1cs_to_sap_auxiliary(constant, input, assignment));
                            }

                        private:
                            static void check_constant_index(std::size_t constant_index) {
                                const auto max_size = std::vector<value_type>().max_size();
                                if (constant_index == 0 || witness_size > max_size ||
                                    constant_index > max_size - witness_size) {
                                    throw std::invalid_argument("modified_sap: invalid constant-recovery witness size");
                                }
                            }
                        };

                    }    // namespace detail

                    /**
                     * Convert ordinary R1CS with one public input into modified squaring constraints.
                     * Recover the source constant one using the nonzero public input.
                     */
                    template<typename FieldType>
                    class r1cs_to_modified_sap {
                    public:
                        using field_type = FieldType;
                        using field_value_type = typename field_type::value_type;
                        using constraint_system_type = modified_sap_constraint_system<field_type>;

                        /**
                         * Return owned, normalized logical rows: binding, constant recovery, then two rows
                         * per source multiplication in source order. No witness values or domain padding
                         * enter the conversion. The source is preserved and no references to it are retained.
                         *
                         * With N source variables excluding the implicit one and M source constraints,
                         * witness entries are the N source values, the recovery block, then M auxiliaries.
                         * Source index 0 maps to recovered one at N; source index i >= 1 maps to i - 1.
                         *
                         * Reject invalid public-input counts, raw indices, size overflows and unsupported
                         * domain dimensions with std::invalid_argument. Unsorted and repeated terms are
                         * accepted; check indices before normalization, including zero and cancelling terms.
                         */
                        static constraint_system_type instance_map(const r1cs_constraint_system<field_type> &cs) {
                            using recovery_type = detail::modified_sap_constant_recovery<field_type>;
                            using constraint_type = typename constraint_system_type::constraint_type;
                            using variable_type = typename constraint_type::variable_type;

                            constraint_system_type result;
                            result.witness_size = validate_source(cs);
                            const std::size_t source_variables = cs.num_variables();
                            const std::size_t first_auxiliary = source_variables + recovery_type::witness_size;
                            const std::size_t logical_rows =
                                1 + recovery_type::constraint_count + 2 * cs.num_constraints();
                            result.constraints.reserve(logical_rows);
                            // binding row
                            result.constraints.push_back({{}, -variable_type(0)});
                            recovery_type::append_constraints(result.constraints, source_variables);

                            const auto minus_one = -field_value_type::one();
                            for (std::size_t i = 0; i < cs.num_constraints(); ++i) {
                                const auto auxiliary = first_auxiliary + i;
                                const auto remap_index = [source_variables, auxiliary](std::size_t index) {
                                    // Source index 0 (implicit one) maps to the recovered one at source_variables;
                                    // source index j >= 1 maps to j - 1, putting the public input at w[0].
                                    // The helper also emits the new auxiliary's index, already in target indexing.
                                    // It lies beyond all validated source indices, so it is left unchanged.
                                    if (index == auxiliary) {
                                        return index;
                                    }
                                    return index == 0 ? source_variables : index - 1;
                                };
                                std::array<constraint_type, 2> rows;
                                const auto &source = cs.constraints[i];
                                detail::r1cs_to_sap_constraint(
                                    source.a, source.b, source.c, auxiliary,
                                    [&](std::size_t row, std::size_t index, const auto &coefficient) {
                                        rows[row].a.add_term(variable_type(remap_index(index)), coefficient);
                                    },
                                    [&](std::size_t row, std::size_t index, const auto &coefficient) {
                                        rows[row].c.add_term(variable_type(remap_index(index)), coefficient);
                                    });
                                for (auto &row : rows) {
                                    row.c.add_term(variable_type(0), minus_one);
                                    result.constraints.push_back(std::move(row));
                                }
                            }
                            return result.normalized();
                        }

                        /**
                         * Return an owned witness: u and the original auxiliaries, the constant-recovery
                         * block, then one (L - R)^2 auxiliary per source constraint in source order.
                         * L and R are the source row's left and right linear combinations.
                         * The source and input vectors are preserved; no references to them are retained.
                         *
                         * Apply the same source validation as instance_map before evaluating any row.
                         * Mismatched input lengths, even canonical(u), or an unsatisfied source throw
                         * std::invalid_argument in both debug and release builds.
                         */
                        static std::vector<field_value_type>
                            witness_map(const r1cs_constraint_system<field_type> &cs,
                                        const r1cs_primary_input<field_type> &primary_input,
                                        const r1cs_auxiliary_input<field_type> &auxiliary_input) {
                            const auto witness_size = validate_source(cs);
                            if (primary_input.size() != 1 || auxiliary_input.size() != cs.auxiliary_input_size) {
                                throw std::invalid_argument("modified_sap: R1CS input size mismatch");
                            }
                            if (!cs.is_satisfied(primary_input, auxiliary_input)) {
                                throw std::invalid_argument("modified_sap: unsatisfied R1CS witness");
                            }

                            std::vector<field_value_type> result;
                            result.reserve(witness_size);
                            result.insert(result.end(), primary_input.begin(), primary_input.end());
                            result.insert(result.end(), auxiliary_input.begin(), auxiliary_input.end());
                            detail::modified_sap_constant_recovery<field_type>::append_witness(result);
                            for (const auto &constraint : cs.constraints) {
                                // Source combinations retain implicit-one indexing and only access the
                                // original prefix, even after recovery values and auxiliaries are appended.
                                result.push_back(detail::r1cs_to_sap_auxiliary(constraint.a, constraint.b, result));
                            }
                            return result;
                        }

                    private:
                        // Validate the source for either map before evaluation or output allocation,
                        // and return the complete modified-SAP witness size.
                        static std::size_t validate_source(const r1cs_constraint_system<field_type> &cs) {
                            using recovery_type = detail::modified_sap_constant_recovery<field_type>;
                            using constraint_type = typename constraint_system_type::constraint_type;
                            using combination_type = typename constraint_type::linear_combination_type;

                            if (cs.primary_input_size != 1) {
                                throw std::invalid_argument("modified_sap: R1CS requires exactly one public input");
                            }
                            const auto max_witness_size = std::vector<field_value_type>().max_size();
                            if (recovery_type::witness_size >= max_witness_size ||
                                cs.auxiliary_input_size > max_witness_size - recovery_type::witness_size - 1) {
                                throw std::invalid_argument("modified_sap: R1CS witness size overflow");
                            }
                            const std::size_t source_variables = cs.num_variables();
                            const std::size_t first_auxiliary = source_variables + recovery_type::witness_size;
                            if (cs.num_constraints() > max_witness_size - first_auxiliary) {
                                throw std::invalid_argument("modified_sap: R1CS witness size overflow");
                            }

                            const auto max_rows = std::vector<constraint_type>().max_size();
                            if (recovery_type::constraint_count >= max_rows ||
                                cs.num_constraints() > (max_rows - recovery_type::constraint_count - 1) / 2) {
                                throw std::invalid_argument("modified_sap: R1CS row size overflow");
                            }
                            const std::size_t logical_rows =
                                1 + recovery_type::constraint_count + 2 * cs.num_constraints();
                            // Discard the size: this only validates domain support. The polynomial reduction
                            // creates the domain and pads the rows later.
                            modified_sap_to_polynomials<field_type>::get_domain_size(logical_rows);

                            // Ordinary R1CS validation also requires sorted, unique terms. Here only raw
                            // index validity is required; instance_map normalizes the generated rows.
                            const auto in_range = [source_variables](const auto &term) {
                                return term.index <= source_variables;
                            };
                            const auto max_terms = combination_type().terms.max_size();
                            for (const auto &constraint : cs.constraints) {
                                if (!std::all_of(constraint.a.begin(), constraint.a.end(), in_range) ||
                                    !std::all_of(constraint.b.begin(), constraint.b.end(), in_range) ||
                                    !std::all_of(constraint.c.begin(), constraint.c.end(), in_range)) {
                                    throw std::invalid_argument("modified_sap: R1CS term index out of range");
                                }
                                // Each A contains all left/right terms. The first C contains all output
                                // terms, plus the auxiliary and -w[0].
                                if (max_terms < 2 ||
                                    constraint.a.terms.size() > max_terms - constraint.b.terms.size() ||
                                    constraint.c.terms.size() > max_terms - 2) {
                                    throw std::invalid_argument("modified_sap: R1CS row term size overflow");
                                }
                            }

                            return first_auxiliary + cs.num_constraints();
                        }
                    };

                }    // namespace reductions
            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_R1CS_TO_MODIFIED_SAP_HPP
