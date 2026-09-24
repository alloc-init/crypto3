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

#include <array>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/multiprecision/integer.hpp>

#include <nil/crypto3/zk/snark/arithmetization/bit_comparison.hpp>
#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/modified_sap.hpp>
#include <nil/crypto3/zk/snark/reductions/detail/r1cs_to_sap.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {
                namespace reductions {
                    namespace detail {

                        /**
                         * Recover constant one from the canonical bits of the odd public input w[0].
                         * The caller supplies the binding row w[0] = u and enforces canonical(u) odd.
                         * Appended entries are bits, prefix markers, then product auxiliaries.
                         * Appended rows are Booleanity, reconstruction, then comparison products.
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

                            constexpr static std::size_t bit_count =
                                boost::multiprecision::msb(field_type::modulus) + 1;
                            constexpr static auto comparison_counts = for_each_bit_comparison(
                                field_type::modulus, bit_count, field_type::modulus, [](std::size_t) { },
                                [](std::size_t) { }, [](std::size_t, std::size_t) { });
                            constexpr static std::size_t marker_count = comparison_counts.marker_count;
                            constexpr static std::size_t product_count = comparison_counts.product_count;
                            constexpr static std::size_t witness_size = bit_count + marker_count + product_count;
                            constexpr static std::size_t constraint_count = bit_count + 1 + 2 * product_count;

                            /**
                             * Append recovery rows with the bits starting at first_bit >= 1.
                             * The caller reserves witness_size fresh entries and normalizes the complete system.
                             * Invalid offsets or unrepresentable sizes throw std::invalid_argument before mutation.
                             */
                            static void append_constraints(std::vector<constraint_type> &constraints,
                                                           std::size_t first_bit) {
                                check_first_bit(first_bit);
                                if (constraint_count > constraints.max_size() ||
                                    constraints.size() > constraints.max_size() - constraint_count) {
                                    throw std::invalid_argument("modified_sap: constant-recovery row size overflow");
                                }
                                constraints.reserve(constraints.size() + constraint_count);
                                const auto minus_one = -value_type::one();
                                for (std::size_t bit = 0; bit < bit_count; ++bit) {
                                    const variable_type variable(first_bit + bit);
                                    constraint_type row {variable, variable};
                                    row.c.add_term(variable_type(0), minus_one);
                                    constraints.push_back(std::move(row));
                                }

                                constraint_type reconstruction;
                                auto coefficient = minus_one;
                                for (std::size_t bit = 0; bit < bit_count; ++bit) {
                                    reconstruction.c.add_term(variable_type(first_bit + bit), coefficient);
                                    coefficient += coefficient;
                                }
                                constraints.push_back(std::move(reconstruction));

                                for_each_product(first_bit,
                                                 [&](std::size_t marker, const combination_type &factor,
                                                     std::optional<std::size_t> next_marker, std::size_t auxiliary) {
                                                     const combination_type left {variable_type(marker)};
                                                     combination_type output;
                                                     if (next_marker) {
                                                         output.add_term(variable_type(*next_marker));
                                                     }
                                                     std::array<constraint_type, 2> rows;
                                                     r1cs_to_sap_constraint(
                                                         left, factor, output, auxiliary,
                                                         [&](std::size_t row, std::size_t index, const auto &value) {
                                                             rows[row].a.add_term(variable_type(index), value);
                                                         },
                                                         [&](std::size_t row, std::size_t index, const auto &value) {
                                                             rows[row].c.add_term(variable_type(index), value);
                                                         });
                                                     for (auto &row : rows) {
                                                         row.c.add_term(variable_type(0), minus_one);
                                                         constraints.push_back(std::move(row));
                                                     }
                                                 });
                            }

                            /**
                             * Append recovery values to an assignment beginning with u, preserving its prefix.
                             * Empty input, even canonical(u), or size overflow throws std::invalid_argument
                             * before mutation. The new first entry is the recovered constant-one bit.
                             */
                            static void append_witness(std::vector<value_type> &assignment) {
                                const auto first_bit = assignment.size();
                                check_first_bit(first_bit);
                                const auto integer = assignment.front().to_integral();
                                if (!boost::multiprecision::bit_test(integer, 0)) {
                                    throw std::invalid_argument("modified_sap: constant recovery requires odd u");
                                }
                                assignment.resize(first_bit + witness_size, value_type::zero());
                                for (std::size_t bit = 0; bit < bit_count; ++bit) {
                                    assignment[first_bit + bit] =
                                        value_type(boost::multiprecision::bit_test(integer, bit));
                                }
                                for_each_product(first_bit, [&](std::size_t marker, const combination_type &factor,
                                                                std::optional<std::size_t> next_marker,
                                                                std::size_t auxiliary) {
                                    const combination_type left {variable_type(marker)};
                                    if (next_marker) {
                                        assignment[*next_marker] = assignment[marker] * factor.evaluate(assignment);
                                    }
                                    assignment[auxiliary] = r1cs_to_sap_auxiliary(left, factor, assignment);
                                });
                            }

                        private:
                            static void check_first_bit(std::size_t first_bit) {
                                const auto max_size = std::vector<value_type>().max_size();
                                if (first_bit == 0 || witness_size > max_size || first_bit > max_size - witness_size) {
                                    throw std::invalid_argument("modified_sap: invalid constant-recovery witness size");
                                }
                            }

                            // Share product order and indices between constraint and witness construction.
                            template<typename Emit>
                            static void for_each_product(std::size_t first_bit, Emit &&emit) {
                                std::size_t marker = first_bit + bit_count - 1;
                                std::size_t next_marker = first_bit + bit_count;
                                std::size_t auxiliary = next_marker + marker_count;
                                // The matching width and odd modulus > 2 reuse the top bit before any product.
                                for_each_bit_comparison(
                                    field_type::modulus, bit_count, field_type::modulus,
                                    [&](std::size_t bit) { marker = first_bit + bit; },
                                    [&](std::size_t bit) {
                                        emit(marker, combination_type(variable_type(first_bit + bit)), next_marker,
                                             auxiliary++);
                                        marker = next_marker++;
                                    },
                                    [&](std::size_t begin, std::size_t end) {
                                        combination_type sum;
                                        for (std::size_t bit = end; bit > begin;) {
                                            sum.add_term(variable_type(first_bit + --bit));
                                        }
                                        emit(marker, sum, std::nullopt, auxiliary++);
                                    });
                            }
                        };

                    }    // namespace detail
                }    // namespace reductions
            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_R1CS_TO_MODIFIED_SAP_HPP
