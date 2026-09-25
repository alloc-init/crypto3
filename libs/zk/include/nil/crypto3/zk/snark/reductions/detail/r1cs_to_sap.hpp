//---------------------------------------------------------------------------//
// Copyright (c) 2018-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
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

#ifndef CRYPTO3_ZK_DETAIL_R1CS_TO_SAP_HPP
#define CRYPTO3_ZK_DETAIL_R1CS_TO_SAP_HPP

#include <cstddef>

#include <nil/crypto3/math/linear_combination.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {
                namespace reductions {
                    namespace detail {

                        template<typename FieldType>
                        typename FieldType::value_type times_four(typename FieldType::value_type x) {
                            const typename FieldType::value_type times_two = x + x;
                            return times_two + times_two;
                        }

                        /**
                         * L, R and O denote the left, right and output linear combinations, respectively.
                         * Convert L*R = O into the two square equations
                         * (L + R)^2 = 4*O + t and (L - R)^2 = t.
                         * auxiliary_index names a fresh variable t, using the caller's index convention.
                         *
                         * add_a(row, index, coefficient) and add_c(row, index, coefficient) receive
                         * contributions to rows 0 and 1. Callers accumulate contributions: repeated,
                         * unsorted and zero terms are preserved. No intermediate rows are allocated.
                         * The caller validates indices and provides an odd-characteristic field.
                         */
                        template<typename VariableType, math::assignment_layout AssignmentLayout, typename AddA,
                                 typename AddC>
                        void r1cs_to_sap_constraint(
                            const math::linear_combination<VariableType, AssignmentLayout> &left,
                            const math::linear_combination<VariableType, AssignmentLayout> &right,
                            const math::linear_combination<VariableType, AssignmentLayout> &output,
                            std::size_t auxiliary_index, AddA &&add_a, AddC &&add_c) {
                            using field_type = typename VariableType::field_type;
                            using value_type = typename field_type::value_type;

                            for (const auto &term : left) {
                                add_a(0, term.index, term.coeff);
                                add_a(1, term.index, term.coeff);
                            }
                            for (const auto &term : right) {
                                add_a(0, term.index, term.coeff);
                                add_a(1, term.index, -term.coeff);
                            }
                            for (const auto &term : output) {
                                add_c(0, term.index, times_four<field_type>(term.coeff));
                            }
                            add_c(0, auxiliary_index, value_type::one());
                            add_c(1, auxiliary_index, value_type::one());
                        }

                        /**
                         * Compute the auxiliary value t = (L(w) - R(w))^2 using the combinations'
                         * assignment layout. The caller supplies valid indices and assignment dimensions.
                         */
                        template<typename VariableType, math::assignment_layout AssignmentLayout,
                                 typename AssignmentType>
                        typename AssignmentType::value_type
                            r1cs_to_sap_auxiliary(const math::linear_combination<VariableType, AssignmentLayout> &left,
                                                  const math::linear_combination<VariableType, AssignmentLayout> &right,
                                                  const AssignmentType &assignment) {
                            const auto difference = left.evaluate(assignment) - right.evaluate(assignment);
                            return difference * difference;
                        }

                    }    // namespace detail
                }    // namespace reductions
            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_DETAIL_R1CS_TO_SAP_HPP
