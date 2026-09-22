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

#ifndef CRYPTO3_ZK_MODIFIED_SQAP_CONSTRAINT_SYSTEM_HPP
#define CRYPTO3_ZK_MODIFIED_SQAP_CONSTRAINT_SYSTEM_HPP

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <nil/crypto3/math/linear_combination.hpp>
#include <nil/crypto3/math/linear_variable.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {

                /**
                 * A constraint a(w)^2 - c(w) = u.
                 * Every term index directly addresses the witness, including u at index zero.
                 */
                template<typename FieldType>
                struct modified_sqap_constraint {
                    using field_type = FieldType;
                    using variable_type = math::linear_variable<field_type>;
                    using linear_combination_type =
                        math::linear_combination<variable_type, math::assignment_layout::explicit_constant>;

                    linear_combination_type a;
                    linear_combination_type c;

                    bool operator==(const modified_sqap_constraint &other) const = default;
                };

                /**
                 * Ordered logical constraints with explicitly indexed witness entries.
                 * Default construction does not produce a valid constraint system.
                 */
                template<typename FieldType>
                struct modified_sqap_constraint_system {
                    using field_type = FieldType;
                    using field_value_type = typename field_type::value_type;
                    using constraint_type = modified_sqap_constraint<field_type>;
                    using linear_combination_type = typename constraint_type::linear_combination_type;

                    // Total number of witness entries, including the public input at index zero.
                    std::size_t witness_size = 0;
                    // Logical rows include the binding row and exclude domain padding.
                    std::vector<constraint_type> constraints;

                    std::size_t num_variables() const {
                        return witness_size;
                    }

                    std::size_t num_constraints() const {
                        return constraints.size();
                    }

                    /**
                     * Check dimensions, every supplied term index, and the normalized binding row.
                     * Unsorted, repeated and zero terms are allowed when their indices are in range.
                     */
                    bool is_valid() const {
                        if (witness_size == 0 || constraints.empty()) {
                            return false;
                        }

                        // Check before merging or dropping terms, including zero and cancelling terms.
                        const auto in_range = [this](const auto &term) { return term.index < witness_size; };
                        for (const auto &constraint : constraints) {
                            if (!std::all_of(constraint.a.begin(), constraint.a.end(), in_range) ||
                                !std::all_of(constraint.c.begin(), constraint.c.end(), in_range)) {
                                return false;
                            }
                        }

                        const auto a = normalize_linear_combination(constraints.front().a);
                        const auto c = normalize_linear_combination(constraints.front().c);
                        return a.terms.empty() && c.terms.size() == 1 && c.terms.front().index == 0 &&
                               c.terms.front().coeff == -field_value_type::one();
                    }

                    /**
                     * Return an owned canonical copy, preserving row order and witness dimension.
                     * Invalid structure throws std::invalid_argument; the source is not modified.
                     */
                    modified_sqap_constraint_system normalized() const {
                        if (!is_valid()) {
                            throw std::invalid_argument("modified_sqap: invalid dimensions, term index or binding row");
                        }

                        modified_sqap_constraint_system result;
                        result.witness_size = witness_size;
                        result.constraints.reserve(constraints.size());
                        for (const auto &constraint : constraints) {
                            result.constraints.push_back({normalize_linear_combination(constraint.a),
                                                          normalize_linear_combination(constraint.c)});
                        }
                        return result;
                    }

                    bool operator==(const modified_sqap_constraint_system &other) const = default;

                private:
                    static linear_combination_type
                        normalize_linear_combination(const linear_combination_type &combination) {
                        // The existing constructor sorts terms and combines repeated indices.
                        linear_combination_type result(combination.terms);
                        std::erase_if(result.terms, [](const auto &term) { return term.coeff.is_zero(); });
                        return result;
                    }
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SQAP_CONSTRAINT_SYSTEM_HPP
