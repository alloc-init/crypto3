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

#ifndef CRYPTO3_ZK_MODIFIED_SQAP_TO_POLYNOMIALS_HPP
#define CRYPTO3_ZK_MODIFIED_SQAP_TO_POLYNOMIALS_HPP

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
#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/modified_sqap.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {
                namespace reductions {

                    /**
                     * Reduction of modified squaring constraints to polynomials on a radix-two domain.
                     */
                    template<typename FieldType>
                    class modified_sqap_to_polynomials {
                    public:
                        using field_type = FieldType;
                        using field_value_type = typename field_type::value_type;
                        using constraint_system_type = modified_sqap_constraint_system<field_type>;
                        using constraint_type = typename constraint_system_type::constraint_type;
                        using domain_type = math::basic_radix2_domain<field_type>;

                        /**
                         * Smallest supported power of two >= max(2, logical_rows), without allocation.
                         * Reject zero rows, unsupported roots and unrepresentable storage sizes.
                         */
                        static std::size_t get_domain_size(std::size_t logical_rows) {
                            if (logical_rows == 0) {
                                throw std::invalid_argument("modified_sqap: expected at least one logical row");
                            }

                            const std::size_t log_m = std::max<std::size_t>(1, std::bit_width(logical_rows - 1));
                            if (log_m >= std::numeric_limits<std::size_t>::digits) {
                                throw std::invalid_argument("modified_sqap: domain size overflow");
                            }
                            // Existing unity-root and radix-two FFT helpers use shifts of 1u.
                            if (log_m > algebra::fields::arithmetic_params<field_type>::s ||
                                log_m >= std::numeric_limits<unsigned int>::digits) {
                                throw std::invalid_argument("modified_sqap: unsupported radix-two domain size");
                            }

                            const std::size_t m = std::size_t(1) << log_m;
                            // Z needs m + 1 coefficients; padding needs m constraint rows.
                            if (m >= std::vector<field_value_type>().max_size() ||
                                m > std::vector<constraint_type>().max_size()) {
                                throw std::invalid_argument("modified_sqap: domain storage size overflow");
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
                                throw std::invalid_argument("modified_sqap: invalid constraint system");
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
                    };

                }    // namespace reductions
            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SQAP_TO_POLYNOMIALS_HPP
