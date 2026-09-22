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

#ifndef CRYPTO3_ZK_MODIFIED_SQAP_PROVING_KEY_HPP
#define CRYPTO3_ZK_MODIFIED_SQAP_PROVING_KEY_HPP

#include <vector>

#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/modified_sqap.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/verification_key.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {

                /**
                 * Owned proving parameters for the modified SQAP SNARK.
                 * Default construction does not produce a valid proving key.
                 */
                template<typename CurveType, typename ConstraintSystem =
                                                 modified_sqap_constraint_system<typename CurveType::scalar_field_type>>
                struct modified_sqap_proving_key {
                    using curve_type = CurveType;
                    using g1_type = typename curve_type::template g1_type<>;
                    using g1_value_type = typename g1_type::value_type;
                    using constraint_system_type = ConstraintSystem;
                    using verification_key_type = modified_sqap_verification_key<curve_type>;

                    // n = verification_key.num_variables; m = verification_key.domain_size.
                    std::vector<g1_value_type> W;          // n mixed witness bases.
                    std::vector<g1_value_type> H_query;    // m - 1 quotient commitment bases.
                    std::vector<g1_value_type> T_A;        // m - 1 opening bases.
                    std::vector<g1_value_type> T_C;        // m - 1 opening bases.
                    std::vector<g1_value_type> T_H;        // m - 2 opening bases; empty when m = 2.
                    std::vector<g1_value_type> T_Z;        // m opening bases.

                    // Canonical logical rows, including the binding row, before padding.
                    constraint_system_type constraint_system;
                    verification_key_type verification_key;

                    bool operator==(const modified_sqap_proving_key &other) const = default;
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SQAP_PROVING_KEY_HPP
