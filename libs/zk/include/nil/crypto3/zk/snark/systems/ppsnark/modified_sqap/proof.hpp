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

#ifndef CRYPTO3_ZK_MODIFIED_SQAP_PROOF_HPP
#define CRYPTO3_ZK_MODIFIED_SQAP_PROOF_HPP

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {

                /**
                 * A modified SQAP proof. The public input is supplied separately.
                 */
                template<typename CurveType>
                struct modified_sqap_proof {
                    using curve_type = CurveType;
                    using scalar_field_type = typename curve_type::scalar_field_type;
                    using scalar_value_type = typename scalar_field_type::value_type;
                    using g1_type = typename curve_type::template g1_type<>;
                    using g1_value_type = typename g1_type::value_type;

                    g1_value_type P = g1_value_type::zero();
                    g1_value_type Q = g1_value_type::zero();
                    scalar_value_type v_A = scalar_value_type::zero();
                    scalar_value_type v_C = scalar_value_type::zero();
                    scalar_value_type v_H = scalar_value_type::zero();
                    scalar_value_type v_Z = scalar_value_type::zero();

                    bool operator==(const modified_sqap_proof &other) const = default;
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SQAP_PROOF_HPP
