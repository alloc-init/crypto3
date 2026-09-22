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

#ifndef CRYPTO3_ZK_MODIFIED_SQAP_VERIFICATION_KEY_HPP
#define CRYPTO3_ZK_MODIFIED_SQAP_VERIFICATION_KEY_HPP

#include <array>
#include <cstddef>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {

                /**
                 * Verification parameters and circuit binding for the modified SQAP SNARK.
                 * Default construction does not produce a valid verification key.
                 */
                template<typename CurveType>
                struct modified_sqap_verification_key {
                    using curve_type = CurveType;
                    using base_field_type = typename curve_type::base_field_type;
                    using base_value_type = typename base_field_type::value_type;
                    using g2_type = typename curve_type::template g2_type<>;
                    using g2_value_type = typename g2_type::value_type;
                    using gt_type = typename curve_type::gt_type;
                    using gt_value_type = typename gt_type::value_type;

                    g2_value_type g2_one = g2_value_type::one();
                    g2_value_type tau_g2 = g2_value_type::zero();
                    g2_value_type gamma_inverse_g2 = g2_value_type::zero();

                    gt_value_type alpha_z_vanishing_gt = gt_value_type::one();
                    // [alpha_U]_T in A, C, H, Z order. GT's identity is field one.
                    std::array<gt_value_type, 4> alpha_gt = {gt_value_type::one(), gt_value_type::one(),
                                                             gt_value_type::one(), gt_value_type::one()};

                    std::size_t num_variables = 0;
                    std::size_t domain_size = 0;
                    base_value_type circuit_digest = base_value_type::zero();

                    bool operator==(const modified_sqap_verification_key &other) const = default;
                };

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_MODIFIED_SQAP_VERIFICATION_KEY_HPP
