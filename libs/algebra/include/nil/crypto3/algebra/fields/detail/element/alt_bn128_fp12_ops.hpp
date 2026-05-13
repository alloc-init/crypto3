//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_ALGEBRA_FIELDS_DETAIL_ELEMENT_ALT_BN128_FP12_OPS_HPP
#define CRYPTO3_ALGEBRA_FIELDS_DETAIL_ELEMENT_ALT_BN128_FP12_OPS_HPP

#include <cstddef>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/fp2.hpp>
#include <nil/crypto3/algebra/fields/fp6_3over2.hpp>

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {
                namespace detail {

                    template<typename BaseField>
                    struct alt_bn128_fp12_ops {
                        constexpr static bool is_supported = false;
                    };

                    template<>
                    struct alt_bn128_fp12_ops<fields::alt_bn128_base_field<254>> {
                        constexpr static bool is_supported = true;

                        using base_field_type = fields::alt_bn128_base_field<254>;
                        using fp_type = typename base_field_type::value_type;

                        using fp2_field_type = fields::fp2<base_field_type>;
                        using fp2_type = typename fp2_field_type::value_type;

                        using fp6_field_type = fields::fp6_3over2<base_field_type>;
                        using fp6_type = typename fp6_field_type::value_type;

                        constexpr static std::size_t fp12_arity = 12;

                        static fp2_type non_residue() {
                            return fp2_type(fp_type(9u), fp_type(1u));
                        }

                        static fp2_type mul_by_xi(fp2_type const& x) {
                            return fp2_type(x.data[0].doubled().doubled().doubled() +
                                                x.data[0] - x.data[1],
                                            x.data[1].doubled().doubled().doubled() +
                                                x.data[1] + x.data[0]);
                        }

                        static fp6_type mul_by_v(fp6_type const& x) {
                            return fp6_type(mul_by_xi(x.data[2]), x.data[0], x.data[1]);
                        }

                        static fp6_type mul_v_add(fp6_type const& x, fp6_type const& y) {
                            return fp6_type(mul_by_xi(x.data[2]) + y.data[0],
                                            x.data[0] + y.data[1],
                                            x.data[1] + y.data[2]);
                        }

                        template<typename Fp12Value>
                        static Fp12Value mul_reduced_mcl_shape(Fp12Value const& x,
                                                               Fp12Value const& y) {
                            fp6_type const& a = x.data[0];
                            fp6_type const& b = x.data[1];
                            fp6_type const& c = y.data[0];
                            fp6_type const& d = y.data[1];

                            fp6_type const ac = a * c;
                            fp6_type const bd = b * d;

                            return Fp12Value(mul_v_add(bd, ac),
                                             (a + b) * (c + d) - ac - bd);
                        }
                    };

                }    // namespace detail
            }        // namespace fields
        }            // namespace algebra
    }                // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ALGEBRA_FIELDS_DETAIL_ELEMENT_ALT_BN128_FP12_OPS_HPP
