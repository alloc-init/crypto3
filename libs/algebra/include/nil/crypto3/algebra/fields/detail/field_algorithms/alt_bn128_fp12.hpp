//---------------------------------------------------------------------------//
// Copyright (c) 2026 Alloc Init Labs Inc.
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

#ifndef CRYPTO3_ALGEBRA_FIELDS_DETAIL_ALT_BN128_FP12_ALGORITHMS_HPP
#define CRYPTO3_ALGEBRA_FIELDS_DETAIL_ALT_BN128_FP12_ALGORITHMS_HPP

#include <cassert>
#include <cstddef>

#include <boost/multiprecision/cpp_int.hpp>

#include <nil/crypto3/algebra/fields/detail/field_algorithms.hpp>
#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/field_order.hpp>
#include <nil/crypto3/algebra/fields/fp2.hpp>
#include <nil/crypto3/algebra/fields/fp6_3over2.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

namespace nil::crypto3::algebra::fields::detail {

    using alt_bn128_254_base_field = ::nil::crypto3::algebra::fields::alt_bn128_base_field<254>;
    using alt_bn128_254_fp2_field = ::nil::crypto3::algebra::fields::fp2<alt_bn128_254_base_field>;
    using alt_bn128_254_fp6_field = ::nil::crypto3::algebra::fields::fp6_3over2<alt_bn128_254_base_field>;
    using alt_bn128_254_fp12_field = ::nil::crypto3::algebra::fields::fp12_2over3over2<alt_bn128_254_base_field>;
    using alt_bn128_254_fp_value = typename alt_bn128_254_base_field::value_type;
    using alt_bn128_254_fp2_value = typename alt_bn128_254_fp2_field::value_type;
    using alt_bn128_254_fp6_value = typename alt_bn128_254_fp6_field::value_type;
    using alt_bn128_254_fp12_value = typename alt_bn128_254_fp12_field::value_type;

    struct alt_bn128_fp4_subfield_value {
        alt_bn128_254_fp2_value c0;
        alt_bn128_254_fp2_value c1;

        static const alt_bn128_254_fp2_value &non_residue() {
            static const alt_bn128_254_fp2_value value(alt_bn128_254_fp_value(9), alt_bn128_254_fp_value(1));
            return value;
        }

        alt_bn128_fp4_subfield_value squared() const {
            const alt_bn128_254_fp2_value aa = c0.squared();
            const alt_bn128_254_fp2_value bb = c1.squared();
            const alt_bn128_254_fp2_value ab = c0 * c1;
            return {aa + non_residue() * bb, ab + ab};
        }
    };

    inline alt_bn128_fp4_subfield_value operator*(const alt_bn128_fp4_subfield_value &a,
                                                  const alt_bn128_fp4_subfield_value &b) {
        const alt_bn128_254_fp2_value a0b0 = a.c0 * b.c0;
        const alt_bn128_254_fp2_value a1b1 = a.c1 * b.c1;
        return {a0b0 + alt_bn128_fp4_subfield_value::non_residue() * a1b1, (a.c0 + a.c1) * (b.c0 + b.c1) - a0b0 - a1b1};
    }

    inline alt_bn128_fp4_subfield_value pow(alt_bn128_fp4_subfield_value x, boost::multiprecision::cpp_int exponent) {
        alt_bn128_fp4_subfield_value result {alt_bn128_254_fp2_value::one(), alt_bn128_254_fp2_value::zero()};
        while (exponent != 0) {
            if (boost::multiprecision::bit_test(exponent, 0)) {
                result = result * x;
            }
            exponent >>= 1;
            if (exponent != 0) {
                x = x.squared();
            }
        }
        return result;
    }

    inline alt_bn128_254_fp12_value fp12_from_fp4_subfield(const alt_bn128_fp4_subfield_value &x) {
        return alt_bn128_254_fp12_value(
            alt_bn128_254_fp6_value(x.c0, alt_bn128_254_fp2_value::zero(), alt_bn128_254_fp2_value::zero()),
            alt_bn128_254_fp6_value(alt_bn128_254_fp2_value::zero(), x.c1, alt_bn128_254_fp2_value::zero()));
    }

    inline alt_bn128_fp4_subfield_value fp12_to_fp4_subfield(const alt_bn128_254_fp12_value &x) {
        assert(x.data[0].data[1].is_zero());
        assert(x.data[0].data[2].is_zero());
        assert(x.data[1].data[0].is_zero());
        assert(x.data[1].data[2].is_zero());
        return {x.data[0].data[0], x.data[1].data[1]};
    }

    inline alt_bn128_254_fp6_value fp12_to_fp6_subfield(const alt_bn128_254_fp12_value &x) {
        assert(x.data[1].is_zero());
        return x.data[0];
    }

    inline alt_bn128_254_fp12_value fp12_norm_one_sqrt_frobenius(const alt_bn128_254_fp12_value &x) {
        static const boost::multiprecision::cpp_int exponent = (field_order<alt_bn128_254_fp2_field>() - 1) >> 2;
        const alt_bn128_254_fp12_value y = x.pow(exponent);
        return x * y * y.Frobenius_map(2) * y.Frobenius_map(4);
    }

    inline alt_bn128_254_fp6_value fp6_fourth_root_odd_frobenius(const alt_bn128_254_fp6_value &x) {
        static const boost::multiprecision::cpp_int q = field_order<alt_bn128_254_fp2_field>();
        static const boost::multiprecision::cpp_int a = (q - 1) >> 4;
        static const boost::multiprecision::cpp_int b = (boost::multiprecision::cpp_int(3) * a - 1) >> 2;

        assert((q - 1) % 16 == 0);
        assert((boost::multiprecision::cpp_int(3) * a - 1) % 4 == 0);

        const alt_bn128_254_fp6_value y = x.pow(b);
        const alt_bn128_254_fp6_value z = x.pow(a);
        const alt_bn128_254_fp6_value z4 = z.squared().squared();
        const alt_bn128_254_fp6_value z8 = z4.squared();
        return y * y.Frobenius_map(2) * y.Frobenius_map(4) * z4.Frobenius_map(2) * x * z8;
    }

    inline alt_bn128_254_fp12_value
        alt_bn128_fp12_two_primary_component(const alt_bn128_254_fp12_value &x,
                                             const boost::multiprecision::cpp_int &odd_order,
                                             std::size_t two_adicity) {
        static const boost::multiprecision::cpp_int expected_odd_order =
            (field_order<alt_bn128_254_fp12_field>() - 1) >> 5;
        if (two_adicity != 5 || odd_order != expected_odd_order) {
            return x.pow(odd_order);
        }

        static const boost::multiprecision::cpp_int fp4_odd_order = []() {
            const boost::multiprecision::cpp_int p = field_characteristic<alt_bn128_254_base_field>();
            boost::multiprecision::cpp_int order = p * p * p * p;
            order -= 1;
            order >>= 5;
            return order;
        }();
        const alt_bn128_254_fp12_value norm = x * x.Frobenius_map(4) * x.Frobenius_map(8);
        return fp12_from_fp4_subfield(pow(fp12_to_fp4_subfield(norm), fp4_odd_order));
    }

    inline alt_bn128_254_fp12_value alt_bn128_fp12_odd_subgroup_sqrt(const alt_bn128_254_fp12_value &x,
                                                                     const boost::multiprecision::cpp_int &odd_order) {
        if (x.is_zero() || x.is_one()) {
            return x;
        }

        static const boost::multiprecision::cpp_int expected_odd_order =
            (field_order<alt_bn128_254_fp12_field>() - 1) >> 5;
        if (odd_order != expected_odd_order) {
            assert((odd_order & 1) != 0);
            return x.pow((odd_order + 1) >> 1);
        }

        const alt_bn128_254_fp6_value norm = fp12_to_fp6_subfield(x * x.Frobenius_map(6));
        const alt_bn128_254_fp6_value fourth_root = fp6_fourth_root_odd_frobenius(norm);
        const alt_bn128_254_fp6_value fp6_component = fourth_root.squared();
        const alt_bn128_254_fp12_value norm_one = x * alt_bn128_254_fp12_value(fp6_component.inversed());
        const alt_bn128_254_fp12_value result =
            alt_bn128_254_fp12_value(fourth_root) * fp12_norm_one_sqrt_frobenius(norm_one);
        assert(result.squared() == x);
        return result;
    }

    inline alt_bn128_254_fp12_value alt_bn128_fp12_primitive_root(std::size_t two_adicity) {
        if (two_adicity != 5) {
            throw std::invalid_argument("BN254 Fp12 primitive root discovery supports two-adicity 5");
        }

        static const alt_bn128_254_fp12_value root = []() {
            const boost::multiprecision::cpp_int exponent = (field_order<alt_bn128_254_fp12_field>() - 1) >> 5;
            for (unsigned candidate = 2; candidate < 64; ++candidate) {
                const alt_bn128_254_fp12_value value(alt_bn128_254_fp6_value::zero(),
                                                     alt_bn128_254_fp6_value(alt_bn128_254_fp_value(candidate)));
                const alt_bn128_254_fp12_value possible_root = value.pow(exponent);
                if (!possible_root.is_one() && possible_root.pow(16) != alt_bn128_254_fp12_value::one() &&
                    possible_root.pow(32) == alt_bn128_254_fp12_value::one()) {
                    return possible_root;
                }
            }
            assert(false && "failed to find a primitive 32nd root of unity");
            return alt_bn128_254_fp12_value::one();
        }();
        return root;
    }

    template<>
    struct field_algorithms<alt_bn128_254_fp12_value> {
        static alt_bn128_254_fp12_value two_primary_component(const alt_bn128_254_fp12_value &x,
                                                              const boost::multiprecision::cpp_int &odd_order,
                                                              std::size_t two_adicity) {
            return alt_bn128_fp12_two_primary_component(x, odd_order, two_adicity);
        }

        static alt_bn128_254_fp12_value odd_subgroup_sqrt(const alt_bn128_254_fp12_value &x,
                                                          const boost::multiprecision::cpp_int &odd_order) {
            return alt_bn128_fp12_odd_subgroup_sqrt(x, odd_order);
        }

        static alt_bn128_254_fp12_value sqrt_known_square(const alt_bn128_254_fp12_value &x) {
            return x.sqrt_known_square();
        }

        static alt_bn128_254_fp12_value primitive_two_power_root_of_unity(std::size_t two_adicity) {
            return alt_bn128_fp12_primitive_root(two_adicity);
        }
    };

}    // namespace nil::crypto3::algebra::fields::detail

#endif    // CRYPTO3_ALGEBRA_FIELDS_DETAIL_ALT_BN128_FP12_ALGORITHMS_HPP
