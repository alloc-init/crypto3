#pragma once

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>
#include <nil/crypto3/algebra/fields/detail/element/fp12_fast_utils/fp12_limb_ops.hpp>

#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_ops_x86.hpp>
#endif

namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    template<class Field>
    inline void add_8_limbs_mod(limb_array &z, const limb_array &x, const limb_array &y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        add_8_limbs_mod_x86<Field>(z, x, y);
#else
        add_limbs_portable<8>(z.data(), x.data(), y.data());
        static const auto p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
        if (ge_modulus<4>(z.data() + 4, p.data())) {
            subtract_limbs_portable<4>(z.data() + 4, z.data() + 4, p.data());
        }
#endif
    }

    template<class Field>
    inline void subtract_8_limbs_mod(limb_array &z, const limb_array &x, const limb_array &y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        subtract_8_limbs_mod_x86<Field>(z, x, y);
#else
        bool borrow = subtract_limbs_portable<8>(z.data(), x.data(), y.data());
        if (borrow) {
            // If the full 8-limb subtraction borrowed, add p to the high half
            // to keep the double value in the p * R residue class.
            static const auto p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
            add_limbs_portable<4>(z.data() + 4, z.data() + 4, p.data());
        }
#endif
    }

    template<class Field>
    inline void mul_8_limbs_by_9(limb_array &dst, const limb_array &src) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        mul_8_limbs_by_9_x86<Field>(dst, src);
#else
        limb_array cpy = src;
        add_8_limbs_mod<Field>(dst, src, src);    // 2x
        add_8_limbs_mod<Field>(dst, dst, dst);    // 4x
        add_8_limbs_mod<Field>(dst, dst, dst);    // 8x
        add_8_limbs_mod<Field>(dst, dst, cpy);    // 9x
#endif
    }

    inline void multiply_4x4(limb *z, const limb *x, const limb *y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        multiply_4x4_x86(z, x, y);
#else
        multiply_portable<4>(z, x, y);
#endif
    }

    template<class Field>
    inline void montgomery_reduce(limb *result, const limb *data) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        montgomery_reduce_x86<Field>(result, data);
#else
        montgomery_reduce_portable<Field, 8>(result, data);
#endif
    }

    // fp2_base values are two contiguous 4-limb coefficients; each output
    // coefficient is normalized modulo p.
    template<class Field>
    inline void fp2_base_add_mod(limb *z, const limb *x, const limb *y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        fp2_base_add_mod_x86<Field>(z, x, y);
#else
        add_low_limbs_mod_portable<Field, 4>(z, x, y);
        add_low_limbs_mod_portable<Field, 4>(z + 4, x + 4, y + 4);
#endif
    }

    // Raw fp2_base coefficient-wise add. BN254 inputs are bounded so each
    // coefficient sum fits in four limbs.
    inline void fp2_base_add_pre(limb *z, const limb *x, const limb *y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        fp2_base_add_pre_x86(z, x, y);
#else
        add_limbs_portable<4>(z, x, y);
        add_limbs_portable<4>(z + 4, x + 4, y + 4);
#endif
    }

    template<class Field>
    inline void fp2_sub_pre(limb_array *data, const limb_array *other) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        fp2_sub_pre_x86<Field>(data, other);
#else
        subtract_8_limbs_mod<Field>(data[0], data[0], other[0]);
        subtract_limbs_portable<8>(data[1].data(), data[1].data(), other[1].data());
#endif
    }

    // Out-of-place dst = src * xi + addend. Callers must use one of the
    // modify_* variants when dst aliases an input.
    template<class Field>
    inline void fp2_mul_xi_add(limb_array *dst, const limb_array *src, const limb_array *addend) {
        mul_8_limbs_by_9<Field>(dst[0], src[0]);
        subtract_8_limbs_mod<Field>(dst[0], dst[0], src[1]);
        mul_8_limbs_by_9<Field>(dst[1], src[1]);
        add_8_limbs_mod<Field>(dst[1], dst[1], src[0]);
        add_8_limbs_mod<Field>(dst[0], dst[0], addend[0]);
        add_8_limbs_mod<Field>(dst[1], dst[1], addend[1]);
    }

    template<class Field>
    inline void fp2_mul_xi_add_modify_src(limb_array *src, const limb_array *addend) {
        limb_array src0 = src[0];
        mul_8_limbs_by_9<Field>(src[0], src0);
        subtract_8_limbs_mod<Field>(src[0], src[0], src[1]);
        mul_8_limbs_by_9<Field>(src[1], src[1]);
        add_8_limbs_mod<Field>(src[1], src[1], src0);
        add_8_limbs_mod<Field>(src[0], src[0], addend[0]);
        add_8_limbs_mod<Field>(src[1], src[1], addend[1]);
    }

    template<class Field>
    inline void fp2_mul_xi_add_modify_addend(limb_array *addend, const limb_array *src) {
        limb_array tmp;
        mul_8_limbs_by_9<Field>(tmp, src[0]);
        subtract_8_limbs_mod<Field>(tmp, tmp, src[1]);
        add_8_limbs_mod<Field>(addend[0], addend[0], tmp);
        mul_8_limbs_by_9<Field>(tmp, src[1]);
        add_8_limbs_mod<Field>(tmp, tmp, src[0]);
        add_8_limbs_mod<Field>(addend[1], addend[1], tmp);
    }

    template<class Field>
    inline void fp2_mul_pre_portable(limb_array *z, const limb *x, const limb *y) {
        // For x = a + bu and y = c + du:
        //   xy = (a + bu) * (c + du)
        //      = ac + adu + bcu + bdu^2
        //      = ac + (ad + bc)u - bd      # since u^2 = -1
        //      = (ac - bd) + (ad + bc)u
        // Karatsuba computes the cross term with one product:
        //   ad + bc = (a + b)(c + d) - ac - bd.
        limb_array ac, bd, a_plus_b, c_plus_d;
        multiply_4x4(ac.data(), x, y);
        multiply_4x4(bd.data(), x + 4, y + 4);
        add_limbs_portable<4>(a_plus_b.data(), x, x + 4);
        add_limbs_portable<4>(c_plus_d.data(), y, y + 4);
        subtract_8_limbs_mod<Field>(z[0], ac, bd);
        multiply_4x4(z[1].data(), a_plus_b.data(), c_plus_d.data());
        subtract_limbs_portable<8>(z[1].data(), z[1].data(), ac.data());
        subtract_limbs_portable<8>(z[1].data(), z[1].data(), bd.data());
    }

    template<class Field>
    inline void fp2_mul_pre(limb_array *z, const limb *x, const limb *y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        fp2_mul_pre_x86<Field>(z, x, y);
#else
        fp2_mul_pre_portable<Field>(z, x, y);
#endif
    }

    template<class Field>
    inline void fp2_add_mul_pre_portable(limb_array *z, const limb *a, const limb *b, const limb *c, const limb *d) {
        // Build the raw fp2 sums in the same packed layout expected by fp2_mul_pre.
        limb x[8], y[8];
        add_limbs_portable<4>(x, a, b);
        add_limbs_portable<4>(x + 4, a + 4, b + 4);
        add_limbs_portable<4>(y, c, d);
        add_limbs_portable<4>(y + 4, c + 4, d + 4);
        fp2_mul_pre_portable<Field>(z, x, y);
    }

    template<class Field>
    inline void fp2_add_mul_pre(limb_array *z, const limb *a, const limb *b, const limb *c, const limb *d) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        fp2_add_mul_pre_x86<Field>(z, a, b, c, d);
#else
        fp2_add_mul_pre_portable<Field>(z, a, b, c, d);
#endif
    }

}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops
