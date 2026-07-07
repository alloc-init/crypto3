#pragma once

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>
#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/limb_ops.hpp>

#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_ops_x86.hpp>
#endif

namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {

    template<class Field>
    inline void mul_8_limbs_by_9(limb_array &dst, const limb_array &src) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        mul_8_limbs_by_9_x86<Field>(dst, src);
#else
        limb_array cpy = src;
        add_limbs_mod<Field, 8>(dst.data(), src.data(), src.data());    // 2x
        add_limbs_mod<Field, 8>(dst.data(), dst.data(), dst.data());    // 4x
        add_limbs_mod<Field, 8>(dst.data(), dst.data(), dst.data());    // 8x
        add_limbs_mod<Field, 8>(dst.data(), dst.data(), cpy.data());    // 9x
#endif
    }

    // Out-of-place dst = src * xi + addend. Callers must use one of the
    // modify_* variants when dst aliases an input.
    template<class Field>
    inline void fp2_mul_xi_add(limb_array *dst, const limb_array *src, const limb_array *addend) {
        mul_8_limbs_by_9<Field>(dst[0], src[0]);
        subtract_limbs_mod<Field, 8>(dst[0].data(), dst[0].data(), src[1].data());
        mul_8_limbs_by_9<Field>(dst[1], src[1]);
        add_limbs_mod<Field, 8>(dst[1].data(), dst[1].data(), src[0].data());
        add_limbs_mod<Field, 8>(dst[0].data(), dst[0].data(), addend[0].data());
        add_limbs_mod<Field, 8>(dst[1].data(), dst[1].data(), addend[1].data());
    }

    template<class Field>
    inline void fp2_mul_xi_add_modify_src(limb_array *src, const limb_array *addend) {
        limb_array src0 = src[0];
        mul_8_limbs_by_9<Field>(src[0], src0);
        subtract_limbs_mod<Field, 8>(src[0].data(), src[0].data(), src[1].data());
        mul_8_limbs_by_9<Field>(src[1], src[1]);
        add_limbs_mod<Field, 8>(src[1].data(), src[1].data(), src0.data());
        add_limbs_mod<Field, 8>(src[0].data(), src[0].data(), addend[0].data());
        add_limbs_mod<Field, 8>(src[1].data(), src[1].data(), addend[1].data());
    }

    template<class Field>
    inline void fp2_mul_xi_add_modify_addend(limb_array *addend, const limb_array *src) {
        limb_array tmp;
        mul_8_limbs_by_9<Field>(tmp, src[0]);
        subtract_limbs_mod<Field, 8>(tmp.data(), tmp.data(), src[1].data());
        add_limbs_mod<Field, 8>(addend[0].data(), addend[0].data(), tmp.data());
        mul_8_limbs_by_9<Field>(tmp, src[1]);
        add_limbs_mod<Field, 8>(tmp.data(), tmp.data(), src[0].data());
        add_limbs_mod<Field, 8>(addend[1].data(), addend[1].data(), tmp.data());
    }

}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops
