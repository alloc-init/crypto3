#pragma once

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>
#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/limb_ops.hpp>

namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    template<class Field>
    inline void fp2_mul_xi_add(limb_array *dst, const limb_array *src, const limb_array *addend) {
        limb_array buf[2];
        mul_8_limbs_by_9<Field>(buf[0].data(), src[0].data());
        subtract_limbs_mod<Field, 8>(buf[0].data(), buf[0].data(), src[1].data());
        mul_8_limbs_by_9<Field>(buf[1].data(), src[1].data());
        add_limbs_mod<Field, 8>(buf[1].data(), buf[1].data(), src[0].data());
        add_limbs_mod<Field, 8>(buf[0].data(), buf[0].data(), addend[0].data());
        add_limbs_mod<Field, 8>(buf[1].data(), buf[1].data(), addend[1].data());
        dst[0] = buf[0];
        dst[1] = buf[1];
    }
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops
