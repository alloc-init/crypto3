#pragma once

#include <array>

#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/types.hpp>

namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    using namespace nil::crypto3::algebra::fields::detail::fp12_fast;
    const size_t base_value_limb_count = 4u;
    const size_t storage_limb_count = 8u;
    using limb_array = std::array<limb, storage_limb_count>;
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops