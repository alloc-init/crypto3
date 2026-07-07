#pragma once

#include <array>

#include <nil/crypto3/algebra/fields/detail/element/fp12/fp12_limb_types.hpp>

namespace nil::crypto3::algebra::fields::detail::bls12_fp12_limb_ops {
    const size_t base_value_limb_count = 6u;
    const size_t storage_limb_count = 12u;
    using limb_array = std::array<fp12::limb, storage_limb_count>;
}    // namespace nil::crypto3::algebra::fields::detail::bls12_fp12_limb_ops