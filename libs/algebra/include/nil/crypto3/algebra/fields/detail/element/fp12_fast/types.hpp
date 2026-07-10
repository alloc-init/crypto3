#pragma once

#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <concepts>
#include <type_traits>

#include <nil/crypto3/algebra/fields/field.hpp>

namespace nil::crypto3::algebra::fields::detail::fp12_fast {
    using limb = std::uint64_t;
#if defined(__SIZEOF_INT128__)
    using wide_limb = unsigned __int128;
#else
#error "bls12 fp12 limb ops require unsigned __int128 support"
#endif
    const size_t limb_bits = sizeof(limb) * CHAR_BIT;

    template<typename Params>
    concept Fp12FastParams =
        requires(std::array<typename Params::limb_array, 2> &dst, const std::array<typename Params::limb_array, 2> &src,
                 const std::array<typename Params::limb_array, 2> &addend) {
            typename Params::base_field_type;     // Fp
            typename Params::non_residue_type;    // Fp2
            typename Params::underlying_type;     // Fp6
            typename Params::limb_array;

            typename std::integral_constant<std::size_t, Params::base_value_limb_count>;
            typename std::integral_constant<std::size_t, Params::storage_limb_count>;

            requires Params::storage_limb_count == Params::base_value_limb_count * 2;

            requires std::derived_from<typename Params::base_field_type, field_base>;
            requires std::same_as<typename Params::limb_array, std::array<limb, Params::storage_limb_count>>;

            { Params::fp2_mul_xi_add(dst, src, addend) } -> std::same_as<void>;
        };
}    // namespace nil::crypto3::algebra::fields::detail::fp12_fast
