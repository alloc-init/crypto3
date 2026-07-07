#pragma once

#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>

namespace nil::crypto3::algebra::fields::detail::fp12_fast {
    using limb = std::uint64_t;
#if defined(__SIZEOF_INT128__)
    using wide_limb = unsigned __int128;
#else
#error "bls12 fp12 limb ops require unsigned __int128 support"
#endif
    const size_t limb_bits = sizeof(limb) * CHAR_BIT;
}
