#pragma once

#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {
                namespace detail {
                    namespace alt_bn128_fp12_limb_ops {

                        using limb = std::uint64_t;
#if defined(__SIZEOF_INT128__)
                        using wide_limb = unsigned __int128;
#else
#error "alt_bn128 fp12 limb ops require unsigned __int128 support"
#endif

                        const size_t limb_bits = sizeof(limb) * CHAR_BIT;
                        const size_t base_value_limb_count = 4u;
                        const size_t storage_limb_count = 8u;

                        using limb_array = std::array<limb, storage_limb_count>;
                    }    // namespace alt_bn128_fp12_limb_ops
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil