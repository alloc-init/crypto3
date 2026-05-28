#pragma once

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {
                namespace detail {
                    namespace alt_bn128_fp12_limb_ops {
                        inline void multiply_4x4_x86_bmi2_adx(limb_array &result, const limb_array &x,
                                                              const limb_array &y) {
                            limb *r = result.data();
                            const limb *xp = x.data();
                            const limb *yp = y.data();

                            // asm volatile(
                            //     // full mulx/adcx/adox Comba kernel here
                            //     :
                            //     : [r] "r"(r), [x] "r"(xp), [y] "r"(yp)
                            //     : "rax", "rdx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "cc",
                            //     "memory");
                        }
                    }    // namespace alt_bn128_fp12_limb_ops
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil
