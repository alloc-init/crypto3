#pragma once

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>

#define bn254_fp12_multiply_partial_x86(X_INDEX, Y_INDEX)   \
    "movq " #X_INDEX "(%[x]), %%rax\n"                      \
    "mulq " #Y_INDEX "(%[y])\n"                             \
    "add %%rax, %[acc0]\n"                                  \
    "adc %%rdx, %[acc1]\n"                                  \
    "adc $0, %[acc2]\n"

#define bn254_fp12_multiply_emit_x86(INDEX) \
    "movq %[acc0], " #INDEX "(%[result])\n" \
    "movq %[acc1], %[acc0]\n"               \
    "movq %[acc2], %[acc1]\n"               \
    "xor %[acc2], %[acc2]\n"

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {
                namespace detail {
                    namespace alt_bn128_fp12_limb_ops {
                        inline void multiply_4x4_x86(limb_array &result, const limb_array &x, const limb_array &y) {
                            limb acc0 = 0;
                            limb acc1 = 0;
                            limb acc2 = 0;

                            asm volatile(
                                // round 1, x0*y0
                                bn254_fp12_multiply_partial_x86(0, 0)
                                bn254_fp12_multiply_emit_x86(0)

                                // round 2, x0 * y1, x1 * y0
                                bn254_fp12_multiply_partial_x86(0, 8)
                                bn254_fp12_multiply_partial_x86(8, 0)
                                bn254_fp12_multiply_emit_x86(8)

                                // round 3, x0*y2, x1*y1, x2*y0
                                bn254_fp12_multiply_partial_x86(0,16)
                                bn254_fp12_multiply_partial_x86(8,8)
                                bn254_fp12_multiply_partial_x86(16,0)
                                bn254_fp12_multiply_emit_x86(16)

                                 // round 4, x0*y3, x1*y2, x2*y1, x3*y0
                                bn254_fp12_multiply_partial_x86(0, 24)
                                bn254_fp12_multiply_partial_x86(8, 16)
                                bn254_fp12_multiply_partial_x86(16, 8)
                                bn254_fp12_multiply_partial_x86(24, 0)
                                bn254_fp12_multiply_emit_x86(24)

                                 // round 5, x1*y3, x2*y2, x3*y1
                                bn254_fp12_multiply_partial_x86(8, 24)
                                bn254_fp12_multiply_partial_x86(16, 16)
                                bn254_fp12_multiply_partial_x86(24, 8)
                                bn254_fp12_multiply_emit_x86(32)

                                 // round 6, x2*y3, x3*y2
                                bn254_fp12_multiply_partial_x86(16, 24)
                                bn254_fp12_multiply_partial_x86(24, 16)
                                bn254_fp12_multiply_emit_x86(40)

                                 // round 7, x3*y3
                                bn254_fp12_multiply_partial_x86(24, 24)
                                bn254_fp12_multiply_emit_x86(48)
                                "movq %[acc0], 56(%[result])\n"
                                "movq %[acc1], 64(%[result])\n"


                                :   [acc0]"+r"(acc0),
                                    [acc1]"+r"(acc1),
                                    [acc2]"+r"(acc2)
                                :   [result]"r"(result.data()),
                                    [x]"r"(x.data()),
                                    [y]"r"(y.data())
                                :   "rax", "rdx", "cc", "memory"
                            );
                        }

                        inline void multiply_5x5_x86(limb_array &result, const limb_array &x, const limb_array &y) {
                            limb acc0 = 0;
                            limb acc1 = 0;
                            limb acc2 = 0;

                            asm volatile(
                                // round 1, x0*y0
                                bn254_fp12_multiply_partial_x86(0, 0)
                                bn254_fp12_multiply_emit_x86(0)

                                // round 2, x0 * y1, x1 * y0
                                bn254_fp12_multiply_partial_x86(0, 8)
                                bn254_fp12_multiply_partial_x86(8, 0)
                                bn254_fp12_multiply_emit_x86(8)

                                // round 3, x0*y2, x1*y1, x2*y0
                                bn254_fp12_multiply_partial_x86(0, 16)
                                bn254_fp12_multiply_partial_x86(8, 8)
                                bn254_fp12_multiply_partial_x86(16, 0)
                                bn254_fp12_multiply_emit_x86(16)

                                 // round 4, x0*y3, x1*y2, x2*y1, x3*y0
                                bn254_fp12_multiply_partial_x86(0, 24)
                                bn254_fp12_multiply_partial_x86(8, 16)
                                bn254_fp12_multiply_partial_x86(16, 8)
                                bn254_fp12_multiply_partial_x86(24, 0)
                                bn254_fp12_multiply_emit_x86(24)

                                 // round 5, x0*y4, x1*y3, x2*y2, x3*y1, x4*y0
                                bn254_fp12_multiply_partial_x86(0, 32)
                                bn254_fp12_multiply_partial_x86(8, 24)
                                bn254_fp12_multiply_partial_x86(16, 16)
                                bn254_fp12_multiply_partial_x86(24, 8)
                                bn254_fp12_multiply_partial_x86(32, 0)
                                bn254_fp12_multiply_emit_x86(32)

                                 // round 6, x1*y4, x2*y3, x3*y2, x4*y1
                                bn254_fp12_multiply_partial_x86(8, 32)
                                bn254_fp12_multiply_partial_x86(16, 24)
                                bn254_fp12_multiply_partial_x86(24, 16)
                                bn254_fp12_multiply_partial_x86(32, 8)
                                bn254_fp12_multiply_emit_x86(40)

                                 // round 7, x2*y4, x3*y3, x4*y2
                                bn254_fp12_multiply_partial_x86(16, 32)
                                bn254_fp12_multiply_partial_x86(24, 24)
                                bn254_fp12_multiply_partial_x86(32, 16)
                                bn254_fp12_multiply_emit_x86(48)

                                 // round 8, x3*y4, x4*y3
                                bn254_fp12_multiply_partial_x86(24, 32)
                                bn254_fp12_multiply_partial_x86(32, 24)
                                bn254_fp12_multiply_emit_x86(56)

                                 // round 9, x4*y4, drop high bits (see note in portable version)
                                "movq 32(%[x]), %%rax\n"
                                "mulq 32(%[y])\n"
                                "add %%rax, %[acc0]\n"
                                "movq %[acc0], 64(%[result])\n"

                                :   [acc0]"+r"(acc0),
                                    [acc1]"+r"(acc1),
                                    [acc2]"+r"(acc2)
                                :   [result]"r"(result.data()),
                                    [x]"r"(x.data()),
                                    [y]"r"(y.data())
                                :   "rax", "rdx", "cc", "memory"
                            );
                        }

                        void montgomery_reduce_x86(limb_array &data, limb_array &p, limb p_dash) {
                            return;
                        }

                    }    // namespace alt_bn128_fp12_limb_ops
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil