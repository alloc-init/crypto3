#pragma once

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>

#define bn254_fp12_multiply_partial_x86(X_INDEX, Y_INDEX) \
    "movq " #X_INDEX "(%%rbx), %%rdx\n"                    \
    "mulx " #Y_INDEX "(%%rcx), %%r11, %%r12\n"          \
    "adcx %%r11, %%r8\n"                              \
    "adcx %%r12, %%r9\n"                             \
    "adcx %[zero], %%r10\n"

#define bn254_fp12_multiply_emit_x86(INDEX) \
    "movq %%r8, " #INDEX "(%%rax)\n" \
    "movq %%r9, %%r8\n"               \
    "movq %%r10, %%r9\n"               \
    "xor %%r10, %%r10\n"                

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {
                namespace detail {
                    namespace alt_bn128_fp12_limb_ops {
                        inline void multiply_4x4_x86_bmi2_adx(limb_array &result, const limb_array &x,
                                                              const limb_array &y) {
                            limb zero = 0;

                            asm volatile(
                                ".byte 0x0f, 0x1f, 0x40, 0x00\n\t"  // 4-byte NOP marker

                                // initialize accumulators
                                "xor %%r8, %%r8\n"
                                "xor %%r9, %%r9\n"
                                "xor %%r10, %%r10\n"

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

                                "movq %%r8, 56(%%rax)\n"
                                "movq %%r9, 64(%%rax)\n"

                                ".byte 0x0f, 0x1f, 0x40, 0x00\n\t"  // 4-byte NOP marker
                                :
                                :   "a"(result.data()), "b"(x.data()), "c"(y.data()), [zero] "r"(zero)
                                :   "r8", // acc0
                                    "r9", // acc1
                                    "r10", // acc2
                                    "r11", // low
                                    "r12", // high
                                    "rdx", // used by mulx 
                                    "cc", "memory");
                        }

                        inline void multiply_5x5_x86_bmi2_adx(limb_array &result, const limb_array &x,
                                                              const limb_array &y) {
                            limb zero = 0;

                            asm volatile(
                                // initialize accumulators
                                "xor %%r8, %%r8\n"
                                "xor %%r9, %%r9\n"
                                "xor %%r10, %%r10\n"

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

                                 // round 9, x4*y4
                                bn254_fp12_multiply_partial_x86(32, 32)
                                bn254_fp12_multiply_emit_x86(64)

                                :
                                :   "a"(result.data()), "b"(x.data()), "c"(y.data()), [zero] "r"(zero)
                                :   "r8", // acc0
                                    "r9", // acc1
                                    "r10", // acc2
                                    "r11", // low
                                    "r12", // high
                                    "rdx", // used by mulx 
                                    "cc", "memory");
                        }
                    }    // namespace alt_bn128_fp12_limb_ops
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil
