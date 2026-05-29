#pragma once

#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__) && (defined(__GNUC__) || defined(__clang__))

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>

#define bn254_fp12_multiply_partial_x86(X_INDEX, Y_INDEX) \
    "movq " #X_INDEX "(%[x]), %%rdx\n"                    \
    "mulx " #Y_INDEX "(%[y]), %[low], %[high]\n"          \
    "adcx %[low], %[acc0]\n"                              \
    "adcx %[high], %[acc1]\n"                             \
    "adcx %%r11, %[acc2]\n"

// two of the multiply_partials with dual carry chains
#define bn254_fp12_multiply_partial_x86_interleaved2(X1_INDEX, Y1_INDEX, X2_INDEX, Y2_INDEX)    \
    "movq " #X1_INDEX "(%[x]), %%rdx\n"                                                         \
    "mulx " #Y1_INDEX "(%[y]), %[low], %[high]\n"                                               \
    "movq " #X2_INDEX "(%[x]), %%rdx\n"                                                         \
    "adcx %[low], %[acc0]\n"                                                                    \
    "mulx " #Y2_INDEX "(%[y]), %[low2], %[high2]\n"                                             \
    "adcx %[high], %[acc1]\n"                                                                   \
    "adox %[low2], %[acc0]\n"                                                                   \
    "adcx %%r11, %[acc2]\n"                                                                     \
    "adox %[high2], %[acc1]\n"                                                                  \
    "adox %%r11, %[acc2]\n"

// three of the multiply_partials with dual carry chains
#define bn254_fp12_multiply_partial_x86_interleaved3(X1_INDEX, Y1_INDEX, X2_INDEX, Y2_INDEX, X3_INDEX, Y3_INDEX)    \
    "movq " #X1_INDEX "(%[x]), %%rdx\n"                                                                             \
    "mulx " #Y1_INDEX "(%[y]), %[low], %[high]\n"                                                                   \
    "movq " #X2_INDEX "(%[x]), %%rdx\n"                                                                             \
    "adcx %[low], %[acc0]\n"                                                                                        \
    "mulx " #Y2_INDEX "(%[y]), %[low2], %[high2]\n"                                                                 \
    "movq " #X3_INDEX "(%[x]), %%rdx\n"                                                                             \
    "adcx %[high], %[acc1]\n"                                                                                       \
    "mulx " #Y3_INDEX "(%[y]), %[low], %[high]\n"                                                                   \
    "adox %[low2], %[acc0]\n"                                                                                       \
    "adcx %%r11, %[acc2]\n"                                                                                         \
    "adcx %[low], %[acc0]\n"                                                                                        \
    "adox %[high2], %[acc1]\n"                                                                                      \
    "adcx %[high], %[acc1]\n"                                                                                       \
    "adox %%r11, %[acc2]\n"                                                                                         \
    "adcx %%r11, %[acc2]\n" 

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
                        inline void multiply_4x4_x86_bmi2_adx(limb_array &result, const limb_array &x,
                                                              const limb_array &y) {
                            limb low, high, low2, high2;
                            limb acc0 = 0;
                            limb acc1 = 0;
                            limb acc2 = 0;

                            asm volatile(
                                "xor %%r11, %%r11\n" // clear zero register

                                // round 1, x0*y0
                                bn254_fp12_multiply_partial_x86(0, 0)
                                bn254_fp12_multiply_emit_x86(0)

                                // round 2, x0 * y1, x1 * y0
                                bn254_fp12_multiply_partial_x86_interleaved2(0, 8, 8, 0)
                                bn254_fp12_multiply_emit_x86(8)

                                // round 3, x0*y2, x1*y1, x2*y0
                                bn254_fp12_multiply_partial_x86_interleaved3(0,16, 8, 8, 16,0)
                                bn254_fp12_multiply_emit_x86(16)

                                 // round 4, x0*y3, x1*y2, x2*y1, x3*y0
                                bn254_fp12_multiply_partial_x86_interleaved2(0, 24, 8, 16)
                                bn254_fp12_multiply_partial_x86_interleaved2(16, 8, 24, 0)
                                bn254_fp12_multiply_emit_x86(24)

                                 // round 5, x1*y3, x2*y2, x3*y1
                                bn254_fp12_multiply_partial_x86_interleaved3(8, 24, 16, 16, 24, 8)
                                bn254_fp12_multiply_emit_x86(32)

                                 // round 6, x2*y3, x3*y2
                                bn254_fp12_multiply_partial_x86_interleaved2(16, 24, 24, 16)
                                bn254_fp12_multiply_emit_x86(40)

                                 // round 7, x3*y3
                                bn254_fp12_multiply_partial_x86(24, 24)
                                bn254_fp12_multiply_emit_x86(48)

                                "movq %[acc0], 56(%[result])\n"
                                "movq %[acc1], 64(%[result])\n"

                                :   [acc0]"+r"(acc0),
                                    [acc1]"+r"(acc1),
                                    [acc2]"+r"(acc2),
                                    [low]"=&r"(low),
                                    [high]"=&r"(high),
                                    [low2]"=&r"(low2),
                                    [high2]"=&r"(high2)
                                :   [result]"r"(result.data()),
                                    [x]"r"(x.data()),
                                    [y]"r"(y.data())
                                :   "rdx", "r11", "cc", "memory"
                            );
                        }

                        inline void multiply_5x5_x86_bmi2_adx(limb_array &result, const limb_array &x,
                                                              const limb_array &y) {
                            limb low, high, low2, high2;
                            limb acc0 = 0;
                            limb acc1 = 0;
                            limb acc2 = 0;

                            asm volatile(
                                "xor %%r11, %%r11\n" // clear zero register

                                // round 1, x0*y0
                                bn254_fp12_multiply_partial_x86(0, 0)
                                bn254_fp12_multiply_emit_x86(0)

                                // round 2, x0 * y1, x1 * y0
                                bn254_fp12_multiply_partial_x86_interleaved2(0, 8, 8, 0)
                                bn254_fp12_multiply_emit_x86(8)

                                // round 3, x0*y2, x1*y1, x2*y0
                                bn254_fp12_multiply_partial_x86_interleaved3(0, 16, 8, 8, 16, 0)
                                bn254_fp12_multiply_emit_x86(16)

                                 // round 4, x0*y3, x1*y2, x2*y1, x3*y0
                                bn254_fp12_multiply_partial_x86_interleaved2(0, 24, 8, 16)
                                bn254_fp12_multiply_partial_x86_interleaved2(16, 8, 24, 0)
                                bn254_fp12_multiply_emit_x86(24)

                                 // round 5, x0*y4, x1*y3, x2*y2, x3*y1, x4*y0
                                bn254_fp12_multiply_partial_x86_interleaved3(0, 32, 8, 24, 16, 16)
                                bn254_fp12_multiply_partial_x86_interleaved2(24, 8, 32, 0)
                                bn254_fp12_multiply_emit_x86(32)

                                 // round 6, x1*y4, x2*y3, x3*y2, x4*y1
                                bn254_fp12_multiply_partial_x86_interleaved2(8, 32, 16, 24)
                                bn254_fp12_multiply_partial_x86_interleaved2(24, 16, 32, 8)
                                bn254_fp12_multiply_emit_x86(40)

                                 // round 7, x2*y4, x3*y3, x4*y2
                                bn254_fp12_multiply_partial_x86_interleaved3(16, 32, 24, 24, 32, 16)
                                bn254_fp12_multiply_emit_x86(48)

                                 // round 8, x3*y4, x4*y3
                                bn254_fp12_multiply_partial_x86_interleaved2(24, 32, 32, 24)
                                bn254_fp12_multiply_emit_x86(56)

                                 // round 9, x4*y4
                                bn254_fp12_multiply_partial_x86(32, 32)
                                bn254_fp12_multiply_emit_x86(64)


                                :   [acc0]"+r"(acc0),
                                    [acc1]"+r"(acc1),
                                    [acc2]"+r"(acc2),
                                    [low]"=&r"(low),
                                    [high]"=&r"(high),
                                    [low2]"=&r"(low2),
                                    [high2]"=&r"(high2)
                                :   [result]"r"(result.data()),
                                    [x]"r"(x.data()),
                                    [y]"r"(y.data())
                                :   "rdx", "r11", "cc", "memory"
                            );
                        }
                    }    // namespace alt_bn128_fp12_limb_ops
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil

#endif