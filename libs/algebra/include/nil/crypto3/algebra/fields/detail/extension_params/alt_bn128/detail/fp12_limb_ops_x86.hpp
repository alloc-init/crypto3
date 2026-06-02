#pragma once

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>
#include <boost/preprocessor.hpp>

#define STR_IMPL(X) #X
#define STR(X) STR_IMPL(X)
#define BYTE_OFFSET(I) STR(BOOST_PP_MUL(I, 8))
#define BYTE_OFFSET2(I, J) STR(BOOST_PP_ADD(BOOST_PP_MUL(I, 8), BOOST_PP_MUL(J, 8)))

#define bn254_fp12_multiply_partial_x86(X, Y)   \
    "movq " BYTE_OFFSET(X) "(%[x]), %%rax\n"                      \
    "mulq " BYTE_OFFSET(Y) "(%[y])\n"                             \
    "add %%rax, %[acc0]\n"                                  \
    "adc %%rdx, %[acc1]\n"                                  \
    "adc $0, %[acc2]\n"

#define bn254_fp12_multiply_emit_x86(I) \
    "movq %[acc0], " BYTE_OFFSET(I) "(%[result])\n" \
    "movq %[acc1], %[acc0]\n"               \
    "movq %[acc2], %[acc1]\n"               \
    "xor %[acc2], %[acc2]\n"

#define bn254_fp12_montgomery_reduce_mul_mp(I, J) \
    /* rax, rdx = m * p[j] */ \
    "movq " BYTE_OFFSET(J) "(%[p]), %%rdx\n" \
    "mulq %[m]\n" \
    /* rcx += data[i+j] // add data to carry accumulator */ \
    "add " BYTE_OFFSET2(I, J) "(%[data]), %%rcx\n" \
    /* add any overflow to high */ \
    "adc $0, %%rdx\n" \
    /* add carry/accumulator to low */ \
    "add %%rcx, %%rax\n" \
    /* add overflow to high */ \
    "adc $0, %%rdx\n" \
    /* data[i,j] = low */ \
    "movq %%rax, " BYTE_OFFSET2(I, J) "(%[data])\n" \
    /* carry = high */ \
    "movq %%rdx, %%rcx\n"

// main body of loop in montgomery reduce
#define bn254_fp12_montgomery_reduce_cancel_low(I) \
    /* m = data[i] * p_dash */ \
    "movq " BYTE_OFFSET(I) "(%[data]), %%rax\n" \
    "mulq %[p_dash]\n" \
    "movq %%rax, %[m]\n" \
    /* main loop, multiply limbs by m*p */ \
    "xor %%rcx, %%rcx\n" /* clear carry */ \
    bn254_fp12_montgomery_reduce_mul_mp(I, 0) \
    bn254_fp12_montgomery_reduce_mul_mp(I, 1) \
    bn254_fp12_montgomery_reduce_mul_mp(I, 2) \
    bn254_fp12_montgomery_reduce_mul_mp(I, 3) \
    "shr $1, %%rcx\n" \
    /* propagate carry to higher limbs */ \
    "adcq $0, " BYTE_OFFSET(4) "(%[data])\n" \
    "adcq $0, " BYTE_OFFSET(5) "(%[data])\n" \
    "adcq $0, " BYTE_OFFSET(6) "(%[data])\n" \
    "adcq $0, " BYTE_OFFSET(7) "(%[data])\n"


#define bn254_fp12_montgomery_reduce_modulus_loop_cmp(I) \
    "movq " BYTE_OFFSET(I) "(%[p]), %%rax\n" \
    "cmpq " BYTE_OFFSET2(I, 4) "(%[data]), %%rax\n" \
    "ja modulus_loop_subtract\n" \
    "jb modulus_loop_end\n" 

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
                            limb carry = 0;

                            asm volatile(
                                // round 1, x0*y0
                                bn254_fp12_multiply_partial_x86(0, 0)
                                bn254_fp12_multiply_emit_x86(0)

                                // round 2, x0 * y1, x1 * y0
                                bn254_fp12_multiply_partial_x86(0, 1)
                                bn254_fp12_multiply_partial_x86(1, 0)
                                bn254_fp12_multiply_emit_x86(1)

                                // round 3, x0*y2, x1*y1, x2*y0
                                bn254_fp12_multiply_partial_x86(0,2)
                                bn254_fp12_multiply_partial_x86(1,1)
                                bn254_fp12_multiply_partial_x86(2,0)
                                bn254_fp12_multiply_emit_x86(2)

                                 // round 4, x0*y3, x1*y2, x2*y1, x3*y0
                                bn254_fp12_multiply_partial_x86(0, 3)
                                bn254_fp12_multiply_partial_x86(1, 2)
                                bn254_fp12_multiply_partial_x86(2, 1)
                                bn254_fp12_multiply_partial_x86(3, 0)
                                bn254_fp12_multiply_emit_x86(3)

                                 // round 5, x1*y3, x2*y2, x3*y1
                                bn254_fp12_multiply_partial_x86(1, 3)
                                bn254_fp12_multiply_partial_x86(2, 2)
                                bn254_fp12_multiply_partial_x86(3, 1)
                                bn254_fp12_multiply_emit_x86(4)

                                 // round 6, x2*y3, x3*y2
                                bn254_fp12_multiply_partial_x86(2, 3)
                                bn254_fp12_multiply_partial_x86(3, 2)
                                bn254_fp12_multiply_emit_x86(5)

                                 // round 7, x3*y3
                                bn254_fp12_multiply_partial_x86(3, 3)
                                bn254_fp12_multiply_emit_x86(6)
                                "movq %[acc0], " BYTE_OFFSET(7) "(%[result])\n"
                                "movq %[acc1], " BYTE_OFFSET(8) "(%[result])\n"


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
                                bn254_fp12_multiply_partial_x86(0, 1)
                                bn254_fp12_multiply_partial_x86(1, 0)
                                bn254_fp12_multiply_emit_x86(1)

                                // round 3, x0*y2, x1*y1, x2*y0
                                bn254_fp12_multiply_partial_x86(0, 2)
                                bn254_fp12_multiply_partial_x86(1, 1)
                                bn254_fp12_multiply_partial_x86(2, 0)
                                bn254_fp12_multiply_emit_x86(2)

                                 // round 4, x0*y3, x1*y2, x2*y1, x3*y0
                                bn254_fp12_multiply_partial_x86(0, 3)
                                bn254_fp12_multiply_partial_x86(1, 2)
                                bn254_fp12_multiply_partial_x86(2, 1)
                                bn254_fp12_multiply_partial_x86(3, 0)
                                bn254_fp12_multiply_emit_x86(3)

                                 // round 5, x0*y4, x1*y3, x2*y2, x3*y1, x4*y0
                                bn254_fp12_multiply_partial_x86(0, 4)
                                bn254_fp12_multiply_partial_x86(1, 3)
                                bn254_fp12_multiply_partial_x86(2, 2)
                                bn254_fp12_multiply_partial_x86(3, 1)
                                bn254_fp12_multiply_partial_x86(4, 0)
                                bn254_fp12_multiply_emit_x86(4)

                                 // round 6, x1*y4, x2*y3, x3*y2, x4*y1
                                bn254_fp12_multiply_partial_x86(1, 4)
                                bn254_fp12_multiply_partial_x86(2, 3)
                                bn254_fp12_multiply_partial_x86(3, 2)
                                bn254_fp12_multiply_partial_x86(4, 1)
                                bn254_fp12_multiply_emit_x86(5)

                                 // round 7, x2*y4, x3*y3, x4*y2
                                bn254_fp12_multiply_partial_x86(2, 4)
                                bn254_fp12_multiply_partial_x86(3, 3)
                                bn254_fp12_multiply_partial_x86(4, 2)
                                bn254_fp12_multiply_emit_x86(6)

                                 // round 8, x3*y4, x4*y3
                                bn254_fp12_multiply_partial_x86(3, 4)
                                bn254_fp12_multiply_partial_x86(4, 3)
                                bn254_fp12_multiply_emit_x86(7)

                                 // round 9, x4*y4, drop high bits (see note in portable version)
                                "movq " BYTE_OFFSET(4) "(%[x]), %%rax\n"
                                "mulq " BYTE_OFFSET(4) "(%[y])\n"
                                "add %%rax, %[acc0]\n"
                                "movq %[acc0], " BYTE_OFFSET(8) "(%[result])\n"

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
                            limb m;

                            asm volatile(
                                // initial loop: for each limb, compute m and multiply each limb by m*p
                                bn254_fp12_montgomery_reduce_cancel_low(0)
                                bn254_fp12_montgomery_reduce_cancel_low(1)
                                bn254_fp12_montgomery_reduce_cancel_low(2)
                                bn254_fp12_montgomery_reduce_cancel_low(3)

                                // subtract modulus
                                "modulus_loop_start:\n"
                                // data[9] is nonzero, its definitely greater than p
                                "cmpq $0, " BYTE_OFFSET(9) "(%[data])\n"
                                "jne modulus_loop_subtract\n"
                                bn254_fp12_montgomery_reduce_modulus_loop_cmp(3)
                                bn254_fp12_montgomery_reduce_modulus_loop_cmp(2)
                                bn254_fp12_montgomery_reduce_modulus_loop_cmp(1)
                                bn254_fp12_montgomery_reduce_modulus_loop_cmp(0)
                                "jmp modulus_loop_end\n"
                                "modulus_loop_subtract:\n"
                                "movq " BYTE_OFFSET(0) "(%[p]), %%rax\n"
                                "subq %%rax, " BYTE_OFFSET(4) "(%[data])\n"
                                "movq " BYTE_OFFSET(1) "(%[p]), %%rax\n"
                                "sbbq %%rax, " BYTE_OFFSET(5) "(%[data])\n"
                                "movq " BYTE_OFFSET(2) "(%[p]), %%rax\n"
                                "sbbq %%rax, " BYTE_OFFSET(6) "(%[data])\n"
                                "movq " BYTE_OFFSET(3) "(%[p]), %%rax\n"
                                "sbbq %%rax, " BYTE_OFFSET(7) "(%[data])\n"
                                "sbbq $0, " BYTE_OFFSET(8) "(%[data])\n"
                                "jmp modulus_loop_start\n"
                                "modulus_loop_end:\n"

                                "movq " BYTE_OFFSET(4) "(%[data]), %%rax\n"
                                "movq %%rax, " BYTE_OFFSET(0) "(%[data])\n"
                                "movq " BYTE_OFFSET(5) "(%[data]), %%rax\n"
                                "movq %%rax, " BYTE_OFFSET(1) "(%[data])\n"
                                "movq " BYTE_OFFSET(6) "(%[data]), %%rax\n"
                                "movq %%rax, " BYTE_OFFSET(2) "(%[data])\n"
                                "movq " BYTE_OFFSET(7) "(%[data]), %%rax\n"
                                "movq %%rax, " BYTE_OFFSET(3) "(%[data])\n"
                                "movq $0, " BYTE_OFFSET(4) "(%[data])\n"
                                "movq $0, " BYTE_OFFSET(5) "(%[data])\n"
                                "movq $0, " BYTE_OFFSET(6) "(%[data])\n"
                                "movq $0, " BYTE_OFFSET(7) "(%[data])\n"
                                "movq $0, " BYTE_OFFSET(8) "(%[data])\n"

                                : [m]"=&r"(m)
                                : [data]"r"(data.data()), [p]"r"(p.data()), [p_dash]"r"(p_dash)
                                : "rax", "rcx", "rdx", "cc", "memory"
                            );
                        }

                    }    // namespace alt_bn128_fp12_limb_ops
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil