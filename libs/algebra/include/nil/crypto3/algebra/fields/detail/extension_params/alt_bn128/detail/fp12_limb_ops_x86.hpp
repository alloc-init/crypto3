#pragma once

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>
#include <boost/preprocessor.hpp>

#define STR_IMPL(X) #X
#define STR(X) STR_IMPL(X)
#define CAT_IMPL(a, b) a##b
#define CAT(a, b) CAT_IMPL(a, b)

#define BYTE_OFFSET(I) STR(BOOST_PP_MUL(I, 8))
#define BYTE_OFFSET2(I, J) STR(BOOST_PP_MUL(BOOST_PP_ADD(I, J), 8))

#define bn254_fp12_multiply_partial_x86(X, Y)   \
    "movq " BYTE_OFFSET(X) "(%[x]), %%rax\n"    \
    "mulq " BYTE_OFFSET(Y) "(%[y])\n"           \
    "add %%rax, %[acc0]\n"                      \
    "adc %%rdx, %[acc1]\n"                      \
    "adc $0, %[acc2]\n"

#define bn254_fp12_multiply_emit_x86(I)             \
    "movq %[acc0], " BYTE_OFFSET(I) "(%[result])\n" \
    "movq %[acc1], %[acc0]\n"                       \
    "movq %[acc2], %[acc1]\n"                       \
    "xor %[acc2], %[acc2]\n"

#define REG_T_I(I) "%[t" STR(BOOST_PP_MOD(I, 5)) "]"
#define REG_T_IJ(I, J) "%[t" STR(BOOST_PP_MOD(BOOST_PP_ADD(I, J), 5)) "]"
#define VAR_T_AT(I) CAT(t, BOOST_PP_MOD(BOOST_PP_ADD(I, 4), 5))

// multiply bottom limb by m*p and propagate carries
// regs:
//      * rdx: m
//      * rax, rsi: low, high
//      * rcx: row carry
#define bn254_fp12_montgomery_reduce_mul_mp(I, J)           \
    /* compute m * p[j], m is in rdx */                     \
    /* rax = low, rsi = high */                             \
    "mulxq %[p" #J "], %%rax, %%rsi \n"                     \
    /* rcx += data[i+j] // add data to carry accumulator */ \
    "add " REG_T_IJ(I, J) ", %%rcx\n"                       \
    /* add any overflow to high */                          \
    "adc $0, %%rsi\n"                                       \
    /* add carry/accumulator to low */                      \
    "add %%rcx, %%rax\n"                                    \
    /* add overflow to high */                              \
    "adc $0, %%rsi\n"                                       \
    /* data[i,j] = low */                                   \
    "movq %%rax, " REG_T_IJ(I, J) "\n"                      \
    /* carry = high */                                      \
    "movq %%rsi, %%rcx\n"

// main body of loop in montgomery reduce
// regs:
//      * rdx: m
//      * rcx: row carry
//      * dil: window carry - needs to be zeroed outside loop
#define bn254_fp12_montgomery_reduce_cancel_low(I)                  \
    /* m = data[i] * p_dash */                                      \
    "movq %[t" #I "], %%rdx\n"                                      \
    /* dont modify rdx in mul_mp */                                 \
    "imulq %[p_dash], %%rdx\n"                                      \
    /* clear carry */                                               \
    "xor %%rcx, %%rcx\n"                                            \
    /* main loop, multiply limbs by m*p */                          \
    bn254_fp12_montgomery_reduce_mul_mp(I, 0)                       \
    bn254_fp12_montgomery_reduce_mul_mp(I, 1)                       \
    bn254_fp12_montgomery_reduce_mul_mp(I, 2)                       \
    bn254_fp12_montgomery_reduce_mul_mp(I, 3)                       \
    /* load next limb */                                            \
    "movq " BYTE_OFFSET2(I, 5) "(%[data]), " REG_T_I(I) "\n"        \
    /* propagate carries */                                         \
    "add %%rcx, " REG_T_IJ(I, 4) "\n"                               \
    "adcq %%rdi, " REG_T_I(I) "\n"                                  \
    "setc %%dil\n"                                                 

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

                        template <class Field>
                        void montgomery_reduce_x86(limb_array &data) {
                            constexpr auto mod_obj = Field::modulus_params.get_mod_obj();
                            static constexpr limb p0 = limb(mod_obj.get_mod().limbs()[0]);
                            static constexpr limb p1 = limb(mod_obj.get_mod().limbs()[1]);
                            static constexpr limb p2 = limb(mod_obj.get_mod().limbs()[2]);
                            static constexpr limb p3 = limb(mod_obj.get_mod().limbs()[3]);
                            static constexpr limb p_dash = limb(mod_obj.get_p_dash());

                            limb t0 = data[0];
                            limb t1 = data[1];
                            limb t2 = data[2];
                            limb t3 = data[3];
                            limb t4 = data[4];

                            asm volatile(
                                // initial loop: for each limb, compute m and multiply each limb by m*p
                                // make sure window carry is initialized
                                "xor %%rdi, %%rdi\n"
                                bn254_fp12_montgomery_reduce_cancel_low(0)
                                bn254_fp12_montgomery_reduce_cancel_low(1)
                                bn254_fp12_montgomery_reduce_cancel_low(2)
                                bn254_fp12_montgomery_reduce_cancel_low(3)

                                // subtract modulus
                                "modulus_loop_start%=:\n"
                                "testq " REG_T_IJ(4, 4) ", " REG_T_IJ(4, 4) "\n"
                                "jnz modulus_loop_subtract%=\n"
                                "cmpq %[p3], " REG_T_IJ(3, 4) "\n"
                                "ja modulus_loop_subtract%=\n"
                                "jb modulus_loop_end%=\n"
                                "cmpq %[p2], " REG_T_IJ(2, 4) "\n"
                                "ja modulus_loop_subtract%=\n"
                                "jb modulus_loop_end%=\n"
                                "cmpq %[p1], " REG_T_IJ(1, 4) "\n"
                                "ja modulus_loop_subtract%=\n"
                                "jb modulus_loop_end%=\n"
                                "cmpq %[p0], " REG_T_IJ(0, 4) "\n"
                                "jb modulus_loop_end%=\n"
                                "modulus_loop_subtract%=:\n"
                                "subq %[p0], " REG_T_IJ(0, 4) "\n"
                                "sbbq %[p1], " REG_T_IJ(1, 4) "\n"
                                "sbbq %[p2], " REG_T_IJ(2, 4) "\n"
                                "sbbq %[p3], " REG_T_IJ(3, 4) "\n"
                                "sbbq $0, " REG_T_IJ(4, 4) "\n"
                                "jmp modulus_loop_start%=\n"
                                "modulus_loop_end%=:\n"

                                : [t0]"+r"(t0),
                                  [t1]"+r"(t1),
                                  [t2]"+r"(t2),
                                  [t3]"+r"(t3),
                                  [t4]"+r"(t4)
                                : [data]"r"(data.data()),
                                  [p0]"m"(p0),
                                  [p1]"m"(p1),
                                  [p2]"m"(p2),
                                  [p3]"m"(p3),
                                  [p_dash]"m"(p_dash)
                                : "rax", "rcx", "rdx", "rsi", "rdi", "cc", "memory"
                            );

                            data[0] = VAR_T_AT(0);
                            data[1] = VAR_T_AT(1);
                            data[2] = VAR_T_AT(2);
                            data[3] = VAR_T_AT(3);
                            data[4] = 0;
                            data[5] = 0;
                            data[6] = 0;
                            data[7] = 0;
                            data[8] = 0;
                        }

                    }    // namespace alt_bn128_fp12_limb_ops
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil
