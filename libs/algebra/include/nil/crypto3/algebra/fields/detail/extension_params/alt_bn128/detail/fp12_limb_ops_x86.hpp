#pragma once

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>
#include <boost/preprocessor.hpp>

// clang-format off

#define STR_IMPL(X) #X
#define STR(X) STR_IMPL(X)

#define PTR(REGNAME, I) STR(BOOST_PP_MUL(I, 8)) "(%[" #REGNAME "])"
#define PTR2(REGNAME, I, J) PTR(REGNAME, BOOST_PP_ADD(I, J))

#define bn254_fp12_multiply_partial_x86(X, Y) \
    "movq " PTR(x, X) ", %%rdx\n"             \
    "mulx " PTR(y, Y) ", %[low], %[high] \n"  \
    "add %[low], %[acc0]\n"                   \
    "adc %[high], %[acc1]\n"                  \
    "adc $0, %[acc2]\n"

#define bn254_fp12_multiply_emit_x86(I)  \
    "movq %[acc0], " PTR(result, I) "\n" \
    "movq %[acc1], %[acc0]\n"            \
    "movq %[acc2], %[acc1]\n"            \
    "xor %[acc2], %[acc2]\n"

// clang-format on
namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    inline void multiply_4x4_x86(limb_array &result, const limb_array &x, const limb_array &y) {
        limb acc0 = 0;
        limb acc1 = 0;
        limb acc2 = 0;
        limb low, high;

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
            "movq %[acc0], " PTR(result, 7) "\n"
            "movq %[acc1], " PTR(result, 8) "\n"

            : [acc0]"+r"(acc0),
              [acc1]"+r"(acc1),
              [acc2]"+r"(acc2),
              [low]"=&r"(low),
              [high]"=&r"(high)
            : [result]"r"(result.data()),
              [x]"r"(x.data()),
              [y]"r"(y.data())
            : "rdx", "cc", "memory"
        );
    }

    inline void multiply_5x5_x86(limb_array &result, const limb_array &x, const limb_array &y) {
        limb acc0 = 0;
        limb acc1 = 0;
        limb acc2 = 0;
        limb low, high;

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
            "movq " PTR(x, 4) ", %%rdx\n"
            "mulx " PTR(y, 4) ", %[low], %[high]\n"
            "add %[low], %[acc0]\n"
            "movq %[acc0], " PTR(result, 8) "\n"

            : [acc0]"+r"(acc0),
              [acc1]"+r"(acc1),
              [acc2]"+r"(acc2),
              [low]"=&r"(low),
              [high]"=&r"(high)
            : [result]"r"(result.data()),
              [x]"r"(x.data()),
              [y]"r"(y.data())
            : "rdx", "cc", "memory"
        );
    }
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops

// get the i+j%5-th "t" register
#define T(I, J) "%[t" STR(BOOST_PP_MOD(BOOST_PP_ADD(I, J), 5)) "]"

// clang-format off

// multiply bottom limb by m*p and propagate carries
// HIGH and CARRY are parameterized so you can alternate them and avoid a copy
#define bn254_fp12_montgomery_reduce_mul_mp(I, J, HIGH, CARRY)  \
    /* compute m * p[j], m is in rdx */                         \
    "mulxq %[p" #J "], %[low], %[" #HIGH "] \n"                 \
    /* carry += data[i+j] // add data to carry accumulator */   \
    "add " T(I, J) ", %[" #CARRY "]\n"                          \
    /* add any overflow to high */                              \
    "adc $0, %[" #HIGH "]\n"                                    \
    /* add carry/accumulator to low */                          \
    "add %[" #CARRY "], %[low]\n"                               \
    /* add overflow to high */                                  \
    "adc $0, %[" #HIGH "]\n"                                    \
    /* data[i,j] = low */                                       \
    "movq %[low], " T(I, J) "\n"

// main body of loop in montgomery reduce
#define bn254_fp12_montgomery_reduce_cancel_low(I)                      \
    /* m = data[i] * p_dash */                                          \
    "movq %[t" #I "], %%rdx\n"                                          \
    /* dont modify rdx in mul_mp */                                     \
    "imulq %[p_dash], %%rdx\n"                                          \
    /* first iteration - avoid generic loop body w carry */             \
    "mulxq %[p0], %[low], %[high]\n"                                    \
    "add %[t" #I "], %[low]\n"                                          \
    /* add overflow to high */                                          \
    "adc $0, %[high]\n"                                                 \
    /* dont have to store low since it t[i] canceled this round */      \
    /* high becomes carry in next round */                              \
    bn254_fp12_montgomery_reduce_mul_mp(I, 1, carry, high)              \
    bn254_fp12_montgomery_reduce_mul_mp(I, 2, high, carry)              \
    bn254_fp12_montgomery_reduce_mul_mp(I, 3, carry, high)              \
    /* load next limb */                                                \
    "movq " PTR2(data, I, 5) ", %[t" #I "]\n"                           \
    "add %[carry], " T(I, 4) "\n"                                       \
    "adcq %[pending], %[t" #I "]\n"                                     \
    "setc %b[pending]\n"

// clang-format on
namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    template<class Field>
    inline void montgomery_reduce_x86(limb_array &data) {
        constexpr auto mod_obj = Field::modulus_params.get_mod_obj();
        static constexpr limb p0 = limb(mod_obj.get_mod().limbs()[0]);
        static constexpr limb p1 = limb(mod_obj.get_mod().limbs()[1]);
        static constexpr limb p2 = limb(mod_obj.get_mod().limbs()[2]);
        static constexpr limb p3 = limb(mod_obj.get_mod().limbs()[3]);
        static constexpr limb p_dash = limb(mod_obj.get_p_dash());

        limb t0, t1, t2, t3, t4, low, high, pending, carry;

        asm volatile(
            "movq " PTR(data, 0) ", %[t0]\n"
            "movq " PTR(data, 1) ", %[t1]\n"
            "movq " PTR(data, 2) ", %[t2]\n"
            "movq " PTR(data, 3) ", %[t3]\n"
            "movq " PTR(data, 4) ", %[t4]\n"

            // initial loop: for each limb, compute m and multiply each limb by m*p
            // make sure window carry is initialized
            "xor %[pending], %[pending]\n"
            bn254_fp12_montgomery_reduce_cancel_low(0)
            bn254_fp12_montgomery_reduce_cancel_low(1)
            bn254_fp12_montgomery_reduce_cancel_low(2)
            bn254_fp12_montgomery_reduce_cancel_low(3)

            // subtract modulus
            "modulus_loop_start%=:\n"
            "testq " T(4, 4) ", " T(4, 4) "\n"
            "jnz modulus_loop_subtract%=\n"
            "cmpq %[p3], " T(3, 4) "\n"
            "ja modulus_loop_subtract%=\n"
            "jb modulus_loop_end%=\n"
            "cmpq %[p2], " T(2, 4) "\n"
            "ja modulus_loop_subtract%=\n"
            "jb modulus_loop_end%=\n"
            "cmpq %[p1], " T(1, 4) "\n"
            "ja modulus_loop_subtract%=\n"
            "jb modulus_loop_end%=\n"
            "cmpq %[p0], " T(0, 4) "\n"
            "jb modulus_loop_end%=\n"
            "modulus_loop_subtract%=:\n"
            "subq %[p0], " T(0, 4) "\n"
            "sbbq %[p1], " T(1, 4) "\n"
            "sbbq %[p2], " T(2, 4) "\n"
            "sbbq %[p3], " T(3, 4) "\n"
            "sbbq $0, " T(4, 4) "\n"
            "jmp modulus_loop_start%=\n"
            "modulus_loop_end%=:\n"

            "movq " T(0, 4) ", " PTR(data, 0) "\n"
            "movq " T(1, 4) ", " PTR(data, 1) "\n"
            "movq " T(2, 4) ", " PTR(data, 2) "\n"
            "movq " T(3, 4) ", " PTR(data, 3) "\n"
            "movq $0, " PTR(data, 4) "\n"
            "movq $0, " PTR(data, 5) "\n"
            "movq $0, " PTR(data, 6) "\n"
            "movq $0, " PTR(data, 7) "\n"
            "movq $0, " PTR(data, 8) "\n"

            : [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3),
              [t4]"=&r"(t4),
              [low]"=&r"(low),
              [high]"=&r"(high),
              [pending]"=&r"(pending),
              [carry]"=&r"(carry)
            : [data]"r"(data.data()),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3),
              [p_dash]"r"(p_dash)
            : "rdx", "cc", "memory"
        );
    }
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops

namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    inline bool subtract_9_limbs_x86(limb *result, const limb *other) {
        bool borrow;
        asm volatile(
            "movq " PTR(result, 0) ", %%rax\n"
            "subq " PTR(other, 0) ", %%rax\n"
            "movq %%rax, " PTR(result, 0) "\n"
            "movq " PTR(result, 1) ", %%rax\n"
            "sbbq " PTR(other, 1) ", %%rax\n"
            "movq %%rax, " PTR(result, 1) "\n"
            "movq " PTR(result, 2) ", %%rax\n"
            "sbbq " PTR(other, 2) ", %%rax\n"
            "movq %%rax, " PTR(result, 2) "\n"
            "movq " PTR(result, 3) ", %%rax\n"
            "sbbq " PTR(other, 3) ", %%rax\n"
            "movq %%rax, " PTR(result, 3) "\n"
            "movq " PTR(result, 4) ", %%rax\n"
            "sbbq " PTR(other, 4) ", %%rax\n"
            "movq %%rax, " PTR(result, 4) "\n"
            "movq " PTR(result, 5) ", %%rax\n"
            "sbbq " PTR(other, 5) ", %%rax\n"
            "movq %%rax, " PTR(result, 5) "\n"
            "movq " PTR(result, 6) ", %%rax\n"
            "sbbq " PTR(other, 6) ", %%rax\n"
            "movq %%rax, " PTR(result, 6) "\n"
            "movq " PTR(result, 7) ", %%rax\n"
            "sbbq " PTR(other, 7) ", %%rax\n"
            "movq %%rax, " PTR(result, 7) "\n"
            "movq " PTR(result, 8) ", %%rax\n"
            "sbbq " PTR(other, 8) ", %%rax\n"
            "movq %%rax, " PTR(result, 8) "\n"
            "setc %[borrow]\n"
            : [borrow]"=r"(borrow)
            : [result]"r"(result),
              [other]"r"(other)
            : "rax", "cc", "memory"
        );
        return borrow;
    }
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops

namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    inline void add_4_limbs_x86(limb *result, const limb *other) {
        asm volatile(
            "movq " PTR(other, 0) ", %%rax\n"
            "addq %%rax, " PTR(result, 0) "\n"
            "movq " PTR(other, 1) ", %%rax\n"
            "adcq %%rax, " PTR(result, 1) "\n"
            "movq " PTR(other, 2) ", %%rax\n"
            "adcq %%rax, " PTR(result, 2) "\n"
            "movq " PTR(other, 3) ", %%rax\n"
            "adcq %%rax, " PTR(result, 3) "\n"
            :
            : [result]"r"(result),
              [other]"r"(other)
            : "rax", "cc", "memory"
        );
    }
    inline void add_9_limbs_x86(limb *result, const limb *other) {
        asm volatile(
            "movq " PTR(other, 0) ", %%rax\n"
            "addq %%rax, " PTR(result, 0) "\n"
            "movq " PTR(other, 1) ", %%rax\n"
            "adcq %%rax, " PTR(result, 1) "\n"
            "movq " PTR(other, 2) ", %%rax\n"
            "adcq %%rax, " PTR(result, 2) "\n"
            "movq " PTR(other, 3) ", %%rax\n"
            "adcq %%rax, " PTR(result, 3) "\n"
            "movq " PTR(other, 4) ", %%rax\n"
            "adcq %%rax, " PTR(result, 4) "\n"
            "movq " PTR(other, 5) ", %%rax\n"
            "adcq %%rax, " PTR(result, 5) "\n"
            "movq " PTR(other, 6) ", %%rax\n"
            "adcq %%rax, " PTR(result, 6) "\n"
            "movq " PTR(other, 7) ", %%rax\n"
            "adcq %%rax, " PTR(result, 7) "\n"
            "movq " PTR(other, 8) ", %%rax\n"
            "adcq %%rax, " PTR(result, 8) "\n"
            :
            : [result]"r"(result),
              [other]"r"(other)
            : "rax", "cc", "memory"
        );
    }
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops

#undef STR_IMPL
#undef STR
#undef PTR
#undef PTR2
#undef bn254_fp12_multiply_partial_x86
#undef bn254_fp12_multiply_emit_x86
#undef T
#undef bn254_fp12_montgomery_reduce_mul_mp
#undef bn254_fp12_montgomery_reduce_cancel_low