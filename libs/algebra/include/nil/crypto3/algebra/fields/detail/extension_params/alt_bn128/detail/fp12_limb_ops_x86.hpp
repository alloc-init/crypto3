#pragma once

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>
#include <boost/preprocessor.hpp>

// clang-format off

#define STR_IMPL(X) #X
#define STR(X) STR_IMPL(X)

#define PTR(REGNAME, I) STR(BOOST_PP_MUL(I, 8)) "(%[" #REGNAME "])"
#define PTR2(REGNAME, I, J) PTR(REGNAME, BOOST_PP_ADD(I, J))

// get the i+j%5-th "d" register
#define D(I, J) "%[d" STR(BOOST_PP_MOD(BOOST_PP_ADD(I, J), 5)) "]"

#define schoolbook_round(I)                     \
    "mov " PTR(y, I) ", %%rdx\n"                \
    "mulx " PTR(x, 0) ", %%rax, %[high]\n"      \
    "adox %%rax, " D(0, I) "\n"                 \
    "adcx %[high], " D(1, I) "\n"               \
    "mulx " PTR(x, 1) ", %%rax, %[high]\n"      \
    "adox %%rax, " D(1, I) "\n"                 \
    "adcx %[high], " D(2, I) "\n"               \
    "mulx " PTR(x, 2) ", %%rax, %[high]\n"      \
    "adox %%rax, " D(2, I) "\n"                 \
    "adcx %[high], " D(3, I) "\n"               \
    "mulx " PTR(x, 3) ", %%rax, %[high]\n"      \
    "adox %%rax, " D(3, I) "\n"                 \
    "adcx %[high], " D(4, I) "\n"               \
    "adox %[zero], " D(4, I) "\n"               \
    "adc %[zero], " D(4, I) "\n"                \
    "mov " D(0, I) ", " PTR(result, I) "\n"     \
    "xor " D(0, I) ", " D(0, I) "\n"

// clang-format on
namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    inline void multiply_4x4_x86(limb_array &result, const limb_array &x, const limb_array &y) {
        limb zero, high;
        limb d0 = 0, d1 = 0, d2 = 0, d3 = 0, d4 = 0;

        asm volatile(
            // initial round
            "mov " PTR(y, 0) ", %%rdx\n"
            "mulx " PTR(x, 0) ", %[d0], %[d1]\n"
            "mulx " PTR(x, 1) ", %%rax, %[d2]\n"
            "add %%rax, %[d1]\n"
            "mulx " PTR(x, 2) ", %%rax, %[d3]\n"
            "adc %%rax, %[d2]\n"
            "mulx " PTR(x, 3) ", %%rax, %[d4]\n"
            "adc %%rax, %[d3]\n"
            "adc $0, %[d4]\n"
            "mov %[d0], " PTR(result, 0) "\n"
            "xor %[d0], %[d0]\n"
            "xor %[zero], %[zero]\n"

            schoolbook_round(1)
            schoolbook_round(2)
            schoolbook_round(3)

            "mov " D(1, 3) ", " PTR(result, 4) "\n"
            "mov " D(2, 3) ", " PTR(result, 5) "\n"
            "mov " D(3, 3) ", " PTR(result, 6) "\n"
            "mov " D(4, 3) ", " PTR(result, 7) "\n"

            : [high]"=&r"(high),
              [zero]"=&r"(zero),
              [d0]"=&r"(d0),
              [d1]"=&r"(d1),
              [d2]"=&r"(d2),
              [d3]"=&r"(d3),
              [d4]"=&r"(d4)
            : [result]"r"(result.data()),
              [y]"r"(y.data()),
              [x]"r"(x.data())
            : "rax", "rdx", "cc", "memory"
        );
    }

    inline void multiply_5x5_x86(limb_array &result, const limb_array &x, const limb_array &y) {
        multiply_4x4_x86(result, x, y);

        limb d0, d1, d2, d3;
        asm volatile(
            // round 8, x4 and y4 can only be carries
            "cmpq $0, " PTR(x, 4) "\n"
            "je done_adding_y%=\n"
            "movq " PTR(y, 0) ", %[d0]\n"
            "movq " PTR(y, 1) ", %[d1]\n"
            "movq " PTR(y, 2) ", %[d2]\n"
            "movq " PTR(y, 3) ", %[d3]\n"
            "addq %[d0], " PTR(result, 4) "\n"
            "adcq %[d1], " PTR(result, 5) "\n"
            "adcq %[d2], " PTR(result, 6) "\n"
            "adcq %[d3], " PTR(result, 7) "\n"
            "adcq $0, " PTR(result, 8) "\n"
            "done_adding_y%=:\n"

            "cmpq $0, " PTR(y, 4) "\n"
            "je done_adding_x%=\n"
            "movq " PTR(x, 0) ", %[d0]\n"
            "movq " PTR(x, 1) ", %[d1]\n"
            "movq " PTR(x, 2) ", %[d2]\n"
            "movq " PTR(x, 3) ", %[d3]\n"
            "addq %[d0], " PTR(result, 4) "\n"
            "adcq %[d1], " PTR(result, 5) "\n"
            "adcq %[d2], " PTR(result, 6) "\n"
            "adcq %[d3], " PTR(result, 7) "\n"
            "adcq $0, " PTR(result, 8) "\n"
            "done_adding_x%=:\n"

            "movq " PTR(x, 4) ", %[d0]\n"
            "andq " PTR(y, 4) ", %[d0]\n"
            "addq %[d0], " PTR(result, 8) "\n"

            : [d0]"=&r"(d0),
              [d1]"=&r"(d1),
              [d2]"=&r"(d2),
              [d3]"=&r"(d3)
            : [result]"r"(result.data()),
              [x]"r"(x.data()),
              [y]"r"(y.data())
            : "cc", "memory"
        );
    }
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops

// get the i+j%5-th "t" register
#define T(I, J) "%[t" STR(BOOST_PP_MOD(BOOST_PP_ADD(I, J), 5)) "]"

// clang-format off

// multiply bottom limb by m*p and propagate carries
// HIGH and CARRY are parameterized so you can alternate them and avoid a copy
#define montgomery_reduce_mul_mp(I, J, HIGH, CARRY)             \
    /* compute m * p[j], m is in rdx */                         \
    "mulx %[p" #J "], %[low], %[" #HIGH "] \n"                  \
    "add %[" #CARRY "], " T(I, J) "\n"                          \
    "adc $0, %[" #HIGH "]\n"                                    \
    "add %[low], " T(I, J) "\n"                                 \
    "adc $0, %[" #HIGH "]\n"

// main body of loop in montgomery reduce
#define montgomery_reduce_cancel_low(I)                                 \
    /* m = data[i] * p_dash */                                          \
    "movq %[t" #I "], %%rdx\n"                                          \
    "imulq %[p_dash], %%rdx\n"                                          \
    /* first iteration - avoid generic loop body w carry */             \
    "mulxq %[p0], %[low], %[carry]\n"                                   \
    /* add just to see if carry */                                      \
    "add %[t" #I "], %[low]\n"                                          \
    /* add overflow to carry */                                         \
    "adc $0, %[carry]\n"                                                \
    /* dont have to store low since it t[i] canceled this round */      \
    montgomery_reduce_mul_mp(I, 1, high, carry)                         \
    montgomery_reduce_mul_mp(I, 2, carry, high)                         \
    montgomery_reduce_mul_mp(I, 3, high, carry)                         \
    /* load next limb */                                                \
    "movq " PTR2(data, I, 5) ", %[t" #I "]\n"                           \
    "add %[high], " T(I, 4) "\n"                                        \
    "adc %[pending], %[t" #I "]\n"                                      \
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
            montgomery_reduce_cancel_low(0)
            montgomery_reduce_cancel_low(1)
            montgomery_reduce_cancel_low(2)
            montgomery_reduce_cancel_low(3)

            "modulus%=:\n"
            // try t =- p until we get a carry flag
            "subq %[p0], " T(0, 4) "\n"
            "sbbq %[p1], " T(1, 4) "\n"
            "sbbq %[p2], " T(2, 4) "\n"
            "sbbq %[p3], " T(3, 4) "\n"
            "sbbq $0, " T(4, 4) "\n"
            "jnc modulus%=\n"

            // add back 1 p
            "add %[p0], " T(0, 4) "\n"
            "adc %[p1], " T(1, 4) "\n"
            "adc %[p2], " T(2, 4) "\n"
            "adc %[p3], " T(3, 4) "\n"

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

    template<class Field>
    inline void add_8_limbs_mod_x86(limb_array &data, const limb_array &other) {
        constexpr auto mod_obj = Field::modulus_params.get_mod_obj();
        static constexpr limb p0 = limb(mod_obj.get_mod().limbs()[0]);
        static constexpr limb p1 = limb(mod_obj.get_mod().limbs()[1]);
        static constexpr limb p2 = limb(mod_obj.get_mod().limbs()[2]);
        static constexpr limb p3 = limb(mod_obj.get_mod().limbs()[3]);
        limb t0 = data[4];
        limb t1 = data[5];
        limb t2 = data[6];
        limb t3 = data[7];
        limb q0, q1, q2, q3;
        bool cond;
        asm volatile(
            "movq " PTR(data, 0) ", %%rax\n"
            "addq " PTR(other, 0) ", %%rax\n"
            "movq %%rax, " PTR(data, 0) "\n"
            "movq " PTR(data, 1) ", %%rax\n"
            "adcq " PTR(other, 1) ", %%rax\n"
            "movq %%rax, " PTR(data, 1) "\n"
            "movq " PTR(data, 2) ", %%rax\n"
            "adcq " PTR(other, 2) ", %%rax\n"
            "movq %%rax, " PTR(data, 2) "\n"
            "movq " PTR(data, 3) ", %%rax\n"
            "adcq " PTR(other, 3) ", %%rax\n"
            "movq %%rax, " PTR(data, 3) "\n"
            "movq " PTR(other, 4) ", %%rax\n"
            "adcq %%rax, %[t0]\n"
            "movq " PTR(other, 5) ", %%rax\n"
            "adcq %%rax, %[t1]\n"
            "movq " PTR(other, 6) ", %%rax\n"
            "adcq %%rax, %[t2]\n"
            "movq " PTR(other, 7) ", %%rax\n"
            "adcq %%rax, %[t3]\n"

            // q = t
            "movq %[t0], %[q0]\n"
            "movq %[t1], %[q1]\n"
            "movq %[t2], %[q2]\n"
            "movq %[t3], %[q3]\n"

            // try q - p
            "subq %[p0], %[q0]\n"
            "sbbq %[p1], %[q1]\n"
            "sbbq %[p2], %[q2]\n"
            "sbbq %[p3], %[q3]\n"

            // if q - p didnt result in a borrrow, set t=q
            "cmovnc %[q0], %[t0]\n"
            "cmovnc %[q1], %[t1]\n"
            "cmovnc %[q2], %[t2]\n"
            "cmovnc %[q3], %[t3]\n"

            : [t0]"+r"(t0),
              [t1]"+r"(t1),
              [t2]"+r"(t2),
              [t3]"+r"(t3),
              [q0]"=&r"(q0),
              [q1]"=&r"(q1),
              [q2]"=&r"(q2),
              [q3]"=&r"(q3)
            : [data]"r"(data.data()),
              [other]"r"(other.data()),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3)
            : "rax", "cc", "memory"
        );
        data[4] = t0;
        data[5] = t1;
        data[6] = t2;
        data[7] = t3;
    }
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops

#undef STR_IMPL
#undef STR
#undef PTR
#undef PTR2
#undef multiply_partial
#undef multiply_emit
#undef D
#undef T
#undef montgomery_reduce_mul_mp
#undef montgomery_reduce_cancel_low