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
    "adcx %[zero], " D(4, I) "\n"                \
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
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops

// clang-format off

#define set_static_modulus_limbs_from_field()                           \
        constexpr auto mod_obj = Field::modulus_params.get_mod_obj();   \
        static constexpr limb p0 = limb(mod_obj.get_mod().limbs()[0]);  \
        static constexpr limb p1 = limb(mod_obj.get_mod().limbs()[1]);  \
        static constexpr limb p2 = limb(mod_obj.get_mod().limbs()[2]);  \
        static constexpr limb p3 = limb(mod_obj.get_mod().limbs()[3]);  \
        static constexpr limb p_dash = limb(mod_obj.get_p_dash());

// get the i+j%5-th "t" register
#define T(I, J) "%[t" STR(BOOST_PP_MOD(BOOST_PP_ADD(I, J), 5)) "]"

#define montgomery_reduce_load_next(I)                  \
    "setc %b[pending]\n"                                \
    "movq " PTR2(data, I, 5) ", " T(I, 5) "\n"

// main body of loop in montgomery reduce
#define montgomery_reduce_cancel_low(I)                 \
    /* m = data[i] * p_dash */                          \
    "movq " T(I, 0) ", %%rdx\n"                         \
    "imulq %[p_dash], %%rdx\n"                          \
    "xor %[zero], %[zero]\n"                            \
    /* multiply m * pi for each i and add it to data */ \
    "mulxq %[p0], %[low], %[high]\n"                    \
    /* use dual carry chains */                         \
    "adcx %[low], " T(I, 0) "\n"                        \
    "adox %[high], " T(I, 1) "\n"                       \
    "mulxq %[p1], %[low], %[high]\n"                    \
    "adcx %[low], " T(I, 1) "\n"                        \
    "adox %[high], " T(I, 2) "\n"                       \
    "mulxq %[p2], %[low], %[high]\n"                    \
    "adcx %[low], " T(I, 2) "\n"                        \
    "adox %[high], " T(I, 3) "\n"                       \
    "mulxq %[p3], %[low], %[high]\n"                    \
    "adcx %[low], " T(I, 3) "\n"                        \
    /* merge carry chains */                            \
    "adox %[zero], %[high]\n"                           \
    "adox %[pending], %[high]\n"                        \
    "adcx %[high], " T(I, 4) "\n"                       \
    BOOST_PP_IF(BOOST_PP_LESS(I, 3), montgomery_reduce_load_next(I), "")

// clang-format on

namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    template<class Field>
    inline void montgomery_reduce_x86(limb_array &data) {
        set_static_modulus_limbs_from_field();

        limb t0, t1, t2, t3, t4;
        limb low, high, pending, zero;

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

            // Reduce high limbs mod p, reuse some registers for it
            // q = t
            "movq " T(0, 4) ", %[low]\n"
            "movq " T(1, 4) ", %[high]\n"
            "movq " T(2, 4) ", %[zero]\n"
            "movq " T(3, 4) ", %%rdx\n"

            // try q - p
            "subq %[p0], %[low]\n"
            "sbbq %[p1], %[high]\n"
            "sbbq %[p2], %[zero]\n"
            "sbbq %[p3], %%rdx\n"

            // if q - p didnt result in a borrrow, set t=q-p
            "cmovnc %[low], " T(0, 4) "\n"
            "cmovnc %[high], " T(1, 4) "\n"
            "cmovnc %[zero], " T(2, 4) "\n"
            "cmovnc %%rdx, " T(3, 4) "\n"

            "movq " T(0, 4) ", " PTR(data, 0) "\n"
            "movq " T(1, 4) ", " PTR(data, 1) "\n"
            "movq " T(2, 4) ", " PTR(data, 2) "\n"
            "movq " T(3, 4) ", " PTR(data, 3) "\n"
            "movq $0, " PTR(data, 4) "\n"
            "movq $0, " PTR(data, 5) "\n"
            "movq $0, " PTR(data, 6) "\n"
            "movq $0, " PTR(data, 7) "\n"

            : [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3),
              [t4]"=&r"(t4),
              [low]"=&r"(low),
              [high]"=&r"(high),
              [pending]"=&r"(pending),
              [zero]"=&r"(zero)
            : [data]"r"(data.data()),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3),
              [p_dash]"m"(p_dash)
            : "rdx", "cc", "memory"
        );
    }
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops

namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    inline void add_8_limbs_x86(limb *result, const limb *other) {
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
            :
            : [result]"r"(result),
              [other]"r"(other)
            : "rax", "cc", "memory"
        );
    }

    template<class Field>
    inline void add_low_4_limbs_mod_x86(limb *data, const limb *other) {
        set_static_modulus_limbs_from_field();
        limb t0 = data[0];
        limb t1 = data[1];
        limb t2 = data[2];
        limb t3 = data[3];
        limb q0, q1, q2, q3;
        bool cond;
        asm volatile(
            "movq " PTR(other, 0) ", %%rax\n"
            "addq %%rax, %[t0]\n"
            "movq " PTR(other, 1) ", %%rax\n"
            "adcq %%rax, %[t1]\n"
            "movq " PTR(other, 2) ", %%rax\n"
            "adcq %%rax, %[t2]\n"
            "movq " PTR(other, 3) ", %%rax\n"
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
            : [other]"r"(other),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3)
            : "rax", "cc", "memory"
        );
        data[0] = t0;
        data[1] = t1;
        data[2] = t2;
        data[3] = t3;
    }

    inline void subtract_8_limbs_x86(limb_array &result, const limb_array &other) {
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
            :
            : [result]"r"(result.data()),
              [other]"r"(other.data())
            : "rax", "cc", "memory"
        );
    }

    template<class Field>
    inline void add_8_limbs_mod_x86(limb_array &data, const limb_array &other) {
        set_static_modulus_limbs_from_field();
        limb t0 = data[4];
        limb t1 = data[5];
        limb t2 = data[6];
        limb t3 = data[7];
        limb q0, q1, q2, q3;
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
            "movq %[t0], %[q0]\n"
            "movq %[t1], %[q1]\n"
            "movq %[t2], %[q2]\n"
            "movq %[t3], %[q3]\n"
            "subq %[p0], %[q0]\n"
            "sbbq %[p1], %[q1]\n"
            "sbbq %[p2], %[q2]\n"
            "sbbq %[p3], %[q3]\n"
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

    template<class Field>
    inline void subtract_8_limbs_mod_x86(limb_array &result, const limb_array &other) {
        set_static_modulus_limbs_from_field();
        limb t0 = result[4];
        limb t1 = result[5];
        limb t2 = result[6];
        limb t3 = result[7];
        limb q0, q1, q2, q3;
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

            "sbbq " PTR(other, 4) ", %[t0]\n"
            "sbbq " PTR(other, 5) ", %[t1]\n"
            "sbbq " PTR(other, 6) ", %[t2]\n"
            "sbbq " PTR(other, 7) ", %[t3]\n"

            "setc %%al\n" // save carry flag

            "movq %[t0], %[q0]\n"
            "movq %[t1], %[q1]\n"
            "movq %[t2], %[q2]\n"
            "movq %[t3], %[q3]\n"

            "addq %[p0], %[q0]\n"
            "adcq %[p1], %[q1]\n"
            "adcq %[p2], %[q2]\n"
            "adcq %[p3], %[q3]\n"

            "add $255, %%al\n" // restore carry flag

            "cmovc %[q0], %[t0]\n"
            "cmovc %[q1], %[t1]\n"
            "cmovc %[q2], %[t2]\n"
            "cmovc %[q3], %[t3]\n"

            : [t0]"+r"(t0),
              [t1]"+r"(t1),
              [t2]"+r"(t2),
              [t3]"+r"(t3),
              [q0]"=&r"(q0),
              [q1]"=&r"(q1),
              [q2]"=&r"(q2),
              [q3]"=&r"(q3)
            : [result]"r"(result.data()),
              [other]"r"(other.data()),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3)
            : "rax", "cc", "memory"
        );
        result[4] = t0;
        result[5] = t1;
        result[6] = t2;
        result[7] = t3;
    }
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops

#define double_mod_p()      \
    "addq %[t0], %[t0]\n"   \
    "adcq %[t1], %[t1]\n"   \
    "adcq %[t2], %[t2]\n"   \
    "adcq %[t3], %[t3]\n"   \
    "adcq %[t4], %[t4]\n"   \
    "adcq %[t5], %[t5]\n"   \
    "adcq %[t6], %[t6]\n"   \
    "adcq %[t7], %[t7]\n"   \
    "movq %[t4], %[q0]\n"   \
    "movq %[t5], %[q1]\n"   \
    "movq %[t6], %[q2]\n"   \
    "movq %[t7], %[q3]\n"   \
    "subq %[p0], %[q0]\n"   \
    "sbbq %[p1], %[q1]\n"   \
    "sbbq %[p2], %[q2]\n"   \
    "sbbq %[p3], %[q3]\n"   \
    "cmovnc %[q0], %[t4]\n" \
    "cmovnc %[q1], %[t5]\n" \
    "cmovnc %[q2], %[t6]\n" \
    "cmovnc %[q3], %[t7]\n"

namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    template<class Field>
    inline void mul_8_limbs_by_9_x86(limb_array &t) {
        set_static_modulus_limbs_from_field();
        limb_array x(t);
        limb t0 = t[0];
        limb t1 = t[1];
        limb t2 = t[2];
        limb t3 = t[3];
        limb t4 = t[4];
        limb t5 = t[5];
        limb t6 = t[6];
        limb t7 = t[7];
        limb q0, q1, q2, q3;

        asm volatile(
            double_mod_p() // 2x
            double_mod_p() // 4x
            double_mod_p() // 8x

            // 9x, final add by original value
            "addq " PTR(x, 0) ", %[t0]\n"
            "adcq " PTR(x, 1) ", %[t1]\n"
            "adcq " PTR(x, 2) ", %[t2]\n"
            "adcq " PTR(x, 3) ", %[t3]\n"
            "adcq " PTR(x, 4) ", %[t4]\n"
            "adcq " PTR(x, 5) ", %[t5]\n"
            "adcq " PTR(x, 6) ", %[t6]\n"
            "adcq " PTR(x, 7) ", %[t7]\n"
            "movq %[t4], %[q0]\n"
            "movq %[t5], %[q1]\n"
            "movq %[t6], %[q2]\n"
            "movq %[t7], %[q3]\n"
            "subq %[p0], %[q0]\n"
            "sbbq %[p1], %[q1]\n"
            "sbbq %[p2], %[q2]\n"
            "sbbq %[p3], %[q3]\n"
            "cmovnc %[q0], %[t4]\n"
            "cmovnc %[q1], %[t5]\n"
            "cmovnc %[q2], %[t6]\n"
            "cmovnc %[q3], %[t7]\n"

            : [t0]"+r"(t0),
              [t1]"+r"(t1),
              [t2]"+r"(t2),
              [t3]"+r"(t3),
              [t4]"+r"(t4),
              [t5]"+r"(t5),
              [t6]"+r"(t6),
              [t7]"+r"(t7),
              [q0]"=&r"(q0),
              [q1]"=&r"(q1),
              [q2]"=&r"(q2),
              [q3]"=&r"(q3)
            : [x]"r"(x.data()),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3)
            : "cc"
        );
        t[0] = t0;
        t[1] = t1;
        t[2] = t2;
        t[3] = t3;
        t[4] = t4;
        t[5] = t5;
        t[6] = t6;
        t[7] = t7;
    }

    template<class Field>
    inline void fp2_mul_pre_x86(limb *z, const limb *x, const limb *y) {

        // For x = a + bu and y = c + du:
        //   xy = (a + bu) * (c + du)
        //      = ac + adu + bcu + bdu^2
        //      = ac + (ad + bc)u - bd      # since u^2 = -1
        //      = (ac - bd) + (ad + bc)u
        // Karatsuba computes the cross term with one product:
        //   ad + bc = (a + b)(c + d) - ac - bd.
        // limb_array ac;
        // limb_array bd;
        // multiply_4x4_x86(ac, a, c);
    
        // multiply_4x4_x86(bd, b, d);
        // limb_array a_plus_b = a;
        // limb_array c_plus_d = c;

        // add_8_limbs(a_plus_b, b); // adcx
        // add_8_limbs(c_plus_d, d); // adox

        // z0 = ac;
        // subtract_8_limbs_mod_x86<Field>(z0, bd);

        // multiply_4x4(z1, a_plus_b, c_plus_d);
        // subtract_8_limbs(z1, ac);
        // subtract_8_limbs(z1, bd);

        // set_static_modulus_limbs_from_field();
        // limb_array ac, bd;
        // limb_array a_plus_b = a;
        // limb_array c_plus_d = c;
        // asm volatile { : : : "rdx", "cc", "memory"};
    }
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops

#undef STR_IMPL
#undef STR
#undef PTR
#undef PTR2
#undef schoolbook_round
#undef set_static_modulus_limbs_from_field
#undef D
#undef T
#undef montgomery_reduce_load_next
#undef montgomery_reduce_cancel_low
#undef double_mod_p