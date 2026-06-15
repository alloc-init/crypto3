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

#define initial_schoolbook_round(X, X_BASE, Y, Y_BASE, Z, Z_BASE)   \
    "mov " PTR2(Y, Y_BASE, 0) ", %%rdx\n"                           \
    "mulx " PTR2(X, X_BASE, 0) ", %[d0], %[d1]\n"                   \
    "mulx " PTR2(X, X_BASE, 1) ", %[low], %[d2]\n"                  \
    "add %[low], %[d1]\n"                                           \
    "mulx " PTR2(X, X_BASE, 2) ", %[low], %[d3]\n"                  \
    "adc %[low], %[d2]\n"                                           \
    "mulx " PTR2(X, X_BASE, 3) ", %[low], %[d4]\n"                  \
    "adc %[low], %[d3]\n"                                           \
    "adc $0, %[d4]\n"                                               \
    "mov %[d0], " PTR2(Z, Z_BASE, 0) "\n"                           \
    "xor %[d0], %[d0]\n"                                            \
    "xor %[zero], %[zero]\n"

#define schoolbook_round(ROUND, X, X_BASE, Y, Y_BASE, Z, Z_BASE)    \
    "mov " PTR2(Y, Y_BASE, ROUND) ", %%rdx\n"                       \
    "mulx " PTR2(X, X_BASE, 0) ", %[low], %[high]\n"                \
    "adox %[low], " D(0, ROUND) "\n"                                \
    "adcx %[high], " D(1, ROUND) "\n"                               \
    "mulx " PTR2(X, X_BASE, 1) ", %[low], %[high]\n"                \
    "adox %[low], " D(1, ROUND) "\n"                                \
    "adcx %[high], " D(2, ROUND) "\n"                               \
    "mulx " PTR2(X, X_BASE, 2) ", %[low], %[high]\n"                \
    "adox %[low], " D(2, ROUND) "\n"                                \
    "adcx %[high], " D(3, ROUND) "\n"                               \
    "mulx " PTR2(X, X_BASE, 3) ", %[low], %[high]\n"                \
    "adox %[low], " D(3, ROUND) "\n"                                \
    "adcx %[high], " D(4, ROUND) "\n"                               \
    "adox %[zero], " D(4, ROUND) "\n"                               \
    "adcx %[zero], " D(4, ROUND) "\n"                               \
    "mov " D(0, ROUND) ", " PTR2(Z, Z_BASE, ROUND) "\n"             \
    "xor " D(0, ROUND) ", " D(0, ROUND) "\n"

#define final_schoolbook_round(Z, Z_BASE)                           \
    "mov " D(1, 3) ", " PTR2(Z, Z_BASE, 4) "\n"                     \
    "mov " D(2, 3) ", " PTR2(Z, Z_BASE, 5) "\n"                     \
    "mov " D(3, 3) ", " PTR2(Z, Z_BASE, 6) "\n"                     \
    "mov " D(4, 3) ", " PTR2(Z, Z_BASE, 7) "\n"

#define schoolbook(X, X_BASE, Y, Y_BASE, Z, Z_BASE)                 \
    initial_schoolbook_round(X, X_BASE, Y, Y_BASE, Z, Z_BASE)       \
    schoolbook_round(1, X, X_BASE, Y, Y_BASE, Z, Z_BASE)            \
    schoolbook_round(2, X, X_BASE, Y, Y_BASE, Z, Z_BASE)            \
    schoolbook_round(3, X, X_BASE, Y, Y_BASE, Z, Z_BASE)            \
    final_schoolbook_round(Z, Z_BASE)

// clang-format on
namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    inline void multiply_4x4_x86(limb_array &result, const limb_array &x, const limb_array &y) {
        limb low, zero, high;
        limb d0, d1, d2, d3, d4;

        asm volatile(
            schoolbook(x, 0, y, 0, result, 0)
            : [low]"=&r"(low),
              [high]"=&r"(high),
              [zero]"=&r"(zero),
              [d0]"=&r"(d0),
              [d1]"=&r"(d1),
              [d2]"=&r"(d2),
              [d3]"=&r"(d3),
              [d4]"=&r"(d4)
            : [result]"r"(result.data()),
              [y]"r"(y.data()),
              [x]"r"(x.data())
            : "rdx", "cc", "memory"
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
}

#define subtract_mod(RESULT, RESULT_BASE, OTHER, OTHER_BASE, SCRATCH, T0, T1, T2, T3, Q0, Q1, Q2, Q3) \
    "movq " PTR2(RESULT, RESULT_BASE, 0) ", %[" #SCRATCH "]\n"                                        \
    "subq " PTR2(OTHER, OTHER_BASE, 0) ", %[" #SCRATCH "]\n"                                          \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 0) "\n"                                        \
    "movq " PTR2(RESULT, RESULT_BASE, 1) ", %[" #SCRATCH "]\n"                                        \
    "sbbq " PTR2(OTHER, OTHER_BASE, 1) ", %[" #SCRATCH "]\n"                                          \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 1) "\n"                                        \
    "movq " PTR2(RESULT, RESULT_BASE, 2) ", %[" #SCRATCH "]\n"                                        \
    "sbbq " PTR2(OTHER, OTHER_BASE, 2) ", %[" #SCRATCH "]\n"                                          \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 2) "\n"                                        \
    "movq " PTR2(RESULT, RESULT_BASE, 3) ", %[" #SCRATCH "]\n"                                        \
    "sbbq " PTR2(OTHER, OTHER_BASE, 3) ", %[" #SCRATCH "]\n"                                          \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 3) "\n"                                        \
    "movq " PTR2(RESULT, RESULT_BASE, 4) ", %[" #T0 "]\n"                                             \
    "movq " PTR2(RESULT, RESULT_BASE, 5) ", %[" #T1 "]\n"                                             \
    "movq " PTR2(RESULT, RESULT_BASE, 6) ", %[" #T2 "]\n"                                             \
    "movq " PTR2(RESULT, RESULT_BASE, 7) ", %[" #T3 "]\n"                                             \
    "sbbq " PTR2(OTHER, OTHER_BASE, 4) ", %[" #T0 "]\n"                                               \
    "sbbq " PTR2(OTHER, OTHER_BASE, 5) ", %[" #T1 "]\n"                                               \
    "sbbq " PTR2(OTHER, OTHER_BASE, 6) ", %[" #T2 "]\n"                                               \
    "sbbq " PTR2(OTHER, OTHER_BASE, 7) ", %[" #T3 "]\n"                                               \
    "setc %b[" #SCRATCH "]\n" /* save carry flag */                                                   \
    "movq %[" #T0 "], %[" #Q0 "]\n"                                                                   \
    "movq %[" #T1 "], %[" #Q1 "]\n"                                                                   \
    "movq %[" #T2 "], %[" #Q2 "]\n"                                                                   \
    "movq %[" #T3 "], %[" #Q3 "]\n"                                                                   \
    "addq %[p0], %[" #Q0 "]\n"                                                                        \
    "adcq %[p1], %[" #Q1 "]\n"                                                                        \
    "adcq %[p2], %[" #Q2 "]\n"                                                                        \
    "adcq %[p3], %[" #Q3 "]\n"                                                                        \
    "add $255, %b[" #SCRATCH "]\n" /* restore carry flag */                                           \
    "cmovc %[" #Q0 "], %[" #T0 "]\n"                                                                  \
    "cmovc %[" #Q1 "], %[" #T1 "]\n"                                                                  \
    "cmovc %[" #Q2 "], %[" #T2 "]\n"                                                                  \
    "cmovc %[" #Q3 "], %[" #T3 "]\n"                                                                  \
    "movq %[" #T0 "], " PTR2(RESULT, RESULT_BASE, 4) "\n"                                             \
    "movq %[" #T1 "], " PTR2(RESULT, RESULT_BASE, 5) "\n"                                             \
    "movq %[" #T2 "], " PTR2(RESULT, RESULT_BASE, 6) "\n"                                             \
    "movq %[" #T3 "], " PTR2(RESULT, RESULT_BASE, 7) "\n"                              

namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops {
    template<class Field>
    inline void subtract_8_limbs_mod_x86(limb_array &result, const limb_array &other) {
        set_static_modulus_limbs_from_field();
        limb scratch, t0, t1, t2, t3, q0, q1, q2, q3;
        asm volatile(
            subtract_mod(result, 0, other, 0, scratch, t0, t1, t2, t3, q0, q1, q2, q3)
            : [scratch]"=&r"(scratch),
              [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3),
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
            : "cc", "memory"
        );
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
    inline void fp2_mul_pre_x86(limb_array *z, const limb_array *x, const limb_array *y) {
        set_static_modulus_limbs_from_field();

        // For x = a + bu and y = c + du:
        //   xy = (a + bu) * (c + du)
        //      = ac + adu + bcu + bdu^2
        //      = ac + (ad + bc)u - bd      # since u^2 = -1
        //      = (ac - bd) + (ad + bc)u
        
        const limb_array &a = x[0];
        const limb_array &b = x[1];
        const limb_array &c = y[0];
        const limb_array &d = y[1];
        limb_array &z0 = z[0];
        limb_array &z1 = z[1];

        limb low, high, zero, d0, d1, d2, d3, d4, tmp;
        limb_array scratch;
        asm volatile(
            schoolbook(x, 0, y, 0, z, 0)
            schoolbook(x, 8, y, 8, scratch, 0)
            subtract_mod(z, 0, scratch, 0, tmp, low, high, zero, d0, d1, d2, d3, d4)

            :   [low]"=&r"(low),
                [high]"=&r"(high),
                [zero]"=&r"(zero),
                [d0]"=&r"(d0),
                [d1]"=&r"(d1),
                [d2]"=&r"(d2),
                [d3]"=&r"(d3),
                [d4]"=&r"(d4),
                [tmp]"=&r"(tmp)
            :   
                [x]"r"(x),
                [y]"r"(y),
                [z]"r"(z),
                [scratch]"r"(scratch.data()),
                [p0]"m"(p0),
                [p1]"m"(p1),
                [p2]"m"(p2),
                [p3]"m"(p3)
            : "rdx", "cc", "memory"
        );

        // limb_array bd;
        // multiply_4x4_x86(z0, a, c);
        // multiply_4x4_x86(bd, b, d);
        // subtract_8_limbs_mod_x86<Field>(z0, bd);

        limb_array bc;
        multiply_4x4_x86(z1, a, d);
        multiply_4x4_x86(bc, b, c);
        add_8_limbs_mod_x86<Field>(z1, bc);

    }
}    // namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops

#undef STR_IMPL
#undef STR
#undef PTR
#undef PTR2
#undef initial_schoolbook_round
#undef schoolbook_round
#undef final_schoolbook_round
#undef schoolbook
#undef subtract_mod
#undef set_static_modulus_limbs_from_field
#undef D
#undef T
#undef montgomery_reduce_load_next
#undef montgomery_reduce_cancel_low
#undef double_mod_p