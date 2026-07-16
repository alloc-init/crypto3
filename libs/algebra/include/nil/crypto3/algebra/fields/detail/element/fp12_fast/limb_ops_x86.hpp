#pragma once

#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/types.hpp>
#include <boost/preprocessor.hpp>

// clang-format off

#define STR_IMPL(X) #X
#define STR(X) STR_IMPL(X)

#define PTR(REGNAME, I) STR(I) "*8(%[" #REGNAME "])"
#define PTR2(REGNAME, I, J) PTR(REGNAME, BOOST_PP_ADD(I, J))

// Get the (i + j) mod 4 d register in the 4-live-limb schoolbook window.
#define D4(I, J) "%[d" STR(BOOST_PP_MOD(BOOST_PP_ADD(I, J), 4)) "]"

#define SCHOOLBOOK_4x4_ROUND(ROUND, Z, Z_BASE, X, X_BASE, Y, Y_BASE) \
    "mov " PTR2(Y, Y_BASE, ROUND) ", %%rdx\n"                        \
    "mulx " PTR2(X, X_BASE, 0) ", %[low], %[high]\n"                 \
    "adox %[low], " D4(0, ROUND) "\n"                                \
    "mov " D4(0, ROUND) ", " PTR2(Z, Z_BASE, ROUND) "\n"             \
    "adcx %[high], " D4(1, ROUND) "\n"                               \
    "mulx " PTR2(X, X_BASE, 1) ", %[low], %[high]\n"                 \
    "adox %[low], " D4(1, ROUND) "\n"                                \
    "adcx %[high], " D4(2, ROUND) "\n"                               \
    "mulx " PTR2(X, X_BASE, 2) ", %[low], %[high]\n"                 \
    "adox %[low], " D4(2, ROUND) "\n"                                \
    "adcx %[high], " D4(3, ROUND) "\n"                               \
    "mulx " PTR2(X, X_BASE, 3) ", %[low], " D4(4, ROUND) "\n"        \
    "adox %[low], " D4(3, ROUND) "\n"                                \
    "adox %[zero], " D4(4, ROUND) "\n"                               \
    "adcx %[zero], " D4(4, ROUND) "\n"

#define SCHOOLBOOK_4x4(Z, Z_BASE, X, X_BASE, Y, Y_BASE)             \
    "xor %[zero], %[zero]\n"                                        \
    "mov " PTR2(Y, Y_BASE, 0) ", %%rdx\n"                           \
    "mulx " PTR2(X, X_BASE, 0) ", %[d0], %[d1]\n"                   \
    "mov %[d0], " PTR2(Z, Z_BASE, 0) "\n"                           \
    "mulx " PTR2(X, X_BASE, 1) ", %[low], %[d2]\n"                  \
    "add %[low], %[d1]\n"                                           \
    "mulx " PTR2(X, X_BASE, 2) ", %[low], %[d3]\n"                  \
    "adc %[low], %[d2]\n"                                           \
    "mulx " PTR2(X, X_BASE, 3) ", %[low], %[d0]\n"                  \
    "adc %[low], %[d3]\n"                                           \
    "adc $0, %[d0]\n"                                               \
    SCHOOLBOOK_4x4_ROUND(1, Z, Z_BASE, X, X_BASE, Y, Y_BASE)        \
    SCHOOLBOOK_4x4_ROUND(2, Z, Z_BASE, X, X_BASE, Y, Y_BASE)        \
    SCHOOLBOOK_4x4_ROUND(3, Z, Z_BASE, X, X_BASE, Y, Y_BASE)        \
    "mov " D4(0, 4) ", " PTR2(Z, Z_BASE, 4) "\n"                    \
    "mov " D4(1, 4) ", " PTR2(Z, Z_BASE, 5) "\n"                    \
    "mov " D4(2, 4) ", " PTR2(Z, Z_BASE, 6) "\n"                    \
    "mov " D4(3, 4) ", " PTR2(Z, Z_BASE, 7) "\n"

// Get the (i + j) mod 6 d register in the 6-live-limb schoolbook window.
#define D6(I, J) "%[d" STR(BOOST_PP_MOD(BOOST_PP_ADD(I, J), 6)) "]"

#define SCHOOLBOOK_6x6_ROUND(ROUND, Z, Z_BASE, X, X_BASE, Y, Y_BASE) \
    "mov " PTR2(Y, Y_BASE, ROUND) ", %%rdx\n"                        \
    "mulx " PTR2(X, X_BASE, 0) ", %[low], %[high]\n"                 \
    "adox %[low], " D6(0, ROUND) "\n"                                \
    "mov " D6(0, ROUND) ", " PTR2(Z, Z_BASE, ROUND) "\n"             \
    "adcx %[high], " D6(1, ROUND) "\n"                               \
    "mulx " PTR2(X, X_BASE, 1) ", %[low], %[high]\n"                 \
    "adox %[low], " D6(1, ROUND) "\n"                                \
    "adcx %[high], " D6(2, ROUND) "\n"                               \
    "mulx " PTR2(X, X_BASE, 2) ", %[low], %[high]\n"                 \
    "adox %[low], " D6(2, ROUND) "\n"                                \
    "adcx %[high], " D6(3, ROUND) "\n"                               \
    "mulx " PTR2(X, X_BASE, 3) ", %[low], %[high]\n"                 \
    "adox %[low], " D6(3, ROUND) "\n"                                \
    "adcx %[high], " D6(4, ROUND) "\n"                               \
    "mulx " PTR2(X, X_BASE, 4) ", %[low], %[high]\n"                 \
    "adox %[low], " D6(4, ROUND) "\n"                                \
    "adcx %[high], " D6(5, ROUND) "\n"                               \
    "mulx " PTR2(X, X_BASE, 5) ", %[low], " D6(6, ROUND) "\n"        \
    "adox %[low], " D6(5, ROUND) "\n"                                \
    "adox %[zero], " D6(6, ROUND) "\n"                               \
    "adcx %[zero], " D6(6, ROUND) "\n"

#define SCHOOLBOOK_6x6(Z, Z_BASE, X, X_BASE, Y, Y_BASE)             \
    "xor %[zero], %[zero]\n"                                        \
    "mov " PTR2(Y, Y_BASE, 0) ", %%rdx\n"                           \
    "mulx " PTR2(X, X_BASE, 0) ", %[d0], %[d1]\n"                   \
    "mov %[d0], " PTR2(Z, Z_BASE, 0) "\n"                           \
    "mulx " PTR2(X, X_BASE, 1) ", %[low], %[d2]\n"                  \
    "add %[low], %[d1]\n"                                           \
    "mulx " PTR2(X, X_BASE, 2) ", %[low], %[d3]\n"                  \
    "adc %[low], %[d2]\n"                                           \
    "mulx " PTR2(X, X_BASE, 3) ", %[low], %[d4]\n"                  \
    "adc %[low], %[d3]\n"                                           \
    "mulx " PTR2(X, X_BASE, 4) ", %[low], %[d5]\n"                  \
    "adc %[low], %[d4]\n"                                           \
    "mulx " PTR2(X, X_BASE, 5) ", %[low], %[d0]\n"                  \
    "adc %[low], %[d5]\n"                                           \
    "adc $0, %[d0]\n"                                               \
    SCHOOLBOOK_6x6_ROUND(1, Z, Z_BASE, X, X_BASE, Y, Y_BASE)        \
    SCHOOLBOOK_6x6_ROUND(2, Z, Z_BASE, X, X_BASE, Y, Y_BASE)        \
    SCHOOLBOOK_6x6_ROUND(3, Z, Z_BASE, X, X_BASE, Y, Y_BASE)        \
    SCHOOLBOOK_6x6_ROUND(4, Z, Z_BASE, X, X_BASE, Y, Y_BASE)        \
    SCHOOLBOOK_6x6_ROUND(5, Z, Z_BASE, X, X_BASE, Y, Y_BASE)        \
    "mov " D6(0, 6) ", " PTR2(Z, Z_BASE, 6) "\n"                    \
    "mov " D6(1, 6) ", " PTR2(Z, Z_BASE, 7) "\n"                    \
    "mov " D6(2, 6) ", " PTR2(Z, Z_BASE, 8) "\n"                    \
    "mov " D6(3, 6) ", " PTR2(Z, Z_BASE, 9) "\n"                    \
    "mov " D6(4, 6) ", " PTR2(Z, Z_BASE, 10) "\n"                   \
    "mov " D6(5, 6) ", " PTR2(Z, Z_BASE, 11) "\n"

#define DOUBLE_MOD_P()      \
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

#define MUL_8_LIMBS_BY_9(DST, DST_BASE, SRC, SRC_BASE) \
    "movq " PTR2(SRC, SRC_BASE, 0) ", %[t0]\n"         \
    "movq " PTR2(SRC, SRC_BASE, 1) ", %[t1]\n"         \
    "movq " PTR2(SRC, SRC_BASE, 2) ", %[t2]\n"         \
    "movq " PTR2(SRC, SRC_BASE, 3) ", %[t3]\n"         \
    "movq " PTR2(SRC, SRC_BASE, 4) ", %[t4]\n"         \
    "movq " PTR2(SRC, SRC_BASE, 5) ", %[t5]\n"         \
    "movq " PTR2(SRC, SRC_BASE, 6) ", %[t6]\n"         \
    "movq " PTR2(SRC, SRC_BASE, 7) ", %[t7]\n"         \
    DOUBLE_MOD_P() /* 2x */                            \
    DOUBLE_MOD_P() /* 4x */                            \
    DOUBLE_MOD_P() /* 8x */                            \
    "addq " PTR2(SRC, SRC_BASE, 0) ", %[t0]\n"         \
    "adcq " PTR2(SRC, SRC_BASE, 1) ", %[t1]\n"         \
    "adcq " PTR2(SRC, SRC_BASE, 2) ", %[t2]\n"         \
    "adcq " PTR2(SRC, SRC_BASE, 3) ", %[t3]\n"         \
    "adcq " PTR2(SRC, SRC_BASE, 4) ", %[t4]\n"         \
    "adcq " PTR2(SRC, SRC_BASE, 5) ", %[t5]\n"         \
    "adcq " PTR2(SRC, SRC_BASE, 6) ", %[t6]\n"         \
    "adcq " PTR2(SRC, SRC_BASE, 7) ", %[t7]\n"         \
    "movq %[t4], %[q0]\n"                              \
    "movq %[t5], %[q1]\n"                              \
    "movq %[t6], %[q2]\n"                              \
    "movq %[t7], %[q3]\n"                              \
    "subq %[p0], %[q0]\n"                              \
    "sbbq %[p1], %[q1]\n"                              \
    "sbbq %[p2], %[q2]\n"                              \
    "sbbq %[p3], %[q3]\n"                              \
    "cmovnc %[q0], %[t4]\n"                            \
    "cmovnc %[q1], %[t5]\n"                            \
    "cmovnc %[q2], %[t6]\n"                            \
    "cmovnc %[q3], %[t7]\n"                            \
    "movq %[t0], " PTR2(DST, DST_BASE, 0) "\n"         \
    "movq %[t1], " PTR2(DST, DST_BASE, 1) "\n"         \
    "movq %[t2], " PTR2(DST, DST_BASE, 2) "\n"         \
    "movq %[t3], " PTR2(DST, DST_BASE, 3) "\n"         \
    "movq %[t4], " PTR2(DST, DST_BASE, 4) "\n"         \
    "movq %[t5], " PTR2(DST, DST_BASE, 5) "\n"         \
    "movq %[t6], " PTR2(DST, DST_BASE, 6) "\n"         \
    "movq %[t7], " PTR2(DST, DST_BASE, 7) "\n"

#define MUL_12_LIMBS_BY_5_LOOP(DST, I) \
    "movq " PTR(DST, I) ", %[d0]\n" \
    "movq $2, %%rcx\n" \
    "shlxq %%rcx, %[d0], %[d2]\n" \
    "leaq (%[d1], %[d2]), %[d2]\n" \
    "movq $62, %%rcx\n" \
    "shrxq %%rcx, %[d0], %[d1]\n" \
    "adcx %[d2], %[d0]\n" \
    "movq %[d0], " PTR(DST, I) "\n"

#define MUL_12_LIMBS_BY_5(DST) \
    "movq " PTR(DST, 0) ", %[d0]\n" \
    "movq %[d0], %[d1]\n" \
    "shlq $2, %[d0]\n" \
    "shrq $62, %[d1]\n" \
    "addq %[d0], " PTR(DST, 0) "\n" \
    MUL_12_LIMBS_BY_5_LOOP(DST, 1) \
    MUL_12_LIMBS_BY_5_LOOP(DST, 2) \
    MUL_12_LIMBS_BY_5_LOOP(DST, 3) \
    MUL_12_LIMBS_BY_5_LOOP(DST, 4) \
    MUL_12_LIMBS_BY_5_LOOP(DST, 5) \
    MUL_12_LIMBS_BY_5_LOOP(DST, 6) \
    MUL_12_LIMBS_BY_5_LOOP(DST, 7) \
    MUL_12_LIMBS_BY_5_LOOP(DST, 8) \
    MUL_12_LIMBS_BY_5_LOOP(DST, 9) \
    MUL_12_LIMBS_BY_5_LOOP(DST, 10) \
    MUL_12_LIMBS_BY_5_LOOP(DST, 11)

#define GET_MODULUS_4_LIMBS(FIELD)                                  \
    constexpr auto mod_obj = FIELD::modulus_params.get_mod_obj();   \
    static constexpr limb p0 = limb(mod_obj.get_mod().limbs()[0]);  \
    static constexpr limb p1 = limb(mod_obj.get_mod().limbs()[1]);  \
    static constexpr limb p2 = limb(mod_obj.get_mod().limbs()[2]);  \
    static constexpr limb p3 = limb(mod_obj.get_mod().limbs()[3]);  \
    static constexpr limb p_dash = limb(mod_obj.get_p_dash());

#define GET_MODULUS_6_LIMBS(FIELD)                                  \
    constexpr auto mod_obj = FIELD::modulus_params.get_mod_obj();   \
    static constexpr limb p0 = limb(mod_obj.get_mod().limbs()[0]);  \
    static constexpr limb p1 = limb(mod_obj.get_mod().limbs()[1]);  \
    static constexpr limb p2 = limb(mod_obj.get_mod().limbs()[2]);  \
    static constexpr limb p3 = limb(mod_obj.get_mod().limbs()[3]);  \
    static constexpr limb p4 = limb(mod_obj.get_mod().limbs()[4]);  \
    static constexpr limb p5 = limb(mod_obj.get_mod().limbs()[5]);  \
    static constexpr limb p_dash = limb(mod_obj.get_p_dash());

// Get the (i + j) mod 5 t register in the REDC rotating window for 8 limb values.
#define T8(I, J) "%[t" STR(BOOST_PP_MOD(BOOST_PP_ADD(I, J), 5)) "]"

#define REDC_8_LIMBS_LOAD_NEXT(I)                \
    "setc %b[pending]\n"                         \
    "movq " PTR2(data, I, 5) ", " T8(I, 5) "\n"

// main body of loop in montgomery reduce
#define REDC_8_LIMBS_CANCEL_LOW(I)                                         \
    /* m = data[i] * p_dash */                                             \
    "movq " T8(I, 0) ", %%rdx\n"                                           \
    "imulq %[p_dash], %%rdx\n"                                             \
    "xor %[zero], %[zero]\n"                                               \
    /* multiply m by each modulus limb and add into the rotating window */ \
    "mulxq %[p0], %[low], %[high]\n"                                       \
    /* use dual carry chains */                                            \
    "adcx %[low], " T8(I, 0) "\n"                                          \
    "adox %[high], " T8(I, 1) "\n"                                         \
    "mulxq %[p1], %[low], %[high]\n"                                       \
    "adcx %[low], " T8(I, 1) "\n"                                          \
    "adox %[high], " T8(I, 2) "\n"                                         \
    "mulxq %[p2], %[low], %[high]\n"                                       \
    "adcx %[low], " T8(I, 2) "\n"                                          \
    "adox %[high], " T8(I, 3) "\n"                                         \
    "mulxq %[p3], %[low], %[high]\n"                                       \
    "adcx %[low], " T8(I, 3) "\n"                                          \
    /* merge OF, prior pending carry, and CF into the next live limb */    \
    "adox %[zero], %[high]\n"                                              \
    "adox %[pending], %[high]\n"                                           \
    "adcx %[high], " T8(I, 4) "\n"                                         \
    BOOST_PP_IF(BOOST_PP_LESS(I, 3), REDC_8_LIMBS_LOAD_NEXT(I), "")

// Get the (i + j) mod 7 t register in the REDC rotating window for 12 limb values.
#define T12(I, J) "%[t" STR(BOOST_PP_MOD(BOOST_PP_ADD(I, J), 7)) "]"

#define REDC_12_LIMBS_LOAD_NEXT(I)               \
    "setc %b[pending]\n"                         \
    "movq " PTR2(data, I, 7) ", " T12(I, 7) "\n"

// main body of loop in montgomery reduce
#define REDC_12_LIMBS_CANCEL_LOW(I)                                         \
    /* m = data[i] * p_dash */                                              \
    "movq " T12(I, 0) ", %%rdx\n"                                           \
    "imulq %[p_dash], %%rdx\n"                                              \
    "xor %[zero], %[zero]\n"                                                \
    /* multiply m by each modulus limb and add into the rotating window */  \
    "mulxq %[p0], %[low], %[high]\n"                                        \
    /* use dual carry chains */                                             \
    "adcx %[low], " T12(I, 0) "\n"                                          \
    "adox %[high], " T12(I, 1) "\n"                                         \
    "mulxq %[p1], %[low], %[high]\n"                                        \
    "adcx %[low], " T12(I, 1) "\n"                                          \
    "adox %[high], " T12(I, 2) "\n"                                         \
    "mulxq %[p2], %[low], %[high]\n"                                        \
    "adcx %[low], " T12(I, 2) "\n"                                          \
    "adox %[high], " T12(I, 3) "\n"                                         \
    "mulxq %[p3], %[low], %[high]\n"                                        \
    "adcx %[low], " T12(I, 3) "\n"                                          \
    "adox %[high], " T12(I, 4) "\n"                                         \
    "mulxq %[p4], %[low], %[high]\n"                                        \
    "adcx %[low], " T12(I, 4) "\n"                                          \
    "adox %[high], " T12(I, 5) "\n"                                         \
    "mulxq %[p5], %[low], %[high]\n"                                        \
    "adcx %[low], " T12(I, 5) "\n"                                          \
    /* merge OF, prior pending carry, and CF into the next live limb */     \
    "adox %[zero], %[high]\n"                                               \
    "adox %[pending], %[high]\n"                                            \
    "adcx %[high], " T12(I, 6) "\n"                                         \
    BOOST_PP_IF(BOOST_PP_LESS(I, 5), REDC_12_LIMBS_LOAD_NEXT(I), "")

#define ADD_4_LIMBS(Z, Z_BASE, X, X_BASE, Y, Y_BASE, TMP) \
    "movq " PTR2(X, X_BASE, 0) ", %[" #TMP "]\n"        \
    "addq " PTR2(Y, Y_BASE, 0) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 0) "\n"        \
    "movq " PTR2(X, X_BASE, 1) ", %[" #TMP "]\n"        \
    "adcq " PTR2(Y, Y_BASE, 1) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 1) "\n"        \
    "movq " PTR2(X, X_BASE, 2) ", %[" #TMP "]\n"        \
    "adcq " PTR2(Y, Y_BASE, 2) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 2) "\n"        \
    "movq " PTR2(X, X_BASE, 3) ", %[" #TMP "]\n"        \
    "adcq " PTR2(Y, Y_BASE, 3) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 3) "\n"

#define ADD_6_LIMBS(Z, Z_BASE, X, X_BASE, Y, Y_BASE, TMP) \
    "movq " PTR2(X, X_BASE, 0) ", %[" #TMP "]\n"        \
    "addq " PTR2(Y, Y_BASE, 0) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 0) "\n"        \
    "movq " PTR2(X, X_BASE, 1) ", %[" #TMP "]\n"        \
    "adcq " PTR2(Y, Y_BASE, 1) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 1) "\n"        \
    "movq " PTR2(X, X_BASE, 2) ", %[" #TMP "]\n"        \
    "adcq " PTR2(Y, Y_BASE, 2) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 2) "\n"        \
    "movq " PTR2(X, X_BASE, 3) ", %[" #TMP "]\n"        \
    "adcq " PTR2(Y, Y_BASE, 3) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 3) "\n"        \
    "movq " PTR2(X, X_BASE, 4) ", %[" #TMP "]\n"        \
    "adcq " PTR2(Y, Y_BASE, 4) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 4) "\n"        \
    "movq " PTR2(X, X_BASE, 5) ", %[" #TMP "]\n"        \
    "adcq " PTR2(Y, Y_BASE, 5) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 5) "\n"

#define ADD_8_LIMBS(Z, Z_BASE, X, X_BASE, Y, Y_BASE, TMP)   \
    "movq " PTR2(X, X_BASE, 0) ", %[" #TMP "]\n"            \
    "addq " PTR2(Y, Y_BASE, 0) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 0) "\n"            \
    "movq " PTR2(X, X_BASE, 1) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 1) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 1) "\n"            \
    "movq " PTR2(X, X_BASE, 2) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 2) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 2) "\n"            \
    "movq " PTR2(X, X_BASE, 3) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 3) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 3) "\n"            \
    "movq " PTR2(X, X_BASE, 4) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 4) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 4) "\n"            \
    "movq " PTR2(X, X_BASE, 5) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 5) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 5) "\n"            \
    "movq " PTR2(X, X_BASE, 6) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 6) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 6) "\n"            \
    "movq " PTR2(X, X_BASE, 7) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 7) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 7) "\n"

#define ADD_12_LIMBS(Z, Z_BASE, X, X_BASE, Y, Y_BASE, TMP)  \
    "movq " PTR2(X, X_BASE, 0) ", %[" #TMP "]\n"            \
    "addq " PTR2(Y, Y_BASE, 0) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 0) "\n"            \
    "movq " PTR2(X, X_BASE, 1) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 1) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 1) "\n"            \
    "movq " PTR2(X, X_BASE, 2) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 2) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 2) "\n"            \
    "movq " PTR2(X, X_BASE, 3) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 3) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 3) "\n"            \
    "movq " PTR2(X, X_BASE, 4) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 4) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 4) "\n"            \
    "movq " PTR2(X, X_BASE, 5) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 5) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 5) "\n"            \
    "movq " PTR2(X, X_BASE, 6) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 6) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 6) "\n"            \
    "movq " PTR2(X, X_BASE, 7) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 7) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 7) "\n"            \
    "movq " PTR2(X, X_BASE, 8) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 8) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 8) "\n"            \
    "movq " PTR2(X, X_BASE, 9) ", %[" #TMP "]\n"            \
    "adcq " PTR2(Y, Y_BASE, 9) ", %[" #TMP "]\n"            \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 9) "\n"            \
    "movq " PTR2(X, X_BASE, 10) ", %[" #TMP "]\n"           \
    "adcq " PTR2(Y, Y_BASE, 10) ", %[" #TMP "]\n"           \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 10) "\n"           \
    "movq " PTR2(X, X_BASE, 11) ", %[" #TMP "]\n"           \
    "adcq " PTR2(Y, Y_BASE, 11) ", %[" #TMP "]\n"           \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 11) "\n"

#define ADD_4_LIMBS_MOD(Z, Z_BASE, X, X_BASE, Y, Y_BASE)    \
    "movq " PTR2(X, X_BASE, 0) ", %[t0]\n"                  \
    "movq " PTR2(X, X_BASE, 1) ", %[t1]\n"                  \
    "movq " PTR2(X, X_BASE, 2) ", %[t2]\n"                  \
    "movq " PTR2(X, X_BASE, 3) ", %[t3]\n"                  \
    "movq " PTR2(Y, Y_BASE, 0) ", %[q0]\n"                  \
    "addq %[q0], %[t0]\n"                                   \
    "movq " PTR2(Y, Y_BASE, 1) ", %[q0]\n"                  \
    "adcq %[q0], %[t1]\n"                                   \
    "movq " PTR2(Y, Y_BASE, 2) ", %[q0]\n"                  \
    "adcq %[q0], %[t2]\n"                                   \
    "movq " PTR2(Y, Y_BASE, 3) ", %[q0]\n"                  \
    "adcq %[q0], %[t3]\n"                                   \
    "movq %[t0], %[q0]\n"                                   \
    "movq %[t1], %[q1]\n"                                   \
    "movq %[t2], %[q2]\n"                                   \
    "movq %[t3], %[q3]\n"                                   \
    "subq %[p0], %[q0]\n"                                   \
    "sbbq %[p1], %[q1]\n"                                   \
    "sbbq %[p2], %[q2]\n"                                   \
    "sbbq %[p3], %[q3]\n"                                   \
    "cmovnc %[q0], %[t0]\n"                                 \
    "cmovnc %[q1], %[t1]\n"                                 \
    "cmovnc %[q2], %[t2]\n"                                 \
    "cmovnc %[q3], %[t3]\n"                                 \
    "movq %[t0], " PTR2(Z, Z_BASE, 0) "\n"                  \
    "movq %[t1], " PTR2(Z, Z_BASE, 1) "\n"                  \
    "movq %[t2], " PTR2(Z, Z_BASE, 2) "\n"                  \
    "movq %[t3], " PTR2(Z, Z_BASE, 3) "\n"

#define ADD_6_LIMBS_MOD(Z, Z_BASE, X, X_BASE, Y, Y_BASE) \
    "movq " PTR2(X, X_BASE, 0) ", %%rdx\n"               \
    "addq " PTR2(Y, Y_BASE, 0) ", %%rdx\n"               \
    "movq %%rdx, " PTR2(Z, Z_BASE, 0) "\n"               \
    "movq " PTR2(X, X_BASE, 1) ", %%rdx\n"               \
    "adcq " PTR2(Y, Y_BASE, 1) ", %%rdx\n"               \
    "movq %%rdx, " PTR2(Z, Z_BASE, 1) "\n"               \
    "movq " PTR2(X, X_BASE, 2) ", %%rdx\n"               \
    "adcq " PTR2(Y, Y_BASE, 2) ", %%rdx\n"               \
    "movq %%rdx, " PTR2(Z, Z_BASE, 2) "\n"               \
    "movq " PTR2(X, X_BASE, 3) ", %%rdx\n"               \
    "adcq " PTR2(Y, Y_BASE, 3) ", %%rdx\n"               \
    "movq %%rdx, " PTR2(Z, Z_BASE, 3) "\n"               \
    "movq " PTR2(X, X_BASE, 4) ", %%rdx\n"               \
    "adcq " PTR2(Y, Y_BASE, 4) ", %%rdx\n"               \
    "movq %%rdx, " PTR2(Z, Z_BASE, 4) "\n"               \
    "movq " PTR2(X, X_BASE, 5) ", %%rdx\n"               \
    "adcq " PTR2(Y, Y_BASE, 5) ", %%rdx\n"               \
    "movq %%rdx, " PTR2(Z, Z_BASE, 5) "\n"               \
    "movq " PTR2(Z, Z_BASE, 0) ", %[t0]\n"               \
    "movq " PTR2(Z, Z_BASE, 1) ", %[t1]\n"               \
    "movq " PTR2(Z, Z_BASE, 2) ", %[t2]\n"               \
    "movq " PTR2(Z, Z_BASE, 3) ", %[t3]\n"               \
    "movq " PTR2(Z, Z_BASE, 4) ", %[t4]\n"               \
    "movq " PTR2(Z, Z_BASE, 5) ", %[t5]\n"               \
    "subq %[p0], %[t0]\n"                                \
    "sbbq %[p1], %[t1]\n"                                \
    "sbbq %[p2], %[t2]\n"                                \
    "sbbq %[p3], %[t3]\n"                                \
    "sbbq %[p4], %[t4]\n"                                \
    "sbbq %[p5], %[t5]\n"                                \
    "jc done" #Z_BASE "%=\n"                             \
    "movq %[t0], " PTR2(Z, Z_BASE, 0) "\n"               \
    "movq %[t1], " PTR2(Z, Z_BASE, 1) "\n"               \
    "movq %[t2], " PTR2(Z, Z_BASE, 2) "\n"               \
    "movq %[t3], " PTR2(Z, Z_BASE, 3) "\n"               \
    "movq %[t4], " PTR2(Z, Z_BASE, 4) "\n"               \
    "movq %[t5], " PTR2(Z, Z_BASE, 5) "\n"               \
    "done" #Z_BASE "%=:\n"

#define ADD_8_LIMBS_MOD(Z, Z_BASE, X, X_BASE, Y, Y_BASE, T0, T1, T2, T3, Q0, Q1, Q2, Q3) \
    "movq " PTR2(X, X_BASE, 4) ", %[" #T0 "]\n"                                          \
    "movq " PTR2(X, X_BASE, 5) ", %[" #T1 "]\n"                                          \
    "movq " PTR2(X, X_BASE, 6) ", %[" #T2 "]\n"                                          \
    "movq " PTR2(X, X_BASE, 7) ", %[" #T3 "]\n"                                          \
    "movq " PTR2(X, X_BASE, 0) ", %[" #Q0 "]\n"                                          \
    "addq " PTR2(Y, Y_BASE, 0) ", %[" #Q0 "]\n"                                          \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 0) "\n"                                          \
    "movq " PTR2(X, X_BASE, 1) ", %[" #Q0 "]\n"                                          \
    "adcq " PTR2(Y, Y_BASE, 1) ", %[" #Q0 "]\n"                                          \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 1) "\n"                                          \
    "movq " PTR2(X, X_BASE, 2) ", %[" #Q0 "]\n"                                          \
    "adcq " PTR2(Y, Y_BASE, 2) ", %[" #Q0 "]\n"                                          \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 2) "\n"                                          \
    "movq " PTR2(X, X_BASE, 3) ", %[" #Q0 "]\n"                                          \
    "adcq " PTR2(Y, Y_BASE, 3) ", %[" #Q0 "]\n"                                          \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 3) "\n"                                          \
    "movq " PTR2(Y, Y_BASE, 4) ", %[" #Q0 "]\n"                                          \
    "adcq %[" #Q0 "], %[" #T0 "]\n"                                                      \
    "movq " PTR2(Y, Y_BASE, 5) ", %[" #Q0 "]\n"                                          \
    "adcq %[" #Q0 "], %[" #T1 "]\n"                                                      \
    "movq " PTR2(Y, Y_BASE, 6) ", %[" #Q0 "]\n"                                          \
    "adcq %[" #Q0 "], %[" #T2 "]\n"                                                      \
    "movq " PTR2(Y, Y_BASE, 7) ", %[" #Q0 "]\n"                                          \
    "adcq %[" #Q0 "], %[" #T3 "]\n"                                                      \
    "movq %[" #T0 "], %[" #Q0 "]\n"                                                      \
    "movq %[" #T1 "], %[" #Q1 "]\n"                                                      \
    "movq %[" #T2 "], %[" #Q2 "]\n"                                                      \
    "movq %[" #T3 "], %[" #Q3 "]\n"                                                      \
    "subq %[p0], %[" #Q0 "]\n"                                                           \
    "sbbq %[p1], %[" #Q1 "]\n"                                                           \
    "sbbq %[p2], %[" #Q2 "]\n"                                                           \
    "sbbq %[p3], %[" #Q3 "]\n"                                                           \
    "cmovnc %[" #Q0 "], %[" #T0 "]\n"                                                    \
    "cmovnc %[" #Q1 "], %[" #T1 "]\n"                                                    \
    "cmovnc %[" #Q2 "], %[" #T2 "]\n"                                                    \
    "cmovnc %[" #Q3 "], %[" #T3 "]\n"                                                    \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 4) "\n"                                          \
    "movq %[" #T1 "], " PTR2(Z, Z_BASE, 5) "\n"                                          \
    "movq %[" #T2 "], " PTR2(Z, Z_BASE, 6) "\n"                                          \
    "movq %[" #T3 "], " PTR2(Z, Z_BASE, 7) "\n"

#define ADD_12_LIMBS_MOD(Z, Z_BASE, X, X_BASE, Y, Y_BASE, T0, T1, T2, T3, T4, T5, Q0, Q1, Q2, Q3, Q4, Q5) \
    "movq " PTR2(X, X_BASE, 6) ", %[" #T0 "]\n"     \
    "movq " PTR2(X, X_BASE, 7) ", %[" #T1 "]\n"     \
    "movq " PTR2(X, X_BASE, 8) ", %[" #T2 "]\n"     \
    "movq " PTR2(X, X_BASE, 9) ", %[" #T3 "]\n"     \
    "movq " PTR2(X, X_BASE, 10) ", %[" #T4 "]\n"    \
    "movq " PTR2(X, X_BASE, 11) ", %[" #T5 "]\n"    \
    "movq " PTR2(X, X_BASE, 0) ", %[" #Q0 "]\n"     \
    "addq " PTR2(Y, Y_BASE, 0) ", %[" #Q0 "]\n"     \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 0) "\n"     \
    "movq " PTR2(X, X_BASE, 1) ", %[" #Q0 "]\n"     \
    "adcq " PTR2(Y, Y_BASE, 1) ", %[" #Q0 "]\n"     \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 1) "\n"     \
    "movq " PTR2(X, X_BASE, 2) ", %[" #Q0 "]\n"     \
    "adcq " PTR2(Y, Y_BASE, 2) ", %[" #Q0 "]\n"     \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 2) "\n"     \
    "movq " PTR2(X, X_BASE, 3) ", %[" #Q0 "]\n"     \
    "adcq " PTR2(Y, Y_BASE, 3) ", %[" #Q0 "]\n"     \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 3) "\n"     \
    "movq " PTR2(X, X_BASE, 4) ", %[" #Q0 "]\n"     \
    "adcq " PTR2(Y, Y_BASE, 4) ", %[" #Q0 "]\n"     \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 4) "\n"     \
    "movq " PTR2(X, X_BASE, 5) ", %[" #Q0 "]\n"     \
    "adcq " PTR2(Y, Y_BASE, 5) ", %[" #Q0 "]\n"     \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 5) "\n"     \
    "movq " PTR2(Y, Y_BASE, 6) ", %[" #Q0 "]\n"     \
    "adcq %[" #Q0 "], %[" #T0 "]\n"                 \
    "movq " PTR2(Y, Y_BASE, 7) ", %[" #Q0 "]\n"     \
    "adcq %[" #Q0 "], %[" #T1 "]\n"                 \
    "movq " PTR2(Y, Y_BASE, 8) ", %[" #Q0 "]\n"     \
    "adcq %[" #Q0 "], %[" #T2 "]\n"                 \
    "movq " PTR2(Y, Y_BASE, 9) ", %[" #Q0 "]\n"     \
    "adcq %[" #Q0 "], %[" #T3 "]\n"                 \
    "movq " PTR2(Y, Y_BASE, 10) ", %[" #Q0 "]\n"    \
    "adcq %[" #Q0 "], %[" #T4 "]\n"                 \
    "movq " PTR2(Y, Y_BASE, 11) ", %[" #Q0 "]\n"    \
    "adcq %[" #Q0 "], %[" #T5 "]\n"                 \
    "movq %[" #T0 "], %[" #Q0 "]\n"                 \
    "movq %[" #T1 "], %[" #Q1 "]\n"                 \
    "movq %[" #T2 "], %[" #Q2 "]\n"                 \
    "movq %[" #T3 "], %[" #Q3 "]\n"                 \
    "movq %[" #T4 "], %[" #Q4 "]\n"                 \
    "movq %[" #T5 "], %[" #Q5 "]\n"                 \
    "subq %[p0], %[" #Q0 "]\n"                      \
    "sbbq %[p1], %[" #Q1 "]\n"                      \
    "sbbq %[p2], %[" #Q2 "]\n"                      \
    "sbbq %[p3], %[" #Q3 "]\n"                      \
    "sbbq %[p4], %[" #Q4 "]\n"                      \
    "sbbq %[p5], %[" #Q5 "]\n"                      \
    "cmovnc %[" #Q0 "], %[" #T0 "]\n"               \
    "cmovnc %[" #Q1 "], %[" #T1 "]\n"               \
    "cmovnc %[" #Q2 "], %[" #T2 "]\n"               \
    "cmovnc %[" #Q3 "], %[" #T3 "]\n"               \
    "cmovnc %[" #Q4 "], %[" #T4 "]\n"               \
    "cmovnc %[" #Q5 "], %[" #T5 "]\n"               \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 6) "\n"     \
    "movq %[" #T1 "], " PTR2(Z, Z_BASE, 7) "\n"     \
    "movq %[" #T2 "], " PTR2(Z, Z_BASE, 8) "\n"     \
    "movq %[" #T3 "], " PTR2(Z, Z_BASE, 9) "\n"     \
    "movq %[" #T4 "], " PTR2(Z, Z_BASE, 10) "\n"    \
    "movq %[" #T5 "], " PTR2(Z, Z_BASE, 11) "\n"

#define SUB_8_LIMBS(RESULT, RESULT_BASE, OTHER, OTHER_BASE, SCRATCH) \
    "movq " PTR2(RESULT, RESULT_BASE, 0) ", %[" #SCRATCH "]\n"     \
    "subq " PTR2(OTHER, OTHER_BASE, 0) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 0) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 1) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 1) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 1) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 2) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 2) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 2) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 3) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 3) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 3) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 4) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 4) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 4) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 5) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 5) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 5) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 6) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 6) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 6) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 7) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 7) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 7) "\n"

#define SUB_12_LIMBS(RESULT, RESULT_BASE, OTHER, OTHER_BASE, SCRATCH) \
    "movq " PTR2(RESULT, RESULT_BASE, 0) ", %[" #SCRATCH "]\n"     \
    "subq " PTR2(OTHER, OTHER_BASE, 0) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 0) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 1) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 1) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 1) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 2) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 2) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 2) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 3) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 3) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 3) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 4) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 4) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 4) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 5) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 5) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 5) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 6) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 6) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 6) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 7) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 7) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 7) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 8) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 8) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 8) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 9) ", %[" #SCRATCH "]\n"     \
    "sbbq " PTR2(OTHER, OTHER_BASE, 9) ", %[" #SCRATCH "]\n"       \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 9) "\n"     \
    "movq " PTR2(RESULT, RESULT_BASE, 10) ", %[" #SCRATCH "]\n"    \
    "sbbq " PTR2(OTHER, OTHER_BASE, 10) ", %[" #SCRATCH "]\n"      \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 10) "\n"    \
    "movq " PTR2(RESULT, RESULT_BASE, 11) ", %[" #SCRATCH "]\n"    \
    "sbbq " PTR2(OTHER, OTHER_BASE, 11) ", %[" #SCRATCH "]\n"      \
    "movq %[" #SCRATCH "], " PTR2(RESULT, RESULT_BASE, 11) "\n"

#define SUB_8_LIMBS_MOD(Z, Z_BASE, X, X_BASE, Y, Y_BASE, T0, T1, T2, T3) \
    "movq " PTR2(X, X_BASE, 0) ", %[" #T0 "]\n"                        \
    "subq " PTR2(Y, Y_BASE, 0) ", %[" #T0 "]\n"                        \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 0)"\n"                         \
    "movq " PTR2(X, X_BASE, 1) ", %[" #T0 "]\n"                        \
    "sbbq " PTR2(Y, Y_BASE, 1) ", %[" #T0 "]\n"                        \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 1)"\n"                         \
    "movq " PTR2(X, X_BASE, 2) ", %[" #T0 "]\n"                        \
    "sbbq " PTR2(Y, Y_BASE, 2) ", %[" #T0 "]\n"                        \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 2)"\n"                         \
    "movq " PTR2(X, X_BASE, 3) ", %[" #T0 "]\n"                        \
    "sbbq " PTR2(Y, Y_BASE, 3) ", %[" #T0 "]\n"                        \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 3) "\n"                        \
    "movq " PTR2(X, X_BASE, 4) ", %[" #T0 "]\n"                        \
    "sbbq " PTR2(Y, Y_BASE, 4) ", %[" #T0 "]\n"                        \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 4)"\n"                         \
    "movq " PTR2(X, X_BASE, 5) ", %[" #T0 "]\n"                        \
    "sbbq " PTR2(Y, Y_BASE, 5) ", %[" #T0 "]\n"                        \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 5)"\n"                         \
    "movq " PTR2(X, X_BASE, 6) ", %[" #T0 "]\n"                        \
    "sbbq " PTR2(Y, Y_BASE, 6) ", %[" #T0 "]\n"                        \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 6) "\n"                        \
    "movq " PTR2(X, X_BASE, 7) ", %[" #T0 "]\n"                        \
    "sbbq " PTR2(Y, Y_BASE, 7) ", %[" #T0 "]\n"                        \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 7)"\n"                         \
    "mov $0, %[" #T0 "]\n"                                             \
    "mov $0, %[" #T1 "]\n"                                             \
    "mov $0, %[" #T2 "]\n"                                             \
    "mov $0, %[" #T3 "]\n"                                             \
    "cmovc %[p0], %[" #T0 "]\n"                                        \
    "cmovc %[p1], %[" #T1 "]\n"                                        \
    "cmovc %[p2], %[" #T2 "]\n"                                        \
    "cmovc %[p3], %[" #T3 "]\n"                                        \
    "addq %[" #T0 "], " PTR2(Z, Z_BASE, 4) "\n"                        \
    "adcq %[" #T1 "], " PTR2(Z, Z_BASE, 5) "\n"                        \
    "adcq %[" #T2 "], " PTR2(Z, Z_BASE, 6) "\n"                        \
    "adcq %[" #T3 "], " PTR2(Z, Z_BASE, 7) "\n"

#define SUB_12_LIMBS_MOD(Z, Z_BASE, X, X_BASE, Y, Y_BASE, T0, T1, T2, T3, T4, T5) \
    "movq " PTR2(X, X_BASE, 0) ", %[" #T0 "]\n"         \
    "subq " PTR2(Y, Y_BASE, 0) ", %[" #T0 "]\n"         \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 0)"\n"          \
    "movq " PTR2(X, X_BASE, 1) ", %[" #T0 "]\n"         \
    "sbbq " PTR2(Y, Y_BASE, 1) ", %[" #T0 "]\n"         \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 1)"\n"          \
    "movq " PTR2(X, X_BASE, 2) ", %[" #T0 "]\n"         \
    "sbbq " PTR2(Y, Y_BASE, 2) ", %[" #T0 "]\n"         \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 2)"\n"          \
    "movq " PTR2(X, X_BASE, 3) ", %[" #T0 "]\n"         \
    "sbbq " PTR2(Y, Y_BASE, 3) ", %[" #T0 "]\n"         \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 3) "\n"         \
    "movq " PTR2(X, X_BASE, 4) ", %[" #T0 "]\n"         \
    "sbbq " PTR2(Y, Y_BASE, 4) ", %[" #T0 "]\n"         \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 4)"\n"          \
    "movq " PTR2(X, X_BASE, 5) ", %[" #T0 "]\n"         \
    "sbbq " PTR2(Y, Y_BASE, 5) ", %[" #T0 "]\n"         \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 5)"\n"          \
    "movq " PTR2(X, X_BASE, 6) ", %[" #T0 "]\n"         \
    "sbbq " PTR2(Y, Y_BASE, 6) ", %[" #T0 "]\n"         \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 6) "\n"         \
    "movq " PTR2(X, X_BASE, 7) ", %[" #T0 "]\n"         \
    "sbbq " PTR2(Y, Y_BASE, 7) ", %[" #T0 "]\n"         \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 7)"\n"          \
    "movq " PTR2(X, X_BASE, 8) ", %[" #T0 "]\n"         \
    "sbbq " PTR2(Y, Y_BASE, 8) ", %[" #T0 "]\n"         \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 8)"\n"          \
    "movq " PTR2(X, X_BASE, 9) ", %[" #T0 "]\n"         \
    "sbbq " PTR2(Y, Y_BASE, 9) ", %[" #T0 "]\n"         \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 9)"\n"          \
    "movq " PTR2(X, X_BASE, 10) ", %[" #T0 "]\n"        \
    "sbbq " PTR2(Y, Y_BASE, 10) ", %[" #T0 "]\n"        \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 10)"\n"         \
    "movq " PTR2(X, X_BASE, 11) ", %[" #T0 "]\n"        \
    "sbbq " PTR2(Y, Y_BASE, 11) ", %[" #T0 "]\n"        \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 11)"\n"         \
    "mov $0, %[" #T0 "]\n"                              \
    "mov $0, %[" #T1 "]\n"                              \
    "mov $0, %[" #T2 "]\n"                              \
    "mov $0, %[" #T3 "]\n"                              \
    "mov $0, %[" #T4 "]\n"                              \
    "mov $0, %[" #T5 "]\n"                              \
    "cmovc %[p0], %[" #T0 "]\n"                         \
    "cmovc %[p1], %[" #T1 "]\n"                         \
    "cmovc %[p2], %[" #T2 "]\n"                         \
    "cmovc %[p3], %[" #T3 "]\n"                         \
    "cmovc %[p4], %[" #T4 "]\n"                         \
    "cmovc %[p5], %[" #T5 "]\n"                         \
    "addq %[" #T0 "], " PTR2(Z, Z_BASE, 6) "\n"         \
    "adcq %[" #T1 "], " PTR2(Z, Z_BASE, 7) "\n"         \
    "adcq %[" #T2 "], " PTR2(Z, Z_BASE, 8) "\n"         \
    "adcq %[" #T3 "], " PTR2(Z, Z_BASE, 9) "\n"         \
    "adcq %[" #T4 "], " PTR2(Z, Z_BASE, 10) "\n"        \
    "adcq %[" #T5 "], " PTR2(Z, Z_BASE, 11) "\n"

namespace nil::crypto3::algebra::fields::detail::fp12_fast {
    template<Fp12FastParams Params>
        requires(Params::storage_limb_count == 8)
    inline void montgomery_reduce_x86(limb *result, const limb *data) {
        GET_MODULUS_4_LIMBS(Params::base_field_type);
        limb t0, t1, t2, t3, t4;
        limb low, high, pending, zero;
        asm volatile(
            "movq " PTR(data, 0) ", %[t0]\n"
            "movq " PTR(data, 1) ", %[t1]\n"
            "movq " PTR(data, 2) ", %[t2]\n"
            "movq " PTR(data, 3) ", %[t3]\n"
            "movq " PTR(data, 4) ", %[t4]\n"

            // Unrolled REDC rounds. Each round cancels the current low limb and
            // rotates in one higher input limb.
            "xor %[pending], %[pending]\n"
            REDC_8_LIMBS_CANCEL_LOW(0)
            REDC_8_LIMBS_CANCEL_LOW(1)
            REDC_8_LIMBS_CANCEL_LOW(2)
            REDC_8_LIMBS_CANCEL_LOW(3)

            // Conditionally reduce the high limbs modulo p, reusing scratch registers.
            // q = t - p
            "movq " T8(0, 4) ", %[low]\n"
            "movq " T8(1, 4) ", %[high]\n"
            "movq " T8(2, 4) ", %[zero]\n"
            "movq " T8(3, 4) ", %%rdx\n"

            // Try q = t - p.
            "subq %[p0], %[low]\n"
            "sbbq %[p1], %[high]\n"
            "sbbq %[p2], %[zero]\n"
            "sbbq %[p3], %%rdx\n"

            // If t - p did not borrow, keep q as the canonical result.
            "cmovnc %[low], " T8(0, 4) "\n"
            "cmovnc %[high], " T8(1, 4) "\n"
            "cmovnc %[zero], " T8(2, 4) "\n"
            "cmovnc %%rdx, " T8(3, 4) "\n"

            "movq " T8(0, 4) ", " PTR(result, 0) "\n"
            "movq " T8(1, 4) ", " PTR(result, 1) "\n"
            "movq " T8(2, 4) ", " PTR(result, 2) "\n"
            "movq " T8(3, 4) ", " PTR(result, 3) "\n"

            : [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3),
              [t4]"=&r"(t4),
              [low]"=&r"(low),
              [high]"=&r"(high),
              [pending]"=&r"(pending),
              [zero]"=&r"(zero)
            : [data]"r"(data),
              [result]"r"(result),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3),
              [p_dash]"m"(p_dash)
            : "rdx", "cc", "memory"
        );
    }


    template<Fp12FastParams Params>
        requires(Params::storage_limb_count == 12)
    inline void montgomery_reduce_x86(limb *result, const limb *data) {
        GET_MODULUS_6_LIMBS(Params::base_field_type); // sets p, p_dash
        limb t0, t1, t2, t3, t4, t5, t6;
        limb low, high, pending, zero;
        asm volatile(
            "movq " PTR(data, 0) ", %[t0]\n"
            "movq " PTR(data, 1) ", %[t1]\n"
            "movq " PTR(data, 2) ", %[t2]\n"
            "movq " PTR(data, 3) ", %[t3]\n"
            "movq " PTR(data, 4) ", %[t4]\n"
            "movq " PTR(data, 5) ", %[t5]\n"
            "movq " PTR(data, 6) ", %[t6]\n"

            // Unrolled REDC rounds. Each round cancels the current low limb and
            // rotates in one higher input limb.
            "xor %[pending], %[pending]\n"
            REDC_12_LIMBS_CANCEL_LOW(0)
            REDC_12_LIMBS_CANCEL_LOW(1)
            REDC_12_LIMBS_CANCEL_LOW(2)
            REDC_12_LIMBS_CANCEL_LOW(3)
            REDC_12_LIMBS_CANCEL_LOW(4)
            REDC_12_LIMBS_CANCEL_LOW(5)

            // Conditionally reduce the high limbs modulo p, reusing scratch registers.
            // q = t - p
            "movq " T12(0, 6) ", %[low]\n"
            "movq " T12(1, 6) ", %[high]\n"
            "movq " T12(2, 6) ", %[zero]\n"
            "movq " T12(3, 6) ", %%rdx\n"
            "movq " T12(4, 6) ", %[pending]\n"
            "movq " T12(5, 6) ", %[data]\n"

            // Try q = t - p.
            "subq %[p0], %[low]\n"
            "sbbq %[p1], %[high]\n"
            "sbbq %[p2], %[zero]\n"
            "sbbq %[p3], %%rdx\n"
            "sbbq %[p4], %[pending]\n"
            "sbbq %[p5], %[data]\n"

            // If t - p did not borrow, keep q as the canonical result.
            "cmovnc %[low], " T12(0, 6) "\n"
            "cmovnc %[high], " T12(1, 6) "\n"
            "cmovnc %[zero], " T12(2, 6) "\n"
            "cmovnc %%rdx, " T12(3, 6) "\n"
            "cmovnc %[pending], " T12(4, 6) "\n"
            "cmovnc %[data], " T12(5, 6) "\n"

            "movq " T12(0, 6) ", " PTR(result, 0) "\n"
            "movq " T12(1, 6) ", " PTR(result, 1) "\n"
            "movq " T12(2, 6) ", " PTR(result, 2) "\n"
            "movq " T12(3, 6) ", " PTR(result, 3) "\n"
            "movq " T12(4, 6) ", " PTR(result, 4) "\n"
            "movq " T12(5, 6) ", " PTR(result, 5) "\n"

            : [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3),
              [t4]"=&r"(t4),
              [t5]"=&r"(t5),
              [t6]"=&r"(t6),
              [low]"=&r"(low),
              [high]"=&r"(high),
              [pending]"=&r"(pending),
              [zero]"=&r"(zero)
            : [data]"r"(data),
              [result]"r"(result),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3),
              [p4]"m"(p4),
              [p5]"m"(p5),
              [p_dash]"m"(p_dash)
            : "rdx", "cc", "memory"
        );
    }

    template<Fp12FastParams Params>
        requires(Params::storage_limb_count == 8)
    inline void add_limbs_mod_x86(limb *z, const limb *x, const limb *y) {
        GET_MODULUS_4_LIMBS(Params::base_field_type);
        limb t0, t1, t2, t3, q0, q1, q2, q3;
        asm volatile(
            ADD_8_LIMBS_MOD(z, 0, x, 0, y, 0, t0, t1, t2, t3, q0, q1, q2, q3)
            : [t0]"+r"(t0),
              [t1]"+r"(t1),
              [t2]"+r"(t2),
              [t3]"+r"(t3),
              [q0]"=&r"(q0),
              [q1]"=&r"(q1),
              [q2]"=&r"(q2),
              [q3]"=&r"(q3)
            : [z]"r"(z),
              [x]"r"(x),
              [y]"r"(y),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3)
            : "cc", "memory"
        );
    }

    template<Fp12FastParams Params>
        requires(Params::storage_limb_count == 12)
    inline void add_limbs_mod_x86(limb *z, const limb *x, const limb *y) {
        GET_MODULUS_6_LIMBS(Params::base_field_type);
        limb t0, t1, t2, t3, t4, t5, q0, q1, q2, q3, q4, q5;
        asm volatile(
            ADD_12_LIMBS_MOD(z, 0, x, 0, y, 0,
                t0, t1, t2, t3, t4, t5,
                q0, q1, q2, q3, q4, q5)
            : [t0]"+r"(t0),
              [t1]"+r"(t1),
              [t2]"+r"(t2),
              [t3]"+r"(t3),
              [t4]"+r"(t4),
              [t5]"+r"(t5),
              [q0]"=&r"(q0),
              [q1]"=&r"(q1),
              [q2]"=&r"(q2),
              [q3]"=&r"(q3),
              [q4]"=&r"(q4),
              [q5]"=&r"(q5)
            : [z]"r"(z),
              [x]"r"(x),
              [y]"r"(y),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3),
              [p4]"m"(p4),
              [p5]"m"(p5)
            : "cc", "memory"
        );
    }

    template<Fp12FastParams Params>
        requires(Params::storage_limb_count == 8)
    inline void subtract_limbs_mod_x86(limb *z, const limb *x, const limb *y) {
        GET_MODULUS_4_LIMBS(Params::base_field_type);
        limb t0, t1, t2, t3;
        asm volatile(
            SUB_8_LIMBS_MOD(z, 0, x, 0, y, 0, t0, t1, t2, t3)
            : [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3)
            : [z]"r"(z),
              [x]"r"(x),
              [y]"r"(y),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3)
            : "cc", "memory"
        );
    }

    template<Fp12FastParams Params>
        requires(Params::storage_limb_count == 12)
    inline void subtract_limbs_mod_x86(limb *z, const limb *x, const limb *y) {
        GET_MODULUS_6_LIMBS(Params::base_field_type);
        limb t0, t1, t2, t3, t4, t5;
        asm volatile(
            SUB_12_LIMBS_MOD(z, 0, x, 0, y, 0, t0, t1, t2, t3, t4, t5)
            : [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3),
              [t4]"=&r"(t4),
              [t5]"=&r"(t5)
            : [z]"r"(z),
              [x]"r"(x),
              [y]"r"(y),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3),
              [p4]"m"(p4),
              [p5]"m"(p5)
            : "cc", "memory"
        );
    }

    template<class Field>
    inline void mul_8_limbs_by_9_x86(limb *dst, const limb *src) {
        GET_MODULUS_4_LIMBS(Field);
        limb t0, t1, t2, t3, t4, t5, t6, t7;
        limb q0, q1, q2, q3;
        asm volatile(
            MUL_8_LIMBS_BY_9(dst, 0, src, 0)
            : [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3),
              [t4]"=&r"(t4),
              [t5]"=&r"(t5),
              [t6]"=&r"(t6),
              [t7]"=&r"(t7),
              [q0]"=&r"(q0),
              [q1]"=&r"(q1),
              [q2]"=&r"(q2),
              [q3]"=&r"(q3)
            : [src]"r"(src),
              [dst]"r"(dst),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3)
            : "cc", "memory"
        );
    }

    template<Fp12FastParams Params>
        requires(Params::base_value_limb_count == 4)
    inline void fp2_base_add_mod_x86(limb *z, const limb *x, const limb *y) {
        GET_MODULUS_4_LIMBS(Params::base_field_type);
        limb t0, t1, t2, t3;
        limb q0, q1, q2, q3;
        asm volatile(
            ADD_4_LIMBS_MOD(z, 0, x, 0, y, 0)
            ADD_4_LIMBS_MOD(z, 4, x, 4, y, 4)
            : [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3),
              [q0]"=&r"(q0),
              [q1]"=&r"(q1),
              [q2]"=&r"(q2),
              [q3]"=&r"(q3)
            : [x]"r"(x),
              [y]"r"(y),
              [z]"r"(z),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3)
            : "cc", "memory"
        );
    }


    template<Fp12FastParams Params>
        requires(Params::base_value_limb_count == 6)
    inline void fp2_base_add_mod_x86(limb *z, const limb *x, const limb *y) {
        GET_MODULUS_6_LIMBS(Params::base_field_type);
        limb t0, t1, t2, t3, t4, t5;
        asm volatile(
            ADD_6_LIMBS_MOD(z, 0, x, 0, y, 0)
            ADD_6_LIMBS_MOD(z, 6, x, 6, y, 6)
            : [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3),
              [t4]"=&r"(t4),
              [t5]"=&r"(t5)
            : [x]"r"(x),
              [y]"r"(y),
              [z]"r"(z),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3),
              [p4]"m"(p4),
              [p5]"m"(p5)
            : "rdx", "cc", "memory"
        );
    }


    template<Fp12FastParams Params>
        requires(Params::storage_limb_count == 8)
    inline void fp2_sub_pre_x86(limb *data, const limb *other) {
        GET_MODULUS_4_LIMBS(Params::base_field_type);
        limb t0, t1, t2, t3;
        asm volatile(
            SUB_8_LIMBS_MOD(data, 0, data, 0, other, 0, t0, t1, t2, t3)
            SUB_8_LIMBS(data, 8, other, 8, t0)
            : [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3)
            : [data]"r"(data),
              [other]"r"(other),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3)
            : "cc", "memory"
        );
    }

    template<Fp12FastParams Params>
        requires(Params::storage_limb_count == 12)
    inline void fp2_sub_pre_x86(limb *data, const limb *other) {
        GET_MODULUS_6_LIMBS(Params::base_field_type);
        limb t0, t1, t2, t3, t4, t5;
        asm volatile(
            SUB_12_LIMBS_MOD(data, 0, data, 0, other, 0, t0, t1, t2, t3, t4, t5)
            SUB_12_LIMBS(data, 12, other, 12, t0)
            : [t0]"=&r"(t0),
              [t1]"=&r"(t1),
              [t2]"=&r"(t2),
              [t3]"=&r"(t3),
              [t4]"=&r"(t4),
              [t5]"=&r"(t5)
            : [data]"r"(data),
              [other]"r"(other),
              [p0]"m"(p0),
              [p1]"m"(p1),
              [p2]"m"(p2),
              [p3]"m"(p3),
              [p4]"m"(p4),
              [p5]"m"(p5)
            : "cc", "memory"
        );
    }

    template<Fp12FastParams Params>
        requires(Params::u_squared == -1 && Params::storage_limb_count == 8)
    inline void fp2_mul_pre_x86(limb *z, const limb *x, const limb *y) {
        GET_MODULUS_4_LIMBS(Params::base_field_type);
        // For x = a + bu and y = c + du:
        //   xy = (a + bu) * (c + du)
        //      = ac + adu + bcu + bdu^2
        //      = ac + (ad + bc)u - bd      # since u^2 = -1
        //      = (ac - bd) + (ad + bc)u
        limb low, high, zero, d0, d1, d2, d3;
        limb scratch[8];
        asm volatile(
            SCHOOLBOOK_4x4(z, 0, x, 0, y, 0)
            SCHOOLBOOK_4x4(scratch, 0, x, 4, y, 4)
            SUB_8_LIMBS_MOD(z, 0, z, 0, scratch, 0, d0, d1, d2, d3)
            SCHOOLBOOK_4x4(z, 8, x, 0, y, 4)
            SCHOOLBOOK_4x4(scratch, 0, x, 4, y, 0)
            ADD_8_LIMBS(z, 8, z, 8, scratch, 0, low)
            :   [low]"=&r"(low),
                [high]"=&r"(high),
                [zero]"=&r"(zero),
                [d0]"=&r"(d0),
                [d1]"=&r"(d1),
                [d2]"=&r"(d2),
                [d3]"=&r"(d3)
            :   [x]"r"(x),
                [y]"r"(y),
                [z]"r"(z),
                [scratch]"r"(scratch),
                [p0]"m"(p0),
                [p1]"m"(p1),
                [p2]"m"(p2),
                [p3]"m"(p3)
            : "rdx", "cc", "memory"
        );
    }

    template<Fp12FastParams Params>
        requires(Params::u_squared == -1 && Params::storage_limb_count == 12)
    inline void fp2_mul_pre_x86(limb *z, const limb *x, const limb *y) {
        GET_MODULUS_6_LIMBS(Params::base_field_type);
        // For x = a + bu and y = c + du:
        //   xy = (a + bu) * (c + du)
        //      = ac + adu + bcu + bdu^2
        //      = ac + (ad + bc)u - bd      # since u^2 = -1
        //      = (ac - bd) + (ad + bc)u
        limb low, high, zero, d0, d1, d2, d3, d4, d5;
        limb scratch[12];
        asm volatile(
            SCHOOLBOOK_6x6(z, 0, x, 0, y, 0) // z[0] = ac
            SCHOOLBOOK_6x6(scratch, 0, x, 6, y, 6) // scratch = bd
            SUB_12_LIMBS_MOD(z, 0, z, 0, scratch, 0, d0, d1, d2, d3, d4, d5) // z[0] -= bd == ac - bd
            SCHOOLBOOK_6x6(z, 12, x, 0, y, 6) // z[1] = ad
            SCHOOLBOOK_6x6(scratch, 0, x, 6, y, 0) // scratch = bc
            ADD_12_LIMBS(z, 12, z, 12, scratch, 0, low) // z[1] += bc == ad + bc
            :   [low]"=&r"(low),
                [high]"=&r"(high),
                [zero]"=&r"(zero),
                [d0]"=&r"(d0),
                [d1]"=&r"(d1),
                [d2]"=&r"(d2),
                [d3]"=&r"(d3),
                [d4]"=&r"(d4),
                [d5]"=&r"(d5)
            :   [x]"r"(x),
                [y]"r"(y),
                [z]"r"(z),
                [scratch]"r"(scratch),
                [p0]"m"(p0),
                [p1]"m"(p1),
                [p2]"m"(p2),
                [p3]"m"(p3),
                [p4]"m"(p4),
                [p5]"m"(p5)
            : "rdx", "cc", "memory"
        );
    }

    template<Fp12FastParams Params>
        requires(Params::u_squared == -5 && Params::storage_limb_count == 12)
    inline void fp2_mul_pre_x86(limb *z, const limb *x, const limb *y) {
        GET_MODULUS_6_LIMBS(Params::base_field_type);
        // For x = a + bu and y = c + du:
        //   xy = (a + bu) * (c + du)
        //      = ac + adu + bcu + bdu^2
        //      = ac + (ad + bc)u - bd      # since u^2 = -1
        //      = (ac - bd) + (ad + bc)u
        limb low, high, zero, d0, d1, d2, d3, d4, d5;
        limb scratch[24];
        asm volatile(
            SCHOOLBOOK_6x6(z, 0, x, 0, y, 0) // z[0] = ac
            SCHOOLBOOK_6x6(scratch, 0, x, 6, y, 6) // scratch = bd
            MUL_12_LIMBS_BY_5(scratch)
            SUB_12_LIMBS_MOD(z, 0, z, 0, scratch, 0, d0, d1, d2, d3, d4, d5) // z[0] -= bd == ac - bd
            SCHOOLBOOK_6x6(z, 12, x, 0, y, 6) // z[1] = ad
            SCHOOLBOOK_6x6(scratch, 0, x, 6, y, 0) // scratch = bc
            ADD_12_LIMBS(z, 12, z, 12, scratch, 0, low) // z[1] += bc == ad + bc
            :   [low]"=&r"(low),
                [high]"=&r"(high),
                [zero]"=&r"(zero),
                [d0]"=&r"(d0),
                [d1]"=&r"(d1),
                [d2]"=&r"(d2),
                [d3]"=&r"(d3),
                [d4]"=&r"(d4),
                [d5]"=&r"(d5)
            :   [x]"r"(x),
                [y]"r"(y),
                [z]"r"(z),
                [scratch]"r"(scratch),
                [p0]"m"(p0),
                [p1]"m"(p1),
                [p2]"m"(p2),
                [p3]"m"(p3),
                [p4]"m"(p4),
                [p5]"m"(p5)
            : "rdx", "cc", "memory"
        );
    }

    template<Fp12FastParams Params>
        requires(Params::u_squared == -1 && Params::storage_limb_count == 8)
    inline void fp2_add_mul_pre_x86(limb *z, const limb  *a, const limb  *b, const limb  *c, const limb  *d) {
        GET_MODULUS_4_LIMBS(Params::base_field_type);
        limb low, high, zero, d0, d1, d2, d3;
        struct {
            limb x[8];
            limb y[8];
            limb tmp[8];
        } data;
        asm volatile(
            ADD_4_LIMBS(data, 0, a, 0, b, 0, low)
            ADD_4_LIMBS(data, 4, a, 4, b, 4, low)
            ADD_4_LIMBS(data, 8, c, 0, d, 0, low)
            ADD_4_LIMBS(data, 12, c, 4, d, 4, low)
            :   [low]"=&r"(low)
            :   [data]"r"(&data),
                [a]"r"(a),
                [b]"r"(b),
                [c]"r"(c),
                [d]"r"(d)
            : "cc", "memory"
        );
        asm volatile(
            SCHOOLBOOK_4x4(z, 0, data, 0, data, 8) // z[0] = x[0] * y[0]
            SCHOOLBOOK_4x4(data, 16, data, 4, data, 12) // tmp = x[1] * y[1]
            SUB_8_LIMBS_MOD(z, 0, z, 0, data, 16, d0, d1, d2, d3)
            SCHOOLBOOK_4x4(z, 8, data, 0, data, 12)
            SCHOOLBOOK_4x4(data, 16, data, 4, data, 8)
            ADD_8_LIMBS(z, 8, z, 8, data, 16, low)
            :   [low]"=&r"(low),
                [high]"=&r"(high),
                [zero]"=&r"(zero),
                [d0]"=&r"(d0),
                [d1]"=&r"(d1),
                [d2]"=&r"(d2),
                [d3]"=&r"(d3)
            :   [data]"r"(&data),
                [z]"r"(z),
                [p0]"m"(p0),
                [p1]"m"(p1),
                [p2]"m"(p2),
                [p3]"m"(p3)
            : "rdx", "cc", "memory"
        );
    }

    template<Fp12FastParams Params>
        requires(Params::u_squared == -1 && Params::storage_limb_count == 12)
    inline void fp2_add_mul_pre_x86(limb *z, const limb  *a, const limb  *b, const limb  *c, const limb  *d) {
        GET_MODULUS_6_LIMBS(Params::base_field_type);
        limb low, high, zero, d0, d1, d2, d3, d4, d5;
        struct {
            limb x[12];
            limb y[12];
            limb tmp[12];
        } data;
        asm volatile(
            ADD_6_LIMBS(data, 0, a, 0, b, 0, low)
            ADD_6_LIMBS(data, 6, a, 6, b, 6, low)
            ADD_6_LIMBS(data, 12, c, 0, d, 0, low)
            ADD_6_LIMBS(data, 18, c, 6, d, 6, low)
            :   [low]"=&r"(low)
            :   [data]"r"(&data),
                [a]"r"(a),
                [b]"r"(b),
                [c]"r"(c),
                [d]"r"(d)
            : "cc", "memory"
        );
        asm volatile(
            SCHOOLBOOK_6x6(z, 0, data, 0, data, 12) // z[0] = x[0] * y[0]
            SCHOOLBOOK_6x6(data, 24, data, 6, data, 18) // tmp = x[1] * y[1]
            SUB_12_LIMBS_MOD(z, 0, z, 0, data, 24, d0, d1, d2, d3, d4, d5)
            SCHOOLBOOK_6x6(z, 12, data, 0, data, 18) // ad
            SCHOOLBOOK_6x6(data, 24, data, 6, data, 12) // bc
            ADD_12_LIMBS(z, 12, z, 12, data, 24, low)
            :   [low]"=&r"(low),
                [high]"=&r"(high),
                [zero]"=&r"(zero),
                [d0]"=&r"(d0),
                [d1]"=&r"(d1),
                [d2]"=&r"(d2),
                [d3]"=&r"(d3),
                [d4]"=&r"(d4),
                [d5]"=&r"(d5)
            :   [data]"r"(&data),
                [z]"r"(z),
                [p0]"m"(p0),
                [p1]"m"(p1),
                [p2]"m"(p2),
                [p3]"m"(p3),
                [p4]"m"(p4),
                [p5]"m"(p5)
            : "rdx", "cc", "memory"
        );
    }
}

#undef STR_IMPL
#undef STR
#undef PTR
#undef PTR2
#undef D4
#undef D6
#undef SCHOOLBOOK_4x4_ROUND
#undef SCHOOLBOOK_4x4
#undef SCHOOLBOOK_6x6_ROUND
#undef SCHOOLBOOK_6x6
#undef DOUBLE_MOD_P
#undef MUL_8_LIMBS_BY_9
#undef GET_MODULUS_4_LIMBS
#undef GET_MODULUS_6_LIMBS
#undef T8
#undef T12
#undef REDC_8_LIMBS_LOAD_NEXT
#undef REDC_8_LIMBS_CANCEL_LOW
#undef REDC_12_LIMBS_LOAD_NEXT
#undef REDC_12_LIMBS_CANCEL_LOW
#undef ADD_4_LIMBS
#undef ADD_4_LIMBS_MOD
#undef ADD_6_LIMBS
#undef ADD_6_LIMBS_MOD
#undef ADD_8_LIMBS
#undef ADD_8_LIMBS_MOD
#undef ADD_12_LIMBS
#undef ADD_12_LIMBS_MOD
#undef SUB_8_LIMBS
#undef SUB_8_LIMBS_MOD
#undef SUB_12_LIMBS
#undef SUB_12_LIMBS_MOD

// clang-format on
