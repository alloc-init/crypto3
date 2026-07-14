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

#define SET_STATIC_MODULUS_FROM_FIELD()                             \
    constexpr auto mod_obj = Field::modulus_params.get_mod_obj();   \
    static constexpr limb p0 = limb(mod_obj.get_mod().limbs()[0]);  \
    static constexpr limb p1 = limb(mod_obj.get_mod().limbs()[1]);  \
    static constexpr limb p2 = limb(mod_obj.get_mod().limbs()[2]);  \
    static constexpr limb p3 = limb(mod_obj.get_mod().limbs()[3]);  \
    static constexpr limb p_dash = limb(mod_obj.get_p_dash());

// Get the (i + j) mod 5 t register in the REDC rotating window.
#define T(I, J) "%[t" STR(BOOST_PP_MOD(BOOST_PP_ADD(I, J), 5)) "]"

#define MONTGOMERY_REDUCE_LOAD_NEXT(I)                  \
    "setc %b[pending]\n"                                \
    "movq " PTR2(data, I, 5) ", " T(I, 5) "\n"

// main body of loop in montgomery reduce
#define MONTGOMERY_REDUCE_CANCEL_LOW(I)                 \
    /* m = data[i] * p_dash */                          \
    "movq " T(I, 0) ", %%rdx\n"                         \
    "imulq %[p_dash], %%rdx\n"                          \
    "xor %[zero], %[zero]\n"                            \
    /* multiply m by each modulus limb and add into the rotating window */ \
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
    /* merge OF, prior pending carry, and CF into the next live limb */ \
    "adox %[zero], %[high]\n"                           \
    "adox %[pending], %[high]\n"                        \
    "adcx %[high], " T(I, 4) "\n"                       \
    BOOST_PP_IF(BOOST_PP_LESS(I, 3), MONTGOMERY_REDUCE_LOAD_NEXT(I), "")

#define ADD_LOW_4_LIMBS(Z, Z_BASE, X, X_BASE, Y, Y_BASE, TMP) \
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

#define ADD_LIMBS(Z, Z_BASE, X, X_BASE, Y, Y_BASE, TMP) \
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
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 5) "\n"        \
    "movq " PTR2(X, X_BASE, 6) ", %[" #TMP "]\n"        \
    "adcq " PTR2(Y, Y_BASE, 6) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 6) "\n"        \
    "movq " PTR2(X, X_BASE, 7) ", %[" #TMP "]\n"        \
    "adcq " PTR2(Y, Y_BASE, 7) ", %[" #TMP "]\n"        \
    "movq %[" #TMP "], " PTR2(Z, Z_BASE, 7) "\n"

#define ADD_LOW_4_LIMBS_MOD(Z, Z_BASE, X, X_BASE, Y, Y_BASE)    \
    "movq " PTR2(X, X_BASE, 0) ", %[t0]\n"                      \
    "movq " PTR2(X, X_BASE, 1) ", %[t1]\n"                      \
    "movq " PTR2(X, X_BASE, 2) ", %[t2]\n"                      \
    "movq " PTR2(X, X_BASE, 3) ", %[t3]\n"                      \
    "movq " PTR2(Y, Y_BASE, 0) ", %[q0]\n"                      \
    "addq %[q0], %[t0]\n"                                       \
    "movq " PTR2(Y, Y_BASE, 1) ", %[q0]\n"                      \
    "adcq %[q0], %[t1]\n"                                       \
    "movq " PTR2(Y, Y_BASE, 2) ", %[q0]\n"                      \
    "adcq %[q0], %[t2]\n"                                       \
    "movq " PTR2(Y, Y_BASE, 3) ", %[q0]\n"                      \
    "adcq %[q0], %[t3]\n"                                       \
    "movq %[t0], %[q0]\n"                                       \
    "movq %[t1], %[q1]\n"                                       \
    "movq %[t2], %[q2]\n"                                       \
    "movq %[t3], %[q3]\n"                                       \
    "subq %[p0], %[q0]\n"                                       \
    "sbbq %[p1], %[q1]\n"                                       \
    "sbbq %[p2], %[q2]\n"                                       \
    "sbbq %[p3], %[q3]\n"                                       \
    "cmovnc %[q0], %[t0]\n"                                     \
    "cmovnc %[q1], %[t1]\n"                                     \
    "cmovnc %[q2], %[t2]\n"                                     \
    "cmovnc %[q3], %[t3]\n"                                     \
    "movq %[t0], " PTR2(Z, Z_BASE, 0) "\n"                      \
    "movq %[t1], " PTR2(Z, Z_BASE, 1) "\n"                      \
    "movq %[t2], " PTR2(Z, Z_BASE, 2) "\n"                      \
    "movq %[t3], " PTR2(Z, Z_BASE, 3) "\n"

#define ADD_LIMBS_MOD(Z, Z_BASE, X, X_BASE, Y, Y_BASE, T0, T1, T2, T3, Q0, Q1, Q2, Q3) \
    "movq " PTR2(X, X_BASE, 4) ", %[" #T0 "]\n"                                        \
    "movq " PTR2(X, X_BASE, 5) ", %[" #T1 "]\n"                                        \
    "movq " PTR2(X, X_BASE, 6) ", %[" #T2 "]\n"                                        \
    "movq " PTR2(X, X_BASE, 7) ", %[" #T3 "]\n"                                        \
    "movq " PTR2(X, X_BASE, 0) ", %[" #Q0 "]\n"                                        \
    "addq " PTR2(Y, Y_BASE, 0) ", %[" #Q0 "]\n"                                        \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 0) "\n"                                        \
    "movq " PTR2(X, X_BASE, 1) ", %[" #Q0 "]\n"                                        \
    "adcq " PTR2(Y, Y_BASE, 1) ", %[" #Q0 "]\n"                                        \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 1) "\n"                                        \
    "movq " PTR2(X, X_BASE, 2) ", %[" #Q0 "]\n"                                        \
    "adcq " PTR2(Y, Y_BASE, 2) ", %[" #Q0 "]\n"                                        \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 2) "\n"                                        \
    "movq " PTR2(X, X_BASE, 3) ", %[" #Q0 "]\n"                                        \
    "adcq " PTR2(Y, Y_BASE, 3) ", %[" #Q0 "]\n"                                        \
    "movq %[" #Q0 "], " PTR2(Z, Z_BASE, 3) "\n"                                        \
    "movq " PTR2(Y, Y_BASE, 4) ", %[" #Q0 "]\n"                                        \
    "adcq %[" #Q0 "], %[" #T0 "]\n"                                                    \
    "movq " PTR2(Y, Y_BASE, 5) ", %[" #Q0 "]\n"                                        \
    "adcq %[" #Q0 "], %[" #T1 "]\n"                                                    \
    "movq " PTR2(Y, Y_BASE, 6) ", %[" #Q0 "]\n"                                        \
    "adcq %[" #Q0 "], %[" #T2 "]\n"                                                    \
    "movq " PTR2(Y, Y_BASE, 7) ", %[" #Q0 "]\n"                                        \
    "adcq %[" #Q0 "], %[" #T3 "]\n"                                                    \
    "movq %[" #T0 "], %[" #Q0 "]\n"                                                    \
    "movq %[" #T1 "], %[" #Q1 "]\n"                                                    \
    "movq %[" #T2 "], %[" #Q2 "]\n"                                                    \
    "movq %[" #T3 "], %[" #Q3 "]\n"                                                    \
    "subq %[p0], %[" #Q0 "]\n"                                                         \
    "sbbq %[p1], %[" #Q1 "]\n"                                                         \
    "sbbq %[p2], %[" #Q2 "]\n"                                                         \
    "sbbq %[p3], %[" #Q3 "]\n"                                                         \
    "cmovnc %[" #Q0 "], %[" #T0 "]\n"                                                  \
    "cmovnc %[" #Q1 "], %[" #T1 "]\n"                                                  \
    "cmovnc %[" #Q2 "], %[" #T2 "]\n"                                                  \
    "cmovnc %[" #Q3 "], %[" #T3 "]\n"                                                  \
    "movq %[" #T0 "], " PTR2(Z, Z_BASE, 4) "\n"                                        \
    "movq %[" #T1 "], " PTR2(Z, Z_BASE, 5) "\n"                                        \
    "movq %[" #T2 "], " PTR2(Z, Z_BASE, 6) "\n"                                        \
    "movq %[" #T3 "], " PTR2(Z, Z_BASE, 7) "\n"

#define SUB_LIMBS(RESULT, RESULT_BASE, OTHER, OTHER_BASE, SCRATCH) \
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

#define SUB_LIMBS_MOD(Z, Z_BASE, X, X_BASE, Y, Y_BASE, T0, T1, T2, T3) \
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

namespace nil::crypto3::algebra::fields::detail::fp12_fast {
    inline void multiply_4x4_x86(limb *z, const limb *x, const limb *y) {
        limb low, zero, high;
        limb d0, d1, d2, d3;
        asm volatile(
            SCHOOLBOOK_4x4(z, 0, x, 0, y, 0)
            : [low]"=&r"(low),
              [high]"=&r"(high),
              [zero]"=&r"(zero),
              [d0]"=&r"(d0),
              [d1]"=&r"(d1),
              [d2]"=&r"(d2),
              [d3]"=&r"(d3)
            : [z]"r"(z),
              [y]"r"(y),
              [x]"r"(x)
            : "rdx", "cc", "memory"
        );
    }

    inline void multiply_6x6_x86(limb *z, const limb *x, const limb *y) {
        limb low, zero, high;
        limb d0, d1, d2, d3;
        asm volatile(
            SCHOOLBOOK_4x4(z, 0, x, 0, y, 0)
            : [low]"=&r"(low),
              [high]"=&r"(high),
              [zero]"=&r"(zero),
              [d0]"=&r"(d0),
              [d1]"=&r"(d1),
              [d2]"=&r"(d2),
              [d3]"=&r"(d3)
            : [z]"r"(z),
              [y]"r"(y),
              [x]"r"(x)
            : "rdx", "cc", "memory"
        );
    }

    template<class Field>
    inline void montgomery_reduce_x86(limb *result, const limb *data) {
        SET_STATIC_MODULUS_FROM_FIELD();
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
            MONTGOMERY_REDUCE_CANCEL_LOW(0)
            MONTGOMERY_REDUCE_CANCEL_LOW(1)
            MONTGOMERY_REDUCE_CANCEL_LOW(2)
            MONTGOMERY_REDUCE_CANCEL_LOW(3)

            // Conditionally reduce the high limbs modulo p, reusing scratch registers.
            // q = t - p
            "movq " T(0, 4) ", %[low]\n"
            "movq " T(1, 4) ", %[high]\n"
            "movq " T(2, 4) ", %[zero]\n"
            "movq " T(3, 4) ", %%rdx\n"

            // Try q = t - p.
            "subq %[p0], %[low]\n"
            "sbbq %[p1], %[high]\n"
            "sbbq %[p2], %[zero]\n"
            "sbbq %[p3], %%rdx\n"

            // If t - p did not borrow, keep q as the canonical result.
            "cmovnc %[low], " T(0, 4) "\n"
            "cmovnc %[high], " T(1, 4) "\n"
            "cmovnc %[zero], " T(2, 4) "\n"
            "cmovnc %%rdx, " T(3, 4) "\n"

            "movq " T(0, 4) ", " PTR(result, 0) "\n"
            "movq " T(1, 4) ", " PTR(result, 1) "\n"
            "movq " T(2, 4) ", " PTR(result, 2) "\n"
            "movq " T(3, 4) ", " PTR(result, 3) "\n"

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

    inline void subtract_8_limbs_x86(limb *result, const limb *other) {
        limb scratch;
        asm volatile(
            SUB_LIMBS(result, 0, other, 0, scratch)
            : [scratch]"=&r"(scratch)
            : [result]"r"(result),
              [other]"r"(other)
            : "cc", "memory"
        );
    }

    template<class Field>
    inline void add_8_limbs_mod_x86(limb *z, const limb *x, const limb *y) {
        SET_STATIC_MODULUS_FROM_FIELD();
        limb t0, t1, t2, t3, q0, q1, q2, q3;
        asm volatile(
            ADD_LIMBS_MOD(z, 0, x, 0, y, 0, t0, t1, t2, t3, q0, q1, q2, q3)
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

    template<class Field>
    inline void subtract_8_limbs_mod_x86(limb *z, const limb *x, const limb *y) {
        SET_STATIC_MODULUS_FROM_FIELD();
        limb t0, t1, t2, t3;
        asm volatile(
            SUB_LIMBS_MOD(z, 0, x, 0, y, 0, t0, t1, t2, t3)
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

    template<class Field>
    inline void mul_8_limbs_by_9_x86(limb *dst, const limb *src) {
        SET_STATIC_MODULUS_FROM_FIELD();
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

    template <class Field>
    inline void fp2_base_add_mod_x86(limb *z, const limb *x, const limb *y) {
        SET_STATIC_MODULUS_FROM_FIELD();
        limb t0, t1, t2, t3;
        limb q0, q1, q2, q3;
        asm volatile(
            ADD_LOW_4_LIMBS_MOD(z, 0, x, 0, y, 0)
            ADD_LOW_4_LIMBS_MOD(z, 4, x, 4, y, 4)
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

    template<class Field>
    inline void fp2_sub_pre_x86(limb *data, const limb *other) {
        SET_STATIC_MODULUS_FROM_FIELD();
        limb t0, t1, t2, t3;
        asm volatile(
            SUB_LIMBS_MOD(data, 0, data, 0, other, 0, t0, t1, t2, t3)
            SUB_LIMBS(data, 8, other, 8, t0)
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

    template<class Field>
    inline void fp2_mul_pre_x86(limb *z, const limb *x, const limb *y) {
        SET_STATIC_MODULUS_FROM_FIELD();
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
            SUB_LIMBS_MOD(z, 0, z, 0, scratch, 0, d0, d1, d2, d3)
            SCHOOLBOOK_4x4(z, 8, x, 0, y, 4)
            SCHOOLBOOK_4x4(scratch, 0, x, 4, y, 0)
            ADD_LIMBS(z, 8, z, 8, scratch, 0, low)
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

    template<class Field>
    inline void fp2_add_mul_pre_x86(limb *z, const limb  *a, const limb  *b, const limb  *c, const limb  *d) {
        SET_STATIC_MODULUS_FROM_FIELD();
        limb low, high, zero, d0, d1, d2, d3;
        struct {
            limb x[8];
            limb y[8];
            limb tmp[8];
        } data;
        asm volatile(
            ADD_LOW_4_LIMBS(data, 0, a, 0, b, 0, low)
            ADD_LOW_4_LIMBS(data, 4, a, 4, b, 4, low)
            ADD_LOW_4_LIMBS(data, 8, c, 0, d, 0, low)
            ADD_LOW_4_LIMBS(data, 12, c, 4, d, 4, low)
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
            SUB_LIMBS_MOD(z, 0, z, 0, data, 16, d0, d1, d2, d3)
            SCHOOLBOOK_4x4(z, 8, data, 0, data, 12)
            SCHOOLBOOK_4x4(data, 16, data, 4, data, 8)
            ADD_LIMBS(z, 8, z, 8, data, 16, low)
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
}

#undef STR_IMPL
#undef STR
#undef PTR
#undef PTR2
#undef D4
#undef SCHOOLBOOK_ROUND
#undef SCHOOLBOOK_4x4
#undef DOUBLE_MOD_P
#undef MUL_8_LIMBS_BY_9
#undef SET_STATIC_MODULUS_FROM_FIELD
#undef T
#undef MONTGOMERY_REDUCE_LOAD_NEXT
#undef MONTGOMERY_REDUCE_CANCEL_LOW
#undef ADD_LOW_4_LIMBS
#undef ADD_LIMBS
#undef ADD_LOW_4_LIMBS_MOD
#undef ADD_LIMBS_MOD
#undef SUB_LIMBS
#undef SUB_LIMBS_MOD
#undef BRANCHY_SUB_LIMBS_MOD

// clang-format on
