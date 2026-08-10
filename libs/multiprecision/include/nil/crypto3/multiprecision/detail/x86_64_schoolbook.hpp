//---------------------------------------------------------------------------//
// Copyright (c) 2026 Riccardo Abbate
//
// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_MULTIPRECISION_DETAIL_X86_64_SCHOOLBOOK_HPP
#define CRYPTO3_MULTIPRECISION_DETAIL_X86_64_SCHOOLBOOK_HPP

#include <cstddef>
#include <type_traits>

#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/mod.hpp>

// clang-format off

// These macros form an assembly fragment rather than a complete asm statement so fp12_fast can combine several
// products in one block. Keep the operands named low, high, zero, and d0..d3 in every caller.
#define CRYPTO3_MP_DETAIL_STR_IMPL(X) #X
#define CRYPTO3_MP_DETAIL_STR(X) CRYPTO3_MP_DETAIL_STR_IMPL(X)
#define CRYPTO3_MP_DETAIL_PTR(REGNAME, I) CRYPTO3_MP_DETAIL_STR(I) "*8(%[" #REGNAME "])"
#define CRYPTO3_MP_DETAIL_PTR2(REGNAME, I, J) CRYPTO3_MP_DETAIL_PTR(REGNAME, BOOST_PP_ADD(I, J))
#define CRYPTO3_MP_DETAIL_D4(I, J)                                                                                   \
    "%[d" CRYPTO3_MP_DETAIL_STR(BOOST_PP_MOD(BOOST_PP_ADD(I, J), 4)) "]"

#define CRYPTO3_MP_DETAIL_SCHOOLBOOK_4X4_ROUND(ROUND, Z, Z_BASE, X, X_BASE, Y, Y_BASE)                              \
    "mov " CRYPTO3_MP_DETAIL_PTR2(Y, Y_BASE, ROUND) ", %%rdx\n"                                                    \
    "mulx " CRYPTO3_MP_DETAIL_PTR2(X, X_BASE, 0) ", %[low], %[high]\n"                                                 \
    "adox %[low], " CRYPTO3_MP_DETAIL_D4(0, ROUND) "\n"                                                          \
    "mov " CRYPTO3_MP_DETAIL_D4(0, ROUND) ", " CRYPTO3_MP_DETAIL_PTR2(Z, Z_BASE, ROUND) "\n"                    \
    "adcx %[high], " CRYPTO3_MP_DETAIL_D4(1, ROUND) "\n"                                                         \
    "mulx " CRYPTO3_MP_DETAIL_PTR2(X, X_BASE, 1) ", %[low], %[high]\n"                                                 \
    "adox %[low], " CRYPTO3_MP_DETAIL_D4(1, ROUND) "\n"                                                          \
    "adcx %[high], " CRYPTO3_MP_DETAIL_D4(2, ROUND) "\n"                                                         \
    "mulx " CRYPTO3_MP_DETAIL_PTR2(X, X_BASE, 2) ", %[low], %[high]\n"                                                 \
    "adox %[low], " CRYPTO3_MP_DETAIL_D4(2, ROUND) "\n"                                                          \
    "adcx %[high], " CRYPTO3_MP_DETAIL_D4(3, ROUND) "\n"                                                         \
    "mulx " CRYPTO3_MP_DETAIL_PTR2(X, X_BASE, 3) ", %[low], " CRYPTO3_MP_DETAIL_D4(4, ROUND) "\n"              \
    "adox %[low], " CRYPTO3_MP_DETAIL_D4(3, ROUND) "\n"                                                          \
    "adox %[zero], " CRYPTO3_MP_DETAIL_D4(4, ROUND) "\n"                                                         \
    "adcx %[zero], " CRYPTO3_MP_DETAIL_D4(4, ROUND) "\n"

#define CRYPTO3_MP_DETAIL_SCHOOLBOOK_4X4(Z, Z_BASE, X, X_BASE, Y, Y_BASE)                                           \
    "xor %[zero], %[zero]\n"                                                                                       \
    "mov " CRYPTO3_MP_DETAIL_PTR2(Y, Y_BASE, 0) ", %%rdx\n"                                                       \
    "mulx " CRYPTO3_MP_DETAIL_PTR2(X, X_BASE, 0) ", %[d0], %[d1]\n"                                                  \
    "mov %[d0], " CRYPTO3_MP_DETAIL_PTR2(Z, Z_BASE, 0) "\n"                                                       \
    "mulx " CRYPTO3_MP_DETAIL_PTR2(X, X_BASE, 1) ", %[low], %[d2]\n"                                                 \
    "add %[low], %[d1]\n"                                                                                          \
    "mulx " CRYPTO3_MP_DETAIL_PTR2(X, X_BASE, 2) ", %[low], %[d3]\n"                                                 \
    "adc %[low], %[d2]\n"                                                                                          \
    "mulx " CRYPTO3_MP_DETAIL_PTR2(X, X_BASE, 3) ", %[low], %[d0]\n"                                                 \
    "adc %[low], %[d3]\n"                                                                                          \
    "adc $0, %[d0]\n"                                                                                              \
    CRYPTO3_MP_DETAIL_SCHOOLBOOK_4X4_ROUND(1, Z, Z_BASE, X, X_BASE, Y, Y_BASE)                                     \
    CRYPTO3_MP_DETAIL_SCHOOLBOOK_4X4_ROUND(2, Z, Z_BASE, X, X_BASE, Y, Y_BASE)                                     \
    CRYPTO3_MP_DETAIL_SCHOOLBOOK_4X4_ROUND(3, Z, Z_BASE, X, X_BASE, Y, Y_BASE)                                     \
    "mov " CRYPTO3_MP_DETAIL_D4(0, 4) ", " CRYPTO3_MP_DETAIL_PTR2(Z, Z_BASE, 4) "\n"                             \
    "mov " CRYPTO3_MP_DETAIL_D4(1, 4) ", " CRYPTO3_MP_DETAIL_PTR2(Z, Z_BASE, 5) "\n"                             \
    "mov " CRYPTO3_MP_DETAIL_D4(2, 4) ", " CRYPTO3_MP_DETAIL_PTR2(Z, Z_BASE, 6) "\n"                             \
    "mov " CRYPTO3_MP_DETAIL_D4(3, 4) ", " CRYPTO3_MP_DETAIL_PTR2(Z, Z_BASE, 7) "\n"

// clang-format on

namespace nil {
    namespace crypto3 {
        namespace multiprecision {
            namespace detail {

                /** Multiply two four-limb integers with the ADX/BMI2 kernel shared with fp12_fast. */
                template<typename Limb>
                inline void schoolbook_4x4_adx_bmi2(Limb *result, const Limb *left, const Limb *right) {
                    static_assert(std::is_unsigned<Limb>::value, "the assembly kernel expects unsigned limbs");
                    static_assert(sizeof(Limb) == 8, "the assembly kernel expects 64-bit limbs");

                    Limb low, high, zero, d0, d1, d2, d3;
                    asm volatile(CRYPTO3_MP_DETAIL_SCHOOLBOOK_4X4(result, 0, left, 0, right, 0)
                                 : [low] "=&r"(low), [high] "=&r"(high), [zero] "=&r"(zero), [d0] "=&r"(d0),
                                   [d1] "=&r"(d1), [d2] "=&r"(d2), [d3] "=&r"(d3)
                                 : [result] "r"(result), [left] "r"(left), [right] "r"(right)
                                 : "rdx", "cc", "memory");
                }

            }    // namespace detail
        }    // namespace multiprecision
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MULTIPRECISION_DETAIL_X86_64_SCHOOLBOOK_HPP
