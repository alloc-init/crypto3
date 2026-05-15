///////////////////////////////////////////////////////////////
//  Copyright (c) 2023 Martun Karapetyan <martun@nil.foundation>
//
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt
//
//  Contains eval_multiply for cpp_int_modular_backend.
//

#ifndef CRYPTO3_MP_CPP_INT_MUL_HPP
#define CRYPTO3_MP_CPP_INT_MUL_HPP

// This implementation header is normally included from cpp_int_modular.hpp after
// cpp_int_modular_backend and its traits are declared. Include the parent when a
// parser opens this file directly, so editor tooling sees the required context.
#ifndef CRYPTO3_CPP_INT_MODULAR_HPP
#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>
#endif

namespace boost {
    namespace multiprecision {
        namespace backends {
            // Functions in this file are used by modular parameter construction, reduction, and
            // fixed-width lazy arithmetic. Keep them on cpp_int_modular_backend limbs instead of
            // converting through cpp_int.

            inline constexpr void mul_4x4_add(
                    limb_type &acc0,
                    limb_type &acc1,
                    limb_type &acc2,
                    limb_type x,
                    limb_type y) noexcept {
                constexpr std::size_t limb_bits =
                    sizeof(limb_type) * CHAR_BIT;
                const double_limb_type product =
                    static_cast<double_limb_type>(x) *
                    static_cast<double_limb_type>(y);
                const double_limb_type sum0 =
                    static_cast<double_limb_type>(acc0) +
                    static_cast<limb_type>(product);
                acc0 = static_cast<limb_type>(sum0);

                const double_limb_type sum1 =
                    static_cast<double_limb_type>(acc1) +
                    static_cast<limb_type>(product >> limb_bits) +
                    (sum0 >> limb_bits);
                acc1 = static_cast<limb_type>(sum1);
                acc2 += static_cast<limb_type>(
                    sum1 >> limb_bits);
            }

            template<unsigned Bits>
            inline constexpr void mul_4x4_emit(
                    cpp_int_modular_backend<Bits> &result,
                    std::size_t idx,
                    limb_type &acc0,
                    limb_type &acc1,
                    limb_type &acc2) noexcept {
                result.limbs()[idx] = acc0;
                acc0 = acc1;
                acc1 = acc2;
                acc2 = 0u;
            }

            template<unsigned Bits>
            inline constexpr void eval_multiply_4x4(
                    cpp_int_modular_backend<Bits> &result,
                    const limb_type *x,
                    const limb_type *y) noexcept {
                result.zero_after(0u);

                limb_type acc0 = 0u;
                limb_type acc1 = 0u;
                limb_type acc2 = 0u;

                mul_4x4_add(acc0, acc1, acc2, x[0], y[0]);
                mul_4x4_emit(result, 0u, acc0, acc1, acc2);

                mul_4x4_add(acc0, acc1, acc2, x[0], y[1]);
                mul_4x4_add(acc0, acc1, acc2, x[1], y[0]);
                mul_4x4_emit(result, 1u, acc0, acc1, acc2);

                mul_4x4_add(acc0, acc1, acc2, x[0], y[2]);
                mul_4x4_add(acc0, acc1, acc2, x[1], y[1]);
                mul_4x4_add(acc0, acc1, acc2, x[2], y[0]);
                mul_4x4_emit(result, 2u, acc0, acc1, acc2);

                mul_4x4_add(acc0, acc1, acc2, x[0], y[3]);
                mul_4x4_add(acc0, acc1, acc2, x[1], y[2]);
                mul_4x4_add(acc0, acc1, acc2, x[2], y[1]);
                mul_4x4_add(acc0, acc1, acc2, x[3], y[0]);
                mul_4x4_emit(result, 3u, acc0, acc1, acc2);

                mul_4x4_add(acc0, acc1, acc2, x[1], y[3]);
                mul_4x4_add(acc0, acc1, acc2, x[2], y[2]);
                mul_4x4_add(acc0, acc1, acc2, x[3], y[1]);
                mul_4x4_emit(result, 4u, acc0, acc1, acc2);

                mul_4x4_add(acc0, acc1, acc2, x[2], y[3]);
                mul_4x4_add(acc0, acc1, acc2, x[3], y[2]);
                mul_4x4_emit(result, 5u, acc0, acc1, acc2);

                mul_4x4_add(acc0, acc1, acc2, x[3], y[3]);
                mul_4x4_emit(result, 6u, acc0, acc1, acc2);
                result.limbs()[7] = acc0;
                result.limbs()[8] = acc1;
                result.zero_after(9u);
                result.set_carry(false);
                result.normalize();
            }

            template<std::size_t InputLimbs, unsigned Bits>
            inline constexpr void eval_multiply_low_limbs(
                    cpp_int_modular_backend<Bits> &result,
                    const limb_type *x,
                    const limb_type *y) noexcept {
                result.zero_after(0u);
                result.set_carry(false);

                for (std::size_t i = 0; i < InputLimbs; ++i) {
                    limb_type carry = 0u;
                    for (std::size_t j = 0; j < InputLimbs && i + j < result.size(); ++j) {
                        const double_limb_type product =
                            static_cast<double_limb_type>(x[i]) *
                                static_cast<double_limb_type>(y[j]) +
                            static_cast<double_limb_type>(result.limbs()[i + j]) +
                            carry;
                        result.limbs()[i + j] = static_cast<limb_type>(product);
                        carry = static_cast<limb_type>(
                            product >> cpp_int_modular_backend<Bits>::limb_bits);
                    }

                    for (std::size_t idx = i + InputLimbs; carry != 0u && idx < result.size(); ++idx) {
                        const double_limb_type sum =
                            static_cast<double_limb_type>(result.limbs()[idx]) +
                            carry;
                        result.limbs()[idx] = static_cast<limb_type>(sum);
                        carry = static_cast<limb_type>(
                            sum >> cpp_int_modular_backend<Bits>::limb_bits);
                    }
                }

                result.normalize();
            }

            template<unsigned Bits1, unsigned Bits2>
            inline constexpr void
            eval_multiply(cpp_int_modular_backend<Bits1 + Bits2> &result,
                          const cpp_int_modular_backend<Bits1> &a,
                          const cpp_int_modular_backend<Bits2> &b) noexcept {
                result = a;
                // Call the lower function, we don't care about speed here.
                eval_multiply(result, b);
            }

            // If the second argument is trivial or not, we still assign it to the first the exact same way.
            template<unsigned Bits1, unsigned Bits2>
            inline constexpr typename std::enable_if<
                !is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits1>>::value
                >::type
            eval_multiply(cpp_int_modular_backend<Bits1> &result,
                          const cpp_int_modular_backend<Bits2> &a,
                          const limb_type &b) noexcept {
                result = a;
                // Call the lower function, we don't care about speed here.
                eval_multiply(result, b);
            }


            template<unsigned Bits1>
            inline constexpr typename std::enable_if<
                !is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits1>>::value>::type
            eval_multiply(cpp_int_modular_backend<Bits1> &result,
                          const limb_type &b) noexcept {
                limb_type carry = 0u;
                for (std::size_t i = 0; i < result.size(); ++i) {
                    const double_limb_type product =
                        static_cast<double_limb_type>(result.limbs()[i]) *
                            static_cast<double_limb_type>(b) +
                        carry;
                    result.limbs()[i] = static_cast<limb_type>(product);
                    carry = static_cast<limb_type>(
                        product >> cpp_int_modular_backend<Bits1>::limb_bits);
                }
                result.set_carry(false);
                result.normalize();
            }

            // Caller is responsible for the result to fit in Bits1 bits, we will NOT throw!!!
            // Covers the case where the second operand is not trivial.
            template<unsigned Bits1, unsigned Bits2>
            inline constexpr typename std::enable_if<
                !is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits1>>::value &&
                !is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits2>>::value>::type
            eval_multiply(cpp_int_modular_backend<Bits1> &result,
                          const cpp_int_modular_backend<Bits2> &a) noexcept {
                if (result.size() >= 9u && a.size() >= 4u) {
                    bool fits_4x4 = true;
                    for (std::size_t i = 4u; fits_4x4 && i < result.size(); ++i) {
                        fits_4x4 = result.limbs()[i] == 0u;
                    }
                    for (std::size_t i = 4u; fits_4x4 && i < a.size(); ++i) {
                        fits_4x4 = a.limbs()[i] == 0u;
                    }

                    if (fits_4x4) {
                        limb_type lhs4[4] = {
                            result.limbs()[0], result.limbs()[1], result.limbs()[2], result.limbs()[3]};
                        limb_type rhs4[4] = {
                            a.limbs()[0], a.limbs()[1], a.limbs()[2], a.limbs()[3]};
                        eval_multiply_4x4(result, lhs4, rhs4);
                        return;
                    }
                }

                limb_type lhs[cpp_int_modular_backend<Bits1>::internal_limb_count] = {};
                limb_type rhs[cpp_int_modular_backend<Bits2>::internal_limb_count] = {};

                for (std::size_t i = 0; i < result.size(); ++i) {
                    lhs[i] = result.limbs()[i];
                }
                for (std::size_t i = 0; i < a.size(); ++i) {
                    rhs[i] = a.limbs()[i];
                }

                std::size_t lhs_size = result.size();
                while (lhs_size > 0u && lhs[lhs_size - 1u] == 0u) {
                    --lhs_size;
                }
                std::size_t rhs_size = a.size();
                while (rhs_size > 0u && rhs[rhs_size - 1u] == 0u) {
                    --rhs_size;
                }

                result.zero_after(0u);
                result.set_carry(false);

                for (std::size_t i = 0; i < lhs_size; ++i) {
                    if (lhs[i] == 0u) {
                        continue;
                    }

                    limb_type carry = 0u;
                    std::size_t j = 0;
                    for (; j < rhs_size && i + j < result.size(); ++j) {
                        const double_limb_type product =
                            static_cast<double_limb_type>(lhs[i]) *
                                static_cast<double_limb_type>(rhs[j]) +
                            static_cast<double_limb_type>(result.limbs()[i + j]) +
                            carry;
                        result.limbs()[i + j] = static_cast<limb_type>(product);
                        carry = static_cast<limb_type>(
                            product >> cpp_int_modular_backend<Bits1>::limb_bits);
                    }

                    for (std::size_t idx = i + j; carry != 0u && idx < result.size(); ++idx) {
                        const double_limb_type sum =
                            static_cast<double_limb_type>(result.limbs()[idx]) +
                            carry;
                        result.limbs()[idx] = static_cast<limb_type>(sum);
                        carry = static_cast<limb_type>(
                            sum >> cpp_int_modular_backend<Bits1>::limb_bits);
                    }
                }

                result.normalize();
            }

            template<unsigned Bits1, unsigned Bits2>
            inline constexpr typename std::enable_if<
                is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits1>>::value &&
                !is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits2>>::value>::type
            eval_multiply(cpp_int_modular_backend<Bits1> &result,
                          const cpp_int_modular_backend<Bits2> &a) noexcept {
                cpp_int_backend<Bits1, Bits1, unsigned_magnitude, unchecked> result_cpp_int = result.to_cpp_int();
                eval_multiply(result_cpp_int, a.to_cpp_int());
                result.from_cpp_int(result_cpp_int);
            }

            // We need to specialize this for goldilock fields, where the second argument is trivial, while the first may not.
            // Caller is responsible for the result to fit in "Bits1" bits, we WILL NOT THROW!!!
            template<unsigned Bits1, unsigned Bits2>
            inline constexpr typename std::enable_if<
                !is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits1>>::value &&
                is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits2>>::value>::type
            eval_multiply(
                cpp_int_modular_backend<Bits1>& result,
                const cpp_int_modular_backend<Bits2>& o) noexcept
            {
                eval_multiply(result, limb_type(*o.limbs()));
            }

            // It looks to me boost has a bug with handling 'eval_multiply' for both arguments being trivial but of
            // different bit length, so we need to specialize this one here.
            // Caller is responsible for the result to fit in "Bits1" bits, we WILL NOT THROW!!!
            template<unsigned Bits1, unsigned Bits2>
            inline constexpr typename std::enable_if<
                is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits1>>::value &&
                is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits2>>::value>::type
            eval_multiply(
                cpp_int_modular_backend<Bits1>& result,
                const cpp_int_modular_backend<Bits2>& o) noexcept
            {
               *result.limbs() *= *o.limbs();
            }

            // Multiplication with an unsigned integral type.
            template<unsigned Bits1, unsigned Bits2>
            inline constexpr typename std::enable_if<
                is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits1>>::value &&
                is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits2>>::value
                >::type
            eval_multiply(cpp_int_modular_backend<Bits1> &result,
                          const cpp_int_modular_backend<Bits2> &a,
                          const limb_type &b) noexcept {
                result = a;
                // Call the lower function, we don't care about speed here.
                eval_multiply(result, b);
            }

            // Multiplication with an unsigned integral type.
            template<unsigned Bits1>
            inline constexpr typename std::enable_if<
                is_trivial_cpp_int_modular<cpp_int_modular_backend<Bits1>>::value>::type
            eval_multiply(cpp_int_modular_backend<Bits1> &result,
                          const limb_type &b) noexcept {
                *result.limbs() *= b;
            }


        }    // namespace backends
    }   // namespace multiprecision
}   // namespace boost


#endif 
