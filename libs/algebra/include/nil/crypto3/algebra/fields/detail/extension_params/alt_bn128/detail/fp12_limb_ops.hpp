#pragma once

#include <climits>
#include <cstddef>

#include <boost/multiprecision/cpp_int.hpp>

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {
                namespace detail {
                    namespace alt_bn128_fp12_limb_ops {

                        typedef boost::multiprecision::limb_type limb_type;
                        typedef boost::multiprecision::double_limb_type double_limb_type;

                        inline constexpr void multiply_4x4_add(limb_type &acc0, limb_type &acc1, limb_type &acc2,
                                                               limb_type x, limb_type y) noexcept {
                            constexpr std::size_t limb_bits = sizeof(limb_type) * CHAR_BIT;
                            const double_limb_type product =
                                static_cast<double_limb_type>(x) * static_cast<double_limb_type>(y);
                            const double_limb_type sum0 =
                                static_cast<double_limb_type>(acc0) + static_cast<limb_type>(product);
                            acc0 = static_cast<limb_type>(sum0);

                            const double_limb_type sum1 = static_cast<double_limb_type>(acc1) +
                                                          static_cast<limb_type>(product >> limb_bits) +
                                                          (sum0 >> limb_bits);
                            acc1 = static_cast<limb_type>(sum1);
                            acc2 += static_cast<limb_type>(sum1 >> limb_bits);
                        }

                        template<typename Backend>
                        inline constexpr void multiply_4x4_emit(Backend &result, std::size_t idx, limb_type &acc0,
                                                                limb_type &acc1, limb_type &acc2) noexcept {
                            result.limbs()[idx] = acc0;
                            acc0 = acc1;
                            acc1 = acc2;
                            acc2 = 0u;
                        }

                        template<typename Backend>
                        inline constexpr void multiply_4x4(Backend &result, const limb_type *x,
                                                           const limb_type *y) noexcept {
                            result.zero_after(0u);

                            limb_type acc0 = 0u;
                            limb_type acc1 = 0u;
                            limb_type acc2 = 0u;

                            multiply_4x4_add(acc0, acc1, acc2, x[0], y[0]);
                            multiply_4x4_emit(result, 0u, acc0, acc1, acc2);

                            multiply_4x4_add(acc0, acc1, acc2, x[0], y[1]);
                            multiply_4x4_add(acc0, acc1, acc2, x[1], y[0]);
                            multiply_4x4_emit(result, 1u, acc0, acc1, acc2);

                            multiply_4x4_add(acc0, acc1, acc2, x[0], y[2]);
                            multiply_4x4_add(acc0, acc1, acc2, x[1], y[1]);
                            multiply_4x4_add(acc0, acc1, acc2, x[2], y[0]);
                            multiply_4x4_emit(result, 2u, acc0, acc1, acc2);

                            multiply_4x4_add(acc0, acc1, acc2, x[0], y[3]);
                            multiply_4x4_add(acc0, acc1, acc2, x[1], y[2]);
                            multiply_4x4_add(acc0, acc1, acc2, x[2], y[1]);
                            multiply_4x4_add(acc0, acc1, acc2, x[3], y[0]);
                            multiply_4x4_emit(result, 3u, acc0, acc1, acc2);

                            multiply_4x4_add(acc0, acc1, acc2, x[1], y[3]);
                            multiply_4x4_add(acc0, acc1, acc2, x[2], y[2]);
                            multiply_4x4_add(acc0, acc1, acc2, x[3], y[1]);
                            multiply_4x4_emit(result, 4u, acc0, acc1, acc2);

                            multiply_4x4_add(acc0, acc1, acc2, x[2], y[3]);
                            multiply_4x4_add(acc0, acc1, acc2, x[3], y[2]);
                            multiply_4x4_emit(result, 5u, acc0, acc1, acc2);

                            multiply_4x4_add(acc0, acc1, acc2, x[3], y[3]);
                            multiply_4x4_emit(result, 6u, acc0, acc1, acc2);
                            result.limbs()[7] = acc0;
                            result.limbs()[8] = acc1;
                            result.zero_after(9u);
                            result.set_carry(false);
                            result.normalize();
                        }

                        template<std::size_t InputLimbs, typename Backend>
                        inline constexpr void multiply_low_limbs(Backend &result, const limb_type *x,
                                                                 const limb_type *y) noexcept {
                            result.zero_after(0u);
                            result.set_carry(false);

                            for (std::size_t i = 0; i < InputLimbs; ++i) {
                                limb_type carry = 0u;
                                for (std::size_t j = 0; j < InputLimbs && i + j < result.size(); ++j) {
                                    const double_limb_type product =
                                        static_cast<double_limb_type>(x[i]) * static_cast<double_limb_type>(y[j]) +
                                        static_cast<double_limb_type>(result.limbs()[i + j]) + carry;
                                    result.limbs()[i + j] = static_cast<limb_type>(product);
                                    carry = static_cast<limb_type>(product >> Backend::limb_bits);
                                }

                                for (std::size_t idx = i + InputLimbs; carry != 0u && idx < result.size(); ++idx) {
                                    const double_limb_type sum =
                                        static_cast<double_limb_type>(result.limbs()[idx]) + carry;
                                    result.limbs()[idx] = static_cast<limb_type>(sum);
                                    carry = static_cast<limb_type>(sum >> Backend::limb_bits);
                                }
                            }

                            result.normalize();
                        }

                        template<typename Backend>
                        inline constexpr void multiply_by_limb(Backend &result, limb_type value) noexcept {
                            limb_type carry = 0u;
                            for (std::size_t i = 0; i < result.size(); ++i) {
                                const double_limb_type product = static_cast<double_limb_type>(result.limbs()[i]) *
                                                                     static_cast<double_limb_type>(value) +
                                                                 carry;
                                result.limbs()[i] = static_cast<limb_type>(product);
                                carry = static_cast<limb_type>(product >> Backend::limb_bits);
                            }
                            result.set_carry(false);
                            result.normalize();
                        }

                        inline constexpr bool limbs_ge_modulus_4(const limb_type *x, const limb_type *p) noexcept {
                            for (std::size_t step = 4u; step > 0u; --step) {
                                const std::size_t idx = step - 1u;
                                if (x[idx] < p[idx]) {
                                    return false;
                                }
                                if (x[idx] > p[idx]) {
                                    return true;
                                }
                            }
                            return true;
                        }

                        inline constexpr bool redc_result_ge_modulus_4(const limb_type *x,
                                                                       const limb_type *p) noexcept {
                            return x[4] != 0u || limbs_ge_modulus_4(x, p);
                        }

                        inline constexpr void subtract_modulus_4(limb_type *x, const limb_type *p) noexcept {
                            limb_type borrow = 0;
                            for (std::size_t i = 0; i < 4u; ++i) {
                                const limb_type subtrahend = p[i] + borrow;
                                const bool subtrahend_carry = subtrahend < p[i];
                                const limb_type current = x[i];
                                x[i] = current - subtrahend;
                                borrow = (subtrahend_carry || current < subtrahend) ? 1u : 0u;
                            }
                            x[4] -= borrow;
                        }

                        template<typename Backend, typename ModObj>
                        inline constexpr bool try_montgomery_reduce_4(Backend &result, const ModObj &mod_obj) {
                            if (result.size() > 9u) {
                                return false;
                            }

                            limb_type t[9] = {};
                            for (std::size_t i = 0; i < result.size(); ++i) {
                                t[i] = result.limbs()[i];
                            }

                            const limb_type *p = mod_obj.get_mod().limbs();
                            const limb_type p_dash = mod_obj.get_p_dash();
                            for (std::size_t i = 0; i < 4u; ++i) {
                                const limb_type m = t[i] * p_dash;
                                limb_type carry = 0;

                                for (std::size_t j = 0; j < 4u; ++j) {
                                    const double_limb_type product =
                                        static_cast<double_limb_type>(m) * static_cast<double_limb_type>(p[j]) +
                                        static_cast<double_limb_type>(t[i + j]) + carry;
                                    t[i + j] = static_cast<limb_type>(product);
                                    carry = static_cast<limb_type>(product >> Backend::limb_bits);
                                }

                                for (std::size_t idx = i + 4u; carry != 0u && idx < 9u; ++idx) {
                                    const double_limb_type sum = static_cast<double_limb_type>(t[idx]) + carry;
                                    t[idx] = static_cast<limb_type>(sum);
                                    carry = static_cast<limb_type>(sum >> Backend::limb_bits);
                                }
                            }

                            for (std::size_t i = 0u; i < 16u && redc_result_ge_modulus_4(t + 4u, p); ++i) {
                                subtract_modulus_4(t + 4u, p);
                            }

                            if (redc_result_ge_modulus_4(t + 4u, p)) {
                                Backend normalized;
                                for (std::size_t i = 0; i < 5u; ++i) {
                                    normalized.limbs()[i] = t[4u + i];
                                }
                                normalized.zero_after(5u);
                                normalized.set_carry(false);
                                normalized.normalize();

                                Backend large_mod = mod_obj.get_mod();
                                boost::multiprecision::backends::eval_modulus(normalized, large_mod);
                                result = normalized;
                                return true;
                            }

                            for (std::size_t i = 0; i < 4u; ++i) {
                                result.limbs()[i] = t[4u + i];
                            }
                            result.zero_after(4u);
                            result.set_carry(false);
                            result.normalize();
                            return true;
                        }
                    }    // namespace alt_bn128_fp12_limb_ops
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil
