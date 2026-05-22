#pragma once

#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {
                namespace detail {
                    namespace alt_bn128_fp12_limb_ops {

                        using limb_type = std::uint64_t;
#if defined(__SIZEOF_INT128__)
                        using double_limb_type = unsigned __int128;
#else
#error "alt_bn128 fp12 limb ops require unsigned __int128 support"
#endif

                        const size_t limb_bits = sizeof(limb_type) * CHAR_BIT;
                        const size_t base_value_limb_count = 4u;
                        const size_t storage_limb_count = 9u;

                        using limb_array = std::array<limb_type, storage_limb_count>;

                        // Loads limbs from a multiprecision backend value
                        template<typename Backend>
                        static limb_array load_limbs(const Backend &backend) {
                            static_assert(Backend::limb_bits == limb_bits,
                                          "alt_bn128 fp12 fast path expects 64-bit field limbs");
                            static_assert(backend.size() == base_value_limb_count,
                                          "alt_bn128 fp12 fast path expects 4 64-bit limbs");
                            limb_array result = {};
                            for (size_t i = 0; i < base_value_limb_count; i++) {
                                result[i] = static_cast<limb_type>(backend.limbs()[i]);
                            }
                            return result;
                        }

                        bool is_zero(const limb_array &x) noexcept {
                            for (size_t i = 0; i < x.size(); i++) {
                                if (x[i] != 0u) {
                                    return false;
                                }
                            }
                            return true;
                        }

                        int compare_limbs(const limb_array &x, const limb_array &y) noexcept {
                            for (int i = x.size() - 1; i >= 0; i--) {
                                if (x[i] < y[i]) {
                                    return -1;
                                }
                                if (x[i] > y[i]) {
                                    return 1;
                                }
                            }
                            return 0;
                        }

                        void add_limbs(limb_array &result, const limb_array &other) noexcept {
                            limb_type carry = 0u;
                            for (size_t i = 0; i < result.size(); i++) {
                                const auto sum = (double_limb_type)result[i] + other[i] + carry;
                                result[i] = (limb_type)sum;
                                carry = (limb_type)(sum >> limb_bits);
                            }
                        }

                        void subtract_limbs(limb_array &result, const limb_array &other) noexcept {
                            limb_type borrow = 0u;
                            for (size_t i = 0; i < result.size(); i++) {
                                const limb_type subtrahend = other[i] + borrow;
                                const bool subtrahend_carry = subtrahend < other[i];
                                const limb_type current = result[i];
                                result[i] = current - subtrahend;
                                borrow = (subtrahend_carry || current < subtrahend) ? 1u : 0u;
                            }
                        }

                        void left_shift_one(limb_array &result) noexcept {
                            limb_type carry = 0u;
                            for (size_t i = 0; i < result.size(); i++) {
                                const limb_type next_carry = result[i] >> (limb_bits - 1u);
                                result[i] = (result[i] << 1u) | carry;
                                carry = next_carry;
                            }
                        }

                        // Add one limb product into the current Comba column accumulator.
                        //
                        // acc0 holds the limb being emitted for the current output column, acc1 holds the next
                        // carry limb, and acc2 collects overflow from acc1. Each x*y product is 128 bits, so adding
                        // its low half to acc0 and high half to acc1 lets a column accumulate several partial
                        // products before multiply_emit advances to the next column.
                        void multiply_partial(limb_type &acc0, limb_type &acc1, limb_type &acc2, limb_type x,
                                              limb_type y) noexcept {
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

                        // Emit the completed low limb for one output column and shift the accumulator forward.
                        //
                        // After this, acc0/acc1 contain the carry state for the next column and acc2 is clear for
                        // new overflow. This is shared by the fixed 4x4 and 5x5 kernels.
                        void multiply_emit(limb_array &result, size_t idx, limb_type &acc0, limb_type &acc1,
                                           limb_type &acc2) noexcept {
                            result[idx] = acc0;
                            acc0 = acc1;
                            acc1 = acc2;
                            acc2 = 0u;
                        }

                        // Multiply two base-width values using their low four limbs.
                        //
                        // The result is the full 8-limb product placed in the 9-limb storage shape. This is the
                        // common path for products of ordinary BN254 Fp Montgomery residues.
                        void multiply_4x4(limb_array &result, const limb_array &x, const limb_array &y) noexcept {
                            result = {};
                            limb_type acc0 = 0u;
                            limb_type acc1 = 0u;
                            limb_type acc2 = 0u;

                            multiply_partial(acc0, acc1, acc2, x[0], y[0]);
                            multiply_emit(result, 0u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[0], y[1]);
                            multiply_partial(acc0, acc1, acc2, x[1], y[0]);
                            multiply_emit(result, 1u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[0], y[2]);
                            multiply_partial(acc0, acc1, acc2, x[1], y[1]);
                            multiply_partial(acc0, acc1, acc2, x[2], y[0]);
                            multiply_emit(result, 2u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[0], y[3]);
                            multiply_partial(acc0, acc1, acc2, x[1], y[2]);
                            multiply_partial(acc0, acc1, acc2, x[2], y[1]);
                            multiply_partial(acc0, acc1, acc2, x[3], y[0]);
                            multiply_emit(result, 3u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[1], y[3]);
                            multiply_partial(acc0, acc1, acc2, x[2], y[2]);
                            multiply_partial(acc0, acc1, acc2, x[3], y[1]);
                            multiply_emit(result, 4u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[2], y[3]);
                            multiply_partial(acc0, acc1, acc2, x[3], y[2]);
                            multiply_emit(result, 5u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[3], y[3]);
                            multiply_emit(result, 6u, acc0, acc1, acc2);
                            result[7] = acc0;
                            result[8] = acc1;
                        }

                        // Multiply two Fp2-sum-width values using their low five limbs.
                        //
                        // Fp2 Karatsuba sums such as (a + b) can carry once past the four-limb base field value,
                        // so the cross-term product needs a 5x5 kernel. These inputs are bounded by the tower
                        // formulas, and their product fits the nine-limb pre-REDC storage used by this fast path.
                        void multiply_5x5(limb_array &result, const limb_array &x, const limb_array &y) noexcept {
                            result = {};
                            limb_type acc0 = 0u;
                            limb_type acc1 = 0u;
                            limb_type acc2 = 0u;

                            multiply_partial(acc0, acc1, acc2, x[0], y[0]);
                            multiply_emit(result, 0u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[0], y[1]);
                            multiply_partial(acc0, acc1, acc2, x[1], y[0]);
                            multiply_emit(result, 1u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[0], y[2]);
                            multiply_partial(acc0, acc1, acc2, x[1], y[1]);
                            multiply_partial(acc0, acc1, acc2, x[2], y[0]);
                            multiply_emit(result, 2u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[0], y[3]);
                            multiply_partial(acc0, acc1, acc2, x[1], y[2]);
                            multiply_partial(acc0, acc1, acc2, x[2], y[1]);
                            multiply_partial(acc0, acc1, acc2, x[3], y[0]);
                            multiply_emit(result, 3u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[0], y[4]);
                            multiply_partial(acc0, acc1, acc2, x[1], y[3]);
                            multiply_partial(acc0, acc1, acc2, x[2], y[2]);
                            multiply_partial(acc0, acc1, acc2, x[3], y[1]);
                            multiply_partial(acc0, acc1, acc2, x[4], y[0]);
                            multiply_emit(result, 4u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[1], y[4]);
                            multiply_partial(acc0, acc1, acc2, x[2], y[3]);
                            multiply_partial(acc0, acc1, acc2, x[3], y[2]);
                            multiply_partial(acc0, acc1, acc2, x[4], y[1]);
                            multiply_emit(result, 5u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[2], y[4]);
                            multiply_partial(acc0, acc1, acc2, x[3], y[3]);
                            multiply_partial(acc0, acc1, acc2, x[4], y[2]);
                            multiply_emit(result, 6u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[3], y[4]);
                            multiply_partial(acc0, acc1, acc2, x[4], y[3]);
                            multiply_emit(result, 7u, acc0, acc1, acc2);

                            multiply_partial(acc0, acc1, acc2, x[4], y[4]);
                            multiply_emit(result, 8u, acc0, acc1, acc2);
                        }

                        void multiply_by_limb(limb_array &result, limb_type value) noexcept {
                            limb_type carry = 0u;
                            for (size_t i = 0; i < result.size(); i++) {
                                const double_limb_type product =
                                    static_cast<double_limb_type>(result[i]) * static_cast<double_limb_type>(value) +
                                    carry;
                                result[i] = static_cast<limb_type>(product);
                                carry = static_cast<limb_type>(product >> limb_bits);
                            }
                        }

                        bool limbs_ge_modulus_4(const limb_type *x, const limb_type *p) noexcept {
                            for (size_t step = 4u; step > 0u; --step) {
                                const size_t idx = step - 1u;
                                if (x[idx] < p[idx]) {
                                    return false;
                                }
                                if (x[idx] > p[idx]) {
                                    return true;
                                }
                            }
                            return true;
                        }

                        bool redc_result_ge_modulus_4(const limb_type *x, const limb_type *p) noexcept {
                            return x[4] != 0u || limbs_ge_modulus_4(x, p);
                        }

                        void subtract_modulus_4(limb_type *x, const limb_type *p) noexcept {
                            limb_type borrow = 0;
                            for (size_t i = 0; i < 4u; i++) {
                                const limb_type subtrahend = p[i] + borrow;
                                const bool subtrahend_carry = subtrahend < p[i];
                                const limb_type current = x[i];
                                x[i] = current - subtrahend;
                                borrow = (subtrahend_carry || current < subtrahend) ? 1u : 0u;
                            }
                            x[4] -= borrow;
                        }

                        template<class Field>
                        void montgomery_reduce_4(limb_array &data) noexcept {
                            static limb_array p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
                            limb_type p_dash = static_cast<limb_type>(Field::modulus_params.get_mod_obj().get_p_dash());

                            limb_array t = data;
                            for (size_t i = 0; i < base_value_limb_count; i++) {
                                const limb_type m = t[i] * p_dash;
                                limb_type carry = 0;

                                for (size_t j = 0; j < base_value_limb_count; ++j) {
                                    const double_limb_type product =
                                        static_cast<double_limb_type>(m) * static_cast<double_limb_type>(p[j]) +
                                        static_cast<double_limb_type>(t[i + j]) + carry;
                                    t[i + j] = static_cast<limb_type>(product);
                                    carry = static_cast<limb_type>(product >> limb_bits);
                                }

                                for (size_t idx = i + base_value_limb_count; carry != 0u && idx < t.size(); idx++) {
                                    const double_limb_type sum = static_cast<double_limb_type>(t[idx]) + carry;
                                    t[idx] = static_cast<limb_type>(sum);
                                    carry = static_cast<limb_type>(sum >> limb_bits);
                                }
                            }

                            for (size_t i = 0u;
                                 i < 16u && redc_result_ge_modulus_4(t.data() + base_value_limb_count, p.data());
                                 i++) {
                                subtract_modulus_4(t.data() + base_value_limb_count, p.data());
                            }

                            for (size_t i = 0; i < base_value_limb_count; i++) {
                                data[i] = t[base_value_limb_count + i];
                            }
                            for (size_t i = base_value_limb_count; i < data.size(); i++) {
                                data[i] = 0;
                            }
                        }
                    }    // namespace alt_bn128_fp12_limb_ops
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil
