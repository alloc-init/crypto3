#pragma once

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_types.hpp>

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_ops_x86.hpp>
#endif

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {
                namespace detail {
                    namespace alt_bn128_fp12_limb_ops {
                        // Loads limbs from a multiprecision backend value
                        template<typename Backend>
                        static limb_array load_limbs(const Backend &backend) {
                            static_assert(Backend::limb_bits == limb_bits,
                                          "alt_bn128 fp12 fast path expects 64-bit field limbs");
                            static_assert(Backend::internal_limb_count == base_value_limb_count,
                                          "alt_bn128 fp12 fast path expects 4 64-bit limbs");
                            limb_array result = {};
                            for (size_t i = 0; i < base_value_limb_count; i++) {
                                result[i] = (limb)backend.limbs()[i];
                            }
                            return result;
                        }

                        inline bool is_zero(const limb_array &x) {
                            for (size_t i = 0; i < x.size(); i++) {
                                if (x[i] != 0u) {
                                    return false;
                                }
                            }
                            return true;
                        }

                        template<size_t N>
                        inline void add_limbs_portable(limb *result, const limb *other) {
                            limb carry = 0u;
                            for (size_t i = 0; i < N; i++) {
                                const auto sum = (wide_limb)result[i] + other[i] + carry;
                                result[i] = (limb)sum;
                                carry = (limb)(sum >> limb_bits);
                            }
                        }

                        template<size_t N>
                        inline void add_limbs(limb *result, const limb *other) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            if constexpr (N == 4)
                                add_4_limbs_x86(result, other);
                            else if constexpr (N == 9)
                                add_9_limbs_x86(result, other);
                            else
                                add_limbs_portable<N>(result, other);
#else
                            add_limbs_portable<N>(result, other);
#endif
                        }

                        template<size_t N>
                        inline void add_limbs(limb_array &result, const limb_array &other) {
                            add_limbs<N>(result.data(), other.data());
                        }

                        template<size_t N>
                        inline bool subtract_limbs_portable(limb *result, const limb *other) {
                            bool borrow = false;
                            for (size_t i = 0; i < N; i++) {
                                const limb subtrahend = other[i] + (limb)borrow;
                                const bool subtrahend_carry = subtrahend < other[i];
                                const limb current = result[i];
                                result[i] = current - subtrahend;
                                borrow = subtrahend_carry || current < subtrahend;
                            }
                            return borrow;
                        }

                        template<size_t N>
                        inline bool subtract_limbs(limb *result, const limb *other) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            if constexpr (N == 9)
                                return subtract_9_limbs_x86(result, other);
                            else
                                return subtract_limbs_portable<N>(result, other);
#else
                            return subtract_limbs_portable<N>(result, other);
#endif
                        }

                        inline bool ge_modulus(const limb *x, const limb *y) {
                            if (x[4] != 0u) {
                                // p has 4 limbs, so if x has a nonzero 5th digit, it is greater
                                return true;
                            } else {
                                for (int i = 4; i >= 0; i--) {
                                    if (x[i] < y[i]) {
                                        return false;
                                    }
                                    if (x[i] > y[i]) {
                                        return true;
                                    }
                                }
                                return true;
                            }
                        }

                        // do one pass of normalization on lower limbs
                        template<class Field>
                        inline void subtract_modulus_lower(limb_array &data) {
                            static const limb_array p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
                            if (ge_modulus(data.data(), p.data())) {
                                subtract_limbs<5>(data.data(), p.data());
                            }
                        }

                        // loop until upper limbs are normalized
                        template<class Field>
                        inline void subtract_modulus_upper(limb_array &data) {
                            static const limb_array p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
                            if (ge_modulus(data.data() + 4, p.data())) {
                                subtract_limbs<5>(data.data() + 4, p.data());
                            }
                        }

                        template<class Field>
                        inline void add_4_limbs_mod(limb_array &data, const limb_array &other) {
                            add_limbs<4>(data, other);
                            subtract_modulus_lower<Field>(data);
                        }

                        template<class Field>
                        inline void add_8_limbs_mod(limb_array &data, const limb_array &other) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            add_8_limbs_mod_x86<Field>(data, other);
#else
                            add_limbs<8>(data, other);
                            subtract_modulus_upper<Field>(data);
#endif
                        }

                        template<class Field>
                        inline void subtract_limbs_mod(limb_array &data, const limb_array &other) {
                            bool borrow = subtract_limbs<8>(data.data(), other.data());
                            if (borrow) {
                                // if we went negative, add p
                                static const limb_array p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
                                add_limbs<4>(data.data() + base_value_limb_count, p.data());
                            }
                        }

                        inline void left_shift_one(limb_array &result) {
                            limb carry = 0u;
                            for (size_t i = 0; i < result.size(); i++) {
                                const limb next_carry = result[i] >> (limb_bits - 1u);
                                result[i] = (result[i] << 1u) | carry;
                                carry = next_carry;
                            }
                        }

                        inline void multiply_by_limb(limb_array &result, limb value) {
                            limb carry = 0u;
                            for (size_t i = 0; i < result.size(); i++) {
                                const wide_limb product = (wide_limb)result[i] * (wide_limb)value + carry;
                                result[i] = (limb)product;
                                carry = (limb)(product >> limb_bits);
                            }
                        }

                        // Add one limb product into the current Comba column accumulator.
                        //
                        // acc0 holds the limb being emitted for the current output column, acc1 holds the next
                        // carry limb, and acc2 collects overflow from acc1. Each x*y product is 128 bits, so adding
                        // its low half to acc0 and high half to acc1 lets a column accumulate several partial
                        // products before multiply_emit advances to the next column.
                        inline void multiply_partial(limb &acc0, limb &acc1, limb &acc2, limb x, limb y) {
                            // compute the base product x * y
                            const wide_limb product = (wide_limb)x * y;
                            // add the low bits of the new product to the lowest accumulator
                            const wide_limb sum0 = (wide_limb)acc0 + (limb)product;
                            // add the high bits of the new product to the mid accumulator, along with any carry from
                            // the low accumulator
                            const wide_limb sum1 = (wide_limb)acc1 + (limb)(product >> limb_bits) + (sum0 >> limb_bits);
                            // set the output accumulators to the low bits of the sums
                            acc0 = (limb)sum0;
                            acc1 = (limb)sum1;
                            // finally, add any overflow from the middle accumulator to the high accumulator
                            acc2 += (limb)(sum1 >> limb_bits);
                        }

                        // Emit the completed low limb for one output column and shift the accumulator forward.
                        //
                        // After this, acc0/acc1 contain the carry state for the next column and acc2 is clear for
                        // new overflow. This is shared by the fixed 4x4 and 5x5 kernels.
                        inline void multiply_emit(limb_array &result, size_t idx, limb &acc0, limb &acc1, limb &acc2) {
                            result[idx] = acc0;
                            acc0 = acc1;
                            acc1 = acc2;
                            acc2 = 0u;
                        }

                        // Multiply two base-width values using their low four limbs.
                        //
                        // The result is the full 8-limb product placed in the 9-limb storage shape. This is the
                        // common path for products of ordinary BN254 Fp Montgomery residues.
                        inline void
                            multiply_4x4_portable(limb_array &result, const limb_array &x, const limb_array &y) {
                            result = {};
                            limb acc0 = 0u;
                            limb acc1 = 0u;
                            limb acc2 = 0u;

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

                        inline void multiply_4x4(limb_array &result, const limb_array &x, const limb_array &y) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            multiply_4x4_x86(result, x, y);
#else
                            multiply_4x4_portable(result, x, y);
#endif
                        }

                        template<class Field>
                        inline void montgomery_reduce_portable(limb_array &data) {
                            // p is the field modulus as 4 limbs
                            static limb_array p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
                            // p_dash is -p^{-1} modulo one limb, B = 2^64.
                            // Multiplying the current low limb by p_dash gives the
                            // one-limb factor m that makes t[i] + m * p[0] == 0 mod B.
                            limb p_dash = Field::modulus_params.get_mod_obj().get_p_dash();

                            // REDC over R = 2^(64 * 4). At step i, choose m so adding
                            // m * p shifted by i limbs makes data[i] zero modulo 2^64.
                            // After four steps the low four limbs have been cancelled,
                            // so the high four limbs contain data * R^-1 modulo p.
                            for (size_t i = 0; i < base_value_limb_count; i++) {
                                // Only the low limb of this product is used. Because
                                // p[0] * p_dash == -1 mod B, this m cancels data[i] when
                                // m * p is added into the current REDC column.
                                const limb m = data[i] * p_dash;
                                limb carry = 0;

                                // Add m * p into data starting at limb i. The low limb of
                                // this sum is constructed to cancel data[i].
                                for (size_t j = 0; j < base_value_limb_count; ++j) {
                                    const wide_limb product =
                                        (wide_limb)m * (wide_limb)p[j] + (wide_limb)data[i + j] + carry;
                                    data[i + j] = (limb)product;
                                    carry = (limb)(product >> limb_bits);
                                }

                                // Propagate any carry beyond the four modulus limbs.
                                for (size_t j = i + base_value_limb_count; carry != 0 && j < data.size(); j++) {
                                    const wide_limb sum = (wide_limb)data[j] + carry;
                                    data[j] = (limb)sum;
                                    carry = (limb)(sum >> limb_bits);
                                }
                            }

                            // The REDC output lives in data[4..8]. Bring it back into the
                            // canonical field range before moving the low four limbs.
                            subtract_modulus_upper<Field>(data);

                            // Keep the reduced 4-limb field value and clear the lazy
                            // extension limbs in the shared storage shape.
                            for (size_t i = 0; i < base_value_limb_count; i++) {
                                data[i] = data[base_value_limb_count + i];
                            }
                            for (size_t i = base_value_limb_count; i < data.size(); i++) {
                                data[i] = 0;
                            }
                        }

                        template<class Field>
                        inline void montgomery_reduce(limb_array &data) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            montgomery_reduce_x86<Field>(data);
#else
                            montgomery_reduce_portable<Field>(data);
#endif
                        }
                    }    // namespace alt_bn128_fp12_limb_ops
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil
