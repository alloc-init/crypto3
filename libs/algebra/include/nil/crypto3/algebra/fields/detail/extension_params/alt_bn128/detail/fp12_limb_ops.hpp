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

                        template<size_t N>
                        inline void add_limbs_portable(limb *z, const limb *x, const limb *y) {
                            limb carry = 0u;
                            for (size_t i = 0; i < N; i++) {
                                const auto sum = (wide_limb)x[i] + y[i] + carry;
                                z[i] = (limb)sum;
                                carry = (limb)(sum >> limb_bits);
                            }
                        }

                        template<size_t N>
                        inline bool subtract_limbs_portable(limb *z, const limb *x, const limb *y) {
                            bool borrow = false;
                            for (size_t i = 0; i < N; i++) {
                                const limb subtrahend = y[i] + (limb)borrow;
                                const bool subtrahend_carry = subtrahend < y[i];
                                const limb current = x[i];
                                z[i] = current - subtrahend;
                                borrow = subtrahend_carry || current < subtrahend;
                            }
                            return borrow;
                        }

                        inline bool ge_modulus_4(const limb *x, const limb *mod) {
                            for (int i = 3; i >= 0; i--) {
                                if (x[i] < mod[i]) {
                                    return false;
                                }
                                if (x[i] > mod[i]) {
                                    return true;
                                }
                            }
                            return true;
                        }

                        inline bool ge_modulus(const limb *x, const limb *mod) {
                            if (x[4] != 0u) {
                                return true;
                            }
                            return ge_modulus_4(x, mod);
                        }

                        template<class Field>
                        inline void add_low_4_limbs_mod_portable(limb *z, const limb *x, const limb *y) {
                            add_limbs_portable<4>(z, x, y);
                            static const limb_array p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
                            // do one pass of normalization on lower limbs
                            if (ge_modulus(z, p.data())) {
                                subtract_limbs_portable<5>(z, z, p.data());
                            }
                        }

                        template<class Field>
                        inline void add_8_limbs_mod(limb_array &z, const limb_array &x, const limb_array &y) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            add_8_limbs_mod_x86<Field>(z, x, y);
#else
                            add_limbs_portable<8>(z.data(), x.data(), y.data());
                            static const limb_array p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
                            if (ge_modulus_4(z.data() + 4, p.data())) {
                                subtract_limbs_portable<4>(z.data() + 4, z.data() + 4, p.data());
                            }
#endif
                        }

                        template<class Field>
                        inline void subtract_8_limbs_mod(limb_array &z, const limb_array &x, const limb_array &y) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            subtract_8_limbs_mod_x86<Field>(z, x, y);
#else
                            bool borrow = subtract_limbs_portable<8>(z.data(), x.data(), y.data());
                            if (borrow) {
                                // if we went negative, add p
                                static const limb_array p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
                                add_limbs_portable<4>(z.data() + 4, z.data() + 4, p.data());
                            }
#endif
                        }

                        template<class Field>
                        inline void mul_8_limbs_by_9(limb_array &dst, const limb_array &src) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            mul_8_limbs_by_9_x86<Field>(dst, src);
#else
                            limb_array cpy = src;
                            add_8_limbs_mod<Field>(dst, src, src);    // 2x
                            add_8_limbs_mod<Field>(dst, dst, dst);    // 4x
                            add_8_limbs_mod<Field>(dst, dst, dst);    // 8x
                            add_8_limbs_mod<Field>(dst, dst, cpy);    // 9x
#endif
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
                        inline void multiply_emit(limb *result, size_t idx, limb &acc0, limb &acc1, limb &acc2) {
                            result[idx] = acc0;
                            acc0 = acc1;
                            acc1 = acc2;
                            acc2 = 0u;
                        }

                        // Multiply two base-width values using their low four limbs.
                        //
                        // The result is the full 8-limb product placed in the 9-limb storage shape. This is the
                        // common path for products of ordinary BN254 Fp Montgomery residues.
                        inline void multiply_4x4_portable(limb *result, const limb *x, const limb *y) {
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
                        }

                        inline void multiply_4x4(limb *z, const limb *x, const limb *y) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            multiply_4x4_x86(z, x, y);
#else
                            multiply_4x4_portable(z, x, y);
#endif
                        }

                        template<class Field>
                        inline void montgomery_reduce_portable(limb *result, const limb_array &data) {
                            // p is the field modulus as 4 limbs
                            static const limb_array p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
                            // p_dash is -p^{-1} modulo one limb, B = 2^64.
                            // Multiplying the current low limb by p_dash gives the
                            // one-limb factor m that makes t[i] + m * p[0] == 0 mod B.
                            limb p_dash = Field::modulus_params.get_mod_obj().get_p_dash();
                            limb_array buf = data;

                            // REDC over R = 2^(64 * 4). At step i, choose m so adding
                            // m * p shifted by i limbs makes buf[i] zero modulo 2^64.
                            // After four steps the low four limbs have been cancelled,
                            // so the high four limbs contain buf * R^-1 modulo p.
                            for (size_t i = 0; i < base_value_limb_count; i++) {
                                // Only the low limb of this product is used. Because
                                // p[0] * p_dash == -1 mod B, this m cancels buf[i] when
                                // m * p is added into the current REDC column.
                                const limb m = buf[i] * p_dash;
                                limb carry = 0;

                                // Add m * p into buf starting at limb i. The low limb of
                                // this sum is constructed to cancel buf[i].
                                for (size_t j = 0; j < base_value_limb_count; ++j) {
                                    const wide_limb product =
                                        (wide_limb)m * (wide_limb)p[j] + (wide_limb)buf[i + j] + carry;
                                    buf[i + j] = (limb)product;
                                    carry = (limb)(product >> limb_bits);
                                }

                                // Propagate any carry beyond the four modulus limbs.
                                for (size_t j = i + base_value_limb_count; carry != 0 && j < storage_limb_count; j++) {
                                    const wide_limb sum = (wide_limb)buf[j] + carry;
                                    buf[j] = (limb)sum;
                                    carry = (limb)(sum >> limb_bits);
                                }
                            }

                            // The REDC output lives in buf[4..7]. Bring it back into the
                            // canonical field range before moving the low four limbs.
                            if (ge_modulus_4(buf.data() + 4, p.data())) {
                                subtract_limbs_portable<4>(buf.data() + 4, buf.data() + 4, p.data());
                            }

                            // Keep the reduced 4-limb field value and clear the lazy
                            // extension limbs in the shared storage shape.
                            for (size_t i = 0; i < base_value_limb_count; i++) {
                                result[i] = buf[base_value_limb_count + i];
                            }
                        }

                        template<class Field>
                        inline void montgomery_reduce(limb *result, const limb_array &data) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            montgomery_reduce_x86<Field>(result, data.data());
#else
                            montgomery_reduce_portable<Field>(result, data);
#endif
                        }

                        // fp2_base ops must support pointers becuase fp2_view doesnt own its data
                        // output z is assumed continuous
                        template<class Field>
                        inline void fp2_base_add_mod(limb *z, const limb *const *x, const limb *const *y) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            fp2_base_add_mod_x86<Field>(z, x, y);
#else
                            add_low_4_limbs_mod_portable<Field>(z, x[0], y[0]);
                            add_low_4_limbs_mod_portable<Field>(z + 8, x[1], y[1]);
#endif
                        }

                        // fp2_base ops must support pointers becuase fp2_view doesnt own its data
                        inline void fp2_base_add_pre(limb *z, const limb *const *x, const limb *const *y) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            fp2_base_add_pre_x86(z, x, y);
#else
                            add_limbs_portable<4>(z, x[0], y[0]);
                            add_limbs_portable<4>(z + 8, x[1], y[1]);
#endif
                        }

                        template<class Field>
                        inline void fp2_sub_pre(limb_array *data, const limb_array *other) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            fp2_sub_pre_x86<Field>(data, other);
#else
                            subtract_8_limbs_mod<Field>(data[0], data[0], other[0]);
                            subtract_limbs_portable<8>(data[1].data(), data[1].data(), other[1].data());
#endif
                        }

                        // caller must handle aliasing
                        template<class Field>
                        inline void fp2_mul_xi_add(limb_array *dst, const limb_array *src, const limb_array *addend) {
                            mul_8_limbs_by_9<Field>(dst[0], src[0]);
                            subtract_8_limbs_mod<Field>(dst[0], dst[0], src[1]);
                            mul_8_limbs_by_9<Field>(dst[1], src[1]);
                            add_8_limbs_mod<Field>(dst[1], dst[1], src[0]);
                            add_8_limbs_mod<Field>(dst[0], dst[0], addend[0]);
                            add_8_limbs_mod<Field>(dst[1], dst[1], addend[1]);
                        }

                        template<class Field>
                        inline void fp2_mul_xi_add_modify_src(limb_array *src, const limb_array *addend) {
                            limb_array src0 = src[0];
                            mul_8_limbs_by_9<Field>(src[0], src0);
                            subtract_8_limbs_mod<Field>(src[0], src[0], src[1]);
                            mul_8_limbs_by_9<Field>(src[1], src[1]);
                            add_8_limbs_mod<Field>(src[1], src[1], src0);
                            add_8_limbs_mod<Field>(src[0], src[0], addend[0]);
                            add_8_limbs_mod<Field>(src[1], src[1], addend[1]);
                        }

                        template<class Field>
                        inline void fp2_mul_xi_add_modify_addend(limb_array *addend, const limb_array *src) {
                            limb_array tmp;
                            mul_8_limbs_by_9<Field>(tmp, src[0]);
                            subtract_8_limbs_mod<Field>(tmp, tmp, src[1]);
                            add_8_limbs_mod<Field>(addend[0], addend[0], tmp);
                            mul_8_limbs_by_9<Field>(tmp, src[1]);
                            add_8_limbs_mod<Field>(tmp, tmp, src[0]);
                            add_8_limbs_mod<Field>(addend[1], addend[1], tmp);
                        }

                        template<class Field>
                        inline void fp2_mul_pre_portable(limb_array *z, const limb *const *x, const limb *const *y) {
                            // For x = a + bu and y = c + du:
                            //   xy = (a + bu) * (c + du)
                            //      = ac + adu + bcu + bdu^2
                            //      = ac + (ad + bc)u - bd      # since u^2 = -1
                            //      = (ac - bd) + (ad + bc)u
                            // Karatsuba computes the cross term with one product:
                            //   ad + bc = (a + b)(c + d) - ac - bd.
                            limb_array ac, bd, a_plus_b, c_plus_d;
                            multiply_4x4(ac.data(), x[0], y[0]);
                            multiply_4x4(bd.data(), x[1], y[1]);
                            add_limbs_portable<4>(a_plus_b.data(), x[0], x[1]);
                            add_limbs_portable<4>(c_plus_d.data(), y[0], y[1]);
                            subtract_8_limbs_mod<Field>(z[0], ac, bd);
                            multiply_4x4(z[1].data(), a_plus_b.data(), c_plus_d.data());
                            subtract_limbs_portable<8>(z[1].data(), z[1].data(), ac.data());
                            subtract_limbs_portable<8>(z[1].data(), z[1].data(), bd.data());
                        }

                        template<class Field>
                        inline void fp2_mul_pre(limb_array *z, const limb *const *x, const limb *const *y) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                            fp2_mul_pre_x86<Field>(z, x, y);
#else
                            fp2_mul_pre_portable<Field>(z, x, y);
#endif
                        }
                    }    // namespace alt_bn128_fp12_limb_ops
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil
