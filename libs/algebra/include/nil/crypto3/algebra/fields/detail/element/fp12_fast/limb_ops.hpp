#pragma once

#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/types.hpp>

#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/limb_ops_x86.hpp>
#endif

namespace nil::crypto3::algebra::fields::detail::fp12_fast {
    template<typename Backend>
    static constexpr std::array<limb, Backend::internal_limb_count> load_limbs(const Backend &backend) {
        static_assert(Backend::limb_bits == limb_bits, "fp12 fast path expects 64-bit field limbs");
        std::array<limb, Backend::internal_limb_count> result = {};
        for (size_t i = 0; i < Backend::internal_limb_count; i++) {
            result[i] = (limb)backend.limbs()[i];
        }
        return result;
    }

    template<Fp12FastParams Params>
    auto modulus_limbs() {
        return load_limbs(Params::base_field_type::modulus_params.get_mod_obj().get_mod());
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

    template<size_t N>
    inline bool ge_modulus(const limb *x, const limb *mod) {
        for (int i = N - 1; i >= 0; i--) {
            if (x[i] < mod[i]) {
                return false;
            }
            if (x[i] > mod[i]) {
                return true;
            }
        }
        return true;
    }

    template<size_t N>
    inline bool ge_modulus_wide(const limb *x, const limb *mod) {
        if (x[N] != 0u) {
            return true;
        }
        return ge_modulus<N>(x, mod);
    }

    template<size_t BaseLimbCount>
    inline void multiply_portable(limb *z, const limb *x, const limb *y) {
        for (size_t i = 0; i < 2 * BaseLimbCount; i++) {
            z[i] = 0;
        }
        for (size_t i = 0; i < BaseLimbCount; i++) {
            limb carry = 0;
            for (size_t j = 0; j < BaseLimbCount; j++) {
                wide_limb product = (wide_limb)x[j] * y[i];
                wide_limb sum = (wide_limb)z[i + j] + product + carry;
                z[i + j] = (limb)sum;
                carry = (limb)(sum >> limb_bits);
            }
            z[i + BaseLimbCount] += carry;
        }
    }

    template<Fp12FastParams Params>
    inline void montgomery_reduce(limb *result, const typename Params::limb_array &data) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        if constexpr (Params::storage_limb_count == 8) {
            return montgomery_reduce_8_limbs_x86<typename Params::base_field_type>(result, data.data());
        } else if constexpr (Params::storage_limb_count == 12) {
            return montgomery_reduce_12_limbs_x86<typename Params::base_field_type>(result, data.data());
        }
#endif
        constexpr size_t N = Params::storage_limb_count;
        // p is the field modulus as 4 limbs
        static const auto p = modulus_limbs<Params>();
        // p_dash is -p^{-1} modulo one limb, B = 2^64.
        // Multiplying the current low limb by p_dash gives the
        // one-limb factor m that makes t[i] + m * p[0] == 0 mod B.
        limb p_dash = Params::base_field_type::modulus_params.get_mod_obj().get_p_dash();
        typename Params::limb_array buf = data;
        // REDC over R = 2^(64 * 4). At step i, choose m so adding
        // m * p shifted by i limbs makes buf[i] zero modulo 2^64.
        // After N/2 steps the low N/2 limbs have been cancelled,
        // so the high N/2 limbs contain buf * R^-1 modulo p.
        for (size_t i = 0; i < N / 2; i++) {
            // Only the low limb of this product is used. Because
            // p[0] * p_dash == -1 mod B, this m cancels buf[i] when
            // m * p is added into the current REDC column.
            const limb m = buf[i] * p_dash;
            limb carry = 0;
            // Add m * p into buf starting at limb i. The low limb of
            // this sum is constructed to cancel buf[i].
            for (size_t j = 0; j < N / 2; ++j) {
                const wide_limb product = (wide_limb)m * (wide_limb)p[j] + (wide_limb)buf[i + j] + carry;
                buf[i + j] = (limb)product;
                carry = (limb)(product >> limb_bits);
            }
            // Propagate any carry beyond the N/2 modulus limbs.
            for (size_t j = i + N / 2; carry != 0 && j < N; j++) {
                const wide_limb sum = (wide_limb)buf[j] + carry;
                buf[j] = (limb)sum;
                carry = (limb)(sum >> limb_bits);
            }
        }
        // The REDC output lives in buf[4..7]. Bring it back into the
        // canonical field range before writing the N/2 output limbs.
        if (ge_modulus<N / 2>(buf.data() + N / 2, p.data())) {
            subtract_limbs_portable<N / 2>(buf.data() + N / 2, buf.data() + N / 2, p.data());
        }
        // Write the reduced N/2-limb field value to the caller-provided Fp storage.
        for (size_t i = 0; i < N / 2; i++) {
            result[i] = buf[N / 2 + i];
        }
    }

    template<Fp12FastParams Params>
    inline void add_limbs_mod(typename Params::limb_array &z, const typename Params::limb_array &x,
                              const typename Params::limb_array &y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        if constexpr (Params::storage_limb_count == 8) {
            return add_8_limbs_mod_x86<typename Params::base_field_type>(z.data(), x.data(), y.data());
        } else if constexpr (Params::storage_limb_count == 12) {
            return add_12_limbs_mod_x86<typename Params::base_field_type>(z.data(), x.data(), y.data());
        }
#endif
        constexpr size_t N = Params::base_value_limb_count;
        add_limbs_portable<Params::storage_limb_count>(z.data(), x.data(), y.data());
        static const auto p = modulus_limbs<Params>();
        if (ge_modulus<N>(z.data() + N, p.data())) {
            subtract_limbs_portable<N>(z.data() + N, z.data() + N, p.data());
        }
    }

    template<Fp12FastParams Params>
    inline void subtract_limbs_mod(typename Params::limb_array &z, const typename Params::limb_array &x,
                                   const typename Params::limb_array &y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        if constexpr (Params::storage_limb_count == 8) {
            return subtract_8_limbs_mod_x86<typename Params::base_field_type>(z.data(), x.data(), y.data());
        } else if constexpr (Params::storage_limb_count == 12) {
            return subtract_12_limbs_mod_x86<typename Params::base_field_type>(z.data(), x.data(), y.data());
        }
#endif
        constexpr size_t N = Params::base_value_limb_count;
        bool borrow = subtract_limbs_portable<Params::storage_limb_count>(z.data(), x.data(), y.data());
        if (borrow) {
            // If the full 8-limb subtraction borrowed, add p to the high half
            // to keep the double value in the p * R residue class.
            static const auto p = modulus_limbs<Params>();
            add_limbs_portable<N>(z.data() + N, z.data() + N, p.data());
        }
    }

    template<Fp12FastParams Params>
    inline void mul_limbs_by_5_unreduced(typename Params::limb_array &dst, const typename Params::limb_array &src) {
        constexpr size_t N = Params::storage_limb_count;
        constexpr size_t M = Params::base_value_limb_count;
        typename Params::limb_array cpy = src;
        // shift by 2
        limb overflow = 0;
        for (size_t i = 0; i < N; i++) {
            dst[i] = src[i] << 2 | overflow;
            overflow = src[i] >> (CHAR_BIT * sizeof(limb) - 2);    // overflow = top 2 bits
        }
        // add 1
        add_limbs_portable<N>(dst.data(), dst.data(), cpy.data());
        // note: we dont have to reduce for BLS12_377, this will always fit within bounds < pR, see fp12 test
    }

    template<Fp12FastParams Params>
    inline void mul_limbs_by_9(typename Params::limb_array &dst, const typename Params::limb_array &src) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        if constexpr (Params::storage_limb_count == 8) {
            return mul_8_limbs_by_9_x86<typename Params::base_field_type>(dst.data(), src.data());
        }
#endif
        typename Params::limb_array cpy = src;
        add_limbs_mod<Params>(dst, src, src);    // 2x
        add_limbs_mod<Params>(dst, dst, dst);    // 4x
        add_limbs_mod<Params>(dst, dst, dst);    // 8x
        add_limbs_mod<Params>(dst, dst, cpy);    // 9x
    }

    template<Fp12FastParams Params>
    inline void add_low_limbs_mod_portable(limb *z, const limb *x, const limb *y) {
        constexpr size_t N = Params::base_value_limb_count;
        limb tmp[N + 1] = {};
        limb carry = 0u;
        for (size_t i = 0; i < N; i++) {
            const auto sum = (wide_limb)x[i] + y[i] + carry;
            tmp[i] = (limb)sum;
            carry = (limb)(sum >> limb_bits);
        }
        tmp[N] = carry;
        static const auto p = modulus_limbs<Params>();
        // Normalize the 5-limb scratch, then copy only this coefficient back.
        // z may point into the middle of a contiguous fp2_base value.
        if (ge_modulus_wide<N>(tmp, p.data())) {
            subtract_limbs_portable<N + 1>(tmp, tmp, p.data());
        }
        for (size_t i = 0; i < N; i++) {
            z[i] = tmp[i];
        }
    }

    // fp2_base values are two contiguous 4-limb coefficients; each output
    // coefficient is normalized modulo p.
    template<Fp12FastParams Params>
    inline void fp2_base_add_mod(typename Params::limb_array &z, const typename Params::limb_array &x,
                                 const typename Params::limb_array &y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        if constexpr (Params::base_value_limb_count == 4) {
            return fp2_base_4_limbs_add_mod_x86<typename Params::base_field_type>(z.data(), x.data(), y.data());
        } else if constexpr (Params::base_value_limb_count == 6) {
            return fp2_base_6_limbs_add_mod_x86<typename Params::base_field_type>(z.data(), x.data(), y.data());
        }
#endif
        constexpr size_t N = Params::base_value_limb_count;
        add_low_limbs_mod_portable<Params>(z.data(), x.data(), y.data());
        add_low_limbs_mod_portable<Params>(z.data() + N, x.data() + N, y.data() + N);
    }

    template<Fp12FastParams Params>
    inline void fp2_sub_pre(std::array<typename Params::limb_array, 2> &data,
                            const std::array<typename Params::limb_array, 2> &other) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        if constexpr (Params::storage_limb_count == 8) {
            return fp2_sub_pre_x86<typename Params::base_field_type>((limb *)&data, (limb *)&other);
        }
#endif
        constexpr size_t N = Params::storage_limb_count;
        subtract_limbs_mod<Params>(data[0], data[0], other[0]);
        subtract_limbs_portable<Params::storage_limb_count>(data[1].data(), data[1].data(), other[1].data());
    }

    // xi = 9 + u, u^2 = -1
    template<Fp12FastParams Params>
        requires(Params::xi == std::array {9, 1} && Params::u_squared == -1)
    inline void fp2_mul_xi_add(std::array<typename Params::limb_array, 2> &dst,
                               const std::array<typename Params::limb_array, 2> &src,
                               const std::array<typename Params::limb_array, 2> &addend) {
        typename Params::limb_array buf[2];    // necessary since inputs might alias dst
        mul_limbs_by_9<Params>(buf[0], src[0]);
        subtract_limbs_mod<Params>(buf[0], buf[0], src[1]);
        mul_limbs_by_9<Params>(buf[1], src[1]);
        add_limbs_mod<Params>(buf[1], buf[1], src[0]);
        add_limbs_mod<Params>(buf[0], buf[0], addend[0]);
        add_limbs_mod<Params>(buf[1], buf[1], addend[1]);
        dst[0] = buf[0];
        dst[1] = buf[1];
    }

    // xi = 1 + u, u^2 = -1
    template<Fp12FastParams Params>
        requires(Params::xi == std::array {1, 1} && Params::u_squared == -1)
    inline void fp2_mul_xi_add(std::array<typename Params::limb_array, 2> &dst,
                               const std::array<typename Params::limb_array, 2> &src,
                               const std::array<typename Params::limb_array, 2> &addend) {
        typename Params::limb_array buf[2];    // necessary since inputs might alias dst
        subtract_limbs_mod<Params>(buf[0], src[0], src[1]);
        add_limbs_mod<Params>(buf[1], src[0], src[1]);
        add_limbs_mod<Params>(buf[0], buf[0], addend[0]);
        add_limbs_mod<Params>(buf[1], buf[1], addend[1]);
        dst[0] = buf[0];
        dst[1] = buf[1];
    }

    // xi = u, u^2 = -5
    template<Fp12FastParams Params>
        requires(Params::xi == std::array {0, 1} && Params::u_squared == -5)
    inline void fp2_mul_xi_add(std::array<typename Params::limb_array, 2> &dst,
                               const std::array<typename Params::limb_array, 2> &src,
                               const std::array<typename Params::limb_array, 2> &addend) {
        // (a + bu)(0 + u)
        //   = 0 + au + 0 + bu^2
        //   = au + bu^2
        //   = -5b + au
        typename Params::limb_array buf[2];    // necessary since inputs might alias dst
        // src[1] is a nonnegative imaginary coefficient assembled from at most
        // three Fp2 products, so src[1] < 6p^2. For BLS12-377, 30p < R
        // (R / p is about 152), hence 5 * src[1] < p * R.
        mul_limbs_by_5_unreduced<Params>(buf[0], src[1]);
        subtract_limbs_mod<Params>(buf[0], {}, buf[0]);
        add_limbs_mod<Params>(buf[0], buf[0], addend[0]);
        add_limbs_mod<Params>(buf[1], src[0], addend[1]);
        dst[0] = buf[0];
        dst[1] = buf[1];
    }

    // fp2 mul pre for fields where u^2 = -1
    template<Fp12FastParams Params>
        requires(Params::u_squared == -1)
    inline void fp2_mul_pre(std::array<typename Params::limb_array, 2> &z, const typename Params::limb_array &x,
                            const typename Params::limb_array &y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        if constexpr (Params::storage_limb_count == 8) {
            return fp2_mul_pre_x86<typename Params::base_field_type>((limb *)&z, x.data(), y.data());
        }
#endif
        constexpr size_t N = Params::base_value_limb_count;
        constexpr size_t M = Params::storage_limb_count;
        // For x = a + bu and y = c + du:
        //   xy = (a + bu) * (c + du)
        //      = ac + adu + bcu + bdu^2
        //      = ac + (ad + bc)u - bd      # since u^2 = -1
        //      = (ac - bd) + (ad + bc)u
        // Karatsuba computes the cross term with one product:
        //   ad + bc = (a + b)(c + d) - ac - bd.
        typename Params::limb_array ac, bd, a_plus_b, c_plus_d;
        multiply_portable<N>(ac.data(), x.data(), y.data());
        multiply_portable<N>(bd.data(), x.data() + N, y.data() + N);
        add_limbs_portable<N>(a_plus_b.data(), x.data(), x.data() + N);
        add_limbs_portable<N>(c_plus_d.data(), y.data(), y.data() + N);
        subtract_limbs_mod<Params>(z[0], ac, bd);
        multiply_portable<N>(z[1].data(), a_plus_b.data(), c_plus_d.data());
        subtract_limbs_portable<M>(z[1].data(), z[1].data(), ac.data());
        subtract_limbs_portable<M>(z[1].data(), z[1].data(), bd.data());
    }

    // fp2 mul pre for fields where u^2 = -5
    template<Fp12FastParams Params>
        requires(Params::u_squared == -5)
    inline void fp2_mul_pre(std::array<typename Params::limb_array, 2> &z, const typename Params::limb_array &x,
                            const typename Params::limb_array &y) {
        constexpr size_t N = Params::base_value_limb_count;
        constexpr size_t M = Params::storage_limb_count;
        // For x = a + bu and y = c + du:
        //   xy = (a + bu) * (c + du)
        //      = ac + adu + bcu + bdu^2
        //      = ac + (ad + bc)u - 5bd      # since u^2 = -5
        //      = (ac - 5bd) + (ad + bc)u
        // Karatsuba computes the cross term with one product:
        //   ad + bc = (a + b)(c + d) - ac - bd.
        typename Params::limb_array ac, bd, bd5, a_plus_b, c_plus_d;
        multiply_portable<N>(ac.data(), x.data(), y.data());
        multiply_portable<N>(bd.data(), x.data() + N, y.data() + N);
        // x and y coefficients are each below 2p when this is reached through
        // fp2_add_mul_pre, so bd < 4p^2. For BLS12-377, 20p < R, and therefore
        // 5 * bd < p * R without normalization.
        mul_limbs_by_5_unreduced<Params>(bd5, bd);
        subtract_limbs_mod<Params>(z[0], ac, bd5);
        add_limbs_portable<N>(a_plus_b.data(), x.data(), x.data() + N);
        add_limbs_portable<N>(c_plus_d.data(), y.data(), y.data() + N);
        multiply_portable<N>(z[1].data(), a_plus_b.data(), c_plus_d.data());
        subtract_limbs_portable<M>(z[1].data(), z[1].data(), ac.data());
        subtract_limbs_portable<M>(z[1].data(), z[1].data(), bd.data());
    }

    // fp2 add mul pre for fields where u^2 = -1
    template<Fp12FastParams Params>
        requires(Params::u_squared == -1)
    inline void fp2_add_mul_pre(std::array<typename Params::limb_array, 2> &z, const typename Params::limb_array &a,
                                const typename Params::limb_array &b, const typename Params::limb_array &c,
                                const typename Params::limb_array &d) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        if constexpr (Params::storage_limb_count == 8) {
            return fp2_add_mul_pre_x86<typename Params::base_field_type>((limb *)&z, a.data(), b.data(), c.data(),
                                                                         d.data());
        }
#endif
        constexpr size_t N = Params::base_value_limb_count;
        // Build the raw fp2 sums in the same packed layout expected by fp2_mul_pre.
        typename Params::limb_array x, y;
        add_limbs_portable<N>(x.data(), a.data(), b.data());
        add_limbs_portable<N>(x.data() + N, a.data() + N, b.data() + N);
        add_limbs_portable<N>(y.data(), c.data(), d.data());
        add_limbs_portable<N>(y.data() + N, c.data() + N, d.data() + N);
        fp2_mul_pre<Params>(z, x, y);
    }

    // Generic fp2 add mul pre
    template<Fp12FastParams Params>
        requires(Params::u_squared != -1)
    inline void fp2_add_mul_pre(std::array<typename Params::limb_array, 2> &z, const typename Params::limb_array &a,
                                const typename Params::limb_array &b, const typename Params::limb_array &c,
                                const typename Params::limb_array &d) {
        constexpr size_t N = Params::base_value_limb_count;
        // Build the raw fp2 sums in the same packed layout expected by fp2_mul_pre.
        typename Params::limb_array x, y;
        add_limbs_portable<N>(x.data(), a.data(), b.data());
        add_limbs_portable<N>(x.data() + N, a.data() + N, b.data() + N);
        add_limbs_portable<N>(y.data(), c.data(), d.data());
        add_limbs_portable<N>(y.data() + N, c.data() + N, d.data() + N);
        fp2_mul_pre<Params>(z, x, y);
    }

}    // namespace nil::crypto3::algebra::fields::detail::fp12_fast
