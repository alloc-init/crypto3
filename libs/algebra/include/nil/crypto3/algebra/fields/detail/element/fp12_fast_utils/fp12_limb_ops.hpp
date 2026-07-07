#pragma once

#include <nil/crypto3/algebra/fields/detail/element/fp12_fast_utils/fp12_limb_types.hpp>

namespace nil::crypto3::algebra::fields::detail::fp12_fast_utils {
    template<typename Backend>
    static std::array<limb, Backend::internal_limb_count> load_limbs(const Backend &backend) {
        static_assert(Backend::limb_bits == limb_bits, "fp12 fast path expects 64-bit field limbs");
        std::array<limb, Backend::internal_limb_count> result = {};
        for (size_t i = 0; i < Backend::internal_limb_count; i++) {
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
        return fp12_fast_utils::ge_modulus<N>(x, mod);
    }

    template<class Field, size_t N>
    inline void add_low_limbs_mod_portable(limb *z, const limb *x, const limb *y) {
        limb tmp[N + 1] = {};
        limb carry = 0u;
        for (size_t i = 0; i < N; i++) {
            const auto sum = (wide_limb)x[i] + y[i] + carry;
            tmp[i] = (limb)sum;
            carry = (limb)(sum >> limb_bits);
        }
        tmp[N] = carry;
        static const auto p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
        // Normalize the 5-limb scratch, then copy only this coefficient back.
        // z may point into the middle of a contiguous fp2_base value.
        if (ge_modulus_wide<4>(tmp, p.data())) {
            subtract_limbs_portable<5>(tmp, tmp, p.data());
        }
        for (size_t i = 0; i < N; i++) {
            z[i] = tmp[i];
        }
    }

    template<size_t N>
    void multiply_portable(limb *result, const limb *x, const limb *y) {
        for (size_t i = 0; i < 2 * N; i++) {
            result[i] = 0;
        }
        for (size_t i = 0; i < N; i++) {
            limb carry = 0;
            for (size_t j = 0; j < N; j++) {
                wide_limb product = (wide_limb)x[j] * y[i];
                wide_limb sum = (wide_limb)result[i + j] + product + carry;
                result[i + j] = (limb)sum;
                carry = (limb)(sum >> limb_bits);
            }
            result[i + N] += carry;
        }
    }

    template<class Field, size_t N>
    inline void montgomery_reduce_portable(limb *result, const limb *data) {
        // p is the field modulus as 4 limbs
        static const auto p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
        // p_dash is -p^{-1} modulo one limb, B = 2^64.
        // Multiplying the current low limb by p_dash gives the
        // one-limb factor m that makes t[i] + m * p[0] == 0 mod B.
        limb p_dash = Field::modulus_params.get_mod_obj().get_p_dash();
        std::array<limb, N> buf;
        std::copy_n(data, N, buf.begin());
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

    template<class Field, size_t N>
    inline void add_limbs_mod(limb *z, const limb *x, const limb *y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        add_8_limbs_mod_x86<Field>(z, x, y);
#else
        add_limbs_portable<N>(z, x, y);
        static const auto p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
        if (ge_modulus<N / 2>(z + N / 2, p.data())) {
            subtract_limbs_portable<N / 2>(z + N / 2, z + N / 2, p.data());
        }
#endif
    }

    template<class Field, size_t N>
    inline void subtract_limbs_mod(limb *z, const limb *x, const limb *y) {
#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__)
        subtract_8_limbs_mod_x86<Field>(z, x, y);
#else
        bool borrow = subtract_limbs_portable<N>(z, x, y);
        if (borrow) {
            // If the full 8-limb subtraction borrowed, add p to the high half
            // to keep the double value in the p * R residue class.
            static const auto p = load_limbs(Field::modulus_params.get_mod_obj().get_mod());
            add_limbs_portable<N / 2>(z + N / 2, z + N / 2, p.data());
        }
#endif
    }
}    // namespace nil::crypto3::algebra::fields::detail::fp12_fast_utils