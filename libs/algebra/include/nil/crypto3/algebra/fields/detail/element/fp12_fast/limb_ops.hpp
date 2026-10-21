#pragma once

#include <algorithm>

#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/types.hpp>

#if (defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)) || defined(CRYPTO3_HAS_AVX512IFMA)
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

    // Sum canonical Montgomery residues, retaining the multiprecision implementation for other limb counts
    // and targets. The four-limb kernel changes only product accumulation; reduction still uses the existing
    // padded backend and canonicalization, so no dot-length-dependent narrowing is introduced.
    template<Fp12FastParams Params, typename Backend, typename LeftIterator, typename RightIterator>
    inline void montgomery_inner_product(Backend &result, LeftIterator first, LeftIterator last, RightIterator right) {
        const auto &modulus = Params::base_field_type::modulus_params.get_mod_obj();
#if defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)
        if constexpr (Params::base_value_limb_count == 4) {
            using policy = typename std::remove_cvref_t<decltype(modulus)>::policy_type;
            using accumulator_type = typename policy::Backend_doubled_padded_limbs;
            static_assert(Backend::internal_limb_count == 4 && Backend::limb_bits == 64);
            static_assert(sizeof(std::size_t) <= sizeof(limb));
            accumulator_type accumulator;
            for (; first != last; ++first, ++right) {
                BOOST_ASSERT(eval_lt(*first, modulus.get_mod()) && eval_lt(*right, modulus.get_mod()));
                multiply_add_4x4_x86(accumulator.limbs(), (*first).limbs(), (*right).limbs());
            }
            // With B=2^64 and R=B^4, at most B-1 canonical products sum to less than (B-1)*p^2.
            // REDC adds less than p*R; because p<R, both accumulation and reduction fit in nine limbs.
            modulus.montgomery_reduce(accumulator);
            const accumulator_type wide_modulus(modulus.get_mod());
            while (!eval_lt(accumulator, wide_modulus)) {
                eval_subtract(accumulator, wide_modulus);
            }
            result = accumulator;
            return;
        }
#endif
        modulus.montgomery_inner_product(result, first, last, right);
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
#if defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)
        if constexpr (requires { montgomery_reduce_x86<Params>(result, data.data()); }) {
            return montgomery_reduce_x86<Params>(result, data.data());
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
#if defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)
        if constexpr (requires { add_limbs_mod_x86<Params>(z.data(), x.data(), y.data()); }) {
            return add_limbs_mod_x86<Params>(z.data(), x.data(), y.data());
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
#if defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)
        if constexpr (requires { subtract_limbs_mod_x86<Params>(z.data(), x.data(), y.data()); }) {
            return subtract_limbs_mod_x86<Params>(z.data(), x.data(), y.data());
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
#if defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)
        return mul_12_limbs_by_5_unreduced_x86(dst.data(), src.data());
#endif
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
#if defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)
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
            const bool borrow = subtract_limbs_portable<N>(tmp, tmp, p.data());
            tmp[N] -= (limb)borrow;
        }
        for (size_t i = 0; i < N; i++) {
            z[i] = tmp[i];
        }
    }

    template<Fp12FastParams Params>
    inline void subtract_low_limbs_mod_portable(limb *z, const limb *x, const limb *y) {
        constexpr size_t N = Params::base_value_limb_count;
        const bool borrow = subtract_limbs_portable<N>(z, x, y);
        if (borrow) {
            static const auto p = modulus_limbs<Params>();
            add_limbs_portable<N>(z, z, p.data());
        }
    }

    template<Fp12FastParams Params>
    inline void mul_low_limbs_by_5_mod_portable(limb *z, const limb *x) {
        constexpr size_t N = Params::base_value_limb_count;
        std::array<limb, N> doubled, quadrupled;
        add_low_limbs_mod_portable<Params>(doubled.data(), x, x);
        add_low_limbs_mod_portable<Params>(quadrupled.data(), doubled.data(), doubled.data());
        add_low_limbs_mod_portable<Params>(z, quadrupled.data(), x);
    }

    template<Fp12FastParams Params>
    inline void mul_low_limbs_by_9_mod_portable(limb *z, const limb *x) {
        constexpr size_t N = Params::base_value_limb_count;
        std::array<limb, N> doubled, quadrupled, octupled;
        add_low_limbs_mod_portable<Params>(doubled.data(), x, x);
        add_low_limbs_mod_portable<Params>(quadrupled.data(), doubled.data(), doubled.data());
        add_low_limbs_mod_portable<Params>(octupled.data(), quadrupled.data(), quadrupled.data());
        add_low_limbs_mod_portable<Params>(z, octupled.data(), x);
    }

    // fp2_base values are two contiguous 4-limb coefficients; each output
    // coefficient is normalized modulo p.
    template<Fp12FastParams Params>
    inline void fp2_base_add_mod(typename Params::limb_array &z, const typename Params::limb_array &x,
                                 const typename Params::limb_array &y) {
#if defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)
        if constexpr (requires { fp2_base_add_mod_x86<Params>(z.data(), x.data(), y.data()); }) {
            return fp2_base_add_mod_x86<Params>(z.data(), x.data(), y.data());
        }
#endif
        constexpr size_t N = Params::base_value_limb_count;
        add_low_limbs_mod_portable<Params>(z.data(), x.data(), y.data());
        add_low_limbs_mod_portable<Params>(z.data() + N, x.data() + N, y.data() + N);
    }

    template<Fp12FastParams Params>
    inline void fp2_sub_pre(std::array<typename Params::limb_array, 2> &data,
                            const std::array<typename Params::limb_array, 2> &other) {
#if defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)
        if constexpr (requires { fp2_8_limbs_sub_pre_x86<Params>((limb *)&data, (limb *)&other); }) {
            return fp2_8_limbs_sub_pre_x86<Params>((limb *)&data, (limb *)&other);
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

    template<Fp12FastParams Params>
        requires(Params::u_squared == -1)
    inline void fp2_mul_pre_portable(std::array<typename Params::limb_array, 2> &z,
                                     const typename Params::limb_array &x, const typename Params::limb_array &y) {
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

    template<Fp12FastParams Params>
        requires(Params::u_squared == -5)
    inline void fp2_mul_pre_portable(std::array<typename Params::limb_array, 2> &z,
                                     const typename Params::limb_array &x, const typename Params::limb_array &y) {
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

    template<Fp12FastParams Params>
    inline void fp2_mul_pre(std::array<typename Params::limb_array, 2> &z, const typename Params::limb_array &x,
                            const typename Params::limb_array &y) {
#if defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)
        if constexpr (requires { fp2_mul_pre_x86<Params>((limb *)&z, x.data(), y.data()); }) {
            return fp2_mul_pre_x86<Params>((limb *)&z, x.data(), y.data());
        }
#endif
        return fp2_mul_pre_portable<Params>(z, x, y);
    }

    template<Fp12FastParams Params>
    inline void fp2_add_mul_pre(std::array<typename Params::limb_array, 2> &z, const typename Params::limb_array &a,
                                const typename Params::limb_array &b, const typename Params::limb_array &c,
                                const typename Params::limb_array &d) {
#if defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)
        if constexpr (requires { fp2_add_mul_pre_x86<Params>((limb *)&z, a.data(), b.data(), c.data(), d.data()); }) {
            return fp2_add_mul_pre_x86<Params>((limb *)&z, a.data(), b.data(), c.data(), d.data());
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

    // Five radix-2^52 digits, with Montgomery factor S=2^260 rather than R=2^256.
    // Restrict this representation to the four-limb fields whose modulus is below 2^254.
    template<Fp12FastParams Params>
        requires(Params::base_value_limb_count == 4 && Params::base_field_type::modulus_bits <= 254)
    struct montgomery_radix52 {
        static constexpr limb mask = (limb(1) << 52) - 1;
        static constexpr std::size_t digits = 5;
        static constexpr std::size_t max_terms = 256;
        using value_type = std::array<limb, digits>;

        static value_type pack(const limb *x) {
            return {x[0] & mask, ((x[0] >> 52) | (x[1] << 12)) & mask, ((x[1] >> 40) | (x[2] << 24)) & mask,
                    ((x[2] >> 28) | (x[3] << 36)) & mask, x[3] >> 16};
        }

        static const value_type &modulus() {
            static const value_type p = [] {
                const auto words = modulus_limbs<Params>();
                return pack(words.data());
            }();
            return p;
        }

        static limb p_dash() {
            return Params::base_field_type::modulus_params.get_mod_obj().get_p_dash() & mask;
        }

        template<typename Backend>
        static value_type from_montgomery64(const Backend &x) {
            auto words = load_limbs(x);
            // S/R=16: four canonical doublings convert the Montgomery scaling exactly.
            for (unsigned i = 0; i < 4; ++i) {
                add_low_limbs_mod_portable<Params>(words.data(), words.data(), words.data());
            }
            return pack(words.data());
        }

        template<typename Backend>
        static void to_montgomery64(Backend &result, value_type x) {
            // Choose m in [0,15] so x+m*p is divisible by 16. Since x<p, (x+m*p)/16<p.
            // This removes the extra four Montgomery bits without a general modular multiplication.
            const limb m = (x[0] * p_dash()) & 15;
            limb carry = 0;
            for (std::size_t i = 0; i < digits; ++i) {
                const wide_limb sum = wide_limb(m) * modulus()[i] + x[i] + carry;
                x[i] = limb(sum) & mask;
                carry = limb(sum >> 52);
            }
            BOOST_ASSERT(carry == 0);
            for (std::size_t i = 0; i < digits - 1; ++i) {
                x[i] = (x[i] >> 4) | ((x[i + 1] & 15) << 48);
            }
            x[4] >>= 4;
            result = Backend();
            result.limbs()[0] = x[0] | (x[1] << 52);
            result.limbs()[1] = (x[1] >> 12) | (x[2] << 40);
            result.limbs()[2] = (x[2] >> 24) | (x[3] << 28);
            result.limbs()[3] = (x[3] >> 36) | (x[4] << 16);
        }

        static void subtract_modulus(value_type &x) {
            value_type difference;
            limb borrow = 0;
            for (std::size_t i = 0; i < digits; ++i) {
                const limb subtrahend = modulus()[i] + borrow;
                difference[i] = (x[i] - subtrahend) & mask;
                borrow = x[i] < subtrahend;
            }
            if (!borrow) {
                x = difference;
            }
        }

        // Scalar reference: full 104-bit products in 128-bit columns, independently of IFMA's split products.
        // At most 256 terms contribute five products to a column, requiring fewer than 115 bits.
        template<typename LeftIterator, typename RightIterator>
        static value_type inner_product(LeftIterator first, LeftIterator last, RightIterator right) {
            value_type result {};
            while (first != last) {
                std::array<wide_limb, 10> t {};
                std::size_t count = 0;
                for (; first != last && count < max_terms; ++first, ++right, ++count) {
                    for (std::size_t i = 0; i < digits; ++i) {
                        for (std::size_t j = 0; j < digits; ++j) {
                            t[i + j] += wide_limb((*first)[i]) * (*right)[j];
                        }
                    }
                }
                for (std::size_t i = 0; i < digits; ++i) {
                    const limb m = (limb(t[i]) * p_dash()) & mask;
                    for (std::size_t j = 0; j < digits; ++j) {
                        t[i + j] += wide_limb(m) * modulus()[j];
                    }
                    t[i + 1] += t[i] >> 52;
                }
                value_type reduced;
                for (std::size_t i = 5; i < 9; ++i) {
                    reduced[i - 5] = limb(t[i]) & mask;
                    t[i + 1] += t[i] >> 52;
                }
                reduced[4] = limb(t[9]);
                // T<256*p^2 and p<2^254 imply REDC(T)<5*p<S. Four subtractions suffice.
                for (unsigned i = 0; i < 4; ++i) {
                    subtract_modulus(reduced);
                }
                limb carry = 0;
                for (std::size_t i = 0; i < digits; ++i) {
                    const limb sum = result[i] + reduced[i] + carry;
                    result[i] = sum & mask;
                    carry = sum >> 52;
                }
                BOOST_ASSERT(carry == 0);
                subtract_modulus(result);
            }
            return result;
        }
    };

    // Lane-wise field arithmetic on normalized radix-2^52 digits. add/subtract take canonical
    // S-Montgomery residues and return canonical residues with the same scaling; no REDC is needed.
    template<Fp12FastParams Params, typename Operations>
    struct montgomery_radix52_simd {
        using radix = montgomery_radix52<Params>;
        using value_type = std::array<typename Operations::value_type, radix::digits>;

        static value_type load(const limb *x) {
            value_type result;
            for (std::size_t i = 0; i < radix::digits; ++i) {
                result[i] = Operations::load(x + i * Operations::lanes);
            }
            return result;
        }

        static void store(limb *result, const value_type &x) {
            for (std::size_t i = 0; i < radix::digits; ++i) {
                Operations::store(result + i * Operations::lanes, x[i]);
            }
        }

        // Also used by dot-product REDC: x may be noncanonical, but must have normalized digits and x<S.
        static void subtract_modulus(value_type &x) {
            const auto zero = Operations::broadcast(0);
            const auto one = Operations::broadcast(1);
            const auto mask = Operations::broadcast(radix::mask);
            value_type difference;
            auto borrow = zero;
            for (std::size_t i = 0; i < radix::digits; ++i) {
                const auto subtrahend = Operations::add(Operations::broadcast(radix::modulus()[i]), borrow);
                difference[i] = Operations::bit_and(Operations::subtract(x[i], subtrahend), mask);
                borrow = Operations::select(Operations::less(x[i], subtrahend), one, zero);
            }
            const auto keep = Operations::less(zero, borrow);
            for (std::size_t i = 0; i < radix::digits; ++i) {
                x[i] = Operations::select(keep, x[i], difference[i]);
            }
        }

        static value_type add(const value_type &a, const value_type &b) {
            value_type result;
            const auto mask = Operations::broadcast(radix::mask);
            auto carry = Operations::broadcast(0);
            for (std::size_t i = 0; i < radix::digits; ++i) {
                const auto sum = Operations::add(Operations::add(a[i], b[i]), carry);
                result[i] = Operations::bit_and(sum, mask);
                carry = Operations::shift52(sum);
            }
            // a+b<2*p<2^255<S, so there is no carry beyond the fifth digit; one subtraction suffices.
            subtract_modulus(result);
            return result;
        }

        static value_type subtract(const value_type &a, const value_type &b) {
            value_type result;
            const auto zero = Operations::broadcast(0);
            const auto one = Operations::broadcast(1);
            const auto mask = Operations::broadcast(radix::mask);
            auto borrow = zero;
            for (std::size_t i = 0; i < radix::digits; ++i) {
                const auto subtrahend = Operations::add(b[i], borrow);
                result[i] = Operations::bit_and(Operations::subtract(a[i], subtrahend), mask);
                borrow = Operations::select(Operations::less(a[i], subtrahend), one, zero);
            }
            const auto negative = Operations::less(zero, borrow);
            auto carry = zero;
            for (std::size_t i = 0; i < radix::digits; ++i) {
                const auto correction = Operations::select(negative, Operations::broadcast(radix::modulus()[i]), zero);
                const auto sum = Operations::add(Operations::add(result[i], correction), carry);
                result[i] = Operations::bit_and(sum, mask);
                carry = Operations::shift52(sum);
            }
            // A borrowed difference is S+a-b. Discarding its carry after adding p gives a-b+p in [0,p).
            return result;
        }
    };

    // Left is count consecutive five-digit scalars. Right is [term][digit][lane]; indices selects
    // nonzero outer terms without repacking the cached panels. All digits are canonical S-Montgomery inputs.
    // Result has [digit][lane] layout. Padded coefficient lanes are zero and never read past a panel.
    template<Fp12FastParams Params, typename Operations>
    inline void montgomery_inner_product52(limb *result, const limb *left, const limb *right,
                                           const std::size_t *indices, std::size_t count) {
        using radix = montgomery_radix52<Params>;
        using vector = typename Operations::value_type;
        constexpr std::size_t lanes = Operations::lanes;
        const auto zero = Operations::broadcast(0);
        const auto mask = Operations::broadcast(radix::mask);
        const auto inverse = Operations::broadcast(radix::p_dash());
        std::array<vector, 5> modulus, total;
        for (std::size_t i = 0; i < 5; ++i) {
            modulus[i] = Operations::broadcast(radix::modulus()[i]);
            total[i] = zero;
        }
        for (std::size_t first = 0; first < count;) {
            const std::size_t terms = std::min(radix::max_terms, count - first);
            std::array<vector, 10> t;
            t.fill(zero);
            // Each column receives at most 9*256 split products, each <2^52: below 2^64.
            // Normalize between chunks, never permitting an unbounded IFMA accumulator.
            for (std::size_t term = first; term < first + terms; ++term) {
                const limb *y = right + indices[term] * 5 * lanes;
                std::array<vector, 5> b;
                for (std::size_t j = 0; j < 5; ++j) {
                    b[j] = Operations::load(y + j * lanes);
                }
                for (std::size_t i = 0; i < 5; ++i) {
                    const auto a = Operations::broadcast(left[5 * term + i]);
                    for (std::size_t j = 0; j < 5; ++j) {
                        t[i + j] = Operations::multiply_add_low(t[i + j], a, b[j]);
                        t[i + j + 1] = Operations::multiply_add_high(t[i + j + 1], a, b[j]);
                    }
                }
            }
            for (std::size_t i = 0; i < 9; ++i) {
                t[i + 1] = Operations::add(t[i + 1], Operations::shift52(t[i]));
                t[i] = Operations::bit_and(t[i], mask);
            }
            // Five cancellation steps divide by S=2^260. IFMA ignores upper operand bits, which is
            // exactly the reduction modulo 2^52 needed for m; the carry of t[i] goes to t[i+1].
            for (std::size_t i = 0; i < 5; ++i) {
                const auto m = Operations::multiply_add_low(zero, t[i], inverse);
                for (std::size_t j = 0; j < 5; ++j) {
                    t[i + j] = Operations::multiply_add_low(t[i + j], m, modulus[j]);
                    t[i + j + 1] = Operations::multiply_add_high(t[i + j + 1], m, modulus[j]);
                }
                t[i + 1] = Operations::add(t[i + 1], Operations::shift52(t[i]));
            }
            std::array<vector, 5> reduced;
            for (std::size_t i = 5; i < 9; ++i) {
                reduced[i - 5] = Operations::bit_and(t[i], mask);
                t[i + 1] = Operations::add(t[i + 1], Operations::shift52(t[i]));
            }
            reduced[4] = t[9];
            for (unsigned i = 0; i < 4; ++i) {
                montgomery_radix52_simd<Params, Operations>::subtract_modulus(reduced);
            }
            auto carry = zero;
            for (std::size_t i = 0; i < 5; ++i) {
                const auto sum = Operations::add(Operations::add(total[i], reduced[i]), carry);
                total[i] = Operations::bit_and(sum, mask);
                carry = Operations::shift52(sum);
            }
            montgomery_radix52_simd<Params, Operations>::subtract_modulus(total);
            first += terms;
        }
        for (std::size_t i = 0; i < 5; ++i) {
            Operations::store(result + i * lanes, total[i]);
        }
    }

}    // namespace nil::crypto3::algebra::fields::detail::fp12_fast
