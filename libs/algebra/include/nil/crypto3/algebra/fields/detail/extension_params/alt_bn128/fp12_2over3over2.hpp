//---------------------------------------------------------------------------//
// Copyright (c) 2020-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
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

#ifndef CRYPTO3_ALGEBRA_FIELDS_ALT_BN128_FP12_2OVER3OVER2_EXTENSION_PARAMS_HPP
#define CRYPTO3_ALGEBRA_FIELDS_ALT_BN128_FP12_2OVER3OVER2_EXTENSION_PARAMS_HPP

#include <array>
#include <climits>

#include <boost/multiprecision/cpp_int.hpp>

#include <nil/crypto3/algebra/fields/params.hpp>
#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/fp6_3over2.hpp>
#include <nil/crypto3/algebra/fields/fp2.hpp>


namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {
                template<typename BaseField>
                class fp12_2over3over2;

                namespace detail {
                    template<typename BaseField>
                    class fp12_2over3over2_extension_params;

                    template<typename BaseField>
                    class fp12_2over3over2;

                    /************************* ALT_BN128 ***********************************/

                    template<std::size_t Version>
                    class fp12_2over3over2_extension_params<fields::alt_bn128<Version>>
                            : public params<fields::alt_bn128<Version>> {
                        typedef fields::alt_bn128<Version> base_field_type;
                        typedef params<base_field_type> policy_type;

                    public:
                        using field_type = fields::fp12_2over3over2<base_field_type>;
                        typedef typename policy_type::integral_type integral_type;

                        constexpr static const integral_type modulus = policy_type::modulus;

                        typedef fields::fp2<base_field_type> non_residue_field_type;
                        typedef typename non_residue_field_type::value_type non_residue_type;
                        typedef fields::fp6_3over2<base_field_type> underlying_field_type;
                        typedef typename underlying_field_type::value_type underlying_type;
                        typedef typename base_field_type::value_type base_value_type;
                        typedef boost::multiprecision::limb_type limb_type;
                        typedef boost::multiprecision::double_limb_type double_limb_type;

                        constexpr static const std::size_t limb_bits = sizeof(limb_type) * CHAR_BIT;
                        constexpr static const std::size_t base_limb_count =
                            (base_field_type::modulus_bits + limb_bits - 1) / limb_bits;
                        constexpr static const std::size_t wide_limb_count =
                            (2 * base_field_type::modulus_bits + 32 + limb_bits - 1) / limb_bits;
                        constexpr static const limb_type max_limb = ~limb_type(0u);

                        typedef std::array<limb_type, base_limb_count> base_limb_array_type;
                        typedef std::array<limb_type, wide_limb_count> wide_limb_array_type;
                        typedef std::array<limb_type, wide_limb_count + 1> redc_limb_array_type;
                        constexpr static const std::size_t montgomery_limb_count =
                            wide_limb_count + 1 - base_limb_count;
                        typedef std::array<limb_type, montgomery_limb_count> montgomery_limb_array_type;
                        typedef std::array<limb_type, montgomery_limb_count + 1> montgomery_division_limb_array_type;

                        constexpr static const std::array<integral_type, 12 * 2> Frobenius_coeffs_c1 = {
                            0x01,
                            0x00,
                            0x1284B71C2865A7DFE8B99FDD76E68B605C521E08292F2176D60B35DADCC9E470_cppui_modular253,
                            0x246996F3B4FAE7E6A6327CFE12150B8E747992778EEEC7E5CA5CF05F80F362AC_cppui_modular254,
                            0x30644E72E131A0295E6DD9E7E0ACCCB0C28F069FBB966E3DE4BD44E5607CFD49_cppui_modular254,
                            0x00,
                            0x19DC81CFCC82E4BBEFE9608CD0ACAA90894CB38DBE55D24AE86F7D391ED4A67F_cppui_modular253,
                            0xABF8B60BE77D7306CBEEE33576139D7F03A5E397D439EC7694AA2BF4C0C101_cppui_modular248,
                            0x30644E72E131A0295E6DD9E7E0ACCCB0C28F069FBB966E3DE4BD44E5607CFD48_cppui_modular254,
                            0x00,
                            0x757CAB3A41D3CDC072FC0AF59C61F302CFA95859526B0D41264475E420AC20F_cppui_modular251,
                            0xCA6B035381E35B618E9B79BA4E2606CA20B7DFD71573C93E85845E34C4A5B9C_cppui_modular252,
                            0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD46_cppui_modular254,
                            0x00,
                            0x1DDF9756B8CBF849CF96A5D90A9ACCFD3B2F4C893F42A9166615563BFBB318D7_cppui_modular253,
                            0xBFAB77F2C36B843121DC8B86F6C4CCF2307D819D98302A771C39BB757899A9B_cppui_modular252,
                            0x59E26BCEA0D48BACD4F263F1ACDB5C4F5763473177FFFFFE_cppui_modular191,
                            0x00,
                            0x1687CCA314AEBB6DC866E529B0D4ADCD0E34B703AA1BF84253B10EDDB9A856C8_cppui_modular253,
                            0x2FB855BCD54A22B6B18456D34C0B44C0187DC4ADD09D90A0C58BE1EAE3BC3C46_cppui_modular254,
                            0x59E26BCEA0D48BACD4F263F1ACDB5C4F5763473177FFFFFF_cppui_modular191,
                            0x00,
                            0x290C83BF3D14634DB120850727BB392D6A86D50BD34B19B929BC44B896723B38_cppui_modular254,
                            0x23BD9E3DA9136A739F668E1ADC9EF7F0F575EC93F71A8DF953C846338C32A1AB_cppui_modular254
                        };

                        constexpr static const non_residue_type non_residue = non_residue_type(0x09, 0x01);

                        static bool is_zero(const wide_limb_array_type &x) {
                            for (std::size_t i = 0; i < wide_limb_count; ++i) {
                                if (x[i] != 0u) {
                                    return false;
                                }
                            }
                            return true;
                        }

                        static bool is_zero(const base_limb_array_type &x) {
                            for (std::size_t i = 0; i < base_limb_count; ++i) {
                                if (x[i] != 0u) {
                                    return false;
                                }
                            }
                            return true;
                        }

                        static int compare_abs(const wide_limb_array_type &x,
                                               const wide_limb_array_type &y) {
                            for (std::size_t i = wide_limb_count; i > 0; --i) {
                                const std::size_t idx = i - 1;
                                if (x[idx] < y[idx]) {
                                    return -1;
                                }
                                if (x[idx] > y[idx]) {
                                    return 1;
                                }
                            }
                            return 0;
                        }

                        static int compare_base(const base_limb_array_type &x,
                                                const base_limb_array_type &y) {
                            for (std::size_t i = base_limb_count; i > 0; --i) {
                                const std::size_t idx = i - 1;
                                if (x[idx] < y[idx]) {
                                    return -1;
                                }
                                if (x[idx] > y[idx]) {
                                    return 1;
                                }
                            }
                            return 0;
                        }

                        static wide_limb_array_type add_abs(const wide_limb_array_type &x,
                                                            const wide_limb_array_type &y) {
                            wide_limb_array_type result = {};
                            double_limb_type carry = 0;
                            for (std::size_t i = 0; i < wide_limb_count; ++i) {
                                const double_limb_type sum =
                                    static_cast<double_limb_type>(x[i]) +
                                    static_cast<double_limb_type>(y[i]) + carry;
                                result[i] = static_cast<limb_type>(sum);
                                carry = sum >> limb_bits;
                            }
                            return result;
                        }

                        static wide_limb_array_type sub_abs(const wide_limb_array_type &x,
                                                            const wide_limb_array_type &y) {
                            wide_limb_array_type result = {};
                            limb_type borrow = 0;
                            for (std::size_t i = 0; i < wide_limb_count; ++i) {
                                const limb_type yi = y[i] + borrow;
                                const bool carry_from_borrow = yi < y[i];
                                const bool needs_borrow = carry_from_borrow || x[i] < yi;
                                result[i] = x[i] - yi;
                                borrow = needs_borrow ? 1u : 0u;
                            }
                            return result;
                        }

                        static void sub_base_inplace(base_limb_array_type &x,
                                                     const base_limb_array_type &y) {
                            limb_type borrow = 0;
                            for (std::size_t i = 0; i < base_limb_count; ++i) {
                                const limb_type yi = y[i] + borrow;
                                const bool carry_from_borrow = yi < y[i];
                                const bool needs_borrow = carry_from_borrow || x[i] < yi;
                                x[i] -= yi;
                                borrow = needs_borrow ? 1u : 0u;
                            }
                        }

                        static base_limb_array_type sub_base(const base_limb_array_type &x,
                                                             const base_limb_array_type &y) {
                            base_limb_array_type result = x;
                            sub_base_inplace(result, y);
                            return result;
                        }

                        static base_limb_array_type add_base(const base_limb_array_type &x,
                                                             const base_limb_array_type &y) {
                            base_limb_array_type result = {};
                            double_limb_type carry = 0;
                            for (std::size_t i = 0; i < base_limb_count; ++i) {
                                const double_limb_type sum =
                                    static_cast<double_limb_type>(x[i]) +
                                    static_cast<double_limb_type>(y[i]) + carry;
                                result[i] = static_cast<limb_type>(sum);
                                carry = sum >> limb_bits;
                            }
                            return result;
                        }

                        static wide_limb_array_type mul_abs_small(const wide_limb_array_type &x,
                                                                  limb_type factor) {
                            wide_limb_array_type result = {};
                            double_limb_type carry = 0;
                            for (std::size_t i = 0; i < wide_limb_count; ++i) {
                                const double_limb_type product =
                                    static_cast<double_limb_type>(x[i]) *
                                    static_cast<double_limb_type>(factor) + carry;
                                result[i] = static_cast<limb_type>(product);
                                carry = product >> limb_bits;
                            }
                            return result;
                        }

                        static wide_limb_array_type mul_base_abs(const base_limb_array_type &x,
                                                                 const base_limb_array_type &y) {
                            wide_limb_array_type result = {};
                            for (std::size_t i = 0; i < base_limb_count; ++i) {
                                double_limb_type carry = 0;
                                for (std::size_t j = 0; j < base_limb_count; ++j) {
                                    const std::size_t idx = i + j;
                                    const double_limb_type product =
                                        static_cast<double_limb_type>(x[i]) *
                                        static_cast<double_limb_type>(y[j]) +
                                        static_cast<double_limb_type>(result[idx]) + carry;
                                    result[idx] = static_cast<limb_type>(product);
                                    carry = product >> limb_bits;
                                }

                                for (std::size_t idx = i + base_limb_count;
                                     carry != 0u && idx < wide_limb_count; ++idx) {
                                    const double_limb_type sum =
                                        static_cast<double_limb_type>(result[idx]) + carry;
                                    result[idx] = static_cast<limb_type>(sum);
                                    carry = sum >> limb_bits;
                                }
                            }
                            return result;
                        }

                        static base_limb_array_type backend_to_base_limbs(
                            const typename integral_type::backend_type &backend) {
                            base_limb_array_type result = {};
                            for (std::size_t i = 0; i < base_limb_count && i < backend.size(); ++i) {
                                result[i] = backend.limbs()[i];
                            }
                            return result;
                        }

                        static const base_limb_array_type &modulus_limbs() {
                            static const base_limb_array_type value = backend_to_base_limbs(modulus.backend());
                            return value;
                        }

                        static base_limb_array_type as_base_limbs(const base_value_type &x) {
                            return backend_to_base_limbs(x.data.backend().base_data());
                        }

                        static unsigned leading_zeroes(limb_type x) {
                            unsigned result = 0;
                            for (std::size_t i = limb_bits; i > 0; --i) {
                                if (((x >> (i - 1)) & 1u) != 0u) {
                                    break;
                                }
                                ++result;
                            }
                            return result;
                        }

                        static base_limb_array_type shift_left_base(const base_limb_array_type &x,
                                                                    unsigned shift) {
                            if (shift == 0u) {
                                return x;
                            }

                            base_limb_array_type result = {};
                            limb_type carry = 0;
                            for (std::size_t i = 0; i < base_limb_count; ++i) {
                                result[i] = (x[i] << shift) | carry;
                                carry = x[i] >> (limb_bits - shift);
                            }
                            return result;
                        }

                        template<typename LimbArray>
                        static void add_divisor_back(LimbArray &u,
                                                     const base_limb_array_type &v,
                                                     std::size_t offset) {
                            double_limb_type carry = 0;
                            for (std::size_t i = 0; i < base_limb_count; ++i) {
                                const double_limb_type sum =
                                    static_cast<double_limb_type>(u[offset + i]) +
                                    static_cast<double_limb_type>(v[i]) + carry;
                                u[offset + i] = static_cast<limb_type>(sum);
                                carry = sum >> limb_bits;
                            }

                            for (std::size_t idx = offset + base_limb_count;
                                 carry != 0u && idx < u.size(); ++idx) {
                                const double_limb_type sum =
                                    static_cast<double_limb_type>(u[idx]) + carry;
                                u[idx] = static_cast<limb_type>(sum);
                                carry = sum >> limb_bits;
                            }
                        }

                        template<typename LimbArray>
                        static bool sub_mul_divisor(LimbArray &u,
                                                    const base_limb_array_type &v,
                                                    limb_type qhat,
                                                    std::size_t offset) {
                            limb_type borrow = 0;
                            double_limb_type carry = 0;

                            for (std::size_t i = 0; i < base_limb_count; ++i) {
                                const double_limb_type product =
                                    static_cast<double_limb_type>(qhat) *
                                    static_cast<double_limb_type>(v[i]) + carry;
                                const limb_type product_low = static_cast<limb_type>(product);
                                carry = product >> limb_bits;

                                const limb_type subtrahend = product_low + borrow;
                                const bool subtrahend_carry = subtrahend < product_low;
                                const limb_type current = u[offset + i];
                                u[offset + i] = current - subtrahend;
                                borrow = (subtrahend_carry || current < subtrahend) ? 1u : 0u;
                            }

                            const double_limb_type final_subtrahend =
                                carry + static_cast<double_limb_type>(borrow);
                            const limb_type final_low = static_cast<limb_type>(final_subtrahend);
                            const bool final_carry = (final_subtrahend >> limb_bits) != 0u;
                            const limb_type current = u[offset + base_limb_count];
                            u[offset + base_limb_count] = current - final_low;
                            return final_carry || current < final_low;
                        }

                        template<typename LimbArray>
                        static base_limb_array_type shift_right_remainder(const LimbArray &u,
                                                                          unsigned shift) {
                            base_limb_array_type result = {};
                            if (shift == 0u) {
                                for (std::size_t i = 0; i < base_limb_count; ++i) {
                                    result[i] = u[i];
                                }
                                return result;
                            }

                            for (std::size_t i = 0; i < base_limb_count; ++i) {
                                const limb_type high = (i + 1 < u.size()) ? u[i + 1] : 0u;
                                result[i] = (u[i] >> shift) | (high << (limb_bits - shift));
                            }
                            return result;
                        }

                        static montgomery_division_limb_array_type shift_left_montgomery(
                            const montgomery_limb_array_type &x,
                            unsigned shift) {
                            montgomery_division_limb_array_type result = {};
                            if (shift == 0u) {
                                for (std::size_t i = 0; i < montgomery_limb_count; ++i) {
                                    result[i] = x[i];
                                }
                                return result;
                            }

                            limb_type carry = 0;
                            for (std::size_t i = 0; i < montgomery_limb_count; ++i) {
                                result[i] = (x[i] << shift) | carry;
                                carry = x[i] >> (limb_bits - shift);
                            }
                            result[montgomery_limb_count] = carry;
                            return result;
                        }

                        static base_limb_array_type reduce_montgomery_result_abs(
                            const montgomery_limb_array_type &x) {
                            const base_limb_array_type &p = modulus_limbs();
                            const unsigned shift = leading_zeroes(p[base_limb_count - 1]);
                            const base_limb_array_type v = shift_left_base(p, shift);
                            montgomery_division_limb_array_type u = shift_left_montgomery(x, shift);

                            constexpr std::size_t quotient_steps =
                                montgomery_limb_count - base_limb_count;
                            for (std::size_t step = quotient_steps + 1; step > 0; --step) {
                                const std::size_t j = step - 1;
                                const double_limb_type numerator =
                                    (static_cast<double_limb_type>(u[j + base_limb_count]) << limb_bits) |
                                    static_cast<double_limb_type>(u[j + base_limb_count - 1]);

                                double_limb_type qhat_wide =
                                    numerator / static_cast<double_limb_type>(v[base_limb_count - 1]);
                                limb_type qhat =
                                    qhat_wide > static_cast<double_limb_type>(max_limb) ?
                                        max_limb :
                                        static_cast<limb_type>(qhat_wide);
                                double_limb_type rhat =
                                    numerator -
                                    static_cast<double_limb_type>(qhat) *
                                        static_cast<double_limb_type>(v[base_limb_count - 1]);

                                while ((rhat >> limb_bits) == 0u &&
                                       static_cast<double_limb_type>(qhat) *
                                               static_cast<double_limb_type>(v[base_limb_count - 2]) >
                                           ((rhat << limb_bits) |
                                            static_cast<double_limb_type>(u[j + base_limb_count - 2]))) {
                                    --qhat;
                                    rhat += static_cast<double_limb_type>(v[base_limb_count - 1]);
                                }

                                if (sub_mul_divisor(u, v, qhat, j)) {
                                    add_divisor_back(u, v, j);
                                }
                            }

                            return shift_right_remainder(u, shift);
                        }

                        static limb_type montgomery_neg_inverse() {
                            static const limb_type value = []() {
                                limb_type inverse = 1u;
                                const limb_type p0 = modulus_limbs()[0];
                                for (std::size_t i = 0; i < limb_bits; ++i) {
                                    inverse *= limb_type(2u) - p0 * inverse;
                                }
                                return limb_type(0u) - inverse;
                            }();
                            return value;
                        }

                        static base_limb_array_type montgomery_reduce_abs(const wide_limb_array_type &x) {
                            redc_limb_array_type t = {};
                            for (std::size_t i = 0; i < wide_limb_count; ++i) {
                                t[i] = x[i];
                            }

                            const base_limb_array_type &p = modulus_limbs();
                            const limb_type p_dash = montgomery_neg_inverse();

                            for (std::size_t i = 0; i < base_limb_count; ++i) {
                                const limb_type m = t[i] * p_dash;
                                double_limb_type carry = 0;

                                for (std::size_t j = 0; j < base_limb_count; ++j) {
                                    const std::size_t idx = i + j;
                                    const double_limb_type product =
                                        static_cast<double_limb_type>(m) *
                                        static_cast<double_limb_type>(p[j]) +
                                        static_cast<double_limb_type>(t[idx]) + carry;
                                    t[idx] = static_cast<limb_type>(product);
                                    carry = product >> limb_bits;
                                }

                                for (std::size_t idx = i + base_limb_count;
                                     carry != 0u && idx < t.size(); ++idx) {
                                    const double_limb_type sum =
                                        static_cast<double_limb_type>(t[idx]) + carry;
                                    t[idx] = static_cast<limb_type>(sum);
                                    carry = sum >> limb_bits;
                                }
                            }

                            montgomery_limb_array_type shifted = {};
                            for (std::size_t i = 0; i < montgomery_limb_count; ++i) {
                                shifted[i] = t[i + base_limb_count];
                            }

                            return reduce_montgomery_result_abs(shifted);
                        }

                        static base_value_type make_montgomery_base_value(const base_limb_array_type &x) {
                            typename integral_type::backend_type backend;
                            for (std::size_t i = 0; i < base_limb_count && i < backend.size(); ++i) {
                                backend.limbs()[i] = x[i];
                            }
                            backend.zero_after(base_limb_count);
                            backend.normalize();

                            typename base_value_type::data_type data;
                            data.backend().base_data() = backend;
                            return base_value_type(data);
                        }

                        struct fp_dbl {
                            wide_limb_array_type data = {};
                            bool negative = false;

                            fp_dbl() = default;
                            explicit fp_dbl(const wide_limb_array_type &in_data,
                                            bool in_negative = false)
                                : data(in_data), negative(in_negative) {
                                normalize();
                            }

                            void normalize() {
                                if (is_zero(data)) {
                                    negative = false;
                                }
                            }

                            fp_dbl operator-() const {
                                fp_dbl result(*this);
                                if (!is_zero(result.data)) {
                                    result.negative = !result.negative;
                                }
                                return result;
                            }

                            fp_dbl operator+(const fp_dbl &other) const {
                                if (negative == other.negative) {
                                    return fp_dbl(add_abs(data, other.data), negative);
                                }

                                const int cmp = compare_abs(data, other.data);
                                if (cmp == 0) {
                                    return fp_dbl();
                                }
                                if (cmp > 0) {
                                    return fp_dbl(sub_abs(data, other.data), negative);
                                }
                                return fp_dbl(sub_abs(other.data, data), other.negative);
                            }

                            fp_dbl operator-(const fp_dbl &other) const {
                                return *this + (-other);
                            }

                            fp_dbl &operator+=(const fp_dbl &other) {
                                *this = *this + other;
                                return *this;
                            }

                            fp_dbl &operator-=(const fp_dbl &other) {
                                *this = *this - other;
                                return *this;
                            }

                            fp_dbl mul_small(limb_type factor) const {
                                if (factor == 0u || is_zero(data)) {
                                    return fp_dbl();
                                }
                                return fp_dbl(mul_abs_small(data, factor), negative);
                            }

                            static fp_dbl mul_pre(const base_limb_array_type &x,
                                                  const base_limb_array_type &y) {
                                return fp_dbl(mul_base_abs(x, y));
                            }

                            static fp_dbl mul_pre(const base_value_type &x,
                                                  const base_value_type &y) {
                                return mul_pre(as_base_limbs(x), as_base_limbs(y));
                            }

                            static fp_dbl sqr_pre(const base_value_type &x) {
                                const base_limb_array_type x_limbs = as_base_limbs(x);
                                return mul_pre(x_limbs, x_limbs);
                            }

                            static base_value_type reduce(const fp_dbl &x) {
                                base_limb_array_type reduced = montgomery_reduce_abs(x.data);
                                if (x.negative && !is_zero(reduced)) {
                                    reduced = sub_base(modulus_limbs(), reduced);
                                }
                                return make_montgomery_base_value(reduced);
                            }
                        };

                        struct fp2_dbl {
                            std::array<fp_dbl, 2> data;

                            fp2_dbl() = default;
                            fp2_dbl(const fp_dbl &c0, const fp_dbl &c1) : data({c0, c1}) {}

                            fp2_dbl operator+(const fp2_dbl &other) const {
                                return fp2_dbl(data[0] + other.data[0], data[1] + other.data[1]);
                            }

                            fp2_dbl operator-(const fp2_dbl &other) const {
                                return fp2_dbl(data[0] - other.data[0], data[1] - other.data[1]);
                            }

                            fp2_dbl &operator+=(const fp2_dbl &other) {
                                data[0] += other.data[0];
                                data[1] += other.data[1];
                                return *this;
                            }

                            fp2_dbl &operator-=(const fp2_dbl &other) {
                                data[0] -= other.data[0];
                                data[1] -= other.data[1];
                                return *this;
                            }

                            static fp2_dbl mul_pre(const non_residue_type &x,
                                                   const non_residue_type &y) {
                                const base_limb_array_type a = as_base_limbs(x.data[0]);
                                const base_limb_array_type b = as_base_limbs(x.data[1]);
                                const base_limb_array_type c = as_base_limbs(y.data[0]);
                                const base_limb_array_type d = as_base_limbs(y.data[1]);
                                const fp_dbl ac = fp_dbl::mul_pre(a, c);
                                const fp_dbl bd = fp_dbl::mul_pre(b, d);

                                return fp2_dbl(ac - bd,
                                               fp_dbl::mul_pre(add_base(a, b), add_base(c, d)) - ac - bd);
                            }

                            static fp2_dbl sqr_pre(const non_residue_type &x) {
                                const base_limb_array_type a = as_base_limbs(x.data[0]);
                                const base_limb_array_type b = as_base_limbs(x.data[1]);
                                const fp_dbl aa = fp_dbl::mul_pre(a, a);
                                const fp_dbl bb = fp_dbl::mul_pre(b, b);
                                return fp2_dbl(aa - bb, fp_dbl::mul_pre(a, b).mul_small(2u));
                            }

                            static fp2_dbl mul_xi(const fp2_dbl &x) {
                                return fp2_dbl(x.data[0].mul_small(9u) - x.data[1],
                                               x.data[0] + x.data[1].mul_small(9u));
                            }

                            static non_residue_type reduce(const fp2_dbl &x) {
                                return non_residue_type(fp_dbl::reduce(x.data[0]),
                                                        fp_dbl::reduce(x.data[1]));
                            }
                        };

                        struct fp6_dbl {
                            std::array<fp2_dbl, 3> data;

                            fp6_dbl() = default;
                            fp6_dbl(const fp2_dbl &c0, const fp2_dbl &c1,
                                    const fp2_dbl &c2)
                                : data({c0, c1, c2}) {}

                            fp6_dbl operator+(const fp6_dbl &other) const {
                                return fp6_dbl(data[0] + other.data[0],
                                               data[1] + other.data[1],
                                               data[2] + other.data[2]);
                            }

                            fp6_dbl operator-(const fp6_dbl &other) const {
                                return fp6_dbl(data[0] - other.data[0],
                                               data[1] - other.data[1],
                                               data[2] - other.data[2]);
                            }

                            fp6_dbl &operator+=(const fp6_dbl &other) {
                                data[0] += other.data[0];
                                data[1] += other.data[1];
                                data[2] += other.data[2];
                                return *this;
                            }

                            fp6_dbl &operator-=(const fp6_dbl &other) {
                                data[0] -= other.data[0];
                                data[1] -= other.data[1];
                                data[2] -= other.data[2];
                                return *this;
                            }

                            static fp6_dbl mul_pre(const underlying_type &x,
                                                   const underlying_type &y) {
                                const non_residue_type &a = x.data[0];
                                const non_residue_type &b = x.data[1];
                                const non_residue_type &c = x.data[2];
                                const non_residue_type &d = y.data[0];
                                const non_residue_type &e = y.data[1];
                                const non_residue_type &f = y.data[2];

                                fp2_dbl za = fp2_dbl::mul_pre(b + c, e + f);
                                fp2_dbl zb = fp2_dbl::mul_pre(a + b, e + d);
                                fp2_dbl zc = fp2_dbl::mul_pre(a + c, d + f);
                                const fp2_dbl be = fp2_dbl::mul_pre(b, e);
                                fp2_dbl cf = fp2_dbl::mul_pre(c, f);
                                const fp2_dbl ad = fp2_dbl::mul_pre(a, d);

                                za -= be;
                                za -= cf;
                                zb -= ad;
                                zb -= be;
                                zc -= ad;
                                zc -= cf;
                                za = fp2_dbl::mul_xi(za);
                                za += ad;
                                cf = fp2_dbl::mul_xi(cf);
                                zb += cf;
                                zc += be;

                                return fp6_dbl(za, zb, zc);
                            }

                            static fp6_dbl sqr_pre(const underlying_type &x) {
                                const non_residue_type &a = x.data[0];
                                const non_residue_type &b = x.data[1];
                                const non_residue_type &c = x.data[2];

                                const non_residue_type two_b = b.doubled();
                                fp2_dbl bc2 = fp2_dbl::mul_pre(two_b, c);
                                const fp2_dbl ab2 = fp2_dbl::mul_pre(two_b, a);
                                const fp2_dbl aa = fp2_dbl::sqr_pre(a);
                                fp2_dbl cc = fp2_dbl::sqr_pre(c);
                                fp2_dbl t = fp2_dbl::sqr_pre(a + b + c);
                                t -= aa;
                                t -= bc2;
                                t -= cc;
                                const fp2_dbl yc = t - ab2;
                                bc2 = fp2_dbl::mul_xi(bc2);
                                const fp2_dbl ya = aa + bc2;
                                cc = fp2_dbl::mul_xi(cc);
                                const fp2_dbl yb = cc + ab2;

                                return fp6_dbl(ya, yb, yc);
                            }

                            static underlying_type reduce(const fp6_dbl &x) {
                                return underlying_type(fp2_dbl::reduce(x.data[0]),
                                                       fp2_dbl::reduce(x.data[1]),
                                                       fp2_dbl::reduce(x.data[2]));
                            }
                        };

                        static non_residue_type mul_by_xi(const non_residue_type &x) {
                            return non_residue_type(
                                x.data[0].doubled().doubled().doubled() + x.data[0] - x.data[1],
                                x.data[1].doubled().doubled().doubled() + x.data[1] + x.data[0]);
                        }

                        static underlying_type mul_by_v(const underlying_type &x) {
                            return underlying_type(mul_by_xi(x.data[2]), x.data[0], x.data[1]);
                        }

                        static underlying_type mul_v_add(const underlying_type &x,
                                                         const underlying_type &y) {
                            return underlying_type(mul_by_xi(x.data[2]) + y.data[0],
                                                   x.data[0] + y.data[1],
                                                   x.data[1] + y.data[2]);
                        }

                        static fp6_dbl mul_v_add(const fp6_dbl &x, const fp6_dbl &y) {
                            return fp6_dbl(fp2_dbl::mul_xi(x.data[2]) + y.data[0],
                                           x.data[0] + y.data[1],
                                           x.data[1] + y.data[2]);
                        }

                        template<typename Fp12Value>
                        static Fp12Value multiply(const Fp12Value &x, const Fp12Value &y) {
                            const underlying_type &a = x.data[0];
                            const underlying_type &b = x.data[1];
                            const underlying_type &c = y.data[0];
                            const underlying_type &d = y.data[1];

                            const fp6_dbl ac = fp6_dbl::mul_pre(a, c);
                            const fp6_dbl bd = fp6_dbl::mul_pre(b, d);
                            const underlying_type z0 = fp6_dbl::reduce(mul_v_add(bd, ac));

                            fp6_dbl z1 = fp6_dbl::mul_pre(a + b, c + d);
                            z1 -= ac;
                            z1 -= bd;

                            return Fp12Value(z0, fp6_dbl::reduce(z1));
                        }
                    };

                    template<std::size_t Version>
                    constexpr typename fp12_2over3over2_extension_params<
                        alt_bn128_base_field<Version>>::non_residue_type const
                    fp12_2over3over2_extension_params<alt_bn128_base_field<Version>>::non_residue;

                    template<std::size_t Version>
                    constexpr std::array<
                        typename fp12_2over3over2_extension_params<alt_bn128_base_field<Version>>::integral_type,
                        12 * 2> const
                    fp12_2over3over2_extension_params<alt_bn128_base_field<Version>>::Frobenius_coeffs_c1;
                } // namespace detail
            } // namespace fields
        } // namespace algebra
    } // namespace crypto3
} // namespace nil

#endif    // CRYPTO3_ALGEBRA_FIELDS_ALT_BN128_FP12_2OVER3OVER2_EXTENSION_PARAMS_HPP
