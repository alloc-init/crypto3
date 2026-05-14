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
                        typedef boost::multiprecision::number<
                            boost::multiprecision::cpp_int_backend<
                                2 * base_field_type::modulus_bits + 32,
                                2 * base_field_type::modulus_bits + 32,
                                boost::multiprecision::signed_magnitude,
                                boost::multiprecision::unchecked,
                                void>>
                            wide_integral_type;
                        typedef typename integral_type::backend_type::cpp_int_type base_cpp_int_type;

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

                        static wide_integral_type modulus_as_wide() {
                            return wide_integral_type(modulus.backend().to_cpp_int());
                        }

                        static wide_integral_type as_wide(const base_value_type &x) {
                            return wide_integral_type(x.data.backend().convert_to_cpp_int());
                        }

                        static base_value_type reduce_base(wide_integral_type x) {
                            const wide_integral_type p = modulus_as_wide();
                            x %= p;
                            if (x < 0) {
                                x += p;
                            }
                            const base_cpp_int_type reduced(x.backend());
                            const integral_type reduced_value{typename integral_type::backend_type(reduced)};
                            return base_value_type(reduced_value);
                        }

                        struct fp_dbl {
                            wide_integral_type data;

                            fp_dbl() = default;
                            explicit fp_dbl(const wide_integral_type &in_data) : data(in_data) {}

                            fp_dbl operator+(const fp_dbl &other) const {
                                return fp_dbl(data + other.data);
                            }

                            fp_dbl operator-(const fp_dbl &other) const {
                                return fp_dbl(data - other.data);
                            }

                            fp_dbl &operator+=(const fp_dbl &other) {
                                data += other.data;
                                return *this;
                            }

                            fp_dbl &operator-=(const fp_dbl &other) {
                                data -= other.data;
                                return *this;
                            }

                            static fp_dbl mul_pre(const base_value_type &x,
                                                  const base_value_type &y) {
                                return fp_dbl(as_wide(x) * as_wide(y));
                            }

                            static fp_dbl sqr_pre(const base_value_type &x) {
                                const wide_integral_type x_wide = as_wide(x);
                                return fp_dbl(x_wide * x_wide);
                            }

                            static base_value_type reduce(const fp_dbl &x) {
                                return reduce_base(x.data);
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
                                const wide_integral_type a = as_wide(x.data[0]);
                                const wide_integral_type b = as_wide(x.data[1]);
                                const wide_integral_type c = as_wide(y.data[0]);
                                const wide_integral_type d = as_wide(y.data[1]);
                                const wide_integral_type ac = a * c;
                                const wide_integral_type bd = b * d;

                                return fp2_dbl(fp_dbl(ac - bd),
                                               fp_dbl((a + b) * (c + d) - ac - bd));
                            }

                            static fp2_dbl sqr_pre(const non_residue_type &x) {
                                const wide_integral_type a = as_wide(x.data[0]);
                                const wide_integral_type b = as_wide(x.data[1]);
                                return fp2_dbl(fp_dbl((a + b) * (a - b)), fp_dbl((a + a) * b));
                            }

                            static fp2_dbl mul_xi(const fp2_dbl &x) {
                                return fp2_dbl(fp_dbl(9 * x.data[0].data - x.data[1].data),
                                               fp_dbl(x.data[0].data + 9 * x.data[1].data));
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
