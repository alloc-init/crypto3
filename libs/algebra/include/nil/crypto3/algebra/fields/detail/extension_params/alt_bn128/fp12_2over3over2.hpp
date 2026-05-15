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
#include <nil/crypto3/multiprecision/modular/modular_policy_fixed.hpp>

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

                    template<size_t Version>
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
                            0x23BD9E3DA9136A739F668E1ADC9EF7F0F575EC93F71A8DF953C846338C32A1AB_cppui_modular254};

                        /////////////////////////////////////////////////////

                        typedef boost::multiprecision::limb_type limb_type;
                        typedef typename base_field_type::modular_backend modular_backend_type;
                        typedef boost::multiprecision::backends::modular_policy<modular_backend_type>
                            modular_policy_type;
                        typedef typename modular_policy_type::Backend_doubled_padded_limbs reduction_backend_type;
                        typedef reduction_backend_type lazy_backend_type;

                        constexpr static const size_t limb_bits = sizeof(limb_type) * CHAR_BIT;
                        constexpr static const size_t base_limb_count =
                            (base_field_type::modulus_bits + limb_bits - 1) / limb_bits;
                        constexpr static const size_t wide_limb_count =
                            (2 * base_field_type::modulus_bits + 32 + limb_bits - 1) / limb_bits;

                        typedef lazy_backend_type base_limb_array_type;
                        typedef lazy_backend_type wide_limb_array_type;
                        static_assert(reduction_backend_type::internal_limb_count >= wide_limb_count,
                                      "Fp12 lazy product must fit the padded Montgomery reduction backend");

                        // Tower: Fp2 = Fp[u]/(u^2 + 1), Fp6 = Fp2[v]/(v^3 - xi),
                        // Fp12 = Fp6[w]/(w^2 - v). Here xi = 9 + u.
                        constexpr static const non_residue_type non_residue = non_residue_type(0x09, 0x01);

                        template<typename Backend>
                        static lazy_backend_type backend_to_base_limbs(const Backend &backend) {
                            lazy_backend_type result;
                            for (size_t i = 0; i < base_limb_count && i < backend.size(); ++i) {
                                result.limbs()[i] = backend.limbs()[i];
                            }
                            result.zero_after(base_limb_count);
                            result.set_carry(false);
                            result.normalize();
                            return result;
                        }

                        static const base_limb_array_type &modulus_limbs() {
                            static const base_limb_array_type value = backend_to_base_limbs(modulus.backend());
                            return value;
                        }

                        static lazy_backend_type as_base_limbs(const base_value_type &x) {
                            return backend_to_base_limbs(x.data.backend().base_data());
                        }

                        static lazy_backend_type montgomery_reduce_abs(const lazy_backend_type &x) {
                            reduction_backend_type backend = x;
                            base_field_type::modulus_params.get_mod_obj().montgomery_reduce(backend);
                            return backend_to_base_limbs(backend);
                        }

                        static base_value_type make_montgomery_base_value(const lazy_backend_type &x) {
                            typename integral_type::backend_type backend;
                            for (size_t i = 0; i < base_limb_count && i < backend.size(); ++i) {
                                backend.limbs()[i] = x.limbs()[i];
                            }
                            backend.zero_after(base_limb_count);

                            typename base_value_type::data_type data;
                            data.backend().base_data() = backend;
                            return base_value_type(data);
                        }

                        struct fp_dbl {
                            wide_limb_array_type data = {};
                            bool negative = false;

                            fp_dbl() = default;
                            explicit fp_dbl(const wide_limb_array_type &in_data, bool in_negative = false) :
                                data(in_data), negative(in_negative) {
                                normalize();
                            }

                            void normalize() {
                                if (data.compare(limb_type(0u)) == 0) {
                                    negative = false;
                                }
                            }

                            fp_dbl operator-() const {
                                fp_dbl result(*this);
                                if (result.data.compare(limb_type(0u)) != 0) {
                                    result.negative = !result.negative;
                                }
                                return result;
                            }

                            fp_dbl operator+(const fp_dbl &other) const {
                                fp_dbl result(*this);
                                result += other;
                                return result;
                            }

                            fp_dbl operator-(const fp_dbl &other) const {
                                fp_dbl result(*this);
                                result -= other;
                                return result;
                            }

                            fp_dbl &operator+=(const fp_dbl &other) {
                                if (negative == other.negative) {
                                    boost::multiprecision::backends::eval_add(data, other.data);
                                    return *this;
                                }
                                const int cmp = data.compare(other.data);
                                if (cmp == 0) {
                                    data = {};
                                    negative = false;
                                    return *this;
                                }
                                if (cmp > 0) {
                                    boost::multiprecision::backends::eval_subtract(data, other.data);
                                    return *this;
                                }
                                wide_limb_array_type magnitude = other.data;
                                boost::multiprecision::backends::eval_subtract(magnitude, data);
                                data = magnitude;
                                negative = other.negative;
                                return *this;
                            }

                            fp_dbl &operator-=(const fp_dbl &other) {
                                if (negative != other.negative) {
                                    boost::multiprecision::backends::eval_add(data, other.data);
                                    return *this;
                                }
                                const int cmp = data.compare(other.data);
                                if (cmp == 0) {
                                    data = {};
                                    negative = false;
                                    return *this;
                                }
                                if (cmp > 0) {
                                    boost::multiprecision::backends::eval_subtract(data, other.data);
                                    return *this;
                                }
                                wide_limb_array_type magnitude = other.data;
                                boost::multiprecision::backends::eval_subtract(magnitude, data);
                                data = magnitude;
                                negative = !other.negative;
                                return *this;
                            }

                            fp_dbl mul_small(limb_type factor) const {
                                if (factor == 0u || data.compare(limb_type(0u)) == 0) {
                                    return fp_dbl();
                                }
                                fp_dbl result(*this);
                                if (factor == 1u) {
                                    return result;
                                }
                                if (factor == 2u) {
                                    boost::multiprecision::backends::eval_left_shift(result.data, 1u);
                                    result.normalize();
                                    return result;
                                }
                                if (factor == 9u) {
                                    return result.mul_by_9_inplace();
                                }
                                boost::multiprecision::backends::eval_multiply(result.data, factor);
                                result.normalize();
                                return result;
                            }

                            fp_dbl &mul_by_9_inplace() {
                                wide_limb_array_type shifted = data;
                                boost::multiprecision::backends::eval_left_shift(shifted, 3u);
                                boost::multiprecision::backends::eval_add(shifted, data);
                                data = shifted;
                                normalize();
                                return *this;
                            }

                            static fp_dbl mul_pre(const base_limb_array_type &x, const base_limb_array_type &y) {
                                wide_limb_array_type product = x;
                                boost::multiprecision::backends::eval_multiply(product, y);
                                return fp_dbl(product);
                            }

                            static fp_dbl mul_pre(const base_value_type &x, const base_value_type &y) {
                                return mul_pre(as_base_limbs(x), as_base_limbs(y));
                            }

                            static base_value_type reduce(const fp_dbl &x) {
                                base_limb_array_type reduced = montgomery_reduce_abs(x.data);
                                if (x.negative && reduced.compare(limb_type(0u)) != 0) {
                                    base_limb_array_type negated = modulus_limbs();
                                    boost::multiprecision::backends::eval_subtract(negated, reduced);
                                    reduced = negated;
                                }
                                return make_montgomery_base_value(reduced);
                            }
                        };

                        struct fp2_base {
                            std::array<base_limb_array_type, 2> data;

                            fp2_base() = default;
                            fp2_base(const base_limb_array_type &c0, const base_limb_array_type &c1) : data({c0, c1}) {
                            }

                            static fp2_base from(const non_residue_type &x) {
                                return fp2_base(as_base_limbs(x.data[0]), as_base_limbs(x.data[1]));
                            }

                            static fp2_base from_sum(const non_residue_type &x, const non_residue_type &y) {
                                fp2_base result = from(x);
                                const base_limb_array_type y0 = as_base_limbs(y.data[0]);
                                const base_limb_array_type y1 = as_base_limbs(y.data[1]);
                                boost::multiprecision::backends::eval_add(result.data[0], y0);
                                boost::multiprecision::backends::eval_add(result.data[1], y1);
                                return result;
                            }

                            fp2_base operator+(const fp2_base &other) const {
                                fp2_base result(*this);
                                boost::multiprecision::backends::eval_add(result.data[0], other.data[0]);
                                boost::multiprecision::backends::eval_add(result.data[1], other.data[1]);
                                return result;
                            }
                        };

                        struct fp2_dbl {
                            std::array<fp_dbl, 2> data;

                            fp2_dbl() = default;
                            fp2_dbl(const fp_dbl &c0, const fp_dbl &c1) : data({c0, c1}) {
                            }

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

                            static void mul_pre(fp2_dbl &result,
                                                const non_residue_type &x,
                                                const non_residue_type &y) {
                                mul_pre(result, fp2_base::from(x), fp2_base::from(y));
                            }

                            static fp2_dbl mul_pre(const non_residue_type &x, const non_residue_type &y) {
                                fp2_dbl result;
                                mul_pre(result, x, y);
                                return result;
                            }

                            static void mul_pre(fp2_dbl &result, const fp2_base &x, const fp2_base &y) {
                                const base_limb_array_type &a = x.data[0];
                                const base_limb_array_type &b = x.data[1];
                                const base_limb_array_type &c = y.data[0];
                                const base_limb_array_type &d = y.data[1];
                                const fp_dbl ac = fp_dbl::mul_pre(a, c);
                                const fp_dbl bd = fp_dbl::mul_pre(b, d);
                                base_limb_array_type ab = a;
                                base_limb_array_type cd = c;
                                boost::multiprecision::backends::eval_add(ab, b);
                                boost::multiprecision::backends::eval_add(cd, d);

                                result.data[0] = ac;
                                result.data[0] -= bd;
                                result.data[1] = fp_dbl::mul_pre(ab, cd);
                                result.data[1] -= ac;
                                result.data[1] -= bd;
                            }

                            static fp2_dbl mul_pre(const fp2_base &x, const fp2_base &y) {
                                fp2_dbl result;
                                mul_pre(result, x, y);
                                return result;
                            }

                            static void mul_pre_loose_sum(fp2_dbl &result,
                                                          const fp2_base &x0,
                                                          const fp2_base &x1,
                                                          const fp2_base &y0,
                                                          const fp2_base &y1) {
                                // Used by the Fp12 cross term: the Fp2 inputs are already sums,
                                // so the inner Karatsuba sum needs one extra limb for range.
                                base_limb_array_type a = x0.data[0];
                                base_limb_array_type b = x0.data[1];
                                base_limb_array_type c = y0.data[0];
                                base_limb_array_type d = y0.data[1];
                                boost::multiprecision::backends::eval_add(a, x1.data[0]);
                                boost::multiprecision::backends::eval_add(b, x1.data[1]);
                                boost::multiprecision::backends::eval_add(c, y1.data[0]);
                                boost::multiprecision::backends::eval_add(d, y1.data[1]);
                                const fp_dbl ac = fp_dbl::mul_pre(a, c);
                                const fp_dbl bd = fp_dbl::mul_pre(b, d);
                                base_limb_array_type ab = a;
                                base_limb_array_type cd = c;
                                boost::multiprecision::backends::eval_add(ab, b);
                                boost::multiprecision::backends::eval_add(cd, d);

                                result.data[0] = ac;
                                result.data[0] -= bd;
                                result.data[1] = fp_dbl::mul_pre(ab, cd);
                                result.data[1] -= ac;
                                result.data[1] -= bd;
                            }

                            static fp2_dbl mul_pre_loose_sum(const fp2_base &x0,
                                                             const fp2_base &x1,
                                                             const fp2_base &y0,
                                                             const fp2_base &y1) {
                                fp2_dbl result;
                                mul_pre_loose_sum(result, x0, x1, y0, y1);
                                return result;
                            }

                            static fp2_dbl sqr_pre(const non_residue_type &x) {
                                return sqr_pre(fp2_base::from(x));
                            }

                            static fp2_dbl sqr_pre(const fp2_base &x) {
                                const base_limb_array_type &a = x.data[0];
                                const base_limb_array_type &b = x.data[1];
                                const fp_dbl aa = fp_dbl::mul_pre(a, a);
                                const fp_dbl bb = fp_dbl::mul_pre(b, b);
                                return fp2_dbl(aa - bb, fp_dbl::mul_pre(a, b).mul_small(2u));
                            }

                            static fp2_dbl mul_xi(const fp2_dbl &x) {
                                fp2_dbl result(x);
                                result.mul_xi_inplace();
                                return result;
                            }

                            void mul_xi_inplace() {
                                // Lazy multiply by xi = 9 + u:
                                // (a + b*u) * xi = (9a - b) + (a + 9b) * u.
                                const fp_dbl c0 = data[0];
                                data[0].mul_by_9_inplace();
                                data[0] -= data[1];
                                data[1].mul_by_9_inplace();
                                data[1] += c0;
                            }

                            static non_residue_type reduce(const fp2_dbl &x) {
                                return non_residue_type(fp_dbl::reduce(x.data[0]), fp_dbl::reduce(x.data[1]));
                            }
                        };

                        struct fp6_dbl {
                            std::array<fp2_dbl, 3> data;

                            fp6_dbl() = default;
                            fp6_dbl(const fp2_dbl &c0, const fp2_dbl &c1, const fp2_dbl &c2) : data({c0, c1, c2}) {
                            }

                            fp6_dbl operator+(const fp6_dbl &other) const {
                                return fp6_dbl(data[0] + other.data[0], data[1] + other.data[1],
                                               data[2] + other.data[2]);
                            }

                            fp6_dbl operator-(const fp6_dbl &other) const {
                                return fp6_dbl(data[0] - other.data[0], data[1] - other.data[1],
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

                            template<typename SumProduct>
                            static void mul_pre_impl(fp6_dbl &result,
                                                     const fp2_base &a,
                                                     const fp2_base &b,
                                                     const fp2_base &c,
                                                     const fp2_base &d,
                                                     const fp2_base &e,
                                                     const fp2_base &f,
                                                     SumProduct sum_product) {
                                // Fp6 is Fp2[v]/(v^3 - xi); multiplying by xi folds v^3
                                // terms back into the constant coefficient.
                                fp2_dbl za;
                                fp2_dbl zb;
                                fp2_dbl zc;
                                fp2_dbl be;
                                fp2_dbl cf;
                                fp2_dbl ad;
                                sum_product(za, b, c, e, f);
                                sum_product(zb, a, b, e, d);
                                sum_product(zc, a, c, d, f);
                                fp2_dbl::mul_pre(be, b, e);
                                fp2_dbl::mul_pre(cf, c, f);
                                fp2_dbl::mul_pre(ad, a, d);

                                za -= be;
                                za -= cf;
                                zb -= ad;
                                zb -= be;
                                zc -= ad;
                                zc -= cf;
                                za.mul_xi_inplace();
                                za += ad;
                                cf.mul_xi_inplace();
                                zb += cf;
                                zc += be;

                                result.data[0] = za;
                                result.data[1] = zb;
                                result.data[2] = zc;
                            }

                            static void mul_pre(fp6_dbl &result, const underlying_type &x, const underlying_type &y) {
                                const fp2_base a = fp2_base::from(x.data[0]);
                                const fp2_base b = fp2_base::from(x.data[1]);
                                const fp2_base c = fp2_base::from(x.data[2]);
                                const fp2_base d = fp2_base::from(y.data[0]);
                                const fp2_base e = fp2_base::from(y.data[1]);
                                const fp2_base f = fp2_base::from(y.data[2]);

                                mul_pre_impl(
                                    result,
                                    a, b, c, d, e, f,
                                    [](fp2_dbl &out,
                                       const fp2_base &x0,
                                       const fp2_base &x1,
                                       const fp2_base &y0,
                                       const fp2_base &y1) {
                                        fp2_base x_sum = x0 + x1;
                                        fp2_base y_sum = y0 + y1;
                                        fp2_dbl::mul_pre(out, x_sum, y_sum);
                                    });
                            }

                            static fp6_dbl mul_pre(const underlying_type &x, const underlying_type &y) {
                                fp6_dbl result;
                                mul_pre(result, x, y);
                                return result;
                            }

                            static void mul_pre_sum(fp6_dbl &result,
                                                    const underlying_type &x0,
                                                    const underlying_type &x1,
                                                    const underlying_type &y0,
                                                    const underlying_type &y1) {
                                // Computes (x0 + x1) * (y0 + y1) for the Fp12 Karatsuba
                                // cross term without constructing generic Fp6 field sums.
                                const fp2_base a = fp2_base::from_sum(x0.data[0], x1.data[0]);
                                const fp2_base b = fp2_base::from_sum(x0.data[1], x1.data[1]);
                                const fp2_base c = fp2_base::from_sum(x0.data[2], x1.data[2]);
                                const fp2_base d = fp2_base::from_sum(y0.data[0], y1.data[0]);
                                const fp2_base e = fp2_base::from_sum(y0.data[1], y1.data[1]);
                                const fp2_base f = fp2_base::from_sum(y0.data[2], y1.data[2]);

                                mul_pre_impl(
                                    result,
                                    a, b, c, d, e, f,
                                    [](fp2_dbl &out,
                                       const fp2_base &a0,
                                       const fp2_base &a1,
                                       const fp2_base &b0,
                                       const fp2_base &b1) {
                                        fp2_dbl::mul_pre_loose_sum(out, a0, a1, b0, b1);
                                    });
                            }

                            static fp6_dbl mul_pre_sum(const underlying_type &x0,
                                                       const underlying_type &x1,
                                                       const underlying_type &y0,
                                                       const underlying_type &y1) {
                                fp6_dbl result;
                                mul_pre_sum(result, x0, x1, y0, y1);
                                return result;
                            }

                            static fp6_dbl sqr_pre(const underlying_type &x) {
                                const fp2_base a = fp2_base::from(x.data[0]);
                                const fp2_base b = fp2_base::from(x.data[1]);
                                const fp2_base c = fp2_base::from(x.data[2]);
                                const fp2_base two_b = b + b;
                                fp2_dbl bc2 = fp2_dbl::mul_pre(two_b, c);
                                const fp2_dbl ab2 = fp2_dbl::mul_pre(two_b, a);
                                const fp2_dbl aa = fp2_dbl::sqr_pre(a);
                                fp2_dbl cc = fp2_dbl::sqr_pre(c);
                                fp2_dbl t = fp2_dbl::sqr_pre(a + b + c);
                                t -= aa;
                                t -= bc2;
                                t -= cc;
                                bc2.mul_xi_inplace();
                                cc.mul_xi_inplace();
                                const fp2_dbl ya = aa + bc2;
                                const fp2_dbl yb = cc + ab2;
                                const fp2_dbl yc = t - ab2;
                                return fp6_dbl(ya, yb, yc);
                            }

                            static underlying_type reduce(const fp6_dbl &x) {
                                return underlying_type(fp2_dbl::reduce(x.data[0]), fp2_dbl::reduce(x.data[1]),
                                                       fp2_dbl::reduce(x.data[2]));
                            }
                        };

                        static void mul_v_add(fp6_dbl &result, const fp6_dbl &x, const fp6_dbl &y) {
                            result.data[0] = x.data[2];
                            result.data[0].mul_xi_inplace();
                            result.data[0] += y.data[0];
                            result.data[1] = x.data[0];
                            result.data[1] += y.data[1];
                            result.data[2] = x.data[1];
                            result.data[2] += y.data[2];
                        }

                        template<typename Fp12Value>
                        static Fp12Value multiply(const Fp12Value &x, const Fp12Value &y) {
                            // For Fp12 = Fp6[w]/(w^2 - v):
                            // (a + b*w)(c + d*w) = (ac + bd*v) + ((a+b)(c+d)-ac-bd)*w.
                            const underlying_type &a = x.data[0];
                            const underlying_type &b = x.data[1];
                            const underlying_type &c = y.data[0];
                            const underlying_type &d = y.data[1];

                            fp6_dbl ac;
                            fp6_dbl bd;
                            fp6_dbl::mul_pre(ac, a, c);
                            fp6_dbl::mul_pre(bd, b, d);
                            fp6_dbl z0_dbl;
                            mul_v_add(z0_dbl, bd, ac);
                            const underlying_type z0 = fp6_dbl::reduce(z0_dbl);

                            fp6_dbl z1;
                            fp6_dbl::mul_pre_sum(z1, a, b, c, d);
                            z1 -= ac;
                            z1 -= bd;

                            return Fp12Value(z0, fp6_dbl::reduce(z1));
                        }
                    };

                    template<size_t Version>
                    constexpr typename fp12_2over3over2_extension_params<
                        alt_bn128_base_field<Version>>::non_residue_type const
                        fp12_2over3over2_extension_params<alt_bn128_base_field<Version>>::non_residue;

                    template<size_t Version>
                    constexpr std::array<
                        typename fp12_2over3over2_extension_params<alt_bn128_base_field<Version>>::integral_type,
                        12 * 2> const
                        fp12_2over3over2_extension_params<alt_bn128_base_field<Version>>::Frobenius_coeffs_c1;
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ALGEBRA_FIELDS_ALT_BN128_FP12_2OVER3OVER2_EXTENSION_PARAMS_HPP
