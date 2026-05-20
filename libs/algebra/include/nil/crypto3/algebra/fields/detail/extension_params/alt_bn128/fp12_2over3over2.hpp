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

                        // Tower: Fp2  = Fp[u]/(u^2 + 1)
                        //        Fp6  = Fp2[v]/(v^3 - xi),
                        //        Fp12 = Fp6[w]/(w^2 - v).
                        // Here xi = 9 + u.
                        constexpr static const non_residue_type non_residue = non_residue_type(0x09, 0x01);

                        /////////////////////////////////////////////////////

                        typedef boost::multiprecision::limb_type limb_type;
                        typedef boost::multiprecision::backends::modular_policy<
                            typename base_field_type::modular_backend>::Backend_doubled_padded_limbs
                            padded_limb_storage_type;

                        // The following aliases are the same padded backend storage.
                        // base_limb_storage_type means only the low base_value_limb_count limbs carry a base Fp value;
                        typedef padded_limb_storage_type base_limb_storage_type;
                        // lazy_limb_storage_type may use up to lazy_product_limb_count limbs for unreduced products and
                        // tower sums before Montgomery reduction.
                        typedef padded_limb_storage_type lazy_limb_storage_type;

                        // For BN254 with 64-bit limbs, base_value_limb_count is 4 and lazy_product_limb_count is 9.
                        constexpr static const size_t limb_bits = sizeof(limb_type) * CHAR_BIT;
                        constexpr static const size_t base_value_limb_count =
                            (base_field_type::modulus_bits + limb_bits - 1) / limb_bits;
                        constexpr static const size_t lazy_product_limb_count =
                            (2 * base_field_type::modulus_bits + 32 + limb_bits - 1) / limb_bits;

                        static_assert(padded_limb_storage_type::internal_limb_count >= lazy_product_limb_count,
                                      "Fp12 lazy product must fit the padded Montgomery reduction backend");

                        template<typename Backend>
                        static base_limb_storage_type backend_to_base_limb_storage(const Backend &backend) {
                            base_limb_storage_type result;
                            for (size_t i = 0; i < base_value_limb_count && i < backend.size(); ++i) {
                                result.limbs()[i] = backend.limbs()[i];
                            }
                            return result;
                        }

                        static base_limb_storage_type as_base_limbs(const base_value_type &x) {
                            return backend_to_base_limb_storage(x.data.backend().base_data());
                        }

                        // Lazy, signed, double-width base-Fp value used inside the tower fast path.
                        //
                        // A normal base_value_type is a canonical Montgomery residue x*R mod p in the
                        // low base_value_limb_count limbs. Multiplying two such residues first produces
                        // a double-width integer:
                        //   (x*R) * (y*R) = x*y*R^2.
                        // REDC/Montgomery reduction removes one factor of R and returns x*y*R mod p.
                        //
                        // The Fp2/Fp6/Fp12 formulas often add or subtract several raw products before
                        // the coefficient is needed as a normal Fp value, for example ac - bd in Fp2.
                        // fp_dbl stores those bounded pre-REDC expressions so the tower can combine
                        // products first and call REDC only for each final coefficient. The underlying
                        // backend stores an unsigned magnitude, so `fp_dbl::negative` records the integer
                        // sign until reduce() maps the signed result back into reduced Montgomery limbs.
                        struct fp_dbl {
                            // Magnitude of a bounded pre-REDC expression. For BN254 this may occupy up to
                            // lazy_product_limb_count limbs, while normal Fp values use only the low four.
                            lazy_limb_storage_type data = {};
                            // Sign of the integer expression represented by data
                            // Zero is non-negative
                            bool negative = false;

                            fp_dbl() = default;
                            explicit fp_dbl(const lazy_limb_storage_type &in_data, bool in_negative = false) :
                                data(in_data), negative(in_negative) {
                                if (data.compare(limb_type(0u)) == 0) {
                                    // Enforce that 0 is positive
                                    negative = false;
                                }
                            }

                            // If this value is not zero, flip the sign
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

                            fp_dbl &add_magnitude(const fp_dbl &other, bool other_negative) {
                                if (negative == other_negative) {
                                    // if they both have the same sign, you can just add
                                    boost::multiprecision::backends::eval_add(data, other.data);
                                    return *this;
                                }
                                // .. they have different signs
                                const int cmp = data.compare(other.data);
                                if (cmp == 0) {
                                    // same value but different signs, result is 0
                                    data = {};
                                    negative = false;
                                    return *this;
                                }
                                if (cmp > 0) {
                                    // 'this' is greater, subtract the other
                                    // if 'this' is negative, result will be negative, and vice versa
                                    boost::multiprecision::backends::eval_subtract(data, other.data);
                                    return *this;
                                }
                                // .. 'this' is less than 'other' in magnitude
                                // if 'other' is negative, result will be negative, and vice versa
                                lazy_limb_storage_type magnitude = other.data;
                                boost::multiprecision::backends::eval_subtract(magnitude, data);
                                data = magnitude;
                                negative = other_negative;
                                return *this;
                            }

                            fp_dbl &operator+=(const fp_dbl &other) {
                                return add_magnitude(other, other.negative);
                            }

                            fp_dbl &operator-=(const fp_dbl &other) {
                                return add_magnitude(other, !other.negative);
                            }

                            fp_dbl doubled() const {
                                fp_dbl result(*this);
                                boost::multiprecision::backends::eval_left_shift(result.data, 1u);
                                return result;
                            }

                            fp_dbl &mul_by_9() {
                                boost::multiprecision::backends::eval_multiply(data, limb_type(9u));
                                return *this;
                            }

                            // "pre" means before Montgomery reduction
                            // x = aR and y = bR are Montgomery residues,
                            // and the result is the raw x*y product still scaled by R^2
                            // ie. xy = abR^2
                            // We leave it in this form because we have wide enough limb type to support extra
                            // additions and subtractions before reducing.
                            static fp_dbl mul_pre_4limb(const base_limb_storage_type &x,
                                                        const base_limb_storage_type &y) {
                                fp_dbl product;
                                if (base_value_limb_count != 4u) {
                                    product.data = x;
                                    boost::multiprecision::backends::eval_multiply(product.data, y);
                                } else {
                                    // Most tower products multiply values below 4p, so they fit
                                    // in the low four 64-bit limbs. The loose Fp12 cross term has
                                    // a separate five-limb path for its wider Karatsuba sum.
                                    boost::multiprecision::backends::eval_multiply_4x4(product.data, x.limbs(),
                                                                                       y.limbs());
                                }
                                return product;
                            }

                            static fp_dbl mul_pre_5limb(const base_limb_storage_type &x,
                                                        const base_limb_storage_type &y) {
                                fp_dbl product;
                                if (base_value_limb_count != 4u) {
                                    product.data = x;
                                    boost::multiprecision::backends::eval_multiply(product.data, y);
                                } else {
                                    boost::multiprecision::backends::eval_multiply_low_limbs<5u>(product.data,
                                                                                                 x.limbs(), y.limbs());
                                }
                                return product;
                            }

                            void reduce() {
                                // Convert a bounded signed pre-REDC expression to reduced Montgomery limbs.
                                // REDC removes one Montgomery factor; a negative integer representative is
                                // then mapped to p - reduced, which is the same value modulo p.
                                base_field_type::modulus_params.get_mod_obj().montgomery_reduce(data);
                                if (negative && data.compare(limb_type(0u)) != 0) {
                                    static const base_limb_storage_type modulus_limbs =
                                        backend_to_base_limb_storage(modulus.backend());
                                    base_limb_storage_type negated = modulus_limbs;
                                    boost::multiprecision::backends::eval_subtract(negated, data);
                                    data = negated;
                                }
                                negative = false;
                            }

                            base_value_type to_base_value() const {
                                // The data limbs must already be reduced Montgomery base-Fp limbs. Construct
                                // base_value_type directly from those limbs to avoid converting them again.
                                base_value_type out;
                                typename integral_type::backend_type &backend = out.data.backend().base_data();
                                for (size_t i = 0; i < base_value_limb_count && i < backend.size(); ++i) {
                                    backend.limbs()[i] = data.limbs()[i];
                                }
                                backend.zero_after(base_value_limb_count);
                                return out;
                            }
                        };

                        struct fp2_base {
                            std::array<base_limb_storage_type, 2> data;

                            fp2_base() = default;
                            fp2_base(const base_limb_storage_type &c0, const base_limb_storage_type &c1) :
                                data({c0, c1}) {
                            }

                            static fp2_base from(const non_residue_type &x) {
                                return fp2_base(as_base_limbs(x.data[0]), as_base_limbs(x.data[1]));
                            }

                            static fp2_base from_sum(const non_residue_type &x, const non_residue_type &y) {
                                fp2_base result = from(x);
                                const base_limb_storage_type y0 = as_base_limbs(y.data[0]);
                                const base_limb_storage_type y1 = as_base_limbs(y.data[1]);
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
                            // Lazy Fp2 value in the same coefficient order as generic fp2:
                            //   data[0] + data[1] * u, with u^2 = -1.
                            // Each coefficient is an unreduced signed double-width Fp value.
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

                            static void mul_pre(fp2_dbl &result, const non_residue_type &x, const non_residue_type &y) {
                                mul_pre(result, fp2_base::from(x), fp2_base::from(y));
                            }

                            static fp2_dbl mul_pre(const non_residue_type &x, const non_residue_type &y) {
                                fp2_dbl result;
                                mul_pre(result, x, y);
                                return result;
                            }

                            static void mul_pre(fp2_dbl &result, const fp2_base &x, const fp2_base &y) {
                                // For x = a + bu and y = c + du:
                                //   xy = (a + bu) * (c + du)
                                //      = ac + adu + bcu + bdu^2
                                //      = ac + (ad + bc)u - bd      # since u^2 = -1
                                //      = (ac - bd) + (ad + bc)u
                                // Karatsuba computes the cross term with one product:
                                //   ad + bc = (a + b)(c + d) - ac - bd.
                                const base_limb_storage_type &a = x.data[0];
                                const base_limb_storage_type &b = x.data[1];
                                const base_limb_storage_type &c = y.data[0];
                                const base_limb_storage_type &d = y.data[1];
                                const fp_dbl ac = fp_dbl::mul_pre_4limb(a, c);
                                const fp_dbl bd = fp_dbl::mul_pre_4limb(b, d);
                                base_limb_storage_type ab = a;
                                base_limb_storage_type cd = c;
                                boost::multiprecision::backends::eval_add(ab, b);
                                boost::multiprecision::backends::eval_add(cd, d);

                                result.data[0] = ac;
                                result.data[0] -= bd;
                                result.data[1] = fp_dbl::mul_pre_4limb(ab, cd);
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
                                // Same Karatsuba Fp2 multiply as mul_pre(), but used when the
                                // logical inputs are (x0 + x1) and (y0 + y1). This is the
                                // Fp12 cross term path, so we form those Fp2 sums here and
                                // avoid constructing generic field elements or reducing them.
                                base_limb_storage_type a = x0.data[0];
                                base_limb_storage_type b = x0.data[1];
                                base_limb_storage_type c = y0.data[0];
                                base_limb_storage_type d = y0.data[1];
                                boost::multiprecision::backends::eval_add(a, x1.data[0]);
                                boost::multiprecision::backends::eval_add(b, x1.data[1]);
                                boost::multiprecision::backends::eval_add(c, y1.data[0]);
                                boost::multiprecision::backends::eval_add(d, y1.data[1]);
                                const fp_dbl ac = fp_dbl::mul_pre_4limb(a, c);
                                const fp_dbl bd = fp_dbl::mul_pre_4limb(b, d);
                                base_limb_storage_type &ab = a;
                                base_limb_storage_type &cd = c;
                                boost::multiprecision::backends::eval_add(ab, b);
                                boost::multiprecision::backends::eval_add(cd, d);

                                result.data[0] = ac;
                                result.data[0] -= bd;
                                // After the outer sums, (a+b) and (c+d) can use a fifth limb,
                                // so this product uses the wider fixed-limb helper.
                                result.data[1] = fp_dbl::mul_pre_5limb(ab, cd);
                                result.data[1] -= ac;
                                result.data[1] -= bd;
                            }

                            static fp2_dbl sqr_pre(const fp2_base &x) {
                                // (a + b*u)^2 = (a^2 - b^2) + 2ab*u, again using u^2 = -1.
                                const base_limb_storage_type &a = x.data[0];
                                const base_limb_storage_type &b = x.data[1];
                                const fp_dbl aa = fp_dbl::mul_pre_4limb(a, a);
                                const fp_dbl bb = fp_dbl::mul_pre_4limb(b, b);
                                return fp2_dbl(aa - bb, fp_dbl::mul_pre_4limb(a, b).doubled());
                            }

                            void mul_by_xi() {
                                // Lazy multiply by xi = 9 + u:
                                // (a + b*u) * xi = (9a - b) + (a + 9b) * u.
                                const fp_dbl tmp_a = data[0];
                                data[0].mul_by_9();
                                data[0] -= data[1];
                                data[1].mul_by_9();
                                data[1] += tmp_a;
                            }

                            void reduce() {
                                data[0].reduce();
                                data[1].reduce();
                            }

                            non_residue_type to_non_residue() const {
                                return non_residue_type(data[0].to_base_value(), data[1].to_base_value());
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
                                za.mul_by_xi();
                                za += ad;
                                cf.mul_by_xi();
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

                                mul_pre_impl(result, a, b, c, d, e, f,
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
                                    result, a, b, c, d, e, f,
                                    [](fp2_dbl &out,
                                       const fp2_base &a0,
                                       const fp2_base &a1,
                                       const fp2_base &b0,
                                       const fp2_base &b1) { fp2_dbl::mul_pre_loose_sum(out, a0, a1, b0, b1); });
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
                                bc2.mul_by_xi();
                                cc.mul_by_xi();
                                const fp2_dbl ya = aa + bc2;
                                const fp2_dbl yb = cc + ab2;
                                const fp2_dbl yc = t - ab2;
                                return fp6_dbl(ya, yb, yc);
                            }

                            void reduce() {
                                data[0].reduce();
                                data[1].reduce();
                                data[2].reduce();
                            }

                            underlying_type to_underlying() const {
                                return underlying_type(data[0].to_non_residue(), data[1].to_non_residue(),
                                                       data[2].to_non_residue());
                            }
                        };

                        static void mul_v_add(fp6_dbl &result, const fp6_dbl &x, const fp6_dbl &y) {
                            result.data[0] = x.data[2];
                            result.data[0].mul_by_xi();
                            result.data[0] += y.data[0];
                            result.data[1] = x.data[0];
                            result.data[1] += y.data[1];
                            result.data[2] = x.data[1];
                            result.data[2] += y.data[2];
                        }

                        template<typename Fp12Value>
                        static Fp12Value multiply(const Fp12Value &x, const Fp12Value &y) {
                            // Tower layout:
                            //   Fp2  = Fp[u]  / (u^2 + 1)
                            //   Fp6  = Fp2[v] / (v^3 - xi), xi = 9 + u
                            //   Fp12 = Fp6[w] / (w^2 - v)
                            // xi is the chosen Fp2 cubic non-residue; adjoining v with
                            // v^3 = xi turns Fp2 into the degree-3 extension Fp6.
                            //
                            // Write x = a + b*w and y = c + d*w with a,b,c,d in Fp6.
                            // Since w^2 = v, the constant coefficient gets bd multiplied by v:
                            //   z0 = ac + v*bd
                            //   z1 = ad + bc = (a + b)(c + d) - ac - bd
                            //
                            // The Fp6 multiply-by-v is just a coefficient rotation with one
                            // xi multiplication because v^3 = xi:
                            //   (r0 + r1*v + r2*v^2) * v = xi*r2 + r0*v + r1*v^2.
                            //
                            // ac, bd, and z1 stay in the lazy doubled representation until
                            // the two final Fp6 reductions. mul_pre_sum computes the
                            // (a+b)(c+d) term directly in the lazy tower, avoiding generic
                            // Fp6 temporaries and intermediate Montgomery reductions.
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
                            z0_dbl.reduce();
                            const underlying_type z0 = z0_dbl.to_underlying();

                            fp6_dbl z1;
                            fp6_dbl::mul_pre_sum(z1, a, b, c, d);
                            z1 -= ac;
                            z1 -= bd;
                            z1.reduce();

                            return Fp12Value(z0, z1.to_underlying());
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
