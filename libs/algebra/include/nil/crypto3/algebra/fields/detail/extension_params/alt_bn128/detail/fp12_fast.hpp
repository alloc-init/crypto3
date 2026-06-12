#pragma once

#include <array>
#include <cstddef>
#include <tuple>

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_ops.hpp>

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {
                namespace detail {

                    template<typename BaseField, typename ExtensionParams>
                    struct alt_bn128_fp12_fast_multiply {
                        typedef BaseField base_field_type;
                        typedef ExtensionParams extension_policy;
                        typedef typename extension_policy::integral_type integral_type;
                        typedef typename extension_policy::base_value_type base_value_type;      // fp
                        typedef typename extension_policy::non_residue_type non_residue_type;    // fp2
                        typedef typename extension_policy::underlying_type underlying_type;      // fp6

                        using limb_array = alt_bn128_fp12_limb_ops::limb_array;

                        // Lazy double-width base-Fp value used inside the tower fast path.
                        //
                        // A normal base_value_type is a canonical Montgomery residue x*R mod p in the
                        // low four limbs. Multiplying two such residues first produces
                        // a double-width integer:
                        //   (x*R) * (y*R) = x*y*R^2.
                        // REDC/Montgomery reduction removes one factor of R and returns x*y*R mod p.
                        //
                        // The Fp2/Fp6/Fp12 formulas often add or subtract several raw products before
                        // the coefficient is needed as a normal Fp value, for example ac - bd in Fp2.
                        // fp_dbl stores those bounded pre-REDC expressions modulo p * R, where
                        // R = 2^(64 * 4). This keeps subtractions non-negative without a separate sign bit.
                        struct fp_dbl {
                            // For BN254 this may occupy up to eight storage limbs. The ninth storage limb is
                            // kept available for product/reducer scratch, but normalized fp_dbl values keep it zero.
                            limb_array data = {};

                            fp_dbl() = default;

                            explicit fp_dbl(const limb_array &in_data) : data(in_data) {
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
                                alt_bn128_fp12_limb_ops::add_8_limbs_mod<base_field_type>(data, other.data);
                                return *this;
                            }

                            fp_dbl &operator-=(const fp_dbl &other) {
                                alt_bn128_fp12_limb_ops::subtract_8_limbs_mod<base_field_type>(data, other.data);
                                return *this;
                            }

                            fp_dbl doubled() const {
                                fp_dbl result(*this);
                                result += *this;
                                return result;
                            }

                            fp_dbl &mul_by_9() {
                                const fp_dbl x(*this);
                                fp_dbl t(*this);
                                t += x;    // 2x
                                t += t;    // 4x
                                t += t;    // 8x
                                t += x;    // 9x
                                *this = t;
                                return *this;
                            }

                            // "pre" means before Montgomery reduction
                            // x = aR and y = bR are Montgomery residues,
                            // and the result is the raw x*y product still scaled by R^2
                            // ie. xy = abR^2
                            // We leave it in this form because we have wide enough limb type to support extra
                            // additions and subtractions before reducing.
                            static fp_dbl mul_pre(const limb_array &x, const limb_array &y) {
                                fp_dbl product;
                                alt_bn128_fp12_limb_ops::multiply_4x4(product.data, x, y);
                                return product;
                            }

                            void reduce() {
                                // Convert a bounded pR-residue pre-REDC expression to reduced four-limb
                                // Montgomery limbs. REDC removes one Montgomery factor.
                                alt_bn128_fp12_limb_ops::montgomery_reduce<base_field_type>(data);
                            }

                            base_value_type to_base_value() const {
                                // The data limbs must already be reduced Montgomery base-Fp limbs. Construct
                                // base_value_type directly from those limbs to avoid converting them again.
                                base_value_type out;
                                typename integral_type::backend_type &backend = out.data.backend().base_data();
                                for (std::size_t i = 0; i < backend.size(); ++i) {
                                    backend.limbs()[i] = data[i];
                                }
                                backend.set_carry(false);
                                backend.normalize();
                                return out;
                            }
                        };

                        // fp2 before multiplying - equivalent to 2 regular Fp's
                        struct fp2_base {
                            std::array<limb_array, 2> data;

                            fp2_base() = default;

                            fp2_base(const limb_array &c0, const limb_array &c1) : data({c0, c1}) {
                            }

                            fp2_base(const non_residue_type &x) :
                                data({alt_bn128_fp12_limb_ops::load_limbs(x.data[0].data.backend().base_data()),
                                      alt_bn128_fp12_limb_ops::load_limbs(x.data[1].data.backend().base_data())}) {
                            }

                            fp2_base &operator+=(const fp2_base &other) {
                                // no fp_base type, do 4 limb addition manually here
                                alt_bn128_fp12_limb_ops::add_low_4_limbs_mod<base_field_type>(data[0], other.data[0]);
                                alt_bn128_fp12_limb_ops::add_low_4_limbs_mod<base_field_type>(data[1], other.data[1]);
                                return *this;
                            }

                            fp2_base operator+(const fp2_base &other) const {
                                fp2_base result(*this);
                                result += other;
                                return result;
                            }
                        };

                        struct fp2_dbl {
                            // Lazy Fp2 value in the same coefficient order as generic fp2:
                            //   data[0] + data[1] * u, with u^2 = -1.
                            // Each coefficient is an unreduced double-width Fp value represented modulo p * R.
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

                            static void mul_pre(fp2_dbl &result, const fp2_base &x, const fp2_base &y) {
                                alt_bn128_fp12_limb_ops::fp2_mul_pre<base_field_type>(
                                    result.data[0].data, result.data[1].data, x.data[0].data, x.data[1].data,
                                    y.data[0].data, y.data[1].data);
                            }

                            static fp2_dbl mul_pre(const fp2_base &x, const fp2_base &y) {
                                fp2_dbl result;
                                fp2_dbl::mul_pre(result, x, y);
                                return result;
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

                        // Raw-limb input view of an Fp6 value:
                        //   data[0] + data[1] * v + data[2] * v^2.
                        // Each coefficient is an fp2_base, so this is still an input-side
                        // representation, not a lazy/pre-REDC result.
                        struct fp6_base {
                            std::array<fp2_base, 3> data;

                            fp6_base() = default;

                            fp6_base(const fp2_base &c0, const fp2_base &c1, const fp2_base &c2) : data({c0, c1, c2}) {
                            }

                            explicit fp6_base(const underlying_type &x) :
                                data({fp2_base(x.data[0]), fp2_base(x.data[1]), fp2_base(x.data[2])}) {
                            }

                            std::tuple<const fp2_base &, const fp2_base &, const fp2_base &> coeffs() const {
                                return {data[0], data[1], data[2]};
                            }

                            fp6_base &operator+=(const fp6_base &other) {
                                data[0] += other.data[0];
                                data[1] += other.data[1];
                                data[2] += other.data[2];
                                return *this;
                            }

                            fp6_base operator+(const fp6_base &other) const {
                                fp6_base result(*this);
                                result += other;
                                return result;
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

                            static fp6_dbl mul_pre(const fp6_base &x, const fp6_base &y) {
                                // Multiply two Fp6 values in the tower Fp6 = Fp2[v]/(v^3 - xi):
                                //   x = a + b*v + c*v^2
                                //   y = d + e*v + f*v^2
                                //
                                // Expanding and replacing v^3 with xi gives:
                                //   z0 = a*d + xi*(b*f + c*e)
                                //   z1 = a*e + b*d + xi*(c*f)
                                //   z2 = a*f + b*e + c*d
                                //
                                // This function computes those cross terms with three Karatsuba-style
                                // sum products, so only six Fp2 products are needed instead of nine:
                                //   za = (b + c)(e + f) - b*e - c*f = b*f + c*e
                                //   zb = (a + b)(d + e) - a*d - b*e = a*e + b*d
                                //   zc = (a + c)(d + f) - a*d - c*f = a*f + c*d
                                const auto &[a, b, c] = x.coeffs();    // a, b, c are fp2_base
                                const auto &[d, e, f] = y.coeffs();
                                fp2_dbl za, zb, zc;
                                fp2_dbl::mul_pre(za, b + c, e + f);
                                fp2_dbl::mul_pre(zb, a + b, e + d);
                                fp2_dbl::mul_pre(zc, a + c, d + f);
                                // Direct products reused by the three Karatsuba corrections.
                                fp2_dbl be = fp2_dbl::mul_pre(b, e);
                                fp2_dbl cf = fp2_dbl::mul_pre(c, f);
                                fp2_dbl ad = fp2_dbl::mul_pre(a, d);
                                // Finish the Karatsuba corrections
                                za -= be;
                                za -= cf;
                                zb -= ad;
                                zb -= be;
                                zc -= ad;
                                zc -= cf;
                                // Fold the v^3 and v^4 terms back into the tower:
                                //   z0 = ad + xi*za
                                //   z1 = zb + xi*cf
                                //   z2 = zc + be
                                za.mul_by_xi();
                                za += ad;
                                cf.mul_by_xi();
                                zb += cf;
                                zc += be;
                                return fp6_dbl(za, zb, zc);
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

                            // Fp6 multiply-by-v is a coefficient rotation with one xi multiplication because v^3 = xi:
                            //   (a + b*v + c*v^2) * v
                            //     = a*v + b*v^2 + c*v^3
                            //     = xi*c + a*v + b*v^2
                            fp6_dbl mul_v() const {
                                fp6_dbl result;
                                result.data[0] = data[2];
                                result.data[0].mul_by_xi();
                                result.data[1] = data[0];
                                result.data[2] = data[1];
                                return result;
                            }
                        };

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
                            //   z = (a + bw)(c + dw)
                            //     = ac + adw + bcw + bdw^2
                            //     = ac + (ad + bc)w + bdv
                            //     = (ac + bdv) + (ad + bc)w
                            //
                            // So,
                            //   z0 = ac + v*bd
                            //   z1 = ad + bc
                            //
                            // And you can avoid computing ad and bc (3 muls instead of 4 total)
                            //   z1 = (a + b)(c + d) - ac - bd
                            //
                            // ac, bd, and z1 stay in the lazy doubled representation until the
                            // two final Fp6 reductions. The input-side fp6 sums are reduced modulo p,
                            // so z1 can use the ordinary pre-REDC multiply path.
                            const fp6_base a(x.data[0]);
                            const fp6_base b(x.data[1]);
                            const fp6_base c(y.data[0]);
                            const fp6_base d(y.data[1]);

                            // false = ordinary Fp6 multiplication; inner Fp2 sums fit mul_pre().
                            fp6_dbl ac = fp6_dbl::mul_pre(a, c);
                            fp6_dbl bd = fp6_dbl::mul_pre(b, d);

                            fp6_dbl z0_dbl = ac + bd.mul_v();

                            fp6_dbl z1_dbl = fp6_dbl::mul_pre(a + b, c + d);
                            z1_dbl -= ac;    // first correction (see above)
                            z1_dbl -= bd;    // second correction

                            // the whole point; delaying reduction until the very end
                            z0_dbl.reduce();
                            z1_dbl.reduce();

                            // convert back to generic crypto3 tower type
                            const underlying_type z0 = z0_dbl.to_underlying();
                            const underlying_type z1 = z1_dbl.to_underlying();

                            return Fp12Value(z0, z1);
                        }
                    };
                }    // namespace detail
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil
