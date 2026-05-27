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

                        // Limb contract with fp12_limb_ops:
                        // - limb_bits == 64.
                        // - base_value_limb_count == 4: a normal BN254 Fp residue or modulus value.
                        // - storage_limb_count == 9: shared storage for base residues and bounded pre-REDC values.
                        constexpr static const std::size_t limb_bits = alt_bn128_fp12_limb_ops::limb_bits;
                        constexpr static const std::size_t base_value_limb_count =
                            alt_bn128_fp12_limb_ops::base_value_limb_count;
                        constexpr static const std::size_t storage_limb_count =
                            alt_bn128_fp12_limb_ops::storage_limb_count;

                        constexpr static const std::size_t field_limb_count =
                            (base_field_type::modulus_bits + limb_bits - 1) / limb_bits;
                        constexpr static const std::size_t bounded_product_limb_count =
                            (2 * base_field_type::modulus_bits + 32 + limb_bits - 1) / limb_bits;

                        static_assert(limb_bits == 64u, "alt_bn128 fp12 limb ops assume 64-bit limbs");
                        static_assert(base_value_limb_count == 4u,
                                      "alt_bn128 fp12 limb ops assume four base-value limbs");
                        static_assert(storage_limb_count == 9u, "alt_bn128 fp12 limb ops assume nine storage limbs");
                        static_assert(field_limb_count == base_value_limb_count,
                                      "alt_bn128 fp12 fast path expects four 64-bit base-field limbs");
                        static_assert(bounded_product_limb_count == storage_limb_count,
                                      "alt_bn128 fp12 fast path expects nine lazy product limbs");

                        // Base Fp values are kept in the nine-limb storage shape so tower additions and lazy products
                        // use one container type. Only the low four limbs hold the actual base-field value.
                        typedef alt_bn128_fp12_limb_ops::limb_array base_limb_storage_type;
                        // Lazy pre-REDC values may use all nine storage limbs before Montgomery reduction.
                        typedef alt_bn128_fp12_limb_ops::limb_array lazy_limb_storage_type;

                        // Lazy, signed, double-width base-Fp value used inside the tower fast path.
                        //
                        // A normal base_value_type is a canonical Montgomery residue x*R mod p in the
                        // low four limbs. Multiplying two such residues first produces
                        // a double-width integer:
                        //   (x*R) * (y*R) = x*y*R^2.
                        // REDC/Montgomery reduction removes one factor of R and returns x*y*R mod p.
                        //
                        // The Fp2/Fp6/Fp12 formulas often add or subtract several raw products before
                        // the coefficient is needed as a normal Fp value, for example ac - bd in Fp2.
                        // fp_dbl stores those bounded pre-REDC expressions so the tower can combine
                        // products first and call REDC only for each final coefficient. The limb storage
                        // stores an unsigned magnitude, so `fp_dbl::negative` records the integer sign
                        // until reduce() maps the signed result back into reduced Montgomery limbs.
                        struct fp_dbl {
                            // Magnitude of a bounded pre-REDC expression. For BN254 this may occupy up to
                            // nine storage limbs, while normal Fp values use only the low four.
                            lazy_limb_storage_type data = {};
                            // Sign of the integer expression represented by data
                            // Zero is non-negative
                            bool negative = false;

                            fp_dbl() = default;

                            explicit fp_dbl(const lazy_limb_storage_type &in_data, bool in_negative = false) :
                                data(in_data), negative(in_negative) {
                                if (alt_bn128_fp12_limb_ops::is_zero(data)) {
                                    // Enforce that 0 is positive
                                    negative = false;
                                }
                            }

                            // If this value is not zero, flip the sign
                            fp_dbl operator-() const {
                                fp_dbl result(*this);
                                if (!alt_bn128_fp12_limb_ops::is_zero(result.data)) {
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
                                    alt_bn128_fp12_limb_ops::add_limbs(data, other.data);
                                    return *this;
                                }
                                // .. they have different signs
                                const int cmp = alt_bn128_fp12_limb_ops::compare_limbs(data, other.data);
                                if (cmp == 0) {
                                    // same value but different signs, result is 0
                                    data = {};
                                    negative = false;
                                    return *this;
                                }
                                if (cmp > 0) {
                                    // 'this' is greater, subtract the other
                                    // if 'this' is negative, result will be negative, and vice versa
                                    alt_bn128_fp12_limb_ops::subtract_limbs(data, data, other.data);
                                    return *this;
                                }
                                // .. 'this' is less than 'other' in magnitude
                                // if 'other' is negative, result will be negative, and vice versa
                                lazy_limb_storage_type magnitude;
                                alt_bn128_fp12_limb_ops::subtract_limbs(magnitude, other.data, data);
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
                                alt_bn128_fp12_limb_ops::left_shift_one(result.data);
                                return result;
                            }

                            fp_dbl &mul_by_9() {
                                alt_bn128_fp12_limb_ops::multiply_by_limb(data, alt_bn128_fp12_limb_ops::limb(9u));
                                return *this;
                            }

                            // "pre" means before Montgomery reduction
                            // x = aR and y = bR are Montgomery residues,
                            // and the result is the raw x*y product still scaled by R^2
                            // ie. xy = abR^2
                            // We leave it in this form because we have wide enough limb type to support extra
                            // additions and subtractions before reducing.
                            template<bool Wide = false>
                            static fp_dbl mul_pre(const base_limb_storage_type &x, const base_limb_storage_type &y) {
                                fp_dbl product;
                                if constexpr (Wide) {
                                    // The Fp2 Karatsuba cross term multiplies coefficient sums; limb_ops reads
                                    // the low five limbs because each sum may carry once past the four-limb field.
                                    alt_bn128_fp12_limb_ops::multiply_5x5(product.data, x, y);
                                } else {
                                    // Most tower products multiply normal base values, so limb_ops reads only
                                    // the low four limbs from each nine-limb storage value.
                                    alt_bn128_fp12_limb_ops::multiply_4x4(product.data, x, y);
                                }
                                return product;
                            }

                            void reduce() {
                                // Convert a bounded signed nine-limb pre-REDC expression to reduced four-limb
                                // Montgomery limbs. REDC removes one Montgomery factor; a negative integer
                                // representative is then mapped to p - reduced, which is the same value modulo p.
                                alt_bn128_fp12_limb_ops::montgomery_reduce<base_field_type>(data);
                                if (negative && !alt_bn128_fp12_limb_ops::is_zero(data)) {
                                    // if this fp_dbl went negative, compute x = p - x
                                    static const base_limb_storage_type mod =
                                        alt_bn128_fp12_limb_ops::load_limbs(extension_policy::modulus.backend());
                                    alt_bn128_fp12_limb_ops::subtract_limbs(data, mod, data);
                                }
                                negative = false;
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

                        struct fp2_base {
                            std::array<base_limb_storage_type, 2> data;

                            fp2_base() = default;

                            fp2_base(const base_limb_storage_type &c0, const base_limb_storage_type &c1) :
                                data({c0, c1}) {
                            }

                            fp2_base(const non_residue_type &x) :
                                data({alt_bn128_fp12_limb_ops::load_limbs(x.data[0].data.backend().base_data()),
                                      alt_bn128_fp12_limb_ops::load_limbs(x.data[1].data.backend().base_data())}) {
                            }

                            fp2_base &operator+=(const fp2_base &other) {
                                alt_bn128_fp12_limb_ops::add_limbs(data[0], other.data[0]);
                                alt_bn128_fp12_limb_ops::add_limbs(data[1], other.data[1]);
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

                            template<bool Wide = false>
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
                                const fp_dbl ac = fp_dbl::mul_pre(a, c);
                                const fp_dbl bd = fp_dbl::mul_pre(b, d);
                                base_limb_storage_type a_plus_b = a;
                                base_limb_storage_type c_plus_d = c;
                                alt_bn128_fp12_limb_ops::add_limbs(a_plus_b, b);
                                alt_bn128_fp12_limb_ops::add_limbs(c_plus_d, d);
                                result.data[0] = ac;
                                result.data[0] -= bd;
                                result.data[1] = fp_dbl::template mul_pre<Wide>(a_plus_b, c_plus_d);
                                result.data[1] -= ac;
                                result.data[1] -= bd;
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

                            template<bool Wide = false>
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
                                // Use Wide when x and y already contain sums, so inner Fp2 cross terms need five limbs.
                                fp2_dbl::template mul_pre<Wide>(za, b + c, e + f);
                                fp2_dbl::template mul_pre<Wide>(zb, a + b, e + d);
                                fp2_dbl::template mul_pre<Wide>(zc, a + c, d + f);
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
                            // two final Fp6 reductions. mul_pre<true> computes (a+b)(c+d)
                            // directly in the lazy tower, avoiding generic Fp6 temporaries and
                            // intermediate Montgomery reductions.
                            const fp6_base a(x.data[0]);
                            const fp6_base b(x.data[1]);
                            const fp6_base c(y.data[0]);
                            const fp6_base d(y.data[1]);

                            // false = ordinary Fp6 multiplication; inner Fp2 sums fit mul_pre().
                            fp6_dbl ac = fp6_dbl::mul_pre(a, c);
                            fp6_dbl bd = fp6_dbl::mul_pre(b, d);

                            fp6_dbl z0_dbl = ac + bd.mul_v();

                            // Inner Fp2 sums need the wide multiplication path because they already include two
                            // rounds of addition.
                            fp6_dbl z1_dbl = fp6_dbl::template mul_pre<true>(a + b, c + d);
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
