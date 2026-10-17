#pragma once

#include <array>
#include <cstddef>
#include <iterator>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/iterator/transform_iterator.hpp>

#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/types.hpp>
#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/limb_ops.hpp>

namespace nil::crypto3::algebra::fields::detail::fp12_fast {

    template<Fp12FastParams Params>
    struct fp12_fast {
        using limb_array = Params::limb_array;
        using base_value_type = typename Params::base_field_type::value_type;
        using base_backend_type =
            std::remove_cvref_t<decltype(std::declval<base_value_type &>().data.backend().base_data())>;

        // Three Fp6 products, six Fp2 products each, and three base-field products per Fp2 product.
        // Every entry is a canonical Montgomery residue. These are reusable linear forms of one operand;
        // multiplying corresponding entries and applying recombine_products gives the Fp12 product.
        static constexpr std::size_t product_coordinate_count = 54;
        using prepared_product_operand = std::array<base_backend_type, product_coordinate_count>;

        // Packed input-side Fp2 value: two canonical Fp coefficients.
        struct fp2_base {
            limb_array data;

            fp2_base() = default;

            fp2_base(const Params::non_residue_type &x) {
                constexpr size_t N = Params::base_value_limb_count;
                for (size_t i = 0; i < N; i++) {
                    data[i] = x.data[0].data.backend().base_data().limbs()[i];
                }
                for (size_t i = 0; i < N; i++) {
                    data[i + N] = x.data[1].data.backend().base_data().limbs()[i];
                }
            }

            static void add_mod(fp2_base &z, const fp2_base &x, const fp2_base &y) {
                fp2_base_add_mod<Params>(z.data, x.data, y.data);
            }

            static void mul_xi(fp2_base &z, const fp2_base &x) {
                constexpr size_t N = Params::base_value_limb_count;
                limb_array result {};
                if constexpr (Params::xi == std::array {9, 1} && Params::u_squared == -1) {
                    std::array<limb, N> a_times_9, b_times_9;
                    mul_low_limbs_by_9_mod_portable<Params>(a_times_9.data(), x.data.data());
                    mul_low_limbs_by_9_mod_portable<Params>(b_times_9.data(), x.data.data() + N);
                    subtract_low_limbs_mod_portable<Params>(result.data(), a_times_9.data(), x.data.data() + N);
                    add_low_limbs_mod_portable<Params>(result.data() + N, x.data.data(), b_times_9.data());
                } else if constexpr (Params::xi == std::array {1, 1} && Params::u_squared == -1) {
                    subtract_low_limbs_mod_portable<Params>(result.data(), x.data.data(), x.data.data() + N);
                    add_low_limbs_mod_portable<Params>(result.data() + N, x.data.data(), x.data.data() + N);
                } else if constexpr (Params::xi == std::array {0, 1} && Params::u_squared == -5) {
                    std::array<limb, N> b_times_5;
                    const std::array<limb, N> zero {};
                    mul_low_limbs_by_5_mod_portable<Params>(b_times_5.data(), x.data.data() + N);
                    subtract_low_limbs_mod_portable<Params>(result.data(), zero.data(), b_times_5.data());
                    for (size_t i = 0; i < N; ++i) {
                        result[i + N] = x.data[i];
                    }
                } else {
                    static_assert(Params::xi != Params::xi, "unsupported Fp2 non-residue");
                }
                z.data = result;
            }

            fp2_base &operator+=(const fp2_base &other) {
                fp2_base::add_mod(*this, *this, other);
                return *this;
            }

            void to_non_residue(Params::non_residue_type &result) const {
                constexpr std::size_t N = Params::base_value_limb_count;
                for (std::size_t coefficient = 0; coefficient < 2; ++coefficient) {
                    for (std::size_t i = 0; i < N; ++i) {
                        result.data[coefficient].data.backend().base_data().limbs()[i] = data[coefficient * N + i];
                    }
                }
            }
        };

        struct fp2_dbl {
            // Lazy Fp2 value in the same coefficient order as generic fp2:
            //   data[0] + data[1] * u, with u^2 = -1.
            // Each coefficient is an unreduced double-width Fp value kept in the p * R residue class.
            std::array<limb_array, 2> data;

            fp2_dbl() = default;

            static void add_mod(fp2_dbl &z, const fp2_dbl &x, const fp2_dbl &y) {
                add_limbs_mod<Params>(z.data[0], x.data[0], y.data[0]);
                add_limbs_mod<Params>(z.data[1], x.data[1], y.data[1]);
            }

            static void sub_mod(fp2_dbl &z, const fp2_dbl &x, const fp2_dbl &y) {
                subtract_limbs_mod<Params>(z.data[0], x.data[0], y.data[0]);
                subtract_limbs_mod<Params>(z.data[1], x.data[1], y.data[1]);
            }

            fp2_dbl &operator+=(const fp2_dbl &other) {
                fp2_dbl::add_mod(*this, *this, other);
                return *this;
            }

            fp2_dbl &operator-=(const fp2_dbl &other) {
                fp2_dbl::sub_mod(*this, *this, other);
                return *this;
            }

            // Karatsuba correction: the real coefficient uses bounded subtraction,
            // while the imaginary coefficient is known nonnegative and can raw-subtract.
            void sub_pre(const fp2_dbl &other) {
                fp2_sub_pre<Params>(data, other.data);
            }

            static void add_mul_pre(fp2_dbl &result, const fp2_base &a, const fp2_base &b, const fp2_base &c,
                                    const fp2_base &d) {
                fp2_add_mul_pre<Params>(result.data, a.data, b.data, c.data, d.data);
            }

            static void mul_pre(fp2_dbl &result, const fp2_base &x, const fp2_base &y) {
                fp2_mul_pre<Params>(result.data, x.data, y.data);
            }

            // dst = src * xi + addend
            static void mul_xi_add(fp2_dbl &dst, const fp2_dbl &src, const fp2_dbl &addend) {
                fp2_mul_xi_add<Params>(dst.data, src.data, addend.data);
            }

            void to_non_residue(Params::non_residue_type &ret) const {
                montgomery_reduce<Params>((limb *)ret.data[0].data.backend().base_data().limbs(), data[0]);
                montgomery_reduce<Params>((limb *)ret.data[1].data.backend().base_data().limbs(), data[1]);
            }
        };

        // Packed input-side Fp6 value:
        //   data[0] + data[1] * v + data[2] * v^2.
        // Each coefficient is an fp2_base, so this is still an input-side
        // representation, not a lazy/pre-REDC result.
        struct fp6_base {
            std::array<fp2_base, 3> data;

            fp6_base() = default;

            fp6_base(const Params::underlying_type &x) : data({x.data[0], x.data[1], x.data[2]}) {
            }

            fp6_base(const fp2_base &c0, const fp2_base &c1, const fp2_base &c2) : data({c0, c1, c2}) {
            }

            std::tuple<const fp2_base &, const fp2_base &, const fp2_base &> coeffs() const {
                return {data[0], data[1], data[2]};
            }

            static void add_mod(fp6_base &z, const fp6_base &x, const fp6_base &y) {
                fp2_base::add_mod(z.data[0], x.data[0], y.data[0]);
                fp2_base::add_mod(z.data[1], x.data[1], y.data[1]);
                fp2_base::add_mod(z.data[2], x.data[2], y.data[2]);
            }

            static fp6_base mul_v(const fp6_base &x) {
                fp2_base c0;
                fp2_base::mul_xi(c0, x.data[2]);
                return fp6_base(c0, x.data[0], x.data[1]);
            }

            fp6_base operator+(const fp6_base &other) const {
                fp6_base result;
                fp6_base::add_mod(result, *this, other);
                return result;
            }
        };

        struct fp6_dbl {
            std::array<fp2_dbl, 3> data;

            fp6_dbl() = default;

            fp6_dbl(const fp2_dbl &c0, const fp2_dbl &c1, const fp2_dbl &c2) : data({c0, c1, c2}) {
            }

            std::tuple<const fp2_dbl &, const fp2_dbl &, const fp2_dbl &> coeffs() const {
                return {data[0], data[1], data[2]};
            }

            fp6_dbl operator+(const fp6_dbl &other) const {
                return fp6_dbl(data[0] + other.data[0], data[1] + other.data[1], data[2] + other.data[2]);
            }

            fp6_dbl operator-(const fp6_dbl &other) const {
                return fp6_dbl(data[0] - other.data[0], data[1] - other.data[1], data[2] - other.data[2]);
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

            static void mul_pre(fp6_dbl &result, const fp6_base &x, const fp6_base &y) {
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
                fp2_dbl &za = result.data[0];
                fp2_dbl &zb = result.data[1];
                fp2_dbl &zc = result.data[2];
                fp2_dbl::add_mul_pre(za, b, c, e, f);
                fp2_dbl::add_mul_pre(zb, a, b, e, d);
                fp2_dbl::add_mul_pre(zc, a, c, d, f);
                // Direct products reused by the three Karatsuba corrections.
                fp2_dbl be, cf, ad;
                fp2_dbl::mul_pre(be, b, e);
                fp2_dbl::mul_pre(cf, c, f);
                fp2_dbl::mul_pre(ad, a, d);
                // Finish the Karatsuba corrections
                za.sub_pre(be);
                za.sub_pre(cf);
                zb.sub_pre(ad);
                zb.sub_pre(be);
                zc.sub_pre(ad);
                zc.sub_pre(cf);
                // Fold the v^3 and v^4 terms back into the tower:
                //   z0 = ad + xi*za
                //   z1 = zb + xi*cf
                //   z2 = zc + be
                fp2_dbl::mul_xi_add(za, za, ad);
                fp2_dbl::mul_xi_add(zb, cf, zb);
                zc += be;
            }

            static void add_mul_pre(fp6_dbl &result, const fp6_base &a, const fp6_base &b, const fp6_base &c,
                                    const fp6_base &d) {
                fp6_base x, y;
                fp6_base::add_mod(x, a, b);
                fp6_base::add_mod(y, c, d);
                fp6_dbl::mul_pre(result, x, y);
            }

            void to_underlying(Params::underlying_type &ret) const {
                data[0].to_non_residue(ret.data[0]);
                data[1].to_non_residue(ret.data[1]);
                data[2].to_non_residue(ret.data[2]);
            }

            // Fp6 multiply-by-v is a coefficient rotation with one xi multiplication because v^3 = xi:
            //   (a + b*v + c*v^2) * v
            //     = a*v + b*v^2 + c*v^3
            //     = xi*c + a*v + b*v^2
            // Computes z = x * v + y and folds x.c * v^3 as xi * x.c.
            static void mul_v_add(fp6_dbl &z, const fp6_dbl &x, const fp6_dbl &y) {
                fp2_dbl::mul_xi_add(z.data[0], x.data[2], y.data[0]);
                z.data[1] = x.data[0];
                z.data[1] += y.data[1];
                z.data[2] = x.data[1];
                z.data[2] += y.data[2];
            }

            static void mul_v(fp6_dbl &z, const fp6_dbl &x) {
                fp2_dbl zero;
                zero.data = {};
                fp2_dbl::mul_xi_add(z.data[0], x.data[2], zero);
                z.data[1] = x.data[0];
                z.data[2] = x.data[1];
            }
        };

        template<typename Fp12Value>
        static prepared_product_operand prepare_product_operand(const Fp12Value &x) {
            using fp2_value_type = typename Params::non_residue_type;
            using fp6_value_type = typename Params::underlying_type;
            const std::array<fp6_value_type, 3> outer = {x.data[0], x.data[1], x.data[0] + x.data[1]};
            prepared_product_operand result;
            for (std::size_t w = 0; w < outer.size(); ++w) {
                const auto &a = outer[w].data;
                const std::array<fp2_value_type, 6> forms = {a[0], a[1], a[2], a[0] + a[1], a[0] + a[2], a[1] + a[2]};
                for (std::size_t v = 0; v < forms.size(); ++v) {
                    const auto &b = forms[v].data;
                    const base_value_type sum = b[0] + b[1];
                    const std::size_t offset = 18 * w + 3 * v;
                    result[offset] = b[0].data.backend().base_data();
                    result[offset + 1] = b[1].data.backend().base_data();
                    result[offset + 2] = sum.data.backend().base_data();
                }
            }
            return result;
        }

        // This map is linear over Fp, so it also recombines sums of the 54 coordinate products. All corrections
        // take place after the dot products, with canonical field arithmetic rather than unbounded lazy values.
        template<typename Fp12Value>
        static Fp12Value recombine_products(const std::array<base_value_type, product_coordinate_count> &products) {
            using fp2_value_type = typename Params::non_residue_type;
            using fp6_value_type = typename Params::underlying_type;
            const auto mul_xi = [](const fp2_value_type &x) {
                fp2_base packed;
                fp2_base::mul_xi(packed, fp2_base(x));
                fp2_value_type result;
                packed.to_non_residue(result);
                return result;
            };
            std::array<fp6_value_type, 3> outer;
            for (std::size_t w = 0; w < outer.size(); ++w) {
                std::array<fp2_value_type, 6> p;
                for (std::size_t v = 0; v < p.size(); ++v) {
                    const std::size_t offset = 18 * w + 3 * v;
                    const auto &ac = products[offset];
                    const auto &bd = products[offset + 1];
                    auto real = ac - bd;
                    if constexpr (Params::u_squared == -5) {
                        real -= bd.doubled().doubled();
                    } else {
                        static_assert(Params::u_squared == -1, "unsupported quadratic non-residue");
                    }
                    p[v] = fp2_value_type(real, products[offset + 2] - ac - bd);
                }
                // p[0..2] are the diagonal Fp2 products; p[3..5] are the pairwise sum products.
                outer[w] = fp6_value_type(p[0] + mul_xi(p[5] - p[1] - p[2]), p[3] - p[0] - p[1] + mul_xi(p[2]),
                                          p[4] - p[0] - p[2] + p[1]);
            }
            const auto &ac = outer[0];
            const auto &bd = outer[1];
            const fp6_value_type v_bd(mul_xi(bd.data[2]), bd.data[0], bd.data[1]);
            return Fp12Value(ac + v_bd, outer[2] - ac - bd);
        }

        // Preconditions: both ranges contain prepared canonical operands; the right range is at least as long
        // as the left. Multipass iterators are required because each coordinate traverses the same operands.
        template<typename Fp12Value, std::forward_iterator LeftIterator, std::forward_iterator RightIterator>
        static Fp12Value inner_product(LeftIterator first, LeftIterator last, RightIterator right) {
            if (first == last) {
                return Fp12Value::zero();
            }
            std::array<base_value_type, product_coordinate_count> products;
            for (std::size_t coordinate = 0; coordinate < product_coordinate_count; ++coordinate) {
                const auto select = [coordinate](const prepared_product_operand &operand) -> const base_backend_type & {
                    return operand[coordinate];
                };
                // Reuse the existing 2*N+1-limb accumulator: canonical inputs are < p < R, so any size_t-length
                // sum fits. Its final canonicalization handles sums larger than p*R, including 123-term dots.
                montgomery_inner_product<Params>(
                    products[coordinate].data.backend().base_data(), boost::make_transform_iterator(first, select),
                    boost::make_transform_iterator(last, select), boost::make_transform_iterator(right, select));
            }
            return recombine_products<Fp12Value>(products);
        }

        template<typename Fp12Value>
        static Fp12Value square(const Fp12Value &x) {
            // Write x = a + b*w in Fp6[w]/(w^2 - v). The quadratic-extension
            // squaring formula
            //
            //   ab = a*b
            //   c0 = (a + b)(a + v*b) - ab - v*ab = a^2 + v*b^2
            //   c1 = 2*ab
            //
            // needs two Fp6 products instead of the three used by generic
            // Fp12 multiplication. Keep both products in the lazy doubled
            // representation until the two final Montgomery reductions.
            const fp6_base a(x.data[0]);
            const fp6_base b(x.data[1]);
            const fp6_base vb = fp6_base::mul_v(b);

            fp6_dbl ab, c0, v_ab;
            fp6_dbl::mul_pre(ab, a, b);
            fp6_dbl::add_mul_pre(c0, a, b, a, vb);
            fp6_dbl::mul_v(v_ab, ab);
            c0 -= ab;
            c0 -= v_ab;
            ab += ab;

            Fp12Value ret;
            c0.to_underlying(ret.data[0]);
            ab.to_underlying(ret.data[1]);
            return ret;
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
            // two final Fp6 reductions. add_mul_pre forms canonical input-side
            // Fp6 sums before the ordinary pre-REDC multiply path.
            const fp6_base a(x.data[0]);
            const fp6_base b(x.data[1]);
            const fp6_base c(y.data[0]);
            const fp6_base d(y.data[1]);
            fp6_dbl ac, bd, z;
            fp6_dbl::mul_pre(ac, a, c);
            fp6_dbl::mul_pre(bd, b, d);
            fp6_dbl::add_mul_pre(z, a, b, c, d);
            z -= ac;    // first correction (see above)
            z -= bd;    // second correction
            Fp12Value ret;
            z.to_underlying(ret.data[1]);
            fp6_dbl::mul_v_add(z, bd, ac);
            z.to_underlying(ret.data[0]);

            return ret;
        }
    };
}    // namespace nil::crypto3::algebra::fields::detail::fp12_fast
