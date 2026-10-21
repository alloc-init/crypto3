//---------------------------------------------------------------------------//
// Copyright (c) 2026
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

#define BOOST_TEST_MODULE algebra_fp12_test

#include <array>
#include <cstddef>
#include <memory>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <vector>

#include <boost/mpl/list.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/algebra/fields/field_algorithms.hpp>
#include <nil/crypto3/algebra/random_element.hpp>

namespace {
    namespace fields = nil::crypto3::algebra::fields;

    using alt_bn128_254_fp12 = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using bls12_377_fp12 = fields::fp12_2over3over2<fields::bls12<377>>;
    using bls12_381_fp12 = fields::fp12_2over3over2<fields::bls12<381>>;

    using fp12_field_types = boost::mpl::list<alt_bn128_254_fp12, bls12_377_fp12, bls12_381_fp12>;

    template<typename FieldType>
    struct fp12_field_name;

    template<>
    struct fp12_field_name<alt_bn128_254_fp12> {
        constexpr static std::string_view value = "alt_bn128_254";
    };

    template<>
    struct fp12_field_name<bls12_377_fp12> {
        constexpr static std::string_view value = "bls12_377";
    };

    template<>
    struct fp12_field_name<bls12_381_fp12> {
        constexpr static std::string_view value = "bls12_381";
    };

    constexpr std::size_t random_samples = 32;

    // element_fp12_2over3over2 uses policy hooks for optimized multiplication
    // and squaring when the policy provides them. This test-only policy reuses
    // the real field parameters but intentionally omits both hooks, so the
    // expected value is produced by the generic implementation rather than by
    // a duplicated reference formula in this file.
    template<typename Fp12Field>
    struct generic_fp12_policy {
        using optimized_policy = typename Fp12Field::extension_policy;

        using field_type = typename optimized_policy::field_type;
        using non_residue_type = typename optimized_policy::non_residue_type;
        using underlying_type = typename optimized_policy::underlying_type;

        inline constexpr static const non_residue_type non_residue = optimized_policy::non_residue;
    };

    template<typename Fp12Field>
    using generic_fp12_value_type = fields::detail::element_fp12_2over3over2<generic_fp12_policy<Fp12Field>>;

    template<typename Fp12Field>
    generic_fp12_value_type<Fp12Field> to_generic_fp12(const typename Fp12Field::value_type &x) {
        return generic_fp12_value_type<Fp12Field>(x.data[0], x.data[1]);
    }

    template<typename Fp12Field>
    typename Fp12Field::value_type from_generic_fp12(const generic_fp12_value_type<Fp12Field> &x) {
        return typename Fp12Field::value_type(x.data[0], x.data[1]);
    }

    template<typename Fp12Field>
    typename Fp12Field::value_type generic_fp12_mul(const typename Fp12Field::value_type &x,
                                                    const typename Fp12Field::value_type &y) {
        return from_generic_fp12<Fp12Field>(to_generic_fp12<Fp12Field>(x) * to_generic_fp12<Fp12Field>(y));
    }

    template<typename Fp12Field, typename Rng>
    typename Fp12Field::value_type random_fp12(Rng &rng) {
        return nil::crypto3::algebra::random_element<Fp12Field>(rng);
    }

    template<typename Fp12Field>
    typename Fp12Field::value_type boundary_fp12(std::size_t pattern) {
        using base_value_type = typename Fp12Field::base_field_type::value_type;
        using integral_type = typename Fp12Field::integral_type;

        const integral_type p_minus_one = Fp12Field::modulus - 1;
        const integral_type p_minus_two = Fp12Field::modulus - 2;
        typename Fp12Field::value_type value;

        for (std::size_t w = 0; w < 2; ++w) {
            for (std::size_t v = 0; v < 3; ++v) {
                for (std::size_t u = 0; u < 2; ++u) {
                    const std::size_t i = (w * 3 + v) * 2 + u;
                    base_value_type coefficient;
                    switch (pattern) {
                        case 0:
                            coefficient = base_value_type(p_minus_one);
                            break;
                        case 1:
                            coefficient = base_value_type(i % 2 == 0 ? p_minus_one : integral_type(0));
                            break;
                        case 2:
                            coefficient = base_value_type(i % 2 == 0 ? integral_type(0) : p_minus_one);
                            break;
                        case 3:
                            coefficient = base_value_type(u == 0 ? p_minus_two : p_minus_one);
                            break;
                        case 4:
                            coefficient = base_value_type(i % 3 == 0 ? p_minus_one : integral_type(1));
                            break;
                        default:
                            coefficient = base_value_type(i % 4 < 2 ? p_minus_two : integral_type(1));
                            break;
                    }
                    value.data[w].data[v].data[u] = coefficient;
                }
            }
        }
        return value;
    }

    // Execute the vector algorithm on scalar lanes as well, independently of the 128-bit-column reference.
    // These operations model IFMA's ignored upper operand bits and reject any accumulator overflow.
    template<std::size_t Lanes>
    struct scalar_ifma52_operations {
        using limb = fields::detail::fp12_fast::limb;
        using wide = fields::detail::fp12_fast::wide_limb;
        static constexpr std::size_t lanes = Lanes;
        static constexpr limb mask = (limb(1) << 52) - 1;
        using value_type = std::array<limb, lanes>;
        using mask_type = std::array<bool, lanes>;
        static value_type broadcast(limb x) {
            value_type r;
            r.fill(x);
            return r;
        }
        static value_type load(const limb *x) {
            value_type r;
            std::copy_n(x, lanes, r.begin());
            return r;
        }
        static void store(limb *x, const value_type &a) {
            std::copy(a.begin(), a.end(), x);
        }
        static value_type add(value_type a, const value_type &b) {
            for (std::size_t i = 0; i < lanes; ++i) {
                if (a[i] > std::numeric_limits<limb>::max() - b[i]) {
                    throw std::overflow_error("IFMA lane overflow");
                }
                a[i] += b[i];
            }
            return a;
        }
        static value_type subtract(value_type a, const value_type &b) {
            for (std::size_t i = 0; i < lanes; ++i)
                a[i] -= b[i];
            return a;
        }
        static value_type bit_and(value_type a, const value_type &b) {
            for (std::size_t i = 0; i < lanes; ++i)
                a[i] &= b[i];
            return a;
        }
        static value_type shift52(value_type a) {
            for (auto &x : a)
                x >>= 52;
            return a;
        }
        static mask_type less(const value_type &a, const value_type &b) {
            mask_type r;
            for (std::size_t i = 0; i < lanes; ++i)
                r[i] = a[i] < b[i];
            return r;
        }
        static value_type select(const mask_type &mask, value_type yes, const value_type &no) {
            for (std::size_t i = 0; i < lanes; ++i)
                if (!mask[i])
                    yes[i] = no[i];
            return yes;
        }
        static value_type multiply_add_low(value_type sum, const value_type &a, const value_type &b) {
            value_type product;
            for (std::size_t i = 0; i < lanes; ++i)
                product[i] = limb(wide(a[i] & mask) * (b[i] & mask)) & mask;
            return add(sum, product);
        }
        static value_type multiply_add_high(value_type sum, const value_type &a, const value_type &b) {
            value_type product;
            for (std::size_t i = 0; i < lanes; ++i)
                product[i] = limb((wide(a[i] & mask) * (b[i] & mask)) >> 52);
            return add(sum, product);
        }
    };

    template<typename Operations>
    void check_radix52_dots() {
        namespace fast = fields::detail::fp12_fast;
        using params = alt_bn128_254_fp12::extension_policy::fast_params;
        using kernel = fast::fp12_fast<params>;
        using radix = fast::montgomery_radix52<params>;
        using backend = kernel::base_backend_type;
        constexpr std::size_t lanes = Operations::lanes, terms = 1025;
        const auto &mod = params::base_field_type::modulus_params.get_mod_obj();
        auto maximum = radix::modulus();
        --maximum[0];
        boost::random::mt19937 rng(0x52F00D);
        for (bool maximal : {true, false}) {
            std::vector<radix::value_type> left(terms);
            std::array<std::vector<radix::value_type>, lanes> right;
            std::vector<fast::limb> packed_left(terms * 5), packed_right(terms * 5 * lanes);
            std::vector<std::size_t> indices(terms);
            for (auto &lane : right)
                lane.resize(terms);
            for (std::size_t j = 0; j < terms; ++j) {
                left[j] =
                    maximal ?
                        maximum :
                        radix::from_montgomery64(
                            random_fp12<alt_bn128_254_fp12>(rng).data[0].data[0].data[0].data.backend().base_data());
                indices[j] = terms - j - 1;
                for (std::size_t lane = 0; lane < lanes; ++lane) {
                    right[lane][j] = lane == 0 ? radix::value_type {} :
                                     maximal   ? maximum :
                                                 radix::from_montgomery64(random_fp12<alt_bn128_254_fp12>(rng)
                                                                              .data[0]
                                                                              .data[0]
                                                                              .data[0]
                                                                              .data.backend()
                                                                              .base_data());
                    for (std::size_t i = 0; i < 5; ++i) {
                        packed_left[j * 5 + i] = left[j][i];
                        packed_right[(indices[j] * 5 + i) * lanes + lane] = right[lane][j][i];
                    }
                }
            }
            for (std::size_t count : {0, 1, 2, 7, 8, 9, 27, 28, 123, 255, 256, 257, 511, 512, 513, 1025}) {
                std::array<fast::limb, 5 * lanes + 2> actual;
                actual.fill(0xAABBCCDDEEFF0011);
                fast::montgomery_inner_product52<params, Operations>(actual.data() + 1, packed_left.data(),
                                                                     packed_right.data(), indices.data(), count);
                BOOST_CHECK_EQUAL(actual.front(), 0xAABBCCDDEEFF0011);
                BOOST_CHECK_EQUAL(actual.back(), 0xAABBCCDDEEFF0011);
                for (std::size_t lane = 0; lane < lanes; ++lane) {
                    radix::value_type result;
                    for (std::size_t i = 0; i < 5; ++i)
                        result[i] = actual[1 + i * lanes + lane];
                    const auto expected = radix::inner_product(left.begin(), left.begin() + count, right[lane].begin());
                    BOOST_CHECK(result == expected);
                    std::vector<backend> l(count), r(count);
                    for (std::size_t j = 0; j < count; ++j) {
                        radix::to_montgomery64(l[j], left[j]);
                        radix::to_montgomery64(r[j], right[lane][j]);
                    }
                    backend original, converted;
                    mod.montgomery_inner_product(original, l.begin(), l.end(), r.begin());
                    radix::to_montgomery64(converted, result);
                    BOOST_CHECK(eval_eq(original, converted));
                    BOOST_CHECK(eval_lt(converted, mod.get_mod()));
                }
            }
        }
    }

    template<typename Operations>
    void check_radix52_add_subtract() {
        namespace fast = fields::detail::fp12_fast;
        using params = alt_bn128_254_fp12::extension_policy::fast_params;
        using radix = fast::montgomery_radix52<params>;
        using arithmetic = fast::montgomery_radix52_simd<params, Operations>;
        using boost::multiprecision::cpp_int;
        constexpr std::size_t lanes = Operations::lanes;
        cpp_int modulus = 0;
        for (auto digit = radix::modulus().rbegin(); digit != radix::modulus().rend(); ++digit) {
            modulus = (modulus << 52) + *digit;
        }
        std::vector<cpp_int> values {0, 1, 2, modulus - 1, modulus - 2};
        for (unsigned digit = 1; digit < 5; ++digit) {
            const cpp_int boundary = cpp_int(1) << (52 * digit);
            values.push_back(boundary - 1);
            values.push_back(boundary);
            values.push_back(boundary + 1);
        }
        // Cross every radix carry/borrow boundary and modulus correction in distinct lanes.
        for (std::size_t a = 0; a < values.size(); ++a) {
            for (std::size_t b = 0; b < values.size(); ++b) {
                std::array<fast::limb, 5 * lanes> left, right, sum, difference;
                for (std::size_t lane = 0; lane < lanes; ++lane) {
                    for (std::size_t digit = 0; digit < 5; ++digit) {
                        left[digit * lanes + lane] =
                            ((values[(a + lane) % values.size()] >> (52 * digit)) & radix::mask)
                                .template convert_to<fast::limb>();
                        right[digit * lanes + lane] =
                            ((values[(b + 3 * lane) % values.size()] >> (52 * digit)) & radix::mask)
                                .template convert_to<fast::limb>();
                    }
                }
                arithmetic::store(sum.data(),
                                  arithmetic::add(arithmetic::load(left.data()), arithmetic::load(right.data())));
                arithmetic::store(difference.data(),
                                  arithmetic::subtract(arithmetic::load(left.data()), arithmetic::load(right.data())));
                for (std::size_t lane = 0; lane < lanes; ++lane) {
                    cpp_int actual_sum = 0, actual_difference = 0;
                    for (std::size_t digit = 5; digit-- > 0;) {
                        BOOST_CHECK_LE(sum[digit * lanes + lane], radix::mask);
                        BOOST_CHECK_LE(difference[digit * lanes + lane], radix::mask);
                        actual_sum = (actual_sum << 52) + sum[digit * lanes + lane];
                        actual_difference = (actual_difference << 52) + difference[digit * lanes + lane];
                    }
                    const auto &x = values[(a + lane) % values.size()];
                    const auto &y = values[(b + 3 * lane) % values.size()];
                    BOOST_CHECK(actual_sum == (x + y) % modulus);
                    BOOST_CHECK(actual_difference == (x - y + modulus) % modulus);
                }
            }
        }
    }

    template<typename Operations>
    void check_radix52_recombination() {
        namespace fast = fields::detail::fp12_fast;
        using value = alt_bn128_254_fp12::value_type;
        using params = alt_bn128_254_fp12::extension_policy::fast_params;
        using kernel = fast::fp12_fast<params>;
        using radix = fast::montgomery_radix52<params>;
        constexpr std::size_t lanes = Operations::lanes, coordinates = kernel::product_coordinate_count;
        auto maximum = radix::modulus();
        --maximum[0];
        boost::random::mt19937 rng(0x521254);
        for (std::size_t pattern = 0; pattern < 58 + random_samples; ++pattern) {
            std::array<fast::limb, coordinates * 5 * lanes> products;
            std::array<std::array<typename kernel::base_value_type, coordinates>, lanes> scalar;
            for (std::size_t lane = 0; lane < lanes; ++lane) {
                for (std::size_t c = 0; c < coordinates; ++c) {
                    typename radix::value_type packed {};
                    if (pattern == 1 || (pattern == 2 && (c + lane) % 2)) {
                        packed = maximum;
                    } else if (pattern == 3 || pattern >= 58) {
                        packed = radix::from_montgomery64(random_fp12<alt_bn128_254_fp12>(rng)
                                                              .coordinate((c + lane) % 12)
                                                              .data.backend()
                                                              .base_data());
                    } else if (pattern >= 4 && c == pattern - 4) {
                        packed = radix::from_montgomery64(
                            typename kernel::base_value_type(lane + 1).data.backend().base_data());
                    }
                    radix::to_montgomery64(scalar[lane][c].data.backend().base_data(), packed);
                    for (std::size_t digit = 0; digit < 5; ++digit) {
                        products[(c * 5 + digit) * lanes + lane] = packed[digit];
                    }
                }
            }
            std::array<fast::limb, 12 * 5 * lanes + 2> combined;
            combined.fill(0xAABBCCDDEEFF0011);
            kernel::template recombine_products52<Operations>(combined.data() + 1, products.data());
            BOOST_CHECK_EQUAL(combined.front(), 0xAABBCCDDEEFF0011);
            BOOST_CHECK_EQUAL(combined.back(), 0xAABBCCDDEEFF0011);
            for (std::size_t lane = 0; lane < lanes; ++lane) {
                value actual;
                for (std::size_t c = 0; c < 12; ++c) {
                    typename radix::value_type packed;
                    for (std::size_t digit = 0; digit < 5; ++digit) {
                        packed[digit] = combined[1 + (c * 5 + digit) * lanes + lane];
                        BOOST_CHECK_LE(packed[digit], radix::mask);
                    }
                    BOOST_CHECK(std::lexicographical_compare(packed.rbegin(), packed.rend(), radix::modulus().rbegin(),
                                                             radix::modulus().rend()));
                    radix::to_montgomery64(actual.coordinate(c).data.backend().base_data(), packed);
                }
                BOOST_CHECK_EQUAL(actual, kernel::template recombine_products<value>(scalar[lane]));
            }
        }
    }

    template<typename Operations>
    void check_rank23_recombination() {
        namespace fast = fields::detail::fp12_fast;
        using field = alt_bn128_254_fp12;
        using value = field::value_type;
        using params = field::extension_policy::fast_params;
        using kernel = typename fast::fp12_fast<params>::rank23;
        using base = typename fast::fp12_fast<params>::base_value_type;
        using radix = fast::montgomery_radix52<params>;
        constexpr std::size_t lanes = Operations::lanes, count = kernel::product_coordinate_count;
        static_assert(count == 23);
        // Independent interpolation oracle: evaluate each Lagrange polynomial directly at
        // the tower element w using generic Fp12 multiplication. No flat-basis reduction or
        // production interpolation weights are used to construct the expected result.
        std::array<base, count> points;
        points[0] = base::zero();
        for (std::size_t i = 1; i < count; ++i) {
            points[i] = base((i + 1) / 2);
            if (i % 2 == 0) {
                points[i] = -points[i];
            }
        }
        std::array<value, count> weights;
        for (std::size_t i = 0; i < count; ++i) {
            weights[i] = value::one();
            auto denominator = base::one();
            for (std::size_t j = 0; j < count; ++j) {
                if (i != j) {
                    auto factor = value::zero();
                    factor.coordinate(6) = base::one();    // w
                    factor.coordinate(0) = -points[j];
                    weights[i] = generic_fp12_mul<field>(weights[i], factor);
                    denominator *= points[i] - points[j];
                }
            }
            auto inverse = value::zero();
            inverse.coordinate(0) = denominator.inversed();
            weights[i] = generic_fp12_mul<field>(weights[i], inverse);
        }
        auto maximum = radix::modulus();
        --maximum[0];
        boost::random::mt19937 rng(0x2312);
        // Every interpolation basis vector, maximum raw residues, all 144 tower basis
        // products, boundary pairs, and random products; each lane carries independent data.
        for (std::size_t pattern = 0; pattern < 26 + 144 + 36 + random_samples; ++pattern) {
            std::array<fast::limb, count * 5 * lanes> products {};
            std::array<value, lanes> expected;
            expected.fill(value::zero());
            for (std::size_t lane = 0; lane < lanes; ++lane) {
                std::array<base, count> scalar;
                if (pattern < 26) {
                    for (std::size_t i = 0; i < count; ++i) {
                        if (pattern == 24 || (pattern == 25 && (i + lane) % 2)) {
                            radix::to_montgomery64(scalar[i].data.backend().base_data(), maximum);
                        } else {
                            scalar[i] = i == pattern ? base(lane + 1) : base::zero();
                        }
                        auto factor = value::zero();
                        factor.coordinate(0) = scalar[i];
                        expected[lane] += generic_fp12_mul<field>(weights[i], factor);
                    }
                } else {
                    value left, right;
                    if (pattern < 26 + 144) {
                        left = value::zero();
                        right = value::zero();
                        left.coordinate((pattern - 26) / 12) = base(lane + 1);
                        right.coordinate((pattern - 26) % 12) = base(lane + 2);
                    } else if (pattern < 26 + 144 + 36) {
                        left = boundary_fp12<field>((pattern - 26 - 144) / 6);
                        right = boundary_fp12<field>((pattern - 26 - 144) % 6);
                    } else {
                        left = random_fp12<field>(rng);
                        right = random_fp12<field>(rng);
                    }
                    const auto a = kernel::prepare_product_operand(left);
                    const auto b = kernel::prepare_product_operand(right);
                    for (std::size_t i = 0; i < count; ++i) {
                        const auto &modulus = base::modulus_params.get_mod_obj().get_mod();
                        BOOST_CHECK(eval_lt(a[i], modulus));
                        BOOST_CHECK(eval_lt(b[i], modulus));
                        base x, y;
                        x.data.backend().base_data() = a[i];
                        y.data.backend().base_data() = b[i];
                        scalar[i] = x * y;
                    }
                    expected[lane] = generic_fp12_mul<field>(left, right);
                }
                for (std::size_t i = 0; i < count; ++i) {
                    const auto digits = radix::from_montgomery64(scalar[i].data.backend().base_data());
                    for (std::size_t digit = 0; digit < 5; ++digit) {
                        products[(i * 5 + digit) * lanes + lane] = digits[digit];
                    }
                }
            }
            std::array<fast::limb, 12 * 5 * lanes + 2> combined;
            combined.fill(0xAABBCCDDEEFF0011);
            kernel::template recombine_products52<Operations>(combined.data() + 1, products.data());
            BOOST_CHECK_EQUAL(combined.front(), 0xAABBCCDDEEFF0011);
            BOOST_CHECK_EQUAL(combined.back(), 0xAABBCCDDEEFF0011);
            for (std::size_t lane = 0; lane < lanes; ++lane) {
                value actual;
                for (std::size_t c = 0; c < 12; ++c) {
                    typename radix::value_type packed;
                    for (std::size_t digit = 0; digit < 5; ++digit) {
                        packed[digit] = combined[1 + (c * 5 + digit) * lanes + lane];
                        BOOST_CHECK_LE(packed[digit], radix::mask);
                    }
                    BOOST_CHECK(std::lexicographical_compare(packed.rbegin(), packed.rend(), radix::modulus().rbegin(),
                                                             radix::modulus().rend()));
                    radix::to_montgomery64(actual.coordinate(c).data.backend().base_data(), packed);
                }
                BOOST_CHECK_EQUAL(actual, expected[lane]);
            }
        }
    }

    template<typename Operations, bool Rank23 = false>
    void check_radix52_fp12_dots() {
        namespace fast = fields::detail::fp12_fast;
        using value = alt_bn128_254_fp12::value_type;
        using params = alt_bn128_254_fp12::extension_policy::fast_params;
        using tower_kernel = fast::fp12_fast<params>;
        using kernel = std::conditional_t<Rank23, typename tower_kernel::rank23, tower_kernel>;
        using radix = fast::montgomery_radix52<params>;
        constexpr std::size_t lanes = Operations::lanes, terms = 514, coordinates = kernel::product_coordinate_count;
        constexpr std::size_t panel_stride = 2 * terms * 5 * lanes;
        boost::random::mt19937 rng(0x5212D07);
        std::vector<value> left(terms);
        std::vector<std::array<value, lanes>> right(terms);
        std::vector<fast::limb> packed_left(coordinates * terms * 5), packed_right(coordinates * panel_stride, 0);
        std::array<std::size_t, terms> indices;
        for (std::size_t term = 0; term < terms; ++term) {
            // The second half cancels the first across the 256-term reduction boundaries.
            left[term] = term >= 257   ? -left[term - 257] :
                         term % 3 == 0 ? boundary_fp12<alt_bn128_254_fp12>(term % 6) :
                                         random_fp12<alt_bn128_254_fp12>(rng);
            indices[term] = 2 * (terms - 1 - term);
            const auto prepared_left = kernel::prepare_product_operand(left[term]);
            for (std::size_t c = 0; c < coordinates; ++c) {
                const auto digits = radix::from_montgomery64(prepared_left[c]);
                std::copy(digits.begin(), digits.end(), packed_left.begin() + (c * terms + term) * 5);
            }
            for (std::size_t lane = 0; lane < lanes; ++lane) {
                // Zero first/last lanes also model empty coefficients and a padded panel tail.
                right[term][lane] = lane == 0 || lane + 1 == lanes ? value::zero() :
                                    term >= 257                    ? right[term - 257][lane] :
                                                                     random_fp12<alt_bn128_254_fp12>(rng);
                const auto prepared_right = kernel::prepare_product_operand(right[term][lane]);
                for (std::size_t c = 0; c < coordinates; ++c) {
                    const auto digits = radix::from_montgomery64(prepared_right[c]);
                    for (std::size_t digit = 0; digit < 5; ++digit) {
                        packed_right[c * panel_stride + (indices[term] * 5 + digit) * lanes + lane] = digits[digit];
                    }
                }
            }
        }
        std::array<value, lanes> expected;
        expected.fill(value::zero());
        std::size_t verified_terms = 0;
        for (std::size_t count : {0, 1, 123, 255, 256, 257, 513, 514}) {
            for (; verified_terms < count; ++verified_terms) {
                for (std::size_t lane = 0; lane < lanes; ++lane) {
                    expected[lane] +=
                        generic_fp12_mul<alt_bn128_254_fp12>(left[verified_terms], right[verified_terms][lane]);
                }
            }
            std::array<fast::limb, coordinates * 5 * lanes> products;
            std::array<fast::limb, 12 * 5 * lanes> combined;
            for (std::size_t c = 0; c < coordinates; ++c) {
                fast::montgomery_inner_product52<params, Operations>(
                    products.data() + c * 5 * lanes, packed_left.data() + c * terms * 5,
                    packed_right.data() + c * panel_stride, indices.data(), count);
            }
            kernel::template recombine_products52<Operations>(combined.data(), products.data());
            for (std::size_t lane = 0; lane < lanes; ++lane) {
                value actual;
                for (std::size_t c = 0; c < 12; ++c) {
                    typename radix::value_type packed;
                    for (std::size_t digit = 0; digit < 5; ++digit) {
                        packed[digit] = combined[(c * 5 + digit) * lanes + lane];
                        BOOST_CHECK_LE(packed[digit], radix::mask);
                    }
                    BOOST_CHECK(std::lexicographical_compare(packed.rbegin(), packed.rend(), radix::modulus().rbegin(),
                                                             radix::modulus().rend()));
                    radix::to_montgomery64(actual.coordinate(c).data.backend().base_data(), packed);
                }
                BOOST_CHECK_EQUAL(actual, expected[lane]);
                if (count == terms) {
                    BOOST_CHECK_EQUAL(actual, value::zero());
                }
            }
        }
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(fp12_tests)

BOOST_AUTO_TEST_CASE_TEMPLATE(tower_values_embed_into_fp12, Fp12Field, fp12_field_types) {
    using fp12_value_type = typename Fp12Field::value_type;
    using fp6_value_type = typename fp12_value_type::underlying_type;
    using fp2_value_type = typename fp6_value_type::underlying_type;
    using fp_value_type = typename fp2_value_type::underlying_type;

    const fp_value_type fp(7);
    const fp2_value_type fp2(fp);
    const fp6_value_type fp6(fp2);

    BOOST_CHECK_EQUAL(fp2, fp2_value_type(fp, fp_value_type::zero()));
    BOOST_CHECK_EQUAL(fp6, fp6_value_type(fp2, fp2_value_type::zero(), fp2_value_type::zero()));
    BOOST_CHECK_EQUAL(fp12_value_type(fp), fp12_value_type(fp6));
    BOOST_CHECK_EQUAL(fp12_value_type(fp2), fp12_value_type(fp6));
    BOOST_CHECK_EQUAL(fp12_value_type(fp6), fp12_value_type(fp6, fp6_value_type::zero()));
    BOOST_CHECK_EQUAL(fp * fp2_value_type::one(), fp2_value_type(fp));
    BOOST_CHECK_EQUAL(fp2_value_type::one() * fp, fp2_value_type(fp));
    BOOST_CHECK_EQUAL(fp * fp6_value_type::one(), fp6_value_type(fp));
    BOOST_CHECK_EQUAL(fp6_value_type::one() * fp, fp6_value_type(fp));
    BOOST_CHECK_EQUAL(fp2 * fp6_value_type::one(), fp6_value_type(fp2));
    BOOST_CHECK_EQUAL(fp6_value_type::one() * fp2, fp6_value_type(fp2));
    BOOST_CHECK_EQUAL(fp * fp12_value_type::one(), fp12_value_type(fp));
    BOOST_CHECK_EQUAL(fp12_value_type::one() * fp, fp12_value_type(fp));
    BOOST_CHECK_EQUAL(fp2 * fp12_value_type::one(), fp12_value_type(fp2));
    BOOST_CHECK_EQUAL(fp12_value_type::one() * fp2, fp12_value_type(fp2));
    BOOST_CHECK_EQUAL(fp6 * fp12_value_type::one(), fp12_value_type(fp6));
    BOOST_CHECK_EQUAL(fp12_value_type::one() * fp6, fp12_value_type(fp6));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(square_roots_round_trip, Fp12Field, fp12_field_types) {
    using fp12_value_type = typename Fp12Field::value_type;
    boost::random::mt19937 rng(0x5a127);

    BOOST_CHECK(fp12_value_type::zero().is_square());
    BOOST_CHECK_EQUAL(fp12_value_type::zero().sqrt(), fp12_value_type::zero());
    BOOST_CHECK(fp12_value_type::one().is_square());
    BOOST_CHECK_EQUAL(fp12_value_type::one().sqrt(), fp12_value_type::one());

    for (std::size_t i = 0; i < 2; ++i) {
        const fp12_value_type square = random_fp12<Fp12Field>(rng).squared();
        BOOST_CHECK(square.is_square());
        BOOST_CHECK_EQUAL(square.sqrt_known_square().squared(), square);
        BOOST_CHECK_EQUAL(square.sqrt().squared(), square);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(is_square_rejects_a_non_square, Fp12Field, fp12_field_types) {
    using fp12_value_type = typename Fp12Field::value_type;
    boost::random::mt19937 rng(0x1105);

    for (std::size_t i = 0; i < 64; ++i) {
        const fp12_value_type candidate = random_fp12<Fp12Field>(rng);
        if (!candidate.is_square()) {
            BOOST_CHECK(!candidate.is_zero());
            return;
        }
    }
    BOOST_FAIL("failed to sample an Fp12 non-square");
}

BOOST_AUTO_TEST_CASE(bn254_subgroup_algorithms_match_direct_exponentiation) {
    using field_type = alt_bn128_254_fp12;
    using value_type = typename field_type::value_type;
    using boost::multiprecision::cpp_int;

    const auto order_decomposition = fields::field_multiplicative_group_decomposition<field_type>();
    const cpp_int &odd_order = order_decomposition.odd_order;
    const std::size_t two_adicity = order_decomposition.two_adicity;
    boost::random::mt19937 rng(0x2a11);

    for (std::size_t i = 0; i < 2; ++i) {
        const value_type x = random_fp12<field_type>(rng);
        BOOST_CHECK_EQUAL(fields::two_primary_component(x, odd_order, two_adicity), x.pow(odd_order));
        BOOST_CHECK_EQUAL(fields::two_primary_component(x, odd_order + 2, two_adicity), x.pow(odd_order + 2));

        const value_type odd_subgroup_value = x.pow(cpp_int(1) << two_adicity);
        const value_type root = fields::odd_subgroup_sqrt(odd_subgroup_value, odd_order);
        BOOST_CHECK_EQUAL(root.squared(), odd_subgroup_value);
        BOOST_CHECK_EQUAL(root, odd_subgroup_value.pow((odd_order + 1) >> 1));

        const cpp_int fallback_order = odd_order + 2;
        BOOST_CHECK_EQUAL(fields::odd_subgroup_sqrt(x, fallback_order), x.pow((fallback_order + 1) >> 1));
    }
}

BOOST_AUTO_TEST_CASE(bn254_roots_of_unity_have_order_32) {
    using value_type = typename alt_bn128_254_fp12::value_type;
    const value_type primitive_root = fields::primitive_two_power_root_of_unity<value_type>(5);
    const std::vector<value_type> roots = fields::roots_of_unity<value_type>(5);
    const std::vector<value_type> roots_from_supplied_root = fields::roots_of_unity(primitive_root, 5);

    BOOST_REQUIRE_EQUAL(roots.size(), 32);
    BOOST_CHECK_EQUAL_COLLECTIONS(roots.begin(), roots.end(), roots_from_supplied_root.begin(),
                                  roots_from_supplied_root.end());
    BOOST_CHECK_EQUAL(roots.front(), value_type::one());
    BOOST_CHECK_EQUAL(primitive_root.pow(32), value_type::one());
    BOOST_CHECK_NE(primitive_root.pow(16), value_type::one());
    for (std::size_t i = 0; i < roots.size(); ++i) {
        BOOST_CHECK_EQUAL(roots[i], primitive_root.pow(i));
    }
}

BOOST_AUTO_TEST_CASE(root_of_unity_argument_validation) {
    using bn254_value_type = typename alt_bn128_254_fp12::value_type;
    using bls12_value_type = typename bls12_381_fp12::value_type;

    BOOST_CHECK_THROW(fields::primitive_two_power_root_of_unity<bn254_value_type>(4), std::invalid_argument);
    BOOST_CHECK_THROW(fields::primitive_two_power_root_of_unity<bls12_value_type>(5), std::invalid_argument);
    BOOST_CHECK_THROW(fields::roots_of_unity(bn254_value_type::one(), std::numeric_limits<std::size_t>::digits),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(generic_subgroup_algorithms_match_direct_formulas) {
    using field_type = bls12_381_fp12;
    using value_type = typename field_type::value_type;
    using boost::multiprecision::cpp_int;

    const auto order_decomposition = fields::field_multiplicative_group_decomposition<field_type>();
    const cpp_int &odd_order = order_decomposition.odd_order;
    const std::size_t two_adicity = order_decomposition.two_adicity;
    boost::random::mt19937 rng(0x3812);
    const value_type x = random_fp12<field_type>(rng);
    BOOST_CHECK_EQUAL(fields::two_primary_component(x, odd_order, two_adicity), x.pow(odd_order));

    const value_type odd_subgroup_value = x.pow(cpp_int(1) << two_adicity);
    BOOST_CHECK_EQUAL(fields::odd_subgroup_sqrt(odd_subgroup_value, odd_order),
                      odd_subgroup_value.pow((odd_order + 1) >> 1));

    const value_type square = x.squared();
    BOOST_CHECK_EQUAL(fields::sqrt_known_square(square).squared(), square);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(coordinates_follow_tower_storage_order, Fp12Field, fp12_field_types) {
    using fp12_value_type = typename Fp12Field::value_type;
    using base_value_type = typename Fp12Field::base_field_type::value_type;

    fp12_value_type value = fp12_value_type::zero();
    for (std::size_t index = 0; index < Fp12Field::arity; ++index) {
        const std::size_t outer_index = index / 6;
        const std::size_t middle_index = (index % 6) / 2;
        const std::size_t inner_index = index % 2;
        value.data[outer_index].data[middle_index].data[inner_index] = base_value_type(index + 1);
    }

    const fp12_value_type &const_value = value;
    for (std::size_t index = 0; index < Fp12Field::arity; ++index) {
        const std::size_t outer_index = index / 6;
        const std::size_t middle_index = (index % 6) / 2;
        const std::size_t inner_index = index % 2;
        const base_value_type &nested_coordinate = const_value.data[outer_index].data[middle_index].data[inner_index];
        BOOST_CHECK(std::addressof(const_value.coordinate(index)) == std::addressof(nested_coordinate));
    }

    for (std::size_t index = 0; index < Fp12Field::arity; ++index) {
        const std::size_t outer_index = index / 6;
        const std::size_t middle_index = (index % 6) / 2;
        const std::size_t inner_index = index % 2;
        const base_value_type replacement(index + Fp12Field::arity + 1);
        value.coordinate(index) = replacement;
        BOOST_CHECK_EQUAL(value.data[outer_index].data[middle_index].data[inner_index], replacement);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multiplication_matches_generic_implementation, Fp12Field, fp12_field_types) {
    using fp12_value_type = typename Fp12Field::value_type;
    boost::random::mt19937 rng(0x25412);

    BOOST_TEST_CONTEXT(fp12_field_name<Fp12Field>::value) {
        for (std::size_t i = 0; i < random_samples; ++i) {
            const fp12_value_type x = random_fp12<Fp12Field>(rng);
            const fp12_value_type y = random_fp12<Fp12Field>(rng);

            BOOST_CHECK_EQUAL(x * y, generic_fp12_mul<Fp12Field>(x, y));

            fp12_value_type inplace = x;
            inplace *= y;
            BOOST_CHECK_EQUAL(inplace, generic_fp12_mul<Fp12Field>(x, y));
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(prepared_product_forms_recombine_exactly, Fp12Field, fp12_field_types) {
    using kernel = fields::detail::fp12_fast::fp12_fast<typename Fp12Field::extension_policy::fast_params>;
    using base_value_type = typename kernel::base_value_type;
    boost::random::mt19937 rng(0x5412);
    const auto check = [&](const auto &x, const auto &y) {
        const auto left = kernel::prepare_product_operand(x);
        const auto right = kernel::prepare_product_operand(y);
        std::array<base_value_type, kernel::product_coordinate_count> products;
        for (std::size_t i = 0; i < products.size(); ++i) {
            const auto &modulus = base_value_type::modulus_params.get_mod_obj().get_mod();
            BOOST_CHECK(eval_lt(left[i], modulus));
            BOOST_CHECK(eval_lt(right[i], modulus));
            base_value_type a, b;
            a.data.backend().base_data() = left[i];
            b.data.backend().base_data() = right[i];
            products[i] = a * b;
        }
        BOOST_CHECK_EQUAL(kernel::template recombine_products<typename Fp12Field::value_type>(products),
                          generic_fp12_mul<Fp12Field>(x, y));
    };
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            check(boundary_fp12<Fp12Field>(i), boundary_fp12<Fp12Field>(j));
        }
    }
    for (std::size_t i = 0; i < random_samples; ++i) {
        check(random_fp12<Fp12Field>(rng), random_fp12<Fp12Field>(rng));
    }
    check(Fp12Field::value_type::zero(), Fp12Field::value_type::one());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(prepared_inner_products_match_generic_arithmetic, Fp12Field, fp12_field_types) {
    using value_type = typename Fp12Field::value_type;
    using kernel = fields::detail::fp12_fast::fp12_fast<typename Fp12Field::extension_policy::fast_params>;
    using base_value_type = typename kernel::base_value_type;
    using prepared_type = typename kernel::prepared_product_operand;
    boost::random::mt19937 rng(0x12354);
    std::vector<value_type> left, right;
    std::vector<prepared_type> prepared_left, prepared_right;
    value_type maximum_residues;
    auto maximum = base_value_type::modulus_params.get_mod_obj().get_mod();
    typename kernel::base_backend_type one;
    one.limbs()[0] = 1;
    eval_subtract(maximum, one);
    for (std::size_t i = 0; i < Fp12Field::arity; ++i) {
        maximum_residues.coordinate(i).data.backend().base_data() = maximum;
    }
    for (std::size_t i = 0; i < 513; ++i) {
        left.push_back(i % 3 == 0 ? maximum_residues : random_fp12<Fp12Field>(rng));
        right.push_back(i % 3 == 0 ? maximum_residues : boundary_fp12<Fp12Field>(i % 6));
        prepared_left.push_back(kernel::prepare_product_operand(left.back()));
        prepared_right.push_back(kernel::prepare_product_operand(right.back()));
    }
    for (std::size_t count : {0, 1, 2, 3, 4, 5, 7, 8, 9, 16, 31, 32, 33, 123, 124, 256, 257, 513}) {
        value_type expected = value_type::zero();
        for (std::size_t i = 0; i < count; ++i) {
            expected += generic_fp12_mul<Fp12Field>(left[i], right[i]);
        }
        const auto actual = kernel::template inner_product<value_type>(
            prepared_left.begin(), prepared_left.begin() + count, prepared_right.begin());
        BOOST_CHECK_EQUAL(actual, expected);
    }

    // Cancellation across a long dot must remain canonical, including the empty result.
    for (std::size_t i = 0; i < 246; i += 2) {
        prepared_left[i + 1] = kernel::prepare_product_operand(-left[i]);
        prepared_right[i + 1] = prepared_right[i];
    }
    BOOST_CHECK_EQUAL(kernel::template inner_product<value_type>(prepared_left.begin(), prepared_left.begin() + 246,
                                                                 prepared_right.begin()),
                      value_type::zero());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(prepared_base_dots_match_portable_with_wide_carries, Fp12Field, fp12_field_types) {
    using params = typename Fp12Field::extension_policy::fast_params;
    using kernel = fields::detail::fp12_fast::fp12_fast<params>;
    using backend = typename kernel::base_backend_type;
    const auto &modulus = params::base_field_type::modulus_params.get_mod_obj();
    backend maximum = modulus.get_mod(), one;
    one.limbs()[0] = 1;
    eval_subtract(maximum, one);
    std::vector<backend> left(1025, maximum), right(1025, maximum);
    // For BN254, 28 products of (p-1)^2 first carry beyond the eight ordinary product limbs.
    // Also exercise several REDC canonicalization subtractions and the six-limb portable fallback.
    for (std::size_t count : {0, 1, 2, 3, 27, 28, 29, 31, 32, 33, 122, 123, 124, 255, 256, 257, 513, 1025}) {
        backend actual, expected;
        modulus.montgomery_inner_product(expected, left.begin(), left.begin() + count, right.begin());
        fields::detail::fp12_fast::montgomery_inner_product<params>(actual, left.begin(), left.begin() + count,
                                                                    right.begin());
        BOOST_CHECK(eval_eq(actual, expected));
        BOOST_CHECK(eval_lt(actual, modulus.get_mod()));
        // Writing the result must wait until all input terms have been consumed.
        fields::detail::fp12_fast::montgomery_inner_product<params>(left[0], left.begin(), left.begin() + count,
                                                                    right.begin());
        BOOST_CHECK(eval_eq(left[0], expected));
        left[0] = maximum;
    }
}

BOOST_AUTO_TEST_CASE(radix52_conversion_preserves_canonical_montgomery_values) {
    namespace fast = fields::detail::fp12_fast;
    using params = alt_bn128_254_fp12::extension_policy::fast_params;
    using radix = fast::montgomery_radix52<params>;
    using backend = fast::fp12_fast<params>::base_backend_type;
    using boost::multiprecision::cpp_int;
    const auto integer = [](const auto &x, unsigned bits) {
        cpp_int result = 0;
        for (auto i = x.rbegin(); i != x.rend(); ++i) {
            result <<= bits;
            result += *i;
        }
        return result;
    };
    const auto &mod = params::base_field_type::modulus_params.get_mod_obj();
    const cpp_int p = integer(fast::load_limbs(mod.get_mod()), 64);
    boost::random::mt19937 rng(0x5264);
    for (std::size_t i = 0; i < 300; ++i) {
        backend original = i < 16 ?
                               mod.get_mod() :
                               random_fp12<alt_bn128_254_fp12>(rng).data[0].data[0].data[0].data.backend().base_data();
        if (i < 16) {
            backend delta;
            delta.limbs()[0] = i + 1;
            eval_subtract(original, delta);
        }
        if (i == 16)
            original = backend();
        const auto packed = radix::from_montgomery64(original);
        BOOST_CHECK(integer(packed, 52) == (integer(fast::load_limbs(original), 64) * 16) % p);
        for (auto digit : packed)
            BOOST_CHECK_LE(digit, radix::mask);
        backend restored;
        radix::to_montgomery64(restored, packed);
        BOOST_CHECK(eval_eq(restored, original));
    }
}

BOOST_AUTO_TEST_CASE(radix52_four_and_eight_lane_arithmetic_matches_reference) {
    check_radix52_dots<scalar_ifma52_operations<4>>();
    check_radix52_dots<scalar_ifma52_operations<8>>();
}

BOOST_AUTO_TEST_CASE(radix52_rank23_matches_generic_fp12_in_four_and_eight_lanes) {
    check_rank23_recombination<scalar_ifma52_operations<4>>();
    check_rank23_recombination<scalar_ifma52_operations<8>>();
    check_radix52_fp12_dots<scalar_ifma52_operations<4>, true>();
    check_radix52_fp12_dots<scalar_ifma52_operations<8>, true>();
}

BOOST_AUTO_TEST_CASE(radix52_recombination_matches_scalar_tower_in_four_and_eight_lanes) {
    check_radix52_add_subtract<scalar_ifma52_operations<4>>();
    check_radix52_add_subtract<scalar_ifma52_operations<8>>();
    check_radix52_recombination<scalar_ifma52_operations<4>>();
    check_radix52_recombination<scalar_ifma52_operations<8>>();
    check_radix52_fp12_dots<scalar_ifma52_operations<4>>();
    check_radix52_fp12_dots<scalar_ifma52_operations<8>>();
}

#if defined(CRYPTO3_HAS_AVX512IFMA)
BOOST_AUTO_TEST_CASE(radix52_ifma_rank23_matches_generic_fp12) {
    check_rank23_recombination<fields::detail::fp12_fast::ifma52_operations<4>>();
    check_rank23_recombination<fields::detail::fp12_fast::ifma52_operations<8>>();
    check_radix52_fp12_dots<fields::detail::fp12_fast::ifma52_operations<4>, true>();
    check_radix52_fp12_dots<fields::detail::fp12_fast::ifma52_operations<8>, true>();
}

BOOST_AUTO_TEST_CASE(radix52_ifma_recombination_matches_generic_fp12) {
    check_radix52_add_subtract<fields::detail::fp12_fast::ifma52_operations<4>>();
    check_radix52_add_subtract<fields::detail::fp12_fast::ifma52_operations<8>>();
    check_radix52_recombination<fields::detail::fp12_fast::ifma52_operations<4>>();
    check_radix52_recombination<fields::detail::fp12_fast::ifma52_operations<8>>();
    check_radix52_fp12_dots<fields::detail::fp12_fast::ifma52_operations<4>>();
    check_radix52_fp12_dots<fields::detail::fp12_fast::ifma52_operations<8>>();
}

BOOST_AUTO_TEST_CASE(radix52_ifma_matches_reference_in_every_lane) {
    check_radix52_dots<fields::detail::fp12_fast::ifma52_operations<4>>();
    check_radix52_dots<fields::detail::fp12_fast::ifma52_operations<8>>();
}
#endif

#if defined(CRYPTO3_HAS_BMI2) && defined(CRYPTO3_HAS_ADX)
BOOST_AUTO_TEST_CASE(four_limb_multiply_accumulate_preserves_all_nine_limbs) {
    using limb = fields::detail::fp12_fast::limb;
    using boost::multiprecision::cpp_int;
    constexpr limb maximum = std::numeric_limits<limb>::max();
    const auto integer = [](const limb *data, std::size_t size) {
        cpp_int result = 0;
        while (size > 0) {
            result <<= 64;
            result += data[--size];
        }
        return result;
    };
    const auto check = [&](const std::array<limb, 4> &left, const std::array<limb, 4> &right,
                           const std::array<limb, 9> &initial, bool square = false) {
        std::array<limb, 11> storage;
        storage.front() = 0xAABBCCDDEEFF0011;
        storage.back() = 0x1122334455667788;
        std::copy(initial.begin(), initial.end(), storage.begin() + 1);
        const auto *right_data = square ? left.data() : right.data();
        const cpp_int expected = integer(initial.data(), initial.size()) +
                                 integer(left.data(), left.size()) * integer(right_data, right.size());
        fields::detail::fp12_fast::multiply_add_4x4_x86(storage.data() + 1, left.data(), right_data);
        BOOST_CHECK_EQUAL(integer(storage.data() + 1, 9), expected);
        BOOST_CHECK_EQUAL(storage.front(), 0xAABBCCDDEEFF0011);
        BOOST_CHECK_EQUAL(storage.back(), 0x1122334455667788);
    };
    const std::array<limb, 4> full {maximum, maximum, maximum, maximum}, zero {};
    std::array<limb, 9> initial;
    initial.fill(maximum);
    check(zero, full, initial);
    initial.back() = maximum - 1;    // Carry into a nearly full ninth limb without overflowing it.
    check(full, full, initial);
    check(full, full, initial, true);
    initial.back() = 0;
    check(full, full, initial);
    check(full, full, {});
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            for (unsigned bit : {0, 1, 31, 32, 63}) {
                std::array<limb, 4> left {}, right {};
                left[i] = limb(1) << bit;
                right[j] = maximum;
                check(left, right, initial);
            }
        }
    }
    boost::random::mt19937 rng(0x4ADC);
    const auto random_limb = [&]() { return (limb(rng()) << 32) | rng(); };
    for (std::size_t sample = 0; sample < 1024; ++sample) {
        std::array<limb, 4> left, right;
        for (auto &value : left)
            value = random_limb();
        for (auto &value : right)
            value = random_limb();
        for (auto &value : initial)
            value = random_limb();
        initial.back() >>= 1;    // Leave room for the product's carry.
        check(left, right, initial, sample % 2 == 0);
    }
}
#endif

BOOST_AUTO_TEST_CASE_TEMPLATE(square_matches_generic_implementation, Fp12Field, fp12_field_types) {
    using fp12_value_type = typename Fp12Field::value_type;
    boost::random::mt19937 rng(0x2545);

    BOOST_TEST_CONTEXT(fp12_field_name<Fp12Field>::value) {
        for (std::size_t i = 0; i < random_samples; ++i) {
            const fp12_value_type x = random_fp12<Fp12Field>(rng);

            BOOST_CHECK_EQUAL(x.squared(), generic_fp12_mul<Fp12Field>(x, x));
        }
    }
}

BOOST_AUTO_TEST_CASE(adversarial_bls12_377_multiplication_and_square_match_generic_implementation) {
    using fp12_field_type = bls12_377_fp12;
    using fp12_value_type = fp12_field_type::value_type;

    std::array<fp12_value_type, 6> values;
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = boundary_fp12<fp12_field_type>(i);
    }

    for (const fp12_value_type &x : values) {
        BOOST_CHECK_EQUAL(x.squared(), generic_fp12_mul<fp12_field_type>(x, x));
        for (const fp12_value_type &y : values) {
            BOOST_CHECK_EQUAL(x * y, generic_fp12_mul<fp12_field_type>(x, y));
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multiplication_identities, Fp12Field, fp12_field_types) {
    using fp12_value_type = typename Fp12Field::value_type;
    boost::random::mt19937 rng(0x2541d);

    BOOST_TEST_CONTEXT(fp12_field_name<Fp12Field>::value) {
        for (std::size_t i = 0; i < random_samples; ++i) {
            const fp12_value_type x = random_fp12<Fp12Field>(rng);

            BOOST_CHECK_EQUAL(x * fp12_value_type::zero(), fp12_value_type::zero());
            BOOST_CHECK_EQUAL(fp12_value_type::zero() * x, fp12_value_type::zero());
            BOOST_CHECK_EQUAL(x * fp12_value_type::one(), x);
            BOOST_CHECK_EQUAL(fp12_value_type::one() * x, x);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(inverse_is_multiplicative_identity, Fp12Field, fp12_field_types) {
    using fp12_value_type = typename Fp12Field::value_type;
    boost::random::mt19937 rng(0x25417);

    BOOST_TEST_CONTEXT(fp12_field_name<Fp12Field>::value) {
        for (std::size_t i = 0; i < random_samples; ++i) {
            fp12_value_type x = random_fp12<Fp12Field>(rng);
            if (x.is_zero()) {
                x = fp12_value_type::one();
            }

            const fp12_value_type x_inv = x.inversed();
            BOOST_CHECK_EQUAL(generic_fp12_mul<Fp12Field>(x, x_inv), fp12_value_type::one());
            BOOST_CHECK_EQUAL(x * x_inv, fp12_value_type::one());
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
