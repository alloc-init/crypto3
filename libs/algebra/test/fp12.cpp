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
#include <string_view>

#include <boost/mpl/list.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/algebra/random_element.hpp>

namespace {
    namespace fields = nil::crypto3::algebra::fields;

    using alt_bn128_254_fp12 = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using bls12_377_fp12 = fields::fp12_2over3over2<fields::bls12<377>>;
    using bls12_381_fp12 = fields::fp12_2over3over2<fields::bls12<381>>;

    using fp12_field_types =
        boost::mpl::list<alt_bn128_254_fp12, bls12_377_fp12, bls12_381_fp12>;

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

    constexpr std::size_t random_samples = 10000;

    // element_fp12_2over3over2::operator* uses policy_type::multiply when the
    // policy provides it; otherwise it falls back to the generic Fp12 tower
    // multiplication in the element class. This test-only policy reuses the
    // real field parameters but intentionally omits multiply, so the expected
    // value is produced by the generic implementation rather than by a
    // duplicated reference formula in this file.
    template<typename Fp12Field>
    struct generic_fp12_policy {
        using optimized_policy = typename Fp12Field::extension_policy;

        using field_type = typename optimized_policy::field_type;
        using non_residue_type = typename optimized_policy::non_residue_type;
        using underlying_type = typename optimized_policy::underlying_type;

        inline constexpr static const non_residue_type non_residue = optimized_policy::non_residue;
    };

    template<typename Fp12Field>
    using generic_fp12_value_type =
        fields::detail::element_fp12_2over3over2<generic_fp12_policy<Fp12Field>>;

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

}    // namespace

BOOST_AUTO_TEST_SUITE(fp12_tests)

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

BOOST_AUTO_TEST_CASE(adversarial_bls12_377_multiplication_matches_generic_implementation) {
    using fp12_field_type = bls12_377_fp12;
    using fp12_value_type = fp12_field_type::value_type;

    std::array<fp12_value_type, 6> values;
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = boundary_fp12<fp12_field_type>(i);
    }

    for (const fp12_value_type &x : values) {
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
