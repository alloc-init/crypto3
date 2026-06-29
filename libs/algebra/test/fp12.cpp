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
    using bn128_254_fp12 = fields::fp12_2over3over2<fields::bn128<254>>;
    using bls12_377_fp12 = fields::fp12_2over3over2<fields::bls12<377>>;
    using bls12_381_fp12 = fields::fp12_2over3over2<fields::bls12<381>>;

    using fp12_field_types =
        boost::mpl::list<alt_bn128_254_fp12, bn128_254_fp12, bls12_377_fp12, bls12_381_fp12>;

    template<typename FieldType>
    struct fp12_field_name;

    template<>
    struct fp12_field_name<alt_bn128_254_fp12> {
        constexpr static std::string_view value = "alt_bn128_254";
    };

    template<>
    struct fp12_field_name<bn128_254_fp12> {
        constexpr static std::string_view value = "bn128_254";
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

    template<typename Fp2Value>
    Fp2Value ref_fp2_add(const Fp2Value &x, const Fp2Value &y) {
        return Fp2Value(x.data[0] + y.data[0], x.data[1] + y.data[1]);
    }

    template<typename Fp2Value>
    Fp2Value ref_fp2_mul(const Fp2Value &x, const Fp2Value &y) {
        using fp_value_type = typename Fp2Value::underlying_type;

        const fp_value_type a0b0 = x.data[0] * y.data[0];
        const fp_value_type a1b1 = x.data[1] * y.data[1];
        const fp_value_type non_residue = Fp2Value::field_type::extension_policy::non_residue;

        return Fp2Value(a0b0 + non_residue * a1b1,
                        (x.data[0] + x.data[1]) * (y.data[0] + y.data[1]) - a0b0 - a1b1);
    }

    template<typename Fp12Value>
    std::array<typename Fp12Value::underlying_type::underlying_type, 6>
        to_w_coefficients(const Fp12Value &x) {
        return {x.data[0].data[0], x.data[1].data[0], x.data[0].data[1],
                x.data[1].data[1], x.data[0].data[2], x.data[1].data[2]};
    }

    template<typename Fp12Value>
    Fp12Value from_w_coefficients(
        const std::array<typename Fp12Value::underlying_type::underlying_type, 6> &x) {
        using fp6_value_type = typename Fp12Value::underlying_type;

        return Fp12Value(fp6_value_type(x[0], x[2], x[4]), fp6_value_type(x[1], x[3], x[5]));
    }

    template<typename Fp12Value>
    Fp12Value ref_fp12_mul(const Fp12Value &x, const Fp12Value &y) {
        using fp6_value_type = typename Fp12Value::underlying_type;
        using fp2_value_type = typename fp6_value_type::underlying_type;

        const auto x_coefficients = to_w_coefficients<Fp12Value>(x);
        const auto y_coefficients = to_w_coefficients<Fp12Value>(y);

        std::array<fp2_value_type, 11> product;
        product.fill(fp2_value_type::zero());

        for (std::size_t i = 0; i < x_coefficients.size(); ++i) {
            for (std::size_t j = 0; j < y_coefficients.size(); ++j) {
                product[i + j] =
                    ref_fp2_add(product[i + j], ref_fp2_mul(x_coefficients[i], y_coefficients[j]));
            }
        }

        const fp2_value_type xi = fp6_value_type::field_type::extension_policy::non_residue;
        for (std::size_t i = product.size(); i-- > 6;) {
            product[i - 6] = ref_fp2_add(product[i - 6], ref_fp2_mul(product[i], xi));
        }

        return from_w_coefficients<Fp12Value>(
            {product[0], product[1], product[2], product[3], product[4], product[5]});
    }

    template<typename Fp12Field, typename Rng>
    typename Fp12Field::value_type random_fp12(Rng &rng) {
        return nil::crypto3::algebra::random_element<Fp12Field>(rng);
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(fp12_tests)

BOOST_AUTO_TEST_CASE_TEMPLATE(multiplication_matches_independent_reference, Fp12Field, fp12_field_types) {
    using fp12_value_type = typename Fp12Field::value_type;
    boost::random::mt19937 rng(0x25412);

    BOOST_TEST_CONTEXT(fp12_field_name<Fp12Field>::value) {
        for (std::size_t i = 0; i < random_samples; ++i) {
            const fp12_value_type x = random_fp12<Fp12Field>(rng);
            const fp12_value_type y = random_fp12<Fp12Field>(rng);

            BOOST_CHECK_EQUAL(x * y, ref_fp12_mul(x, y));

            fp12_value_type inplace = x;
            inplace *= y;
            BOOST_CHECK_EQUAL(inplace, ref_fp12_mul(x, y));
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(square_matches_independent_reference, Fp12Field, fp12_field_types) {
    using fp12_value_type = typename Fp12Field::value_type;
    boost::random::mt19937 rng(0x2545);

    BOOST_TEST_CONTEXT(fp12_field_name<Fp12Field>::value) {
        for (std::size_t i = 0; i < random_samples; ++i) {
            const fp12_value_type x = random_fp12<Fp12Field>(rng);

            BOOST_CHECK_EQUAL(x.squared(), ref_fp12_mul(x, x));
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
            BOOST_CHECK_EQUAL(ref_fp12_mul(x, x_inv), fp12_value_type::one());
            BOOST_CHECK_EQUAL(x * x_inv, fp12_value_type::one());
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
