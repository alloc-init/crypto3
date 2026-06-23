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

#define BOOST_TEST_MODULE algebra_bn254_fp12_test

#include <array>
#include <cstddef>

#include <boost/random/mersenne_twister.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/random_element.hpp>

using curve_type = nil::crypto3::algebra::curves::alt_bn128<254>;
using fp12_field_type = typename curve_type::gt_type;
using fp12_value_type = typename fp12_field_type::value_type;
using fp6_value_type = typename fp12_value_type::underlying_type;
using fp2_value_type = typename fp6_value_type::underlying_type;
using fp_value_type = typename fp2_value_type::underlying_type;

namespace {

constexpr std::size_t random_samples = 32;

fp2_value_type ref_add(const fp2_value_type &x, const fp2_value_type &y) {
    return fp2_value_type(x.data[0] + y.data[0], x.data[1] + y.data[1]);
}

fp2_value_type ref_mul(const fp2_value_type &x, const fp2_value_type &y) {
    const fp_value_type a0b0 = x.data[0] * y.data[0];
    const fp_value_type a1b1 = x.data[1] * y.data[1];
    const fp_value_type non_residue = fp2_value_type::field_type::extension_policy::non_residue;

    return fp2_value_type(a0b0 + non_residue * a1b1,
                          (x.data[0] + x.data[1]) * (y.data[0] + y.data[1]) -
                              a0b0 - a1b1);
}

std::array<fp2_value_type, 6> to_w_coefficients(const fp12_value_type &x) {
    return {x.data[0].data[0], x.data[1].data[0],
            x.data[0].data[1], x.data[1].data[1],
            x.data[0].data[2], x.data[1].data[2]};
}

fp12_value_type from_w_coefficients(const std::array<fp2_value_type, 6> &x) {
    return fp12_value_type(fp6_value_type(x[0], x[2], x[4]),
                           fp6_value_type(x[1], x[3], x[5]));
}

fp12_value_type ref_mul(const fp12_value_type &x, const fp12_value_type &y) {
    const auto x_coefficients = to_w_coefficients(x);
    const auto y_coefficients = to_w_coefficients(y);

    std::array<fp2_value_type, 11> product;
    product.fill(fp2_value_type::zero());

    for (std::size_t i = 0; i < x_coefficients.size(); ++i) {
        for (std::size_t j = 0; j < y_coefficients.size(); ++j) {
            product[i + j] = ref_add(product[i + j],
                                     ref_mul(x_coefficients[i], y_coefficients[j]));
        }
    }

    const fp2_value_type xi = fp6_value_type::field_type::extension_policy::non_residue;
    for (std::size_t i = product.size(); i-- > 6;) {
        product[i - 6] = ref_add(product[i - 6], ref_mul(product[i], xi));
    }

    return from_w_coefficients(
        {product[0], product[1], product[2], product[3], product[4], product[5]});
}

template<typename Rng>
fp12_value_type random_fp12(Rng &rng) {
    return nil::crypto3::algebra::random_element<fp12_field_type>(rng);
}

}    // namespace

BOOST_AUTO_TEST_SUITE(bn254_fp12_tests)

BOOST_AUTO_TEST_CASE(multiplication_matches_independent_reference) {
    boost::random::mt19937 rng(0x25412);

    for (std::size_t i = 0; i < random_samples; ++i) {
        const fp12_value_type x = random_fp12(rng);
        const fp12_value_type y = random_fp12(rng);

        BOOST_CHECK_EQUAL(x * y, ref_mul(x, y));

        fp12_value_type inplace = x;
        inplace *= y;
        BOOST_CHECK_EQUAL(inplace, ref_mul(x, y));
    }
}

BOOST_AUTO_TEST_CASE(square_matches_independent_reference) {
    boost::random::mt19937 rng(0x2545);

    for (std::size_t i = 0; i < random_samples; ++i) {
        const fp12_value_type x = random_fp12(rng);

        BOOST_CHECK_EQUAL(x.squared(), ref_mul(x, x));
    }
}

BOOST_AUTO_TEST_CASE(multiplication_identities) {
    boost::random::mt19937 rng(0x2541d);

    for (std::size_t i = 0; i < random_samples; ++i) {
        const fp12_value_type x = random_fp12(rng);

        BOOST_CHECK_EQUAL(x * fp12_value_type::zero(), fp12_value_type::zero());
        BOOST_CHECK_EQUAL(fp12_value_type::zero() * x, fp12_value_type::zero());
        BOOST_CHECK_EQUAL(x * fp12_value_type::one(), x);
        BOOST_CHECK_EQUAL(fp12_value_type::one() * x, x);
    }
}

BOOST_AUTO_TEST_CASE(inverse_is_multiplicative_identity) {
    boost::random::mt19937 rng(0x25417);

    for (std::size_t i = 0; i < random_samples; ++i) {
        const fp12_value_type x = random_fp12(rng);
        BOOST_TEST_REQUIRE(!x.is_zero());

        const fp12_value_type x_inv = x.inversed();
        BOOST_CHECK_EQUAL(ref_mul(x, x_inv), fp12_value_type::one());
        BOOST_CHECK_EQUAL(x * x_inv, fp12_value_type::one());
    }
}

BOOST_AUTO_TEST_SUITE_END()
