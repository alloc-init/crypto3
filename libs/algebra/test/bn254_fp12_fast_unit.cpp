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

#define BOOST_TEST_MODULE algebra_bn254_fp12_fast_unit_test

#include <array>
#include <cstddef>

#include <boost/test/unit_test.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_ops.hpp>

using namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops;

BOOST_AUTO_TEST_SUITE(bn254_fp12_fast_unit_tests)

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
namespace {
    using field_type = nil::crypto3::algebra::fields::alt_bn128<254>;

    limb_array random_limbs(boost::random::mt19937 &rng,
                            boost::random::uniform_int_distribution<limb> &distribution,
                            size_t count) {
        limb_array result = {};
        for (size_t i = 0; i < count; i++) {
            result[i] = distribution(rng);
        }
        return result;
    }

    void require_reduce_x86_matches_portable(const limb_array &input) {
        limb_array expected = input;
        limb_array actual = input;
        montgomery_reduce_portable<field_type>(expected);
        montgomery_reduce_x86<field_type>(actual);
        BOOST_REQUIRE_EQUAL_COLLECTIONS(expected.begin(), expected.end(), actual.begin(), actual.end());
    }

    void require_reduce_product_4x4_matches(const limb_array &x, const limb_array &y) {
        limb_array product = {};
        multiply_4x4_portable(product, x, y);
        require_reduce_x86_matches_portable(product);
    }

    void require_reduce_product_5x5_matches(const limb_array &x, const limb_array &y) {
        limb_array product = {};
        multiply_5x5_portable(product, x, y);
        require_reduce_x86_matches_portable(product);
    }
}    // namespace

BOOST_AUTO_TEST_CASE(mul_4x4_x86_random) {
    boost::random::mt19937 rng(0x2545);
    boost::random::uniform_int_distribution<limb> d;
    for (size_t i = 0; i < 100; i++) {
        limb_array x, y;
        limb_array z0 = {};
        limb_array z1 = {};
        for (size_t i = 0; i < 4; i++) {
            x[i] = d(rng);
            y[i] = d(rng);
        }
        multiply_4x4_portable(z0, x, y);
        multiply_4x4_x86(z1, x, y);
        BOOST_REQUIRE_EQUAL_COLLECTIONS(z0.begin(), z0.end(), z1.begin(), z1.end());
    }
}

BOOST_AUTO_TEST_CASE(mul_5x5_x86_random) {
    boost::random::mt19937 rng(0x2545);
    boost::random::uniform_int_distribution<limb> d;
    for (size_t i = 0; i < 100; i++) {
        limb_array x, y;
        limb_array z0 = {};
        limb_array z1 = {};
        for (size_t i = 0; i < 5; i++) {
            x[i] = d(rng);
            y[i] = d(rng);
        }
        x[4] &= 1u;
        y[4] &= 1u;
        multiply_5x5_portable(z0, x, y);
        multiply_5x5_x86(z1, x, y);
        BOOST_REQUIRE_EQUAL_COLLECTIONS(z0.begin(), z0.end(), z1.begin(), z1.end());
    }
}

BOOST_AUTO_TEST_CASE(reduce_x86_random_products) {
    boost::random::mt19937 rng(0x2545);
    boost::random::uniform_int_distribution<limb> d;
    for (size_t i = 0; i < 100; i++) {
        const limb_array x = random_limbs(rng, d, 4);
        const limb_array y = random_limbs(rng, d, 4);
        require_reduce_product_4x4_matches(x, y);
    }
    for (size_t i = 0; i < 100; i++) {
        limb_array x = random_limbs(rng, d, 5);
        limb_array y = random_limbs(rng, d, 5);
        x[4] &= 1u;
        y[4] &= 1u;
        require_reduce_product_5x5_matches(x, y);
    }
}

BOOST_AUTO_TEST_CASE(reduce_x86_edge_products) {
    const limb max = ~limb(0u);

    limb_array all_max_4 = {};
    for (size_t i = 0; i < 4; i++) {
        all_max_4[i] = max;
    }

    limb_array alternating_a = {};
    limb_array alternating_b = {};
    for (size_t i = 0; i < 4; i++) {
        alternating_a[i] = (i % 2 == 0) ? max : 0u;
        alternating_b[i] = (i % 2 == 0) ? 0u : max;
    }

    limb_array high_5_a = all_max_4;
    limb_array high_5_b = all_max_4;
    high_5_a[4] = 1u;
    high_5_b[4] = 1u;

    require_reduce_product_4x4_matches(all_max_4, all_max_4);
    require_reduce_product_4x4_matches(all_max_4, alternating_a);
    require_reduce_product_4x4_matches(alternating_a, alternating_b);
    require_reduce_product_5x5_matches(high_5_a, high_5_b);
    require_reduce_product_5x5_matches(high_5_a, all_max_4);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
