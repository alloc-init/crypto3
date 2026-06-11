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

    limb_array modulus_limbs() {
        return load_limbs(field_type::modulus_params.get_mod_obj().get_mod());
    }

    bool is_pR_bounded(const limb_array &input) {
        if (input[8] != 0u) {
            return false;
        }

        const limb_array p = modulus_limbs();
        for (size_t i = base_value_limb_count; i > 0; --i) {
            const size_t index = i - 1u;
            const limb upper = input[base_value_limb_count + index];
            if (upper < p[index]) {
                return true;
            }
            if (upper > p[index]) {
                return false;
            }
        }
        return false;
    }

    void force_pR_bound(limb_array &input) {
        const limb_array p = modulus_limbs();
        input[8] = 0u;
        while (!is_pR_bounded(input)) {
            subtract_limbs_portable<base_value_limb_count>(input.data() + base_value_limb_count, p.data());
        }
    }

    limb_array random_pR_bounded_limbs(boost::random::mt19937 &rng,
                                       boost::random::uniform_int_distribution<limb> &distribution) {
        limb_array result = random_limbs(rng, distribution, 8);
        force_pR_bound(result);
        return result;
    }

    void decrement_upper(limb_array &result) {
        limb borrow = 1u;
        for (size_t i = base_value_limb_count; borrow != 0u && i < 2u * base_value_limb_count; i++) {
            const limb current = result[i];
            result[i] = current - borrow;
            borrow = current == 0u;
        }
    }

    limb_array pR_minus_one() {
        const limb max = ~limb(0u);
        const limb_array p = modulus_limbs();
        limb_array result = {};
        for (size_t i = 0; i < base_value_limb_count; i++) {
            result[i] = max;
            result[base_value_limb_count + i] = p[i];
        }
        decrement_upper(result);
        return result;
    }

    void require_reduce_x86_matches_portable(const limb_array &input) {
        BOOST_REQUIRE(is_pR_bounded(input));
        limb_array expected = input;
        limb_array actual = input;
        montgomery_reduce_portable<field_type>(expected);
        montgomery_reduce_x86<field_type>(actual);
        BOOST_REQUIRE_EQUAL_COLLECTIONS(expected.begin(), expected.end(), actual.begin(), actual.end());
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

BOOST_AUTO_TEST_CASE(reduce_x86_random_pR_bounded_inputs) {
    boost::random::mt19937 rng(0x2545);
    boost::random::uniform_int_distribution<limb> d;
    for (size_t i = 0; i < 100; i++) {
        require_reduce_x86_matches_portable(random_pR_bounded_limbs(rng, d));
    }
}

BOOST_AUTO_TEST_CASE(reduce_x86_edge_pR_bounded_inputs) {
    const limb max = ~limb(0u);
    const limb_array p = modulus_limbs();

    limb_array zero = {};

    limb_array low_max = {};
    for (size_t i = 0; i < 4; i++) {
        low_max[i] = max;
    }

    limb_array high_modulus_minus_one = {};
    for (size_t i = 0; i < 4; i++) {
        high_modulus_minus_one[base_value_limb_count + i] = p[i];
    }
    decrement_upper(high_modulus_minus_one);

    limb_array alternating = {};
    for (size_t i = 0; i < 8; i++) {
        alternating[i] = (i % 2 == 0) ? max : 0u;
    }
    force_pR_bound(alternating);

    require_reduce_x86_matches_portable(zero);
    require_reduce_x86_matches_portable(low_max);
    require_reduce_x86_matches_portable(high_modulus_minus_one);
    require_reduce_x86_matches_portable(pR_minus_one());
    require_reduce_x86_matches_portable(alternating);
}
#else
BOOST_AUTO_TEST_CASE(dummy_test) {
    BOOST_REQUIRE(true);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
