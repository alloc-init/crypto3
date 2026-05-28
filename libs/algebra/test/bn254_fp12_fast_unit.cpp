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

#include <nil/crypto3/algebra/fields/detail/extension_params/alt_bn128/detail/fp12_limb_ops.hpp>

using namespace nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops;

BOOST_AUTO_TEST_SUITE(bn254_fp12_fast_unit_tests)

#if defined(__x86_64__) && defined(__BMI2__) && defined(__ADX__) && (defined(__GNUC__) || defined(__clang__))
BOOST_AUTO_TEST_CASE(mul_4x4_x86_correct) {
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
        multiply_4x4_x86_bmi2_adx(z1, x, y);
        BOOST_REQUIRE_EQUAL_COLLECTIONS(z0.begin(), z0.end(), z1.begin(), z1.end());
    }
}
#endif

BOOST_AUTO_TEST_SUITE_END()
