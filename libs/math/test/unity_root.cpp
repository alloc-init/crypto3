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

#define BOOST_TEST_MODULE unity_root_test

#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/bls12.hpp>
#include <nil/crypto3/math/algorithms/unity_root.hpp>

using namespace nil::crypto3;

namespace {

    using bn254_fq = algebra::fields::alt_bn128_base_field<254>;

    constexpr std::size_t odd_smooth_order = 3 * 3 * 29 * 67;
    constexpr std::size_t even_smooth_order = 2 * odd_smooth_order;

    template<typename FieldType>
    void check_exact_order(std::size_t order, const std::initializer_list<std::size_t> &prime_factors) {
        using value_type = typename FieldType::value_type;
        using integral_type = typename FieldType::integral_type;

        const value_type omega = math::unity_root<FieldType>(order);
        const integral_type exponent = (FieldType::modulus - 1) / order;
        const value_type generator =
            value_type(algebra::fields::arithmetic_params<FieldType>::multiplicative_generator);

        BOOST_CHECK_EQUAL(omega, generator.pow(exponent));
        BOOST_CHECK_EQUAL(omega.pow(order), value_type::one());
        for (const std::size_t prime_factor : prime_factors) {
            BOOST_CHECK_NE(omega.pow(order / prime_factor), value_type::one());
        }
    }

    template<typename FieldType>
    typename FieldType::value_type previous_radix2_unity_root(std::size_t order) {
        using value_type = typename FieldType::value_type;

        const std::size_t log_order = std::ceil(std::log2(order));
        value_type omega = value_type(algebra::fields::arithmetic_params<FieldType>::root_of_unity);
        for (std::size_t i = algebra::fields::arithmetic_params<FieldType>::s; i > log_order; --i) {
            omega *= omega;
        }
        return omega;
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(unity_root_test_suite)

BOOST_AUTO_TEST_CASE(bn254_fq_supports_prime_power_and_smooth_composite_orders) {
    check_exact_order<bn254_fq>(3, {3});
    check_exact_order<bn254_fq>(9, {3});
    check_exact_order<bn254_fq>(odd_smooth_order, {3, 29, 67});
    check_exact_order<bn254_fq>(even_smooth_order, {2, 3, 29, 67});
}

BOOST_AUTO_TEST_CASE(bn254_fq_smooth_composite_roots_are_nested) {
    const auto odd_order_root = math::unity_root<bn254_fq>(odd_smooth_order);
    const auto even_order_root = math::unity_root<bn254_fq>(even_smooth_order);

    BOOST_CHECK_EQUAL(even_order_root.squared(), odd_order_root);
}

BOOST_AUTO_TEST_CASE(unsupported_orders_are_rejected) {
    BOOST_CHECK_THROW(math::unity_root<bn254_fq>(0), std::invalid_argument);
    BOOST_CHECK_THROW(math::unity_root<bn254_fq>(4), std::invalid_argument);
    BOOST_CHECK_THROW(math::unity_root<bn254_fq>(5), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(order_one_returns_one) {
    BOOST_CHECK_EQUAL(math::unity_root<bn254_fq>(1), bn254_fq::value_type::one());
}

BOOST_AUTO_TEST_CASE(radix2_roots_are_unchanged) {
    using radix2_field = algebra::fields::bls12_fr<381>;
    constexpr std::array<std::size_t, 6> orders = {1, 2, 4, 16, 256, 65536};

    for (const std::size_t order : orders) {
        BOOST_CHECK_EQUAL(math::unity_root<radix2_field>(order), previous_radix2_unity_root<radix2_field>(order));
    }
}

BOOST_AUTO_TEST_SUITE_END()
