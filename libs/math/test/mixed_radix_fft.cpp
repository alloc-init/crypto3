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

#define BOOST_TEST_MODULE mixed_radix_fft_test

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/algebra/random_element.hpp>
#include <nil/crypto3/math/algorithms/mixed_radix_fft.hpp>
#include <nil/crypto3/math/polynomial/basic_operations.hpp>

using namespace nil::crypto3;

namespace {

    using bn254_fq = algebra::fields::alt_bn128_base_field<254>;
    using bn254_fq_value = bn254_fq::value_type;
    using bn254_fq12 = algebra::fields::fp12_2over3over2<algebra::fields::alt_bn128<254>>;
    using bn254_fq12_value = bn254_fq12::value_type;

    constexpr std::size_t odd_smooth_order = 3 * 3 * 29 * 67;
    constexpr std::size_t even_smooth_order = 2 * odd_smooth_order;

    template<typename ValueType, typename RootType>
    std::vector<ValueType> naive_dft(const std::vector<ValueType> &coefficients, const RootType &omega) {
        std::vector<ValueType> evaluations(coefficients.size(), ValueType::zero());
        for (std::size_t evaluation_index = 0; evaluation_index < coefficients.size(); ++evaluation_index) {
            const RootType evaluation_point = omega.pow(evaluation_index);
            RootType point_power = RootType::one();
            for (const ValueType &coefficient : coefficients) {
                evaluations[evaluation_index] += coefficient * point_power;
                point_power = point_power * evaluation_point;
            }
        }
        return evaluations;
    }

    std::vector<bn254_fq_value> fq_values(std::size_t size) {
        std::vector<bn254_fq_value> values;
        values.reserve(size);
        for (std::size_t i = 0; i < size; ++i) {
            values.emplace_back(17 * i + i % 11 + 1);
        }
        return values;
    }

    std::vector<bn254_fq12_value> fq12_values(std::size_t size, std::uint32_t seed) {
        boost::random::mt19937 rng(seed);
        std::vector<bn254_fq12_value> values;
        values.reserve(size);
        for (std::size_t i = 0; i < size; ++i) {
            values.push_back(algebra::random_element<bn254_fq12>(rng));
        }
        return values;
    }

    template<typename ValueType>
    std::vector<ValueType> schoolbook_multiply(const std::vector<ValueType> &left,
                                               const std::vector<ValueType> &right) {
        std::vector<ValueType> result(left.size() + right.size() - 1, ValueType::zero());
        for (std::size_t i = 0; i < left.size(); ++i) {
            if (left[i].is_zero()) {
                continue;
            }
            for (std::size_t j = 0; j < right.size(); ++j) {
                if (!right[j].is_zero()) {
                    result[i + j] += left[i] * right[j];
                }
            }
        }
        return result;
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(mixed_radix_fft_test_suite)

BOOST_AUTO_TEST_CASE(plan_exposes_exact_size_and_prime_radices) {
    const math::mixed_radix_fft_plan<bn254_fq> plan(even_smooth_order);
    const std::vector<std::size_t> expected_radices = {2, 3, 3, 29, 67};

    BOOST_CHECK_EQUAL(plan.size(), even_smooth_order);
    BOOST_CHECK(plan.radices() == expected_radices);
    BOOST_CHECK_EQUAL(plan.omega().pow(plan.size()), bn254_fq_value::one());
}

BOOST_AUTO_TEST_CASE(plan_rejects_zero_and_unsupported_sizes) {
    BOOST_CHECK_THROW(math::mixed_radix_fft_plan<bn254_fq>(0), std::invalid_argument);
    BOOST_CHECK_THROW(math::mixed_radix_fft_plan<bn254_fq>(5), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(small_fq_transforms_match_naive_dft_and_round_trip) {
    constexpr std::array<std::size_t, 8> sizes = {1, 2, 3, 6, 9, 29, 67, 87};

    for (const std::size_t size : sizes) {
        BOOST_TEST_CONTEXT("size = " << size) {
            const math::mixed_radix_fft_plan<bn254_fq> plan(size);
            const std::vector<bn254_fq_value> coefficients = fq_values(size);
            std::vector<bn254_fq_value> actual = coefficients;

            plan.fft(actual);
            BOOST_CHECK(actual == naive_dft(coefficients, plan.omega()));

            plan.inverse_fft(actual);
            BOOST_CHECK(actual == coefficients);
        }
    }
}

BOOST_AUTO_TEST_CASE(fq_transforms_round_trip_at_large_smooth_orders) {
    constexpr std::array<std::size_t, 2> sizes = {odd_smooth_order, even_smooth_order};

    for (const std::size_t size : sizes) {
        BOOST_TEST_CONTEXT("size = " << size) {
            const math::mixed_radix_fft_plan<bn254_fq> plan(size);
            const std::vector<bn254_fq_value> coefficients = fq_values(size);
            std::vector<bn254_fq_value> actual = coefficients;

            plan.fft(actual);
            plan.inverse_fft(actual);
            BOOST_CHECK(actual == coefficients);
        }
    }
}

BOOST_AUTO_TEST_CASE(fq12_transforms_match_naive_dft_and_round_trip) {
    const math::mixed_radix_fft_plan<bn254_fq> plan(6);
    const std::vector<bn254_fq12_value> coefficients = fq12_values(6, 0xF012);
    std::vector<bn254_fq12_value> actual = coefficients;

    plan.fft(actual);
    BOOST_CHECK(actual == naive_dft(coefficients, plan.omega()));

    plan.inverse_fft(actual);
    BOOST_CHECK(actual == coefficients);
}

BOOST_AUTO_TEST_CASE(medium_fq12_transform_round_trips) {
    constexpr std::size_t size = 3 * 3 * 67;
    const math::mixed_radix_fft_plan<bn254_fq> plan(size);
    const std::vector<bn254_fq12_value> coefficients = fq12_values(size, 0xF01267);
    std::vector<bn254_fq12_value> actual = coefficients;

    plan.fft(actual);
    plan.inverse_fft(actual);
    BOOST_CHECK(actual == coefficients);
}

BOOST_AUTO_TEST_CASE(short_inputs_are_zero_padded) {
    const math::mixed_radix_fft_plan<bn254_fq> plan(6);
    const std::vector<bn254_fq_value> short_coefficients = {bn254_fq_value(3), bn254_fq_value(5), bn254_fq_value(7)};
    std::vector<bn254_fq_value> padded_coefficients = short_coefficients;
    padded_coefficients.resize(plan.size(), bn254_fq_value::zero());
    std::vector<bn254_fq_value> actual = short_coefficients;

    plan.fft(actual);
    BOOST_CHECK(actual == naive_dft(padded_coefficients, plan.omega()));

    plan.inverse_fft(actual);
    BOOST_CHECK(actual == padded_coefficients);

    const std::vector<bn254_fq_value> short_evaluations = {bn254_fq_value(11), bn254_fq_value(13), bn254_fq_value(17)};
    std::vector<bn254_fq_value> padded_evaluations = short_evaluations;
    padded_evaluations.resize(plan.size(), bn254_fq_value::zero());
    actual = short_evaluations;

    plan.inverse_fft(actual);
    plan.fft(actual);
    BOOST_CHECK(actual == padded_evaluations);
}

BOOST_AUTO_TEST_CASE(oversized_inputs_are_rejected) {
    const math::mixed_radix_fft_plan<bn254_fq> plan(6);
    std::vector<bn254_fq_value> values(7, bn254_fq_value::one());

    BOOST_CHECK_THROW(plan.fft(values), std::invalid_argument);
    BOOST_CHECK_THROW(plan.inverse_fft(values), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(zero_constant_and_sparse_polynomials_transform_correctly) {
    const math::mixed_radix_fft_plan<bn254_fq> plan(6);

    std::vector<bn254_fq_value> zero(plan.size(), bn254_fq_value::zero());
    plan.fft(zero);
    BOOST_CHECK(zero == std::vector<bn254_fq_value>(plan.size(), bn254_fq_value::zero()));
    plan.inverse_fft(zero);
    BOOST_CHECK(zero == std::vector<bn254_fq_value>(plan.size(), bn254_fq_value::zero()));

    const bn254_fq_value constant_value(19);
    std::vector<bn254_fq_value> constant = {constant_value};
    plan.fft(constant);
    BOOST_CHECK(constant == std::vector<bn254_fq_value>(plan.size(), constant_value));
    plan.inverse_fft(constant);
    std::vector<bn254_fq_value> expected_constant(plan.size(), bn254_fq_value::zero());
    expected_constant[0] = constant_value;
    BOOST_CHECK(constant == expected_constant);

    std::vector<bn254_fq_value> sparse(plan.size(), bn254_fq_value::zero());
    sparse[1] = bn254_fq_value(23);
    sparse[5] = bn254_fq_value(31);
    const std::vector<bn254_fq_value> sparse_coefficients = sparse;
    plan.fft(sparse);
    BOOST_CHECK(sparse == naive_dft(sparse_coefficients, plan.omega()));
    plan.inverse_fft(sparse);
    BOOST_CHECK(sparse == sparse_coefficients);
}

BOOST_AUTO_TEST_CASE(fq_polynomial_multiplication_matches_schoolbook) {
    const math::mixed_radix_fft_plan<bn254_fq> plan(18);
    const std::vector<bn254_fq_value> left = fq_values(8);
    const std::vector<bn254_fq_value> right = fq_values(7);
    std::vector<bn254_fq_value> result;

    math::multiplication(result, left, right, plan);

    BOOST_CHECK(result == schoolbook_multiply(left, right));
}

BOOST_AUTO_TEST_CASE(fq12_polynomial_multiplication_matches_schoolbook) {
    const math::mixed_radix_fft_plan<bn254_fq> plan(18);
    const std::vector<bn254_fq12_value> left = fq12_values(8, 0xA012);
    const std::vector<bn254_fq12_value> right = fq12_values(7, 0xB012);
    std::vector<bn254_fq12_value> result;

    math::multiplication(result, left, right, plan);

    BOOST_CHECK(result == schoolbook_multiply(left, right));
}

BOOST_AUTO_TEST_CASE(polynomial_multiplication_supports_output_aliasing) {
    const math::mixed_radix_fft_plan<bn254_fq> plan(9);
    const std::vector<bn254_fq_value> left = fq_values(4);
    const std::vector<bn254_fq_value> right = fq_values(4);
    const std::vector<bn254_fq_value> expected = schoolbook_multiply(left, right);

    std::vector<bn254_fq_value> output = left;
    math::multiplication(output, output, right, plan);
    BOOST_CHECK(output == expected);

    output = right;
    math::multiplication(output, left, output, plan);
    BOOST_CHECK(output == expected);
}

BOOST_AUTO_TEST_CASE(polynomial_multiplication_preserves_exact_output_length) {
    const math::mixed_radix_fft_plan<bn254_fq> plan(6);
    const std::vector<bn254_fq_value> left = {bn254_fq_value(1), bn254_fq_value(2), bn254_fq_value::zero()};
    const std::vector<bn254_fq_value> right = {bn254_fq_value(3), bn254_fq_value(4), bn254_fq_value::zero()};
    std::vector<bn254_fq_value> result;

    math::multiplication(result, left, right, plan);

    BOOST_CHECK_EQUAL(result.size(), left.size() + right.size() - 1);
    BOOST_CHECK(result == schoolbook_multiply(left, right));
}

BOOST_AUTO_TEST_CASE(polynomial_multiplication_rejects_invalid_inputs) {
    const math::mixed_radix_fft_plan<bn254_fq> plan(6);
    const std::vector<bn254_fq_value> coefficients = fq_values(4);
    const std::vector<bn254_fq_value> empty;
    std::vector<bn254_fq_value> result;

    BOOST_CHECK_THROW(math::multiplication(result, empty, coefficients, plan), std::invalid_argument);
    BOOST_CHECK_THROW(math::multiplication(result, coefficients, empty, plan), std::invalid_argument);
    BOOST_CHECK_THROW(math::multiplication(result, coefficients, coefficients, plan), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(fq12_polynomial_multiplication_uses_the_larger_smooth_order) {
    constexpr std::size_t operand_size = odd_smooth_order / 2 + 2;
    const std::vector<bn254_fq12_value> samples = fq12_values(6, 0x34974);
    std::vector<bn254_fq12_value> left(operand_size, bn254_fq12_value::zero());
    std::vector<bn254_fq12_value> right(operand_size, bn254_fq12_value::zero());
    left[0] = samples[0];
    left[113] = samples[1];
    left.back() = samples[2];
    right[0] = samples[3];
    right[257] = samples[4];
    right.back() = samples[5];

    const std::vector<bn254_fq12_value> expected = schoolbook_multiply(left, right);
    const math::mixed_radix_fft_plan<bn254_fq> plan(even_smooth_order);
    std::vector<bn254_fq12_value> result;

    BOOST_REQUIRE_GT(expected.size(), odd_smooth_order);
    math::multiplication(result, left, right, plan);
    BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_SUITE_END()
