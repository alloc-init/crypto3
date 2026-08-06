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

#define BOOST_TEST_MODULE geometric_sequence_domain_test

#include <cstddef>
#include <stdexcept>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/mersenne31.hpp>
#include <nil/crypto3/algebra/fields/mersenne31.hpp>
#include <nil/crypto3/math/domains/basic_radix2_domain.hpp>
#include <nil/crypto3/math/domains/geometric_sequence_domain.hpp>

using namespace nil::crypto3;

namespace {

    using bn254_fq = algebra::fields::alt_bn128<254>;

    template<typename Coefficients, typename ValueType>
    ValueType evaluate(const Coefficients &coefficients, const ValueType &point) {
        ValueType result = ValueType::zero();
        for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
            result = result * point + *it;
        }
        return result;
    }

    template<typename DomainType>
    auto construct_vanishing_polynomial_directly(DomainType &domain) {
        using value_type = decltype(domain.get_domain_element(0));

        std::vector<value_type> result(domain.size() + 1, value_type::zero());
        result[0] = value_type::one();
        for (std::size_t i = 0; i < domain.size(); ++i) {
            const value_type point = domain.get_domain_element(i);
            for (std::size_t degree = i + 1; degree > 0; --degree) {
                result[degree] = result[degree - 1] - point * result[degree];
            }
            result[0] = -(point * result[0]);
        }
        return result;
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(geometric_sequence_domain_test_suite)

BOOST_AUTO_TEST_CASE(rejects_invalid_sizes_and_repeated_points) {
    BOOST_CHECK_THROW((math::geometric_sequence_domain<bn254_fq>(0)), std::invalid_argument);
    BOOST_CHECK_THROW((math::geometric_sequence_domain<bn254_fq>(1)), std::invalid_argument);

    // The order of 2 in F_(2^31-1) is 31, so the point at index 31 repeats the point at index zero.
    BOOST_CHECK_NO_THROW((math::geometric_sequence_domain<algebra::fields::mersenne31>(31)));
    BOOST_CHECK_THROW((math::geometric_sequence_domain<algebra::fields::mersenne31>(32)), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(large_bn254_domain_supports_lagrange_evaluation) {
    constexpr std::size_t domain_size = 14923;
    math::geometric_sequence_domain<bn254_fq> domain(domain_size);

    BOOST_CHECK_EQUAL(domain.size(), domain_size);
    BOOST_CHECK_EQUAL(domain.get_domain_element(0), bn254_fq::value_type::one());
    BOOST_CHECK_EQUAL(domain.get_domain_element(domain_size - 1),
                      domain.get_domain_element(domain_size - 2) * domain.get_domain_element(1));

    const bn254_fq::value_type t = bn254_fq::value_type(1234567u);
    bn254_fq::value_type vanishing_at_t;
    math::evaluation_domain<bn254_fq> &abstract_domain = domain;
    const std::vector<bn254_fq::value_type> weights =
        abstract_domain.evaluate_all_lagrange_polynomials(t, vanishing_at_t);
    bn254_fq::value_type weight_sum = bn254_fq::value_type::zero();
    bn254_fq::value_type weighted_points = bn254_fq::value_type::zero();
    for (std::size_t i = 0; i < domain_size; ++i) {
        weight_sum += weights[i];
        weighted_points += weights[i] * domain.get_domain_element(i);
    }
    BOOST_CHECK_EQUAL(weight_sum, bn254_fq::value_type::one());
    BOOST_CHECK_EQUAL(weighted_points, t);
    BOOST_CHECK_EQUAL(vanishing_at_t, domain.compute_vanishing_polynomial(t));

    const math::polynomial<bn254_fq::value_type> vanishing_polynomial = domain.get_vanishing_polynomial();
    BOOST_CHECK_EQUAL(vanishing_polynomial.size(), domain_size + 1);
    BOOST_CHECK_EQUAL(vanishing_polynomial[domain_size], bn254_fq::value_type::one());
    BOOST_CHECK_EQUAL(evaluate(vanishing_polynomial, t), vanishing_at_t);
}

BOOST_AUTO_TEST_CASE(vanishing_polynomial_matches_direct_construction) {
    using value_type = bn254_fq::value_type;

    for (std::size_t domain_size = 2; domain_size <= 64; ++domain_size) {
        math::geometric_sequence_domain<bn254_fq> domain(domain_size);
        const math::polynomial<value_type> actual = domain.get_vanishing_polynomial();
        const std::vector<value_type> expected = construct_vanishing_polynomial_directly(domain);

        BOOST_REQUIRE_EQUAL(actual.size(), expected.size());
        for (std::size_t i = 0; i < actual.size(); ++i) {
            BOOST_CHECK_EQUAL(actual[i], expected[i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(vanishing_polynomial_handles_an_exact_order_generator) {
    using field_type = algebra::fields::mersenne31;
    using value_type = field_type::value_type;

    constexpr std::size_t domain_size = 31;
    math::geometric_sequence_domain<field_type> domain(domain_size);
    const math::polynomial<value_type> vanishing_polynomial = domain.get_vanishing_polynomial();

    BOOST_REQUIRE_EQUAL(vanishing_polynomial.size(), domain_size + 1);
    BOOST_CHECK_EQUAL(vanishing_polynomial[0], -value_type::one());
    BOOST_CHECK_EQUAL(vanishing_polynomial[domain_size], value_type::one());
    for (std::size_t i = 1; i < domain_size; ++i) {
        BOOST_CHECK_EQUAL(vanishing_polynomial[i], value_type::zero());
    }
}

BOOST_AUTO_TEST_CASE(combined_lagrange_api_has_a_default_domain_implementation) {
    using value_type = bn254_fq::value_type;

    math::basic_radix2_domain<bn254_fq> domain(2);
    const value_type t = value_type(19u);
    value_type vanishing_at_t;
    const std::vector<value_type> combined_weights = domain.evaluate_all_lagrange_polynomials(t, vanishing_at_t);

    BOOST_CHECK(combined_weights == domain.evaluate_all_lagrange_polynomials(t));
    BOOST_CHECK_EQUAL(vanishing_at_t, domain.compute_vanishing_polynomial(t));
}

BOOST_AUTO_TEST_CASE(lagrange_weights_match_the_definition) {
    using value_type = bn254_fq::value_type;

    constexpr std::size_t domain_size = 9;
    math::geometric_sequence_domain<bn254_fq> domain(domain_size);
    const value_type t = value_type(19u);
    value_type returned_vanishing_at_t;
    const std::vector<value_type> weights = domain.evaluate_all_lagrange_polynomials(t, returned_vanishing_at_t);

    value_type vanishing_at_t = value_type::one();
    for (std::size_t j = 0; j < domain_size; ++j) {
        vanishing_at_t *= t - domain.get_domain_element(j);
    }
    BOOST_CHECK_EQUAL(returned_vanishing_at_t, vanishing_at_t);

    for (std::size_t i = 0; i < domain_size; ++i) {
        value_type derivative = value_type::one();
        for (std::size_t j = 0; j < domain_size; ++j) {
            if (i != j) {
                derivative *= domain.get_domain_element(i) - domain.get_domain_element(j);
            }
        }
        const value_type expected =
            vanishing_at_t * (t - domain.get_domain_element(i)).inversed() * derivative.inversed();
        BOOST_CHECK_EQUAL(weights[i], expected);
    }
}

BOOST_AUTO_TEST_CASE(lagrange_weights_interpolate_and_are_unit_vectors_on_the_domain) {
    using value_type = bn254_fq::value_type;

    constexpr std::size_t domain_size = 17;
    math::geometric_sequence_domain<bn254_fq> domain(domain_size);
    std::vector<value_type> coefficients(domain_size, value_type::zero());
    std::vector<value_type> evaluations(domain_size, value_type::zero());
    for (std::size_t i = 0; i < domain_size; ++i) {
        coefficients[i] = value_type(static_cast<unsigned>(3 * i + 1));
    }
    for (std::size_t i = 0; i < domain_size; ++i) {
        evaluations[i] = evaluate(coefficients, domain.get_domain_element(i));
    }

    const value_type t = value_type(23u);
    const std::vector<value_type> weights = domain.evaluate_all_lagrange_polynomials(t);
    value_type interpolated = value_type::zero();
    for (std::size_t i = 0; i < domain_size; ++i) {
        interpolated += evaluations[i] * weights[i];
    }
    BOOST_CHECK_EQUAL(interpolated, evaluate(coefficients, t));

    for (std::size_t point_index = 0; point_index < domain_size; ++point_index) {
        value_type vanishing_at_domain_point;
        const std::vector<value_type> unit_weights =
            domain.evaluate_all_lagrange_polynomials(domain.get_domain_element(point_index), vanishing_at_domain_point);
        BOOST_CHECK_EQUAL(vanishing_at_domain_point, value_type::zero());
        for (std::size_t i = 0; i < domain_size; ++i) {
            BOOST_CHECK_EQUAL(unit_weights[i], i == point_index ? value_type::one() : value_type::zero());
        }
    }
}

BOOST_AUTO_TEST_CASE(add_poly_z_adds_the_scaled_vanishing_polynomial) {
    using value_type = bn254_fq::value_type;

    constexpr std::size_t domain_size = 9;
    math::geometric_sequence_domain<bn254_fq> domain(domain_size);
    std::vector<value_type> polynomial(domain_size + 1, value_type::zero());
    for (std::size_t i = 0; i < polynomial.size(); ++i) {
        polynomial[i] = value_type(static_cast<unsigned>(i + 1));
    }

    const value_type coefficient = value_type(7u);
    const std::vector<value_type> original = polynomial;

    domain.add_poly_z(coefficient, polynomial);

    std::vector<value_type> added_polynomial(polynomial.size(), value_type::zero());
    for (std::size_t i = 0; i < polynomial.size(); ++i) {
        added_polynomial[i] = polynomial[i] - original[i];
    }
    BOOST_CHECK_EQUAL(added_polynomial[domain_size], coefficient);
    for (std::size_t i = 0; i < domain_size; ++i) {
        BOOST_CHECK_EQUAL(evaluate(added_polynomial, domain.get_domain_element(i)), value_type::zero());
    }
}

BOOST_AUTO_TEST_SUITE_END()
