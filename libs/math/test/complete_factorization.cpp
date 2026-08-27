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

#define BOOST_TEST_MODULE complete_factorization_test

#include <array>
#include <cstddef>
#include <cstdint>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/factorization/complete_factorization.hpp>
#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/random/algebraic_engine.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using value_type = fq_field_type::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;
    using backend_type = polynomial_arithmetic::schoolbook_backend<value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    fq12_value_type fq12_value(std::size_t first_coordinate) {
        fq12_value_type value = fq12_value_type::zero();
        for (std::size_t coordinate = 0; coordinate < fq12_field_type::arity; ++coordinate) {
            value.coordinate(coordinate) = value_type(first_coordinate + coordinate);
        }
        return value;
    }

    template<polynomial_arithmetic::PolynomialBackend Backend>
    typename Backend::polynomial_type multiply(Backend &backend, const typename Backend::polynomial_type &left,
                                               const typename Backend::polynomial_type &right) {
        typename Backend::polynomial_type result;
        backend.multiply(result, left, right);
        return result;
    }

    template<polynomial_arithmetic::PolynomialBackend Backend>
    typename Backend::polynomial_type
        reconstruct(Backend &backend,
                    const math::polynomial_factorization_result<typename Backend::polynomial_type> &factorization) {
        typename Backend::polynomial_type result = {factorization.leading_coefficient};
        for (const auto &factor : factorization.factors) {
            for (std::size_t copy = 0; copy < factor.multiplicity; ++copy) {
                result = multiply(backend, result, factor.polynomial);
            }
        }
        return result;
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(complete_factorization_test_suite)

BOOST_AUTO_TEST_CASE(complete_factorization_recovers_irreducible_factors_and_multiplicities) {
    backend_type backend;
    const polynomial_type first_linear_factor = {value_type(1), value_type::one()};
    const polynomial_type second_linear_factor = {value_type(2), value_type::one()};
    // Three is a multiplicative generator of BN254 Fq and is therefore not a square, so X^2 - 3 is irreducible.
    const polynomial_type quadratic_factor = {-value_type(3), value_type::zero(), value_type::one()};

    polynomial_type input = multiply(backend, first_linear_factor, second_linear_factor);
    input = multiply(backend, input, quadratic_factor);
    input = multiply(backend, input, quadratic_factor);
    const value_type leading_coefficient(11);
    for (value_type &coefficient : input) {
        coefficient = coefficient * leading_coefficient;
    }

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    nil::crypto3::random::algebraic_engine<fq_field_type> generator(17);
    const auto result = math::complete_factorization<backend_type>(input, arithmetic_context, generator);

    BOOST_CHECK(result.complete);
    BOOST_CHECK(result.leading_coefficient == leading_coefficient);
    BOOST_REQUIRE_EQUAL(result.factors.size(), 3);
    BOOST_CHECK(reconstruct(backend, result) == input);

    std::size_t matched_factors = 0;
    for (const auto &factor : result.factors) {
        if (factor.polynomial == first_linear_factor || factor.polynomial == second_linear_factor) {
            BOOST_CHECK_EQUAL(factor.multiplicity, 1);
            ++matched_factors;
        } else if (factor.polynomial == quadratic_factor) {
            BOOST_CHECK_EQUAL(factor.multiplicity, 2);
            ++matched_factors;
        }
    }
    BOOST_CHECK_EQUAL(matched_factors, 3);
}

BOOST_AUTO_TEST_CASE(staged_complete_factorization_stops_after_the_reported_factor) {
    backend_type backend;
    const polynomial_type first_factor = {value_type(1), value_type::one()};
    const polynomial_type second_factor = {value_type(2), value_type::one()};
    const polynomial_type third_factor = {value_type(3), value_type::one()};
    polynomial_type input = multiply(backend, first_factor, second_factor);
    input = multiply(backend, input, third_factor);

    std::size_t callback_count = 0;
    auto callback = [&callback_count](const math::polynomial_factor<polynomial_type> &) {
        ++callback_count;
        return math::factorization_control::stop_factorization;
    };

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    nil::crypto3::random::algebraic_engine<fq_field_type> generator(29);
    const auto result = math::complete_factorization<backend_type>(input, arithmetic_context, generator, callback);

    BOOST_CHECK(!result.complete);
    BOOST_CHECK_EQUAL(callback_count, 1);
    BOOST_REQUIRE_EQUAL(result.factors.size(), 1);
    BOOST_CHECK_EQUAL(result.factors.front().polynomial.degree(), 1);
    BOOST_CHECK_EQUAL(result.factors.front().multiplicity, 1);
}

BOOST_AUTO_TEST_CASE(zero_and_constant_inputs_preserve_their_scalar_values) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    nil::crypto3::random::algebraic_engine<fq_field_type> generator(31);

    const auto zero =
        math::complete_factorization<backend_type>(polynomial_type {value_type::zero()}, arithmetic_context, generator);
    BOOST_CHECK(zero.complete);
    BOOST_CHECK(zero.leading_coefficient == value_type::zero());
    BOOST_CHECK(zero.factors.empty());

    const auto constant =
        math::complete_factorization<backend_type>(polynomial_type {value_type(13)}, arithmetic_context, generator);
    BOOST_CHECK(constant.complete);
    BOOST_CHECK(constant.leading_coefficient == value_type(13));
    BOOST_CHECK(constant.factors.empty());
}

BOOST_AUTO_TEST_CASE(complete_factorization_is_deterministic_for_a_seeded_generator) {
    backend_type backend;
    const polynomial_type first_factor = {value_type(1), value_type::one()};
    const polynomial_type second_factor = {value_type(2), value_type::one()};
    const polynomial_type third_factor = {value_type(3), value_type::one()};
    polynomial_type input = multiply(backend, first_factor, second_factor);
    input = multiply(backend, input, third_factor);

    auto factor_with_seed = [&input](std::uint32_t seed) {
        polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
        nil::crypto3::random::algebraic_engine<fq_field_type> generator(seed);
        return math::complete_factorization<backend_type>(input, arithmetic_context, generator);
    };

    const auto first_result = factor_with_seed(37);
    const auto second_result = factor_with_seed(37);
    BOOST_CHECK(first_result == second_result);
}

BOOST_AUTO_TEST_CASE(fq12_factorization_reconstructs_input_and_matches_the_mixed_radix_backend) {
    using extension_backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend_type = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using extension_polynomial_type = typename extension_backend_type::polynomial_type;

    extension_backend_type backend;
    const extension_polynomial_type first_factor = {-fq12_value(1), fq12_value_type::one()};
    const extension_polynomial_type second_factor = {-fq12_value(20), fq12_value_type::one()};
    const extension_polynomial_type repeated_factor = {-fq12_value(40), fq12_value_type::one()};
    extension_polynomial_type input = multiply(backend, first_factor, second_factor);
    input = multiply(backend, input, repeated_factor);
    input = multiply(backend, input, repeated_factor);
    const fq12_value_type leading_coefficient = fq12_value(60);
    math::scalar_multiplication(input, input, leading_coefficient);

    // These coefficients make the sampled polynomial equal first_factor. Cantor-Zassenhaus therefore finds a proper
    // factor in its initial GCD, which exercises generator injection and splitting without making this Fq12 unit test
    // pay for a full modular exponentiation.
    const std::array random_coefficients = {-fq12_value(1), fq12_value_type::one()};
    std::size_t next_coefficient = 0;
    auto generator = [&] { return random_coefficients[next_coefficient++]; };
    polynomial_arithmetic::polynomial_context<extension_backend_type> arithmetic_context;
    const auto result = math::complete_factorization<extension_backend_type>(input, arithmetic_context, generator);

    BOOST_CHECK(result.complete);
    BOOST_CHECK(result.leading_coefficient == leading_coefficient);
    BOOST_REQUIRE_EQUAL(result.factors.size(), 3);
    BOOST_CHECK(reconstruct(backend, result) == input);
    BOOST_CHECK_EQUAL(next_coefficient, random_coefficients.size());

    polynomial_arithmetic::polynomial_context_options options;
    // Force fast division so the comparison exercises mixed-radix arithmetic throughout the factorization pipeline.
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<mixed_radix_backend_type> mixed_radix_context(
        mixed_radix_backend_type(18), options);
    std::size_t next_mixed_radix_coefficient = 0;
    auto mixed_radix_generator = [&] { return random_coefficients[next_mixed_radix_coefficient++]; };
    const auto mixed_radix_result =
        math::complete_factorization<mixed_radix_backend_type>(input, mixed_radix_context, mixed_radix_generator);

    BOOST_CHECK(mixed_radix_result == result);
    BOOST_CHECK_EQUAL(next_mixed_radix_coefficient, random_coefficients.size());
}

BOOST_AUTO_TEST_SUITE_END()
