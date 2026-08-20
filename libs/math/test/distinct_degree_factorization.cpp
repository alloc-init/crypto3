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

#define BOOST_TEST_MODULE distinct_degree_factorization_test

#include <cstddef>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/distinct_degree_factorization.hpp>
#include <nil/crypto3/math/polynomial/kaltofen_shoup_distinct_degree_factorization.hpp>
#include <nil/crypto3/math/polynomial/schoolbook_backend.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq_value_type = fq_field_type::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;

    fq12_value_type fq12_value(std::size_t first_coordinate) {
        fq12_value_type value = fq12_value_type::zero();
        for (std::size_t i = 0; i < fq12_field_type::arity; ++i) {
            value.coordinate(i) = fq_value_type(first_coordinate + i);
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
                    const math::distinct_degree_factorization_result<typename Backend::polynomial_type> &result) {
        typename Backend::polynomial_type reconstructed = {result.leading_coefficient};
        for (const auto &factor : result.factors) {
            reconstructed = multiply(backend, reconstructed, factor.polynomial);
        }
        return reconstructed;
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(distinct_degree_factorization_test_suite)

BOOST_AUTO_TEST_CASE(reference_factorization_separates_known_irreducible_degrees) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type linear_one = {fq_value_type(1), fq_value_type::one()};
    const polynomial_type linear_two = {fq_value_type(2), fq_value_type::one()};
    const polynomial_type degree_one_factors = multiply(backend, linear_one, linear_two);

    // Three is a multiplicative generator of BN254 Fq, hence it is neither a square nor a cube. Therefore X^2 - 3
    // and X^3 - 3 are irreducible over this field.
    const polynomial_type degree_two_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type degree_three_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                                 fq_value_type::one()};
    const fq_value_type leading_coefficient(11);

    polynomial_type input = multiply(backend, degree_one_factors, degree_two_factor);
    input = multiply(backend, input, degree_three_factor);
    for (fq_value_type &coefficient : input) {
        coefficient *= leading_coefficient;
    }

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const auto result = math::distinct_degree_factorization_reference<backend_type>(input, arithmetic_context);

    BOOST_CHECK(result.complete);
    BOOST_CHECK(result.leading_coefficient == leading_coefficient);
    BOOST_REQUIRE_EQUAL(result.factors.size(), 3);
    BOOST_CHECK(result.factors[0] == math::distinct_degree_factor<polynomial_type>({degree_one_factors, 1}));
    BOOST_CHECK(result.factors[1] == math::distinct_degree_factor<polynomial_type>({degree_two_factor, 2}));
    BOOST_CHECK(result.factors[2] == math::distinct_degree_factor<polynomial_type>({degree_three_factor, 3}));
    BOOST_CHECK(reconstruct(backend, result) == input);
}

BOOST_AUTO_TEST_CASE(staged_factorization_stops_after_the_reported_degree) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type degree_one_factor = {fq_value_type(1), fq_value_type::one()};
    const polynomial_type degree_two_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type degree_three_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                                 fq_value_type::one()};
    polynomial_type input = multiply(backend, degree_one_factor, degree_two_factor);
    input = multiply(backend, input, degree_three_factor);

    std::size_t callback_count = 0;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const auto result = math::distinct_degree_factorization_reference<backend_type>(
        input, arithmetic_context, [&](const math::distinct_degree_factor<polynomial_type> &factor) {
            ++callback_count;
            return factor.irreducible_factor_degree == 2 ? math::factorization_control::stop_factorization :
                                                           math::factorization_control::continue_factorization;
        });

    BOOST_CHECK(!result.complete);
    BOOST_CHECK_EQUAL(callback_count, 2);
    BOOST_REQUIRE_EQUAL(result.factors.size(), 2);
    BOOST_CHECK_EQUAL(result.factors[0].irreducible_factor_degree, 1);
    BOOST_CHECK_EQUAL(result.factors[1].irreducible_factor_degree, 2);
    BOOST_CHECK(result.factors[1].polynomial == degree_two_factor);
}

BOOST_AUTO_TEST_CASE(reference_factorization_rejects_a_repeated_factor) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type factor = {fq_value_type(1), fq_value_type::one()};
    const polynomial_type input = multiply(backend, factor, factor);
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;

    BOOST_CHECK_THROW(math::distinct_degree_factorization_reference<backend_type>(input, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(zero_and_constant_polynomials_have_no_distinct_degree_factors) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const auto zero = math::distinct_degree_factorization_reference<backend_type>(
        polynomial_type {fq_value_type::zero()}, arithmetic_context);
    BOOST_CHECK(zero.complete);
    BOOST_CHECK(zero.leading_coefficient == fq_value_type::zero());
    BOOST_CHECK(zero.factors.empty());

    const auto constant = math::distinct_degree_factorization_reference<backend_type>(
        polynomial_type {fq_value_type(13)}, arithmetic_context);
    BOOST_CHECK(constant.complete);
    BOOST_CHECK(constant.leading_coefficient == fq_value_type(13));
    BOOST_CHECK(constant.factors.empty());
}

BOOST_AUTO_TEST_CASE(reference_factorization_supports_native_extension_field_coefficients) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type first_linear_factor = {fq12_value(1), fq12_value_type::one()};
    const polynomial_type second_linear_factor = {fq12_value(20), fq12_value_type::one()};
    const polynomial_type input = multiply(backend, first_linear_factor, second_linear_factor);

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const auto result = math::distinct_degree_factorization_reference<backend_type>(input, arithmetic_context);

    BOOST_CHECK(result.complete);
    BOOST_REQUIRE_EQUAL(result.factors.size(), 1);
    BOOST_CHECK_EQUAL(result.factors.front().irreducible_factor_degree, 1);
    BOOST_CHECK(result.factors.front().polynomial == input);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_precomputation_builds_baby_and_giant_frobenius_steps) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    constexpr std::size_t block_size = 3;
    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type x = {fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(block_size, frobenius_context,
                                                                                       arithmetic_context);

    BOOST_CHECK_EQUAL(precomputation.block_size(), block_size);
    polynomial_type expected_baby_step(x);
    for (std::size_t index = 0; index <= block_size; ++index) {
        BOOST_CHECK(precomputation.baby_step(index) == expected_baby_step);
        if (index != block_size) {
            math::frobenius_map(expected_baby_step, expected_baby_step, frobenius_context, arithmetic_context);
        }
    }

    polynomial_type expected_giant_step(x);
    polynomial_type giant_step = precomputation.baby_step(block_size);
    for (std::size_t block = 1; block <= 3; ++block) {
        math::frobenius_map(expected_giant_step, expected_giant_step, block_size, frobenius_context,
                            arithmetic_context);
        BOOST_CHECK(giant_step == expected_giant_step);
        if (block != 3) {
            precomputation.apply_block_frobenius(giant_step, giant_step, frobenius_context, arithmetic_context);
        }
    }
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_precomputation_rejects_a_zero_block_size) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    BOOST_CHECK_THROW(
        math::detail::kaltofen_shoup_frobenius_precomputation<backend_type>(0, frobenius_context, arithmetic_context),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_precomputation_reduces_the_initial_baby_step) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type linear_divisor = {fq_value_type(5), fq_value_type::one()};
    math::polynomial_frobenius_context<backend_type> linear_frobenius_context(linear_divisor, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> linear_precomputation(
        1, linear_frobenius_context, arithmetic_context);
    BOOST_CHECK(linear_precomputation.baby_step(0) == polynomial_type({-fq_value_type(5)}));

    const polynomial_type constant_divisor = {fq_value_type(7)};
    math::polynomial_frobenius_context<backend_type> constant_frobenius_context(constant_divisor, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> constant_precomputation(
        1, constant_frobenius_context, arithmetic_context);
    BOOST_CHECK(constant_precomputation.baby_step(0) == polynomial_type({fq_value_type::zero()}));
}

BOOST_AUTO_TEST_SUITE_END()
