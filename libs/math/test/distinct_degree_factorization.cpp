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
#include <vector>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/distinct_degree_factorization.hpp>
#include <nil/crypto3/math/polynomial/kaltofen_shoup_distinct_degree_factorization.hpp>
#include <nil/crypto3/math/polynomial/mixed_radix_backend.hpp>
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

    template<polynomial_arithmetic::PolynomialBackend Backend>
    void check_reference_and_kaltofen_shoup_agree(
        const typename Backend::polynomial_type &input,
        polynomial_arithmetic::polynomial_context<Backend> &reference_context,
        polynomial_arithmetic::polynomial_context<Backend> &kaltofen_shoup_context) {
        const auto reference = math::distinct_degree_factorization_reference<Backend>(input, reference_context);
        const auto kaltofen_shoup =
            math::distinct_degree_factorization_kaltofen_shoup<Backend>(input, kaltofen_shoup_context);
        BOOST_CHECK(kaltofen_shoup == reference);
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(distinct_degree_factorization_test_suite)

BOOST_AUTO_TEST_CASE(reference_and_kaltofen_shoup_factorizations_separate_known_irreducible_degrees) {
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
    const auto kaltofen_shoup_result =
        math::distinct_degree_factorization_kaltofen_shoup<backend_type>(input, arithmetic_context);

    BOOST_CHECK(result.complete);
    BOOST_CHECK(result.leading_coefficient == leading_coefficient);
    BOOST_REQUIRE_EQUAL(result.factors.size(), 3);
    BOOST_CHECK(result.factors[0] == math::distinct_degree_factor<polynomial_type>({degree_one_factors, 1}));
    BOOST_CHECK(result.factors[1] == math::distinct_degree_factor<polynomial_type>({degree_two_factor, 2}));
    BOOST_CHECK(result.factors[2] == math::distinct_degree_factor<polynomial_type>({degree_three_factor, 3}));
    BOOST_CHECK(reconstruct(backend, result) == input);
    BOOST_CHECK(kaltofen_shoup_result == result);
}

BOOST_AUTO_TEST_CASE(staged_reference_and_kaltofen_shoup_factorizations_stop_after_the_reported_degree) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type degree_one_factor = {fq_value_type(1), fq_value_type::one()};
    const polynomial_type degree_two_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type degree_three_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                                 fq_value_type::one()};
    polynomial_type input = multiply(backend, degree_one_factor, degree_two_factor);
    input = multiply(backend, input, degree_three_factor);

    auto make_callback = [](std::size_t &callback_count) {
        return [&callback_count](const math::distinct_degree_factor<polynomial_type> &factor) {
            ++callback_count;
            return factor.irreducible_factor_degree == 2 ? math::factorization_control::stop_factorization :
                                                           math::factorization_control::continue_factorization;
        };
    };

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::size_t reference_callback_count = 0;
    auto reference_callback = make_callback(reference_callback_count);
    const auto reference_result =
        math::distinct_degree_factorization_reference<backend_type>(input, arithmetic_context, reference_callback);

    std::size_t kaltofen_shoup_callback_count = 0;
    auto kaltofen_shoup_callback = make_callback(kaltofen_shoup_callback_count);
    const auto kaltofen_shoup_result = math::distinct_degree_factorization_kaltofen_shoup<backend_type>(
        input, arithmetic_context, kaltofen_shoup_callback);

    BOOST_CHECK(!reference_result.complete);
    BOOST_CHECK_EQUAL(reference_callback_count, 2);
    BOOST_CHECK_EQUAL(kaltofen_shoup_callback_count, 2);
    BOOST_CHECK(kaltofen_shoup_result == reference_result);
    BOOST_REQUIRE_EQUAL(reference_result.factors.size(), 2);
    BOOST_CHECK(reference_result.factors[0] == math::distinct_degree_factor<polynomial_type>({degree_one_factor, 1}));
    BOOST_CHECK(reference_result.factors[1] == math::distinct_degree_factor<polynomial_type>({degree_two_factor, 2}));
}

BOOST_AUTO_TEST_CASE(reference_and_kaltofen_shoup_factorizations_reject_a_repeated_factor) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type factor = {fq_value_type(1), fq_value_type::one()};
    const polynomial_type input = multiply(backend, factor, factor);
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;

    BOOST_CHECK_THROW(math::distinct_degree_factorization_reference<backend_type>(input, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::distinct_degree_factorization_kaltofen_shoup<backend_type>(input, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(zero_and_constant_polynomials_have_no_distinct_degree_factors) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const auto zero = math::distinct_degree_factorization_reference<backend_type>(
        polynomial_type {fq_value_type::zero()}, arithmetic_context);
    const auto kaltofen_shoup_zero = math::distinct_degree_factorization_kaltofen_shoup<backend_type>(
        polynomial_type {fq_value_type::zero()}, arithmetic_context);
    BOOST_CHECK(zero.complete);
    BOOST_CHECK(zero.leading_coefficient == fq_value_type::zero());
    BOOST_CHECK(zero.factors.empty());
    BOOST_CHECK(kaltofen_shoup_zero == zero);

    const auto constant = math::distinct_degree_factorization_reference<backend_type>(
        polynomial_type {fq_value_type(13)}, arithmetic_context);
    const auto kaltofen_shoup_constant = math::distinct_degree_factorization_kaltofen_shoup<backend_type>(
        polynomial_type {fq_value_type(13)}, arithmetic_context);
    BOOST_CHECK(constant.complete);
    BOOST_CHECK(constant.leading_coefficient == fq_value_type(13));
    BOOST_CHECK(constant.factors.empty());
    BOOST_CHECK(kaltofen_shoup_constant == constant);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_matches_the_reference_across_square_free_degrees) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> reference_context;
    polynomial_arithmetic::polynomial_context<backend_type> kaltofen_shoup_context;
    // Each polynomial is c*X^n - 3 with nonzero c and n below the field characteristic. Any common root of the
    // polynomial and its derivative would have to be zero, but zero is not a root, so every input is square-free.
    for (const std::size_t degree : {1U, 2U, 3U, 4U, 5U, 8U, 9U, 12U, 20U}) {
        BOOST_TEST_CONTEXT("degree = " << degree) {
            polynomial_type input(degree + 1, fq_value_type::zero());
            input.front() = -fq_value_type(3);
            input.back() = fq_value_type(degree + 1);
            check_reference_and_kaltofen_shoup_agree(input, reference_context, kaltofen_shoup_context);
        }
    }
}

BOOST_AUTO_TEST_CASE(reference_and_kaltofen_shoup_factorizations_support_native_extension_field_coefficients) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    schoolbook_backend backend;
    const polynomial_type first_linear_factor = {fq12_value(1), fq12_value_type::one()};
    const polynomial_type second_linear_factor = {fq12_value(20), fq12_value_type::one()};
    const polynomial_type input = multiply(backend, first_linear_factor, second_linear_factor);

    polynomial_arithmetic::polynomial_context<schoolbook_backend> reference_context;
    polynomial_arithmetic::polynomial_context<schoolbook_backend> kaltofen_shoup_context;
    const auto reference = math::distinct_degree_factorization_reference<schoolbook_backend>(input, reference_context);
    const auto kaltofen_shoup =
        math::distinct_degree_factorization_kaltofen_shoup<schoolbook_backend>(input, kaltofen_shoup_context);

    // Products of representatives modulo this degree-two polynomial contain at most three coefficients.
    polynomial_arithmetic::polynomial_context<mixed_radix_backend> mixed_radix_context {mixed_radix_backend(3)};
    const auto mixed_radix_kaltofen_shoup =
        math::distinct_degree_factorization_kaltofen_shoup<mixed_radix_backend>(input, mixed_radix_context);

    BOOST_CHECK(reference.complete);
    BOOST_REQUIRE_EQUAL(reference.factors.size(), 1);
    BOOST_CHECK_EQUAL(reference.factors.front().irreducible_factor_degree, 1);
    BOOST_CHECK(reference.factors.front().polynomial == input);
    BOOST_CHECK(kaltofen_shoup == reference);
    BOOST_CHECK(mixed_radix_kaltofen_shoup == reference);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_block_size_balances_baby_and_giant_steps) {
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(0), 1);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(1), 1);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(2), 1);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(8), 2);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(9), 3);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(18), 3);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(19), 4);
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

BOOST_AUTO_TEST_CASE(kaltofen_shoup_coarse_blocks_extract_degree_intervals) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    // Three generates the multiplicative group of BN254 Fq, so 3 is not a square and neither 3 nor 9 = 3^2 is a
    // cube. Consequently X^2 - 3, X^3 - 3, and X^3 - 9 are irreducible over this field.
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type first_cubic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                                fq_value_type::one()};
    const polynomial_type second_cubic_factor = {-fq_value_type(9), fq_value_type::zero(), fq_value_type::zero(),
                                                 fq_value_type::one()};
    const polynomial_type first_block = multiply(backend, linear_factor, quadratic_factor);
    const polynomial_type second_block = multiply(backend, first_cubic_factor, second_cubic_factor);
    const polynomial_type input = multiply(backend, first_block, second_block);

    constexpr std::size_t block_size = 2;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(input, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(block_size, frobenius_context,
                                                                                       arithmetic_context);

    polynomial_type giant_step = precomputation.baby_step(block_size);
    polynomial_type factor;
    math::detail::kaltofen_shoup_coarse_block_factor(factor, input, giant_step, block_size, precomputation,
                                                     frobenius_context, arithmetic_context);
    BOOST_CHECK(factor == first_block);

    precomputation.apply_block_frobenius(giant_step, giant_step, frobenius_context, arithmetic_context);
    math::detail::kaltofen_shoup_coarse_block_factor(factor, second_block, giant_step, 1, precomputation,
                                                     frobenius_context, arithmetic_context);
    BOOST_CHECK(factor == second_block);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_coarse_block_rejects_an_invalid_degree_count) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(2, frobenius_context,
                                                                                       arithmetic_context);
    const polynomial_type giant_step = precomputation.baby_step(2);
    polynomial_type output;

    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_coarse_block_factor(output, divisor, giant_step, 0, precomputation,
                                                                       frobenius_context, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_coarse_block_factor(output, divisor, giant_step, 3, precomputation,
                                                                       frobenius_context, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_fine_splitting_separates_exact_degrees) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    // Three generates the multiplicative group of BN254 Fq, so 3 is not a square and neither 3 nor 9 = 3^2 is a
    // cube. Consequently X^2 - 3, X^3 - 3, and X^3 - 9 are irreducible over this field.
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type first_cubic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                                fq_value_type::one()};
    const polynomial_type second_cubic_factor = {-fq_value_type(9), fq_value_type::zero(), fq_value_type::zero(),
                                                 fq_value_type::one()};
    const polynomial_type first_block = multiply(backend, linear_factor, quadratic_factor);
    const polynomial_type second_block = multiply(backend, first_cubic_factor, second_cubic_factor);
    const polynomial_type input = multiply(backend, first_block, second_block);

    constexpr std::size_t block_size = 2;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(input, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(block_size, frobenius_context,
                                                                                       arithmetic_context);

    std::vector<factor_type> factors;
    std::size_t callback_count = 0;
    auto record_factor = [&](const factor_type &) {
        ++callback_count;
        return math::factorization_control::continue_factorization;
    };

    polynomial_type giant_step = precomputation.baby_step(block_size);
    BOOST_CHECK(math::detail::kaltofen_shoup_split_coarse_block(factors, first_block, giant_step, 1, block_size,
                                                                precomputation, arithmetic_context, record_factor) ==
                math::factorization_control::continue_factorization);

    precomputation.apply_block_frobenius(giant_step, giant_step, frobenius_context, arithmetic_context);
    BOOST_CHECK(math::detail::kaltofen_shoup_split_coarse_block(factors, second_block, giant_step, 3, 1, precomputation,
                                                                arithmetic_context, record_factor) ==
                math::factorization_control::continue_factorization);

    BOOST_CHECK_EQUAL(callback_count, 3);
    BOOST_REQUIRE_EQUAL(factors.size(), 3);
    BOOST_CHECK(factors[0] == factor_type({linear_factor, 1}));
    BOOST_CHECK(factors[1] == factor_type({quadratic_factor, 2}));
    BOOST_CHECK(factors[2] == factor_type({second_block, 3}));
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_fine_splitting_honors_an_early_stop) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type input = multiply(backend, linear_factor, quadratic_factor);

    constexpr std::size_t block_size = 2;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(input, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(block_size, frobenius_context,
                                                                                       arithmetic_context);
    const polynomial_type giant_step = precomputation.baby_step(block_size);
    std::vector<factor_type> factors;

    const auto control = math::detail::kaltofen_shoup_split_coarse_block(
        factors, input, giant_step, 1, block_size, precomputation, arithmetic_context,
        [](const factor_type &) { return math::factorization_control::stop_factorization; });

    BOOST_CHECK(control == math::factorization_control::stop_factorization);
    BOOST_REQUIRE_EQUAL(factors.size(), 1);
    BOOST_CHECK(factors.front() == factor_type({linear_factor, 1}));
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_fine_splitting_rejects_invalid_degree_intervals) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type input = multiply(backend, linear_factor, quadratic_factor);

    constexpr std::size_t block_size = 2;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(input, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(block_size, frobenius_context,
                                                                                       arithmetic_context);
    const polynomial_type giant_step = precomputation.baby_step(block_size);
    std::vector<factor_type> factors;
    auto continue_factorization = [](const factor_type &) {
        return math::factorization_control::continue_factorization;
    };

    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_split_coarse_block(factors, input, giant_step, 0, 1, precomputation,
                                                                      arithmetic_context, continue_factorization),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_split_coarse_block(factors, input, giant_step, 1, 0, precomputation,
                                                                      arithmetic_context, continue_factorization),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_split_coarse_block(factors, input, giant_step, 1, 3, precomputation,
                                                                      arithmetic_context, continue_factorization),
                      std::invalid_argument);

    // The stated interval contains only degree one, so the unclassified quadratic must be rejected.
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_split_coarse_block(factors, input, giant_step, 1, 1, precomputation,
                                                                      arithmetic_context, continue_factorization),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_internal_driver_processes_multiple_blocks) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type first_cubic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                                fq_value_type::one()};
    const polynomial_type second_cubic_factor = {-fq_value_type(9), fq_value_type::zero(), fq_value_type::zero(),
                                                 fq_value_type::one()};
    const polynomial_type degree_three_factors = multiply(backend, first_cubic_factor, second_cubic_factor);
    polynomial_type input = multiply(backend, linear_factor, quadratic_factor);
    input = multiply(backend, input, degree_three_factors);

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::vector<factor_type> factors;
    std::size_t callback_count = 0;
    // A block size of two forces a full degree-1-to-2 block followed by a partial block containing degree three.
    const auto control = math::detail::kaltofen_shoup_factor_monic_square_free(
        factors, input, 2, arithmetic_context, [&](const factor_type &) {
            ++callback_count;
            return math::factorization_control::continue_factorization;
        });

    BOOST_CHECK(control == math::factorization_control::continue_factorization);
    BOOST_CHECK_EQUAL(callback_count, 3);
    BOOST_REQUIRE_EQUAL(factors.size(), 3);
    BOOST_CHECK(factors[0] == factor_type({linear_factor, 1}));
    BOOST_CHECK(factors[1] == factor_type({quadratic_factor, 2}));
    BOOST_CHECK(factors[2] == factor_type({degree_three_factors, 3}));
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_internal_driver_emits_the_final_irreducible_factor) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    const polynomial_type cubic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                          fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::vector<factor_type> factors;
    const auto control = math::detail::kaltofen_shoup_factor_monic_square_free(
        factors, cubic_factor, 2, arithmetic_context,
        [](const factor_type &) { return math::factorization_control::continue_factorization; });

    BOOST_CHECK(control == math::factorization_control::continue_factorization);
    BOOST_REQUIRE_EQUAL(factors.size(), 1);
    BOOST_CHECK(factors.front() == factor_type({cubic_factor, 3}));
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_internal_driver_rejects_a_zero_block_size) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    const polynomial_type input = {fq_value_type::one(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::vector<factor_type> factors;

    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_factor_monic_square_free(
                          factors, input, 0, arithmetic_context,
                          [](const factor_type &) { return math::factorization_control::continue_factorization; }),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_internal_driver_propagates_an_early_stop) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type cubic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                          fq_value_type::one()};
    polynomial_type input = multiply(backend, linear_factor, quadratic_factor);
    input = multiply(backend, input, cubic_factor);

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::vector<factor_type> factors;
    const auto control = math::detail::kaltofen_shoup_factor_monic_square_free(
        factors, input, 2, arithmetic_context, [](const factor_type &factor) {
            return factor.irreducible_factor_degree == 2 ? math::factorization_control::stop_factorization :
                                                           math::factorization_control::continue_factorization;
        });

    BOOST_CHECK(control == math::factorization_control::stop_factorization);
    BOOST_REQUIRE_EQUAL(factors.size(), 2);
    BOOST_CHECK(factors[0] == factor_type({linear_factor, 1}));
    BOOST_CHECK(factors[1] == factor_type({quadratic_factor, 2}));
}

BOOST_AUTO_TEST_SUITE_END()
