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

#define BOOST_TEST_MODULE polynomial_gcd_test

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <tuple>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/arithmetic/gcd.hpp>
#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = nil::crypto3::math::polynomial_arithmetic;
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

    std::uint64_t next_value(std::uint64_t &state) {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return state % 1000003 + 1;
    }

    math::polynomial<fq_value_type> fq_polynomial(std::size_t coefficient_count, std::uint64_t seed) {
        math::polynomial<fq_value_type> result(coefficient_count, fq_value_type::zero());
        for (fq_value_type &coefficient : result) {
            coefficient = fq_value_type(next_value(seed));
        }
        return result;
    }

    math::polynomial<fq12_value_type> fq12_polynomial(std::size_t coefficient_count, std::uint64_t seed) {
        math::polynomial<fq12_value_type> result(coefficient_count, fq12_value_type::zero());
        for (fq12_value_type &coefficient : result) {
            for (std::size_t i = 0; i < fq12_field_type::arity; ++i) {
                coefficient.coordinate(i) = fq_value_type(next_value(seed));
            }
        }
        return result;
    }

    template<polynomial_arithmetic::PolynomialBackend Backend>
    typename Backend::polynomial_type product(Backend backend, const typename Backend::polynomial_type &left,
                                              const typename Backend::polynomial_type &right) {
        typename Backend::polynomial_type output;
        backend.multiply(output, left, right);
        return output;
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(polynomial_gcd_test_suite)

BOOST_AUTO_TEST_CASE(gcd_recovers_the_monic_common_factor_and_is_alias_safe) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type common_factor = {fq_value_type(2), fq_value_type(4), fq_value_type(2)};
    const polynomial_type expected_gcd = {fq_value_type::one(), fq_value_type(2), fq_value_type::one()};
    const polynomial_type left =
        product(backend_type {}, common_factor, polynomial_type({fq_value_type(3), fq_value_type::one()}));
    const polynomial_type right =
        product(backend_type {}, common_factor,
                polynomial_type({fq_value_type(5), fq_value_type::zero(), fq_value_type::one()}));

    polynomial_arithmetic::polynomial_context<backend_type> context;
    polynomial_type output;
    math::gcd(output, left, right, context);
    BOOST_CHECK(output == expected_gcd);

    polynomial_type left_alias = left;
    math::gcd(left_alias, left_alias, right, context);
    BOOST_CHECK(left_alias == expected_gcd);

    polynomial_type right_alias = right;
    math::gcd(right_alias, left, right_alias, context);
    BOOST_CHECK(right_alias == expected_gcd);
}

BOOST_AUTO_TEST_CASE(gcd_handles_zero_and_relatively_prime_inputs) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type zero = {fq_value_type::zero()};
    const polynomial_type input = {fq_value_type(2), fq_value_type(4), fq_value_type(2)};
    const polynomial_type monic_input = {fq_value_type::one(), fq_value_type(2), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> context;
    polynomial_type output;

    math::gcd(output, zero, zero, context);
    BOOST_CHECK(output == zero);
    math::gcd(output, input, zero, context);
    BOOST_CHECK(output == monic_input);
    math::gcd(output, zero, input, context);
    BOOST_CHECK(output == monic_input);

    const polynomial_type left = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type right = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    math::gcd(output, left, right, context);
    BOOST_CHECK(output == polynomial_type({fq_value_type::one()}));
}

BOOST_AUTO_TEST_CASE(gcd_divrem_step_rejects_invalid_operands) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> context;
    polynomial_type quotient;
    polynomial_type remainder;
    const polynomial_type smaller = {fq_value_type::one()};
    const polynomial_type larger = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type zero = {fq_value_type::zero()};

    BOOST_CHECK_THROW(math::detail::gcd_divrem_step(quotient, remainder, smaller, larger, context),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::gcd_divrem_step(quotient, remainder, larger, zero, context), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(gcd_supports_extension_coefficients_and_mixed_radix_multiplication) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    const polynomial_type common_factor = {fq12_value(1), fq12_value(13), fq12_value(25)};
    const polynomial_type left =
        product(schoolbook_backend {}, common_factor, polynomial_type({fq12_value(37), fq12_value(49)}));
    const polynomial_type right = product(schoolbook_backend {}, common_factor,
                                          polynomial_type({fq12_value(61), fq12_value(73), fq12_value(85)}));
    polynomial_type expected_gcd;
    math::make_monic(expected_gcd, common_factor);

    polynomial_arithmetic::polynomial_context<schoolbook_backend> schoolbook_context;
    polynomial_type schoolbook_result;
    math::gcd(schoolbook_result, left, right, schoolbook_context);
    BOOST_CHECK(schoolbook_result == expected_gcd);

    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<mixed_radix_backend> mixed_radix_context(mixed_radix_backend(18),
                                                                                       options);
    polynomial_type mixed_radix_result;
    math::gcd(mixed_radix_result, left, right, mixed_radix_context);
    BOOST_CHECK(mixed_radix_result == expected_gcd);
}

BOOST_AUTO_TEST_CASE(half_gcd_matrix_preserves_the_pair_with_and_without_an_iterative_basecase) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type first = fq_polynomial(121, 0x123456789abcdefULL);
    const polynomial_type second = fq_polynomial(120, 0xfedcba987654321ULL);
    polynomial_arithmetic::polynomial_context<backend_type> context;

    for (const std::size_t basecase_cutoff : {0, 8}) {
        math::detail::half_gcd_matrix<backend_type> transformation;
        polynomial_type reduced_first;
        polynomial_type reduced_second;
        math::detail::half_gcd_reduce(transformation, reduced_first, reduced_second, first, second, basecase_cutoff,
                                      context);

        BOOST_CHECK(math::is_zero(reduced_second) || reduced_second.size() <= first.size() / 2);
        polynomial_type reconstructed_first;
        polynomial_type reconstructed_second;
        math::detail::apply_half_gcd_matrix(reconstructed_first, reconstructed_second, transformation, first, second,
                                            context);
        BOOST_CHECK(reconstructed_first == reduced_first);
        BOOST_CHECK(reconstructed_second == reduced_second);
    }
}

BOOST_AUTO_TEST_CASE(half_gcd_rejects_invalid_entry_pairs) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> context;
    polynomial_type first_output;
    polynomial_type second_output;
    const polynomial_type zero = {fq_value_type::zero()};
    const polynomial_type constant = {fq_value_type::one()};
    const polynomial_type linear = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type noncanonical = {fq_value_type::one(), fq_value_type::zero()};

    BOOST_CHECK_THROW(math::detail::half_gcd_reduce(first_output, second_output, zero, zero, 8, context),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::half_gcd_reduce(first_output, second_output, constant, linear, 8, context),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::half_gcd_reduce(first_output, second_output, linear, linear, 8, context),
                      std::invalid_argument);

    math::detail::half_gcd_matrix<backend_type> transformation;
    BOOST_CHECK_THROW(
        math::detail::half_gcd_reduce(transformation, first_output, second_output, noncanonical, zero, 8, context),
        std::invalid_argument);

    math::detail::half_gcd_reduce(first_output, second_output, constant, zero, 8, context);
    BOOST_CHECK(first_output == constant);
    BOOST_CHECK(second_output == zero);
}

BOOST_AUTO_TEST_CASE(half_gcd_matches_euclidean_gcd_for_extension_coefficients) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    const polynomial_type common_factor = fq12_polynomial(5, 0x123456789abcdefULL);
    const polynomial_type left =
        product(schoolbook_backend {}, common_factor, fq12_polynomial(26, 0x314159265358979ULL));
    const polynomial_type right =
        product(schoolbook_backend {}, common_factor, fq12_polynomial(25, 0xfedcba987654321ULL));

    polynomial_arithmetic::polynomial_context_options euclidean_options;
    euclidean_options.gcd_half_gcd_cutoff = 0;
    polynomial_arithmetic::polynomial_context<schoolbook_backend> euclidean_context(schoolbook_backend {},
                                                                                    euclidean_options);
    polynomial_type expected;
    math::gcd(expected, left, right, euclidean_context);

    polynomial_arithmetic::polynomial_context_options half_gcd_options;
    half_gcd_options.half_gcd_basecase_cutoff = 4;
    half_gcd_options.gcd_half_gcd_cutoff = 1;
    polynomial_arithmetic::polynomial_context<mixed_radix_backend> half_gcd_context(mixed_radix_backend(174),
                                                                                    half_gcd_options);
    polynomial_type result;
    math::gcd(result, left, right, half_gcd_context);
    BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_CASE(half_gcd_matches_euclidean_gcd_across_odd_and_even_sizes) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context_options euclidean_options;
    euclidean_options.gcd_half_gcd_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> euclidean_context(backend_type {}, euclidean_options);

    polynomial_arithmetic::polynomial_context_options half_gcd_options;
    half_gcd_options.half_gcd_basecase_cutoff = 4;
    half_gcd_options.gcd_half_gcd_cutoff = 1;
    polynomial_arithmetic::polynomial_context<backend_type> half_gcd_context(backend_type {}, half_gcd_options);

    for (const auto [common_size, left_cofactor_size, right_cofactor_size] :
         {std::tuple<std::size_t, std::size_t, std::size_t> {2, 6, 5},
          {3, 15, 14},
          {5, 27, 25},
          {8, 56, 53},
          {13, 109, 108}}) {
        const polynomial_type common = fq_polynomial(common_size, common_size);
        const polynomial_type left =
            product(backend_type {}, common, fq_polynomial(left_cofactor_size, left_cofactor_size + 100));
        const polynomial_type right =
            product(backend_type {}, common, fq_polynomial(right_cofactor_size, right_cofactor_size + 200));

        polynomial_type expected;
        polynomial_type result;
        math::gcd(expected, left, right, euclidean_context);
        math::gcd(result, left, right, half_gcd_context);
        BOOST_CHECK(result == expected);
    }
}

BOOST_AUTO_TEST_CASE(gcd_is_alias_safe_when_half_gcd_is_forced) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type common = fq_polynomial(5, 0x123456789abcdefULL);
    const polynomial_type left = product(backend_type {}, common, fq_polynomial(56, 0x314159265358979ULL));
    const polynomial_type right = product(backend_type {}, common, fq_polynomial(53, 0xfedcba987654321ULL));

    polynomial_arithmetic::polynomial_context_options euclidean_options;
    euclidean_options.gcd_half_gcd_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> euclidean_context(backend_type {}, euclidean_options);
    polynomial_type expected;
    math::gcd(expected, left, right, euclidean_context);

    polynomial_arithmetic::polynomial_context_options half_gcd_options;
    half_gcd_options.half_gcd_basecase_cutoff = 4;
    half_gcd_options.gcd_half_gcd_cutoff = 1;
    polynomial_arithmetic::polynomial_context<backend_type> half_gcd_context(backend_type {}, half_gcd_options);

    polynomial_type left_alias = left;
    math::gcd(left_alias, left_alias, right, half_gcd_context);
    BOOST_CHECK(left_alias == expected);

    polynomial_type right_alias = right;
    math::gcd(right_alias, left, right_alias, half_gcd_context);
    BOOST_CHECK(right_alias == expected);
}

BOOST_AUTO_TEST_SUITE_END()
