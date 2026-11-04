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

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/algebra/fields/koalabear.hpp>

#include <nil/crypto3/math/polynomial/arithmetic/gcd.hpp>
#include <nil/crypto3/math/polynomial/backends/bn254_small_prime_fft_backend.hpp>
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

    struct matrix_product_counts {
        std::array<std::size_t, 2> accepted {};
        std::size_t declined = 0;
    };

    struct counting_matrix_backend : polynomial_arithmetic::bn254_small_prime_fft_backend {
        using base_type = polynomial_arithmetic::bn254_small_prime_fft_backend;

        counting_matrix_backend(std::size_t max_size, matrix_product_counts &counts, std::size_t byte_limit) :
            base_type(max_size, std::numeric_limits<std::size_t>::max(), byte_limit), counts(counts) {
        }

        template<std::size_t Columns>
            requires(Columns == 1 || Columns == 2)
        bool try_multiply_matrix(const std::array<polynomial_type *, 2 * Columns> &output,
                                 const std::array<const polynomial_type *, 4> &left,
                                 const std::array<const polynomial_type *, 2 * Columns> &right) {
            const bool accepted = base_type::try_multiply_matrix<Columns>(output, left, right);
            if (accepted) {
                ++counts.accepted[Columns - 1];
            } else {
                ++counts.declined;
            }
            return accepted;
        }

        matrix_product_counts &counts;
    };

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

    // Independent of the accelerated division helper and its certificate: inspect the complete remainders.
    template<typename Polynomial>
    std::pair<Polynomial, bool> classical_chain(Polynomial first, Polynomial second) {
        bool normal = first.size() > 1 && first.size() == second.size() + 1;
        if (math::is_zero(first)) {
            return {second, false};
        }
        if (first.size() < second.size()) {
            std::swap(first, second);
        }
        while (!math::is_zero(second)) {
            Polynomial quotient, remainder;
            math::division(quotient, remainder, first, second);
            normal = normal && first.size() == second.size() + 1 &&
                     (math::is_zero(remainder) ? second.size() == 1 : second.size() == remainder.size() + 1);
            first = std::move(second);
            second = std::move(remainder);
        }
        return {first, normal && first.size() == 1 && !math::is_zero(first)};
    }

    // Test-only resultant oracle using a Sylvester determinant, not Euclidean recurrence. The existing static
    // matrix elimination orders pivots by magnitude, which is not available for extension-field elements.
    template<typename Polynomial>
    typename Polynomial::value_type sylvester_resultant(const Polynomial &left, const Polynomial &right) {
        using value_type = typename Polynomial::value_type;
        const std::size_t m = left.size() - 1, n = right.size() - 1, size = m + n;
        std::vector<std::vector<value_type>> matrix(size, std::vector<value_type>(size, value_type::zero()));
        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t j = 0; j <= m; ++j) {
                matrix[row][row + j] = left[m - j];
            }
        }
        for (std::size_t row = 0; row < m; ++row) {
            for (std::size_t j = 0; j <= n; ++j) {
                matrix[n + row][row + j] = right[n - j];
            }
        }
        value_type determinant = value_type::one();
        for (std::size_t column = 0; column < size; ++column) {
            std::size_t pivot = column;
            while (pivot < size && matrix[pivot][column].is_zero()) {
                ++pivot;
            }
            if (pivot == size) {
                return value_type::zero();
            }
            if (pivot != column) {
                std::swap(matrix[pivot], matrix[column]);
                determinant = -determinant;
            }
            determinant *= matrix[column][column];
            const value_type inverse = matrix[column][column].inversed();
            for (std::size_t row = column + 1; row < size; ++row) {
                const value_type scale = matrix[row][column] * inverse;
                for (std::size_t j = column + 1; j < size; ++j) {
                    matrix[row][j] -= scale * matrix[column][j];
                }
            }
        }
        return determinant;
    }

    template<typename Backend>
    void check_chain(const typename Backend::polynomial_type &left, const typename Backend::polynomial_type &right,
                     polynomial_arithmetic::polynomial_context<Backend> &context) {
        using polynomial_type = typename Backend::polynomial_type;
        const auto [last, normal] = classical_chain(left, right);
        polynomial_type expected, actual;
        if (math::is_zero(last)) {
            expected = last;
        } else {
            math::make_monic(expected, last);
        }
        const auto certificate = math::gcd_with_normal_chain_certificate(actual, left, right, context);
        BOOST_CHECK(actual == expected);
        BOOST_REQUIRE_EQUAL(certificate.has_value(), normal);
        if (certificate) {
            BOOST_CHECK(*certificate == last[0]);
        }
        math::gcd(actual, left, right, context);
        BOOST_CHECK(actual == expected);
    }

    // Build a known remainder sequence backwards. A chosen quadratic quotient creates a degree gap at an exact
    // position, including positions that half-GCD visits only through truncated high parts.
    std::pair<math::polynomial<fq_value_type>, math::polynomial<fq_value_type>> constructed_chain(std::size_t divisions,
                                                                                                  std::size_t gap) {
        using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
        using polynomial_type = backend_type::polynomial_type;
        polynomial_type first {fq_value_type(7)}, second {fq_value_type::zero()};
        for (std::size_t step = 1; step <= divisions; ++step) {
            polynomial_type quotient {fq_value_type(step + 2), fq_value_type(step + 3)};
            if (step == gap) {
                quotient.push_back(fq_value_type(step + 4));
            }
            polynomial_type next;
            math::addition(next, product(backend_type {}, quotient, first), second);
            second = std::move(first);
            first = std::move(next);
        }
        return {first, second};
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

BOOST_AUTO_TEST_CASE(small_prime_half_gcd_uses_fused_matrices_and_preserves_the_fallback) {
    using reference_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = reference_backend::polynomial_type;
    polynomial_arithmetic::polynomial_context<reference_backend> euclidean_context;
    polynomial_arithmetic::polynomial_context_options options;
    options.gcd_half_gcd_cutoff = 1;
    options.half_gcd_basecase_cutoff = 8;
    matrix_product_counts fused_counts, fallback_counts;
    polynomial_arithmetic::polynomial_context<counting_matrix_backend> fused_context(
        counting_matrix_backend(1024, fused_counts, 32 * 1024 * 1024), options);
    polynomial_arithmetic::polynomial_context<counting_matrix_backend> fallback_context(
        counting_matrix_backend(1024, fallback_counts, 0), options);

    for (const std::size_t common_size : {1, 5}) {
        const polynomial_type common = fq12_polynomial(common_size, 0x123456789abcdefULL);
        const polynomial_type left = product(reference_backend {}, common, fq12_polynomial(129, 0x31415926));
        const polynomial_type right = product(reference_backend {}, common, fq12_polynomial(128, 0xfedcba98));
        polynomial_type expected, actual;
        math::gcd(expected, left, right, euclidean_context);
        math::gcd(actual, left, right, fused_context);
        BOOST_CHECK(actual == expected);
        math::gcd(actual, left, right, fallback_context);
        BOOST_CHECK(actual == expected);

        auto left_alias = left;
        auto right_alias = right;
        math::gcd(left_alias, left_alias, right, fused_context);
        math::gcd(right_alias, left, right_alias, fused_context);
        BOOST_CHECK(left_alias == expected);
        BOOST_CHECK(right_alias == expected);

        math::detail::half_gcd_matrix<counting_matrix_backend> matrix;
        polynomial_type reduced_first, reduced_second;
        math::detail::half_gcd_reduce(matrix, reduced_first, reduced_second, left, right,
                                      options.half_gcd_basecase_cutoff, fused_context);
        BOOST_CHECK(math::is_zero(reduced_second) || reduced_second.size() <= left.size() / 2);
        polynomial_type expected_first, expected_second;
        math::addition(expected_first, product(reference_backend {}, matrix.entry_00, left),
                       product(reference_backend {}, matrix.entry_01, right));
        math::addition(expected_second, product(reference_backend {}, matrix.entry_10, left),
                       product(reference_backend {}, matrix.entry_11, right));
        BOOST_CHECK(reduced_first == expected_first);
        BOOST_CHECK(reduced_second == expected_second);
    }
    BOOST_CHECK_GT(fused_counts.accepted[0], 0);
    BOOST_CHECK_GT(fused_counts.accepted[1], 0);
    BOOST_CHECK_GT(fused_counts.declined, 0);
    BOOST_CHECK_EQUAL(fallback_counts.accepted[0], 0);
    BOOST_CHECK_EQUAL(fallback_counts.accepted[1], 0);
    BOOST_CHECK_GT(fallback_counts.declined, 0);
}

BOOST_AUTO_TEST_CASE(normal_chain_certificate_matches_classical_remainders_across_recursive_boundaries) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    for (const std::size_t divisions : {1, 2, 7, 8, 9, 31, 32, 33, 63, 64, 65, 127, 128, 129}) {
        for (const std::size_t gap : {std::size_t(0), std::size_t(1), divisions / 2, divisions}) {
            const auto [left, right] = constructed_chain(divisions, gap);
            for (const std::size_t basecase : {0, 4, 8, 60}) {
                polynomial_arithmetic::polynomial_context_options options;
                options.gcd_half_gcd_cutoff = 1;
                options.half_gcd_basecase_cutoff = basecase;
                polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
                BOOST_TEST_CONTEXT("divisions=" << divisions << " gap=" << gap << " basecase=" << basecase) {
                    check_chain(left, right, context);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(normal_chain_certificate_agrees_with_independent_sylvester_resultants) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = backend_type::polynomial_type;
    polynomial_arithmetic::polynomial_context_options options;
    options.gcd_half_gcd_cutoff = 1;
    options.half_gcd_basecase_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
    for (std::size_t degree = 2; degree <= 9; ++degree) {
        for (std::uint64_t seed = 1; seed <= 8; ++seed) {
            const auto left = fq_polynomial(degree + 1, seed);
            const auto right = fq_polynomial(degree, seed + 137);
            const auto expected = sylvester_resultant(left, right);
            polynomial_type common;
            const auto certificate = math::gcd_with_normal_chain_certificate(common, left, right, context);
            BOOST_REQUIRE(certificate.has_value());
            BOOST_CHECK(!expected.is_zero());
            BOOST_CHECK((expected * certificate->inversed()).is_square());
            check_chain(left, right, context);
        }
    }
}

BOOST_AUTO_TEST_CASE(normal_chain_certificate_preserves_the_raw_constant_and_is_alias_safe) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = backend_type::polynomial_type;
    const auto [left, right] = constructed_chain(65, 0);
    for (const std::size_t cutoff : {0, 1, 128}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.gcd_half_gcd_cutoff = cutoff;
        options.half_gcd_basecase_cutoff = 4;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        polynomial_type left_alias = left, right_alias = right;
        const auto left_certificate = math::gcd_with_normal_chain_certificate(left_alias, left_alias, right, context);
        const auto right_certificate = math::gcd_with_normal_chain_certificate(right_alias, left, right_alias, context);
        BOOST_REQUIRE(left_certificate.has_value());
        BOOST_REQUIRE(right_certificate.has_value());
        BOOST_CHECK(*left_certificate == fq_value_type(7));
        BOOST_CHECK(*right_certificate == fq_value_type(7));
        BOOST_CHECK(left_alias == polynomial_type({fq_value_type::one()}));
        BOOST_CHECK(right_alias == left_alias);
    }
}

BOOST_AUTO_TEST_CASE(normal_chain_certificate_declines_degenerate_and_repeated_factor_inputs) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = backend_type::polynomial_type;
    const polynomial_type zero {fq_value_type::zero()}, constant {fq_value_type(3)};
    const polynomial_type linear {fq_value_type(2), fq_value_type::one()};
    const polynomial_type repeated = product(backend_type {}, linear, linear);
    const polynomial_type cubic {fq_value_type::one(), fq_value_type::zero(), fq_value_type::zero(),
                                 fq_value_type::one()};
    polynomial_type repeated_derivative, cubic_derivative;
    math::derivative(repeated_derivative, repeated);
    math::derivative(cubic_derivative, cubic);
    const std::array<std::pair<polynomial_type, polynomial_type>, 10> cases {std::pair {zero, zero},
                                                                             std::pair {zero, linear},
                                                                             std::pair {linear, zero},
                                                                             std::pair {constant, constant},
                                                                             std::pair {linear, linear},
                                                                             std::pair {constant, linear},
                                                                             std::pair {cubic, linear},
                                                                             std::pair {repeated, repeated_derivative},
                                                                             std::pair {cubic, cubic_derivative},
                                                                             std::pair {linear, repeated}};
    for (const std::size_t cutoff : {0, 1}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.gcd_half_gcd_cutoff = cutoff;
        options.half_gcd_basecase_cutoff = 0;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        for (const auto &[left, right] : cases) {
            polynomial_type common;
            BOOST_CHECK(!math::gcd_with_normal_chain_certificate(common, left, right, context).has_value());
            check_chain(left, right, context);
        }
        check_chain(linear, constant, context);
    }
}

BOOST_AUTO_TEST_CASE(irregular_chain_constant_can_give_the_wrong_factor_parity) {
    using value_type = fields::koalabear::value_type;
    using backend_type = polynomial_arithmetic::schoolbook_backend<value_type>;
    using polynomial_type = backend_type::polynomial_type;
    // This field is 1 mod 4 and 2 mod 3. X^3+1 has two irreducible factors; its raw last constant is a square.
    const polynomial_type input {value_type::one(), value_type::zero(), value_type::zero(), value_type::one()};
    polynomial_type derivative, common;
    math::derivative(derivative, input);
    const auto [last, normal] = classical_chain(input, derivative);
    BOOST_CHECK(!normal);
    BOOST_CHECK(last == polynomial_type({value_type::one()}));
    const auto discriminant = -sylvester_resultant(input, derivative);
    BOOST_CHECK(discriminant == -value_type(27));
    BOOST_CHECK(!discriminant.is_square());
    polynomial_arithmetic::polynomial_context_options options;
    options.gcd_half_gcd_cutoff = 1;
    options.half_gcd_basecase_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
    BOOST_CHECK(!math::gcd_with_normal_chain_certificate(common, input, derivative, context).has_value());
    BOOST_CHECK(common == polynomial_type({value_type::one()}));
}

BOOST_AUTO_TEST_CASE(prime_field_factor_parity_needs_the_discriminant_sign) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = backend_type::polynomial_type;
    // (X-1)(X-2) has two factors and discriminant 1, but resultant(F,F') = -1 is not a square in BN254 Fp.
    const polynomial_type input {fq_value_type(2), -fq_value_type(3), fq_value_type::one()};
    polynomial_type derivative, common;
    math::derivative(derivative, input);
    polynomial_arithmetic::polynomial_context<backend_type> context;
    const auto certificate = math::gcd_with_normal_chain_certificate(common, input, derivative, context);
    BOOST_REQUIRE(certificate.has_value());
    BOOST_CHECK(!certificate->is_square());
    BOOST_CHECK((-*certificate).is_square());
    BOOST_CHECK(sylvester_resultant(input, derivative) == -fq_value_type::one());
}

BOOST_AUTO_TEST_CASE(fp12_certificate_recovers_both_factor_count_parities) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = backend_type::polynomial_type;
    fq12_value_type nonsquare = fq12_value_type::zero();
    for (std::size_t i = 1; i <= 32; ++i) {
        const auto candidate = fq12_value(i);
        if (!candidate.is_square()) {
            nonsquare = candidate;
            break;
        }
    }
    BOOST_REQUIRE(!nonsquare.is_zero());
    const polynomial_type irreducible_quadratic {(fq12_value_type::one() - nonsquare) * fq12_value_type(4).inversed(),
                                                 fq12_value_type::one(), fq12_value_type::one()};
    polynomial_arithmetic::polynomial_context_options options;
    options.gcd_half_gcd_cutoff = 1;
    options.half_gcd_basecase_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
    for (const bool use_quadratic : {false, true}) {
        polynomial_type input = use_quadratic ? irreducible_quadratic : polynomial_type {fq12_value_type::one()};
        for (std::size_t linear_factors = 1; linear_factors <= 4; ++linear_factors) {
            input = product(backend_type {}, input,
                            polynomial_type {fq12_value_type(linear_factors + 5), fq12_value_type::one()});
            polynomial_type derivative, common;
            math::derivative(derivative, input);
            const auto certificate = math::gcd_with_normal_chain_certificate(common, input, derivative, context);
            BOOST_REQUIRE(certificate.has_value());
            const auto resultant = sylvester_resultant(input, derivative);
            BOOST_CHECK((resultant * certificate->inversed()).is_square());
            const bool odd = ((input.size() - 1) % 2 != 0) != !certificate->is_square();
            BOOST_CHECK_EQUAL(odd, (linear_factors + std::size_t(use_quadratic)) % 2 != 0);
            check_chain(input, derivative, context);
        }
    }
}

BOOST_AUTO_TEST_CASE(normal_chain_certificate_survives_fused_and_fallback_matrix_products) {
    polynomial_arithmetic::polynomial_context_options options;
    options.gcd_half_gcd_cutoff = 1;
    options.half_gcd_basecase_cutoff = 8;
    matrix_product_counts fused_counts, fallback_counts;
    polynomial_arithmetic::polynomial_context<counting_matrix_backend> fused_context(
        counting_matrix_backend(1024, fused_counts, 32 * 1024 * 1024), options);
    polynomial_arithmetic::polynomial_context<counting_matrix_backend> fallback_context(
        counting_matrix_backend(1024, fallback_counts, 0), options);
    const auto left = fq12_polynomial(130, 0x123456789abcdefULL);
    const auto right = fq12_polynomial(129, 0x31415926);
    check_chain(left, right, fused_context);
    check_chain(left, right, fallback_context);
    BOOST_CHECK_GT(fused_counts.accepted[0], 0);
    BOOST_CHECK_GT(fused_counts.accepted[1], 0);
    BOOST_CHECK_EQUAL(fallback_counts.accepted[0], 0);
    BOOST_CHECK_EQUAL(fallback_counts.accepted[1], 0);
    BOOST_CHECK_GT(fallback_counts.declined, 0);
}

BOOST_AUTO_TEST_SUITE_END()
