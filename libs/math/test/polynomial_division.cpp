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

#define BOOST_TEST_MODULE polynomial_division_test

#include <cstddef>
#include <utility>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/arithmetic/polynomial_division.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_modular_arithmetic.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = nil::crypto3::math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq_value_type = fq_field_type::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;

    template<polynomial_arithmetic::PolynomialBackend Backend>
    typename Backend::polynomial_type build_and_check_context(Backend backend,
                                                              const typename Backend::polynomial_type &divisor,
                                                              std::size_t inverse_precision) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        polynomial_arithmetic::polynomial_context<Backend> arithmetic_context(std::move(backend));
        math::polynomial_divisor_context<Backend> divisor_context(divisor, inverse_precision, arithmetic_context);

        polynomial_type canonical_divisor(divisor);
        math::condense(canonical_divisor);
        BOOST_CHECK(divisor_context.divisor() == canonical_divisor);
        BOOST_CHECK_EQUAL(divisor_context.degree(), canonical_divisor.size() - 1);
        BOOST_CHECK_EQUAL(divisor_context.inverse_precision(), inverse_precision);

        polynomial_type reversed_divisor(canonical_divisor);
        math::reverse(reversed_divisor, reversed_divisor.size());
        polynomial_type product;
        math::multiply_low(product, reversed_divisor, divisor_context.reversed_divisor_inverse(), inverse_precision,
                           arithmetic_context);
        BOOST_CHECK(product == polynomial_type({value_type::one()}));
        return divisor_context.reversed_divisor_inverse();
    }

    template<polynomial_arithmetic::PolynomialBackend Backend>
    std::pair<typename Backend::polynomial_type, typename Backend::polynomial_type>
        compute_divrem(Backend backend, const typename Backend::polynomial_type &dividend,
                       const typename Backend::polynomial_type &divisor, std::size_t inverse_precision) {
        using polynomial_type = typename Backend::polynomial_type;

        polynomial_arithmetic::polynomial_context_options options;
        options.basecase_divisor_coefficient_cutoff = 0;
        options.basecase_quotient_coefficient_cutoff = 0;
        polynomial_arithmetic::polynomial_context<Backend> arithmetic_context(std::move(backend), options);
        math::polynomial_divisor_context<Backend> divisor_context(divisor, inverse_precision, arithmetic_context);
        polynomial_type quotient;
        polynomial_type remainder;
        math::divrem(quotient, remainder, dividend, divisor_context, arithmetic_context);
        return {std::move(quotient), std::move(remainder)};
    }

    fq12_value_type fq12_value(std::size_t first_coordinate) {
        fq12_value_type value = fq12_value_type::zero();
        for (std::size_t i = 0; i < fq12_field_type::arity; ++i) {
            value.coordinate(i) = fq_value_type(first_coordinate + i);
        }
        return value;
    }

    struct prepared_product_counts {
        std::size_t preparations = 0;
        std::size_t prepared_products = 0;
        std::size_t ordinary_products = 0;
        std::size_t squares = 0;
    };

    // Count calls at the backend boundary, as in the exponentiation tests; all arithmetic uses the real backend.
    template<typename ValueType = fq_value_type>
    struct counting_mixed_radix_backend : polynomial_arithmetic::mixed_radix_backend<fq_field_type, ValueType> {
        using base_type = polynomial_arithmetic::mixed_radix_backend<fq_field_type, ValueType>;
        using polynomial_type = typename base_type::polynomial_type;
        using prepared_type = typename base_type::prepared_low_product_type;

        counting_mixed_radix_backend(std::size_t order, prepared_product_counts &counts) :
            base_type(order), counts(counts) {
        }

        auto prepare_low_product(const polynomial_type &fixed, std::size_t variable_count, std::size_t count) {
            ++counts.preparations;
            return base_type::prepare_low_product(fixed, variable_count, count);
        }

        bool try_multiply_low_prepared(polynomial_type &output, const polynomial_type &input,
                                       const prepared_type &prepared, std::size_t count) {
            const bool used = base_type::try_multiply_low_prepared(output, input, prepared, count);
            counts.prepared_products += used;
            return used;
        }

        void multiply_low(polynomial_type &output, const polynomial_type &left, const polynomial_type &right,
                          std::size_t count) {
            ++counts.ordinary_products;
            base_type::multiply_low(output, left, right, count);
        }

        void square(polynomial_type &output, const polynomial_type &input) {
            ++counts.squares;
            base_type::square(output, input);
        }

        prepared_product_counts &counts;
    };

    template<typename Backend>
    void check_reused_division(const typename Backend::polynomial_type &dividend,
                               const math::polynomial_divisor_context<Backend> &divisor,
                               polynomial_arithmetic::polynomial_context<Backend> &context) {
        using polynomial_type = typename Backend::polynomial_type;
        polynomial_type expected_quotient, expected_remainder, quotient, remainder, reconstructed;
        math::division(expected_quotient, expected_remainder, dividend, divisor.divisor());
        math::divrem(quotient, remainder, dividend, divisor, context);
        BOOST_CHECK(quotient == expected_quotient);
        BOOST_CHECK(remainder == expected_remainder);
        BOOST_CHECK(remainder.size() < divisor.divisor().size());
        polynomial_arithmetic::schoolbook_backend<typename Backend::value_type> {}.multiply(reconstructed, quotient,
                                                                                            divisor.divisor());
        math::addition(reconstructed, reconstructed, remainder);
        BOOST_CHECK(reconstructed == dividend);
    }

    template<typename ValueType>
    void check_cached_squaremod(const math::polynomial<ValueType> &divisor, math::polynomial<ValueType> input) {
        using backend_type = counting_mixed_radix_backend<ValueType>;
        using polynomial_type = typename backend_type::polynomial_type;
        prepared_product_counts counts;
        polynomial_arithmetic::polynomial_context<backend_type> context {backend_type(58, counts)};
        math::polynomial_divisor_context<backend_type> fixed(divisor, divisor.size() - 2, context);
        BOOST_CHECK_EQUAL(counts.preparations, 2);

        for (std::size_t repeat = 0; repeat < 3; ++repeat) {
            polynomial_type square, quotient, expected, result;
            polynomial_arithmetic::schoolbook_backend<ValueType> {}.square(square, input);
            math::division(quotient, expected, square, divisor);
            counts = {};
            math::squaremod(result, input, fixed, context);
            BOOST_CHECK(result == expected);
            BOOST_CHECK_EQUAL(counts.squares, 1);
            BOOST_CHECK_EQUAL(counts.prepared_products, 2);
            BOOST_CHECK_EQUAL(counts.ordinary_products, 0);
            BOOST_CHECK_EQUAL(counts.preparations, 0);
            // All three operations select order 29 for these degree-11 moduli. One square and two prepared
            // products therefore use six transforms; ordinary reduction would use eight.
            polynomial_type alias = input;
            math::squaremod(alias, alias, fixed, context);
            BOOST_CHECK(alias == expected);
            check_reused_division(square, fixed, context);
            input[0] += ValueType::one();
        }
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(polynomial_division_test_suite)

BOOST_AUTO_TEST_CASE(divisor_is_canonical_and_reversed_inverse_is_correct) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type(4), fq_value_type::zero()};
    build_and_check_context(backend_type {}, divisor, 7);
}

BOOST_AUTO_TEST_CASE(nonzero_constant_divisor_is_supported) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq_value_type(7)};
    build_and_check_context(backend_type {}, divisor, 5);
}

BOOST_AUTO_TEST_CASE(schoolbook_and_mixed_radix_precomputation_agree) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    constexpr std::size_t inverse_precision = 9;
    // The final Newton step multiplies prefixes of lengths 8 and 9. A size-18 transform is the smallest supported
    // mixed-radix transform covering their 16-coefficient product.
    constexpr std::size_t transform_size = 18;
    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type(5), fq_value_type(7)};
    const polynomial_type schoolbook_inverse =
        build_and_check_context(schoolbook_backend {}, divisor, inverse_precision);
    const polynomial_type mixed_radix_inverse =
        build_and_check_context(mixed_radix_backend(transform_size), divisor, inverse_precision);
    BOOST_CHECK(mixed_radix_inverse == schoolbook_inverse);
}

BOOST_AUTO_TEST_CASE(extension_field_divisor_is_supported) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    constexpr std::size_t inverse_precision = 7;
    // The final Newton step multiplies prefixes of lengths 4 and 6, whose product has 9 coefficients.
    constexpr std::size_t transform_size = 9;
    const polynomial_type divisor = {fq12_value(1), fq12_value(13), fq12_value(25)};
    const polynomial_type schoolbook_inverse =
        build_and_check_context(schoolbook_backend {}, divisor, inverse_precision);
    const polynomial_type mixed_radix_inverse =
        build_and_check_context(mixed_radix_backend(transform_size), divisor, inverse_precision);
    BOOST_CHECK(mixed_radix_inverse == schoolbook_inverse);
}

BOOST_AUTO_TEST_CASE(zero_divisor_is_rejected) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type zero = {fq_value_type::zero(), fq_value_type::zero()};
    BOOST_CHECK_THROW(math::polynomial_divisor_context<backend_type>(zero, 3, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(zero_inverse_precision_is_rejected) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::one()};
    BOOST_CHECK_THROW(math::polynomial_divisor_context<backend_type>(divisor, 0, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(fast_divrem_recovers_known_quotient_and_remainder) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type dividend = {fq_value_type(7), fq_value_type(11), fq_value_type(9), fq_value_type(7),
                                      fq_value_type(4)};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::one(), fq_value_type::one()};
    const polynomial_type expected_quotient = {fq_value_type(2), fq_value_type(3), fq_value_type(4)};
    const polynomial_type expected_remainder = {fq_value_type(5), fq_value_type(6)};

    const auto [quotient, remainder] = compute_divrem(backend_type {}, dividend, divisor, 3);
    BOOST_CHECK(quotient == expected_quotient);
    BOOST_CHECK(remainder == expected_remainder);
}

BOOST_AUTO_TEST_CASE(fast_divrem_preserves_zero_low_quotient_coefficients) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type dividend = {fq_value_type(5), fq_value_type::zero(), fq_value_type::one(),
                                      fq_value_type::one()};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type expected_quotient = {fq_value_type::zero(), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type expected_remainder = {fq_value_type(5)};

    const auto [quotient, remainder] = compute_divrem(backend_type {}, dividend, divisor, 3);
    BOOST_CHECK(quotient == expected_quotient);
    BOOST_CHECK(remainder == expected_remainder);
}

BOOST_AUTO_TEST_CASE(fast_divrem_handles_small_dividends_and_constant_divisors) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type small_dividend = {fq_value_type(4), fq_value_type(5)};
    const polynomial_type larger_divisor = {fq_value_type::one(), fq_value_type(2), fq_value_type(3)};
    const auto [zero_quotient, unchanged_remainder] =
        compute_divrem(backend_type {}, small_dividend, larger_divisor, 1);
    BOOST_CHECK(zero_quotient == polynomial_type({fq_value_type::zero()}));
    BOOST_CHECK(unchanged_remainder == small_dividend);

    const polynomial_type dividend = {fq_value_type(2), fq_value_type(4), fq_value_type(6)};
    const polynomial_type constant_divisor = {fq_value_type(2)};
    const auto [quotient, zero_remainder] = compute_divrem(backend_type {}, dividend, constant_divisor, 3);
    BOOST_CHECK(quotient == polynomial_type({fq_value_type::one(), fq_value_type(2), fq_value_type(3)}));
    BOOST_CHECK(zero_remainder == polynomial_type({fq_value_type::zero()}));
}

BOOST_AUTO_TEST_CASE(fast_divrem_outputs_may_alias_the_dividend) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type dividend = {fq_value_type(7), fq_value_type(11), fq_value_type(9), fq_value_type(7),
                                      fq_value_type(4)};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::one(), fq_value_type::one()};
    const polynomial_type expected_quotient = {fq_value_type(2), fq_value_type(3), fq_value_type(4)};
    const polynomial_type expected_remainder = {fq_value_type(5), fq_value_type(6)};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 3, arithmetic_context);

    polynomial_type quotient_alias = dividend;
    polynomial_type remainder;
    math::divrem(quotient_alias, remainder, quotient_alias, divisor_context, arithmetic_context);
    BOOST_CHECK(quotient_alias == expected_quotient);
    BOOST_CHECK(remainder == expected_remainder);

    polynomial_type quotient;
    polynomial_type remainder_alias = dividend;
    math::divrem(quotient, remainder_alias, remainder_alias, divisor_context, arithmetic_context);
    BOOST_CHECK(quotient == expected_quotient);
    BOOST_CHECK(remainder_alias == expected_remainder);
}

BOOST_AUTO_TEST_CASE(fast_remainder_discards_the_quotient_and_may_alias_the_dividend) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type dividend = {fq_value_type(7), fq_value_type(11), fq_value_type(9), fq_value_type(7),
                                      fq_value_type(4)};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::one(), fq_value_type::one()};
    const polynomial_type expected_remainder = {fq_value_type(5), fq_value_type(6)};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 3, arithmetic_context);

    polynomial_type result;
    math::remainder(result, dividend, divisor_context, arithmetic_context);
    BOOST_CHECK(result == expected_remainder);

    polynomial_type aliased_result = dividend;
    math::remainder(aliased_result, aliased_result, divisor_context, arithmetic_context);
    BOOST_CHECK(aliased_result == expected_remainder);

    const polynomial_type small_dividend = {fq_value_type(2), fq_value_type(3)};
    math::remainder(result, small_dividend, divisor_context, arithmetic_context);
    BOOST_CHECK(result == small_dividend);
}

BOOST_AUTO_TEST_CASE(fast_exact_division_rejects_a_nonzero_remainder_and_may_alias_the_dividend) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::one(), fq_value_type::one()};
    const polynomial_type expected_quotient = {fq_value_type(2), fq_value_type(3), fq_value_type(4)};
    const polynomial_type exact_dividend = {fq_value_type(2), fq_value_type(5), fq_value_type(9), fq_value_type(7),
                                            fq_value_type(4)};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 3, arithmetic_context);

    polynomial_type result;
    math::exact_division(result, exact_dividend, divisor_context, arithmetic_context);
    BOOST_CHECK(result == expected_quotient);

    polynomial_type aliased_result = exact_dividend;
    math::exact_division(aliased_result, aliased_result, divisor_context, arithmetic_context);
    BOOST_CHECK(aliased_result == expected_quotient);

    const polynomial_type inexact_dividend = {fq_value_type(7), fq_value_type(11), fq_value_type(9), fq_value_type(7),
                                              fq_value_type(4)};
    result = {fq_value_type(17)};
    BOOST_CHECK_THROW(math::exact_division(result, inexact_dividend, divisor_context, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK(result == polynomial_type({fq_value_type(17)}));

    polynomial_type inexact_alias = inexact_dividend;
    BOOST_CHECK_THROW(math::exact_division(inexact_alias, inexact_alias, divisor_context, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK(inexact_alias == inexact_dividend);
}

BOOST_AUTO_TEST_CASE(fast_mulmod_multiplies_reduces_and_may_alias_either_input) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type left = {fq_value_type::one(), fq_value_type(2), fq_value_type(3)};
    const polynomial_type right = {fq_value_type(4), fq_value_type(5), fq_value_type(6)};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::one(), fq_value_type::one()};
    const polynomial_type expected_remainder = {fq_value_type(3), fq_value_type(3)};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 3, arithmetic_context);

    polynomial_type result;
    math::mulmod(result, left, right, divisor_context, arithmetic_context);
    BOOST_CHECK(result == expected_remainder);

    polynomial_type left_alias = left;
    math::mulmod(left_alias, left_alias, right, divisor_context, arithmetic_context);
    BOOST_CHECK(left_alias == expected_remainder);

    polynomial_type right_alias = right;
    math::mulmod(right_alias, left, right_alias, divisor_context, arithmetic_context);
    BOOST_CHECK(right_alias == expected_remainder);

    const polynomial_type constant_divisor = {fq_value_type(7)};
    math::polynomial_divisor_context<backend_type> constant_context(constant_divisor, 1, arithmetic_context);
    math::mulmod(result, left, right, constant_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));
}

BOOST_AUTO_TEST_CASE(fast_mulmod_supports_extension_coefficients_and_base_field_roots) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    const polynomial_type left = {fq12_value(1), fq12_value(13), fq12_value(25)};
    const polynomial_type right = {fq12_value(37), fq12_value(49), fq12_value(61)};
    const polynomial_type divisor = {fq12_value(73), fq12_value(85), fq12_value(97)};

    polynomial_type product;
    schoolbook_backend {}.multiply(product, left, right);
    const auto expected = compute_divrem(schoolbook_backend {}, product, divisor, 3);

    constexpr std::size_t transform_size = 9;
    polynomial_arithmetic::polynomial_context<mixed_radix_backend> arithmetic_context {
        mixed_radix_backend(transform_size)};
    math::polynomial_divisor_context<mixed_radix_backend> divisor_context(divisor, 3, arithmetic_context);
    polynomial_type result;
    math::mulmod(result, left, right, divisor_context, arithmetic_context);
    BOOST_CHECK(result == expected.second);
}

BOOST_AUTO_TEST_CASE(fast_divrem_rejects_shared_outputs_and_insufficient_precision) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type dividend = {fq_value_type(7), fq_value_type(11), fq_value_type(9), fq_value_type(7),
                                      fq_value_type(4)};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::one(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type {}, options);

    math::polynomial_divisor_context<backend_type> precise_context(divisor, 3, arithmetic_context);
    polynomial_type shared_output;
    BOOST_CHECK_THROW(math::divrem(shared_output, shared_output, dividend, precise_context, arithmetic_context),
                      std::invalid_argument);

    math::polynomial_divisor_context<backend_type> imprecise_context(divisor, 2, arithmetic_context);
    polynomial_type quotient;
    polynomial_type remainder;
    BOOST_CHECK_THROW(math::divrem(quotient, remainder, dividend, imprecise_context, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(divrem_selects_the_configurable_basecase_fallback) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type small_divisor = {fq_value_type::one(), fq_value_type::one(), fq_value_type::one()};
    const polynomial_type dividend = {fq_value_type(7), fq_value_type(11), fq_value_type(9), fq_value_type(7),
                                      fq_value_type(4)};
    const polynomial_type expected_quotient = {fq_value_type(2), fq_value_type(3), fq_value_type(4)};
    const polynomial_type expected_remainder = {fq_value_type(5), fq_value_type(6)};

    polynomial_arithmetic::polynomial_context<backend_type> default_context;
    math::polynomial_divisor_context<backend_type> small_divisor_context(small_divisor, 1, default_context);
    polynomial_type quotient;
    polynomial_type remainder;
    math::divrem(quotient, remainder, dividend, small_divisor_context, default_context);
    BOOST_CHECK(quotient == expected_quotient);
    BOOST_CHECK(remainder == expected_remainder);

    const polynomial_type large_divisor = {fq_value_type(1), fq_value_type(2),  fq_value_type(3), fq_value_type(4),
                                           fq_value_type(5), fq_value_type(6),  fq_value_type(7), fq_value_type(8),
                                           fq_value_type(9), fq_value_type(10), fq_value_type(11)};
    const polynomial_type short_quotient = {fq_value_type(12), fq_value_type(13)};
    const polynomial_type short_remainder = {fq_value_type(14), fq_value_type(15)};
    polynomial_type large_dividend;
    backend_type {}.multiply(large_dividend, large_divisor, short_quotient);
    math::addition(large_dividend, large_dividend, short_remainder);

    math::polynomial_divisor_context<backend_type> large_divisor_context(large_divisor, 1, default_context);
    math::divrem(quotient, remainder, large_dividend, large_divisor_context, default_context);
    BOOST_CHECK(quotient == short_quotient);
    BOOST_CHECK(remainder == short_remainder);

    const polynomial_type long_quotient = {fq_value_type(12), fq_value_type(13), fq_value_type(14)};
    backend_type {}.multiply(large_dividend, large_divisor, long_quotient);
    math::polynomial_divisor_context<backend_type> default_insufficient_context(large_divisor, 2, default_context);
    BOOST_CHECK_THROW(math::divrem(quotient, remainder, large_dividend, default_insufficient_context, default_context),
                      std::invalid_argument);

    polynomial_arithmetic::polynomial_context_options no_basecase_options;
    no_basecase_options.basecase_divisor_coefficient_cutoff = 0;
    no_basecase_options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> no_basecase_context(backend_type {}, no_basecase_options);
    math::polynomial_divisor_context<backend_type> insufficient_context(small_divisor, 1, no_basecase_context);
    BOOST_CHECK_THROW(math::divrem(quotient, remainder, dividend, insufficient_context, no_basecase_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(fast_divrem_supports_extension_coefficients_and_base_field_roots) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    const polynomial_type divisor = {fq12_value(1), fq12_value(13), fq12_value(25)};
    const polynomial_type expected_quotient = {fq12_value(37), fq12_value(49), fq12_value(61), fq12_value(73)};
    const polynomial_type expected_remainder = {fq12_value(85), fq12_value(97)};
    polynomial_type product;
    polynomial_type dividend;
    schoolbook_backend {}.multiply(product, expected_quotient, divisor);
    math::addition(dividend, product, expected_remainder);

    const auto [schoolbook_quotient, schoolbook_remainder] =
        compute_divrem(schoolbook_backend {}, dividend, divisor, 4);
    constexpr std::size_t transform_size = 9;
    const auto [mixed_radix_quotient, mixed_radix_remainder] =
        compute_divrem(mixed_radix_backend(transform_size), dividend, divisor, 4);

    BOOST_CHECK(schoolbook_quotient == expected_quotient);
    BOOST_CHECK(schoolbook_remainder == expected_remainder);
    BOOST_CHECK(mixed_radix_quotient == expected_quotient);
    BOOST_CHECK(mixed_radix_remainder == expected_remainder);
}

BOOST_AUTO_TEST_CASE(squaremod_reuses_fixed_products_for_prime_and_extension_coefficients) {
    using polynomial_type = math::polynomial<fq_value_type>;
    polynomial_type factor(11, fq_value_type::one()), input(11, fq_value_type::one()), divisor;
    for (std::size_t i = 0; i < 10; ++i) {
        factor[i] = fq_value_type(i + 2);
        input[i] = fq_value_type(2 * i + 3);
    }
    // A reducible degree-11 modulus, first monic and then nonmonic.
    polynomial_arithmetic::schoolbook_backend<fq_value_type> {}.multiply(divisor, factor, polynomial_type {1, 1});
    check_cached_squaremod(divisor, input);
    for (auto &coefficient : divisor) {
        coefficient *= fq_value_type(7);
    }
    check_cached_squaremod(divisor, input);

    math::polynomial<fq12_value_type> extension_divisor(12, fq12_value_type::one());
    math::polynomial<fq12_value_type> extension_input(11, fq12_value_type::one());
    for (std::size_t i = 0; i < extension_input.size(); ++i) {
        extension_divisor[i] = fq12_value(12 * i + 1);
        extension_input[i] = fq12_value(12 * i + 145);
    }
    check_cached_squaremod(extension_divisor, extension_input);
    extension_divisor.back() = fq12_value(289);
    check_cached_squaremod(extension_divisor, extension_input);
}

BOOST_AUTO_TEST_CASE(prepared_division_falls_back_for_changed_quotient_lengths_without_growing_the_cache) {
    using backend_type = counting_mixed_radix_backend<>;
    using polynomial_type = typename backend_type::polynomial_type;
    const polynomial_type divisor {2, 3, 5, 7, 11, 13, 1};
    prepared_product_counts counts;
    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type(18, counts), options);
    const math::polynomial_divisor_context<backend_type> fixed(divisor, 5, context);
    BOOST_CHECK_EQUAL(counts.preparations, 2);

    for (const auto &remainder : {polynomial_type {fq_value_type::zero()}, polynomial_type {19, 23}}) {
        // Quotient size five uses both spectra. Shorter quotients have a different inverse precision and a
        // smaller remainder-product transform. Returning to size five must reuse the original preparations.
        for (std::size_t size : {5, 4, 2, 1, 5}) {
            polynomial_type quotient(size, fq_value_type::zero()), dividend;
            quotient.back() = fq_value_type(17);
            polynomial_arithmetic::schoolbook_backend<fq_value_type> {}.multiply(dividend, quotient, divisor);
            math::addition(dividend, dividend, remainder);
            counts = {};
            check_reused_division(dividend, fixed, context);
            BOOST_CHECK_EQUAL(counts.prepared_products, size == 5 ? 2 : 0);
            BOOST_CHECK_EQUAL(counts.ordinary_products, size == 5 ? 0 : 2);
            BOOST_CHECK_EQUAL(counts.preparations, 0);

            polynomial_type expected_quotient, expected_remainder, aliased = dividend;
            math::divrem(aliased, expected_remainder, aliased, fixed, context);
            BOOST_CHECK(aliased == quotient);
            BOOST_CHECK(expected_remainder == remainder);
            aliased = dividend;
            math::divrem(expected_quotient, aliased, aliased, fixed, context);
            BOOST_CHECK(expected_quotient == quotient);
            BOOST_CHECK(aliased == remainder);
        }
    }
    for (const auto &small : {polynomial_type {fq_value_type::zero()}, polynomial_type {0, 0, 1}}) {
        counts = {};
        check_reused_division(small, fixed, context);
        BOOST_CHECK_EQUAL(counts.prepared_products + counts.ordinary_products + counts.preparations, 0);
    }
}

BOOST_AUTO_TEST_CASE(prepared_divisor_snapshots_survive_copy_move_and_arithmetic_context_changes) {
    using backend_type = counting_mixed_radix_backend<>;
    using polynomial_type = typename backend_type::polynomial_type;
    const polynomial_type original {2, 3, 5, 7, 11, 13, 1};
    const polynomial_type dividend {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    prepared_product_counts counts;
    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    const auto fixed = [&] {
        polynomial_arithmetic::polynomial_context<backend_type> temporary(backend_type(18, counts), options);
        polynomial_type source = original;
        math::polynomial_divisor_context<backend_type> result(source, 5, temporary);
        source[0] = fq_value_type(99);
        return result;
    }();
    auto copy = fixed;
    const auto moved = std::move(copy);
    for (std::size_t order : {522, 67}) {
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type(order, counts), options);
        for (const auto *snapshot : {&fixed, &moved}) {
            counts = {};
            check_reused_division(dividend, *snapshot, context);
            BOOST_CHECK(snapshot->divisor() == original);
            BOOST_CHECK_EQUAL(counts.prepared_products, order == 522 ? 2 : 0);
            BOOST_CHECK_EQUAL(counts.ordinary_products, order == 522 ? 0 : 2);
            BOOST_CHECK_EQUAL(counts.preparations, 0);
        }
        polynomial_type different = original;
        different[0] = fq_value_type(99);
        const math::polynomial_divisor_context<backend_type> changed(different, 5, context);
        check_reused_division(dividend, changed, context);
        check_reused_division(dividend, fixed, context);
    }
}

BOOST_AUTO_TEST_CASE(prepared_division_handles_sparse_reversals_and_zero_fixed_prefixes) {
    using backend_type = counting_mixed_radix_backend<>;
    using polynomial_type = typename backend_type::polynomial_type;
    prepared_product_counts counts;
    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type(18, counts), options);
    // Reversing the top five coefficients of X^10 + X^4 produces [1, 0, 0, 0, 0]. Preserve this padded
    // prefix length when selecting the quotient FFT. The inverse is constant for both moduli.
    const polynomial_type dividend {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1};
    for (const auto &divisor : {polynomial_type {1, 0, 0, 0, 0, 0, 1}, polynomial_type {0, 0, 0, 0, 0, 0, 1}}) {
        const math::polynomial_divisor_context<backend_type> fixed(divisor, 5, context);
        counts = {};
        check_reused_division(dividend, fixed, context);
        BOOST_CHECK_EQUAL(counts.prepared_products, 2);
        BOOST_CHECK_EQUAL(counts.ordinary_products, 0);
        BOOST_CHECK_EQUAL(counts.preparations, 0);
    }
}

BOOST_AUTO_TEST_CASE(unavailable_or_unused_preparations_preserve_existing_division_behavior) {
    using backend_type = counting_mixed_radix_backend<>;
    using polynomial_type = typename backend_type::polynomial_type;
    const polynomial_type divisor {2, 3, 5};
    const polynomial_type dividend {7, 11, 13, 17, 19};
    prepared_product_counts counts;
    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type(9, counts), options);
    // Newton inversion to precision seven fits order nine, but preparing its full length-seven inverse product
    // would need order thirteen or larger. The optional preparation must not make construction fail.
    const math::polynomial_divisor_context<backend_type> partial(divisor, 7, context);
    BOOST_CHECK_EQUAL(counts.preparations, 2);
    counts = {};
    check_reused_division(dividend, partial, context);
    BOOST_CHECK_EQUAL(counts.ordinary_products, 1);
    BOOST_CHECK_EQUAL(counts.prepared_products, 1);
    BOOST_CHECK_EQUAL(counts.preparations, 0);

    polynomial_arithmetic::polynomial_context<backend_type> basecase(backend_type(9, counts));
    counts = {};
    const math::polynomial_divisor_context<backend_type> unprepared(divisor, 3, basecase);
    BOOST_CHECK_EQUAL(counts.preparations, 0);
    counts = {};
    check_reused_division(dividend, unprepared, basecase);
    BOOST_CHECK_EQUAL(counts.ordinary_products + counts.prepared_products, 0);
    // A context constructed under different cutoffs still works when Newton division is subsequently requested.
    check_reused_division(dividend, unprepared, context);
    BOOST_CHECK_EQUAL(counts.ordinary_products, 2);
    BOOST_CHECK_EQUAL(counts.preparations, 0);
    counts = {};
    const math::polynomial_divisor_context<backend_type> constant(polynomial_type {fq_value_type(7)}, 3, context);
    BOOST_CHECK_EQUAL(counts.preparations, 0);
}

BOOST_AUTO_TEST_SUITE_END()
