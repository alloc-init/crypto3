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

#include <nil/crypto3/math/polynomial/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/polynomial_division.hpp>
#include <nil/crypto3/math/polynomial/schoolbook_backend.hpp>

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

        polynomial_arithmetic::polynomial_context<Backend> arithmetic_context(std::move(backend));
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

BOOST_AUTO_TEST_CASE(fast_divrem_rejects_shared_outputs_and_insufficient_precision) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type dividend = {fq_value_type(7), fq_value_type(11), fq_value_type(9), fq_value_type(7),
                                      fq_value_type(4)};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::one(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;

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

BOOST_AUTO_TEST_SUITE_END()
