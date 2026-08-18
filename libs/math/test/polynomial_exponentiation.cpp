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

#define BOOST_TEST_MODULE polynomial_exponentiation_test

#include <cstddef>
#include <utility>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/polynomial_exponentiation.hpp>
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
    typename Backend::polynomial_type
        repeated_powmod(const typename Backend::polynomial_type &base, std::size_t exponent,
                        const math::polynomial_divisor_context<Backend> &divisor_context,
                        polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        polynomial_type result = {value_type::one()};
        polynomial_type reduced_base;
        math::remainder(reduced_base, base, divisor_context, arithmetic_context);
        for (std::size_t i = 0; i < exponent; ++i) {
            math::mulmod(result, result, reduced_base, divisor_context, arithmetic_context);
        }
        return result;
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(polynomial_exponentiation_test_suite)

BOOST_AUTO_TEST_CASE(squaremod_uses_modular_squaring_and_may_alias_the_input) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type input = {fq_value_type(2), fq_value_type(3), fq_value_type(5), fq_value_type(7)};
    const polynomial_type divisor = {fq_value_type(3), fq_value_type(2), fq_value_type(4)};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 5, arithmetic_context);

    polynomial_type square;
    arithmetic_context.square(square, input);
    polynomial_type unused_quotient;
    polynomial_type expected;
    math::division(unused_quotient, expected, square, divisor);

    polynomial_type result;
    math::squaremod(result, input, divisor_context, arithmetic_context);
    BOOST_CHECK(result == expected);

    polynomial_type alias = input;
    math::squaremod(alias, alias, divisor_context, arithmetic_context);
    BOOST_CHECK(alias == expected);
}

BOOST_AUTO_TEST_CASE(powmod_matches_repeated_multiplication_and_may_alias_the_base) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type base = {fq_value_type(2), fq_value_type(3), fq_value_type(5), fq_value_type(7)};
    const polynomial_type divisor = {fq_value_type(3), fq_value_type(2), fq_value_type(4)};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 2, arithmetic_context);

    for (std::size_t exponent = 0; exponent <= 12; ++exponent) {
        const polynomial_type expected = repeated_powmod(base, exponent, divisor_context, arithmetic_context);
        polynomial_type result;
        math::powmod(result, base, exponent, divisor_context, arithmetic_context);
        BOOST_CHECK(result == expected);
    }

    const polynomial_type expected = repeated_powmod(base, 11, divisor_context, arithmetic_context);
    polynomial_type alias = base;
    math::powmod(alias, alias, std::size_t(11), divisor_context, arithmetic_context);
    BOOST_CHECK(alias == expected);
}

BOOST_AUTO_TEST_CASE(powmod_supports_large_multiprecision_exponents_and_rejects_negative_ones) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using boost::multiprecision::cpp_int;

    const polynomial_type base = {fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);

    const cpp_int exponent = (cpp_int(1) << 130) + 3;
    polynomial_type result;
    math::powmod(result, base, exponent, divisor_context, arithmetic_context);
    const polynomial_type expected = {fq_value_type::zero(), -fq_value_type::one()};
    BOOST_CHECK(result == expected);

    result = {fq_value_type(17)};
    BOOST_CHECK_THROW(math::powmod(result, base, cpp_int(-1), divisor_context, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK(result == polynomial_type({fq_value_type(17)}));
}

BOOST_AUTO_TEST_CASE(powmod_returns_zero_modulo_a_nonzero_constant) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type base = {fq_value_type(2), fq_value_type(3)};
    const polynomial_type divisor = {fq_value_type(7)};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);

    for (const std::size_t exponent : {0, 1, 9}) {
        polynomial_type result;
        math::powmod(result, base, exponent, divisor_context, arithmetic_context);
        BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));
    }
}

BOOST_AUTO_TEST_CASE(powmod_rejects_insufficient_inverse_precision_without_changing_output) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type {}, options);
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::zero(),
                                     fq_value_type::one()};
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);
    const polynomial_type base = {fq_value_type::one(), fq_value_type(2), fq_value_type(3)};

    polynomial_type result = {fq_value_type(17)};
    BOOST_CHECK_THROW(math::powmod(result, base, std::size_t(2), divisor_context, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK(result == polynomial_type({fq_value_type(17)}));
}

BOOST_AUTO_TEST_CASE(powmod_supports_extension_coefficients_and_base_field_roots) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    const polynomial_type base = {fq12_value(1), fq12_value(13), fq12_value(25)};
    const polynomial_type divisor = {fq12_value(37), fq12_value(49), fq12_value(61), fq12_value(73)};

    polynomial_arithmetic::polynomial_context<schoolbook_backend> schoolbook_context;
    math::polynomial_divisor_context<schoolbook_backend> schoolbook_divisor(divisor, 2, schoolbook_context);
    const polynomial_type expected = repeated_powmod(base, 11, schoolbook_divisor, schoolbook_context);

    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<mixed_radix_backend> mixed_radix_context(mixed_radix_backend(18),
                                                                                       options);
    math::polynomial_divisor_context<mixed_radix_backend> mixed_radix_divisor(divisor, 2, mixed_radix_context);
    polynomial_type result;
    math::powmod(result, base, std::size_t(11), mixed_radix_divisor, mixed_radix_context);
    BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_SUITE_END()
