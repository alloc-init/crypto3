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

#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/field_order.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_exponentiation.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>

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

    // Count calls at the backend boundary; the schoolbook implementation still performs the actual arithmetic.
    struct counting_schoolbook_backend : polynomial_arithmetic::schoolbook_backend<fq_value_type> {
        using base_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;

        counting_schoolbook_backend(std::size_t &squares, std::size_t &products) :
            squares(squares), products(products) {
        }

        void square(polynomial_type &output, const polynomial_type &input) const {
            ++squares;
            base_type::square(output, input);
        }

        void multiply(polynomial_type &output, const polynomial_type &left, const polynomial_type &right) const {
            ++products;
            base_type::multiply(output, left, right);
        }

        std::size_t &squares;
        std::size_t &products;
    };
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

BOOST_AUTO_TEST_CASE(powmod_accepts_an_extension_field_order_as_its_exponent) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type base = {fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);

    polynomial_type result;
    math::powmod(result, base, fields::field_order<fq12_field_type>(), divisor_context, arithmetic_context);

    // X^2 = -1 modulo X^2 + 1. The twelfth power of an odd characteristic is 1 modulo 4, so X^(q^12) = X.
    BOOST_CHECK(result == base);
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

BOOST_AUTO_TEST_CASE(powmod_x_matches_generic_exponentiation_across_moduli_and_reduction_paths) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const auto zero = fq_value_type::zero();
    const auto one = fq_value_type::one();
    const polynomial_type x = {zero, one};
    const polynomial_type divisors[] = {
        {fq_value_type(7)},                                             // Nonzero constant: the zero ring.
        {zero, one},                                                    // X: the initial residue is zero.
        {fq_value_type(3), fq_value_type(2)},                           // Nonmonic linear modulus.
        {one, zero, one},                                               // Monic quadratic modulus.
        {-one, zero, one},                                              // X^2 - 1 is reducible.
        {zero, zero, zero, one},                                        // X^3: positive powers can become zero.
        {fq_value_type(3), fq_value_type(2), fq_value_type(4)},         // Nonmonic quadratic modulus.
        {fq_value_type(3), fq_value_type(2), zero, fq_value_type(4)}    // X^3 reduces to degree one.
    };

    for (const bool force_newton : {false, true}) {
        polynomial_arithmetic::polynomial_context_options options;
        if (force_newton) {
            options.basecase_divisor_coefficient_cutoff = 0;
            options.basecase_quotient_coefficient_cutoff = 0;
        }
        polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type {}, options);
        for (const polynomial_type &divisor : divisors) {
            const std::size_t precision = divisor.size() > 2 ? divisor.size() - 2 : 1;
            math::polynomial_divisor_context<backend_type> divisor_context(divisor, precision, arithmetic_context);
            for (const std::size_t exponent : {0, 1, 2, 3, 4, 5, 7, 8, 9, 11, 16, 21, 31}) {
                BOOST_TEST_CONTEXT("degree = " << divisor.degree() << ", exponent = " << exponent
                                               << ", Newton = " << force_newton) {
                    polynomial_type expected;
                    math::powmod(expected, x, exponent, divisor_context, arithmetic_context);
                    polynomial_type result = {fq_value_type(17), fq_value_type(19)};
                    math::powmod_x(result, exponent, divisor_context, arithmetic_context);
                    BOOST_CHECK(result == expected);
                    BOOST_REQUIRE(!result.empty());
                    BOOST_CHECK(result.size() == 1 || !result.back().is_zero());
                    BOOST_CHECK(result.size() <= (divisor_context.degree() == 0 ? 1 : divisor_context.degree()));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(powmod_x_preserves_high_bits_of_multiprecision_and_native_exponents) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using boost::multiprecision::cpp_int;
    using boost::multiprecision::uint256_t;

    const auto zero = fq_value_type::zero();
    const auto one = fq_value_type::one();
    const polynomial_type x = {zero, one};
    const polynomial_type divisor = {-one, zero, zero, one};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 2, arithmetic_context);

    const auto check_exponent = [&](const auto &exponent, const polynomial_type &expected) {
        polynomial_type reference;
        math::powmod(reference, x, exponent, divisor_context, arithmetic_context);
        polynomial_type result;
        math::powmod_x(result, exponent, divisor_context, arithmetic_context);
        BOOST_CHECK(result == reference);
        BOOST_CHECK(result == expected);
    };

    // X^3 = 1. The high bits change the exponent modulo 3, so narrowing to the low machine word gives a wrong result.
    const cpp_int dynamic_exponent = (cpp_int(1) << 130) + 3;      // 1 modulo 3; the low word is 0 modulo 3.
    const uint256_t fixed_exponent = (uint256_t(1) << 200) + 5;    // 0 modulo 3; the low word is 2 modulo 3.
    check_exponent(dynamic_exponent, x);
    check_exponent(fixed_exponent, polynomial_type {one});
    check_exponent(std::uint64_t(1) << 63, polynomial_type {zero, zero, one});
    check_exponent(std::numeric_limits<std::uint64_t>::max(), polynomial_type {one});
}

BOOST_AUTO_TEST_CASE(powmod_x_accepts_an_extension_field_order_as_its_exponent) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type x = {fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);

    polynomial_type result;
    const auto exponent = fields::field_order<fq12_field_type>();
    math::powmod_x(result, exponent, divisor_context, arithmetic_context);
    polynomial_type reference;
    math::powmod(reference, x, exponent, divisor_context, arithmetic_context);
    BOOST_CHECK(result == reference);
    BOOST_CHECK(result == x);    // X^2 = -1 and the odd characteristic raised to the twelfth power is 1 modulo 4.
}

BOOST_AUTO_TEST_CASE(powmod_x_rejects_negative_exponents_without_changing_output_even_in_the_zero_ring) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type divisors[] = {{fq_value_type(7)}, {fq_value_type(3), fq_value_type(2)}};
    for (const polynomial_type &divisor : divisors) {
        math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);
        polynomial_type result = {fq_value_type(17)};
        BOOST_CHECK_THROW(math::powmod_x(result, -1, divisor_context, arithmetic_context), std::invalid_argument);
        BOOST_CHECK(result == polynomial_type({fq_value_type(17)}));
        BOOST_CHECK_THROW(
            math::powmod_x(result, boost::multiprecision::cpp_int(-1), divisor_context, arithmetic_context),
            std::invalid_argument);
        BOOST_CHECK(result == polynomial_type({fq_value_type(17)}));
    }
}

BOOST_AUTO_TEST_CASE(powmod_x_uses_the_existing_squaring_inverse_precision_checks) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type {}, options);
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::zero(),
                                     fq_value_type::one()};
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);

    polynomial_type result = {fq_value_type(17)};
    // Squaring X^2 produces X^4, whose reduction modulo this cubic needs two inverse coefficients.
    BOOST_CHECK_THROW(math::powmod_x(result, 4, divisor_context, arithmetic_context), std::invalid_argument);
    BOOST_CHECK(result == polynomial_type({fq_value_type(17)}));

    // X^3 only squares X and cancels the leading term of X*X^2; the stored constant inverse coefficient suffices.
    math::powmod_x(result, 3, divisor_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({-fq_value_type::one()}));
    math::powmod_x(result, 0, divisor_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::one()}));
}

BOOST_AUTO_TEST_CASE(powmod_x_supports_extension_coefficients_and_both_backends) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    const polynomial_type x = {fq12_value_type::zero(), fq12_value_type::one()};
    const polynomial_type divisor = {fq12_value(37), fq12_value(49), fq12_value(61), fq12_value(73)};
    polynomial_arithmetic::polynomial_context<schoolbook_backend> schoolbook_context;
    math::polynomial_divisor_context<schoolbook_backend> schoolbook_divisor(divisor, 2, schoolbook_context);

    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<mixed_radix_backend> mixed_radix_context(mixed_radix_backend(18),
                                                                                       options);
    math::polynomial_divisor_context<mixed_radix_backend> mixed_radix_divisor(divisor, 2, mixed_radix_context);

    for (const std::size_t exponent : {0, 1, 3, 8, 11, 21}) {
        const polynomial_type expected = repeated_powmod(x, exponent, schoolbook_divisor, schoolbook_context);
        polynomial_type reference, result;
        math::powmod(reference, x, exponent, mixed_radix_divisor, mixed_radix_context);
        BOOST_CHECK(reference == expected);
        math::powmod_x(result, exponent, schoolbook_divisor, schoolbook_context);
        BOOST_CHECK(result == expected);
        math::powmod_x(result, exponent, mixed_radix_divisor, mixed_radix_context);
        BOOST_CHECK(result == expected);
    }
}

BOOST_AUTO_TEST_CASE(powmod_x_replaces_backend_products_with_coefficient_steps) {
    using backend_type = counting_schoolbook_backend;
    using polynomial_type = typename backend_type::polynomial_type;

    std::size_t squares = 0, products = 0;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type(squares, products));
    const polynomial_type x = {fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type divisor = {fq_value_type(3), fq_value_type(2), fq_value_type::zero(), fq_value_type(4)};
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 2, arithmetic_context);

    // This small modulus uses basecase division: these counters isolate exponentiation's backend operations.
    // Fast reduction may itself multiply polynomials, which powmod_x deliberately leaves to squaremod.
    for (const unsigned exponent : {1, 2, 8, 11, 31}) {
        squares = products = 0;
        polynomial_type reference;
        math::powmod(reference, x, exponent, divisor_context, arithmetic_context);
        BOOST_CHECK_EQUAL(squares, std::bit_width(exponent) - 1);
        BOOST_CHECK_EQUAL(products, std::popcount(exponent) - 1);

        squares = products = 0;
        polynomial_type result;
        math::powmod_x(result, exponent, divisor_context, arithmetic_context);
        BOOST_CHECK(result == reference);
        BOOST_CHECK_EQUAL(squares, std::bit_width(exponent) - 1);
        BOOST_CHECK_EQUAL(products, 0);
    }
}

BOOST_AUTO_TEST_SUITE_END()
