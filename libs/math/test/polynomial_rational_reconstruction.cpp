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

#define BOOST_TEST_MODULE polynomial_rational_reconstruction_test

#include <array>
#include <limits>
#include <random>
#include <stdexcept>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/babybear.hpp>
#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/reconstruction/polynomial_rational_reconstruction.hpp>
#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>

namespace {
    namespace fields = nil::crypto3::algebra::fields;
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = math::polynomial_arithmetic;

    using value_type = fields::babybear::value_type;
    using backend_type = polynomial_arithmetic::schoolbook_backend<value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    struct counting_backend : backend_type {
        explicit counting_backend(std::size_t &matrix_products) : matrix_products(matrix_products) {
        }

        template<std::size_t Columns>
        bool try_multiply_matrix(const std::array<polynomial_type *, 2 * Columns> &,
                                 const std::array<const polynomial_type *, 4> &,
                                 const std::array<const polynomial_type *, 2 * Columns> &) {
            ++matrix_products;
            return false;
        }

        std::size_t &matrix_products;
    };

    using fq_value_type = fields::alt_bn128_base_field<254>::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;
    using fq12_backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using fq12_polynomial_type = typename fq12_backend_type::polynomial_type;

    fq12_value_type fq12_value(std::size_t first_coordinate) {
        fq12_value_type result = fq12_value_type::zero();
        for (std::size_t i = 0; i < fq12_field_type::arity; ++i) {
            result.coordinate(i) = fq_value_type(first_coordinate + i);
        }
        return result;
    }

    const polynomial_type modulus = {value_type(3), value_type(3), value_type(3), value_type(2), value_type::one()};
    const polynomial_type residue = {value_type(2), value_type::one(), value_type::one(), value_type::one()};
    const polynomial_type expected_numerator = {value_type::one()};
    const polynomial_type expected_denominator = {value_type(2), value_type(2), value_type::one()};
}    // namespace

BOOST_AUTO_TEST_SUITE(polynomial_rational_reconstruction_test_suite)

BOOST_AUTO_TEST_CASE(reconstruction_stops_at_the_degree_bound_and_normalizes_the_denominator) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    polynomial_type numerator;
    polynomial_type denominator;

    BOOST_REQUIRE(math::rational_reconstruct(numerator, denominator, residue, modulus, 0, 2, arithmetic_context));
    BOOST_CHECK(numerator == expected_numerator);
    BOOST_CHECK(denominator == expected_denominator);

    polynomial_type residue_times_denominator;
    arithmetic_context.multiply(residue_times_denominator, residue, denominator);
    polynomial_type divisor_quotient;
    polynomial_type congruence_remainder;
    math::division(divisor_quotient, congruence_remainder, residue_times_denominator, modulus);
    BOOST_CHECK(congruence_remainder == numerator);
}

BOOST_AUTO_TEST_CASE(reconstruction_failure_preserves_the_outputs) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    polynomial_type numerator = {value_type(7)};
    polynomial_type denominator = {value_type(11)};

    BOOST_CHECK(!math::rational_reconstruct(numerator, denominator, residue, modulus, 0, 1, arithmetic_context));
    BOOST_CHECK(numerator == polynomial_type({value_type(7)}));
    BOOST_CHECK(denominator == polynomial_type({value_type(11)}));
}

BOOST_AUTO_TEST_CASE(reconstruction_supports_bn254_fq12_coefficients) {
    polynomial_arithmetic::polynomial_context<fq12_backend_type> arithmetic_context;

    const fq12_polynomial_type first_quotient = {fq12_value(1), fq12_value(13)};
    const fq12_polynomial_type second_quotient = {fq12_value(25), fq12_value(37)};
    const fq12_polynomial_type second_remainder = {fq12_value(49), fq12_value(61), fq12_value(73)};
    const fq12_polynomial_type expected_numerator_before_normalization = {fq12_value(85)};

    // Construct two Euclidean steps backwards:
    //
    //     residue = second_quotient * second_remainder + expected_numerator,
    //     modulus = first_quotient * residue + second_remainder.
    fq12_polynomial_type product;
    fq12_polynomial_type extension_residue;
    arithmetic_context.multiply(product, second_quotient, second_remainder);
    math::addition(extension_residue, product, expected_numerator_before_normalization);
    fq12_polynomial_type extension_modulus;
    arithmetic_context.multiply(product, first_quotient, extension_residue);
    math::addition(extension_modulus, product, second_remainder);

    // The matching residue coefficient is 1 + second_quotient * first_quotient. Normalize it and the numerator by
    // the same Fq12 scalar to obtain the function's canonical result.
    fq12_polynomial_type extension_expected_denominator;
    arithmetic_context.multiply(extension_expected_denominator, second_quotient, first_quotient);
    math::addition(extension_expected_denominator, extension_expected_denominator,
                   fq12_polynomial_type {fq12_value_type::one()});
    const fq12_value_type normalization = extension_expected_denominator.back().inversed();
    fq12_polynomial_type extension_expected_numerator;
    math::scalar_multiplication(extension_expected_numerator, expected_numerator_before_normalization, normalization);
    math::scalar_multiplication(extension_expected_denominator, extension_expected_denominator, normalization);

    fq12_polynomial_type numerator;
    fq12_polynomial_type denominator;
    BOOST_REQUIRE(math::rational_reconstruct(numerator, denominator, extension_residue, extension_modulus, 0, 2,
                                             arithmetic_context));
    BOOST_CHECK(numerator == extension_expected_numerator);
    BOOST_CHECK(denominator == extension_expected_denominator);
}

BOOST_AUTO_TEST_CASE(reconstruction_handles_zero_and_input_aliasing) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    polynomial_type numerator;
    polynomial_type denominator;
    const polynomial_type zero = {value_type::zero()};

    BOOST_REQUIRE(math::rational_reconstruct(numerator, denominator, zero, modulus, 0, 0, arithmetic_context));
    BOOST_CHECK(numerator == zero);
    BOOST_CHECK(denominator == polynomial_type({value_type::one()}));

    polynomial_type residue_alias = residue;
    polynomial_type modulus_alias = modulus;
    BOOST_REQUIRE(math::rational_reconstruct(residue_alias, modulus_alias, residue_alias, modulus_alias, 0, 2,
                                             arithmetic_context));
    BOOST_CHECK(residue_alias == expected_numerator);
    BOOST_CHECK(modulus_alias == expected_denominator);
}

BOOST_AUTO_TEST_CASE(reconstruction_rejects_invalid_inputs_and_bounds) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    polynomial_type numerator;
    polynomial_type denominator;

    BOOST_CHECK_THROW(math::rational_reconstruct(numerator, numerator, residue, modulus, 0, 2, arithmetic_context),
                      std::invalid_argument);

    polynomial_type empty;
    empty.get_storage().clear();
    BOOST_CHECK_THROW(math::rational_reconstruct(numerator, denominator, empty, modulus, 0, 2, arithmetic_context),
                      std::invalid_argument);

    polynomial_type noncanonical(2);
    noncanonical[0] = value_type::one();
    noncanonical[1] = value_type::zero();
    BOOST_CHECK_THROW(
        math::rational_reconstruct(numerator, denominator, noncanonical, modulus, 0, 2, arithmetic_context),
        std::invalid_argument);

    const polynomial_type constant_modulus = {value_type::one()};
    BOOST_CHECK_THROW(
        math::rational_reconstruct(numerator, denominator, residue, constant_modulus, 0, 0, arithmetic_context),
        std::invalid_argument);

    const polynomial_type quadratic_modulus = {value_type::one(), value_type::zero(), value_type::one()};
    const polynomial_type unreduced_residue = {value_type::one(), value_type::one(), value_type::one()};
    BOOST_CHECK_THROW(math::rational_reconstruct(numerator, denominator, unreduced_residue, quadratic_modulus, 0, 0,
                                                 arithmetic_context),
                      std::invalid_argument);

    BOOST_CHECK_THROW(math::rational_reconstruct(numerator, denominator, residue, modulus, 2, 2, arithmetic_context),
                      std::invalid_argument);

    const std::size_t maximum_size = std::numeric_limits<std::size_t>::max();
    BOOST_CHECK_THROW(
        math::rational_reconstruct(numerator, denominator, residue, modulus, maximum_size, 0, arithmetic_context),
        std::invalid_argument);
    BOOST_CHECK_THROW(
        math::rational_reconstruct(numerator, denominator, residue, modulus, 0, maximum_size, arithmetic_context),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(balanced_reconstruction_matches_individual_steps_across_the_half_gcd_cutoff) {
    std::size_t matrix_products = 0;
    polynomial_arithmetic::polynomial_context<counting_backend> context {counting_backend(matrix_products)};
    // A 128-coefficient FFT fits every classical denominator update on these inputs, including degree 255.
    // The larger matrix products would exceed this plan.
    using fft_backend = polynomial_arithmetic::mixed_radix_backend<fields::babybear>;
    polynomial_arithmetic::polynomial_context<fft_backend> classical {fft_backend(128)};
    std::mt19937 random(1824);

    for (const std::size_t degree : {127, 128, 129, 255}) {
        BOOST_TEST_CONTEXT("modulus degree " << degree) {
            polynomial_type f(degree + 1);
            polynomial_type r(degree);
            for (auto *input : {&f, &r}) {
                for (auto &coefficient : *input) {
                    coefficient = value_type(random());
                }
                input->back() = value_type::one();
            }
            polynomial_type p, q, expected_p, expected_q;
            const auto numerator_bound = degree / 2;
            const auto denominator_bound = (degree - 1) / 2;
            BOOST_REQUIRE(math::rational_reconstruct(expected_p, expected_q, r, f, numerator_bound, denominator_bound,
                                                     classical));
            matrix_products = 0;
            BOOST_REQUIRE(math::rational_reconstruct(p, q, r, f, numerator_bound, denominator_bound, context));
            BOOST_CHECK(p == expected_p);
            BOOST_CHECK(q == expected_q);
            BOOST_CHECK_EQUAL(matrix_products != 0, degree >= 128);

            matrix_products = 0;
            BOOST_CHECK(!math::rational_reconstruct(p, q, r, f, numerator_bound, 0, context));
            BOOST_CHECK(p == expected_p);
            BOOST_CHECK(q == expected_q);
            BOOST_CHECK_EQUAL(matrix_products, 0U);
        }
    }
}

BOOST_AUTO_TEST_CASE(half_gcd_selects_the_first_remainder_within_the_inclusive_bound) {
    polynomial_arithmetic::polynomial_context<backend_type> context;
    const polynomial_type x(value_type::one(), 1);
    for (const std::size_t degree : {128, 129}) {
        const std::size_t bound = degree / 2;
        const polynomial_type r(value_type::one(), degree - 1);
        polynomial_type f(value_type(2), degree);
        polynomial_type p, q;

        // With f = 2 X^degree + 3 X^e, multiplying r by X gives -(3/2) X^e modulo f.
        for (const auto exponent : {bound, bound - 1}) {
            f[exponent] = value_type(3);
            p = f;
            q = r;
            BOOST_REQUIRE(math::rational_reconstruct(p, q, q, p, bound, (degree - 1) / 2, context));
            BOOST_CHECK(p == polynomial_type(-value_type(3) * value_type(2).inversed(), exponent));
            BOOST_CHECK(q == x);
            f[exponent] = value_type::zero();
        }

        BOOST_REQUIRE(math::rational_reconstruct(p, q, r, f, bound, (degree - 1) / 2, context));
        BOOST_CHECK(p == polynomial_type({value_type::zero()}));
        BOOST_CHECK(q == x);

        q = polynomial_type(value_type(3), bound);
        BOOST_REQUIRE(math::rational_reconstruct(p, q, q, f, bound, (degree - 1) / 2, context));
        BOOST_CHECK(p == polynomial_type(value_type(3), bound));
        BOOST_CHECK(q == polynomial_type({value_type::one()}));
    }
}

BOOST_AUTO_TEST_SUITE_END()
