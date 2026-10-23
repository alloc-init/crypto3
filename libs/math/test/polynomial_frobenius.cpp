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

#define BOOST_TEST_MODULE polynomial_frobenius_test

#include <cstddef>
#include <limits>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/algebra/fields/field_order.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_frobenius.hpp>
#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/bn254_small_prime_fft_backend.hpp>
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

    struct counting_product_backend : polynomial_arithmetic::schoolbook_backend<fields::babybear::value_type> {
        using base_type = polynomial_arithmetic::schoolbook_backend<fields::babybear::value_type>;
        explicit counting_product_backend(std::size_t &operations) : operations(operations) {
        }
        void square(polynomial_type &output, const polynomial_type &input) const {
            ++operations;
            base_type::square(output, input);
        }
        void multiply(polynomial_type &output, const polynomial_type &left, const polynomial_type &right) const {
            ++operations;
            base_type::multiply(output, left, right);
        }
        std::size_t &operations;
    };
}    // namespace

BOOST_AUTO_TEST_SUITE(polynomial_frobenius_test_suite)

BOOST_AUTO_TEST_CASE(precomputation_matches_direct_exponentiation) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type x = {fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    polynomial_type expected;
    math::powmod(expected, x, fields::field_order<fq_field_type>(), frobenius_context.divisor_context(),
                 arithmetic_context);
    BOOST_CHECK(frobenius_context.x_to_field_order() == expected);
}

BOOST_AUTO_TEST_CASE(frobenius_matches_direct_exponentiation_and_may_alias_the_input) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type input = {fq_value_type(5), fq_value_type(7), fq_value_type(11)};
    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    polynomial_type expected;
    math::powmod(expected, input, fields::field_order<fq_field_type>(), frobenius_context.divisor_context(),
                 arithmetic_context);

    polynomial_type result;
    math::frobenius_map(result, input, frobenius_context, arithmetic_context);
    BOOST_CHECK(result == expected);

    polynomial_type aliased = input;
    math::frobenius_map(aliased, aliased, frobenius_context, arithmetic_context);
    BOOST_CHECK(aliased == expected);

    const polynomial_type unreduced = {fq_value_type(13), fq_value_type(17), fq_value_type(19), fq_value_type(23),
                                       fq_value_type(29)};
    math::powmod(expected, unreduced, fields::field_order<fq_field_type>(), frobenius_context.divisor_context(),
                 arithmetic_context);

    math::frobenius_map(result, unreduced, frobenius_context, arithmetic_context);
    BOOST_CHECK(result == expected);

    aliased = unreduced;
    math::frobenius_map(aliased, aliased, frobenius_context, arithmetic_context);
    BOOST_CHECK(aliased == expected);
}

BOOST_AUTO_TEST_CASE(iterated_frobenius_matches_a_power_of_the_field_order) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type input = {fq_value_type(5), fq_value_type(7), fq_value_type(11)};
    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    const boost::multiprecision::cpp_int order = fields::field_order<fq_field_type>();
    const boost::multiprecision::cpp_int order_squared = order * order;
    polynomial_type expected;
    math::powmod(expected, input, order_squared, frobenius_context.divisor_context(), arithmetic_context);

    polynomial_type result;
    math::frobenius_map(result, input, 2, frobenius_context, arithmetic_context);
    BOOST_CHECK(result == expected);

    polynomial_type aliased = input;
    math::frobenius_map(aliased, aliased, 2, frobenius_context, arithmetic_context);
    BOOST_CHECK(aliased == expected);

    const polynomial_type unreduced = {fq_value_type(5), fq_value_type(7), fq_value_type(11), fq_value_type(13),
                                       fq_value_type(17)};
    polynomial_type reduced;
    math::remainder(reduced, unreduced, frobenius_context.divisor_context(), arithmetic_context);

    polynomial_type zero_iteration_alias = unreduced;
    math::frobenius_map(zero_iteration_alias, zero_iteration_alias, 0, frobenius_context, arithmetic_context);
    BOOST_CHECK(zero_iteration_alias == reduced);
}

BOOST_AUTO_TEST_CASE(frobenius_uses_the_extension_coefficient_field_order) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq12_value_type::one(), fq12_value_type::zero(), fq12_value_type::one()};
    const polynomial_type input = {fq12_value(5), fq12_value(17)};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    // X^2 = -1 and Q = q^12 is 1 modulo 4, so X^Q = X. The Q-power map also fixes every Fq12 coefficient.
    const polynomial_type x = {fq12_value_type::zero(), fq12_value_type::one()};
    BOOST_CHECK(frobenius_context.x_to_field_order() == x);

    polynomial_type result;
    math::frobenius_map(result, input, frobenius_context, arithmetic_context);
    BOOST_CHECK(result == input);
}

BOOST_AUTO_TEST_CASE(frobenius_modulo_a_nonzero_constant_is_zero) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq_value_type(7)};
    const polynomial_type input = {fq_value_type(2), fq_value_type(3)};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    polynomial_type result;
    math::frobenius_map(result, input, frobenius_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));
}

BOOST_AUTO_TEST_CASE(mixed_radix_frobenius_reuses_cyclic_reduction_and_cached_composition) {
    using mixed_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using reference_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename mixed_backend::polynomial_type;
    polynomial_type divisor(7, fq12_value_type::one()), input(6, fq12_value_type::one());
    for (std::size_t i = 0; i < divisor.size(); ++i) {
        divisor[i] = fq12_value(12 * i + 1);
    }
    for (std::size_t i = 0; i < input.size(); ++i) {
        input[i] = fq12_value(12 * i + 85);
    }
    // Force inverse-based division at this bounded size. For degree-six H and a dense reduced product,
    // the quotient uses order nine, while reconstruction uses order six instead of the low product's 18.
    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<mixed_backend> mixed(mixed_backend(18), options);
    polynomial_arithmetic::polynomial_context<reference_backend> reference;
    const math::polynomial_divisor_context<reference_backend> reference_divisor(divisor, 5, reference);
    const math::polynomial_frobenius_context<mixed_backend> frobenius(divisor, mixed);
    const auto order = fields::field_order<fq12_field_type>();
    const polynomial_type x {fq12_value_type::zero(), fq12_value_type::one()};
    polynomial_type expected, result;
    math::powmod(expected, x, order, reference_divisor, reference);
    BOOST_CHECK(frobenius.x_to_field_order() == expected);
    math::powmod(expected, input, order, reference_divisor, reference);
    math::frobenius_map(result, input, frobenius, mixed);
    BOOST_CHECK(result == expected);
    auto alias = input;
    math::frobenius_map(alias, alias, frobenius, mixed);
    BOOST_CHECK(alias == expected);
    // Reuse the same composition cache with another input. Frobenius fixes the added constant one.
    input[0] += fq12_value_type::one();
    expected[0] += fq12_value_type::one();
    math::frobenius_map(result, input, frobenius, mixed);
    BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_CASE(frobenius_context_rejects_invalid_precomputation_inputs) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type zero_divisor = {fq_value_type::zero()};
    BOOST_CHECK_THROW(math::polynomial_frobenius_context<backend_type>(zero_divisor, arithmetic_context),
                      std::invalid_argument);

    polynomial_arithmetic::polynomial_context_options invalid_options;
    invalid_options.modular_composition_cached_power_limit = 0;
    polynomial_arithmetic::polynomial_context<backend_type> invalid_arithmetic_context(backend_type {},
                                                                                       invalid_options);
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    BOOST_CHECK_THROW(math::polynomial_frobenius_context<backend_type>(divisor, invalid_arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(frobenius_products_match_geometric_exponents_and_support_aliasing) {
    using value = fields::babybear::value_type;
    using backend = polynomial_arithmetic::schoolbook_backend<value>;
    using polynomial = backend::polynomial_type;
    polynomial_arithmetic::polynomial_context<backend> context;
    const polynomial modulus {value(5), value(7), value(3), value(2)};
    const math::polynomial_divisor_context<backend> divisor(modulus, 4, context);
    const polynomial input {value(2), value(5), value(7), value(11), value(13), value(17)};
    const boost::multiprecision::cpp_int q = fields::field_order<fields::babybear>();
    for (const boost::multiprecision::cpp_int &step :
         {boost::multiprecision::cpp_int(q), boost::multiprecision::cpp_int(q * q)}) {
        polynomial image;
        math::powmod_x(image, step, divisor, context);
        for (const std::size_t count : {0, 1, 2, 3, 4, 7, 8, 9, 15, 16, 17, 31, 32, 33}) {
            boost::multiprecision::cpp_int exponent = 0, power = 1;
            for (std::size_t i = 0; i < count; ++i) {
                exponent += power;
                power *= step;
            }
            polynomial expected, actual;
            math::powmod(expected, input, exponent, divisor, context);
            math::frobenius_product(actual, input, count, image, divisor, context);
            BOOST_CHECK(actual == expected);
            auto alias = input;
            math::frobenius_product(alias, alias, count, image, divisor, context);
            BOOST_CHECK(alias == expected);
            alias = image;
            math::frobenius_product(alias, input, count, alias, divisor, context);
            BOOST_CHECK(alias == expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(frobenius_products_use_bounded_work_for_maximum_term_count) {
    using backend = counting_product_backend;
    using polynomial = backend::polynomial_type;
    using value = polynomial::value_type;
    std::size_t operations = 0;
    polynomial_arithmetic::polynomial_context<backend> context {backend(operations)};
    const polynomial modulus {-value::one(), value::zero(), value::one()};
    const math::polynomial_divisor_context<backend> divisor(modulus, 1, context);
    const polynomial x {value::zero(), value::one()};
    operations = 0;
    polynomial actual;
    // Q is odd and X^2 = 1, so F fixes X and the product is X for every odd count.
    math::frobenius_product(actual, x, std::numeric_limits<std::size_t>::max(), x, divisor, context);
    BOOST_CHECK(actual == x);
    BOOST_CHECK_LT(operations, 8 * std::numeric_limits<std::size_t>::digits);
}

BOOST_AUTO_TEST_CASE(frobenius_products_over_fp12_reuse_both_q_and_q_squared_images) {
    using backend = polynomial_arithmetic::bn254_small_prime_fft_backend;
    using reference_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial = backend::polynomial_type;
    using value = fq12_value_type;
    value nonsquare = fq12_value(1);
    while (nonsquare.is_square()) {
        nonsquare += value::one();
    }
    const value scale = fq12_value(31);
    const polynomial modulus {-scale * nonsquare, value::zero(), scale};
    const polynomial input {fq12_value(5), fq12_value(17)};
    const polynomial x {value::zero(), value::one()}, minus_x {value::zero(), -value::one()};
    polynomial_arithmetic::polynomial_context<backend> context {backend(8)};
    polynomial_arithmetic::polynomial_context<reference_backend> reference;
    const math::polynomial_divisor_context<backend> divisor(modulus, 1, context);
    const math::polynomial_divisor_context<reference_backend> reference_divisor(modulus, 1, reference);
    const value norm = input[0].squared() - nonsquare * input[1].squared();
    // X^2 = nonsquare: Q-power conjugation negates X but fixes all twelve coefficient coordinates.
    for (const std::size_t count : {0, 1, 2, 3, 7, 8, 9}) {
        polynomial expected, actual;
        math::scalar_multiplication(expected, count % 2 == 0 ? polynomial {value::one()} : input, norm.pow(count / 2));
        math::frobenius_product(actual, input, count, minus_x, divisor, context);
        BOOST_CHECK(actual == expected);
        // The Q^2-power map fixes every element in this quadratic extension.
        math::powmod(expected, input, count, reference_divisor, reference);
        math::frobenius_product(actual, input, count, x, divisor, context);
        BOOST_CHECK(actual == expected);
    }
}

BOOST_AUTO_TEST_CASE(frobenius_products_handle_constants_zero_rings_and_invalid_inputs) {
    using value = fields::babybear::value_type;
    using backend = polynomial_arithmetic::schoolbook_backend<value>;
    using polynomial = backend::polynomial_type;
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_cached_power_limit = 0;
    polynomial_arithmetic::polynomial_context<backend> context(backend {}, options);
    const polynomial modulus {value(5), value(7), value(3), value(2)};
    const math::polynomial_divisor_context<backend> divisor(modulus, 4, context);
    const polynomial x {value::zero(), value::one()};
    polynomial image;
    math::powmod_x(image, fields::field_order<fields::babybear>(), divisor, context);
    for (const auto &constant : {value::zero(), value::one(), value(7)}) {
        for (const std::size_t count : {0, 1, 2, 17}) {
            polynomial result;
            math::frobenius_product(result, polynomial {constant}, count, image, divisor, context);
            BOOST_CHECK(result == polynomial({constant.pow(count)}));
        }
    }
    polynomial output {value(19)};
    math::frobenius_product(output, x, 1, image, divisor, context);
    BOOST_CHECK(output == x);
    BOOST_CHECK_THROW(math::frobenius_product(output, x, 2, image, divisor, context), std::invalid_argument);
    BOOST_CHECK(output == x);
    polynomial empty;
    empty.get_storage().clear();
    polynomial padded(2);
    padded[0] = value::one();
    padded[1] = value::zero();
    BOOST_CHECK_THROW(math::frobenius_product(output, empty, 1, image, divisor, context), std::invalid_argument);
    BOOST_CHECK_THROW(math::frobenius_product(output, padded, 1, image, divisor, context), std::invalid_argument);
    for (const auto &bad_image : {empty, padded, modulus}) {
        BOOST_CHECK_THROW(math::frobenius_product(output, x, 1, bad_image, divisor, context), std::invalid_argument);
    }
    const math::polynomial_divisor_context<backend> constant_divisor(polynomial {value(7)}, 1, context);
    for (const std::size_t count : {0, 1, 9}) {
        math::frobenius_product(output, x, count, polynomial {value::zero()}, constant_divisor, context);
        BOOST_CHECK(output == polynomial({value::zero()}));
    }
}

BOOST_AUTO_TEST_CASE(characteristic_product_twists_coefficients_and_matches_direct_powers) {
    using backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial = backend::polynomial_type;
    const polynomial modulus {fq12_value(3), fq12_value(19), fq12_value(37), fq12_value(53)};
    const polynomial input {fq12_value(71), fq12_value(89), fq12_value(107)};
    BOOST_REQUIRE(input[0].Frobenius_map(1) != input[0]);
    const boost::multiprecision::cpp_int p = fields::field_characteristic<fq12_field_type>();
    for (const std::size_t cache_limit : {1, 3}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = cache_limit;
        polynomial_arithmetic::polynomial_context<backend> context(backend {}, options);
        const math::polynomial_divisor_context<backend> divisor(modulus, 2, context);
        polynomial image;
        math::powmod_x(image, p, divisor, context);
        for (const std::size_t count : {0, 1, 2, 3, 6, 12, 13}) {
            boost::multiprecision::cpp_int exponent = 0;
            for (std::size_t i = 0; i < count; ++i) {
                exponent = exponent * p + 1;
            }
            // Constants also need the coefficient twist: they are not fixed by absolute Frobenius.
            for (const auto &base : {input, polynomial {input[0]}, polynomial {fq12_value_type::zero()}}) {
                polynomial expected, result;
                math::powmod(expected, base, exponent, divisor, context);
                math::detail::frobenius_product_impl<true>(result, base, count, image, divisor, context);
                BOOST_CHECK(result == expected);
                if (count == 12) {
                    result = base;
                    math::detail::frobenius_product_impl<true>(result, result, count, image, divisor, context);
                    BOOST_CHECK(result == expected);
                    result = image;
                    math::detail::frobenius_product_impl<true>(result, base, count, result, divisor, context);
                    BOOST_CHECK(result == expected);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
