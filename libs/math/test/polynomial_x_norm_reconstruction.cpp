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

#define BOOST_TEST_MODULE polynomial_x_norm_reconstruction_test

#include <cstddef>
#include <stdexcept>

#include <boost/random/mersenne_twister.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/algebra/random_element.hpp>

#include <nil/crypto3/math/polynomial/backends/bn254_small_prime_fft_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/reconstruction/polynomial_x_norm_reconstruction.hpp>

namespace {
    namespace fields = nil::crypto3::algebra::fields;
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = math::polynomial_arithmetic;

    using field_type = fields::babybear;
    using value_type = field_type::value_type;
    using backend_type = polynomial_arithmetic::schoolbook_backend<value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq_value_type = fq_field_type::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;
    using fq12_schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using fq12_mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;

    value_type first_quadratic_non_residue() {
        value_type candidate(2);
        while (candidate.is_square()) {
            candidate = candidate + value_type::one();
        }
        return candidate;
    }

    value_type next_quadratic_non_residue(value_type candidate) {
        do {
            candidate = candidate + value_type::one();
        } while (candidate.is_square());
        return candidate;
    }

    template<typename Backend>
    typename Backend::polynomial_type
        context_multiply(const typename Backend::polynomial_type &left, const typename Backend::polynomial_type &right,
                         polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        typename Backend::polynomial_type product;
        arithmetic_context.multiply(product, left, right);
        return product;
    }

    polynomial_type representable_irreducible_quadratic() {
        return {value_type::one(), value_type::zero() - value_type(79), value_type::one()};
    }

    fq12_value_type fq12_scalar(std::size_t value) {
        fq12_value_type result = fq12_value_type::zero();
        result.coordinate(0) = fq_value_type(value);
        return result;
    }

    template<typename Backend>
    void check_bn254_fq12_linear_recovery(polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                          std::size_t seed) {
        using extension_polynomial_type = typename Backend::polynomial_type;

        const extension_polynomial_type g = {fq12_scalar(9), -fq12_scalar(4)};
        boost::random::mt19937 rng(seed);
        auto coefficient_generator = [&] { return nil::crypto3::algebra::random_element<fq12_field_type>(rng); };

        const auto result = math::recover_irreducible_polynomial_x_norm_representation<Backend>(g, arithmetic_context,
                                                                                                coefficient_generator);
        BOOST_REQUIRE(result.has_value());
        BOOST_CHECK_EQUAL(result->p.size(), 1);
        BOOST_CHECK_EQUAL(result->q.size(), 1);
        BOOST_CHECK(result->p[0] * result->p[0] == fq12_scalar(9));
        BOOST_CHECK(result->q[0] * result->q[0] == fq12_scalar(4));
        BOOST_CHECK(math::evaluate_polynomial_x_norm<Backend>(*result, arithmetic_context) == g);
    }

    template<typename Backend>
    void check_bn254_fq12_high_level_recovery(polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                              std::size_t seed) {
        using extension_polynomial_type = typename Backend::polynomial_type;

        const extension_polynomial_type odd_factor = {-fq12_scalar(9), fq12_value_type::one()};
        const extension_polynomial_type even_multiplicity_factor = {-fq12_scalar(4), fq12_value_type::one()};
        extension_polynomial_type even_factor_squared;
        arithmetic_context.square(even_factor_squared, even_multiplicity_factor);
        extension_polynomial_type h = context_multiply<Backend>(odd_factor, even_factor_squared, arithmetic_context);
        math::scalar_multiplication(h, h, fq12_scalar(25));

        boost::random::mt19937 rng(seed);
        auto coefficient_generator = [&] { return nil::crypto3::algebra::random_element<fq12_field_type>(rng); };
        const auto result =
            math::recover_polynomial_x_norm_representation<Backend>(h, arithmetic_context, coefficient_generator);

        BOOST_REQUIRE(result.has_value());
        BOOST_CHECK(math::evaluate_polynomial_x_norm<Backend>(*result, arithmetic_context) == h);
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(polynomial_x_norm_reconstruction_test_suite)

BOOST_AUTO_TEST_CASE(multiplies_x_norm_representations_and_preserves_the_exact_norm_product) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const math::polynomial_x_norm_representation<polynomial_type> left = {
        polynomial_type {value_type::one(), value_type(2)}, polynomial_type {value_type(3)}};
    const math::polynomial_x_norm_representation<polynomial_type> right = {
        polynomial_type {value_type(4), value_type(5)}, polynomial_type {value_type(6), value_type(7)}};

    const auto product =
        math::multiply_polynomial_x_norm_representations<backend_type>(left, right, arithmetic_context);
    BOOST_CHECK(product.p == polynomial_type({value_type(4), value_type(31), value_type(31)}));
    BOOST_CHECK(product.q == polynomial_type({value_type(18), value_type(34), value_type(14)}));

    const polynomial_type left_norm = math::evaluate_polynomial_x_norm<backend_type>(left, arithmetic_context);
    const polynomial_type right_norm = math::evaluate_polynomial_x_norm<backend_type>(right, arithmetic_context);
    polynomial_type expected_norm_product;
    arithmetic_context.multiply(expected_norm_product, left_norm, right_norm);
    BOOST_CHECK(math::evaluate_polynomial_x_norm<backend_type>(product, arithmetic_context) == expected_norm_product);

    const math::polynomial_x_norm_representation<polynomial_type> zero = {polynomial_type {value_type::zero()},
                                                                          polynomial_type {value_type::zero()}};
    const auto zero_product =
        math::multiply_polynomial_x_norm_representations<backend_type>(zero, right, arithmetic_context);
    BOOST_CHECK(zero_product.p == polynomial_type({value_type::zero()}));
    BOOST_CHECK(zero_product.q == polynomial_type({value_type::zero()}));
}

BOOST_AUTO_TEST_CASE(x_norm_representation_multiplication_rejects_malformed_inputs) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const math::polynomial_x_norm_representation<polynomial_type> valid = {polynomial_type {value_type::one()},
                                                                           polynomial_type {value_type::zero()}};
    math::polynomial_x_norm_representation<polynomial_type> malformed = valid;
    malformed.q.get_storage().clear();

    BOOST_CHECK_THROW(
        math::multiply_polynomial_x_norm_representations<backend_type>(valid, malformed, arithmetic_context),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(recovers_and_normalizes_an_irreducible_quadratic_with_the_stated_degree_bounds) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const math::polynomial_x_norm_representation<polynomial_type> original = {
        polynomial_type {value_type::one(), value_type::one()}, polynomial_type {value_type(9)}};
    const polynomial_type g = math::evaluate_polynomial_x_norm<backend_type>(original, arithmetic_context);

    // The discriminant of g is 81 * 77. Since 81 is square and 77 is nonsquare, g is irreducible.
    BOOST_REQUIRE(!value_type(77).is_square());
    BOOST_REQUIRE(g == polynomial_type({value_type::one(), value_type::zero() - value_type(79), value_type::one()}));

    boost::random::mt19937 rng(0x584E4F52);
    std::size_t generator_calls = 0;
    auto coefficient_generator = [&] {
        ++generator_calls;
        return nil::crypto3::algebra::random_element<field_type>(rng);
    };
    const auto result = math::recover_irreducible_polynomial_x_norm_representation<backend_type>(g, arithmetic_context,
                                                                                                 coefficient_generator);

    BOOST_REQUIRE(result.has_value());
    BOOST_CHECK_EQUAL(generator_calls, 0);
    BOOST_CHECK_LE(result->p.size() - 1, 1);
    BOOST_CHECK_LE(result->q.size() - 1, 0);
    // Rational reconstruction makes Q monic, initially producing lambda = 1/81. Exact normalization restores Q^2 = 81.
    BOOST_CHECK(result->q[0] != value_type::one());
    BOOST_CHECK(result->q[0] * result->q[0] == value_type(81));
    BOOST_CHECK(math::evaluate_polynomial_x_norm<backend_type>(*result, arithmetic_context) == g);
}

BOOST_AUTO_TEST_CASE(returns_no_value_when_x_is_nonsquare_modulo_the_irreducible_polynomial) {
    const value_type non_residue = first_quadratic_non_residue();
    BOOST_REQUIRE((value_type::zero() - value_type::one()).is_square());
    const polynomial_type g = {value_type::zero() - non_residue, value_type::zero(), value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::size_t generator_calls = 0;
    auto coefficient_generator = [&] {
        ++generator_calls;
        return value_type::zero();
    };

    const auto result = math::recover_irreducible_polynomial_x_norm_representation<backend_type>(g, arithmetic_context,
                                                                                                 coefficient_generator);
    BOOST_CHECK(!result.has_value());
    BOOST_CHECK_EQUAL(generator_calls, 0);
}

BOOST_AUTO_TEST_CASE(returns_no_value_when_the_scalar_multiple_cannot_be_normalized) {
    const value_type non_residue = first_quadratic_non_residue();
    const value_type minus_one = value_type::zero() - value_type::one();
    BOOST_REQUIRE(minus_one.is_square());

    // X is one modulo non_residue * (X - 1), but the reconstructed norm is multiplied by
    // lambda = -non_residue^-1. This lambda is nonsquare and cannot be removed by scaling P and Q.
    const polynomial_type g = {value_type::zero() - non_residue, non_residue};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    auto coefficient_generator = [&] { return non_residue; };

    const auto result = math::recover_irreducible_polynomial_x_norm_representation<backend_type>(g, arithmetic_context,
                                                                                                 coefficient_generator);
    BOOST_CHECK(!result.has_value());
}

BOOST_AUTO_TEST_CASE(odd_degree_norm_recovery_verifies_exact_identities_without_sampling) {
    polynomial_arithmetic::polynomial_context<backend_type> context;
    value_type c(2);
    const auto noncube_exponent = (field_type::modulus - 1u) / 3u;
    while (c.pow(noncube_exponent) == value_type::one()) {
        c += value_type::one();
    }
    c = c.squared();
    std::size_t generator_calls = 0;
    auto generator = [&]() -> value_type {
        ++generator_calls;
        throw std::logic_error("odd-degree recovery must not sample a polynomial nonresidue");
    };
    for (const std::size_t degree : {3, 9, 27}) {
        // X^degree - c is irreducible by the binomial criterion; translation preserves irreducibility.
        value_type shift = value_type::one();
        while (!(c - shift.pow(degree)).is_square()) {
            shift += value_type::one();
        }
        polynomial_type monic {value_type::one()};
        const polynomial_type linear {shift, value_type::one()};
        for (std::size_t i = 0; i < degree; ++i) {
            context.multiply(monic, monic, linear);
        }
        monic[0] -= c;
        for (const auto &scale : {value_type::one(), value_type(49), first_quadratic_non_residue()}) {
            polynomial_type g;
            math::scalar_multiplication(g, monic, scale);
            const auto result = math::recover_irreducible_polynomial_x_norm_representation(g, context, generator);
            // -1 is square in BabyBear. The leading coefficient of an odd-degree norm is minus a square.
            BOOST_REQUIRE_EQUAL(result.has_value(), scale.is_square());
            if (result) {
                BOOST_CHECK(math::evaluate_polynomial_x_norm(*result, context) == g);
                BOOST_CHECK_LE(result->p.size() - 1, degree / 2);
                BOOST_CHECK_LE(result->q.size() - 1, (degree - 1) / 2);
            }
        }
    }
    const polynomial_type x {value_type::zero(), value_type::one()};
    const auto zero_root = math::recover_irreducible_polynomial_x_norm_representation(x, context, generator);
    BOOST_REQUIRE(zero_root.has_value());
    BOOST_CHECK(math::evaluate_polynomial_x_norm(*zero_root, context) == x);
    const polynomial_type nonsquare_linear {-first_quadratic_non_residue(), value_type::one()};
    BOOST_CHECK(!math::recover_irreducible_polynomial_x_norm_representation(nonsquare_linear, context, generator));
    BOOST_CHECK_EQUAL(generator_calls, 0);
}

BOOST_AUTO_TEST_CASE(odd_degree_recovery_preserves_tonelli_shanks_when_composition_is_disabled) {
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_cached_power_limit = 0;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
    value_type c(2);
    const auto noncube_exponent = (field_type::modulus - 1u) / 3u;
    while (c.pow(noncube_exponent) == value_type::one()) {
        c += value_type::one();
    }
    const polynomial_type g {-c.squared(), value_type::zero(), value_type::zero(), value_type::one()};
    const auto nonresidue = first_quadratic_non_residue();
    std::size_t generator_calls = 0;
    auto generator = [&] {
        // A coefficient-field nonsquare remains nonsquare in an odd-degree extension.
        return generator_calls++ % 3 == 0 ? nonresidue : value_type::zero();
    };
    const auto result = math::recover_irreducible_polynomial_x_norm_representation(g, context, generator);
    BOOST_REQUIRE(result.has_value());
    BOOST_CHECK(math::evaluate_polynomial_x_norm(*result, context) == g);
    BOOST_CHECK_EQUAL(generator_calls, 3);

    const polynomial_type linear {value_type(9), -value_type(4)};
    const auto linear_result = math::recover_irreducible_polynomial_x_norm_representation(linear, context, generator);
    BOOST_REQUIRE(linear_result.has_value());
    BOOST_CHECK(math::evaluate_polynomial_x_norm(*linear_result, context) == linear);
    BOOST_CHECK_EQUAL(generator_calls, 3);
}

BOOST_AUTO_TEST_CASE(odd_degree_bn254_fp12_norm_recovery_agrees_across_backends_without_sampling) {
    using backend = polynomial_arithmetic::bn254_small_prime_fft_backend;
    using polynomial = backend::polynomial_type;
    boost::random::mt19937 rng(0xF0120DD);
    const boost::multiprecision::cpp_int noncube_exponent = (fields::field_order<fq12_field_type>() - 1) / 3;
    fq12_value_type c = fq12_value_type::one();
    for (std::size_t i = 0; i < 64 && c.pow(noncube_exponent) == fq12_value_type::one(); ++i) {
        c = nil::crypto3::algebra::random_element<fq12_field_type>(rng);
    }
    BOOST_REQUIRE(!c.is_zero());
    BOOST_REQUIRE(c.pow(noncube_exponent) != fq12_value_type::one());
    c = c.squared();
    fq12_value_type shift = fq12_value_type::one();
    for (std::size_t i = 0; i < 64 && !(c - shift.pow(3)).is_square(); ++i) {
        shift = nil::crypto3::algebra::random_element<fq12_field_type>(rng);
    }
    BOOST_REQUIRE((c - shift.pow(3)).is_square());
    polynomial g {shift.pow(3) - c, fq12_scalar(3) * shift.squared(), fq12_scalar(3) * shift, fq12_value_type::one()};
    const auto scale = nil::crypto3::algebra::random_element<fq12_field_type>(rng).squared();
    BOOST_REQUIRE(!scale.is_zero());
    math::scalar_multiplication(g, g, scale);
    auto generator = []() -> fq12_value_type {
        throw std::logic_error("odd-degree recovery must not sample a polynomial nonresidue");
    };
    polynomial_arithmetic::polynomial_context<backend> context {backend(16)};
    polynomial_arithmetic::polynomial_context<fq12_schoolbook_backend> reference;
    const auto result = math::recover_irreducible_polynomial_x_norm_representation(g, context, generator);
    const auto expected = math::recover_irreducible_polynomial_x_norm_representation(g, reference, generator);
    BOOST_REQUIRE(result.has_value());
    BOOST_REQUIRE(expected.has_value());
    BOOST_CHECK(math::evaluate_polynomial_x_norm(*result, reference) == g);
    BOOST_CHECK(result->p == expected->p);
    BOOST_CHECK(result->q == expected->q);
}

BOOST_AUTO_TEST_CASE(recovers_a_linear_bn254_fq12_norm_with_the_schoolbook_backend) {
    polynomial_arithmetic::polynomial_context<fq12_schoolbook_backend> arithmetic_context;
    check_bn254_fq12_linear_recovery(arithmetic_context, 0xF0125001);
}

BOOST_AUTO_TEST_CASE(recovers_a_linear_bn254_fq12_norm_with_the_mixed_radix_backend) {
    polynomial_arithmetic::polynomial_context<fq12_mixed_radix_backend> arithmetic_context {
        fq12_mixed_radix_backend(3)};
    check_bn254_fq12_linear_recovery(arithmetic_context, 0xF0125002);
}

BOOST_AUTO_TEST_CASE(rejects_malformed_zero_and_constant_inputs) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    auto coefficient_generator = [] { return value_type::one(); };

    polynomial_type empty;
    empty.get_storage().clear();
    BOOST_CHECK_THROW(math::recover_irreducible_polynomial_x_norm_representation<backend_type>(
                          empty, arithmetic_context, coefficient_generator),
                      std::invalid_argument);

    polynomial_type noncanonical(2);
    noncanonical[0] = value_type::one();
    noncanonical[1] = value_type::zero();
    BOOST_CHECK_THROW(math::recover_irreducible_polynomial_x_norm_representation<backend_type>(
                          noncanonical, arithmetic_context, coefficient_generator),
                      std::invalid_argument);

    BOOST_CHECK_THROW(math::recover_irreducible_polynomial_x_norm_representation<backend_type>(
                          polynomial_type {value_type::zero()}, arithmetic_context, coefficient_generator),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::recover_irreducible_polynomial_x_norm_representation<backend_type>(
                          polynomial_type {value_type::one()}, arithmetic_context, coefficient_generator),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(high_level_recovery_handles_zero_and_constant_polynomials_without_factorization) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::size_t generator_calls = 0;
    auto coefficient_generator = [&] {
        ++generator_calls;
        return value_type::one();
    };

    const auto zero = math::recover_polynomial_x_norm_representation<backend_type>(
        polynomial_type {value_type::zero()}, arithmetic_context, coefficient_generator);
    BOOST_REQUIRE(zero.has_value());
    BOOST_CHECK(zero->p == polynomial_type({value_type::zero()}));
    BOOST_CHECK(zero->q == polynomial_type({value_type::zero()}));
    BOOST_CHECK(math::evaluate_polynomial_x_norm<backend_type>(*zero, arithmetic_context) ==
                polynomial_type({value_type::zero()}));

    const polynomial_type square_constant = {value_type(9)};
    const auto square = math::recover_polynomial_x_norm_representation<backend_type>(
        square_constant, arithmetic_context, coefficient_generator);
    BOOST_REQUIRE(square.has_value());
    BOOST_CHECK(square->p.size() == 1);
    BOOST_CHECK(square->p[0] * square->p[0] == square_constant[0]);
    BOOST_CHECK(square->q == polynomial_type({value_type::zero()}));
    BOOST_CHECK(math::evaluate_polynomial_x_norm<backend_type>(*square, arithmetic_context) == square_constant);

    const polynomial_type nonsquare_constant = {first_quadratic_non_residue()};
    const auto nonsquare = math::recover_polynomial_x_norm_representation<backend_type>(
        nonsquare_constant, arithmetic_context, coefficient_generator);
    BOOST_CHECK(!nonsquare.has_value());
    BOOST_CHECK_EQUAL(generator_calls, 0);
}

BOOST_AUTO_TEST_CASE(high_level_recovery_rejects_nonconstant_inputs_that_fail_necessary_square_filters) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const value_type non_residue = first_quadratic_non_residue();
    const value_type minus_one = value_type::zero() - value_type::one();
    BOOST_REQUIRE(minus_one.is_square());
    BOOST_REQUIRE(!(value_type::zero() - non_residue).is_square());

    std::size_t generator_calls = 0;
    boost::random::mt19937 rng(0x584E1000);
    auto coefficient_generator = [&] {
        ++generator_calls;
        return nil::crypto3::algebra::random_element<field_type>(rng);
    };

    const polynomial_type nonsquare_constant_coefficient = {non_residue, value_type::zero(), value_type::one()};
    const auto constant_filter_result = math::recover_polynomial_x_norm_representation<backend_type>(
        nonsquare_constant_coefficient, arithmetic_context, coefficient_generator);
    BOOST_CHECK(!constant_filter_result.has_value());

    const polynomial_type nonsquare_even_degree_leading_coefficient = {value_type::one(), value_type::zero(),
                                                                       non_residue};
    const auto even_leading_filter_result = math::recover_polynomial_x_norm_representation<backend_type>(
        nonsquare_even_degree_leading_coefficient, arithmetic_context, coefficient_generator);
    BOOST_CHECK(!even_leading_filter_result.has_value());

    const polynomial_type nonsquare_odd_degree_signed_leading_coefficient = {value_type::one(), value_type::zero(),
                                                                             value_type::zero(), non_residue};
    const auto odd_leading_filter_result = math::recover_polynomial_x_norm_representation<backend_type>(
        nonsquare_odd_degree_signed_leading_coefficient, arithmetic_context, coefficient_generator);
    BOOST_CHECK(!odd_leading_filter_result.has_value());

    BOOST_CHECK_EQUAL(generator_calls, 0);
}

BOOST_AUTO_TEST_CASE(high_level_recovery_rejects_empty_and_noncanonical_polynomials) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    auto coefficient_generator = [] { return value_type::one(); };

    polynomial_type empty;
    empty.get_storage().clear();
    BOOST_CHECK_THROW(
        math::recover_polynomial_x_norm_representation<backend_type>(empty, arithmetic_context, coefficient_generator),
        std::invalid_argument);

    polynomial_type noncanonical(2);
    noncanonical[0] = value_type::one();
    noncanonical[1] = value_type::zero();
    BOOST_CHECK_THROW(math::recover_polynomial_x_norm_representation<backend_type>(noncanonical, arithmetic_context,
                                                                                   coefficient_generator),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(high_level_recovery_lifts_an_even_multiplicity_without_recovering_the_factor) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const value_type non_residue = first_quadratic_non_residue();
    const polynomial_type unrepresentable_factor = {value_type::zero() - non_residue, value_type::one()};
    polynomial_type h;
    arithmetic_context.square(h, unrepresentable_factor);

    boost::random::mt19937 rng(0x584E1001);
    auto coefficient_generator = [&] { return nil::crypto3::algebra::random_element<field_type>(rng); };
    const auto result =
        math::recover_polynomial_x_norm_representation<backend_type>(h, arithmetic_context, coefficient_generator);

    BOOST_REQUIRE(result.has_value());
    BOOST_CHECK(result->q == polynomial_type({value_type::zero()}));
    BOOST_CHECK(math::evaluate_polynomial_x_norm<backend_type>(*result, arithmetic_context) == h);
}

BOOST_AUTO_TEST_CASE(high_level_recovery_recovers_an_odd_multiplicity_representable_factor) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type h = representable_irreducible_quadratic();

    boost::random::mt19937 rng(0x584E1002);
    auto coefficient_generator = [&] { return nil::crypto3::algebra::random_element<field_type>(rng); };
    const auto result =
        math::recover_polynomial_x_norm_representation<backend_type>(h, arithmetic_context, coefficient_generator);

    BOOST_REQUIRE(result.has_value());
    BOOST_CHECK(math::evaluate_polynomial_x_norm<backend_type>(*result, arithmetic_context) == h);
}

BOOST_AUTO_TEST_CASE(high_level_recovery_combines_distinct_even_and_odd_factor_multiplicities) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const value_type non_residue = first_quadratic_non_residue();
    const polynomial_type even_factor = {value_type::zero() - non_residue, value_type::one()};
    polynomial_type even_factor_squared;
    arithmetic_context.square(even_factor_squared, even_factor);
    const polynomial_type odd_factor = representable_irreducible_quadratic();
    const polynomial_type h = context_multiply<backend_type>(even_factor_squared, odd_factor, arithmetic_context);

    boost::random::mt19937 rng(0x584E1003);
    auto coefficient_generator = [&] { return nil::crypto3::algebra::random_element<field_type>(rng); };
    const auto result =
        math::recover_polynomial_x_norm_representation<backend_type>(h, arithmetic_context, coefficient_generator);

    BOOST_REQUIRE(result.has_value());
    BOOST_CHECK(math::evaluate_polynomial_x_norm<backend_type>(*result, arithmetic_context) == h);
}

BOOST_AUTO_TEST_CASE(high_level_recovery_incorporates_a_square_leading_scalar) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type monic = representable_irreducible_quadratic();
    polynomial_type h;
    math::scalar_multiplication(h, monic, value_type(9));

    boost::random::mt19937 rng(0x584E1004);
    auto coefficient_generator = [&] { return nil::crypto3::algebra::random_element<field_type>(rng); };
    const auto result =
        math::recover_polynomial_x_norm_representation<backend_type>(h, arithmetic_context, coefficient_generator);

    BOOST_REQUIRE(result.has_value());
    BOOST_CHECK(math::evaluate_polynomial_x_norm<backend_type>(*result, arithmetic_context) == h);
}

BOOST_AUTO_TEST_CASE(high_level_recovery_rejects_an_unrepresentable_odd_factor_after_necessary_filters_pass) {
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const value_type first_non_residue = first_quadratic_non_residue();
    const value_type second_non_residue = next_quadratic_non_residue(first_non_residue);
    const polynomial_type first_factor = {value_type::zero() - first_non_residue, value_type::one()};
    const polynomial_type second_factor = {value_type::zero() - second_non_residue, value_type::one()};
    const polynomial_type h = context_multiply<backend_type>(first_factor, second_factor, arithmetic_context);

    BOOST_REQUIRE(h[0].is_square());
    BOOST_REQUIRE(h[h.size() - 1].is_square());
    boost::random::mt19937 rng(0x584E1005);
    auto coefficient_generator = [&] { return nil::crypto3::algebra::random_element<field_type>(rng); };
    const auto result =
        math::recover_polynomial_x_norm_representation<backend_type>(h, arithmetic_context, coefficient_generator);

    BOOST_CHECK(!result.has_value());
}

BOOST_AUTO_TEST_CASE(high_level_recovery_supports_bn254_fq12_with_the_schoolbook_backend) {
    polynomial_arithmetic::polynomial_context<fq12_schoolbook_backend> arithmetic_context;
    check_bn254_fq12_high_level_recovery(arithmetic_context, 0xF0125101);
}

BOOST_AUTO_TEST_CASE(high_level_recovery_supports_bn254_fq12_with_the_mixed_radix_backend) {
    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<fq12_mixed_radix_backend> arithmetic_context(fq12_mixed_radix_backend(18),
                                                                                           options);
    check_bn254_fq12_high_level_recovery(arithmetic_context, 0xF0125102);
}

BOOST_AUTO_TEST_CASE(screening_preserves_even_multiplicities_and_rejects_odd_obstructions) {
    polynomial_arithmetic::polynomial_context<backend_type> preparation;
    const auto non_residue = first_quadratic_non_residue();
    const auto other_non_residue = next_quadratic_non_residue(non_residue);
    const polynomial_type first = {-non_residue, value_type::one()};
    const polynomial_type second = {-other_non_residue, value_type::one()};
    const auto unsuitable_group = context_multiply<backend_type>(first, second, preparation);
    const auto good_quadratic = representable_irreducible_quadratic();
    const auto rejected_input = context_multiply<backend_type>(unsuitable_group, good_quadratic, preparation);
    polynomial_type even_power;
    preparation.square(even_power, unsuitable_group);
    auto suitable_input = context_multiply<backend_type>(even_power, good_quadratic, preparation);
    suitable_input = context_multiply<backend_type>(
        suitable_input, polynomial_type {value_type::zero(), value_type::one()}, preparation);

    // The product of two unsuitable linear factors passes the scalar filters. They must be classified individually;
    // their even powers remain admissible, including when combined with the zero-root factor X.
    BOOST_REQUIRE(rejected_input[0].is_square());
    BOOST_REQUIRE(rejected_input.back().is_square());
    for (std::size_t limit : {0, 1, 2, 100}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.distinct_degree_screening_limit = limit;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        boost::random::mt19937 rng(0x584E3001);
        auto generator = [&] { return nil::crypto3::algebra::random_element<field_type>(rng); };
        const auto rejected =
            math::recover_polynomial_x_norm_representation<backend_type>(rejected_input, context, generator);
        BOOST_CHECK(!rejected.has_value());
        const auto recovered =
            math::recover_polynomial_x_norm_representation<backend_type>(suitable_input, context, generator);
        BOOST_REQUIRE(recovered.has_value());
        BOOST_CHECK(math::evaluate_polynomial_x_norm<backend_type>(*recovered, context) == suitable_input);
    }
}

BOOST_AUTO_TEST_CASE(even_degree_trace_recovery_normalizes_scalars_and_consumes_only_retry_samples) {
    polynomial_arithmetic::polynomial_context<backend_type> context;
    // This irreducible quartic's root has trace zero; multiplying it by X gives trace -3.
    // The square-root tests independently check irreducibility and the forced retry.
    const polynomial_type monic {value_type(4), -value_type::one(), value_type(4), value_type::zero(),
                                 value_type::one()};
    for (const auto &scale : {value_type::one(), value_type(49), first_quadratic_non_residue()}) {
        polynomial_type g;
        math::scalar_multiplication(g, monic, scale);
        std::size_t calls = 0;
        auto generator = [&] {
            const std::size_t index = calls++;
            BOOST_REQUIRE_LT(index, 12);
            // Three representatives: zero, one (both have zero weighted trace), then X.
            return index == 4 || index == 9 ? value_type::one() : value_type::zero();
        };
        const auto result = math::recover_irreducible_polynomial_x_norm_representation(g, context, generator);
        BOOST_REQUIRE_EQUAL(result.has_value(), scale.is_square());
        BOOST_CHECK_EQUAL(calls, 12);
        if (result) {
            BOOST_CHECK(math::evaluate_polynomial_x_norm(*result, context) == g);
            BOOST_CHECK_LE(result->p.size() - 1, 2);
            BOOST_CHECK_LE(result->q.size() - 1, 1);
        }
    }
    std::size_t calls = 0;
    auto failing_generator = [&]() -> value_type {
        ++calls;
        throw std::runtime_error("retry coefficient generation failed");
    };
    BOOST_CHECK_THROW(math::recover_irreducible_polynomial_x_norm_representation(monic, context, failing_generator),
                      std::runtime_error);
    BOOST_CHECK_EQUAL(calls, 1);
}

BOOST_AUTO_TEST_CASE(even_degree_recovery_preserves_tonelli_shanks_when_composition_is_disabled) {
    const polynomial_type g = representable_irreducible_quadratic();
    for (const std::size_t cache_limit : {0, 1, 16}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = cache_limit;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        std::size_t calls = 0;
        boost::random::mt19937 rng(0xE7E021);
        auto generator = [&] {
            BOOST_REQUIRE_LT(++calls, 128);
            return nil::crypto3::algebra::random_element<field_type>(rng);
        };
        const auto result = math::recover_irreducible_polynomial_x_norm_representation(g, context, generator);
        BOOST_REQUIRE(result.has_value());
        BOOST_CHECK(math::evaluate_polynomial_x_norm(*result, context) == g);
        if (cache_limit == 0) {
            BOOST_CHECK_GT(calls, 0);
            BOOST_CHECK_EQUAL(calls % 2, 0);
        } else {
            BOOST_CHECK_EQUAL(calls, 0);
        }
        const std::size_t calls_before_rejection = calls;
        const polynomial_type nonsquare {-first_quadratic_non_residue(), value_type::zero(), value_type::one()};
        BOOST_CHECK(!math::recover_irreducible_polynomial_x_norm_representation(nonsquare, context, generator));
        BOOST_CHECK_EQUAL(calls, calls_before_rejection);
    }
}

BOOST_AUTO_TEST_CASE(even_degree_bn254_fp12_norm_recovery_agrees_across_backends) {
    using value = fq12_value_type;
    using polynomial = fq12_schoolbook_backend::polynomial_type;
    using fft_backend = polynomial_arithmetic::bn254_small_prime_fft_backend;
    boost::random::mt19937 rng(0xE7E022);
    value c = value::one();
    for (std::size_t i = 0; i < 64 && c.is_square(); ++i) {
        c = nil::crypto3::algebra::random_element<fq12_field_type>(rng);
    }
    BOOST_REQUIRE(!c.is_zero());
    BOOST_REQUIRE(!c.is_square());
    polynomial_arithmetic::polynomial_context<fq12_schoolbook_backend> reference;
    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<fft_backend> optimized(fft_backend(16), options);
    for (const std::size_t degree : {2, 4}) {
        value shift = value::one();
        for (std::size_t i = 0; i < 64 && !(shift.pow(degree) - c).is_square(); ++i) {
            shift = nil::crypto3::algebra::random_element<fq12_field_type>(rng);
        }
        BOOST_REQUIRE((shift.pow(degree) - c).is_square());
        polynomial g {value::one()};
        const polynomial linear {shift, value::one()};
        for (std::size_t i = 0; i < degree; ++i) {
            reference.multiply(g, g, linear);
        }
        g[0] -= c;
        const value scale = nil::crypto3::algebra::random_element<fq12_field_type>(rng).squared();
        BOOST_REQUIRE(!scale.is_zero());
        math::scalar_multiplication(g, g, scale);
        // Fix independent RNG states for equivalent fallback sampling, if a zero trace occurs.
        boost::random::mt19937 rng_reference(0xE7E023), rng_optimized(0xE7E023);
        auto reference_generator = [&] {
            return nil::crypto3::algebra::random_element<fq12_field_type>(rng_reference);
        };
        auto optimized_generator = [&] {
            return nil::crypto3::algebra::random_element<fq12_field_type>(rng_optimized);
        };
        const auto expected =
            math::recover_irreducible_polynomial_x_norm_representation(g, reference, reference_generator);
        const auto actual =
            math::recover_irreducible_polynomial_x_norm_representation(g, optimized, optimized_generator);
        BOOST_REQUIRE(expected.has_value());
        BOOST_REQUIRE(actual.has_value());
        BOOST_CHECK(math::evaluate_polynomial_x_norm(*actual, reference) == g);
        BOOST_CHECK(actual->p == expected->p);
        BOOST_CHECK(actual->q == expected->q);
        BOOST_CHECK_LE(actual->p.size() - 1, degree / 2);
        BOOST_CHECK_LE(actual->q.size() - 1, (degree - 1) / 2);
        BOOST_CHECK(rng_reference == rng_optimized);
    }
}

BOOST_AUTO_TEST_SUITE_END()
