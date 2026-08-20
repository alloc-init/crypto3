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

#define BOOST_TEST_MODULE square_free_factorization_test

#include <cstddef>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/detail/element/fp.hpp>
#include <nil/crypto3/algebra/fields/field.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/algebra/fields/params.hpp>

#include <nil/crypto3/math/polynomial/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/square_free_factorization.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq_value_type = fq_field_type::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;

    using namespace boost::multiprecision::literals;

    class test_prime_field_7 : public fields::field<4> {
    public:
        using policy_type = fields::field<4>;
        using integral_type = policy_type::integral_type;

        constexpr static std::size_t value_bits = modulus_bits;
        constexpr static std::size_t arity = 1;
        constexpr static integral_type modulus = 0x7_cppui_modular4;
        constexpr static integral_type group_order_minus_one_half = (modulus - 1u) / 2;
        constexpr static const modular_params_type modulus_params = modulus.backend();

        using modular_type = boost::multiprecision::number<boost::multiprecision::backends::modular_adaptor<
            modular_backend, boost::multiprecision::backends::modular_params_ct<modular_backend, modulus_params>>>;
        using value_type = fields::detail::element_fp<fields::params<test_prime_field_7>>;
    };

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
    typename Backend::polynomial_type power(Backend &backend, const typename Backend::polynomial_type &base,
                                            std::size_t exponent) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        polynomial_type result = {value_type::one()};
        for (std::size_t i = 0; i < exponent; ++i) {
            result = multiply(backend, result, base);
        }
        return result;
    }

    template<polynomial_arithmetic::PolynomialBackend Backend>
    typename Backend::polynomial_type
        reconstruct(Backend &backend,
                    const math::polynomial_factorization_result<typename Backend::polynomial_type> &result) {
        using polynomial_type = typename Backend::polynomial_type;

        polynomial_type reconstructed = {result.leading_coefficient};
        for (const auto &factor : result.factors) {
            for (std::size_t i = 0; i < factor.multiplicity; ++i) {
                reconstructed = multiply(backend, reconstructed, factor.polynomial);
            }
        }
        return reconstructed;
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(square_free_factorization_test_suite)

BOOST_AUTO_TEST_CASE(fq_factorization_preserves_the_leading_coefficient_and_multiplicities) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type multiplicity_one = {fq_value_type(4), fq_value_type::one()};
    const polynomial_type multiplicity_two = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type multiplicity_three = {fq_value_type(2), fq_value_type::zero(), fq_value_type::one()};

    polynomial_type input = multiply(backend, multiplicity_one, power(backend, multiplicity_two, 2));
    input = multiply(backend, input, power(backend, multiplicity_three, 3));
    math::scalar_multiplication(input, input, fq_value_type(7));

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const auto result = math::square_free_factorization<backend_type>(input, arithmetic_context);

    BOOST_CHECK(result.complete);
    BOOST_CHECK(result.leading_coefficient == fq_value_type(7));
    BOOST_REQUIRE_EQUAL(result.factors.size(), 3);
    BOOST_CHECK(result.factors[0] == math::polynomial_factor<polynomial_type>({multiplicity_one, 1}));
    BOOST_CHECK(result.factors[1] == math::polynomial_factor<polynomial_type>({multiplicity_two, 2}));
    BOOST_CHECK(result.factors[2] == math::polynomial_factor<polynomial_type>({multiplicity_three, 3}));
    BOOST_CHECK(reconstruct(backend, result) == input);
}

BOOST_AUTO_TEST_CASE(staged_factorization_stops_after_the_reported_factor) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type multiplicity_one = {fq_value_type(4), fq_value_type::one()};
    const polynomial_type multiplicity_two = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type input = multiply(backend, multiplicity_one, power(backend, multiplicity_two, 2));

    std::size_t callback_count = 0;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const auto result = math::square_free_factorization<backend_type>(
        input, arithmetic_context, [&](const math::polynomial_factor<polynomial_type> &factor) {
            ++callback_count;
            BOOST_CHECK(factor.polynomial == multiplicity_one);
            BOOST_CHECK_EQUAL(factor.multiplicity, 1);
            return math::factorization_control::stop_factorization;
        });

    BOOST_CHECK(!result.complete);
    BOOST_CHECK_EQUAL(callback_count, 1);
    BOOST_REQUIRE_EQUAL(result.factors.size(), 1);
    BOOST_CHECK(result.factors.front().polynomial == multiplicity_one);
    BOOST_CHECK_EQUAL(result.factors.front().multiplicity, 1);
}

BOOST_AUTO_TEST_CASE(zero_and_constant_polynomials_are_represented_by_the_leading_coefficient) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const auto zero =
        math::square_free_factorization<backend_type>(polynomial_type {fq_value_type::zero()}, arithmetic_context);
    BOOST_CHECK(zero.complete);
    BOOST_CHECK(zero.leading_coefficient == fq_value_type::zero());
    BOOST_CHECK(zero.factors.empty());

    const auto constant =
        math::square_free_factorization<backend_type>(polynomial_type {fq_value_type(13)}, arithmetic_context);
    BOOST_CHECK(constant.complete);
    BOOST_CHECK(constant.leading_coefficient == fq_value_type(13));
    BOOST_CHECK(constant.factors.empty());
}

BOOST_AUTO_TEST_CASE(fq12_factorization_uses_native_extension_field_coefficients) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type multiplicity_two = {fq12_value(1), fq12_value_type::one()};
    const polynomial_type multiplicity_three = {fq12_value(20), fq12_value_type::one()};
    const fq12_value_type leading_coefficient = fq12_value(40);

    polynomial_type input =
        multiply(backend, power(backend, multiplicity_two, 2), power(backend, multiplicity_three, 3));
    math::scalar_multiplication(input, input, leading_coefficient);

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const auto result = math::square_free_factorization<backend_type>(input, arithmetic_context);

    BOOST_CHECK(result.complete);
    BOOST_CHECK(result.leading_coefficient == leading_coefficient);
    BOOST_REQUIRE_EQUAL(result.factors.size(), 2);
    BOOST_CHECK(result.factors[0] == math::polynomial_factor<polynomial_type>({multiplicity_two, 2}));
    BOOST_CHECK(result.factors[1] == math::polynomial_factor<polynomial_type>({multiplicity_three, 3}));
    BOOST_CHECK(reconstruct(backend, result) == input);

    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    polynomial_arithmetic::polynomial_context_options options;
    // Force fast division so the test exercises mixed-radix arithmetic inside the factorization algorithm.
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<mixed_radix_backend> mixed_radix_context(mixed_radix_backend(18),
                                                                                       options);
    const auto mixed_radix_result = math::square_free_factorization<mixed_radix_backend>(input, mixed_radix_context);
    BOOST_CHECK(mixed_radix_result == result);
}

BOOST_AUTO_TEST_CASE(characteristic_must_be_greater_than_the_polynomial_degree) {
    using value_type = test_prime_field_7::value_type;
    using backend_type = polynomial_arithmetic::schoolbook_backend<value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    polynomial_type degree_six(7, value_type::zero());
    degree_six[6] = value_type::one();
    const auto supported = math::square_free_factorization<backend_type>(degree_six, arithmetic_context);
    BOOST_REQUIRE_EQUAL(supported.factors.size(), 1);
    BOOST_CHECK(supported.factors.front().polynomial == polynomial_type({value_type::zero(), value_type::one()}));
    BOOST_CHECK_EQUAL(supported.factors.front().multiplicity, 6);

    polynomial_type degree_seven(8, value_type::zero());
    degree_seven[0] = value_type::one();
    degree_seven[7] = value_type::one();
    BOOST_CHECK_THROW(math::square_free_factorization<backend_type>(degree_seven, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
