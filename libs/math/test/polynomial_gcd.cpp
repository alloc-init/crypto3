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

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/gcd.hpp>
#include <nil/crypto3/math/polynomial/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/schoolbook_backend.hpp>

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

BOOST_AUTO_TEST_SUITE_END()
