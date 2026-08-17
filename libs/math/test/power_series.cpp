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

#define BOOST_TEST_MODULE power_series_test

#include <cstddef>
#include <utility>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/power_series.hpp>
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
    typename Backend::polynomial_type compute_and_check_inverse(Backend backend,
                                                                const typename Backend::polynomial_type &input,
                                                                std::size_t coefficient_count) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        polynomial_arithmetic::polynomial_context<Backend> context(std::move(backend));
        polynomial_type inverse;
        math::inverse_series(inverse, input, coefficient_count, context);

        polynomial_type product;
        math::multiply_low(product, input, inverse, coefficient_count, context);
        BOOST_CHECK(product == polynomial_type({value_type::one()}));

        polynomial_type alias(input);
        math::inverse_series(alias, alias, coefficient_count, context);
        BOOST_CHECK(alias == inverse);
        return inverse;
    }

    fq12_value_type fq12_value(std::size_t first_coordinate) {
        fq12_value_type value = fq12_value_type::zero();
        for (std::size_t i = 0; i < fq12_field_type::arity; ++i) {
            value.coordinate(i) = fq_value_type(first_coordinate + i);
        }
        return value;
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(power_series_test_suite)

BOOST_AUTO_TEST_CASE(inverse_of_one_minus_x_is_geometric_series) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type input = {fq_value_type::one(), -fq_value_type::one()};
    for (const std::size_t coefficient_count : {std::size_t(1), std::size_t(2), std::size_t(3), std::size_t(9)}) {
        const polynomial_type inverse = compute_and_check_inverse(backend_type {}, input, coefficient_count);
        BOOST_CHECK(inverse == polynomial_type(coefficient_count, fq_value_type::one()));
    }
}

BOOST_AUTO_TEST_CASE(schoolbook_and_mixed_radix_backends_agree) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    constexpr std::size_t coefficient_count = 9;
    // The final Newton step multiplies prefixes of lengths 8 and 9. A size-18 transform is the smallest supported
    // mixed-radix transform covering their 16-coefficient product.
    constexpr std::size_t transform_size = 18;
    const polynomial_type input = {fq_value_type(3), fq_value_type(5), fq_value_type(7), fq_value_type(11)};
    const polynomial_type schoolbook_inverse =
        compute_and_check_inverse(schoolbook_backend {}, input, coefficient_count);
    const polynomial_type mixed_radix_inverse =
        compute_and_check_inverse(mixed_radix_backend(transform_size), input, coefficient_count);
    BOOST_CHECK(mixed_radix_inverse == schoolbook_inverse);
}

BOOST_AUTO_TEST_CASE(extension_field_coefficients_are_supported) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    constexpr std::size_t coefficient_count = 7;
    // The final Newton step multiplies prefixes of lengths 4 and 6, whose product has 9 coefficients.
    constexpr std::size_t transform_size = 9;
    const polynomial_type input = {fq12_value(1), fq12_value(13), fq12_value(25)};
    const polynomial_type schoolbook_inverse =
        compute_and_check_inverse(schoolbook_backend {}, input, coefficient_count);
    const polynomial_type mixed_radix_inverse =
        compute_and_check_inverse(mixed_radix_backend(transform_size), input, coefficient_count);
    BOOST_CHECK(mixed_radix_inverse == schoolbook_inverse);
}

BOOST_AUTO_TEST_CASE(zero_precision_and_noninvertible_inputs) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> context;
    polynomial_type output = {fq_value_type::one()};
    const polynomial_type zero = {fq_value_type::zero()};

    math::inverse_series(output, zero, 0, context);
    BOOST_CHECK(output == zero);

    const polynomial_type zero_constant = {fq_value_type::zero(), fq_value_type::one()};
    BOOST_CHECK_THROW(math::inverse_series(output, zero_constant, 3, context), std::invalid_argument);
    BOOST_CHECK_THROW(math::inverse_series(output, zero, 3, context), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
