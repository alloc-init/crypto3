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
                                                              const typename Backend::polynomial_type &modulus,
                                                              std::size_t inverse_precision) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        polynomial_arithmetic::polynomial_context<Backend> arithmetic_context(std::move(backend));
        math::polynomial_modulus_context<Backend> modulus_context(modulus, inverse_precision, arithmetic_context);

        polynomial_type canonical_modulus(modulus);
        math::condense(canonical_modulus);
        BOOST_CHECK(modulus_context.modulus() == canonical_modulus);
        BOOST_CHECK_EQUAL(modulus_context.degree(), canonical_modulus.size() - 1);
        BOOST_CHECK_EQUAL(modulus_context.inverse_precision(), inverse_precision);

        polynomial_type reversed_modulus(canonical_modulus);
        math::reverse(reversed_modulus, reversed_modulus.size());
        polynomial_type product;
        math::multiply_low(product, reversed_modulus, modulus_context.reversed_modulus_inverse(), inverse_precision,
                           arithmetic_context);
        BOOST_CHECK(product == polynomial_type({value_type::one()}));
        return modulus_context.reversed_modulus_inverse();
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

BOOST_AUTO_TEST_CASE(modulus_is_canonical_and_reversed_inverse_is_correct) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type modulus = {fq_value_type(2), fq_value_type(3), fq_value_type(4), fq_value_type::zero()};
    build_and_check_context(backend_type {}, modulus, 7);
}

BOOST_AUTO_TEST_CASE(nonzero_constant_modulus_is_supported) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type modulus = {fq_value_type(7)};
    build_and_check_context(backend_type {}, modulus, 5);
}

BOOST_AUTO_TEST_CASE(schoolbook_and_mixed_radix_precomputation_agree) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    constexpr std::size_t inverse_precision = 9;
    // The final Newton step multiplies prefixes of lengths 8 and 9. A size-18 transform is the smallest supported
    // mixed-radix transform covering their 16-coefficient product.
    constexpr std::size_t transform_size = 18;
    const polynomial_type modulus = {fq_value_type(2), fq_value_type(3), fq_value_type(5), fq_value_type(7)};
    const polynomial_type schoolbook_inverse =
        build_and_check_context(schoolbook_backend {}, modulus, inverse_precision);
    const polynomial_type mixed_radix_inverse =
        build_and_check_context(mixed_radix_backend(transform_size), modulus, inverse_precision);
    BOOST_CHECK(mixed_radix_inverse == schoolbook_inverse);
}

BOOST_AUTO_TEST_CASE(extension_field_modulus_is_supported) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    constexpr std::size_t inverse_precision = 7;
    // The final Newton step multiplies prefixes of lengths 4 and 6, whose product has 9 coefficients.
    constexpr std::size_t transform_size = 9;
    const polynomial_type modulus = {fq12_value(1), fq12_value(13), fq12_value(25)};
    const polynomial_type schoolbook_inverse =
        build_and_check_context(schoolbook_backend {}, modulus, inverse_precision);
    const polynomial_type mixed_radix_inverse =
        build_and_check_context(mixed_radix_backend(transform_size), modulus, inverse_precision);
    BOOST_CHECK(mixed_radix_inverse == schoolbook_inverse);
}

BOOST_AUTO_TEST_CASE(zero_modulus_is_rejected) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type zero = {fq_value_type::zero(), fq_value_type::zero()};
    BOOST_CHECK_THROW(math::polynomial_modulus_context<backend_type>(zero, 3, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(zero_inverse_precision_is_rejected) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type modulus = {fq_value_type::one(), fq_value_type::one()};
    BOOST_CHECK_THROW(math::polynomial_modulus_context<backend_type>(modulus, 0, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
