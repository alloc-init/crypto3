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

#define BOOST_TEST_MODULE polynomial_backend_test

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/coefficient_view.hpp>
#include <nil/crypto3/math/polynomial/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/polynomial_backend.hpp>
#include <nil/crypto3/math/polynomial/schoolbook_backend.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = nil::crypto3::math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq_value_type = fq_field_type::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;

    struct incomplete_backend {
        using polynomial_type = math::polynomial<fq_value_type>;
    };

    using fq_mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type>;
    using fq12_mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;

    constexpr std::size_t odd_smooth_order = 3 * 3 * 29 * 67;
    constexpr std::size_t even_smooth_order = 2 * odd_smooth_order;

    static_assert(!polynomial_arithmetic::PolynomialBackend<incomplete_backend>);
    static_assert(polynomial_arithmetic::PolynomialBackend<polynomial_arithmetic::schoolbook_backend<fq_value_type>>);
    static_assert(polynomial_arithmetic::PolynomialBackend<polynomial_arithmetic::schoolbook_backend<fq12_value_type>>);
    static_assert(polynomial_arithmetic::PolynomialBackend<fq_mixed_radix_backend>);
    static_assert(polynomial_arithmetic::PolynomialBackend<fq12_mixed_radix_backend>);

    template<typename ValueType>
    math::polynomial<ValueType> expected_low_product(const math::polynomial<ValueType> &product,
                                                     std::size_t coefficient_count) {
        if (coefficient_count == 0) {
            return {ValueType::zero()};
        }

        math::polynomial<ValueType> result = product;
        if (result.size() > coefficient_count) {
            result.resize(coefficient_count);
        }
        nil::crypto3::math::condense(result);
        return result;
    }

    template<polynomial_arithmetic::PolynomialBackend Backend>
    void check_backend_conformance(Backend backend, const typename Backend::polynomial_type &left,
                                   const typename Backend::polynomial_type &right,
                                   const typename Backend::polynomial_type &expected_product,
                                   const typename Backend::polynomial_type &expected_square) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        polynomial_arithmetic::polynomial_context<Backend> context(std::move(backend));
        polynomial_type output = {value_type::one()};

        math::multiplication(output, left, right, context);
        BOOST_CHECK(output == expected_product);

        polynomial_type left_storage(left.size() + 2, value_type::zero());
        polynomial_type right_storage(right.size() + 2, value_type::zero());
        std::copy(left.begin(), left.end(), left_storage.begin() + 1);
        std::copy(right.begin(), right.end(), right_storage.begin() + 1);
        const math::coefficient_view<value_type> left_view =
            math::coefficient_view<value_type>(left_storage).subview(1, left.size());
        const math::coefficient_view<value_type> right_view =
            math::coefficient_view<value_type>(right_storage).subview(1, right.size());
        math::multiplication(output, left_view, right_view, context);
        BOOST_CHECK(output == expected_product);

        polynomial_type view_alias = left;
        math::multiplication(view_alias, math::coefficient_view<value_type>(view_alias), right_view, context);
        BOOST_CHECK(view_alias == expected_product);

        math::square(output, left, context);
        BOOST_CHECK(output == expected_square);

        polynomial_type left_alias = left;
        math::multiplication(left_alias, left_alias, right, context);
        BOOST_CHECK(left_alias == expected_product);

        polynomial_type right_alias = right;
        math::multiplication(right_alias, left, right_alias, context);
        BOOST_CHECK(right_alias == expected_product);

        polynomial_type both_aliases = left;
        math::multiplication(both_aliases, both_aliases, both_aliases, context);
        BOOST_CHECK(both_aliases == expected_square);

        polynomial_type square_alias = left;
        math::square(square_alias, square_alias, context);
        BOOST_CHECK(square_alias == expected_square);

        for (std::size_t coefficient_count = 0; coefficient_count <= expected_product.size() + 2; ++coefficient_count) {
            math::multiply_low(output, left, right, coefficient_count, context);
            BOOST_CHECK(output == expected_low_product(expected_product, coefficient_count));
        }

        polynomial_type low_alias = left;
        math::multiply_low(low_alias, low_alias, right, 2, context);
        BOOST_CHECK(low_alias == expected_low_product(expected_product, 2));

        const polynomial_type zero = {value_type::zero()};
        const math::coefficient_view<value_type> empty_view;
        math::multiplication(output, empty_view, right_view, context);
        BOOST_CHECK(output == zero);
        math::multiplication(output, left_view, empty_view, context);
        BOOST_CHECK(output == zero);
        math::multiplication(output, zero, right, context);
        BOOST_CHECK(output == zero);
        math::multiplication(output, left, zero, context);
        BOOST_CHECK(output == zero);
        math::multiply_low(output, zero, right, 2, context);
        BOOST_CHECK(output == zero);
        math::multiply_low(output, left, zero, 2, context);
        BOOST_CHECK(output == zero);
        math::square(output, zero, context);
        BOOST_CHECK(output == zero);
    }

    fq12_value_type fq12_value(std::size_t first_coordinate) {
        fq12_value_type value = fq12_value_type::zero();
        for (std::size_t i = 0; i < fq12_field_type::arity; ++i) {
            value.coordinate(i) = fq_value_type(first_coordinate + i);
        }
        return value;
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(polynomial_backend_test_suite)

BOOST_AUTO_TEST_CASE(fq_backends_conform) {
    using polynomial_type = math::polynomial<fq_value_type>;

    const polynomial_type left = {1, 2, 0, 3};
    const polynomial_type right = {4, 0, 5};
    const polynomial_type expected_product = {4, 8, 5, 22, 0, 15};
    const polynomial_type expected_square = {1, 4, 4, 6, 12, 0, 9};

    BOOST_TEST_CONTEXT("schoolbook") {
        check_backend_conformance(polynomial_arithmetic::schoolbook_backend<fq_value_type> {}, left, right,
                                  expected_product, expected_square);
    }
    BOOST_TEST_CONTEXT("mixed radix") {
        check_backend_conformance(fq_mixed_radix_backend(9), left, right, expected_product, expected_square);
    }
}

BOOST_AUTO_TEST_CASE(fq_backends_accept_noncanonical_coefficient_views) {
    using polynomial_type = math::polynomial<fq_value_type>;

    const polynomial_type left_storage = {1, 2, 0};
    const polynomial_type right_storage = {3, 0};
    const polynomial_type zero_storage = {0, 0, 0};
    const math::coefficient_view<fq_value_type> left(left_storage);
    const math::coefficient_view<fq_value_type> right(right_storage);
    const math::coefficient_view<fq_value_type> zero(zero_storage);
    const polynomial_type expected = {3, 6};
    const polynomial_type expected_zero = {0};

    const auto check_backend = [&](auto backend) {
        using backend_type = decltype(backend);
        polynomial_arithmetic::polynomial_context<backend_type> context(std::move(backend));
        polynomial_type output;
        math::multiplication(output, left, right, context);
        BOOST_CHECK(output == expected);
        math::multiplication(output, zero, right, context);
        BOOST_CHECK(output == expected_zero);
        math::multiplication(output, left, zero, context);
        BOOST_CHECK(output == expected_zero);
    };

    BOOST_TEST_CONTEXT("schoolbook") {
        check_backend(polynomial_arithmetic::schoolbook_backend<fq_value_type> {});
    }
    BOOST_TEST_CONTEXT("mixed radix") {
        check_backend(fq_mixed_radix_backend(9));
    }
}

BOOST_AUTO_TEST_CASE(fq12_backends_conform) {
    using polynomial_type = math::polynomial<fq12_value_type>;

    const fq12_value_type x = fq12_value(1);
    const fq12_value_type y = fq12_value(13);
    const fq12_value_type z = fq12_value(25);
    const fq12_value_type w = fq12_value(37);
    const polynomial_type left = {x, y};
    const polynomial_type right = {z, w};
    const polynomial_type expected_product = {x * z, x * w + y * z, y * w};
    const polynomial_type expected_square = {x * x, x * y + y * x, y * y};

    BOOST_TEST_CONTEXT("schoolbook") {
        check_backend_conformance(polynomial_arithmetic::schoolbook_backend<fq12_value_type> {}, left, right,
                                  expected_product, expected_square);
    }
    BOOST_TEST_CONTEXT("mixed radix") {
        check_backend_conformance(fq12_mixed_radix_backend(3), left, right, expected_product, expected_square);
    }
}

BOOST_AUTO_TEST_CASE(scalar_multiplication_supports_base_field_scalars) {
    using polynomial_type = math::polynomial<fq12_value_type>;

    const polynomial_type input = {fq12_value(1), fq12_value(13)};
    const fq_value_type scalar(7);
    polynomial_type output;

    math::scalar_multiplication(output, input, scalar);
    const polynomial_type expected = {input[0] * scalar, input[1] * scalar};
    BOOST_CHECK(output == expected);
}

BOOST_AUTO_TEST_CASE(derivative_supports_extension_field_coefficients) {
    using polynomial_type = math::polynomial<fq12_value_type>;

    const polynomial_type input = {fq12_value(1), fq12_value(13), fq12_value(25)};
    polynomial_type output;

    math::derivative(output, input);
    const polynomial_type expected = {input[1], input[2] * std::size_t(2)};
    BOOST_CHECK(output == expected);
}

BOOST_AUTO_TEST_CASE(make_monic_supports_extension_field_coefficients) {
    using polynomial_type = math::polynomial<fq12_value_type>;

    const fq12_value_type leading_coefficient = fq12_value(7);
    const polynomial_type input = {fq12_value(3), leading_coefficient};
    polynomial_type output;

    math::make_monic(output, input);
    const polynomial_type expected = {input[0] * leading_coefficient.inversed(), fq12_value_type::one()};
    BOOST_CHECK(output == expected);
}

BOOST_AUTO_TEST_CASE(mixed_radix_backend_uses_only_the_prefix_needed_by_multiply_low) {
    using polynomial_type = math::polynomial<fq_value_type>;

    fq_mixed_radix_backend backend(3);
    const polynomial_type input = {fq_value_type(1), fq_value_type(2), fq_value_type(3)};
    polynomial_type output;

    BOOST_CHECK_THROW(backend.multiply(output, math::coefficient_view(input), math::coefficient_view(input)),
                      std::invalid_argument);
    BOOST_CHECK_THROW(backend.square(output, input), std::invalid_argument);

    backend.multiply_low(output, input, input, 1);
    BOOST_CHECK(output == polynomial_type({fq_value_type(1)}));
    backend.multiply_low(output, input, input, 2);
    BOOST_CHECK(output == polynomial_type({fq_value_type(1), fq_value_type(4)}));
    BOOST_CHECK_THROW(backend.multiply_low(output, input, input, 3), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(transpose_multiplication_uses_the_selected_backend) {
    using polynomial_type = math::polynomial<fq12_value_type>;

    const polynomial_type input = {fq12_value(1), fq12_value(13)};
    const std::vector<fq_value_type> field_coefficients = {fq_value_type(7)};
    polynomial_arithmetic::polynomial_context<fq12_mixed_radix_backend> context {fq12_mixed_radix_backend(3)};

    const polynomial_type result = math::transpose_multiplication(2, input, field_coefficients, context);
    const polynomial_type expected = {input[0] * field_coefficients[0], fq12_value_type::zero(),
                                      fq12_value_type::zero()};
    BOOST_REQUIRE_EQUAL(result.size(), expected.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
        BOOST_CHECK(result[i] == expected[i]);
    }
}

BOOST_AUTO_TEST_CASE(fq12_mixed_radix_backend_uses_the_larger_smooth_order) {
    using polynomial_type = math::polynomial<fq12_value_type>;
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;

    constexpr std::size_t operand_size = odd_smooth_order / 2 + 2;
    const std::array<fq12_value_type, 6> samples = {fq12_value(1),  fq12_value(13), fq12_value(25),
                                                    fq12_value(37), fq12_value(49), fq12_value(61)};
    polynomial_type left(operand_size, fq12_value_type::zero());
    polynomial_type right(operand_size, fq12_value_type::zero());
    left[0] = samples[0];
    left[113] = samples[1];
    left.back() = samples[2];
    right[0] = samples[3];
    right[257] = samples[4];
    right.back() = samples[5];

    polynomial_arithmetic::polynomial_context<schoolbook_backend> schoolbook_context;
    polynomial_type expected;
    math::multiplication(expected, left, right, schoolbook_context);
    BOOST_REQUIRE_GT(expected.size(), odd_smooth_order);

    polynomial_arithmetic::polynomial_context<fq12_mixed_radix_backend> mixed_radix_context {
        fq12_mixed_radix_backend(even_smooth_order)};
    polynomial_type result;
    math::multiplication(result, left, right, mixed_radix_context);
    BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_SUITE_END()
