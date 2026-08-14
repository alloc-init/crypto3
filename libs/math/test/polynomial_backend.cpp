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

#include <array>
#include <cstddef>
#include <utility>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

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
        using value_type = fq_value_type;
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
    polynomial_arithmetic::coefficient_vector<ValueType>
        expected_low_product(const polynomial_arithmetic::coefficient_vector<ValueType> &product,
                             std::size_t coefficient_count) {
        if (coefficient_count == 0) {
            return {ValueType::zero()};
        }

        polynomial_arithmetic::coefficient_vector<ValueType> result = product;
        if (result.size() > coefficient_count) {
            result.resize(coefficient_count);
        }
        nil::crypto3::math::condense(result);
        return result;
    }

    template<polynomial_arithmetic::PolynomialBackend Backend>
    void check_backend_conformance(
        Backend backend,
        const polynomial_arithmetic::coefficient_vector<typename Backend::value_type> &left,
        const polynomial_arithmetic::coefficient_vector<typename Backend::value_type> &right,
        const polynomial_arithmetic::coefficient_vector<typename Backend::value_type> &expected_product,
        const polynomial_arithmetic::coefficient_vector<typename Backend::value_type> &expected_square) {
        using value_type = typename Backend::value_type;
        using coefficient_vector = polynomial_arithmetic::coefficient_vector<value_type>;

        polynomial_arithmetic::polynomial_context<Backend> context(std::move(backend));
        coefficient_vector output = {value_type::one()};

        math::multiplication(output, left, right, context);
        BOOST_CHECK(output == expected_product);

        math::square(output, left, context);
        BOOST_CHECK(output == expected_square);

        coefficient_vector left_alias = left;
        math::multiplication(left_alias, left_alias, right, context);
        BOOST_CHECK(left_alias == expected_product);

        coefficient_vector right_alias = right;
        math::multiplication(right_alias, left, right_alias, context);
        BOOST_CHECK(right_alias == expected_product);

        coefficient_vector both_aliases = left;
        math::multiplication(both_aliases, both_aliases, both_aliases, context);
        BOOST_CHECK(both_aliases == expected_square);

        coefficient_vector square_alias = left;
        math::square(square_alias, square_alias, context);
        BOOST_CHECK(square_alias == expected_square);

        for (std::size_t coefficient_count = 0; coefficient_count <= expected_product.size() + 2; ++coefficient_count) {
            math::multiply_low(output, left, right, coefficient_count, context);
            BOOST_CHECK(output == expected_low_product(expected_product, coefficient_count));
        }

        coefficient_vector low_alias = left;
        math::multiply_low(low_alias, low_alias, right, 2, context);
        BOOST_CHECK(low_alias == expected_low_product(expected_product, 2));

        const coefficient_vector zero = {value_type::zero()};
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
    using coefficient_vector = polynomial_arithmetic::coefficient_vector<fq_value_type>;

    const coefficient_vector left = {1, 2, 0, 3};
    const coefficient_vector right = {4, 0, 5};
    const coefficient_vector expected_product = {4, 8, 5, 22, 0, 15};
    const coefficient_vector expected_square = {1, 4, 4, 6, 12, 0, 9};

    BOOST_TEST_CONTEXT("schoolbook") {
        check_backend_conformance(polynomial_arithmetic::schoolbook_backend<fq_value_type> {}, left, right,
                                  expected_product, expected_square);
    }
    BOOST_TEST_CONTEXT("mixed radix") {
        check_backend_conformance(fq_mixed_radix_backend(9), left, right, expected_product, expected_square);
    }
}

BOOST_AUTO_TEST_CASE(fq12_backends_conform) {
    using coefficient_vector = polynomial_arithmetic::coefficient_vector<fq12_value_type>;

    const fq12_value_type x = fq12_value(1);
    const fq12_value_type y = fq12_value(13);
    const fq12_value_type z = fq12_value(25);
    const fq12_value_type w = fq12_value(37);
    const coefficient_vector left = {x, y};
    const coefficient_vector right = {z, w};
    const coefficient_vector expected_product = {x * z, x * w + y * z, y * w};
    const coefficient_vector expected_square = {x * x, x * y + y * x, y * y};

    BOOST_TEST_CONTEXT("schoolbook") {
        check_backend_conformance(polynomial_arithmetic::schoolbook_backend<fq12_value_type> {}, left, right,
                                  expected_product, expected_square);
    }
    BOOST_TEST_CONTEXT("mixed radix") {
        check_backend_conformance(fq12_mixed_radix_backend(3), left, right, expected_product, expected_square);
    }
}

BOOST_AUTO_TEST_CASE(mixed_radix_backend_uses_only_the_prefix_needed_by_multiply_low) {
    using coefficient_vector = polynomial_arithmetic::coefficient_vector<fq_value_type>;

    fq_mixed_radix_backend backend(3);
    const coefficient_vector input = {fq_value_type(1), fq_value_type(2), fq_value_type(3)};
    coefficient_vector output;

    BOOST_CHECK_THROW(backend.multiply(output, input, input), std::invalid_argument);
    BOOST_CHECK_THROW(backend.square(output, input), std::invalid_argument);

    backend.multiply_low(output, input, input, 1);
    BOOST_CHECK(output == coefficient_vector({fq_value_type(1)}));
    backend.multiply_low(output, input, input, 2);
    BOOST_CHECK(output == coefficient_vector({fq_value_type(1), fq_value_type(4)}));
    BOOST_CHECK_THROW(backend.multiply_low(output, input, input, 3), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(transpose_multiplication_uses_the_selected_backend) {
    using coefficient_vector = polynomial_arithmetic::coefficient_vector<fq12_value_type>;

    const coefficient_vector input = {fq12_value(1), fq12_value(13)};
    const std::vector<fq_value_type> field_coefficients = {fq_value_type(7)};
    polynomial_arithmetic::polynomial_context<fq12_mixed_radix_backend> context {fq12_mixed_radix_backend(3)};

    const coefficient_vector result = math::transpose_multiplication(2, input, field_coefficients, context);
    const coefficient_vector expected = {input[0] * field_coefficients[0], fq12_value_type::zero(),
                                         fq12_value_type::zero()};
    BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_CASE(fq12_mixed_radix_backend_uses_the_larger_smooth_order) {
    using coefficient_vector = polynomial_arithmetic::coefficient_vector<fq12_value_type>;
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;

    constexpr std::size_t operand_size = odd_smooth_order / 2 + 2;
    const std::array<fq12_value_type, 6> samples = {fq12_value(1),  fq12_value(13), fq12_value(25),
                                                    fq12_value(37), fq12_value(49), fq12_value(61)};
    coefficient_vector left(operand_size, fq12_value_type::zero());
    coefficient_vector right(operand_size, fq12_value_type::zero());
    left[0] = samples[0];
    left[113] = samples[1];
    left.back() = samples[2];
    right[0] = samples[3];
    right[257] = samples[4];
    right.back() = samples[5];

    polynomial_arithmetic::polynomial_context<schoolbook_backend> schoolbook_context;
    coefficient_vector expected;
    math::multiplication(expected, left, right, schoolbook_context);
    BOOST_REQUIRE_GT(expected.size(), odd_smooth_order);

    polynomial_arithmetic::polynomial_context<fq12_mixed_radix_backend> mixed_radix_context {
        fq12_mixed_radix_backend(even_smooth_order)};
    coefficient_vector result;
    math::multiplication(result, left, right, mixed_radix_context);
    BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_SUITE_END()
