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

#include <cstddef>
#include <utility>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/polynomial_backend.hpp>
#include <nil/crypto3/math/polynomial/schoolbook_backend.hpp>

namespace {
    namespace polynomial_arithmetic = nil::crypto3::math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq_value_type = fq_field_type::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;

    struct incomplete_backend {
        using value_type = fq_value_type;
    };

    static_assert(!polynomial_arithmetic::PolynomialBackend<incomplete_backend>);
    static_assert(polynomial_arithmetic::PolynomialBackend<polynomial_arithmetic::schoolbook_backend<fq_value_type>>);
    static_assert(polynomial_arithmetic::PolynomialBackend<polynomial_arithmetic::schoolbook_backend<fq12_value_type>>);

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

        context.multiply(output, left, right);
        BOOST_CHECK(output == expected_product);

        context.square(output, left);
        BOOST_CHECK(output == expected_square);

        coefficient_vector left_alias = left;
        context.multiply(left_alias, left_alias, right);
        BOOST_CHECK(left_alias == expected_product);

        coefficient_vector right_alias = right;
        context.multiply(right_alias, left, right_alias);
        BOOST_CHECK(right_alias == expected_product);

        coefficient_vector both_aliases = left;
        context.multiply(both_aliases, both_aliases, both_aliases);
        BOOST_CHECK(both_aliases == expected_square);

        coefficient_vector square_alias = left;
        context.square(square_alias, square_alias);
        BOOST_CHECK(square_alias == expected_square);

        for (std::size_t coefficient_count = 0; coefficient_count <= expected_product.size() + 2; ++coefficient_count) {
            context.multiply_low(output, left, right, coefficient_count);
            BOOST_CHECK(output == expected_low_product(expected_product, coefficient_count));
        }

        coefficient_vector low_alias = left;
        context.multiply_low(low_alias, low_alias, right, 2);
        BOOST_CHECK(low_alias == expected_low_product(expected_product, 2));

        const coefficient_vector zero = {value_type::zero()};
        context.multiply(output, zero, right);
        BOOST_CHECK(output == zero);
        context.multiply(output, left, zero);
        BOOST_CHECK(output == zero);
        context.multiply_low(output, zero, right, 2);
        BOOST_CHECK(output == zero);
        context.multiply_low(output, left, zero, 2);
        BOOST_CHECK(output == zero);
        context.square(output, zero);
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

BOOST_AUTO_TEST_CASE(schoolbook_fq_backend_conforms) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using coefficient_vector = polynomial_arithmetic::coefficient_vector<fq_value_type>;

    const coefficient_vector left = {1, 2, 0, 3};
    const coefficient_vector right = {4, 0, 5};
    const coefficient_vector expected_product = {4, 8, 5, 22, 0, 15};
    const coefficient_vector expected_square = {1, 4, 4, 6, 12, 0, 9};

    check_backend_conformance(backend_type {}, left, right, expected_product, expected_square);
}

BOOST_AUTO_TEST_CASE(schoolbook_fq12_backend_conforms) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using coefficient_vector = polynomial_arithmetic::coefficient_vector<fq12_value_type>;

    const fq12_value_type x = fq12_value(1);
    const fq12_value_type y = fq12_value(13);
    const fq12_value_type z = fq12_value(25);
    const fq12_value_type w = fq12_value(37);
    const coefficient_vector left = {x, y};
    const coefficient_vector right = {z, w};
    const coefficient_vector expected_product = {x * z, x * w + y * z, y * w};
    const coefficient_vector expected_square = {x * x, x * y + y * x, y * y};

    check_backend_conformance(backend_type {}, left, right, expected_product, expected_square);
}

BOOST_AUTO_TEST_SUITE_END()
