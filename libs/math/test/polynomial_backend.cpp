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
#include <limits>
#include <utility>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/polynomial_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>

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

    template<typename Context>
    concept HasPreparedLowProduct = requires(Context &context, const typename Context::polynomial_type &fixed) {
        context.prepare_low_product(fixed, std::size_t(2), std::size_t(3));
    };

    template<typename Context, typename Prepared>
    concept UsesPreparedLowProduct =
        requires(Context &context, typename Context::polynomial_type &output,
                 const typename Context::polynomial_type &input, const Prepared &prepared) {
            { context.try_multiply_low_prepared(output, input, prepared, std::size_t(3)) } -> std::same_as<bool>;
        };

    using fq_context = polynomial_arithmetic::polynomial_context<fq_mixed_radix_backend>;
    using schoolbook_context =
        polynomial_arithmetic::polynomial_context<polynomial_arithmetic::schoolbook_backend<fq_value_type>>;
    static_assert(HasPreparedLowProduct<fq_context>);
    static_assert(!HasPreparedLowProduct<schoolbook_context>);
    static_assert(UsesPreparedLowProduct<fq_context, fq_mixed_radix_backend::prepared_low_product_type>);
    static_assert(!UsesPreparedLowProduct<schoolbook_context, fq_mixed_radix_backend::prepared_low_product_type>);
    static_assert(!UsesPreparedLowProduct<fq_context, int>);

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

    template<typename Backend>
    void check_prepared_low_products(Backend backend, const typename Backend::polynomial_type &left,
                                     const typename Backend::polynomial_type &fixed) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename Backend::value_type;
        polynomial_arithmetic::polynomial_context<Backend> context(std::move(backend));
        polynomial_arithmetic::schoolbook_backend<value_type> reference;
        for (std::size_t count = 0; count <= left.size() + fixed.size() + 1; ++count) {
            const auto prepared = context.prepare_low_product(fixed, left.size(), count);
            BOOST_REQUIRE(prepared.has_value());
            BOOST_CHECK_EQUAL(prepared->coefficient_count(), count);
            polynomial_type variable = left;
            for (int repeat = 0; repeat < 3; ++repeat) {
                polynomial_type full_product, ordinary, result;
                reference.multiply(full_product, variable, fixed);
                const auto expected = expected_low_product(full_product, count);
                context.multiply_low(ordinary, variable, fixed, count);
                BOOST_CHECK(ordinary == expected);
                BOOST_REQUIRE(context.try_multiply_low_prepared(result, variable, *prepared, count));
                BOOST_CHECK(result == expected);
                result = variable;
                BOOST_REQUIRE(context.try_multiply_low_prepared(result, result, *prepared, count));
                BOOST_CHECK(result == expected);
                // Different values with the same prefix length reuse the same immutable spectrum.
                variable[0] += value_type::one();
            }
        }
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

BOOST_AUTO_TEST_CASE(prepared_low_products_match_schoolbook_and_ordinary_products) {
    check_prepared_low_products(fq_mixed_radix_backend(18), fq_mixed_radix_backend::polynomial_type {1, 2, 0, 3},
                                fq_mixed_radix_backend::polynomial_type {0, 4, 5});
    check_prepared_low_products(
        fq12_mixed_radix_backend(18),
        fq12_mixed_radix_backend::polynomial_type {fq12_value(1), fq12_value_type::zero(), fq12_value(13)},
        fq12_mixed_radix_backend::polynomial_type {fq12_value_type::zero(), fq12_value(25), fq12_value(37)});
}

BOOST_AUTO_TEST_CASE(prepared_low_products_handle_zero_and_truncated_prefixes) {
    using polynomial_type = fq_mixed_radix_backend::polynomial_type;
    fq_context context {fq_mixed_radix_backend(3)};
    const polynomial_type input {1, 2, 3};
    const polynomial_type zero {fq_value_type::zero()};
    const polynomial_type zero_prefix {0, 0, 7};
    polynomial_type result {fq_value_type(9)};

    // The full product needs five coefficients, but the selected prefixes need only a size-three FFT.
    auto prepared = context.prepare_low_product(input, input.size(), 2);
    BOOST_REQUIRE(prepared.has_value());
    BOOST_CHECK_EQUAL(prepared->transform_size(), 3);
    BOOST_REQUIRE(context.try_multiply_low_prepared(result, input, *prepared, 2));
    BOOST_CHECK(result == polynomial_type({1, 4}));
    BOOST_REQUIRE(context.try_multiply_low_prepared(result, zero, *prepared, 2));
    BOOST_CHECK(result == zero);

    for (const auto &fixed : {zero, zero_prefix}) {
        const auto zero_prepared = context.prepare_low_product(fixed, input.size(), 2);
        BOOST_REQUIRE(zero_prepared.has_value());
        BOOST_CHECK_EQUAL(zero_prepared->transform_size(), 0);
        BOOST_CHECK_EQUAL(zero_prepared->storage_bytes(), 0);
        BOOST_REQUIRE(context.try_multiply_low_prepared(result, input, *zero_prepared, 2));
        BOOST_CHECK(result == zero);
    }
    prepared = context.prepare_low_product(input, input.size(), 0);
    BOOST_REQUIRE(prepared.has_value());
    BOOST_CHECK_EQUAL(prepared->transform_size(), 0);
    BOOST_REQUIRE(context.try_multiply_low_prepared(result, input, *prepared, 0));
    BOOST_CHECK(result == zero);

    // A nonzero prefix may end in zero coefficients; final output must still be canonical.
    const auto trailing_zero = context.prepare_low_product(polynomial_type {2, 0, 7}, 2, 2);
    BOOST_REQUIRE(trailing_zero.has_value());
    BOOST_REQUIRE(context.try_multiply_low_prepared(result, polynomial_type {3, 0, 1}, *trailing_zero, 2));
    BOOST_CHECK(result == polynomial_type({fq_value_type(6)}));
}

BOOST_AUTO_TEST_CASE(prepared_low_products_validate_precision_and_transform_shape) {
    using polynomial_type = fq_mixed_radix_backend::polynomial_type;
    fq_context context {fq_mixed_radix_backend(18)};
    const polynomial_type fixed {2, 0, 1};
    const auto prepared = context.prepare_low_product(fixed, 2, 10);
    BOOST_REQUIRE(prepared.has_value());
    BOOST_CHECK_EQUAL(prepared->transform_size(), 6);

    // Several variable-prefix lengths select the same transform and may share this spectrum.
    polynomial_arithmetic::schoolbook_backend<fq_value_type> reference;
    for (const auto &left : {polynomial_type {1, 3}, polynomial_type {1, 3, 4}, polynomial_type {1, 3, 4, 5}}) {
        polynomial_type result, expected;
        reference.multiply(expected, left, fixed);
        BOOST_REQUIRE(context.try_multiply_low_prepared(result, left, *prepared, 10));
        BOOST_CHECK(result == expected);
    }

    // A different precision cannot reuse a spectrum containing a different fixed prefix. A different
    // selected length also cannot reuse it, even if the configured maximum can hold the product.
    polynomial_type result {fq_value_type(11)};
    BOOST_CHECK(!context.try_multiply_low_prepared(result, polynomial_type {1, 3}, *prepared, 9));
    BOOST_CHECK(result == polynomial_type({fq_value_type(11)}));
    for (const auto &input : {polynomial_type {fq_value_type::one()}, polynomial_type {1, 2, 3, 4, 5}}) {
        result = input;
        BOOST_CHECK(!context.try_multiply_low_prepared(result, result, *prepared, 10));
        BOOST_CHECK(result == input);
    }
    fq_context different_plan {fq_mixed_radix_backend(9)};
    BOOST_CHECK(!different_plan.try_multiply_low_prepared(result, polynomial_type {1, 3}, *prepared, 10));
    fq_context too_small {fq_mixed_radix_backend(3)};
    BOOST_CHECK(!too_small.try_multiply_low_prepared(result, polynomial_type {1, 3}, *prepared, 10));
    BOOST_CHECK(!too_small.prepare_low_product(fixed, 3, 3).has_value());
    BOOST_CHECK(!context.prepare_low_product(fixed, 0, 3).has_value());
    BOOST_CHECK(!context
                     .prepare_low_product(fixed, std::numeric_limits<std::size_t>::max(),
                                          std::numeric_limits<std::size_t>::max())
                     .has_value());
}

BOOST_AUTO_TEST_CASE(prepared_low_products_own_the_fixed_operand_and_survive_context_changes) {
    using polynomial_type = fq_mixed_radix_backend::polynomial_type;
    const polynomial_type original {2, 0, 1};
    const polynomial_type input {1, 3};
    const auto prepared = [&] {
        fq_context temporary {fq_mixed_radix_backend(6)};
        polynomial_type source = original;
        auto result = temporary.prepare_low_product(source, input.size(), 10);
        source[0] = fq_value_type(99);
        return result;
    }();
    BOOST_REQUIRE(prepared.has_value());
    auto copy = *prepared;
    const auto moved = std::move(copy);
    fq_context context {fq_mixed_radix_backend(18)};
    for (const auto *fixed : {&*prepared, &moved}) {
        polynomial_type result, expected;
        context.multiply_low(expected, input, original, 10);
        BOOST_REQUIRE(context.try_multiply_low_prepared(result, input, *fixed, 10));
        BOOST_CHECK(result == expected);
    }
    const auto changed = context.prepare_low_product(polynomial_type {99, 0, 1}, input.size(), 10);
    BOOST_REQUIRE(changed.has_value());
    polynomial_type result;
    BOOST_REQUIRE(context.try_multiply_low_prepared(result, input, *changed, 10));
    BOOST_CHECK(result == polynomial_type({99, 297, 1, 3}));
    BOOST_REQUIRE(context.try_multiply_low_prepared(result, input, *prepared, 10));
    BOOST_CHECK(result == polynomial_type({2, 6, 1, 3}));
}

BOOST_AUTO_TEST_CASE(prepared_low_products_do_not_retain_larger_workspace_allocations) {
    using polynomial_type = fq_mixed_radix_backend::polynomial_type;
    fq_mixed_radix_backend backend(134);
    polynomial_type large(60, fq_value_type::one()), result;
    backend.square(result, large);    // Populate the size-134 workspace first.
    const polynomial_type small(16, fq_value_type::one());
    const auto prepared = backend.prepare_low_product(small, small.size(), small.size());
    BOOST_REQUIRE(prepared.has_value());
    BOOST_CHECK_EQUAL(prepared->transform_size(), 67);
    // The long-lived spectrum owns a fresh size-67 allocation instead of the size-134 scratch buffer.
    BOOST_CHECK_LT(prepared->storage_bytes(), 134 * sizeof(fq_value_type));
    BOOST_REQUIRE(backend.try_multiply_low_prepared(result, small, *prepared, small.size()));
    BOOST_REQUIRE_EQUAL(result.size(), small.size());
    for (std::size_t i = 0; i < small.size(); ++i) {
        BOOST_CHECK(result[i] == fq_value_type(i + 1));
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

    BOOST_CHECK_THROW(backend.multiply(output, input, input), std::invalid_argument);
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
