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

#define BOOST_TEST_MODULE coefficient_view_test

#include <concepts>
#include <span>
#include <utility>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/math/polynomial/coefficient_view.hpp>
#include <nil/crypto3/math/polynomial/polynomial.hpp>

namespace math = nil::crypto3::math;
namespace fields = nil::crypto3::algebra::fields;

using value_type = fields::babybear::value_type;
using polynomial_type = math::polynomial<value_type>;
using view_type = math::coefficient_view<value_type>;

static_assert(std::constructible_from<view_type, polynomial_type &>);
static_assert(std::constructible_from<view_type, const polynomial_type &>);
static_assert(!std::constructible_from<view_type, polynomial_type &&>);
static_assert(std::constructible_from<view_type, std::span<value_type>>);
static_assert(std::constructible_from<view_type, std::span<const value_type>>);
static_assert(std::constructible_from<view_type, std::span<value_type, 3>>);
static_assert(std::constructible_from<view_type, std::span<const value_type, 3>>);
static_assert(std::same_as<decltype(std::declval<const view_type &>()[0]), const value_type &>);
static_assert(!math::CoefficientPolynomial<view_type>);

BOOST_AUTO_TEST_SUITE(coefficient_view_test_suite)

BOOST_AUTO_TEST_CASE(views_an_owning_polynomial_without_copying) {
    polynomial_type polynomial = {1, 2, 3, 4, 5};
    math::coefficient_view view(polynomial);

    BOOST_CHECK_EQUAL(view.data(), polynomial.data());
    BOOST_CHECK_EQUAL(view.size(), polynomial.size());
    for (std::size_t i = 0; i < polynomial.size(); ++i) {
        BOOST_CHECK_EQUAL(view[i], polynomial[i]);
    }

    polynomial[2] = value_type(17);
    BOOST_CHECK_EQUAL(view[2], value_type(17));
}

BOOST_AUTO_TEST_CASE(creates_constant_time_subviews) {
    const polynomial_type polynomial = {1, 2, 3, 4, 5};
    const view_type view(polynomial);

    const view_type first = view.first(3);
    BOOST_CHECK_EQUAL(first.data(), polynomial.data());
    BOOST_CHECK_EQUAL(first.size(), 3);
    BOOST_CHECK_EQUAL(first[2], value_type(3));

    const view_type last = view.last(2);
    BOOST_CHECK_EQUAL(last.data(), polynomial.data() + 3);
    BOOST_CHECK_EQUAL(last.size(), 2);
    BOOST_CHECK_EQUAL(last[0], value_type(4));

    const view_type middle = view.subview(1, 3);
    BOOST_CHECK_EQUAL(middle.data(), polynomial.data() + 1);
    BOOST_CHECK_EQUAL(middle.size(), 3);
    BOOST_CHECK_EQUAL(middle[0], value_type(2));
    BOOST_CHECK_EQUAL(middle[2], value_type(4));
}

BOOST_AUTO_TEST_CASE(represents_empty_coefficient_ranges) {
    const view_type default_view;
    BOOST_CHECK(default_view.empty());
    BOOST_CHECK_EQUAL(default_view.size(), 0);

    const polynomial_type polynomial = {1, 2, 3};
    const view_type view(polynomial);

    const view_type empty_prefix = view.first(0);
    BOOST_CHECK(empty_prefix.empty());
    BOOST_CHECK_EQUAL(empty_prefix.size(), 0);

    const view_type empty_suffix = view.subview(view.size());
    BOOST_CHECK(empty_suffix.empty());
    BOOST_CHECK_EQUAL(empty_suffix.size(), 0);
}

BOOST_AUTO_TEST_CASE(constructs_from_an_existing_span) {
    value_type coefficients[] = {5, 6, 7};
    const std::span<value_type, 3> span(coefficients);
    math::coefficient_view view(span);

    BOOST_CHECK_EQUAL(view.data(), coefficients);
    BOOST_CHECK_EQUAL(view.size(), 3);
    BOOST_CHECK_EQUAL(view[1], value_type(6));
}

BOOST_AUTO_TEST_SUITE_END()
