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

#include <array>
#include <list>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/math/polynomial/basic_operations.hpp>
#include <nil/crypto3/math/polynomial/basis_change.hpp>
#include <nil/crypto3/math/polynomial/evaluate.hpp>
#include <nil/crypto3/math/polynomial/lagrange_interpolation.hpp>
#include <nil/crypto3/math/polynomial/shift.hpp>
#include <nil/crypto3/math/polynomial/xgcd.hpp>

namespace math = nil::crypto3::math;
namespace fields = nil::crypto3::algebra::fields;

using value_type = fields::babybear::value_type;
using coefficients = std::vector<value_type>;
using coefficient_polynomial = math::polynomial<value_type>;
using evaluation_polynomial = math::polynomial_dfs<value_type>;
using interpolation_points = std::vector<std::pair<value_type, value_type>>;

template<typename Range>
concept SupportsPolynomialMultiplication =
    requires(Range &output, const Range &input) { math::multiplication(output, input, input); };

template<typename Range>
concept SupportsMutableCoefficientOperations = requires(Range &output, const Range &input) {
    math::reverse(output, std::size_t(1));
    math::condense(output);
    math::addition(output, input, input);
    math::subtraction(output, input, input);
};

template<typename Range>
concept SupportsPolynomialDivision = requires(Range &quotient, Range &remainder, const Range &input) {
    math::division(quotient, remainder, input, input);
};

template<typename Range>
concept SupportsExtendedEuclidean =
    requires(Range &g, Range &u, Range &v, const Range &input) { math::extended_euclidean(input, input, g, u, v); };

template<typename Range, typename Point>
concept SupportsPolynomialEvaluation = requires(const Range &coefficients, const Point &point) {
    math::evaluate_polynomial(coefficients, point, coefficients.size());
};

template<typename Range, typename Point>
concept SupportsLagrangeEvaluation = requires(const Range &domain, const Point &point) {
    math::evaluate_lagrange_polynomial(domain, point, domain.size(), std::size_t(0));
};

template<typename Range>
concept SupportsLagrangeInterpolation = requires(const Range &points) { math::lagrange_interpolation(points); };

template<typename Range>
concept SupportsGeometricBasisChange = requires(Range &values, const Range &field_values) {
    math::monomial_to_newton_basis_geometric<fields::babybear>(values, field_values, field_values, std::size_t(1));
    math::newton_to_monomial_basis_geometric<fields::babybear>(values, field_values, field_values, std::size_t(1));
};

template<typename FieldType>
concept SupportsSubproductTree = requires(std::vector<std::vector<std::vector<typename FieldType::value_type>>> &tree) {
    math::compute_subproduct_tree<FieldType>(tree, std::size_t(1));
};

template<typename ValueType>
concept SupportsPolynomialShift = requires(const math::polynomial<ValueType> &polynomial, const ValueType &value) {
    math::polynomial_shift(polynomial, value);
};

static_assert(math::detail::PolynomialCoefficientRange<coefficients>);
static_assert(math::detail::MutablePolynomialCoefficientRange<coefficients>);
static_assert(math::detail::MutableNormalizableCoefficientPolynomial<coefficient_polynomial>);
static_assert(!math::detail::MutableNormalizableCoefficientPolynomial<evaluation_polynomial>);
static_assert(!math::detail::MutableNormalizableCoefficientPolynomial<coefficients>);
static_assert(math::detail::PolynomialCoefficientRange<std::array<value_type, 2>>);
static_assert(!math::detail::MutablePolynomialCoefficientRange<std::array<value_type, 2>>);
static_assert(!math::detail::PolynomialCoefficientRange<std::list<value_type>>);

static_assert(SupportsPolynomialMultiplication<coefficients>);
static_assert(!SupportsPolynomialMultiplication<std::vector<int>>);

static_assert(SupportsMutableCoefficientOperations<coefficients>);
static_assert(SupportsMutableCoefficientOperations<std::vector<int>>);
static_assert(!SupportsMutableCoefficientOperations<std::list<value_type>>);

static_assert(SupportsPolynomialDivision<coefficients>);
static_assert(!SupportsPolynomialDivision<std::vector<int>>);

static_assert(SupportsExtendedEuclidean<coefficients>);
static_assert(!SupportsExtendedEuclidean<std::vector<int>>);

static_assert(SupportsPolynomialEvaluation<coefficients, value_type>);
static_assert(!SupportsPolynomialEvaluation<std::list<value_type>, value_type>);
static_assert(!SupportsPolynomialEvaluation<std::vector<int>, value_type>);

static_assert(SupportsLagrangeEvaluation<coefficients, value_type>);
static_assert(!SupportsLagrangeEvaluation<std::vector<int>, value_type>);
static_assert(!SupportsLagrangeEvaluation<std::list<value_type>, value_type>);

static_assert(SupportsLagrangeInterpolation<interpolation_points>);
static_assert(!SupportsLagrangeInterpolation<std::vector<value_type>>);
static_assert(!SupportsLagrangeInterpolation<std::list<std::pair<value_type, value_type>>>);

static_assert(SupportsGeometricBasisChange<coefficients>);
static_assert(!SupportsGeometricBasisChange<std::list<value_type>>);
static_assert(SupportsSubproductTree<fields::babybear>);

static_assert(SupportsPolynomialShift<value_type>);
static_assert(!SupportsPolynomialShift<int>);
