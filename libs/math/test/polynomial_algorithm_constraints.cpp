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

#include <list>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/math/polynomial/basic_operations.hpp>
#include <nil/crypto3/math/polynomial/evaluate.hpp>
#include <nil/crypto3/math/polynomial/lagrange_interpolation.hpp>

namespace math = nil::crypto3::math;
namespace fields = nil::crypto3::algebra::fields;

using value_type = fields::babybear::value_type;
using coefficients = std::vector<value_type>;
using interpolation_points = std::vector<std::pair<value_type, value_type>>;

template<typename Range>
concept SupportsPolynomialPower = requires(const Range &input) { math::power(input, std::size_t(2)); };

template<typename Range>
concept SupportsPolynomialMultiplication =
    requires(Range &output, const Range &input) { math::multiplication(output, input, input); };

template<typename Range, typename Point>
concept SupportsLagrangeEvaluation = requires(const Range &domain, const Point &point) {
    math::evaluate_lagrange_polynomial(domain, point, domain.size(), std::size_t(0));
};

template<typename Range>
concept SupportsLagrangeInterpolation = requires(const Range &points) { math::lagrange_interpolation(points); };

static_assert(SupportsPolynomialPower<coefficients>);
static_assert(!SupportsPolynomialPower<std::vector<int>>);

static_assert(SupportsPolynomialMultiplication<coefficients>);
static_assert(!SupportsPolynomialMultiplication<std::vector<int>>);

static_assert(SupportsLagrangeEvaluation<coefficients, value_type>);
static_assert(!SupportsLagrangeEvaluation<std::vector<int>, value_type>);
static_assert(!SupportsLagrangeEvaluation<std::list<value_type>, value_type>);

static_assert(SupportsLagrangeInterpolation<interpolation_points>);
static_assert(!SupportsLagrangeInterpolation<std::vector<value_type>>);
static_assert(!SupportsLagrangeInterpolation<std::list<std::pair<value_type, value_type>>>);
