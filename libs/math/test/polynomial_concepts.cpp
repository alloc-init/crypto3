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

#include <ostream>
#include <vector>

#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/math/polynomial/concepts.hpp>
#include <nil/crypto3/math/polynomial/types/polymorphic_polynomial.hpp>
#include <nil/crypto3/math/polynomial/types/polymorphic_polynomial_dfs.hpp>

namespace math = nil::crypto3::math;
namespace fields = nil::crypto3::algebra::fields;

using value_type = fields::babybear::value_type;
using coefficient_polynomial = math::polynomial<value_type>;
using evaluation_polynomial = math::polynomial_dfs<value_type>;
using polymorphic_coefficient_polynomial = math::polymorphic_polynomial<fields::babybear_fp4>;
using polymorphic_evaluation_polynomial = math::polymorphic_polynomial_dfs<fields::babybear_fp4>;

struct non_field_type;
struct small_non_field_type;

struct non_field_value {
    using field_type = non_field_type;
};
struct small_non_field_value {
    using field_type = small_non_field_type;
};

struct non_field_type {
    using value_type = non_field_value;
};

struct small_non_field_type {
    using value_type = small_non_field_value;
};

struct non_field {
    using value_type = non_field_value;

    struct small_subfield {
        using value_type = small_non_field_value;
    };
};

template<typename ValueType>
concept SupportsCoefficientPolynomialScalarOperators =
    requires(std::ostream &stream, const math::polynomial<ValueType> &polynomial, const ValueType &scalar) {
        { math::operator+(polynomial, scalar) } -> std::same_as<math::polynomial<ValueType>>;
        { math::operator+(scalar, polynomial) } -> std::same_as<math::polynomial<ValueType>>;
        { math::operator-(polynomial, scalar) } -> std::same_as<math::polynomial<ValueType>>;
        { math::operator-(scalar, polynomial) } -> std::same_as<math::polynomial<ValueType>>;
        { math::operator*(polynomial, scalar) } -> std::same_as<math::polynomial<ValueType>>;
        { math::operator*(scalar, polynomial) } -> std::same_as<math::polynomial<ValueType>>;
        { math::operator/(polynomial, scalar) } -> std::same_as<math::polynomial<ValueType>>;
        { math::operator/(scalar, polynomial) } -> std::same_as<math::polynomial<ValueType>>;
        { math::operator<<(stream, polynomial) } -> std::same_as<std::ostream &>;
    };

template<typename ValueType>
concept SupportsEvaluationPolynomialScalarOperators =
    requires(std::ostream &stream, const math::polynomial_dfs<ValueType> &polynomial, const ValueType &scalar) {
        { math::operator+(polynomial, scalar) } -> std::same_as<math::polynomial_dfs<ValueType>>;
        { math::operator+(scalar, polynomial) } -> std::same_as<math::polynomial_dfs<ValueType>>;
        { math::operator-(polynomial, scalar) } -> std::same_as<math::polynomial_dfs<ValueType>>;
        { math::operator-(scalar, polynomial) } -> std::same_as<math::polynomial_dfs<ValueType>>;
        { math::operator*(polynomial, scalar) } -> std::same_as<math::polynomial_dfs<ValueType>>;
        { math::operator*(scalar, polynomial) } -> std::same_as<math::polynomial_dfs<ValueType>>;
        { math::operator/(polynomial, scalar) } -> std::same_as<math::polynomial_dfs<ValueType>>;
        { math::operator/(scalar, polynomial) } -> std::same_as<math::polynomial_dfs<ValueType>>;
        { math::operator<<(stream, polynomial) } -> std::same_as<std::ostream &>;
    };

template<typename FieldType>
concept SupportsPolymorphicPolynomialStreams =
    requires(std::ostream &stream, const math::polymorphic_polynomial<FieldType> &coefficient_polynomial,
             const math::polymorphic_polynomial_dfs<FieldType> &evaluation_polynomial) {
        { math::operator<<(stream, coefficient_polynomial) } -> std::same_as<std::ostream &>;
        { math::operator<<(stream, evaluation_polynomial) } -> std::same_as<std::ostream &>;
    };

static_assert(math::CoefficientPolynomial<coefficient_polynomial>);
static_assert(math::CoefficientPolynomial<const coefficient_polynomial &>);
static_assert(!math::EvaluationPolynomial<coefficient_polynomial>);

static_assert(math::EvaluationPolynomial<evaluation_polynomial>);
static_assert(math::EvaluationPolynomial<const evaluation_polynomial &>);
static_assert(!math::CoefficientPolynomial<evaluation_polynomial>);

static_assert(math::CoefficientPolynomial<polymorphic_coefficient_polynomial>);
static_assert(!math::EvaluationPolynomial<polymorphic_coefficient_polynomial>);

static_assert(math::EvaluationPolynomial<polymorphic_evaluation_polynomial>);
static_assert(!math::CoefficientPolynomial<polymorphic_evaluation_polynomial>);

static_assert(!math::CoefficientPolynomial<std::vector<value_type>>);
static_assert(!math::EvaluationPolynomial<std::vector<value_type>>);

static_assert(SupportsCoefficientPolynomialScalarOperators<value_type>);
static_assert(!SupportsCoefficientPolynomialScalarOperators<non_field_value>);
static_assert(SupportsEvaluationPolynomialScalarOperators<value_type>);
static_assert(!SupportsEvaluationPolynomialScalarOperators<non_field_value>);
static_assert(SupportsPolymorphicPolynomialStreams<fields::babybear_fp4>);
static_assert(!SupportsPolymorphicPolynomialStreams<non_field>);

void instantiate_polymorphic_stream_operators(std::ostream &stream,
                                              const polymorphic_coefficient_polynomial &coefficient_polynomial,
                                              const polymorphic_evaluation_polynomial &evaluation_polynomial) {
    math::operator<<(stream, coefficient_polynomial);
    math::operator<<(stream, evaluation_polynomial);
}

void instantiate_scalar_operators(std::ostream &stream, const coefficient_polynomial &coefficient_polynomial,
                                  const evaluation_polynomial &evaluation_polynomial, const value_type &scalar) {
    math::operator+(coefficient_polynomial, scalar);
    math::operator+(scalar, coefficient_polynomial);
    math::operator-(coefficient_polynomial, scalar);
    math::operator-(scalar, coefficient_polynomial);
    math::operator*(coefficient_polynomial, scalar);
    math::operator*(scalar, coefficient_polynomial);
    math::operator/(coefficient_polynomial, scalar);
    math::operator/(scalar, coefficient_polynomial);
    math::operator<<(stream, coefficient_polynomial);

    math::operator+(evaluation_polynomial, scalar);
    math::operator+(scalar, evaluation_polynomial);
    math::operator-(evaluation_polynomial, scalar);
    math::operator-(scalar, evaluation_polynomial);
    math::operator*(evaluation_polynomial, scalar);
    math::operator*(scalar, evaluation_polynomial);
    math::operator/(evaluation_polynomial, scalar);
    math::operator<<(stream, evaluation_polynomial);
}
