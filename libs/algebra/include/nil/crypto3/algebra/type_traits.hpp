//---------------------------------------------------------------------------//
// Copyright (c) 2018-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
// Copyright (c) 2024 Vasiliy Olekhov <vasiliy.olekhov@nil.foundation>
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

#ifndef CRYPTO3_ALGEBRA_TYPE_TRAITS_HPP
#define CRYPTO3_ALGEBRA_TYPE_TRAITS_HPP

#include <concepts>
#include <cstddef>
#include <type_traits>

#include <boost/multiprecision/cpp_int.hpp>

namespace nil {
    namespace crypto3 {
        namespace algebra {

            template<typename T>
            concept FieldValue = requires(const T &a, const T &b, const boost::multiprecision::cpp_int &exponent) {
                typename T::field_type;
                { T::zero() } -> std::convertible_to<T>;
                { T::one() } -> std::convertible_to<T>;
                { a + b } -> std::same_as<T>;
                { a - b } -> std::same_as<T>;
                { a * b } -> std::same_as<T>;
                { a == b } -> std::convertible_to<bool>;
                { a.is_zero() } -> std::convertible_to<bool>;
                { a.is_one() } -> std::convertible_to<bool>;
                { a.inversed() } -> std::same_as<T>;
                { a.pow(exponent) } -> std::same_as<T>;
            };

            template<typename T>
            concept Field = requires {
                typename T::value_type;
                typename T::integral_type;
                typename T::modular_type;
                { T::value_bits } -> std::convertible_to<std::size_t>;
                { T::modulus_bits } -> std::convertible_to<std::size_t>;
                { T::arity } -> std::convertible_to<std::size_t>;
                requires FieldValue<typename T::value_type>;
            };

            template<typename T>
            concept ExtendedField = Field<T> && requires { typename T::extension_policy; };

            template<typename T>
            concept ExtendedFieldValue = FieldValue<T> && requires { typename T::underlying_type; };

            template<typename T>
            concept FieldElementWithCoordinates =
                FieldValue<T> && requires(T &value, const T &const_value, std::size_t index) {
                    value.coordinate(index);
                    const_value.coordinate(index);
                    requires std::is_lvalue_reference_v<decltype(value.coordinate(index))>;
                    requires std::is_lvalue_reference_v<decltype(const_value.coordinate(index))>;
                    requires(!std::is_const_v<std::remove_reference_t<decltype(value.coordinate(index))>>);
                    requires std::is_const_v<std::remove_reference_t<decltype(const_value.coordinate(index))>>;
                    requires std::same_as<std::remove_cvref_t<decltype(value.coordinate(index))>,
                                          std::remove_cvref_t<decltype(const_value.coordinate(index))>>;
                    requires FieldValue<std::remove_cvref_t<decltype(value.coordinate(index))>>;
                };

            template<typename T>
            concept Curve = requires {
                typename T::base_field_type;
                typename T::scalar_field_type;
                typename T::template g1_type<>;
                requires Field<typename T::base_field_type>;
                requires Field<typename T::scalar_field_type>;
            };

            template<typename T>
            concept CurveWithG2 = Curve<T> && requires { typename T::template g2_type<>; };

            template<typename T>
            concept CurveWithTargetGroup = Curve<T> && requires { typename T::gt_type; };

            template<typename T>
            concept CurveElement = requires(T &value, const T &const_value) {
                typename T::field_type;
                typename T::group_type;
                { T::zero() } -> std::convertible_to<T>;
                { T::one() } -> std::convertible_to<T>;
                { const_value.is_zero() } -> std::convertible_to<bool>;
                { const_value.is_well_formed() } -> std::convertible_to<bool>;
                { value.double_inplace() } -> std::same_as<void>;
            };

            template<typename T>
            concept CurveGroup = requires {
                typename T::params_type;
                typename T::curve_type;
                typename T::field_type;
                typename T::value_type;
                requires Curve<typename T::curve_type>;
                requires Field<typename T::field_type>;
                requires CurveElement<typename T::value_type>;
                requires std::same_as<typename T::value_type::group_type, T>;
            };

        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ALGEBRA_TYPE_TRAITS_HPP
