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
#pragma once

#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace nil::crypto3::math {
    // An expression is deliberately weaker than a backend: uBLAS arithmetic
    // returns non-owning expression-template objects that are readable but not
    // necessarily default constructible or mutable.
    template<typename T>
    concept VectorExpression = requires(const std::remove_cvref_t<T> &value,
                                        typename std::remove_cvref_t<T>::size_type i) {
        typename std::remove_cvref_t<T>::value_type;
        typename std::remove_cvref_t<T>::size_type;
        { value.size() } -> std::convertible_to<typename std::remove_cvref_t<T>::size_type>;
        { value(i) } -> std::convertible_to<typename std::remove_cvref_t<T>::value_type>;
    };

    template<typename T>
    concept MatrixExpression = requires(const std::remove_cvref_t<T> &value,
                                        typename std::remove_cvref_t<T>::size_type i) {
        typename std::remove_cvref_t<T>::value_type;
        typename std::remove_cvref_t<T>::size_type;
        { value(i, i) } -> std::convertible_to<typename std::remove_cvref_t<T>::value_type>;
    } && (requires(const std::remove_cvref_t<T> &value) {
        { value.rows() } -> std::convertible_to<typename std::remove_cvref_t<T>::size_type>;
        { value.columns() } -> std::convertible_to<typename std::remove_cvref_t<T>::size_type>;
    } || requires(const std::remove_cvref_t<T> &value) {
        { value.size1() } -> std::convertible_to<typename std::remove_cvref_t<T>::size_type>;
        { value.size2() } -> std::convertible_to<typename std::remove_cvref_t<T>::size_type>;
    });

    template<MatrixExpression Expression>
    auto rows(const Expression &value) {
        if constexpr (requires { value.rows(); }) {
            return value.rows();
        } else {
            return value.size1();
        }
    }

    template<MatrixExpression Expression>
    auto columns(const Expression &value) {
        if constexpr (requires { value.columns(); }) {
            return value.columns();
        } else {
            return value.size2();
        }
    }

    template<typename T>
    concept VectorBackend = VectorExpression<T> && std::semiregular<T> &&
        requires(T &value, typename T::size_type i, typename T::value_type element) {
            value(i) = element;
        };

    template<typename T>
    concept MatrixBackend = MatrixExpression<T> && std::semiregular<T> &&
        requires(T &value, typename T::size_type i, typename T::value_type element) {
            value(i, i) = element;
        };

    template<typename T>
    concept ResizableVectorBackend = VectorBackend<T> && requires(T &value, typename T::size_type i) {
        value.resize(i);
    };

    template<typename T>
    concept ResizableMatrixBackend = MatrixBackend<T> && requires(T &value, typename T::size_type i) {
        value.resize(i, i);
    };
}    // namespace nil::crypto3::math
