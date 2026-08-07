//---------------------------------------------------------------------------//
// Copyright (c) 2020-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
// Copyright (c) 2020-2021 Ilias Khairullin <ilias@nil.foundation>
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

#pragma once

#include <nil/crypto3/math/static_matrix/matrix.hpp>
#include <nil/crypto3/math/static_matrix/vector.hpp>

namespace nil::crypto3::math {

    /** \addtogroup vector
     *  @{
     */

    /** @brief checks equality of two vectors
     *  @param a an N-vector of type T
     *  @param b an N-vector of type T
     *  @return true if and only if \f$ \textbf{a}_i = \textbf{b}_i\ \forall i \in 1\ .. N \f$
     *
     *  Checks the equality of two vectors.
     */
    template<typename T, std::size_t N>
    constexpr bool operator==(const vector<T, N> &a, const vector<T, N> &b) {
        for (std::size_t i = 0; i < vector<T, N>::size; ++i) {
            if (a[i] != b[i])
                return false;
        }
        return true;
    }

    /** @brief checks inequality of two vectors
     *  @param a an N-vector of type T
     *  @param b an N-vector of type T
     *  @return false if and only if \f$ \textbf{a}_i = \textbf{b}_i\ \forall i \in 1\ .. N \f$
     *
     *  Checks the inequality of two vectors.
     */
    template<typename T, std::size_t N>
    constexpr bool operator!=(const vector<T, N> &a, const vector<T, N> &b) {
        return !(a == b);
    }

    /** @brief computes the sum of a vector and a scalar
     *  @param v an N-vector of type T
     *  @param a a scalar of type T
     *  @return \f$ \textbf{v} + a \f$ such that \f$ \left(\textbf{v} + a\right)_i = \textbf{v}_i + a \f$
     *
     *  Computes the sum of a vector and a scalar.
     */
    template<typename T, std::size_t N>
    constexpr vector<T, N> operator+(const vector<T, N> &v, T a) {
        return elementwise([a](T x) { return x + a; }, v);
    }

    /** @brief computes the sum of a vector and a scalar
     *  @param a a scalar of type T
     *  @param v an N-vector of type T
     *  @return \f$ a + \textbf{v} \f$ such that \f$ \left(a + \textbf{v}\right)_i = a + \textbf{v}_i \f$
     *
     *  Computes the sum of a vector and a scalar.
     */
    template<typename T, std::size_t N>
    constexpr vector<T, N> operator+(T a, const vector<T, N> &v) {
        return elementwise([a](T x) { return a + x; }, v);
    }

    /** @brief computes the vector sum
     *  @param a an N-vector of type T
     *  @param b an N-vector of type T
     *  @return \f$ \textbf{a} + \textbf{b} \f$ such that \f$ \left(\textbf{a} + \textbf{b}\right)_i =
     * \textbf{a}_i
     * + \textbf{b}_i \f$
     *
     *  Computes the vector sum.
     */
    template<typename T, std::size_t N>
    constexpr vector<T, N> operator+(const vector<T, N> &a, const vector<T, N> &b) {
        return elementwise(std::plus<T>(), a, b);
    }

    /** @brief computes the product of a vector and a scalar
     *  @param v an N-vector of type T
     *  @param a a scalar of type T
     *  @return \f$ \textbf{v}a \f$ such that \f$ \left(\textbf{v}a\right)_i = \textbf{v}_i a \f$
     *
     *  Computes the product of a vector and a scalar.
     */
    template<typename T, std::size_t N>
    constexpr vector<T, N> operator*(const vector<T, N> &v, T a) {
        return elementwise([a](T x) { return x * a; }, v);
    }

    /** @brief computes the product of a vector and a scalar
     *  @param a a scalar of type T
     *  @param v an N-vector of type T
     *  @return \f$ a\textbf{v} \f$ such that \f$ \left(a\textbf{v}\right)_i = a\textbf{v}_i \f$
     *
     *  Computes the product of a vector and a scalar.
     */
    template<typename T, std::size_t N>
    constexpr vector<T, N> operator*(T a, const vector<T, N> &v) {
        return elementwise([a](T x) { return a * x; }, v);
    }

    /** @brief computes the Hadamard product
     *  @param a an N-vector of type T
     *  @param b an N-vector of type T
     *  @return \f$ \textbf{a} \circ \textbf{b} \f$ such that \f$ \left(\textbf{a} \circ \textbf{b}\right)_i =
     * \textbf{a}_i \textbf{b}_i \f$
     *
     *  Computes the Hadamard, or elementwise, product of two vectors.
     */
    template<typename T, std::size_t N>
    constexpr vector<T, N> operator*(const vector<T, N> &a, const vector<T, N> &b) {
        return elementwise(std::multiplies<T>(), a, b);
    }

    /** @brief computes the quotient between a vector and a scalar
     *  @param v an N-vector of type T
     *  @param a a scalar of type T
     *  @return \f$ \textbf{v}/a \f$ such that \f$ \left(\textbf{v}/a\right)_i = \frac{\textbf{v}_i}{a} \f$
     *
     *  Computes division between a vector and a scalar.
     */
    template<typename T, std::size_t N>
    constexpr vector<T, N> operator/(T a, const vector<T, N> &v) {
        return elementwise([a](T x) { return a / x; }, v);
    }

    /** @brief computes the elementwise vector quotient
     *  @param a an N-vector of type T
     *  @param b an N-vector of type T
     *  @return \f$ \textbf{a} \circ \textbf{b}' \f$ such that \f$ {\textbf{b}_i}' =
     * \left(\textbf{b}_i\right)^{-1}\f$ and \f$ \left(\textbf{a} \circ \textbf{b}'\right)_i = \textbf{a}_i
     * {\textbf{b}'}_i \f$
     *
     *  Computes elementwise division between two vectors
     */
    template<typename T, std::size_t N>
    constexpr vector<T, N> operator/(const vector<T, N> &a, const vector<T, N> &b) {
        return elementwise(std::divides<T>(), a, b);
    }

    /** }@*/

    /** \addtogroup matrix
     *  @{
     */

    /** @brief checks equality of two matrices
     *  @param a an \f$ N \times M \f$ matrix of type T
     *  @param b an \f$ N \times M \f$ matrix of type T
     *  @return true if and only if \f$ \textbf{a}_{ij} = \textbf{b}_{ij}\ \forall i,j \in 1\ .. N \f$
     *
     *  Checks the equality of two matrices.
     */
    template<typename T, std::size_t N, std::size_t M>
    constexpr bool operator==(const matrix<T, N, M> &a, const matrix<T, N, M> &b) {
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < M; ++j) {
                if (a[i][j] != b[i][j])
                    return false;
            }
        }
        return true;
    }

    /** @brief checks inequality of two matrices
     *  @param a an \f$ N \times M \f$ matrix of type T
     *  @param b an \f$ N \times M \f$ matrix of type T
     *  @return false if and only if \f$ \textbf{a}_{ij} = \textbf{b}_{ij}\ \forall i,j \in 1\ .. N \f$
     *
     *  Checks the inequality of two matrices.
     */
    template<typename T, std::size_t N, std::size_t M>
    constexpr bool operator!=(const matrix<T, N, M> &a, const matrix<T, N, M> &b) {
        return !(a == b);
    }

    /** @brief computes the sum of a matrix and a scalar
     *  @param m an \f$ N \times M \f$ matrix of type T
     *  @param a a scalar of type T
     *  @return \f$ \textbf{m} + a \f$ such that \f$ \left(\textbf{m} + a\right)_{ij} = \textbf{m}_{ij} + a
     * \f$
     *
     *  Computes the sum of a matrix and a scalar.
     */
    template<typename T, std::size_t N, std::size_t M>
    constexpr matrix<T, N, M> operator+(const matrix<T, N, M> &m, T a) {
        return elementwise([a](T x) { return x + a; }, m);
    }

    /** @brief computes the sum of a matrix and a scalar
     *  @param a a scalar of type T
     *  @param m an \f$ N \times M \f$ matrix of type T
     *  @return \f$ a + \textbf{m} \f$ such that \f$ \left(a + \textbf{m}\right)_{ij} = a + \textbf{m}_{ij}
     * \f$
     *
     *  Computes the sum of a matrix and a scalar.
     */
    template<typename T, std::size_t N, std::size_t M>
    constexpr matrix<T, N, M> operator+(T a, const matrix<T, N, M> &m) {
        return elementwise([a](T x) { return a + x; }, m);
    }

    /** @brief computes the matrix sum
     *  @param a an \f$ N \times M \f$ matrix of type T
     *  @param b an \f$ N \times M \f$ matrix of type T
     *  @return \f$ \textbf{a} + \textbf{b} \f$ such that \f$ \left(\textbf{a} + \textbf{b}\right)_{ij} =
     * \textbf{a}_{ij} + \textbf{b}_{ij} \f$
     *
     *  Computes the vector sum.
     */
    template<typename T, std::size_t N, std::size_t M>
    constexpr matrix<T, N, M> operator+(const matrix<T, N, M> &a, const matrix<T, N, M> &b) {
        return elementwise(std::plus<T>(), a, b);
    }

    /** @brief computes the product of a matrix and a scalar
     *  @param m an \f$ N \times M \f$ matrix of type T
     *  @param a a scalar of type T
     *  @return \f$ \textbf{m}a \f$ such that \f$ \left(\textbf{m} a\right)_{ij} = \textbf{m}_{ij} a \f$
     *
     *  Computes the sum of a matrix and a scalar.
     */
    template<typename T, std::size_t N, std::size_t M>
    constexpr matrix<T, N, M> operator*(const matrix<T, N, M> &m, T a) {
        return elementwise([a](T x) { return x * a; }, m);
    }

    /** @brief computes the product of a matrix and a scalar
     *  @param a a scalar of type T
     *  @param m an \f$ N \times M \f$ matrix of type T
     *  @return \f$ a\textbf{m} \f$ such that \f$ \left(a\textbf{m}\right)_{ij} = a\textbf{m}_{ij} \f$
     *
     *  Computes the sum of a matrix and a scalar.
     */
    template<typename T, std::size_t N, std::size_t M>
    constexpr matrix<T, N, M> operator*(T a, const matrix<T, N, M> &m) {
        return elementwise([a](T x) { return a * x; }, m);
    }

    /** @brief computes the Hadamard product
     *  @param a an \f$ N \times M \f$ matrix of type T
     *  @param b an \f$ N \times M \f$ matrix of type T
     *  @return \f$ \textbf{a} \circ \textbf{b} \f$ such that \f$ \left(\textbf{a} \circ
     * \textbf{b}\right)_{ij} = \textbf{a}_{ij} \textbf{b}_{ij} \f$
     *
     *  Computes the Hadamard, or elementwise, product of two vectors.
     */
    template<typename T, std::size_t N, std::size_t M>
    constexpr matrix<T, N, M> operator*(const matrix<T, N, M> &a, const matrix<T, N, M> &b) {
        return elementwise(std::multiplies<T>(), a, b);
    }

    /** @brief computes the quotient between a matrix and a scalar
     *  @param m an \f$ N \times M \f$ matrix of type T
     *  @param a a scalar of type T
     *  @return \f$ \textbf{m}/a \f$ such that \f$ \left(\textbf{m}/a\right)_{ij} =
     * \frac{\textbf{m}_{ij}}{a} \f$
     *
     *  Computes division between a matrix and a scalar.
     */
    template<typename T, std::size_t N, std::size_t M>
    constexpr matrix<T, N, M> operator/(T a, const matrix<T, N, M> &m) {
        return elementwise([a](T x) { return a / x; }, m);
    }

    /** @brief computes the elementwise matrix quotient
     *  @param a an \f$ N \times M \f$ matrix of type T
     *  @param b an \f$ N \times M \f$ matrix of type T
     *  @return \f$ \textbf{a} \circ \textbf{b}' \f$ such that \f$ {\textbf{b}_{ij}}' =
     * \left(\textbf{b}_{ij}\right)^{-1}\f$ and \f$ \left(\textbf{a} \circ \textbf{b}'\right)_{ij} =
     * \textbf{a}_{ij}
     * {\textbf{b}'}_{ij} \f$
     *
     *  Computes elementwise division between two matrices
     */
    template<typename T, std::size_t N, std::size_t M>
    constexpr matrix<T, N, M> operator/(const matrix<T, N, M> &a, const matrix<T, N, M> &b) {
        return elementwise(std::divides<T>(), a, b);
    }

    /** }@*/
}    // namespace nil::crypto3::math
