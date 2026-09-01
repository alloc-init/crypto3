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

#include <algorithm>

#include <nil/crypto3/math/static_matrix/static_vector.hpp>
#include <nil/crypto3/math/static_matrix/static_matrix.hpp>
#include <nil/crypto3/math/static_matrix/utility.hpp>

namespace nil::crypto3::math {
    /** \addtogroup static_vector
     *  @{
     */

    /** @brief computes the elementwise square root
     *  @param v an N-static_vector of type T
     *  @return an N-static_vector \f$ \begin{bmatrix} \sqrt{v_1} & \ldots &\sqrt{v_N} \end{bmatrix} \f$ of type T
     *
     *  Computes the elementwise square root of a static_vector.
     */
    template<typename T, std::size_t N>
    constexpr static_vector<T, N> sqrt(const static_vector<T, N> &v) {
        return elementwise(static_cast<T (*)(T)>(sqrt), v);
    }

    /** @brief computes the dot product
     *  @param a an N-static_vector of type T
     *  @param b an N-static_vector of type T
     *  @return a scalar \f$ \textbf{a} \cdot \textbf{b} \f$ of type T such that
     *  \f$ \left(\textbf{a}\cdot\textbf{b}\right)_i = a_i \overline{b_i} \f$
     *
     *  Computes the dot (inner) product of two vectors.
     */
    template<typename T, std::size_t N>
    constexpr T dot(const static_vector<T, N> &a, const static_vector<T, N> &b) {
        T r = 0;
        for (std::size_t i = 0; i < static_vector<T, N>::size; ++i)
            r += a[i] * conj(b[i]);
        return r;
    }

    /** @brief computes the sum of elements
     *  @param v an N-static_vector of type T
     *  @return a scalar \f$ \sum\limits_{i} v_i \f$ of type T
     *
     *  Computes the sum of the elements of a static_vector.
     */
    template<typename T, std::size_t N>
    constexpr T sum(const static_vector<T, N> &v) {
        return accumulate(v, T(0u), std::plus<T>());
    }

    /** @brief computes the minimum valued element
     *  @param v an N-static_vector of type T
     *  @return a scalar \f$ v_i \f$ of type T where \f$ v_i \leq v_j,\ \forall j \f$
     *
     *  Computes the minimum valued element of a static_vector.
     */
    template<typename T, std::size_t N>
    constexpr T min(const static_vector<T, N> &v) {
        return accumulate(v, v[0], [](T a, T b) { return std::min(a, b); });
    }

    /** @brief computes the maximum valued element
     *  @param v an N-static_vector of type T
     *  @return a scalar \f$ v_i \f$ of type T where \f$ v_i \geq v_j,\ \forall j \f$
     *
     *  Computes the maximum valued element of a static_vector.
     */
    template<typename T, std::size_t N>
    constexpr T max(const static_vector<T, N> &v) {
        return accumulate(v, v[0], [](T a, T b) { return std::max(a, b); });
    }

    /** @brief computes the index of the minimum valued element
     *  @param v an N-static_vector of type T
     *  @return an index \f$ i \f$ where \f$ v_i \leq v_j,\ \forall j \f$
     *
     *  Computes the index of the minimum valued element of a static_vector.
     *  Note: the return value is zero-indexed.
     */
    template<typename T, std::size_t N>
    constexpr std::size_t min_index(const static_vector<T, N> &v) {
        T min = v[0];
        std::size_t index = 0;
        for (std::size_t i = 0; i < static_vector<T, N>::size; ++i)
            if (v[i] < min) {
                index = i;
                min = v[i];
            }
        return index;
    }

    /** @brief computes the index of the maximum valued element
     *  @param v an N-static_vector of type T
     *  @return an index \f$ i \f$ where \f$ v_i \geq v_j,\ \forall j \f$
     *
     *  Computes the index of the maximum valued element of a static_vector.
     *  Note: the return value is zero-indexed.
     */
    template<typename T, std::size_t N>
    constexpr std::size_t max_index(const static_vector<T, N> &v) {
        T max = v[0];
        std::size_t index = 0;
        for (std::size_t i = 0; i < static_vector<T, N>::size; ++i)
            if (v[i] > max) {
                index = i;
                max = v[i];
            }
        return index;
    }

    /** \addtogroup static_matrix
     *  @{
     */

    /** @brief computes the transpose
     *  @param m an \f$ M \times N \f$ static_matrix of type T
     *  @return an \f$ N \times M \f$ static_matrix \f$ \textbf{m}^{\mathrm{T}} \f$ of type T such that
     *  \f$ \left(\textbf{m}^{\mathrm{T}}\right)_{ij} = \textbf{m}_{ji},\ \forall i,j \f$
     *
     *  Computes the static_matrix transpose.
     */
    template<typename T, std::size_t M, std::size_t N>
    constexpr static_matrix<T, N, M> transpose(const static_matrix<T, M, N> &m) {
        return generate<N, M>([&m](auto i, auto j) { return m[j][i]; });
    }

    /** @brief computes the static_matrix product
     *  @param a an \f$M \times N\f$ static_matrix
     *  @param b an \f$N \times P\f$ static_matrix
     *  @return an \f$ M \times P \f$ static_matrix \f$ \textbf{a}\textbf{b} \f$ of type T such that
     *  \f$ \left(\textbf{ab}\right)_{ij} = \sum\limits_{k=1}^{N}\textbf{a}_{ik}\textbf{b}_{kj} \f$
     *
     *  Computes the product of two matrices.
     */
    template<typename T, std::size_t M, std::size_t N, std::size_t P>
    constexpr static_matrix<T, M, P> matmul(const static_matrix<T, M, N> &a, const static_matrix<T, N, P> &b) {
        return generate<M, P>([&a, &b](auto i, auto j) { return sum(a.row(i) * b.column(j)); });
    }

    /*!
     * @brief computes the product of static_vector and static_matrix
     * @param v an M-static_vector
     * @param m an \f$M \times N\f$ static_matrix
     * @return an N-static_vector of type T
     */
    template<typename T, std::size_t M, std::size_t N>
    constexpr static_vector<T, N> vectmatmul(const static_vector<T, M> &v, const static_matrix<T, M, N> &m) {
        return generate<N>([&v, &m](auto i) { return sum(v * m.column(i)); });
    }

    /*!
     * @brief computes the product of static_matrix and static_vector
     * @param m an \f$M \times N\f$ static_matrix
     * @param v an N-static_vector
     * @return an M-static_vector of type T
     */
    template<typename T, std::size_t M, std::size_t N>
    constexpr static_vector<T, M> matvectmul(const static_matrix<T, M, N> &m, const static_vector<T, N> &v) {
        return generate<M>([&v, &m](auto i) { return sum(m.row(i) * v); });
    }

    /** @brief Computes the kronecker tensor product
     *  @param a an \f$M \times N\f$ static_matrix
     *  @param b an \f$P \times Q\f$ static_matrix
     *  @return An \f$ MP \times NQ \f$ static_matrix \f$ \textbf{a}\otimes\textbf{b} \f$ of type T such that
     *  \f$ \left(\textbf{a}\otimes\textbf{b}\right)_{ij} = \textbf{a}_{\lfloor i/P \rfloor,\lfloor j/Q
     * \rfloor}\textbf{b}_{i\textrm{%}P,j\textrm{%}Q} \f$ where \f$ i \textrm{%} P \f$ is the remainder of \f$
     * i/P \f$
     *
     * Computes the kronecker tensor product of two matrices.
     */
    template<typename T, std::size_t M, std::size_t N, std::size_t P, std::size_t Q>
    constexpr static_matrix<T, M * P, N * Q> kron(const static_matrix<T, M, N> &a, const static_matrix<T, P, Q> &b) {
        return generate<M * P, N * Q>([&a, &b](auto i, auto j) { return a[i / P][j / Q] * b[i % P][j % Q]; });
    }

    /** @brief Computes the maximum absolute column sum norm
     *  @param m an \f$M \times N\f$ static_matrix
     *  @return a scalar \f$ {\left\lVert \textbf{m} \right\rVert}_1 \f$ of type T
     * such that \f$ {\left\lVert \textbf{m} \right\rVert}_1 = \max\limits_j
     * \sum\limits_{i=1}^M \left\lvert \textbf{m}_{ij} \right\rvert \f$
     *
     *  Computes the maximum absolute column sum norm of a static_matrix.
     */
    template<typename T, std::size_t M, std::size_t N>
    constexpr T macs(const static_matrix<T, M, N> &m) {
        return max(generate<N>([&m](std::size_t i) { return sum(abs(m.column(i))); }));
    }

    /** @brief Computes the maximum absolute row sum norm
     *  @param m an \f$M \times N\f$ static_matrix
     *  @return a scalar \f$ {\left\lVert \textbf{m} \right\rVert}_\infty \f$ of
     * type T such that \f$ {\left\lVert \textbf{m} \right\rVert}_\infty = \max\limits_i
     * \sum\limits_{j=1}^N \left\lvert \textbf{m}_{ij} \right\rvert \f$
     *
     *  Computes the maximum absolute row sum norm of a static_matrix.
     */
    template<typename T, std::size_t M, std::size_t N>
    constexpr T mars(const static_matrix<T, M, N> &m) {
        return max(generate<M>([&m](std::size_t i) { return sum(abs(m.row(i))); }));
    }

    /// @private
    template<typename T, std::size_t M, std::size_t N>
    constexpr std::tuple<static_matrix<T, M, N>, std::size_t, T> gauss_jordan_impl(static_matrix<T, M, N> m) {

        auto negligible = [](const T &v) { return v == T::zero(); };

        T det = 1;
        std::size_t rank = 0;
        std::size_t i = 0, j = 0;
        while (i < M && j < N) {
            // Choose largest magnitude as pivot to avoid adding different magnitudes
            for (std::size_t ip = i + 1; ip < M; ++ip) {
                if (m[ip][j] > m[i][j]) {
                    for (std::size_t jp = 0; jp < N; ++jp) {
                        auto tmp = m[ip][jp];
                        m[ip][jp] = m[i][jp];
                        m[i][jp] = tmp;
                    }
                    det = -det;
                    break;
                }
            }

            // If m_ij is still 0, continue to the next column
            if (!negligible(m[i][j])) {
                // Scale m_ij to 1
                auto s = m[i][j];
                for (std::size_t jp = 0; jp < N; ++jp) {
                    m[i][jp] /= s;
                }
                det /= s;

                // Eliminate other values in the column
                for (std::size_t ip = 0; ip < M; ++ip) {
                    if (ip == i) {
                        continue;
                    }
                    if (!negligible(m[ip][j])) {
                        auto s = m[ip][j];
                        [&]() {    // wrap this in a lambda to get around a gcc bug
                            for (std::size_t jp = 0; jp < N; ++jp) {
                                m[ip][jp] -= s * m[i][jp];
                            }
                        }();
                    }
                }

                // Increment rank
                ++rank;

                // Select next row
                ++i;
            }
            ++j;
        }
        det = (rank == M) ? det : 0;
        return {m, rank, det};
    }

    /** @brief Compute the reduced row echelon form
     *  @param m an \f$ M \times N \f$ static_matrix of type T
     *  @return an \f$ M \times N \f$ static_matrix of type T, the reduced row echelon form
     * of \f$ \textbf{m} \f$
     *
     *  Computes the reduced row echelon form of a static_matrix using Gauss-Jordan
     * elimination.  The tolerance for determining negligible elements is \f$
     * \max\left(N, M\right) \cdot \epsilon \cdot {\left\lVert \textbf{m}
     * \right\rVert}_\infty \f$.
     */
    template<typename T, std::size_t M, std::size_t N>
    constexpr static_matrix<T, M, N> rref(const static_matrix<T, M, N> &m) {
        return std::get<0>(gauss_jordan_impl(m));
    }

    /** @brief Compute the reduced row echelon form
     *  @param m an \f$ M \times N \f$ static_matrix of type T
     *  @param tolerance the tolerance used to determine when an element is
     * negligible (near zero)
     *  @return an \f$ M \times N \f$ static_matrix of type T, the reduced row echelon form
     * of \f$ \textbf{m} \f$
     *
     *  Computes the reduced row echelon form of a static_matrix using Gauss-Jordan
     * elimination.
     */
    template<typename T, std::size_t M, std::size_t N>
    constexpr static_matrix<T, M, N> rref(const static_matrix<T, M, N> &m, T tolerance) {
        return std::get<0>(gauss_jordan_impl(m), tolerance);
    }

    /** @brief Compute the rank
     *  @param m \f$ M \times N \f$ static_matrix of type T
     *  @return a scalar \f$ \textrm{rank}\left(\textbf{m}\right) \f$
     *
     *  Computes the rank using the reduced row echelon form.
     */
    template<typename T, std::size_t M, std::size_t N>
    constexpr std::size_t rank(const static_matrix<T, M, N> &m) {
        return std::get<1>(gauss_jordan_impl(m));
    }

    /** @brief Compute the determinant
     *  @param m \f$ M \times M \f$ static_matrix of type T
     *  @return a scalar \f$ \left\lvert \textbf{m} \right\rvert \f$ of type T
     *
     *  Computes the determinant using the reduced row echelon form.
     */
    template<typename T, std::size_t M>
    constexpr T det(const static_matrix<T, M, M> &m) {
        return std::get<2>(gauss_jordan_impl(m));
    }

    /** @brief computes the static_matrix inverse
     *  @param m an \f$ M \times M \f$ static_matrix of type T
     *  @return The inverse of \f$ \textbf{m} \f$, \f$ \textbf{m}^{-1}\f$ such that
     *  \f$ \textbf{m}\textbf{m}^{-1} = \textbf{m}^{-1}\textbf{m} = \textbf{I}_{M}
     * \f$
     *
     *  Computes the inverse of a static_matrix using the reduced row echelon form.
     */
    template<typename T, std::size_t M>
    constexpr static_matrix<T, M, M> inverse(const static_matrix<T, M, M> &m) {
        if (rank(m) < M)
            throw "static_matrix is not invertible";
        return submat<M, M>(rref(horzcat(m, get_identity<T, M>())), 0, M);
    }

    /** @brief computes the trace
     *  @param m an \f$ M \times M \f$ static_matrix of type T
     *  @return the trace of \f$ \textbf{m} \f$, \f$ \textrm{tr}\left(\textbf{m}\right) \f$
     *  such that \f$ \textrm{tr}\left(\textbf{m}\right) = \sum\limits_{n=1}^{M} \textbf{m}_{nn} \f$
     *
     *  Computes the trace of a static_matrix.
     */
    template<typename T, std::size_t M>
    constexpr T trace(const static_matrix<T, M, M> &m) {
        return sum(generate<M>([&m](std::size_t i) { return m[i][i]; }));
    }

    /** }@*/
}    // namespace nil::crypto3::math
