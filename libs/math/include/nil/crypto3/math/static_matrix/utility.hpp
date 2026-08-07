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
#include <nil/crypto3/math/detail/same_value.hpp>

#include <tuple>

namespace nil::crypto3::math {
    template<typename T, std::size_t N, std::size_t M>
    struct matrix;

    template<typename T, std::size_t N>
    struct vector;

    /** \addtogroup vector
     *  @{
     */

    /** @brief applies a function elementwise between many vectors
     *  @param f a function of type F that operates on many scalars of type T and returns a scalar of type U
     *  @param v an N-vector of type T
     *  @param vectors additional N-vectors of type T
     *  @return an N-vector of type T with elements described by \f$ f\left(\textbf{v}_i, \ldots\right) \f$
     *
     *  Applies a function elementwise between many vectors.
     */
    template<typename F, typename T, typename... Vectors,
             typename U = typename std::invoke_result<F, T, typename Vectors::value_type...>::type,
             std::size_t N = detail::all_same_value<std::size_t, Vectors::size...>::value>
    constexpr vector<U, N> elementwise(F f, const vector<T, N> &v, const Vectors &...vectors) {
        vector<U, N> op_applied = {};
        for (std::size_t i = 0; i < N; ++i) {
            op_applied[i] = std::apply(f, std::forward_as_tuple(v[i], vectors[i]...));
        }
        return op_applied;
    }

    /** @brief accumulates an operation across a vector
     *  @param v an N-vector of type T
     *  @param init the initial value
     *  @param f a function of type F that operates between U and vector elements of type T
     *  @return \f$ f\left(f\left(\ldots f\left(\textrm{init}, \textbf{v}_1\right), \ldots\right), \textbf{v}_N
     * \right) \f$
     *
     *  Accumulates an operation over the elements.  This is equivalent to a functional fold.
     */
    template<typename T, std::size_t N, typename F, typename U>
    constexpr U accumulate(const vector<T, N> &v, U init, F &&f) {
        U r = init;
        for (std::size_t i = 0; i < vector<T, N>::size; ++i)
            r = std::apply(std::forward<F>(f), std::forward_as_tuple(r, v[i]));
        return r;
    }

    /** @brief casts a vector to another type
     *  @param v an N-vector of type U
     *  @return an N-vector of type T containing the casted elements of \f$ \textbf{v} \f$
     *
     *  Casts a vector to another type by `static_cast`ing each element.
     */
    template<typename T, typename U, std::size_t N>
    constexpr vector<T, N> cast(const vector<U, N> &v) {
        return elementwise([](const U u) { return static_cast<T>(u); }, v);
    }

    /** @brief generates a vector containing consecutive elements
     *  @param value the value of the first element of the vector
     *  @return an N-vector of type T such that \f$ \textbf{v}_i = \textrm{start} + i - 1 \f$
     *
     *  Generates a vector containing consecutive elements spaced by 1.
     */
    template<std::size_t N, typename T>
    constexpr vector<T, N> iota(T value = T()) {
        vector<T, N> seq = {};
        for (auto &x : seq) {
            x = value;
            value += 1;    // equivalent to value++, see GCC Bug 91705
        }
        return seq;
    }

    /** @brief generates a vector of equally spaced elements
     *  @param min the value of the first element of the vector
     *  @param max the value of the last element of the vector
     *  @return an N-vector of type T with elements spaced by \f$ \frac{\textbf{v}_N - \textbf{v}_1}{N-1} \f$
     *
     *  Generates a vector of equally spaced elements.
     */
    template<std::size_t N, typename T>
    constexpr vector<T, N> linspace(T min, T max) {
        return ((max - min) / (N - 1)) * iota<N, T>() + min;
    }

    /** @brief generates a vector containing a single value
     *  @param value the scalar value of all elements
     *  @return an N-vector of type T such that \f$ \textbf{v}_i = \textrm{value}\ \forall i \in 1\ .. N\f$
     *
     *  Generates a vector with all elements equal to a single value.
     */
    template<std::size_t N, typename T>
    constexpr vector<T, N> fill(T value) {
        vector<T, N> filled = {};
        for (auto &x : filled)
            x = value;
        return filled;
    }

    /** @brief generates a vector as a function of its index
     *  @param f a function that operates on an integer index
     *  @return an N-vector with type matching the return type of f such that \f$ \textbf{v}_i = f(i) \f$
     *
     *  Generates a vector as a function of its index.
     */
    template<std::size_t N, typename F>
    constexpr decltype(auto) generate(F &&f) {
        return elementwise(f, iota<N, std::size_t>());
    }

    /** @brief shifts vector elements
     *  @param v an N-vector of type T
     *  @param n the amount to shift each element
     *  @return an N-vector of type T \f$ \textbf{v} \gg n \f$ such that
     *  \f$ \left(\textbf{v} \gg n\right)_i = \textbf{v}_{(i + n)\ \textrm{mod}\ N} \f$
     *
     *  Rotates a vector by shifting its elements.
     */
    template<std::size_t N, typename T>
    constexpr vector<T, N> rotate(vector<T, N> v, int n) {
        vector<T, N> rotated = {};
        // add N (the modulus) to n until it is positive
        while (n < 0) {
            n += N;
        }
        for (std::size_t i = 0; i < N; ++i) {
            rotated[i] = v[(i + n) % N];
        }
        return rotated;
    }

    /** @brief slices a vector into a subvector
     *  @param v an N-vector of type T
     *  @param start the first index of the subvector
     *  @return an M-vector \f$ \textbf{v}_{\textrm{start}:\left(\textrm{start} + M - 1\right)} \f$
     *  such that  \f$ \left(\textbf{v}_{\textrm{start}:\left(\textrm{start} + M - 1\right)}\right)_i =
     * \textbf{v}_{\textrm{start} + i} \f$
     *
     *  Slices a vector into a subvector.
     */
    template<std::size_t M, typename T, std::size_t N>
    constexpr vector<T, M> slice(vector<T, N> v, std::size_t start = 0) {

        vector<T, M> sliced = {};
        for (std::size_t i = 0; i < M; ++i)
            sliced[i] = v[i + start];

        return sliced;
    }

    /** @brief concatenates two vectors
     *  @param a an N-vector of type T
     *  @param b an M-vector of type T
     *  @return an N+M-vector \f$ \left[\textbf{a} \textbf{b}\right] \f$ such that
     *  \f$ \left(\left[\textbf{a} \textbf{b}\right]\right)_i = \begin{cases} \textbf{a}_i & i < N\\
     * \textbf{b}_{i - N} & i \ge N \end{cases} \f$
     *
     *  Slices a vector into a subvector.
     */
    template<typename T, std::size_t N, std::size_t M>
    constexpr vector<T, N + M> concat(vector<T, N> a, vector<T, M> b) {
        vector<T, N + M> concatted = {};
        for (std::size_t i = 0; i < N; ++i) {
            concatted[i] = a[i];
        }
        for (std::size_t i = 0; i < M; ++i) {
            concatted[i + N] = b[i];
        }
        return concatted;
    }

    /** }@*/

    /** \addtogroup matrix
     *  @{
     */

    /** @brief applies a function elementwise between many matrices
     *  @param f a function of type F that operates on many scalars of type T and returns a scalar of type U
     *  @param m an \f$ N \times M \f$ matrix of type T
     *  @param matrices additional \f$ N \times M \f$ matrices of type T
     *  @return an \f$ N \times M \f$ matrix of type T with elements described by \f$ f\left(\textbf{m}_{ij},
     * \ldots\right) \f$
     *
     *  Applies a function elementwise between many matrices.
     */
    template<typename F, typename T, typename... Matrices,
             typename U = std::invoke_result_t<F, T, typename Matrices::value_type...>,
             std::size_t N = detail::all_same_value<std::size_t, Matrices::column_size...>::value,
             std::size_t M = detail::all_same_value<std::size_t, Matrices::row_size...>::value>
    constexpr matrix<U, N, M> elementwise(F f, const matrix<T, N, M> &m, const Matrices &...matrices) {
        matrix<U, N, M> op_applied = {};
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < M; ++j) {
                op_applied[i][j] = std::apply(f, std::forward_as_tuple(m[i][j], matrices[i][j]...));
            }
        }
        return op_applied;
    }

    /** @brief casts a vector to another type
     *  @param m an \f$ N \times M \f$ matrix of type U
     *  @return an \f$ N \times M \f$ matrix of type T containing the casted elements of \f$ \textbf{M} \f$
     *
     *  Casts a matrix to another type by `static_cast`ing each element.
     */
    template<typename T, typename U, std::size_t N, std::size_t M>
    constexpr matrix<T, N, M> cast(const matrix<U, N, M> &m) {
        return elementwise([](const U u) { return static_cast<T>(u); }, m);
    }

    /** @brief generates a matrix as a function of its indices
     *  @param f a function that operates on two integer indices
     *  @return an \f$ N \times M \f$ matrix with type matching the return type of f such that \f$
     * \textbf{m}_{ij} = f(i, j) \f$
     *
     *  Generates a matrix as a function of its indices.
     */
    template<std::size_t N, std::size_t M, typename F>
    constexpr decltype(auto) generate(F &&f) {
        matrix<typename std::invoke_result<F, std::size_t, std::size_t>::type, N, M> generated = {};
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < M; ++j) {
                generated[i][j] = std::apply(f, std::forward_as_tuple(i, j));
            }
        }
        return generated;
    }

    /** @brief generates a matrix containing a single value
     *  @param value the scalar value of all elements
     *  @return an \f$ N \times M \f$ matrix of type T such that \f$ \textbf{m}_{ij} = \textrm{value}\ \forall
     * i,j \f$
     *
     *  Generates a matrix with all elements equal to a single value.
     */
    template<std::size_t N, std::size_t M, typename T>
    constexpr matrix<T, N, M> fill(T value) {
        return generate<N, M>([value](std::size_t, std::size_t) { return value; });
    }

    /** @brief the matrix identity
     *
     *  The matrix identity \f$ I_N \f$.
     */
    template<typename T, std::size_t N>
    constexpr matrix<T, N, N> identity = generate<N, N>([](std::size_t i, std::size_t j) { return T(i == j ? 1 : 0); });

    template<typename T, std::size_t N>
    constexpr matrix<T, N, N> get_identity() {
        return generate<N, N>([](std::size_t i, std::size_t j) { return T(i == j ? 1 : 0); });
    }

    /** @brief repeats a matrix
     *  @tparam Row the number of times to repeat in the row direction
     *  @tparam Col the number of times to repeat in the column direction
     *  @param m an \f$ M \times N \f$ matrix of type T
     *  @return an \f$ \left(M\times\textrm{Row}\right) \times \left(N\times\textrm{Col}\right) \f$
     *  matrix \f$ \textbf{m}' \f$ of type T such that
     * \f$ \textbf{m}'_{ij} = \textbf{m}_{\left(i\ \textrm{mod}\ M\right),\ \left(j\ \textrm{mod}\ N \right)}
     * \f$
     *
     *  Repeats copies of a matrix.
     */
    template<std::size_t Row, std::size_t Col, std::size_t M, std::size_t N, typename T>
    constexpr matrix<T, M * Row, N * Col> repmat(const matrix<T, M, N> &m) {
        return generate<M * Row, N * Col>([&m](std::size_t i, std::size_t j) { return m[i % M][j % N]; });
    }

    /** @brief swaps rows of a matrix
     *  @param m an \f$ M \times N \f$ matrix of type T
     *  @param a the index of a row to swap
     *  @param b the index of a row to swap
     *  @return an \f$ M \times N \f$ matrix \f$ \textbf{m}' \f$ of type T such that
     *  \f$ {\textbf{m}'}_{ij} = \begin{cases} \textbf{m}_{bj} & i = a\\ \textbf{m}_{aj} & i = b\\
     * \textbf{m}_{ij} & \textrm{otherwise} \end{cases} \f$
     *
     *  Swap two rows of a matrix.
     */
    template<std::size_t M, std::size_t N, typename T>
    constexpr matrix<T, M, N> swaprow(matrix<T, M, N> m, std::size_t a, std::size_t b) {
        for (int i = 0; i < N; i++) {
            T tmp = m[a][i];
            m[a][i] = m[b][i];
            m[b][i] = tmp;
        }
        return m;
    }

    /** @brief swaps columns of a matrix
     *  @param m an \f$ M \times N \f$ matrix of type T
     *  @param a the index of a column to swap
     *  @param b the index of a column to swap
     *  @return an \f$ M \times N \f$ matrix \f$ \textbf{m}' \f$ of type T such that
     *  \f$ {\textbf{m}'}_{ij} = \begin{cases} \textbf{m}_{ib} & j = a\\ \textbf{m}_{ia} & j = b\\
     * \textbf{m}_{ij} & \textrm{otherwise} \end{cases} \f$
     *
     *  Swap two rows of a matrix.
     */
    template<std::size_t M, std::size_t N, typename T>
    constexpr matrix<T, M, N> swapcol(matrix<T, M, N> m, std::size_t a, std::size_t b) {
        for (int i = 0; i < N; i++) {
            T tmp = m[i][a];
            m[i][a] = m[i][b];
            m[i][b] = tmp;
        }
        return m;
    }

    /** @brief horizontally concatenates two matrices
     *  @param a an \f$ M \times N \f$ matrix of type T
     *  @param b an \f$ M \times P \f$ matrix of type T
     *  @return an \f$ M \times \left(N+P\right) \f$ matrix of type T \f$ \left[\textbf{a} \textbf{b}\right] \f$
     * such that \f$ \left(\left[\textbf{a} \textbf{b}\right]\right)_{ij} = \begin{cases} \textbf{a}_{ij} & j <
     * N\\ \textbf{b}_{i,\ \left(j - N\right)} & j \ge N \end{cases} \f$
     *
     *  Horizontally concatenates two matrices.
     */
    template<std::size_t M, std::size_t N, std::size_t P, typename T>
    constexpr matrix<T, M, N + P> horzcat(const matrix<T, M, N> &a, const matrix<T, M, P> &b) {
        return generate<M, N + P>([&a, &b](std::size_t i, std::size_t j) { return j < N ? a[i][j] : b[i][j - N]; });
    }

    /** @brief vertically concatenates two matrices
     *  @param a an \f$ M \times P \f$ matrix of type T
     *  @param b an \f$ N \times P \f$ matrix of type T
     *  @return an \f$ \left(M+N\right) \times P \f$ matrix of type T \f$ \left[\textbf{a}; \textbf{b}\right]
     * \f$
     * such that \f$ \left(\left[\textbf{a}; \textbf{b}\right]\right)_{ij} = \begin{cases} \textbf{a}_{ij} & i <
     * M\\ \textbf{b}_{\left(i - M\right),\ j} & i \ge M \end{cases} \f$
     *
     *  Vertically concatenates two matrices.
     */
    template<std::size_t M, std::size_t N, std::size_t P, typename T>
    constexpr matrix<T, M + N, P> vertcat(const matrix<T, M, P> &a, const matrix<T, N, P> &b) {
        return generate<M + N, P>([&a, &b](std::size_t i, std::size_t j) { return i < M ? a[i][j] : b[i - M][j]; });
    }

    /** @brief extracts the submatrix of a matrix
     *  @param m an \f$ M \times N \f$ matrix of type T
     *  @param a the starting index into the rows
     *  @param b the starting index into the columns
     *  @return an \f$ P \times Q \f$ submatrix \f$ \textbf{m}' \f$ of type T such that
     *  \f$ {\textbf{m}'}_{ij} = \textbf{m}_{\left(a + i\right),\ \left(b + j\right)} \f$
     *
     *  Extracts the submatrix of a matrix.
     */
    template<std::size_t P, std::size_t Q, std::size_t M, std::size_t N, typename T>
    constexpr matrix<T, P, Q> submat(const matrix<T, M, N> &m, std::size_t a, std::size_t b) {
        if ((a + P > M) || (b + Q > N))
            throw "index out of range";
        return generate<P, Q>([&m, &a, &b](std::size_t i, std::size_t j) { return m[a + i][b + j]; });
    }

    /** @brief reshapes a matrix
     *  @param m an \f$ M \times N \f$ matrix of type T
     *  @return a \f$ P \times Q \f$ matrix of type T
     *
     *  Reshapes a matrix without changing the order of the elements in memory (in row-major order).
     */
    template<std::size_t P, std::size_t Q, std::size_t M, std::size_t N, typename T>
    constexpr matrix<T, P, Q> reshape(const matrix<T, M, N> &m) {
        static_assert(P * Q == M * N, "Reshaped matrix must preserve size! P*Q != M*N");
        return generate<P, Q>([&m](std::size_t i, std::size_t j) { return m[(i * Q + j) / N][(i * Q + j) % N]; });
    }

    /** @brief converts a vector into a column vector
     *  @param v an N-vector of type T
     *  @returns an \f$ N \times 1 \f$ matrix of type T
     *
     *  Converts a vector into a column vector.
     */
    template<typename T, std::size_t N>
    constexpr matrix<T, N, 1> as_column(const std::vector<T> &v) {
        return generate<N, 1>([&v](auto i, auto) { return v[i]; });
    }

    /** @brief converts a vector into a row vector
     *  @param v an N-vector of type T
     *  @returns a \f$ 1 \times N \f$ matrix of type T
     *
     *  Converts a vector into a row vector.
     */
    template<typename T, std::size_t N>
    constexpr matrix<T, 1, N> as_row(const std::vector<T> &v) {
        return generate<1, N>([&v](auto, auto j) { return v[j]; });
    }

    /** }@*/

}    // namespace nil::crypto3::math
