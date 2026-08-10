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

namespace nil::crypto3::math {
    template<typename T, std::size_t N, std::size_t M>
    struct static_matrix {
        static_assert(N != 0 && M != 0, "static_matrix must have have positive dimensions");

        constexpr static_matrix() : arrays {} {
        }

        constexpr static_matrix(const T (&array)[N][M]) {
            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < M; ++j) {
                    arrays[i][j] = array[i][j];
                }
            }
        }

        template<typename... Args>
        constexpr static_matrix(Args... args) : arrays {std::forward<Args>(args)...} {
            static_assert(sizeof...(args) == N * M, "Number of arguments must match the static_matrix size");
        }

        using value_type = T;
        using size_type = std::size_t;
        constexpr static const size_type column_size = N;
        constexpr static const size_type row_size = M;

        constexpr static_vector<T, M> row(size_type i) const {
            if (i >= N) {
                throw "index out of range";
            }
            return generate<M>([i, this](size_type j) { return arrays[i][j]; });
        }

        constexpr static_vector<T, N> column(size_type i) const {
            if (i >= M)
                throw "index out of range";
            return generate<N>([i, this](size_type j) { return arrays[j][i]; });
        }

        constexpr T *operator[](size_type i) {
            return arrays[i];
        }

        constexpr T const *operator[](size_type i) const {
            return arrays[i];
        }

        constexpr bool operator==(const static_matrix &other) const {
            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < M; ++j) {
                    if (arrays[i][j] != other.arrays[i][j]) {
                        return false;
                    }
                }
            }
            return true;
        }

        constexpr const T &operator()(size_type i, size_type j) const {
            return arrays[i][j];
        }

        T arrays[N][M];    ///< @private
    };

    template<typename T, std::size_t M, std::size_t N>
    static_matrix(const T (&)[M][N]) -> static_matrix<T, M, N>;
}    // namespace nil::crypto3::math
