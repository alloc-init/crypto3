// MIT License
#pragma once

#include <cassert>
#include <cstddef>
#include <optional>
#include <type_traits>
#include <utility>

#include <nil/crypto3/math/matrix/matrix.hpp>
#include <nil/crypto3/math/matrix/vector.hpp>

namespace nil::crypto3::math {
    namespace detail {
        template<typename T>
        decltype(auto) matrix_inverse(const T &value) {
            if constexpr (requires { value.inversed(); }) {
                return value.inversed();
            } else {
                return value.inv();
            }
        }
    }    // namespace detail

    template<MatrixExpression Expression>
    bool is_square(const Expression &value) {
        return rows(value) == columns(value);
    }

    template<MatrixExpression Expression>
    bool is_upper_triangular(const Expression &value) {
        using value_type = typename std::remove_cvref_t<Expression>::value_type;
        if (!is_square(value)) {
            return false;
        }
        for (std::size_t i = 0; i < rows(value); ++i) {
            if (value(i, i) == value_type {}) {
                return false;
            }
            for (std::size_t j = 0; j < i; ++j) {
                if (value(i, j) != value_type {}) {
                    return false;
                }
            }
        }
        return true;
    }

    template<MatrixBackend MatrixBackendType, VectorBackend VectorBackendType>
        requires std::same_as<typename MatrixBackendType::value_type, typename VectorBackendType::value_type> &&
                 std::constructible_from<vector<VectorBackendType>, typename VectorBackendType::size_type>
    vector<VectorBackendType> back_substitute(const matrix<MatrixBackendType> &upper,
                                              const vector<VectorBackendType> &right_hand_side) {
        using value_type = typename MatrixBackendType::value_type;
        const std::size_t size = rows(upper);
        assert(is_upper_triangular(upper));
        assert(right_hand_side.size() == size);

        vector<VectorBackendType> result(size);
        for (std::size_t reverse_i = size; reverse_i > 0; --reverse_i) {
            const std::size_t i = reverse_i - 1;
            value_type sum = right_hand_side(i);
            for (std::size_t j = i + 1; j < size; ++j) {
                sum -= upper(i, j) * result(j);
            }
            result(i) = sum * detail::matrix_inverse(upper(i, i));
        }
        return result;
    }

    template<MatrixBackend MatrixBackendType, VectorBackend VectorBackendType>
        requires std::same_as<typename MatrixBackendType::value_type, typename VectorBackendType::value_type>
    std::optional<vector<VectorBackendType>>
        solve(matrix<MatrixBackendType> matrix, vector<VectorBackendType> right_hand_side) {
        using value_type = typename MatrixBackendType::value_type;

        const std::size_t row_count = matrix.rows();
        const std::size_t column_count = matrix.columns();
        if (right_hand_side.size() != row_count) {
            return std::nullopt;
        }

        std::size_t pivot_row = 0;
        std::size_t pivot_column = 0;
        while (pivot_row < row_count && pivot_column < column_count) {
            std::size_t nonzero_row = row_count;
            for (std::size_t i = pivot_row; i < row_count; ++i) {
                if (matrix(i, pivot_column) != value_type {}) {
                    nonzero_row = i;
                    break;
                }
            }
            if (nonzero_row == row_count) {
                ++pivot_column;
                continue;
            }

            if (nonzero_row != pivot_row) {
                for (std::size_t j = 0; j < column_count; ++j) {
                    std::swap(matrix(pivot_row, j), matrix(nonzero_row, j));
                }
                std::swap(right_hand_side(pivot_row), right_hand_side(nonzero_row));
            }

            const value_type pivot_inverse = detail::matrix_inverse(matrix(pivot_row, pivot_column));
            const value_type pivot_rhs = right_hand_side(pivot_row);
            for (std::size_t i = pivot_row + 1; i < row_count; ++i) {
                const value_type factor = matrix(i, pivot_column) * pivot_inverse;
                matrix(i, pivot_column) = value_type {};
                for (std::size_t j = pivot_column + 1; j < column_count; ++j) {
                    matrix(i, j) -= matrix(pivot_row, j) * factor;
                }
                right_hand_side(i) -= pivot_rhs * factor;
            }
            ++pivot_row;
            ++pivot_column;
        }

        for (std::size_t i = 0; i < row_count; ++i) {
            bool zero_row = true;
            for (std::size_t j = 0; j < column_count; ++j) {
                if (matrix(i, j) != value_type {}) {
                    zero_row = false;
                    break;
                }
            }
            if (zero_row) {
                return std::nullopt;
            }
        }

        if (!is_square(matrix)) {
            return std::nullopt;
        }
        return back_substitute(matrix, right_hand_side);
    }
}    // namespace nil::crypto3::math
