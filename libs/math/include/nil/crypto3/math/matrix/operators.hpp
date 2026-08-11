// MIT License
#pragma once

#include <nil/crypto3/math/matrix/matrix.hpp>
#include <nil/crypto3/math/matrix/vector.hpp>

namespace nil::crypto3::math {
    // Arithmetic returns the backend expression itself and therefore remains
    // lazy. The expression borrows both operands, so temporary frontends are
    // explicitly rejected below.
    template<MatrixBackend Left, MatrixBackend Right>
        requires requires(const Left &left, const Right &right) {
            requires MatrixExpression<decltype(left + right)>;
        }
    decltype(auto) operator+(const matrix<Left> &left, const matrix<Right> &right) {
        return left.backend() + right.backend();
    }

    template<MatrixBackend Left, MatrixBackend Right>
        requires requires(const Left &left, const Right &right) {
            requires MatrixExpression<decltype(left - right)>;
        }
    decltype(auto) operator-(const matrix<Left> &left, const matrix<Right> &right) {
        return left.backend() - right.backend();
    }

    template<VectorBackend Left, VectorBackend Right>
        requires requires(const Left &left, const Right &right) {
            requires VectorExpression<decltype(left + right)>;
        }
    decltype(auto) operator+(const vector<Left> &left, const vector<Right> &right) {
        return left.backend() + right.backend();
    }

    template<VectorBackend Left, VectorBackend Right>
        requires requires(const Left &left, const Right &right) {
            requires VectorExpression<decltype(left - right)>;
        }
    decltype(auto) operator-(const vector<Left> &left, const vector<Right> &right) {
        return left.backend() - right.backend();
    }

    template<typename Left, typename Right>
    void operator+(matrix<Left> &&, const matrix<Right> &) = delete;

    template<typename Left, typename Right>
    void operator+(const matrix<Left> &, matrix<Right> &&) = delete;

    template<typename Left, typename Right>
    void operator+(matrix<Left> &&, matrix<Right> &&) = delete;

    template<typename Left, typename Right>
    void operator-(matrix<Left> &&, const matrix<Right> &) = delete;

    template<typename Left, typename Right>
    void operator-(const matrix<Left> &, matrix<Right> &&) = delete;

    template<typename Left, typename Right>
    void operator-(matrix<Left> &&, matrix<Right> &&) = delete;

    template<typename Left, typename Right>
    void operator+(vector<Left> &&, const vector<Right> &) = delete;

    template<typename Left, typename Right>
    void operator+(const vector<Left> &, vector<Right> &&) = delete;

    template<typename Left, typename Right>
    void operator+(vector<Left> &&, vector<Right> &&) = delete;

    template<typename Left, typename Right>
    void operator-(vector<Left> &&, const vector<Right> &) = delete;

    template<typename Left, typename Right>
    void operator-(const vector<Left> &, vector<Right> &&) = delete;

    template<typename Left, typename Right>
    void operator-(vector<Left> &&, vector<Right> &&) = delete;
}    // namespace nil::crypto3::math
