// MIT License
#pragma once

#include <cstddef>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <nil/crypto3/math/matrix/matrix.hpp>
#include <nil/crypto3/math/matrix/vector.hpp>

namespace nil::crypto3::math {
    // Arithmetic returns the backend expression itself and therefore remains
    // lazy. The expression borrows both operands, so temporary frontends are
    // explicitly rejected below.
    template<MatrixBackend Left, MatrixBackend Right>
        requires requires(const Left &left, const Right &right) { requires MatrixExpression<decltype(left + right)>; }
    decltype(auto) operator+(const matrix<Left> &left, const matrix<Right> &right) {
        return left.backend() + right.backend();
    }

    template<MatrixBackend Left, MatrixBackend Right>
        requires requires(const Left &left, const Right &right) { requires MatrixExpression<decltype(left - right)>; }
    decltype(auto) operator-(const matrix<Left> &left, const matrix<Right> &right) {
        return left.backend() - right.backend();
    }

    template<VectorBackend Left, VectorBackend Right>
        requires requires(const Left &left, const Right &right) { requires VectorExpression<decltype(left + right)>; }
    decltype(auto) operator+(const vector<Left> &left, const vector<Right> &right) {
        return left.backend() + right.backend();
    }

    template<VectorBackend Left, VectorBackend Right>
        requires requires(const Left &left, const Right &right) { requires VectorExpression<decltype(left - right)>; }
    decltype(auto) operator-(const vector<Left> &left, const vector<Right> &right) {
        return left.backend() - right.backend();
    }

    template<MatrixBackend Left, MatrixBackend Right>
        requires requires(const Left &left, const Right &right) {
            requires MatrixExpression<decltype(boost::numeric::ublas::prod(left, right))>;
        }
    decltype(auto) product(const matrix<Left> &left, const matrix<Right> &right) {
        return boost::numeric::ublas::prod(left.backend(), right.backend());
    }

    template<MatrixBackend Left, VectorBackend Right>
        requires requires(const Left &left, const Right &right) {
            requires VectorExpression<decltype(boost::numeric::ublas::prod(left, right))>;
        }
    decltype(auto) product(const matrix<Left> &left, const vector<Right> &right) {
        return boost::numeric::ublas::prod(left.backend(), right.backend());
    }

    template<MatrixExpression Left, VectorBackend Right>
        requires requires(const Left &left, const Right &right) {
            requires VectorExpression<decltype(boost::numeric::ublas::prod(left, right))>;
        }
    decltype(auto) product(const Left &left, const vector<Right> &right) {
        return boost::numeric::ublas::prod(left, right.backend());
    }

    template<VectorBackend Left, VectorBackend Right>
        requires requires(const Left &left, const Right &right) { boost::numeric::ublas::inner_prod(left, right); }
    decltype(auto) inner_product(const vector<Left> &left, const vector<Right> &right) {
        return boost::numeric::ublas::inner_prod(left.backend(), right.backend());
    }

    template<VectorExpression Left, VectorBackend Right>
        requires requires(const Left &left, const Right &right) { boost::numeric::ublas::inner_prod(left, right); }
    decltype(auto) inner_product(const Left &left, const vector<Right> &right) {
        return boost::numeric::ublas::inner_prod(left, right.backend());
    }

    template<VectorBackend Left, VectorBackend Right>
        requires requires(const Left &left, const Right &right) {
            requires VectorExpression<decltype(boost::numeric::ublas::element_prod(left, right))>;
        }
    decltype(auto) element_product(const vector<Left> &left, const vector<Right> &right) {
        return boost::numeric::ublas::element_prod(left.backend(), right.backend());
    }

    template<MatrixBackend Backend>
    decltype(auto) row(matrix<Backend> &value, std::size_t index) {
        return boost::numeric::ublas::row(value.backend(), index);
    }

    template<MatrixBackend Backend>
    decltype(auto) row(const matrix<Backend> &value, std::size_t index) {
        return boost::numeric::ublas::row(value.backend(), index);
    }

    template<MatrixBackend Backend>
    decltype(auto) column(matrix<Backend> &value, std::size_t index) {
        return boost::numeric::ublas::column(value.backend(), index);
    }

    template<MatrixBackend Backend>
    decltype(auto) column(const matrix<Backend> &value, std::size_t index) {
        return boost::numeric::ublas::column(value.backend(), index);
    }

    template<MatrixBackend Backend>
    decltype(auto) submatrix(matrix<Backend> &value,
                             std::size_t row_begin,
                             std::size_t row_end,
                             std::size_t column_begin,
                             std::size_t column_end) {
        using boost::numeric::ublas::range;
        return boost::numeric::ublas::project(
            value.backend(), range(row_begin, row_end), range(column_begin, column_end));
    }

    template<MatrixBackend Backend>
    decltype(auto) submatrix(const matrix<Backend> &value,
                             std::size_t row_begin,
                             std::size_t row_end,
                             std::size_t column_begin,
                             std::size_t column_end) {
        using boost::numeric::ublas::range;
        return boost::numeric::ublas::project(
            value.backend(), range(row_begin, row_end), range(column_begin, column_end));
    }

    template<VectorExpression Expression>
    decltype(auto) subvector(Expression &value, std::size_t begin, std::size_t end) {
        return boost::numeric::ublas::subrange(value, begin, end);
    }

    template<VectorExpression Expression>
    decltype(auto) subvector(const Expression &value, std::size_t begin, std::size_t end) {
        return boost::numeric::ublas::subrange(value, begin, end);
    }

    template<MatrixBackend Backend>
    void swap_rows(matrix<Backend> &value, std::size_t left, std::size_t right) {
        row(value, left).swap(row(value, right));
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
