#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <nil/crypto3/math/matrix/matrix.hpp>
#include <nil/crypto3/math/matrix/operators.hpp>
#include <nil/crypto3/math/matrix/vector.hpp>

namespace nil::crypto3::math {
    template<typename T, typename Layout = boost::numeric::ublas::row_major,
             typename Storage = boost::numeric::ublas::unbounded_array<T>>
    using regular_matrix = matrix<boost::numeric::ublas::matrix<T, Layout, Storage>>;

    template<typename T, typename Storage = boost::numeric::ublas::unbounded_array<T>>
    using regular_vector = vector<boost::numeric::ublas::vector<T, Storage>>;

    template<typename T>
    regular_matrix<T> identity_matrix(std::size_t size) {
        return regular_matrix<T>(boost::numeric::ublas::identity_matrix<T>(size));
    }
}    // namespace nil::crypto3::math
