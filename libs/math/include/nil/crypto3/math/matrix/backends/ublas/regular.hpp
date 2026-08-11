#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <nil/crypto3/math/matrix/matrix.hpp>
#include <nil/crypto3/math/matrix/operators.hpp>
#include <nil/crypto3/math/matrix/vector.hpp>

namespace nil::crypto3::math {
    namespace backends::ublas {
        template<typename T,
                 typename Layout = boost::numeric::ublas::row_major,
                 typename Storage = boost::numeric::ublas::unbounded_array<T>>
        using regular_matrix = boost::numeric::ublas::matrix<T, Layout, Storage>;

        template<typename T, typename Storage = boost::numeric::ublas::unbounded_array<T>>
        using regular_vector = boost::numeric::ublas::vector<T, Storage>;
    }    // namespace backends::ublas

    template<typename T,
             typename Layout = boost::numeric::ublas::row_major,
             typename Storage = boost::numeric::ublas::unbounded_array<T>>
    using regular_matrix = matrix<backends::ublas::regular_matrix<T, Layout, Storage>>;

    template<typename T, typename Storage = boost::numeric::ublas::unbounded_array<T>>
    using regular_vector = vector<backends::ublas::regular_vector<T, Storage>>;
}    // namespace nil::crypto3::math
