#pragma once

#include <cstddef>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>

#include <nil/crypto3/math/matrix/matrix.hpp>
#include <nil/crypto3/math/matrix/operators.hpp>
#include <nil/crypto3/math/matrix/vector.hpp>

namespace nil::crypto3::math {
    namespace backends::ublas {
        template<typename T,
                 typename Layout = boost::numeric::ublas::row_major,
                 std::size_t IndexBase = 0,
                 typename IndexStorage = boost::numeric::ublas::unbounded_array<std::size_t>,
                 typename ValueStorage = boost::numeric::ublas::unbounded_array<T>>
        using compressed_matrix =
            boost::numeric::ublas::compressed_matrix<T, Layout, IndexBase, IndexStorage, ValueStorage>;

        template<typename T,
                 std::size_t IndexBase = 0,
                 typename IndexStorage = boost::numeric::ublas::unbounded_array<std::size_t>,
                 typename ValueStorage = boost::numeric::ublas::unbounded_array<T>>
        using compressed_vector = boost::numeric::ublas::compressed_vector<T, IndexBase, IndexStorage, ValueStorage>;
    }    // namespace backends::ublas

    template<typename T,
             typename Layout = boost::numeric::ublas::row_major,
             std::size_t IndexBase = 0,
             typename IndexStorage = boost::numeric::ublas::unbounded_array<std::size_t>,
             typename ValueStorage = boost::numeric::ublas::unbounded_array<T>>
    using compressed_matrix = matrix<backends::ublas::compressed_matrix<T, Layout, IndexBase, IndexStorage, ValueStorage>>;

    template<typename T,
             std::size_t IndexBase = 0,
             typename IndexStorage = boost::numeric::ublas::unbounded_array<std::size_t>,
             typename ValueStorage = boost::numeric::ublas::unbounded_array<T>>
    using compressed_vector = vector<backends::ublas::compressed_vector<T, IndexBase, IndexStorage, ValueStorage>>;
}    // namespace nil::crypto3::math
