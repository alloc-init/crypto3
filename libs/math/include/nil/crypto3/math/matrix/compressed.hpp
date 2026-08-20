#pragma once

#include <cstddef>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>

#include <nil/crypto3/math/matrix/matrix.hpp>
#include <nil/crypto3/math/matrix/operators.hpp>
#include <nil/crypto3/math/matrix/vector.hpp>

namespace nil::crypto3::math {
    template<typename T, typename Layout = boost::numeric::ublas::row_major, std::size_t IndexBase = 0,
             typename IndexStorage = boost::numeric::ublas::unbounded_array<std::size_t>,
             typename ValueStorage = boost::numeric::ublas::unbounded_array<T>>
    using compressed_matrix =
        matrix<boost::numeric::ublas::compressed_matrix<T, Layout, IndexBase, IndexStorage, ValueStorage>>;

    template<typename T, std::size_t IndexBase = 0,
             typename IndexStorage = boost::numeric::ublas::unbounded_array<std::size_t>,
             typename ValueStorage = boost::numeric::ublas::unbounded_array<T>>
    using compressed_vector =
        vector<boost::numeric::ublas::compressed_vector<T, IndexBase, IndexStorage, ValueStorage>>;

    template<typename Backend>
        requires requires(const Backend &backend, std::size_t row, std::size_t column) {
            backend.find_element(row, column);
        }
    const typename Backend::value_type *find_element(const matrix<Backend> &value, std::size_t row,
                                                     std::size_t column) {
        return value.backend().find_element(row, column);
    }

    // An absent sparse source element leaves the destination unchanged.
    template<typename SourceBackend, typename DestinationBackend>
        requires requires(const SourceBackend &source, DestinationBackend &destination, std::size_t row,
                          std::size_t column) {
            source.find_element(row, column);
            destination(row, column) = *source.find_element(row, column);
        }
    void assign_if_stored(const matrix<SourceBackend> &source, std::size_t source_row,
                          std::size_t source_column, matrix<DestinationBackend> &destination,
                          std::size_t destination_row, std::size_t destination_column) {
        if (auto element = source.backend().find_element(source_row, source_column)) {
            destination(destination_row, destination_column) = *element;
        }
    }

    template<typename Backend>
        requires requires(const Backend &backend, std::size_t index) { backend.find_element(index); }
    const typename Backend::value_type *find_element(const vector<Backend> &value, std::size_t index) {
        return value.backend().find_element(index);
    }

    template<typename Backend, typename Function>
        requires requires(const Backend &backend) {
            backend.begin1();
            backend.end1();
        }
    void for_each_nonzero(const matrix<Backend> &value, Function &&function) {
        for (auto outer = value.backend().begin1(); outer != value.backend().end1(); ++outer) {
            for (auto element = outer.begin(); element != outer.end(); ++element) {
                function(element.index1(), element.index2(), *element);
            }
        }
    }

    template<typename Backend, typename Function>
        requires requires(const Backend &backend) {
            backend.begin();
            backend.end();
        }
    void for_each_nonzero(const vector<Backend> &value, Function &&function) {
        for (auto element = value.backend().begin(); element != value.backend().end(); ++element) {
            function(element.index(), *element);
        }
    }
}    // namespace nil::crypto3::math
