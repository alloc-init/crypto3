//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// MIT License
//---------------------------------------------------------------------------//
#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>

#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/types/array_list.hpp>
#include <nil/marshalling/types/bundle.hpp>
#include <nil/marshalling/types/integral.hpp>

#include <nil/crypto3/marshalling/algebra/types/field_element.hpp>
#include <nil/crypto3/math/matrix/compressed.hpp>
#include <nil/crypto3/math/matrix/regular.hpp>

namespace nil::crypto3::marshalling::types {

    template<typename TTypeBase, typename MatrixType>
    struct regular_matrix {
        using index_type = nil::marshalling::types::integral<TTypeBase, std::size_t>;
        using value_type = field_element<TTypeBase, typename MatrixType::value_type>;
        using values_type = nil::marshalling::types::standard_array_list<TTypeBase, value_type>;
        using type = nil::marshalling::types::bundle<
            TTypeBase, std::tuple<index_type, index_type, values_type>>;
    };

    template<typename Endianness, typename MatrixType>
    typename regular_matrix<nil::marshalling::field_type<Endianness>, MatrixType>::type
        fill_regular_matrix(const MatrixType &matrix) {
        using TTypeBase = nil::marshalling::field_type<Endianness>;
        using marshalling_type = regular_matrix<TTypeBase, MatrixType>;

        typename marshalling_type::values_type values;
        values.value().reserve(matrix.rows() * matrix.columns());
        for (std::size_t row = 0; row < matrix.rows(); ++row) {
            for (std::size_t column = 0; column < matrix.columns(); ++column) {
                values.value().emplace_back(matrix(row, column));
            }
        }

        return typename marshalling_type::type(std::make_tuple(
            typename marshalling_type::index_type(matrix.rows()),
            typename marshalling_type::index_type(matrix.columns()),
            std::move(values)));
    }

    template<typename Endianness, typename MatrixType>
    MatrixType make_regular_matrix(
        const typename regular_matrix<nil::marshalling::field_type<Endianness>, MatrixType>::type &filled_matrix) {
        const auto &fields = filled_matrix.value();
        const std::size_t rows = std::get<0>(fields).value();
        const std::size_t columns = std::get<1>(fields).value();
        const auto &values = std::get<2>(fields).value();
        if (columns != 0 && rows > std::numeric_limits<std::size_t>::max() / columns) {
            throw std::invalid_argument("regular matrix dimensions overflow");
        }
        if (values.size() != rows * columns) {
            throw std::invalid_argument("regular matrix element count does not match its dimensions");
        }

        MatrixType matrix(rows, columns);
        auto value = values.begin();
        for (std::size_t row = 0; row < rows; ++row) {
            for (std::size_t column = 0; column < columns; ++column, ++value) {
                matrix(row, column) = value->value();
            }
        }
        return matrix;
    }

    template<typename TTypeBase, typename VectorType>
    struct regular_vector {
        using value_type = field_element<TTypeBase, typename VectorType::value_type>;
        using type = nil::marshalling::types::standard_array_list<TTypeBase, value_type>;
    };

    template<typename Endianness, typename VectorType>
    typename regular_vector<nil::marshalling::field_type<Endianness>, VectorType>::type
        fill_regular_vector(const VectorType &vector) {
        using TTypeBase = nil::marshalling::field_type<Endianness>;
        using marshalling_type = regular_vector<TTypeBase, VectorType>;

        typename marshalling_type::type values;
        values.value().reserve(vector.size());
        for (std::size_t index = 0; index < vector.size(); ++index) {
            values.value().emplace_back(vector[index]);
        }
        return values;
    }

    template<typename Endianness, typename VectorType>
    VectorType make_regular_vector(
        const typename regular_vector<nil::marshalling::field_type<Endianness>, VectorType>::type &filled_vector) {
        const auto &values = filled_vector.value();
        VectorType vector(values.size());
        for (std::size_t index = 0; index < values.size(); ++index) {
            vector[index] = values[index].value();
        }
        return vector;
    }

    template<typename TTypeBase, typename MatrixType>
    struct compressed_matrix {
        using index_type = nil::marshalling::types::integral<TTypeBase, std::size_t>;
        using value_type = field_element<TTypeBase, typename MatrixType::value_type>;
        using entry_type = nil::marshalling::types::bundle<
            TTypeBase, std::tuple<index_type, index_type, value_type>>;
        using entries_type = nil::marshalling::types::standard_array_list<TTypeBase, entry_type>;
        using type = nil::marshalling::types::bundle<
            TTypeBase, std::tuple<index_type, index_type, entries_type>>;
    };

    template<typename Endianness, typename MatrixType>
    typename compressed_matrix<nil::marshalling::field_type<Endianness>, MatrixType>::type
        fill_compressed_matrix(const MatrixType &matrix) {
        using TTypeBase = nil::marshalling::field_type<Endianness>;
        using marshalling_type = compressed_matrix<TTypeBase, MatrixType>;

        typename marshalling_type::entries_type entries;
        nil::crypto3::math::for_each_nonzero(
            matrix, [&](std::size_t row, std::size_t column, const typename MatrixType::value_type &value) {
                entries.value().emplace_back(std::make_tuple(
                    typename marshalling_type::index_type(row),
                    typename marshalling_type::index_type(column),
                    typename marshalling_type::value_type(value)));
            });

        return typename marshalling_type::type(std::make_tuple(
            typename marshalling_type::index_type(matrix.rows()),
            typename marshalling_type::index_type(matrix.columns()),
            std::move(entries)));
    }

    template<typename Endianness, typename MatrixType>
    MatrixType make_compressed_matrix(
        const typename compressed_matrix<nil::marshalling::field_type<Endianness>, MatrixType>::type &filled_matrix) {
        const auto &fields = filled_matrix.value();
        const std::size_t rows = std::get<0>(fields).value();
        const std::size_t columns = std::get<1>(fields).value();
        MatrixType matrix(rows, columns);

        for (const auto &filled_entry : std::get<2>(fields).value()) {
            const auto &entry = filled_entry.value();
            const std::size_t row = std::get<0>(entry).value();
            const std::size_t column = std::get<1>(entry).value();
            if (row >= rows || column >= columns) {
                throw std::invalid_argument("compressed matrix entry is out of bounds");
            }
            if (nil::crypto3::math::find_element(matrix, row, column) != nullptr) {
                throw std::invalid_argument("compressed matrix contains a duplicate entry");
            }
            matrix(row, column) = std::get<2>(entry).value();
        }
        return matrix;
    }

    template<typename TTypeBase, typename VectorType>
    struct compressed_vector {
        using index_type = nil::marshalling::types::integral<TTypeBase, std::size_t>;
        using value_type = field_element<TTypeBase, typename VectorType::value_type>;
        using entry_type = nil::marshalling::types::bundle<TTypeBase, std::tuple<index_type, value_type>>;
        using entries_type = nil::marshalling::types::standard_array_list<TTypeBase, entry_type>;
        using type = nil::marshalling::types::bundle<TTypeBase, std::tuple<index_type, entries_type>>;
    };

    template<typename Endianness, typename VectorType>
    typename compressed_vector<nil::marshalling::field_type<Endianness>, VectorType>::type
        fill_compressed_vector(const VectorType &vector) {
        using TTypeBase = nil::marshalling::field_type<Endianness>;
        using marshalling_type = compressed_vector<TTypeBase, VectorType>;

        typename marshalling_type::entries_type entries;
        nil::crypto3::math::for_each_nonzero(
            vector, [&](std::size_t index, const typename VectorType::value_type &value) {
                entries.value().emplace_back(std::make_tuple(
                    typename marshalling_type::index_type(index),
                    typename marshalling_type::value_type(value)));
            });

        return typename marshalling_type::type(std::make_tuple(
            typename marshalling_type::index_type(vector.size()), std::move(entries)));
    }

    template<typename Endianness, typename VectorType>
    VectorType make_compressed_vector(
        const typename compressed_vector<nil::marshalling::field_type<Endianness>, VectorType>::type &filled_vector) {
        const auto &fields = filled_vector.value();
        const std::size_t size = std::get<0>(fields).value();
        VectorType vector(size);

        for (const auto &filled_entry : std::get<1>(fields).value()) {
            const auto &entry = filled_entry.value();
            const std::size_t index = std::get<0>(entry).value();
            if (index >= size) {
                throw std::invalid_argument("compressed vector entry is out of bounds");
            }
            if (nil::crypto3::math::find_element(vector, index) != nullptr) {
                throw std::invalid_argument("compressed vector contains a duplicate entry");
            }
            vector[index] = std::get<1>(entry).value();
        }
        return vector;
    }

}    // namespace nil::crypto3::marshalling::types
