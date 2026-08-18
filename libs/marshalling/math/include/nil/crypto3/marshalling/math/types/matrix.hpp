//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// MIT License
//---------------------------------------------------------------------------//
#pragma once

#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <utility>

#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/types/array_list.hpp>
#include <nil/marshalling/types/bundle.hpp>
#include <nil/marshalling/types/integral.hpp>

#include <nil/crypto3/marshalling/algebra/types/field_element.hpp>
#include <nil/crypto3/math/matrix/compressed.hpp>

namespace nil::crypto3::marshalling::types {

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

}    // namespace nil::crypto3::marshalling::types
