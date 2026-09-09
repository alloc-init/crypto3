//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// MIT License
//---------------------------------------------------------------------------//
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/options.hpp>
#include <nil/marshalling/status_type.hpp>
#include <nil/marshalling/types/array_list.hpp>
#include <nil/marshalling/types/bundle.hpp>
#include <nil/marshalling/types/integral.hpp>

#include <nil/crypto3/marshalling/algebra/types/field_element.hpp>
#include <nil/crypto3/math/matrix/compressed.hpp>
#include <nil/crypto3/math/matrix/regular.hpp>

namespace nil::crypto3::marshalling::types {

    namespace detail {

        /** Validate raw base-field components before conversion can reduce them modulo the field modulus. */
        template<typename FieldValueType>
        struct canonical_field_element_validator {
            template<typename MarshalledField>
            bool operator()(const MarshalledField &field) const {
                if constexpr (algebra::ExtendedFieldValue<FieldValueType>) {
                    const auto &components = field.value();
                    if (components.size() != FieldValueType::field_type::arity) {
                        return false;
                    }
                    for (const auto &component : components) {
                        if (component.value() >= FieldValueType::field_type::modulus) {
                            return false;
                        }
                    }
                    return true;
                } else {
                    return field.value() < FieldValueType::field_type::modulus;
                }
            }
        };

        template<typename TTypeBase, typename FieldValueType>
        using validated_field_element = field_element<
            TTypeBase, FieldValueType,
            nil::marshalling::option::contents_validator<canonical_field_element_validator<FieldValueType>>>;

        /**
         * Existing marshalling value writers require a random-access output buffer, but incremental serialization must
         * also support sequential output iterators such as stream and back-insert iterators. Serialize one value into a
         * local fixed-size byte array, then copy its bytes sequentially to the caller's output iterator.
         */
        template<typename MarshalledType, typename OutputIterator>
        nil::marshalling::status_type write_marshaled_value(const MarshalledType &value, OutputIterator &output) {
            std::array<std::uint8_t, MarshalledType::max_length()> scratch;
            auto scratch_output = scratch.begin();
            const std::size_t length = value.length();
            const nil::marshalling::status_type status = value.write(scratch_output, length);
            if (status != nil::marshalling::status_type::success) {
                return status;
            }

            for (std::size_t index = 0; index < length; ++index) {
                *output = scratch[index];
                ++output;
            }
            return nil::marshalling::status_type::success;
        }

        /**
         * Read one fixed-length value directly when the input is random-access. For a sequential input iterator, copy
         * one value into a local fixed-size array required by the existing marshalling reader.
         */
        template<typename MarshalledType, typename InputIterator>
        nil::marshalling::status_type read_marshaled_value(MarshalledType &value, InputIterator &input,
                                                           std::size_t &available_size) {
            constexpr std::size_t length = MarshalledType::max_length();
            if (available_size < length) {
                return nil::marshalling::status_type::not_enough_data;
            }

            using iterator_category = typename std::iterator_traits<std::decay_t<InputIterator>>::iterator_category;
            nil::marshalling::status_type read_status;
            if constexpr (std::is_base_of_v<std::random_access_iterator_tag, iterator_category>) {
                read_status = value.read(input, length);
            } else {
                std::array<std::uint8_t, length> scratch;
                for (std::size_t index = 0; index < length; ++index) {
                    scratch[index] = static_cast<std::uint8_t>(*input);
                    ++input;
                }
                auto scratch_input = scratch.begin();
                read_status = value.read(scratch_input, length);
            }
            if (read_status != nil::marshalling::status_type::success) {
                return read_status;
            }
            available_size -= length;
            if (!value.valid() || value.length() != length) {
                return nil::marshalling::status_type::invalid_msg_data;
            }
            return nil::marshalling::status_type::success;
        }

    }    // namespace detail

    template<typename TTypeBase, typename MatrixType>
    struct regular_matrix {
        using index_type = nil::marshalling::types::integral<TTypeBase, std::size_t>;
        using value_type = field_element<TTypeBase, typename MatrixType::value_type>;
        using values_type = nil::marshalling::types::standard_array_list<TTypeBase, value_type>;
        using element_count_type = typename values_type::parsed_options_type::sequence_size_field_prefix;
        using type = nil::marshalling::types::bundle<TTypeBase, std::tuple<index_type, index_type, values_type>>;
    };

    /**
     * Compute the fixed encoded size of a regular matrix without constructing its marshalled element list.
     * encoded_size is modified only on success.
     *
     *     encoded size = 2 * index length + element-count-prefix length
     *                    + rows * columns * field-element length
     *
     * @return invalid_msg_data if the element count or encoded size overflows std::size_t; success otherwise.
     */
    template<typename Endianness, typename MatrixType>
    nil::marshalling::status_type regular_matrix_encoded_size(std::size_t rows, std::size_t columns,
                                                              std::size_t &encoded_size) {
        using TTypeBase = nil::marshalling::field_type<Endianness>;
        using marshalling_type = regular_matrix<TTypeBase, MatrixType>;

        constexpr std::size_t index_length = marshalling_type::index_type::max_length();
        constexpr std::size_t element_count_length = marshalling_type::element_count_type::max_length();
        constexpr std::size_t element_length = marshalling_type::value_type::max_length();
        constexpr std::size_t maximum_size = std::numeric_limits<std::size_t>::max();

        if (columns != 0 && rows > maximum_size / columns) {
            return nil::marshalling::status_type::invalid_msg_data;
        }
        const std::size_t element_count = rows * columns;

        if (index_length > maximum_size - index_length) {
            return nil::marshalling::status_type::invalid_msg_data;
        }
        const std::size_t dimension_length = index_length + index_length;
        if (element_count_length > maximum_size - dimension_length) {
            return nil::marshalling::status_type::invalid_msg_data;
        }
        const std::size_t header_length = dimension_length + element_count_length;
        if (element_count != 0 && element_length > (maximum_size - header_length) / element_count) {
            return nil::marshalling::status_type::invalid_msg_data;
        }

        encoded_size = header_length + element_count * element_length;
        return nil::marshalling::status_type::success;
    }

    /**
     * Write a regular matrix directly to an output iterator without constructing a marshalled element list.
     * The wire format is identical to regular_matrix::type: rows, columns, element count, then row-major values.
     *
     * @return invalid_msg_data if dimensions or encoded size overflow, buffer_overflow if available_size is too small,
     *         or the status returned by an individual marshalling field write.
     */
    template<typename Endianness, typename MatrixType, typename OutputIterator>
    nil::marshalling::status_type write_regular_matrix(const MatrixType &matrix, OutputIterator &output,
                                                       std::size_t available_size) {
        using TTypeBase = nil::marshalling::field_type<Endianness>;
        using marshalling_type = regular_matrix<TTypeBase, MatrixType>;

        const std::size_t rows = matrix.rows();
        const std::size_t columns = matrix.columns();
        std::size_t encoded_size = 0;
        const nil::marshalling::status_type size_status =
            regular_matrix_encoded_size<Endianness, MatrixType>(rows, columns, encoded_size);
        if (size_status != nil::marshalling::status_type::success) {
            return size_status;
        }
        if (available_size < encoded_size) {
            return nil::marshalling::status_type::buffer_overflow;
        }

        const typename marshalling_type::index_type row_count(rows);
        nil::marshalling::status_type status = detail::write_marshaled_value(row_count, output);
        if (status != nil::marshalling::status_type::success) {
            return status;
        }
        const typename marshalling_type::index_type column_count(columns);
        status = detail::write_marshaled_value(column_count, output);
        if (status != nil::marshalling::status_type::success) {
            return status;
        }
        const typename marshalling_type::element_count_type element_count(rows * columns);
        status = detail::write_marshaled_value(element_count, output);
        if (status != nil::marshalling::status_type::success) {
            return status;
        }

        for (std::size_t row = 0; row < rows; ++row) {
            for (std::size_t column = 0; column < columns; ++column) {
                const typename marshalling_type::value_type value(matrix(row, column));
                status = detail::write_marshaled_value(value, output);
                if (status != nil::marshalling::status_type::success) {
                    return status;
                }
            }
        }
        return nil::marshalling::status_type::success;
    }

    /**
     * Read a regular matrix incrementally and pass each decoded field value to visitor(row, column, value).
     * The wire format is identical to regular_matrix::type: rows, columns, element count, then row-major values.
     * No callbacks are made until the complete header and payload size have been validated.
     *
     * @return not_enough_data for truncated input, invalid_msg_data for overflowing dimensions, an incorrect element
     *         count, or a noncanonical field element; otherwise the status returned by a marshalling field read.
     */
    template<typename Endianness, typename MatrixType, typename InputIterator, typename Visitor>
    nil::marshalling::status_type read_regular_matrix(InputIterator &input, std::size_t available_size,
                                                      std::size_t &rows, std::size_t &columns, Visitor &&visitor) {
        using TTypeBase = nil::marshalling::field_type<Endianness>;
        using marshalling_type = regular_matrix<TTypeBase, MatrixType>;

        std::size_t remaining_size = available_size;
        typename marshalling_type::index_type encoded_rows;
        nil::marshalling::status_type status = detail::read_marshaled_value(encoded_rows, input, remaining_size);
        if (status != nil::marshalling::status_type::success) {
            return status;
        }

        typename marshalling_type::index_type encoded_columns;
        status = detail::read_marshaled_value(encoded_columns, input, remaining_size);
        if (status != nil::marshalling::status_type::success) {
            return status;
        }

        typename marshalling_type::element_count_type encoded_element_count;
        status = detail::read_marshaled_value(encoded_element_count, input, remaining_size);
        if (status != nil::marshalling::status_type::success) {
            return status;
        }

        const std::size_t decoded_rows = encoded_rows.value();
        const std::size_t decoded_columns = encoded_columns.value();
        std::size_t encoded_size = 0;
        status = regular_matrix_encoded_size<Endianness, MatrixType>(decoded_rows, decoded_columns, encoded_size);
        if (status != nil::marshalling::status_type::success) {
            return status;
        }

        const std::size_t expected_element_count = decoded_rows * decoded_columns;
        if (encoded_element_count.value() != expected_element_count) {
            return nil::marshalling::status_type::invalid_msg_data;
        }
        if (available_size < encoded_size) {
            return nil::marshalling::status_type::not_enough_data;
        }

        rows = decoded_rows;
        columns = decoded_columns;
        for (std::size_t row = 0; row < decoded_rows; ++row) {
            for (std::size_t column = 0; column < decoded_columns; ++column) {
                detail::validated_field_element<TTypeBase, typename MatrixType::value_type> value;
                status = detail::read_marshaled_value(value, input, remaining_size);
                if (status != nil::marshalling::status_type::success) {
                    return status;
                }
                visitor(row, column, value.value());
            }
        }
        return nil::marshalling::status_type::success;
    }

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

        return typename marshalling_type::type(std::make_tuple(typename marshalling_type::index_type(matrix.rows()),
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
        using entry_type = nil::marshalling::types::bundle<TTypeBase, std::tuple<index_type, index_type, value_type>>;
        using entries_type = nil::marshalling::types::standard_array_list<TTypeBase, entry_type>;
        using type = nil::marshalling::types::bundle<TTypeBase, std::tuple<index_type, index_type, entries_type>>;
    };

    template<typename Endianness, typename MatrixType>
    typename compressed_matrix<nil::marshalling::field_type<Endianness>, MatrixType>::type
        fill_compressed_matrix(const MatrixType &matrix) {
        using TTypeBase = nil::marshalling::field_type<Endianness>;
        using marshalling_type = compressed_matrix<TTypeBase, MatrixType>;

        typename marshalling_type::entries_type entries;
        nil::crypto3::math::for_each_nonzero(
            matrix, [&](std::size_t row, std::size_t column, const typename MatrixType::value_type &value) {
                entries.value().emplace_back(std::make_tuple(typename marshalling_type::index_type(row),
                                                             typename marshalling_type::index_type(column),
                                                             typename marshalling_type::value_type(value)));
            });

        return typename marshalling_type::type(std::make_tuple(typename marshalling_type::index_type(matrix.rows()),
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
                entries.value().emplace_back(std::make_tuple(typename marshalling_type::index_type(index),
                                                             typename marshalling_type::value_type(value)));
            });

        return typename marshalling_type::type(
            std::make_tuple(typename marshalling_type::index_type(vector.size()), std::move(entries)));
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
