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
#include <type_traits>

#include <nil/marshalling/endianness.hpp>
#include <nil/marshalling/options.hpp>
#include <nil/marshalling/status_type.hpp>

#include <nil/crypto3/marshalling/algebra/types/field_element.hpp>

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

        /**
         * Check raw components through valid(). For checked decoding, use read_marshaled_value with
         * canonical_field_element_encoding_validator to reject unused high bits as well.
         */
        template<typename TTypeBase, typename FieldValueType>
        using validated_field_element = field_element<
            TTypeBase, FieldValueType,
            nil::marshalling::option::contents_validator<canonical_field_element_validator<FieldValueType>>>;

        struct accept_any_marshaled_encoding {
            template<typename InputIterator>
            constexpr bool operator()(const InputIterator &) const {
                return true;
            }
        };

        /** Reject nonzero high padding bits before a fixed-width integer decoder can discard them. */
        template<typename TTypeBase, typename FieldValueType>
        struct canonical_field_element_encoding_validator {
            template<typename RandomAccessIterator>
            bool operator()(RandomAccessIterator input) const {
                constexpr std::size_t component_bits = FieldValueType::field_type::modulus_bits;
                constexpr std::size_t component_length = (component_bits + 7) / 8;
                constexpr std::size_t unused_high_bits = component_length * 8 - component_bits;
                if constexpr (unused_high_bits == 0) {
                    return true;
                } else {
                    constexpr std::uint8_t unused_high_bits_mask =
                        static_cast<std::uint8_t>(0xffU << (8 - unused_high_bits));
                    for (std::size_t component = 0; component < FieldValueType::field_type::arity; ++component) {
                        std::size_t high_byte_offset = component * component_length;
                        if constexpr (std::is_same_v<typename TTypeBase::endian_type,
                                                     nil::marshalling::endian::little_endian>) {
                            high_byte_offset += component_length - 1;
                        }
                        if ((static_cast<std::uint8_t>(input[high_byte_offset]) & unused_high_bits_mask) != 0) {
                            return false;
                        }
                    }
                    return true;
                }
            }
        };

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
        template<typename MarshalledType, typename InputIterator,
                 typename EncodingValidator = accept_any_marshaled_encoding>
        nil::marshalling::status_type read_marshaled_value(MarshalledType &value, InputIterator &input,
                                                           std::size_t &available_size,
                                                           EncodingValidator encoding_validator = {}) {
            constexpr std::size_t length = MarshalledType::max_length();
            if (available_size < length) {
                return nil::marshalling::status_type::not_enough_data;
            }

            using iterator_category = typename std::iterator_traits<std::decay_t<InputIterator>>::iterator_category;
            nil::marshalling::status_type read_status;
            bool encoding_is_valid = true;
            if constexpr (std::is_base_of_v<std::random_access_iterator_tag, iterator_category>) {
                const InputIterator encoded_value = input;
                read_status = value.read(input, length);
                encoding_is_valid = encoding_validator(encoded_value);
            } else {
                std::array<std::uint8_t, length> scratch;
                for (std::size_t index = 0; index < length; ++index) {
                    scratch[index] = static_cast<std::uint8_t>(*input);
                    ++input;
                }
                auto scratch_input = scratch.begin();
                read_status = value.read(scratch_input, length);
                encoding_is_valid = encoding_validator(scratch.begin());
            }
            if (read_status != nil::marshalling::status_type::success) {
                return read_status;
            }
            available_size -= length;
            if (!encoding_is_valid || !value.valid() || value.length() != length) {
                return nil::marshalling::status_type::invalid_msg_data;
            }
            return nil::marshalling::status_type::success;
        }

    }    // namespace detail

}    // namespace nil::crypto3::marshalling::types
