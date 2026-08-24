//---------------------------------------------------------------------------//
// Copyright (c) 2026 Alloc Init Labs Inc.
//
// MIT License
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_MARSHALLING_PROCESSING_FIELD_ELEMENT_HPP
#define CRYPTO3_MARSHALLING_PROCESSING_FIELD_ELEMENT_HPP

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <nil/crypto3/marshalling/algebra/types/field_element.hpp>
#include <nil/marshalling/algorithms/pack.hpp>

namespace nil::crypto3::marshalling::processing {

    template<algebra::FieldValue FieldValue>
    constexpr std::size_t field_element_bit_length() {
        if constexpr (requires { FieldValue::bitsize(); }) {
            return FieldValue::bitsize();
        } else {
            constexpr std::size_t value_bits = FieldValue::field_type::value_bits;
            return ((value_bits + 7) / 8) * 8;
        }
    }

    template<algebra::FieldValue FieldValue>
    std::vector<std::uint8_t> field_element_to_le_bytes(const FieldValue &value) {
        if constexpr (requires { value.to_le_bytes(); }) {
            return value.to_le_bytes();
        } else {
            nil::marshalling::status_type status;
            std::vector<std::uint8_t> bytes =
                nil::marshalling::pack<nil::marshalling::option::little_endian>(value, status);
            if (status != nil::marshalling::status_type::success) {
                throw std::runtime_error("error encoding field element");
            }
            return bytes;
        }
    }

    template<algebra::FieldValue FieldValue>
    FieldValue field_element_from_le_bytes(const std::vector<std::uint8_t> &bytes) {
        if constexpr (requires { FieldValue::from_le_bytes(bytes); }) {
            return FieldValue::from_le_bytes(bytes);
        } else {
            nil::marshalling::status_type status;
            FieldValue value = nil::marshalling::pack<nil::marshalling::option::little_endian>(bytes, status);
            if (status != nil::marshalling::status_type::success) {
                throw std::runtime_error("error decoding field element");
            }
            return value;
        }
    }

    template<algebra::FieldValue FieldValue>
    std::vector<bool> field_element_to_le_bits(const FieldValue &value) {
        std::vector<bool> bits(field_element_bit_length<FieldValue>());
        std::size_t current = 0;
        for (std::uint8_t byte : field_element_to_le_bytes(value)) {
            for (std::size_t i = 0; i < 8 && current < bits.size(); ++i) {
                bits[current++] = (byte >> i) & 1;
            }
        }
        return bits;
    }

}    // namespace nil::crypto3::marshalling::processing

#endif    // CRYPTO3_MARSHALLING_PROCESSING_FIELD_ELEMENT_HPP
