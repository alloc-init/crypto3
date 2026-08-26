//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_ALGEBRA_FIELDS_FIELD_ORDER_HPP
#define CRYPTO3_ALGEBRA_FIELDS_FIELD_ORDER_HPP

#include <cstddef>
#include <stdexcept>
#include <utility>

#include <boost/multiprecision/cpp_int.hpp>

#include <nil/crypto3/algebra/type_traits.hpp>

namespace nil::crypto3::algebra::fields {

    struct multiplicative_group_decomposition {
        boost::multiprecision::cpp_int odd_order;
        std::size_t two_adicity = 0;
    };

    /** Decompose a positive multiplicative-group order into odd_order * 2^two_adicity. */
    inline multiplicative_group_decomposition
        decompose_multiplicative_group_order(boost::multiprecision::cpp_int group_order) {
        if (group_order <= 0) {
            throw std::invalid_argument("a multiplicative-group order must be positive");
        }
        std::size_t two_adicity = 0;
        while ((group_order & 1) == 0) {
            group_order >>= 1;
            ++two_adicity;
        }
        return {std::move(group_order), two_adicity};
    }

    template<typename FieldOrValue, bool = FieldValue<FieldOrValue>>
    struct field_type {
        using type = FieldOrValue;
    };

    template<typename FieldValueType>
    struct field_type<FieldValueType, true> {
        using type = typename FieldValueType::field_type;
    };

    template<typename FieldOrValue>
    using field_type_t = typename field_type<FieldOrValue>::type;

    /** Return the characteristic of FieldType. Extension fields have the same characteristic as their prime field. */
    template<typename FieldType>
    boost::multiprecision::cpp_int field_characteristic() {
        using field_type = field_type_t<FieldType>;
        return boost::multiprecision::cpp_int(field_type::modulus.backend().to_cpp_int());
    }

    /**
     * Return the number of elements in FieldType. Crypto3 field types describe extensions with an arity equal to their
     * degree over the prime field, so a field of characteristic p and arity d has p^d elements.
     */
    template<typename FieldType>
    boost::multiprecision::cpp_int field_order() {
        using field_type = field_type_t<FieldType>;
        const boost::multiprecision::cpp_int characteristic = field_characteristic<field_type>();
        boost::multiprecision::cpp_int order = 1;
        for (std::size_t i = 0; i < field_type::arity; ++i) {
            order *= characteristic;
        }
        return order;
    }

    /** Return the number of elements in a positive-degree extension of FieldType. */
    template<typename FieldType>
    boost::multiprecision::cpp_int extension_field_order(std::size_t extension_degree) {
        if (extension_degree == 0) {
            throw std::invalid_argument("a field extension must have positive degree");
        }
        const boost::multiprecision::cpp_int base_field_order = field_order<FieldType>();
        boost::multiprecision::cpp_int order = 1;
        for (std::size_t i = 0; i < extension_degree; ++i) {
            order *= base_field_order;
        }
        return order;
    }

    /** Decompose the multiplicative-group order of FieldType into odd_order * 2^two_adicity. */
    template<typename FieldType>
    multiplicative_group_decomposition field_multiplicative_group_decomposition() {
        return decompose_multiplicative_group_order(field_order<FieldType>() - 1);
    }

    /** Decompose the multiplicative-group order of a positive-degree extension of FieldType. */
    template<typename FieldType>
    multiplicative_group_decomposition
        extension_field_multiplicative_group_decomposition(std::size_t extension_degree) {
        return decompose_multiplicative_group_order(extension_field_order<FieldType>(extension_degree) - 1);
    }

}    // namespace nil::crypto3::algebra::fields

#endif    // CRYPTO3_ALGEBRA_FIELDS_FIELD_ORDER_HPP
