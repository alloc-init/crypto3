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

#include <boost/multiprecision/cpp_int.hpp>

namespace nil::crypto3::algebra::fields {

    /** Return the characteristic of FieldType. Extension fields have the same characteristic as their prime field. */
    template<typename FieldType>
    boost::multiprecision::cpp_int field_characteristic() {
        return boost::multiprecision::cpp_int(FieldType::modulus.backend().to_cpp_int());
    }

    /**
     * Return the number of elements in FieldType. Crypto3 field types describe extensions with an arity equal to their
     * degree over the prime field, so a field of characteristic p and arity d has p^d elements.
     */
    template<typename FieldType>
    boost::multiprecision::cpp_int field_order() {
        const boost::multiprecision::cpp_int characteristic = field_characteristic<FieldType>();
        boost::multiprecision::cpp_int order = 1;
        for (std::size_t i = 0; i < FieldType::arity; ++i) {
            order *= characteristic;
        }
        return order;
    }

}    // namespace nil::crypto3::algebra::fields

#endif    // CRYPTO3_ALGEBRA_FIELDS_FIELD_ORDER_HPP
