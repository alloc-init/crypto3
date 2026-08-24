//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// MIT License
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_ALGEBRA_FIELDS_DETAIL_FIELD_ALGORITHMS_HPP
#define CRYPTO3_ALGEBRA_FIELDS_DETAIL_FIELD_ALGORITHMS_HPP

#include <cassert>
#include <cstddef>
#include <stdexcept>

#include <boost/multiprecision/cpp_int.hpp>

#include <nil/crypto3/algebra/type_traits.hpp>

namespace nil::crypto3::algebra::fields::detail {

    // Customization point used by the public free functions in field_algorithms.hpp.
    // This primary template provides formulas that work for any FieldValue. A field
    // may specialize this class in a separate detail header when its tower structure
    // permits a faster implementation. Keeping dispatch here gives callers one stable
    // public API while keeping field-specific types and algorithms private.
    template<FieldValue Value>
    struct field_algorithms {
        static Value two_primary_component(const Value &x,
                                           const boost::multiprecision::cpp_int &odd_order,
                                           std::size_t) {
            return x.pow(odd_order);
        }

        static Value odd_subgroup_sqrt(const Value &x, const boost::multiprecision::cpp_int &odd_order) {
            assert((odd_order & 1) != 0);
            return x.is_zero() || x.is_one() ? x : x.pow((odd_order + 1) >> 1);
        }

        static Value sqrt_known_square(const Value &x) {
            return x.sqrt();
        }

        static Value primitive_two_power_root_of_unity(std::size_t) {
            throw std::invalid_argument("primitive root discovery is not available for this field");
        }
    };

}    // namespace nil::crypto3::algebra::fields::detail

#endif    // CRYPTO3_ALGEBRA_FIELDS_DETAIL_FIELD_ALGORITHMS_HPP
