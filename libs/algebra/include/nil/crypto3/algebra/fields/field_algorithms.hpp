//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// MIT License
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_ALGEBRA_FIELDS_FIELD_ALGORITHMS_HPP
#define CRYPTO3_ALGEBRA_FIELDS_FIELD_ALGORITHMS_HPP

#include <cassert>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#include <nil/crypto3/algebra/fields/field_order.hpp>
#include <nil/crypto3/algebra/type_traits.hpp>
#include <nil/crypto3/algebra/fields/detail/field_algorithms/alt_bn128_fp12.hpp>

namespace nil::crypto3::algebra::fields {

    template<FieldValue Value>
    Value two_primary_component(const Value &x,
                                const boost::multiprecision::cpp_int &odd_order,
                                std::size_t two_adicity) {
        return detail::field_algorithms<Value>::two_primary_component(x, odd_order, two_adicity);
    }

    template<FieldValue Value>
    Value odd_subgroup_sqrt(const Value &x, const boost::multiprecision::cpp_int &odd_order) {
        return detail::field_algorithms<Value>::odd_subgroup_sqrt(x, odd_order);
    }

    template<FieldValue Value>
    Value sqrt_known_square(const Value &x) {
        return detail::field_algorithms<Value>::sqrt_known_square(x);
    }

    template<FieldValue Value>
    Value primitive_two_power_root_of_unity(std::size_t two_adicity) {
        return detail::field_algorithms<Value>::primitive_two_power_root_of_unity(two_adicity);
    }

    template<FieldValue Value>
    std::vector<Value> roots_of_unity(const Value &primitive_root, std::size_t two_adicity) {
        if (two_adicity >= std::numeric_limits<std::size_t>::digits) {
            throw std::invalid_argument("two-adicity is too large");
        }

        const std::size_t count = std::size_t(1) << two_adicity;
        assert(primitive_root.pow(count) == Value::one());
        assert(two_adicity == 0 || primitive_root.pow(count >> 1) != Value::one());

        std::vector<Value> roots(count);
        roots[0] = Value::one();
        for (std::size_t i = 1; i < count; ++i) {
            roots[i] = roots[i - 1] * primitive_root;
        }
        return roots;
    }

    template<FieldValue Value>
    std::vector<Value> roots_of_unity(std::size_t two_adicity) {
        return roots_of_unity(primitive_two_power_root_of_unity<Value>(two_adicity), two_adicity);
    }

}    // namespace nil::crypto3::algebra::fields

#endif    // CRYPTO3_ALGEBRA_FIELDS_FIELD_ALGORITHMS_HPP
