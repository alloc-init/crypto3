//---------------------------------------------------------------------------//
// Copyright (c) 2026 Alloc Init
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

#ifndef CRYPTO3_MARSHALLING_MODIFIED_SAP_CONSTRAINT_SYSTEM_HPP
#define CRYPTO3_MARSHALLING_MODIFIED_SAP_CONSTRAINT_SYSTEM_HPP

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/options.hpp>
#include <nil/marshalling/status_type.hpp>
#include <nil/marshalling/types/array_list.hpp>
#include <nil/marshalling/types/bundle.hpp>
#include <nil/marshalling/types/integral.hpp>

#include <nil/crypto3/marshalling/algebra/types/detail/checked_field_element.hpp>
#include <nil/crypto3/marshalling/zk/types/r1cs_gg_ppzksnark/r1cs.hpp>
#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/modified_sap.hpp>

namespace nil::crypto3::marshalling::types {

    namespace detail {

        template<typename TTypeBase, typename ConstraintSystem>
        struct modified_sap_constraint_system_reader {
            template<typename MarshalledSystem, typename InputIterator>
            nil::marshalling::status_type operator()(MarshalledSystem &system, InputIterator &input,
                                                     std::size_t available_size) const {
                using status_type = nil::marshalling::status_type;
                using size_type = nil::marshalling::types::integral<TTypeBase, std::size_t>;
                using field_value_type = typename ConstraintSystem::field_value_type;
                using native_combination_type = typename ConstraintSystem::linear_combination_type;
                // The first bundle field stores witness_size, the number of witness entries.
                auto status = read_marshaled_value(std::get<0>(system.value()), input, available_size);
                if (status != status_type::success) {
                    return status;
                }
                size_type row_count;
                status = read_marshaled_value(row_count, input, available_size);
                if (status != status_type::success) {
                    return status;
                }
                auto &rows = std::get<1>(system.value()).value();
                using row_type = typename std::decay_t<decltype(rows)>::value_type;
                using combination_type = std::tuple_element_t<0, typename row_type::value_type>;
                using term_type = typename combination_type::element_type;
                const std::size_t count = row_count.value();
                if (count > rows.max_size() ||
                    count > std::vector<typename ConstraintSystem::constraint_type>().max_size()) {
                    return status_type::invalid_msg_data;
                }
                // Every remaining row needs two list prefixes, even when both combinations are empty.
                // Division checks the bound without multiplying an untrusted count first.
                constexpr std::size_t prefix_size = size_type::max_length();
                if (count > available_size / (2 * prefix_size)) {
                    return status_type::not_enough_data;
                }
                std::size_t remaining_prefixes = 2 * count;    // Safe after the byte-budget check above.
                rows.clear();
                for (std::size_t i = 0; i < count; ++i) {
                    row_type row;
                    for (auto *combination : {&std::get<0>(row.value()), &std::get<1>(row.value())}) {
                        size_type term_count;
                        status = read_marshaled_value(term_count, input, available_size);
                        if (status != status_type::success) {
                            return status;
                        }
                        --remaining_prefixes;
                        auto &terms = combination->value();
                        if (term_count.value() > terms.max_size() ||
                            term_count.value() > native_combination_type().terms.max_size()) {
                            return status_type::invalid_msg_data;
                        }
                        // Terms must leave space for the other combination and all later rows' prefixes.
                        const std::size_t reserved_bytes = remaining_prefixes * prefix_size;
                        if (available_size < reserved_bytes ||
                            term_count.value() > (available_size - reserved_bytes) / term_type::max_length()) {
                            return status_type::not_enough_data;
                        }
                        for (std::size_t j = 0; j < term_count.value(); ++j) {
                            term_type term;
                            status = read_marshaled_value(std::get<0>(term.value()), input, available_size);
                            if (status != status_type::success) {
                                return status;
                            }
                            status = read_marshaled_value(
                                std::get<1>(term.value()), input, available_size,
                                canonical_field_element_encoding_validator<TTypeBase, field_value_type>());
                            if (status != status_type::success) {
                                return status;
                            }
                            // Grow only after checking the count and reading a complete, canonical term.
                            terms.push_back(std::move(term));
                        }
                    }
                    rows.push_back(std::move(row));
                }
                return status_type::success;
            }
        };

    }    // namespace detail

    // Each logical row stores a and c, using the existing explicitly indexed linear-combination codec.
    template<typename TTypeBase, typename Constraint>
    using modified_sap_constraint = nil::marshalling::types::bundle<
        TTypeBase,
        std::tuple<
            linear_combination<TTypeBase, typename Constraint::linear_combination_type, bool,
                               detail::validated_field_element<TTypeBase, typename Constraint::field_type::value_type>>,
            linear_combination<
                TTypeBase, typename Constraint::linear_combination_type, bool,
                detail::validated_field_element<TTypeBase, typename Constraint::field_type::value_type>>>>;

    // witness_size, then the size-prefixed list of logical rows. Domain padding is not serialized.
    template<typename TTypeBase, typename ConstraintSystem>
    using modified_sap_constraint_system = nil::marshalling::types::bundle<
        TTypeBase,
        std::tuple<nil::marshalling::types::integral<TTypeBase, std::size_t>,
                   nil::marshalling::types::standard_array_list<
                       TTypeBase, modified_sap_constraint<TTypeBase, typename ConstraintSystem::constraint_type>>>,
        nil::marshalling::option::custom_value_reader<
            detail::modified_sap_constraint_system_reader<TTypeBase, ConstraintSystem>>>;

    template<typename ConstraintSystem, typename Endianness>
    modified_sap_constraint_system<nil::marshalling::field_type<Endianness>, ConstraintSystem>
        fill_modified_sap_constraint_system(const ConstraintSystem &system) {
        // normalized() checks every raw index before combining or dropping terms, and owns its result.
        const auto normalized = system.normalized();
        using type_base = nil::marshalling::field_type<Endianness>;
        using row_type = typename ConstraintSystem::constraint_type;
        using combination_type = typename ConstraintSystem::linear_combination_type;
        using coefficient_type =
            detail::validated_field_element<type_base, typename ConstraintSystem::field_value_type>;
        using marshalled_row = modified_sap_constraint<type_base, row_type>;
        using size_type = nil::marshalling::types::integral<type_base, std::size_t>;
        auto rows = nil::marshalling::types::fill_standard_array_list<type_base, marshalled_row>(
            normalized.constraints, [](const row_type &row) {
                return marshalled_row(
                    std::make_tuple(fill_linear_combination<combination_type, Endianness, coefficient_type>(row.a),
                                    fill_linear_combination<combination_type, Endianness, coefficient_type>(row.c)));
            });
        return modified_sap_constraint_system<type_base, ConstraintSystem>(
            std::make_tuple(size_type(normalized.witness_size), std::move(rows)));
    }

    template<typename ConstraintSystem, typename Endianness>
    ConstraintSystem make_modified_sap_constraint_system(
        const modified_sap_constraint_system<nil::marshalling::field_type<Endianness>, ConstraintSystem>
            &filled_system) {
        // Checked coefficients must be validated before native field conversion can reduce their raw integers.
        if (!filled_system.valid()) {
            throw std::invalid_argument("invalid modified SAP constraint system representation");
        }
        using type_base = nil::marshalling::field_type<Endianness>;
        using row_type = typename ConstraintSystem::constraint_type;
        using combination_type = typename ConstraintSystem::linear_combination_type;
        using coefficient_type =
            detail::validated_field_element<type_base, typename ConstraintSystem::field_value_type>;
        using marshalled_row = modified_sap_constraint<type_base, row_type>;
        ConstraintSystem system;
        system.witness_size = std::get<0>(filled_system.value()).value();
        system.constraints = nil::marshalling::types::make_standard_array_list<type_base, row_type, marshalled_row>(
            std::get<1>(filled_system.value()), [](const marshalled_row &row) {
                return row_type {
                    make_linear_combination<combination_type, Endianness, coefficient_type>(std::get<0>(row.value())),
                    make_linear_combination<combination_type, Endianness, coefficient_type>(std::get<1>(row.value()))};
            });
        for (const auto &row : system.constraints) {
            for (const auto *combination : {&row.a, &row.c}) {
                // is_valid() checks strict index order and bounds; the encoding also excludes zero coefficients.
                if (!combination->is_valid(system.witness_size) ||
                    std::any_of(combination->terms.begin(), combination->terms.end(),
                                [](const auto &term) { return term.coeff.is_zero(); })) {
                    throw std::invalid_argument("noncanonical modified SAP linear combination");
                }
            }
        }
        if (!system.is_valid()) {
            throw std::invalid_argument("invalid modified SAP constraint system");
        }
        return system;
    }

}    // namespace nil::crypto3::marshalling::types

#endif    // CRYPTO3_MARSHALLING_MODIFIED_SAP_CONSTRAINT_SYSTEM_HPP
