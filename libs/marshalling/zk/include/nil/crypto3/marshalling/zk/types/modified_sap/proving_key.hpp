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

#ifndef CRYPTO3_MARSHALLING_MODIFIED_SAP_PROVING_KEY_HPP
#define CRYPTO3_MARSHALLING_MODIFIED_SAP_PROVING_KEY_HPP

#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/options.hpp>
#include <nil/marshalling/processing/tuple.hpp>
#include <nil/marshalling/status_type.hpp>
#include <nil/marshalling/types/array_list.hpp>
#include <nil/marshalling/types/bundle.hpp>
#include <nil/marshalling/types/integral.hpp>

#include <nil/crypto3/algebra/curves/detail/scalar_mul.hpp>
#include <nil/crypto3/marshalling/algebra/types/curve_element.hpp>
#include <nil/crypto3/marshalling/algebra/types/detail/checked_field_element.hpp>
#include <nil/crypto3/marshalling/zk/types/modified_sap/constraint_system.hpp>
#include <nil/crypto3/marshalling/zk/types/modified_sap/verification_key.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/prover.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/proving_key.hpp>

namespace nil::crypto3::marshalling::types {

    namespace detail {

        template<typename TTypeBase, typename ProvingKey>
        struct modified_sap_proving_key_reader {
            template<typename MarshalledKey, typename InputIterator>
            nil::marshalling::status_type operator()(MarshalledKey &key, InputIterator &input,
                                                     std::size_t available_size) const {
                using status_type = nil::marshalling::status_type;
                using size_type = nil::marshalling::types::integral<TTypeBase, std::size_t>;
                using point_type = curve_element<TTypeBase, typename ProvingKey::g1_type>;
                using system_type = std::tuple_element_t<6, typename MarshalledKey::value_type>;
                using verification_key_type = std::tuple_element_t<7, typename MarshalledKey::value_type>;
                constexpr std::size_t prefix_size = size_type::max_length();
                constexpr std::size_t verification_key_size = verification_key_type::max_length();
                constexpr std::size_t tail_min_size = system_type::min_length() + verification_key_size;

                status_type status = status_type::success;
                std::size_t remaining_prefixes = 6;
                nil::marshalling::processing::tuple_for_each_until<6>(key.value(), [&](auto &query) {
                    if (status != status_type::success) {
                        return;
                    }
                    size_type count;
                    status = read_marshaled_value(count, input, available_size);
                    if (status != status_type::success) {
                        return;
                    }
                    --remaining_prefixes;
                    auto &points = query.value();
                    // Both the marshalled list and its eventual native query must be able to hold the count.
                    if (count.value() > points.max_size() ||
                        count.value() > std::vector<typename ProvingKey::g1_value_type>().max_size()) {
                        status = status_type::invalid_msg_data;
                        return;
                    }
                    // Leave room for later query prefixes and both nested objects before allocating points.
                    // The prefix count is at most five; division avoids multiplying an untrusted point count.
                    const std::size_t reserved_bytes = remaining_prefixes * prefix_size + tail_min_size;
                    if (available_size < reserved_bytes ||
                        count.value() > (available_size - reserved_bytes) / point_type::max_length()) {
                        status = status_type::not_enough_data;
                        return;
                    }
                    points.clear();
                    for (std::size_t i = 0; i < count.value(); ++i) {
                        point_type point;
                        // The existing point reader rejects noncanonical encodings and invalid subgroup points.
                        status = read_marshaled_value(point, input, available_size);
                        if (status != status_type::success) {
                            return;
                        }
                        points.push_back(std::move(point));
                    }
                });
                if (status != status_type::success) {
                    return status;
                }

                // Bound the variable-size system reader so its counts cannot consume the fixed verification key.
                auto &system = std::get<6>(key.value());
                status = system.read(input, available_size - verification_key_size);
                if (status != status_type::success) {
                    return status;
                }
                // A successful read consumed exactly length() bytes, also for sequential iterators.
                available_size -= system.length();
                return read_marshaled_value(std::get<7>(key.value()), input, available_size);
            }
        };

        struct modified_sap_proving_key_queries_validator {
            template<typename MarshalledKey>
            bool operator()(const MarshalledKey &key) const {
                bool valid = true;
                nil::marshalling::processing::tuple_for_each_until<6>(key.value(), [&](const auto &query) {
                    for (const auto &point : query.value()) {
                        // fill/make also accept assembled objects that have not passed through the byte reader.
                        if (!point.value().is_well_formed() ||
                            !algebra::curves::detail::subgroup_check(point.value())) {
                            valid = false;
                            return;
                        }
                    }
                });
                return valid;
            }
        };

    }    // namespace detail

    // W, H_query, T_A, T_C, T_H, T_Z, logical constraint system, and the nested verification key.
    template<typename TTypeBase, typename ProvingKey>
    using modified_sap_proving_key = nil::marshalling::types::bundle<
        TTypeBase,
        std::tuple<nil::marshalling::types::standard_array_list<TTypeBase,
                                                                curve_element<TTypeBase, typename ProvingKey::g1_type>>,
                   nil::marshalling::types::standard_array_list<TTypeBase,
                                                                curve_element<TTypeBase, typename ProvingKey::g1_type>>,
                   nil::marshalling::types::standard_array_list<TTypeBase,
                                                                curve_element<TTypeBase, typename ProvingKey::g1_type>>,
                   nil::marshalling::types::standard_array_list<TTypeBase,
                                                                curve_element<TTypeBase, typename ProvingKey::g1_type>>,
                   nil::marshalling::types::standard_array_list<TTypeBase,
                                                                curve_element<TTypeBase, typename ProvingKey::g1_type>>,
                   nil::marshalling::types::standard_array_list<TTypeBase,
                                                                curve_element<TTypeBase, typename ProvingKey::g1_type>>,
                   modified_sap_constraint_system<TTypeBase, typename ProvingKey::constraint_system_type>,
                   modified_sap_verification_key<TTypeBase, typename ProvingKey::verification_key_type>>,
        nil::marshalling::option::custom_value_reader<detail::modified_sap_proving_key_reader<TTypeBase, ProvingKey>>,
        nil::marshalling::option::contents_validator<detail::modified_sap_proving_key_queries_validator>>;

    template<typename Policy, typename Endianness>
    modified_sap_proving_key<nil::marshalling::field_type<Endianness>, typename Policy::proving_key_type>
        fill_modified_sap_proving_key(const typename Policy::proving_key_type &key) {
        // Reuse the prover's circuit, dimension, query-length, and circuit-digest checks.
        zk::snark::modified_sap_prover<Policy>::validate(key);
        using type_base = nil::marshalling::field_type<Endianness>;
        using key_type = typename Policy::proving_key_type;
        using g1_type = typename key_type::g1_type;
        modified_sap_proving_key<type_base, key_type> result(
            std::make_tuple(fill_curve_element_vector<g1_type, Endianness>(key.W),
                            fill_curve_element_vector<g1_type, Endianness>(key.H_query),
                            fill_curve_element_vector<g1_type, Endianness>(key.T_A),
                            fill_curve_element_vector<g1_type, Endianness>(key.T_C),
                            fill_curve_element_vector<g1_type, Endianness>(key.T_H),
                            fill_curve_element_vector<g1_type, Endianness>(key.T_Z),
                            fill_modified_sap_constraint_system<typename key_type::constraint_system_type, Endianness>(
                                key.constraint_system),
                            fill_modified_sap_verification_key<Policy, Endianness>(key.verification_key)));
        // The nested converters validate their objects; this also checks all six queries' G1 points.
        if (!result.valid()) {
            throw std::invalid_argument("invalid modified SAP proving key representation");
        }
        return result;
    }

    template<typename Policy, typename Endianness>
    typename Policy::proving_key_type
        make_modified_sap_proving_key(const modified_sap_proving_key<nil::marshalling::field_type<Endianness>,
                                                                     typename Policy::proving_key_type> &filled_key) {
        // Check assembled points and raw nested fields before conversion can reduce their integers.
        if (!filled_key.valid()) {
            throw std::invalid_argument("invalid modified SAP proving key representation");
        }
        using key_type = typename Policy::proving_key_type;
        using g1_type = typename key_type::g1_type;
        const auto &members = filled_key.value();
        key_type key;
        key.constraint_system =
            make_modified_sap_constraint_system<typename key_type::constraint_system_type, Endianness>(
                std::get<6>(members));
        key.verification_key = make_modified_sap_verification_key<Policy, Endianness>(std::get<7>(members));
        key.W = make_curve_element_vector<g1_type, Endianness>(std::get<0>(members));
        key.H_query = make_curve_element_vector<g1_type, Endianness>(std::get<1>(members));
        key.T_A = make_curve_element_vector<g1_type, Endianness>(std::get<2>(members));
        key.T_C = make_curve_element_vector<g1_type, Endianness>(std::get<3>(members));
        key.T_H = make_curve_element_vector<g1_type, Endianness>(std::get<4>(members));
        key.T_Z = make_curve_element_vector<g1_type, Endianness>(std::get<5>(members));
        // Publish only a key whose query lengths, dimensions and digest match its decoded circuit.
        zk::snark::modified_sap_prover<Policy>::validate(key);
        return key;
    }

}    // namespace nil::crypto3::marshalling::types

#endif    // CRYPTO3_MARSHALLING_MODIFIED_SAP_PROVING_KEY_HPP
