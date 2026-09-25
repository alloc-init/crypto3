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

#ifndef CRYPTO3_MARSHALLING_MODIFIED_SAP_PROOF_HPP
#define CRYPTO3_MARSHALLING_MODIFIED_SAP_PROOF_HPP

#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <type_traits>

#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/options.hpp>
#include <nil/marshalling/processing/tuple.hpp>
#include <nil/marshalling/status_type.hpp>
#include <nil/marshalling/types/bundle.hpp>

#include <nil/crypto3/algebra/curves/detail/scalar_mul.hpp>
#include <nil/crypto3/marshalling/algebra/types/curve_element.hpp>
#include <nil/crypto3/marshalling/algebra/types/detail/checked_field_element.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/proof.hpp>

namespace nil::crypto3::marshalling::types {

    namespace detail {

        template<typename TTypeBase, typename Proof>
        struct modified_sap_proof_reader {
            template<typename MarshalledProof, typename InputIterator>
            nil::marshalling::status_type operator()(MarshalledProof &proof, InputIterator &input,
                                                     std::size_t available_size) const {
                using status_type = nil::marshalling::status_type;
                using scalar_type = validated_field_element<TTypeBase, typename Proof::scalar_value_type>;
                if (available_size < MarshalledProof::max_length()) {
                    return status_type::not_enough_data;
                }

                status_type status = status_type::success;
                nil::marshalling::processing::tuple_for_each(proof.value(), [&](auto &member) {
                    if (status != status_type::success) {
                        return;
                    }
                    if constexpr (std::is_same_v<std::decay_t<decltype(member)>, scalar_type>) {
                        // Check raw scalars and their padding before conversion can reduce or discard bits.
                        status = read_marshaled_value(
                            member, input, available_size,
                            canonical_field_element_encoding_validator<TTypeBase, typename Proof::scalar_value_type>());
                    } else {
                        // The point codec checks canonical coordinates, flags, and subgroup membership.
                        status = read_marshaled_value(member, input, available_size);
                    }
                });
                return status;
            }
        };

        struct modified_sap_proof_validator {
            template<typename MarshalledProof>
            bool operator()(const MarshalledProof &proof) const {
                // fill/make can receive directly assembled points without going through the byte reader.
                const auto &P = std::get<0>(proof.value()).value();
                const auto &Q = std::get<1>(proof.value()).value();
                return P.is_well_formed() && algebra::curves::detail::subgroup_check(P) && Q.is_well_formed() &&
                       algebra::curves::detail::subgroup_check(Q);
            }
        };

    }    // namespace detail

    // Ordered payload: P, Q, v_A, v_C, v_H, v_Z. The public input is supplied separately.
    template<typename TTypeBase, typename Proof>
    using modified_sap_proof = nil::marshalling::types::bundle<
        TTypeBase,
        std::tuple<curve_element<TTypeBase, typename Proof::g1_type>, curve_element<TTypeBase, typename Proof::g1_type>,
                   detail::validated_field_element<TTypeBase, typename Proof::scalar_value_type>,
                   detail::validated_field_element<TTypeBase, typename Proof::scalar_value_type>,
                   detail::validated_field_element<TTypeBase, typename Proof::scalar_value_type>,
                   detail::validated_field_element<TTypeBase, typename Proof::scalar_value_type>>,
        nil::marshalling::option::custom_value_reader<detail::modified_sap_proof_reader<TTypeBase, Proof>>,
        nil::marshalling::option::contents_validator<detail::modified_sap_proof_validator>>;

    template<typename Proof, typename Endianness>
    modified_sap_proof<nil::marshalling::field_type<Endianness>, Proof> fill_modified_sap_proof(const Proof &proof) {
        using type_base = nil::marshalling::field_type<Endianness>;
        using point_type = curve_element<type_base, typename Proof::g1_type>;
        using scalar_type = detail::validated_field_element<type_base, typename Proof::scalar_value_type>;
        modified_sap_proof<type_base, Proof> result(std::make_tuple(point_type(proof.P), point_type(proof.Q),
                                                                    scalar_type(proof.v_A), scalar_type(proof.v_C),
                                                                    scalar_type(proof.v_H), scalar_type(proof.v_Z)));
        if (!result.valid()) {
            throw std::invalid_argument("invalid modified SAP proof representation");
        }
        return result;
    }

    template<typename Proof, typename Endianness>
    Proof make_modified_sap_proof(
        const modified_sap_proof<nil::marshalling::field_type<Endianness>, Proof> &filled_proof) {
        // Validate raw scalar storage before value() converts it to the native field type.
        if (!filled_proof.valid()) {
            throw std::invalid_argument("invalid modified SAP proof representation");
        }
        const auto &members = filled_proof.value();
        return {std::get<0>(members).value(), std::get<1>(members).value(), std::get<2>(members).value(),
                std::get<3>(members).value(), std::get<4>(members).value(), std::get<5>(members).value()};
    }

}    // namespace nil::crypto3::marshalling::types

#endif    // CRYPTO3_MARSHALLING_MODIFIED_SAP_PROOF_HPP
