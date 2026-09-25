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

#ifndef CRYPTO3_MARSHALLING_MODIFIED_SAP_VERIFICATION_KEY_HPP
#define CRYPTO3_MARSHALLING_MODIFIED_SAP_VERIFICATION_KEY_HPP

#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/options.hpp>
#include <nil/marshalling/processing/tuple.hpp>
#include <nil/marshalling/status_type.hpp>
#include <nil/marshalling/types/array_list.hpp>
#include <nil/marshalling/types/bundle.hpp>
#include <nil/marshalling/types/integral.hpp>

#include <nil/crypto3/marshalling/algebra/types/curve_element.hpp>
#include <nil/crypto3/marshalling/algebra/types/detail/checked_field_element.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/verification_key.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/verifier.hpp>

namespace nil::crypto3::marshalling::types {

    namespace detail {

        template<typename TTypeBase, typename VerificationKey>
        struct modified_sap_verification_key_reader {
            template<typename MarshalledKey, typename InputIterator>
            nil::marshalling::status_type operator()(MarshalledKey &key, InputIterator &input,
                                                     std::size_t available_size) const {
                using status_type = nil::marshalling::status_type;
                using gt_type = validated_field_element<TTypeBase, typename VerificationKey::gt_value_type>;
                using digest_type = validated_field_element<TTypeBase, typename VerificationKey::base_value_type>;
                using alpha_array_type = std::tuple_element_t<4, typename MarshalledKey::value_type>;
                if (available_size < MarshalledKey::max_length()) {
                    return status_type::not_enough_data;
                }

                // A GT field's generic reader uses the supplied length. Bound each read to its twelve components
                // so it cannot consume the next GT value, dimensions, or digest.
                auto read_gt = [&](auto &member) {
                    return read_marshaled_value(
                        member, input, available_size,
                        canonical_field_element_encoding_validator<TTypeBase,
                                                                   typename VerificationKey::gt_value_type>());
                };
                status_type status = status_type::success;
                nil::marshalling::processing::tuple_for_each(key.value(), [&](auto &member) {
                    if (status != status_type::success) {
                        return;
                    }
                    using member_type = std::decay_t<decltype(member)>;
                    if constexpr (std::is_same_v<member_type, gt_type>) {
                        status = read_gt(member);
                    } else if constexpr (std::is_same_v<member_type, alpha_array_type>) {
                        // The four GT entries are fixed by the format; there is no count prefix to decode.
                        member.value().resize(4);
                        for (auto &alpha : member.value()) {
                            status = read_gt(alpha);
                            if (status != status_type::success) {
                                break;
                            }
                        }
                    } else if constexpr (std::is_same_v<member_type, digest_type>) {
                        // The digest is an Fp value, with the same canonical-component checks as GT.
                        status = read_marshaled_value(
                            member, input, available_size,
                            canonical_field_element_encoding_validator<TTypeBase,
                                                                       typename VerificationKey::base_value_type>());
                    } else {
                        status = read_marshaled_value(member, input, available_size);
                    }
                });
                return status;
            }
        };

        struct modified_sap_verification_key_shape_validator {
            template<typename MarshalledKey>
            bool operator()(const MarshalledKey &key) const {
                // sequence_fixed_size also permits shorter native lists; received keys require all four entries.
                return std::get<4>(key.value()).value().size() == 4;
            }
        };

    }    // namespace detail

    // The reader checks canonical encodings. make_modified_sap_verification_key also applies native key invariants.
    template<typename TTypeBase, typename VerificationKey>
    using modified_sap_verification_key = nil::marshalling::types::bundle<
        TTypeBase,
        std::tuple<curve_element<TTypeBase, typename VerificationKey::g2_type>,
                   curve_element<TTypeBase, typename VerificationKey::g2_type>,
                   curve_element<TTypeBase, typename VerificationKey::g2_type>,
                   detail::validated_field_element<TTypeBase, typename VerificationKey::gt_value_type>,
                   nil::marshalling::types::array_list<
                       TTypeBase, detail::validated_field_element<TTypeBase, typename VerificationKey::gt_value_type>,
                       nil::marshalling::option::sequence_fixed_size<4>,
                       nil::marshalling::option::sequence_fixed_size_use_fixed_size_storage>,
                   nil::marshalling::types::integral<TTypeBase, std::size_t>,
                   nil::marshalling::types::integral<TTypeBase, std::size_t>,
                   detail::validated_field_element<TTypeBase, typename VerificationKey::base_value_type>>,
        nil::marshalling::option::custom_value_reader<
            detail::modified_sap_verification_key_reader<TTypeBase, VerificationKey>>,
        nil::marshalling::option::contents_validator<detail::modified_sap_verification_key_shape_validator>>;

    template<typename Policy, typename Endianness>
    modified_sap_verification_key<nil::marshalling::field_type<Endianness>, typename Policy::verification_key_type>
        fill_modified_sap_verification_key(const typename Policy::verification_key_type &key) {
        if (!zk::snark::modified_sap_verifier<Policy>::validate_verification_key(key)) {
            throw std::invalid_argument("invalid modified SAP verification key");
        }
        using type_base = nil::marshalling::field_type<Endianness>;
        using key_type = typename Policy::verification_key_type;
        using point_type = curve_element<type_base, typename key_type::g2_type>;
        using gt_type = detail::validated_field_element<type_base, typename key_type::gt_value_type>;
        using size_type = nil::marshalling::types::integral<type_base, std::size_t>;
        using digest_type = detail::validated_field_element<type_base, typename key_type::base_value_type>;
        using result_type = modified_sap_verification_key<type_base, key_type>;
        std::tuple_element_t<4, typename result_type::value_type> alphas;
        for (const auto &alpha : key.alpha_gt) {
            alphas.value().emplace_back(alpha);
        }
        return result_type(std::make_tuple(point_type(key.g2_one), point_type(key.tau_g2),
                                           point_type(key.gamma_inverse_g2), gt_type(key.alpha_z_vanishing_gt),
                                           std::move(alphas), size_type(key.num_variables), size_type(key.domain_size),
                                           digest_type(key.circuit_digest)));
    }

    template<typename Policy, typename Endianness>
    typename Policy::verification_key_type make_modified_sap_verification_key(
        const modified_sap_verification_key<nil::marshalling::field_type<Endianness>,
                                            typename Policy::verification_key_type> &filled_key) {
        // Validate all raw components and fixed counts before value() can reduce integers or index missing entries.
        if (!filled_key.valid()) {
            throw std::invalid_argument("invalid modified SAP verification key representation");
        }
        const auto &members = filled_key.value();
        typename Policy::verification_key_type key;
        key.g2_one = std::get<0>(members).value();
        key.tau_g2 = std::get<1>(members).value();
        key.gamma_inverse_g2 = std::get<2>(members).value();
        key.alpha_z_vanishing_gt = std::get<3>(members).value();
        for (std::size_t i = 0; i < key.alpha_gt.size(); ++i) {
            key.alpha_gt[i] = std::get<4>(members).value()[i].value();
        }
        key.num_variables = std::get<5>(members).value();
        key.domain_size = std::get<6>(members).value();
        key.circuit_digest = std::get<7>(members).value();
        if (!zk::snark::modified_sap_verifier<Policy>::validate_verification_key(key)) {
            throw std::invalid_argument("invalid modified SAP verification key");
        }
        // The circuit is absent from a standalone key, so its digest cannot be recomputed here.
        return key;
    }

}    // namespace nil::crypto3::marshalling::types

#endif    // CRYPTO3_MARSHALLING_MODIFIED_SAP_VERIFICATION_KEY_HPP
