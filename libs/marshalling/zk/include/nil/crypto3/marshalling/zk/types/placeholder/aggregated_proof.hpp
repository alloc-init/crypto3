//---------------------------------------------------------------------------//
// Copyright (c) 2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2021 Nikita Kaskov <nbering@nil.foundation>
// Copyright (c) 2021 Ilias Khairullin <ilias@nil.foundation>
// Copyright (c) 2022-2023 Elena Tatuzova <e.tatuzova@nil.foundation>
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

#ifndef CRYPTO3_MARSHALLING_PLACEHOLDER_AGGREGATED_PROOF_HPP
#define CRYPTO3_MARSHALLING_PLACEHOLDER_AGGREGATED_PROOF_HPP

#include <nil/crypto3/marshalling/zk/types/commitments/lpc.hpp>
#include <nil/crypto3/marshalling/zk/types/placeholder/proof.hpp>

namespace nil {
    namespace crypto3 {
        namespace marshalling {
            namespace types {
                template<typename TTypeBase, typename Proof>
                using placeholder_aggregated_proof_type = nil::marshalling::types::bundle<
                    TTypeBase,
                    std::tuple<
                        nil::marshalling::types::
                            standard_array_list<TTypeBase, placeholder_partial_evaluation_proof<TTypeBase, Proof>>,
                        nil::crypto3::marshalling::types::aggregated_proof<TTypeBase,
                                                                           typename Proof::commitment_scheme_type>>>;

                template<typename Endianness, typename AggregatedProof, typename Proof>
                placeholder_aggregated_proof_type<nil::marshalling::field_type<Endianness>, Proof>
                    fill_placeholder_aggregated_proof(const AggregatedProof &proof) {
                    using TTypeBase = nil::marshalling::field_type<Endianness>;

                    nil::marshalling::types::standard_array_list<TTypeBase,
                                                                 placeholder_partial_evaluation_proof<TTypeBase, Proof>>
                        filled_partial_proofs;
                    for (const auto &it : proof.partial_proofs) {
                        filled_partial_proofs.value().push_back(
                            fill_placeholder_partial_evaluation_proof<Endianness, Proof>(it));
                    }

                    return placeholder_aggregated_proof_type<TTypeBase, Proof>(
                        std::make_tuple(filled_partial_proofs,
                                        fill_aggregated_proof<Endianness, typename Proof::commitment_scheme_type>(
                                            proof.aggregated_proof)));
                }

                template<typename Endianness, typename AggregatedProof, typename Proof>
                AggregatedProof make_placeholder_aggregated_proof(
                    const placeholder_aggregated_proof_type<nil::marshalling::field_type<Endianness>, Proof>
                        &filled_proof) {
                    AggregatedProof proof;

                    auto filled_partial_proofs = std::get<0>(filled_proof.value()).value();
                    for (const auto &it : filled_partial_proofs) {
                        proof.partial_proofs.push_back(
                            make_placeholder_partial_evaluation_proof<Endianness, Proof>(it));
                    }

                    proof.aggregated_proof =
                        make_aggregated_proof<Endianness, typename AggregatedProof::commitment_scheme_type>(
                            std::get<1>(filled_proof.value()));

                    return proof;
                }
            }    // namespace types
        }    // namespace marshalling
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MARSHALLING_PLACEHOLDER_AGGREGATED_PROOF_HPP
