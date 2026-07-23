//---------------------------------------------------------------------------//
// Copyright (c) 2020 Ilias Khairullin <ilias@nil.foundation>
// Copyright (c) 2020 Mikhail Komarov <nemo@nil.foundation>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON_HPP
#define CRYPTO3_HASH_POSEIDON_HPP

#include <nil/crypto3/hash/accumulators/hash.hpp>
#include <nil/crypto3/hash/detail/poseidon_common/poseidon_sponge.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_optimized_permutation.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_permutation.hpp>
#include <nil/crypto3/hash/detail/poseidon2/poseidon2_permutation.hpp>
#include <nil/crypto3/hash/detail/stream_processors/stream_processors_enum.hpp>
namespace nil {
    namespace crypto3 {
        namespace hashes {
            template<typename PolicyType, typename PermutationType, detail::poseidon_sponge_padding_mode PaddingMode>
            class basic_poseidon1 {
            public:
                using policy_type = PolicyType;
                using permutation_type = PermutationType;

                using word_type = typename policy_type::word_type;
                constexpr static const std::size_t word_bits = policy_type::word_bits;

                constexpr static const std::size_t block_words = policy_type::block_words;
                using block_type = typename policy_type::block_type;

                constexpr static const std::size_t digest_bits = policy_type::digest_bits;
                using digest_type = typename policy_type::digest_type;

                struct construction {
                    struct params_type {
                        // This is required by the hash concept.
                    };

                    using type = detail::poseidon_sponge_construction<
                        policy_type, permutation_type, detail::poseidon_sponge_absorb_mode::overwrite, PaddingMode>;
                };

                constexpr static detail::stream_processor_type stream_processor = detail::stream_processor_type::raw;
                using accumulator_tag =
                    accumulators::tag::algebraic_hash<basic_poseidon1<policy_type, permutation_type, PaddingMode>>;
            };

            template<typename PolicyType>
            class poseidon1 : public basic_poseidon1<PolicyType, detail::poseidon1_optimized_permutation<PolicyType>,
                                                     detail::poseidon_sponge_padding_mode::pad10> {
            public:
                using policy_type = PolicyType;
                using accumulator_tag = accumulators::tag::algebraic_hash<poseidon1<policy_type>>;
            };

            template<typename PolicyType>
            class poseidon : public poseidon1<PolicyType> {
            public:
                using policy_type = PolicyType;
                using accumulator_tag = accumulators::tag::algebraic_hash<poseidon<policy_type>>;
            };

            template<typename PolicyType>
            class poseidon1_dense : public basic_poseidon1<PolicyType, detail::poseidon1_permutation<PolicyType>,
                                                           detail::poseidon_sponge_padding_mode::pad10> {
            public:
                using policy_type = PolicyType;
                using accumulator_tag = accumulators::tag::algebraic_hash<poseidon1_dense<policy_type>>;
            };

            template<typename PolicyType>
            class poseidon1_padding_free
                : public basic_poseidon1<PolicyType, detail::poseidon1_optimized_permutation<PolicyType>,
                                         detail::poseidon_sponge_padding_mode::padding_free> {
            public:
                using policy_type = PolicyType;
                using accumulator_tag = accumulators::tag::algebraic_hash<poseidon1_padding_free<policy_type>>;
            };

            template<typename PolicyType, detail::poseidon_sponge_padding_mode PaddingMode>
            class basic_poseidon2 {
            public:
                using policy_type = PolicyType;
                using permutation_type = detail::poseidon2_permutation<policy_type>;

                using word_type = typename policy_type::word_type;
                constexpr static const std::size_t word_bits = policy_type::word_bits;

                constexpr static const std::size_t block_words = policy_type::block_words;
                using block_type = typename policy_type::block_type;

                constexpr static const std::size_t digest_bits = policy_type::digest_bits;
                using digest_type = typename policy_type::digest_type;

                struct construction {
                    struct params_type {
                        // This is required by the hash concept.
                    };

                    using type = detail::poseidon_sponge_construction<
                        policy_type, permutation_type, detail::poseidon_sponge_absorb_mode::overwrite, PaddingMode>;
                };

                constexpr static detail::stream_processor_type stream_processor = detail::stream_processor_type::raw;
                using accumulator_tag = accumulators::tag::algebraic_hash<basic_poseidon2<policy_type, PaddingMode>>;
            };

            template<typename PolicyType>
            class poseidon2 : public basic_poseidon2<PolicyType, detail::poseidon_sponge_padding_mode::pad10> {
            public:
                using policy_type = PolicyType;
                using accumulator_tag = accumulators::tag::algebraic_hash<poseidon2<policy_type>>;
            };

            template<typename PolicyType>
            class poseidon2_padding_free
                : public basic_poseidon2<PolicyType, detail::poseidon_sponge_padding_mode::padding_free> {
            public:
                using policy_type = PolicyType;
                using accumulator_tag = accumulators::tag::algebraic_hash<poseidon2_padding_free<policy_type>>;
            };
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif
