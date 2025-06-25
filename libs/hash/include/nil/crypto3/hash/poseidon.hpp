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
#include <nil/crypto3/hash/detail/poseidon/poseidon_sponge.hpp>
#include <nil/crypto3/hash/detail/poseidon/poseidon_functions.hpp>
#include <nil/crypto3/hash/detail/poseidon/poseidon_permutation.hpp>
#include <nil/crypto3/hash/detail/sponge_construction.hpp>
#include <nil/crypto3/hash/detail/stream_processors/stream_processors_enum.hpp>
namespace nil {
    namespace crypto3 {
        namespace hashes {
#ifdef __ZKLLVM__
            class poseidon {
            public:
                typedef typename algebra::curves::pallas::base_field_type::value_type block_type;

                struct process {
                    block_type operator()(block_type first_input_block, block_type second_input_block) {
                        return __builtin_assigner_poseidon_pallas_base({0, first_input_block, second_input_block})[2];
                    }
                };
            };
#else
            template<typename PolicyType>
            class poseidon {
            public:
                typedef PolicyType policy_type;

                typedef typename policy_type::word_type word_type;
                constexpr static const std::size_t word_bits = policy_type::word_bits;

                constexpr static const std::size_t block_words = policy_type::block_words;
                typedef typename policy_type::block_type block_type;

                constexpr static const std::size_t digest_bits = policy_type::digest_bits;
                using digest_type = typename policy_type::digest_type;

                struct construction {
                    struct params_type {
                        // This is required by 'is_hash' concept.
                    };

                    using type = detail::poseidon_sponge_construction_custom<policy_type>;
                };

                constexpr static detail::stream_processor_type stream_processor = detail::stream_processor_type::raw;
                using accumulator_tag = accumulators::tag::algebraic_hash<poseidon<policy_type>>;
            };

            template<typename PolicyType>
            class original_poseidon {
            public:
                typedef PolicyType policy_type;

                typedef typename policy_type::word_type word_type;
                constexpr static const std::size_t word_bits = policy_type::word_bits;

                constexpr static const std::size_t block_words = policy_type::block_words;
                typedef typename policy_type::block_type block_type;

                constexpr static const std::size_t digest_bits = policy_type::digest_bits;
                using digest_type = typename policy_type::digest_type;

                struct construction {
                    struct params_type {
                        // This is required by 'is_hash' concept.
                    };

                    typedef algebraic_sponge_construction<policy_type,
                        typename policy_type::iv_generator,
                        detail::poseidon_functions<policy_type>,
                        detail::poseidon_functions<policy_type>,
                        detail::poseidon_functions<policy_type>> type;
                };

                constexpr static detail::stream_processor_type stream_processor = detail::stream_processor_type::raw;
                using accumulator_tag = accumulators::tag::algebraic_hash<original_poseidon<PolicyType>>;
            };
#endif
        } // namespace hashes
    } // namespace crypto3
} // namespace nil

#endif
