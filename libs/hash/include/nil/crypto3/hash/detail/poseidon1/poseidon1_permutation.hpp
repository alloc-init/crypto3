//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON1_PERMUTATION_HPP
#define CRYPTO3_HASH_POSEIDON1_PERMUTATION_HPP

#include <cstddef>

#include <boost/assert.hpp>

#include <nil/crypto3/hash/detail/poseidon1/poseidon1_constants.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_policy.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                /*!
                 * @brief Textbook Poseidon1 permutation.
                 *
                 * This is the unoptimized standard Poseidon1 round structure:
                 *   RF/2 full rounds -> RP partial rounds -> RF/2 full rounds.
                 * Every round uses AddRoundConstants -> S-box -> MDS. The optimized sparse partial-round
                 * backend should be added as a separate permutation and proven equivalent to this one.
                 */
                template<poseidon1_policy_type PolicyType>
                class poseidon1_permutation {
                    using policy_type = PolicyType;
                    using constants_type = poseidon1_constants<policy_type>;

                public:
                    using element_type = typename policy_type::word_type;
                    using state_vector_type = typename constants_type::state_vector_type;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    using state_type = typename policy_type::state_type;

                    constexpr static const std::size_t block_words = policy_type::block_words;
                    using block_type = typename policy_type::block_type;

                    constexpr static const std::size_t full_rounds = policy_type::full_rounds;
                    constexpr static const std::size_t half_full_rounds = policy_type::half_full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;
                    constexpr static const std::size_t sbox_power = policy_type::sbox_power;

                    constexpr static const std::size_t word_bits = policy_type::word_bits;
                    using word_type = typename policy_type::word_type;

                    static void permute(state_type &state) {
                        std::size_t round_number = 0;

                        state_vector_type state_vector;
                        for (std::size_t i = 0; i < state_words; ++i) {
                            state_vector[i] = state[i];
                        }

                        for (std::size_t i = 0; i < half_full_rounds; ++i) {
                            full_round(state_vector, round_number++);
                        }

                        for (std::size_t i = 0; i < part_rounds; ++i) {
                            partial_round(state_vector, round_number++);
                        }

                        for (std::size_t i = half_full_rounds; i < full_rounds; ++i) {
                            full_round(state_vector, round_number++);
                        }

                        for (std::size_t i = 0; i < state_words; ++i) {
                            state[i] = state_vector[i];
                        }
                    }

                private:
                    static void full_round(state_vector_type &state, std::size_t round_number) {
                        BOOST_ASSERT_MSG(round_number < half_full_rounds ||
                                             round_number >= half_full_rounds + part_rounds,
                                         "Wrong usage of the full round function of Poseidon1.");

                        const constants_type &constants = get_constants();
                        for (std::size_t i = 0; i < state_words; ++i) {
                            state[i] += constants.get_round_constant(round_number, i);
                            state[i] = state[i].pow(sbox_power);
                        }
                        constants.product_with_mds_matrix(state);
                    }

                    static void partial_round(state_vector_type &state, std::size_t round_number) {
                        BOOST_ASSERT_MSG(round_number >= half_full_rounds &&
                                             round_number < half_full_rounds + part_rounds,
                                         "Wrong usage of the partial round function of Poseidon1.");

                        const constants_type &constants = get_constants();
                        for (std::size_t i = 0; i < state_words; ++i) {
                            state[i] += constants.get_round_constant(round_number, i);
                        }
                        state[0] = state[0].pow(sbox_power);
                        constants.product_with_mds_matrix(state);
                    }

                    static const constants_type &get_constants() {
                        static const constants_type constants;
                        return constants;
                    }
                };

            }    // namespace detail
        }    // namespace hashes
    }        // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON1_PERMUTATION_HPP
