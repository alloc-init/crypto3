//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON2_PERMUTATION_HPP
#define CRYPTO3_HASH_POSEIDON2_PERMUTATION_HPP

#include <cstddef>

#include <nil/crypto3/hash/detail/poseidon2/poseidon2_constants.hpp>
#include <nil/crypto3/hash/detail/poseidon2/poseidon2_policy.hpp>
#include <nil/crypto3/hash/detail/poseidon2/poseidon2_round_functions.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                /*!
                 * @brief Standard Poseidon2 permutation.
                 *
                 * The schedule is:
                 *   initial external linear layer,
                 *   RF/2 initial external rounds,
                 *   RP internal rounds,
                 *   RF/2 terminal external rounds.
                 */
                template<poseidon2_policy_type PolicyType>
                class poseidon2_permutation {
                    using policy_type = PolicyType;
                    using constants_type = poseidon2_constants<policy_type>;
                    using round_functions_type = poseidon2_round_functions<policy_type>;

                public:
                    using element_type = typename policy_type::word_type;

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
                        const constants_type &constants = get_constants();

                        round_functions_type::external_linear_layer(state);

                        for (std::size_t i = 0; i < half_full_rounds; ++i) {
                            round_functions_type::initial_external_round(state, constants, i);
                        }

                        for (std::size_t i = 0; i < part_rounds; ++i) {
                            round_functions_type::internal_round(state, constants, i);
                        }

                        for (std::size_t i = 0; i < half_full_rounds; ++i) {
                            round_functions_type::terminal_external_round(state, constants, i);
                        }
                    }

                private:
                    static const constants_type &get_constants() {
                        static const constants_type constants;
                        return constants;
                    }
                };

            }    // namespace detail
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON2_PERMUTATION_HPP
