//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON1_OPTIMIZED_PERMUTATION_HPP
#define CRYPTO3_HASH_POSEIDON1_OPTIMIZED_PERMUTATION_HPP

#include <cstddef>

#include <nil/crypto3/hash/detail/poseidon1/poseidon1_optimized_constants.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_policy.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_round_functions.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                /*!
                 * @brief Optimized Poseidon1 permutation.
                 *
                 * This permutation is algebraically equivalent to poseidon1_permutation, but evaluates the
                 * partial-round section with optimized constants and sparse matrices.
                 */
                template<poseidon1_policy_type PolicyType>
                class poseidon1_optimized_permutation {
                    using policy_type = PolicyType;
                    using constants_type = poseidon1_optimized_constants<policy_type>;
                    using round_functions_type = poseidon1_round_functions<policy_type>;

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
                        std::size_t round_number = 0;
                        const constants_type &constants = get_constants();

                        // The first RF/2 full rounds are identical to the dense Poseidon1 permutation:
                        // add a full round-constant vector, apply the S-box to every state element, then
                        // multiply by the dense MDS matrix.
                        for (std::size_t i = 0; i < half_full_rounds; ++i) {
                            round_functions_type::full_round(state, constants, round_number++);
                        }

                        // The middle RP partial rounds are the only optimized part. The constants and MDS
                        // products in this section have been rewritten so that each round only needs one
                        // S-box, one scalar round constant, and one sparse matrix product.
                        optimized_partial_rounds(state, constants);
                        round_number += part_rounds;

                        // The final RF/2 full rounds return to the ordinary dense Poseidon1 schedule.
                        for (std::size_t i = half_full_rounds; i < full_rounds; ++i) {
                            round_functions_type::full_round(state, constants, round_number++);
                        }
                    }

                private:
                    static void optimized_partial_rounds(state_type &state, const constants_type &constants) {
                        // All full-state constants that survive the backwards-constant rewrite are added
                        // once, before entering the sparse matrix chain.
                        constants.add_first_partial_round_constants(state);

                        // This dense bridge matrix accounts for the transition from the first rewritten
                        // constant addition into the sequence of sparse partial-round matrices.
                        constants.product_with_pre_sparse_matrix(state);

                        // For the first RP-1 optimized partial rounds, only state[0] is nonlinear. The
                        // rewritten scalar constant is added after that S-box, then the sparse matrix
                        // advances the linear layer.
                        for (std::size_t i = 0; i < part_rounds - 1; ++i) {
                            round_functions_type::apply_partial_sbox(state);
                            state[0] += constants.get_partial_round_constant(i);
                            constants.product_with_sparse_matrix(state, i);
                        }

                        // The last optimized partial round has no scalar constant left after the backwards
                        // rewrite. It is just the final partial S-box followed by the last sparse product.
                        round_functions_type::apply_partial_sbox(state);
                        constants.product_with_sparse_matrix(state, part_rounds - 1);
                    }

                    static const constants_type &get_constants() {
                        static const constants_type constants;
                        return constants;
                    }
                };

            }    // namespace detail
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON1_OPTIMIZED_PERMUTATION_HPP
