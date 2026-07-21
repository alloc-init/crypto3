//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON1_ROUND_FUNCTIONS_HPP
#define CRYPTO3_HASH_POSEIDON1_ROUND_FUNCTIONS_HPP

#include <cstddef>

#include <boost/assert.hpp>

#include <nil/crypto3/hash/detail/poseidon1/poseidon1_policy.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                /*!
                 * @brief Shared Poseidon1 round operations.
                 *
                 * The standard and optimized permutations differ in the middle partial-round schedule, but
                 * they use the same state representation, S-box, round-constant addition, and dense MDS
                 * multiplication for full rounds.
                 */
                template<poseidon1_policy_type PolicyType>
                struct poseidon1_round_functions {
                    using policy_type = PolicyType;
                    using element_type = typename policy_type::word_type;
                    using state_type = typename policy_type::state_type;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t half_full_rounds = policy_type::half_full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;
                    constexpr static const std::size_t sbox_power = policy_type::sbox_power;

                    template<typename MatrixType>
                    static state_type matrix_vector_mul(const MatrixType &matrix, const state_type &vector) {
                        state_type result;
                        result.fill(zero());
                        for (std::size_t i = 0; i < state_words; ++i) {
                            for (std::size_t j = 0; j < state_words; ++j) {
                                result[i] += matrix[i][j] * vector[j];
                            }
                        }
                        return result;
                    }

                    template<typename MatrixType>
                    static void product_with_matrix(state_type &state, const MatrixType &matrix) {
                        state = matrix_vector_mul(matrix, state);
                    }

                    template<typename ConstantsType>
                    static void add_round_constants(state_type &state, const ConstantsType &constants,
                                                    std::size_t round_number) {
                        for (std::size_t i = 0; i < state_words; ++i) {
                            state[i] += constants.get_round_constant(round_number, i);
                        }
                    }

                    static void apply_full_sbox(state_type &state) {
                        for (std::size_t i = 0; i < state_words; ++i) {
                            state[i] = state[i].pow(sbox_power);
                        }
                    }

                    static void apply_partial_sbox(state_type &state) {
                        state[0] = state[0].pow(sbox_power);
                    }

                    template<typename ConstantsType>
                    static void full_round(state_type &state, const ConstantsType &constants,
                                           std::size_t round_number) {
                        BOOST_ASSERT_MSG(round_number < half_full_rounds ||
                                             round_number >= half_full_rounds + part_rounds,
                                         "Wrong usage of the full round function of Poseidon1.");

                        add_round_constants(state, constants, round_number);
                        apply_full_sbox(state);
                        constants.product_with_mds_matrix(state);
                    }

                    template<typename ConstantsType>
                    static void partial_round(state_type &state, const ConstantsType &constants,
                                              std::size_t round_number) {
                        BOOST_ASSERT_MSG(round_number >= half_full_rounds &&
                                             round_number < half_full_rounds + part_rounds,
                                         "Wrong usage of the partial round function of Poseidon1.");

                        add_round_constants(state, constants, round_number);
                        apply_partial_sbox(state);
                        constants.product_with_mds_matrix(state);
                    }

                private:
                    static constexpr element_type zero() {
                        return element_type(0u);
                    }
                };

            }    // namespace detail
        }    // namespace hashes
    }        // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON1_ROUND_FUNCTIONS_HPP
