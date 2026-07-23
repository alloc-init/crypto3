//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON2_ROUND_FUNCTIONS_HPP
#define CRYPTO3_HASH_POSEIDON2_ROUND_FUNCTIONS_HPP

#include <array>
#include <cstddef>

#include <boost/assert.hpp>

#include <nil/crypto3/hash/detail/poseidon2/poseidon2_policy.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                /*!
                 * @brief Shared Poseidon2 round operations.
                 *
                 * Poseidon2 uses external full rounds and internal partial rounds. The external
                 * linear layer is the light MDS construction for widths 2, 3, and multiples of 4.
                 * The internal linear layer is multiplication by `1 + Diag(v)`, where `v` is
                 * supplied by the concrete constants.
                 */
                template<poseidon2_policy_type PolicyType>
                struct poseidon2_round_functions {
                    using policy_type = PolicyType;
                    using element_type = typename policy_type::word_type;
                    using state_type = typename policy_type::state_type;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t half_full_rounds = policy_type::half_full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;
                    constexpr static const std::size_t sbox_power = policy_type::sbox_power;

                    static_assert(state_words == 2 || state_words == 3 || state_words % 4 == 0,
                                  "Poseidon2 external linear layer supports width 2, width 3, "
                                  "and multiples of 4.");

                    static void external_linear_layer(state_type &state) {
                        if constexpr (state_words == 2 || state_words == 3) {
                            element_type sum = zero();
                            for (std::size_t i = 0; i < state_words; ++i) {
                                sum += state[i];
                            }
                            for (std::size_t i = 0; i < state_words; ++i) {
                                state[i] += sum;
                            }
                        } else {
                            for (std::size_t i = 0; i < state_words; i += 4) {
                                apply_external_mds_4(state, i);
                            }

                            std::array<element_type, 4> sums = {zero(), zero(), zero(), zero()};
                            for (std::size_t i = 0; i < state_words; i += 4) {
                                for (std::size_t j = 0; j < 4; ++j) {
                                    sums[j] += state[i + j];
                                }
                            }

                            for (std::size_t i = 0; i < state_words; ++i) {
                                state[i] += sums[i % 4];
                            }
                        }
                    }

                    template<typename ConstantsType>
                    static void internal_linear_layer(state_type &state, const ConstantsType &constants) {
                        element_type sum = zero();
                        for (std::size_t i = 0; i < state_words; ++i) {
                            sum += state[i];
                        }

                        for (std::size_t i = 0; i < state_words; ++i) {
                            state[i] *= constants.get_internal_diagonal_minus_one(i);
                            state[i] += sum;
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
                    static void add_initial_external_round_constants(state_type &state, const ConstantsType &constants,
                                                                     std::size_t round_number) {
                        BOOST_ASSERT_MSG(round_number < half_full_rounds,
                                         "Wrong usage of Poseidon2 initial external constants.");
                        for (std::size_t i = 0; i < state_words; ++i) {
                            state[i] += constants.get_initial_external_round_constant(round_number, i);
                        }
                    }

                    template<typename ConstantsType>
                    static void add_terminal_external_round_constants(state_type &state,
                                                                      const ConstantsType &constants,
                                                                      std::size_t round_number) {
                        BOOST_ASSERT_MSG(round_number < half_full_rounds,
                                         "Wrong usage of Poseidon2 terminal external constants.");
                        for (std::size_t i = 0; i < state_words; ++i) {
                            state[i] += constants.get_terminal_external_round_constant(round_number, i);
                        }
                    }

                    template<typename ConstantsType>
                    static void initial_external_round(state_type &state, const ConstantsType &constants,
                                                       std::size_t round_number) {
                        add_initial_external_round_constants(state, constants, round_number);
                        apply_full_sbox(state);
                        external_linear_layer(state);
                    }

                    template<typename ConstantsType>
                    static void internal_round(state_type &state, const ConstantsType &constants,
                                               std::size_t round_number) {
                        BOOST_ASSERT_MSG(round_number < part_rounds,
                                         "Wrong usage of Poseidon2 internal round constants.");
                        state[0] += constants.get_internal_round_constant(round_number);
                        apply_partial_sbox(state);
                        internal_linear_layer(state, constants);
                    }

                    template<typename ConstantsType>
                    static void terminal_external_round(state_type &state, const ConstantsType &constants,
                                                        std::size_t round_number) {
                        add_terminal_external_round_constants(state, constants, round_number);
                        apply_full_sbox(state);
                        external_linear_layer(state);
                    }

                private:
                    static void apply_external_mds_4(state_type &state, std::size_t offset) {
                        const element_type t0 = state[offset] + state[offset + 1];
                        const element_type t1 = state[offset + 2] + state[offset + 3];
                        const element_type t2 = mul_by_2(state[offset + 1]) + t1;
                        const element_type t3 = mul_by_2(state[offset + 3]) + t0;
                        const element_type t4 = mul_by_2(mul_by_2(t1)) + t3;
                        const element_type t5 = mul_by_2(mul_by_2(t0)) + t2;
                        const element_type t6 = t3 + t5;
                        const element_type t7 = t2 + t4;

                        state[offset] = t6;
                        state[offset + 1] = t5;
                        state[offset + 2] = t7;
                        state[offset + 3] = t4;
                    }

                    static element_type mul_by_2(const element_type &value) {
                        return value + value;
                    }

                    static constexpr element_type zero() {
                        return element_type(0u);
                    }
                };

            }    // namespace detail
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON2_ROUND_FUNCTIONS_HPP
