//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON_COMMON_POLICY_HPP
#define CRYPTO3_HASH_POSEIDON_COMMON_POLICY_HPP

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                /*!
                 * @brief Shared Poseidon policy surface.
                 *
                 * Poseidon1 and Poseidon2 have different permutation schedules and constant layouts, but
                 * they share the same sponge-facing shape: field words, rate/capacity, state, block, digest,
                 * security level, S-box power, and round counts.
                 */
                template<typename FieldType,
                         std::size_t Security,
                         std::size_t Rate,
                         std::size_t Capacity,
                         std::size_t SBoxPower,
                         std::size_t FullRounds,
                         std::size_t PartRounds,
                         std::size_t DigestBits>
                struct base_standard_poseidon_policy {
                    static_assert(FullRounds % 2 == 0, "Poseidon requires an even number of full rounds.");
                    static_assert(Rate > 0, "Poseidon rate must be positive.");
                    static_assert(Capacity > 0, "Poseidon capacity must be positive.");
                    static_assert(SBoxPower >= 3, "Poseidon S-box power must be at least 3.");
                    static_assert(DigestBits > 0, "Poseidon digest must contain at least one bit.");

                    using field_type = FieldType;

                    constexpr static const std::size_t word_bits = field_type::value_bits;
                    using word_type = typename field_type::value_type;

                    constexpr static const std::size_t block_words = Rate;
                    constexpr static const std::size_t block_bits = block_words * word_bits;
                    using block_type = std::array<word_type, block_words>;

                    constexpr static const std::size_t state_words = Rate + Capacity;
                    constexpr static const std::size_t state_bits = state_words * word_bits;
                    using state_type = std::array<word_type, state_words>;

                    constexpr static const std::size_t digest_bits = DigestBits;
                    constexpr static const std::size_t digest_words = (digest_bits + word_bits - 1) / word_bits;
                    static_assert(digest_words <= Rate,
                                  "Poseidon intentionally supports only OUT <= RATE for now; "
                                  "extend the sponge to squeeze multiple blocks before raising this limit.");
                    using digest_type =
                        std::conditional_t<digest_words == 1, word_type, std::array<word_type, digest_words>>;

                    constexpr static const std::size_t full_rounds = FullRounds;
                    constexpr static const std::size_t half_full_rounds = FullRounds / 2;
                    constexpr static const std::size_t part_rounds = PartRounds;

                    constexpr static const std::size_t security = Security;
                    constexpr static const std::size_t rate = Rate;
                    constexpr static const std::size_t capacity = Capacity;
                    constexpr static const std::size_t sbox_power = SBoxPower;

                    struct iv_generator {
                        static state_type generate() {
                            state_type h;
                            h.fill(word_type(0u));
                            return h;
                        }
                    };
                };

                template<typename PolicyType>
                concept poseidon_common_policy_type =
                    requires(typename PolicyType::word_type word) {
                        typename PolicyType::field_type;
                        typename PolicyType::word_type;
                        typename PolicyType::block_type;
                        typename PolicyType::state_type;
                        typename PolicyType::digest_type;
                        typename PolicyType::iv_generator;

                        { PolicyType::word_bits } -> std::convertible_to<std::size_t>;
                        { PolicyType::block_words } -> std::convertible_to<std::size_t>;
                        { PolicyType::block_bits } -> std::convertible_to<std::size_t>;
                        { PolicyType::state_words } -> std::convertible_to<std::size_t>;
                        { PolicyType::state_bits } -> std::convertible_to<std::size_t>;
                        { PolicyType::digest_words } -> std::convertible_to<std::size_t>;
                        { PolicyType::digest_bits } -> std::convertible_to<std::size_t>;
                        { PolicyType::full_rounds } -> std::convertible_to<std::size_t>;
                        { PolicyType::half_full_rounds } -> std::convertible_to<std::size_t>;
                        { PolicyType::part_rounds } -> std::convertible_to<std::size_t>;
                        { PolicyType::security } -> std::convertible_to<std::size_t>;
                        { PolicyType::rate } -> std::convertible_to<std::size_t>;
                        { PolicyType::capacity } -> std::convertible_to<std::size_t>;
                        { PolicyType::sbox_power } -> std::convertible_to<std::size_t>;

                        { word.pow(PolicyType::sbox_power) } -> std::same_as<typename PolicyType::word_type>;
                        { PolicyType::iv_generator::generate() } -> std::same_as<typename PolicyType::state_type>;
                    } && (PolicyType::rate > 0) && (PolicyType::capacity > 0) &&
                    (PolicyType::state_words == PolicyType::rate + PolicyType::capacity) &&
                    (PolicyType::block_words == PolicyType::rate) && (PolicyType::digest_bits > 0) &&
                    (PolicyType::digest_words > 0) && (PolicyType::digest_words <= PolicyType::block_words) &&
                    (PolicyType::digest_words ==
                     (PolicyType::digest_bits + PolicyType::word_bits - 1) / PolicyType::word_bits) &&
                    (PolicyType::full_rounds >= 6) && (PolicyType::full_rounds % 2 == 0) &&
                    (PolicyType::half_full_rounds == PolicyType::full_rounds / 2) && (PolicyType::sbox_power >= 3);

                struct poseidon1_policy_tag { };
                struct poseidon2_policy_tag { };

                template<typename PolicyType, typename AlgorithmTag>
                concept poseidon_algorithm_policy_type = poseidon_common_policy_type<PolicyType> && requires {
                    typename PolicyType::poseidon_policy_tag;
                } && std::same_as<typename PolicyType::poseidon_policy_tag, AlgorithmTag>;

            }    // namespace detail
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON_COMMON_POLICY_HPP
