//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON1_POLICY_HPP
#define CRYPTO3_HASH_POSEIDON1_POLICY_HPP

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                /*!
                 * @brief Standard Poseidon1 permutation parameters.
                 *
                 * This policy family is separate from the legacy crypto3 Poseidon policies. It models the
                 * original Poseidon1 round structure: RF/2 full rounds, RP partial rounds, RF/2 full rounds,
                 * with AddRoundConstants -> S-box -> MDS inside each round.
                 */
                template<typename FieldType, std::size_t Security, std::size_t Rate, std::size_t Capacity,
                         std::size_t SBoxPower, std::size_t FullRounds, std::size_t PartRounds>
                struct base_poseidon1_policy {
                    static_assert(FullRounds % 2 == 0, "Poseidon1 requires an even number of full rounds.");
                    static_assert(Rate > 0, "Poseidon1 rate must be positive.");
                    static_assert(Capacity > 0, "Poseidon1 capacity must be positive.");
                    static_assert(SBoxPower >= 3, "Poseidon1 S-box power must be at least 3.");

                    using field_type = FieldType;

                    constexpr static const std::size_t word_bits = field_type::value_bits;
                    using word_type = typename field_type::value_type;

                    constexpr static const std::size_t block_words = Rate;
                    constexpr static const std::size_t block_bits = block_words * word_bits;
                    using block_type = std::array<word_type, block_words>;

                    constexpr static const std::size_t state_words = Rate + Capacity;
                    constexpr static const std::size_t state_bits = state_words * word_bits;
                    using state_type = std::array<word_type, state_words>;

                    constexpr static const std::size_t digest_words = 1;
                    constexpr static const std::size_t digest_bits = word_bits;
                    using digest_type = word_type;

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
                concept poseidon1_policy_type = requires(typename PolicyType::word_type word) {
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
                    (PolicyType::block_words == PolicyType::rate) &&
                    (PolicyType::full_rounds >= 6) &&
                    (PolicyType::full_rounds % 2 == 0) &&
                    (PolicyType::half_full_rounds == PolicyType::full_rounds / 2) &&
                    (PolicyType::sbox_power >= 3);

                /*!
                 * @brief Standard Poseidon1 policy with the round counts used by the original Poseidon paper
                 * parameter sets already present in crypto3.
                 *
                 * Only capacity=1 is modeled here because the existing checked-in constants are capacity-1
                 * parameter sets. Add a new specialization when adding constants for another capacity.
                 */
                template<typename FieldType, std::size_t Security, std::size_t Rate, std::size_t Capacity = 1,
                         typename Enable = void>
                struct poseidon1_policy;

                template<typename FieldType, std::size_t Rate>
                struct poseidon1_policy<FieldType, 80, Rate, 1, std::enable_if_t<Rate == 1 || Rate == 2>>
                    : base_poseidon1_policy<FieldType, 80, Rate, 1, 5, 8, 33> { };

                template<typename FieldType>
                struct poseidon1_policy<FieldType, 80, 4, 1, void>
                    : base_poseidon1_policy<FieldType, 80, 4, 1, 5, 8, 35> { };

                template<typename FieldType, std::size_t Rate>
                struct poseidon1_policy<FieldType, 128, Rate, 1, std::enable_if_t<Rate == 1 || Rate == 2>>
                    : base_poseidon1_policy<FieldType, 128, Rate, 1, 5, 8, 57> { };

                template<typename FieldType>
                struct poseidon1_policy<FieldType, 128, 4, 1, void>
                    : base_poseidon1_policy<FieldType, 128, 4, 1, 5, 8, 60> { };

                template<typename FieldType, std::size_t Rate>
                struct poseidon1_policy<FieldType, 256, Rate, 1, std::enable_if_t<(Rate <= 4)>>
                    : base_poseidon1_policy<FieldType, 256, Rate, 1, 5, 8, 120> { };

            }    // namespace detail
        }    // namespace hashes
    }        // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON1_POLICY_HPP
