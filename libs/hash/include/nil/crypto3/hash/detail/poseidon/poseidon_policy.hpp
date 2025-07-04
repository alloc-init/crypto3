//---------------------------------------------------------------------------//
// Copyright (c) 2020 Ilias Khairullin <ilias@nil.foundation>
// Copyright (c) 2020 Mikhail Komarov <nemo@nil.foundation>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON_POLICY_HPP
#define CRYPTO3_HASH_POSEIDON_POLICY_HPP

#include <array>
#include <type_traits>

#include <nil/crypto3/detail/stream_endian.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {
                /*!
                 * @brief Poseidon internal parameters
                 * @tparam FieldType type of field
                 * @tparam Rate Rate of input block for Poseidon permutation in field elements
                 * @tparam Capacity Capacity or inner part of Poseidon permutation in field elements
                 * @tparam Security mode of Poseidon permutation
                 */
                // base_poseidon_policy<FieldType, 128, 4, 1, 5, 8, 60, false> {};
                template<typename FieldType,
                    std::size_t Security, std::size_t Rate, std::size_t Capacity,
                    std::size_t SBoxPower, std::size_t FullRounds, std::size_t PartRounds,
                    bool PastaVersion>
                struct base_poseidon_policy {
                    using field_type = FieldType;

                    constexpr static const std::size_t word_bits = field_type::value_bits;
                    using word_type = typename field_type::value_type;

                    constexpr static const std::size_t block_words = Rate;
                    constexpr static const std::size_t block_bits = block_words * word_bits;

                    using block_type = std::array<word_type, block_words>;
                
                    constexpr static const std::size_t state_words = Rate + Capacity;
                    using state_type = std::array<word_type, state_words>;
                    constexpr static const std::size_t state_bits = state_words * state_bits;

                    constexpr static const std::size_t digest_bits = word_bits;
                    using digest_type = word_type;

                    constexpr static const std::size_t full_rounds = FullRounds;
                    constexpr static const std::size_t half_full_rounds = FullRounds >> 1;
                    constexpr static const std::size_t part_rounds = PartRounds;

                    constexpr static const std::size_t security = Security;
                    constexpr static const std::size_t rate = Rate;
                    constexpr static const std::size_t capacity = Capacity;
                    constexpr static const std::size_t sbox_power = SBoxPower;

                    constexpr static const bool pasta_version = PastaVersion;

                    struct iv_generator {
                        // TODO: maybe it would be done in constexpr way
                        static state_type generate() {
                            static const state_type H0 = []() {
                                state_type H;
                                H.fill(word_type(0u));
                                return H;
                            }();
                            return H0;
                        }
                    };
                };

                /*!
                 * @brief Policy class for the original implementation.
                 * @tparam FieldType Type of the field.
                 * @tparam Security The bit strength of hash, one of [80, 128, 256].
                 * @tparam Rate Rate of input block for Poseidon permutation in field elements. Values of 2 or 4 are used with Merkle Trees.
                 * @tparam Capacity Capacity or inner part of Poseidon permutation in field elements.
                 */
                template<typename FieldType, std::size_t Security, std::size_t Rate, typename Enable = void>
                struct poseidon_policy;

                template<typename FieldType, std::size_t Rate>
                struct poseidon_policy<FieldType, 80, Rate,
                        std::enable_if_t<Rate == 1 || Rate == 2 >> :
                    base_poseidon_policy<FieldType, 80, Rate, 1, 5, 8, 33, false> {};

                template<typename FieldType>
                struct poseidon_policy<FieldType, 80, 4> :
                    base_poseidon_policy<FieldType, 80, 4, 1, 5, 8, 35, false> {};

                template<typename FieldType, std::size_t Rate>
                struct poseidon_policy<FieldType, 128, Rate,
                        std::enable_if_t<Rate == 1 || Rate == 2 >> :
                    base_poseidon_policy<FieldType, 128, Rate, 1, 5, 8, 57, false> {};

                template<typename FieldType>
                struct poseidon_policy<FieldType, 128, 4> :
                    base_poseidon_policy<FieldType, 128, 4, 1, 5, 8, 60, false> {};

                template<typename FieldType>
                struct poseidon_policy<FieldType, 128, 8> :
                    base_poseidon_policy<FieldType, 128, 4, 1, 5, 8, 63, false> {};

                template<typename FieldType, std::size_t Rate>
                struct poseidon_policy<FieldType, 256, Rate,
                        std::enable_if_t<Rate <= 4 >> :
                    base_poseidon_policy<FieldType, 256, Rate, 1, 5, 8, 120, false> {};

                /*!
                 * @brief Policy class for special Poseidon implementation over Pasta curves
                 * This implementation uses X^7 S-boxes,
                 *      changes the order of arc, s-box and mds operations and
                 *      they don't use partial rounds.
                 * Only 1 options is supported, with Rate=2, Capacity=1, 55 full rounds and security of 128 bits.
                 * @tparam FieldType Type of the field.
                 */
                template<typename FieldType>
                struct pasta_poseidon_policy : base_poseidon_policy<FieldType, 128, 2, 1, 7, 55, 0, true> {};

            }    // namespace detail
        }        // namespace hashes
    }            // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON_POLICY_HPP
