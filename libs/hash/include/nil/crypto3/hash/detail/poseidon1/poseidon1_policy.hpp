//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON1_POLICY_HPP
#define CRYPTO3_HASH_POSEIDON1_POLICY_HPP

#include <cstddef>
#include <type_traits>

#include <nil/crypto3/hash/detail/poseidon_common/poseidon_policy.hpp>

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
                         std::size_t SBoxPower, std::size_t FullRounds, std::size_t PartRounds, std::size_t DigestBits>
                struct base_poseidon1_policy
                    : base_standard_poseidon_policy<FieldType, Security, Rate, Capacity, SBoxPower, FullRounds,
                                                    PartRounds, DigestBits> {
                    using poseidon_policy_tag = poseidon1_policy_tag;
                };

                template<typename PolicyType>
                concept poseidon1_policy_type = poseidon_algorithm_policy_type<PolicyType, poseidon1_policy_tag>;

                /*!
                 * @brief Standard Poseidon1 policy with the round counts used by the original Poseidon paper
                 * parameter sets already present in crypto3.
                 *
                 * Only capacity=1 is modeled here because the existing checked-in constants are capacity-1
                 * parameter sets. Add a new specialization when adding constants for another capacity.
                 */
                template<typename FieldType, std::size_t Security, std::size_t Rate, std::size_t Capacity = 1,
                         std::size_t DigestBits = FieldType::value_bits, typename Enable = void>
                struct poseidon1_policy;

                template<typename FieldType, std::size_t Rate, std::size_t DigestBits>
                struct poseidon1_policy<FieldType, 80, Rate, 1, DigestBits, std::enable_if_t<Rate == 1 || Rate == 2>>
                    : base_poseidon1_policy<FieldType, 80, Rate, 1, 5, 8, 33, DigestBits> { };

                template<typename FieldType, std::size_t DigestBits>
                struct poseidon1_policy<FieldType, 80, 4, 1, DigestBits, void>
                    : base_poseidon1_policy<FieldType, 80, 4, 1, 5, 8, 35, DigestBits> { };

                template<typename FieldType, std::size_t Rate, std::size_t DigestBits>
                struct poseidon1_policy<FieldType, 128, Rate, 1, DigestBits, std::enable_if_t<Rate == 1 || Rate == 2>>
                    : base_poseidon1_policy<FieldType, 128, Rate, 1, 5, 8, 57, DigestBits> { };

                template<typename FieldType, std::size_t DigestBits>
                struct poseidon1_policy<FieldType, 128, 4, 1, DigestBits, void>
                    : base_poseidon1_policy<FieldType, 128, 4, 1, 5, 8, 60, DigestBits> { };

                template<typename FieldType, std::size_t Rate, std::size_t DigestBits>
                struct poseidon1_policy<FieldType, 256, Rate, 1, DigestBits, std::enable_if_t<(Rate <= 4)>>
                    : base_poseidon1_policy<FieldType, 256, Rate, 1, 5, 8, 120, DigestBits> { };

            }    // namespace detail
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON1_POLICY_HPP
