//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON2_POLICY_HPP
#define CRYPTO3_HASH_POSEIDON2_POLICY_HPP

#include <cstddef>
#include <type_traits>

#include <nil/crypto3/algebra/fields/alt_bn128/scalar_field.hpp>
#include <nil/crypto3/hash/detail/poseidon_common/poseidon_policy.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                /*!
                 * @brief Standard Poseidon2 permutation parameters.
                 *
                 * This policy family is separate from the legacy crypto3 Poseidon policies and from
                 * the Poseidon1 policies because Poseidon2 has a different linear-layer structure:
                 * external full rounds around internal partial rounds.
                 */
                template<typename FieldType, std::size_t Security, std::size_t Rate, std::size_t Capacity,
                         std::size_t SBoxPower, std::size_t FullRounds, std::size_t PartRounds, std::size_t DigestBits>
                struct base_poseidon2_policy
                    : base_standard_poseidon_policy<FieldType, Security, Rate, Capacity, SBoxPower, FullRounds,
                                                    PartRounds, DigestBits> {
                    using poseidon_policy_tag = poseidon2_policy_tag;
                };

                template<typename PolicyType>
                concept poseidon2_policy_type = poseidon_algorithm_policy_type<PolicyType, poseidon2_policy_tag>;

                /*!
                 * @brief Standard Poseidon2 policy.
                 *
                 * The first supported instance is BN254 width 3, matching the Poseidon2 paper's
                 * `(n, t, d) = (256, 3, 5)` parameter set.
                 */
                template<typename FieldType, std::size_t Security, std::size_t Rate, std::size_t Capacity = 1,
                         std::size_t DigestBits = FieldType::value_bits, typename Enable = void>
                struct poseidon2_policy;

                template<std::size_t DigestBits>
                struct poseidon2_policy<algebra::fields::alt_bn128_scalar_field<254>, 128, 2, 1, DigestBits, void>
                    : base_poseidon2_policy<algebra::fields::alt_bn128_scalar_field<254>, 128, 2, 1, 5, 8, 56,
                                            DigestBits> { };

            }    // namespace detail
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON2_POLICY_HPP
