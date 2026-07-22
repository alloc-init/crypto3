//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON1_CONSTANTS_HPP
#define CRYPTO3_HASH_POSEIDON1_CONSTANTS_HPP

#include <cstddef>

#include <nil/crypto3/algebra/matrix/matrix.hpp>

#include <nil/crypto3/hash/detail/poseidon/original_constants.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_policy.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_round_functions.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                template<poseidon1_policy_type PolicyType>
                class poseidon1_constants {
                public:
                    using policy_type = PolicyType;
                    using element_type = typename policy_type::word_type;
                    using state_type = typename policy_type::state_type;
                    using round_functions_type = poseidon1_round_functions<policy_type>;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t full_rounds = policy_type::full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;

                    using mds_matrix_type = algebra::matrix<element_type, state_words, state_words>;
                    using round_constants_type = algebra::matrix<element_type, full_rounds + part_rounds, state_words>;

                    using constants_data_type = poseidon_original_constants_data<policy_type>;

                    poseidon1_constants() {
                        // Store the raw MDS table and use column-vector multiplication. This is the same
                        // linear map as the old dense row-vector path with a transposed MDS, but it matches
                        // the convention used by the optimized sparse derivation.
                        for (std::size_t i = 0; i < state_words; ++i) {
                            for (std::size_t j = 0; j < state_words; ++j) {
                                mds_matrix[i][j] = constants_data_type::mds_matrix[i][j];
                            }
                        }
                    }

                    inline const element_type &get_round_constant(std::size_t round, std::size_t index) const {
                        return constants_data_type::round_constants[round][index];
                    }

                    inline void product_with_mds_matrix(state_type &state) const {
                        round_functions_type::product_with_matrix(state, mds_matrix);
                    }

                    inline const mds_matrix_type &get_mds_matrix() const {
                        return mds_matrix;
                    }

                private:
                    mds_matrix_type mds_matrix;
                };

            }    // namespace detail
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON1_CONSTANTS_HPP
