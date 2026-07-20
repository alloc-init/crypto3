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
#include <nil/crypto3/algebra/matrix/math.hpp>
#include <nil/crypto3/algebra/vector/vector.hpp>
#include <nil/crypto3/algebra/vector/math.hpp>

#include <nil/crypto3/hash/detail/poseidon/original_constants.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_policy.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                template<poseidon1_policy_type PolicyType>
                class poseidon1_constants {
                public:
                    using policy_type = PolicyType;
                    using element_type = typename policy_type::word_type;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t full_rounds = policy_type::full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;

                    using state_vector_type = algebra::vector<element_type, state_words>;
                    using mds_matrix_type = algebra::matrix<element_type, state_words, state_words>;
                    using round_constants_type = algebra::matrix<element_type, full_rounds + part_rounds, state_words>;

                    using constants_data_type = poseidon_original_constants_data<policy_type>;

                    poseidon1_constants() {
                        // The legacy constants table is stored in the orientation expected by the old
                        // row-vector vectmatmul path after this transpose. Keep that convention so the
                        // new standard Poseidon1 permutation can be checked against the existing vectors.
                        for (std::size_t i = 0; i < state_words; ++i) {
                            for (std::size_t j = 0; j < state_words; ++j) {
                                mds_matrix[i][j] = constants_data_type::mds_matrix[j][i];
                            }
                        }
                    }

                    inline const element_type &get_round_constant(std::size_t round, std::size_t index) const {
                        return constants_data_type::round_constants[round][index];
                    }

                    inline void product_with_mds_matrix(state_vector_type &state) const {
                        state = algebra::vectmatmul(state, mds_matrix);
                    }

                private:
                    mds_matrix_type mds_matrix;
                };

            }    // namespace detail
        }    // namespace hashes
    }        // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON1_CONSTANTS_HPP
