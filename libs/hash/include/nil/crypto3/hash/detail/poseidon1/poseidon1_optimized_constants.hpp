//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON1_OPTIMIZED_CONSTANTS_HPP
#define CRYPTO3_HASH_POSEIDON1_OPTIMIZED_CONSTANTS_HPP

#include <array>
#include <cstddef>

#include <boost/assert.hpp>

#include <nil/crypto3/algebra/matrix/math.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_constants.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_policy.hpp>
#include <nil/crypto3/hash/detail/poseidon1/poseidon1_round_functions.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                /*!
                 * @brief Optimized constants for the Poseidon1 sparse partial-round schedule.
                 *
                 * These constants are derived from the same dense MDS and round constants used by
                 * poseidon1_permutation. The derivation follows the standard Poseidon1 optimization:
                 * move most partial-round constants backwards through MDS^{-1}, then replace the dense
                 * MDS applications in the partial-round section with one pre-sparse matrix and sparse
                 * matrices.
                 */
                template<poseidon1_policy_type PolicyType>
                class poseidon1_optimized_constants : public poseidon1_constants<PolicyType> {
                public:
                    using policy_type = PolicyType;
                    using base_type = poseidon1_constants<policy_type>;
                    using element_type = typename policy_type::word_type;
                    using state_type = typename policy_type::state_type;
                    using round_functions_type = poseidon1_round_functions<policy_type>;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t full_rounds = policy_type::full_rounds;
                    constexpr static const std::size_t half_full_rounds = policy_type::half_full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;

                    static_assert(state_words > 1, "Optimized Poseidon1 requires at least two state elements.");
                    static_assert(part_rounds > 0, "Optimized Poseidon1 requires at least one partial round.");

                    using matrix_type = typename base_type::mds_matrix_type;
                    using partial_matrix_type = algebra::matrix<element_type, state_words - 1, state_words - 1>;
                    using partial_vector_type = std::array<element_type, state_words - 1>;
                    using partial_round_constants_type = std::array<element_type, part_rounds - 1>;
                    using sparse_vectors_type = std::array<partial_vector_type, part_rounds>;
                    using sparse_first_rows_type = std::array<state_type, part_rounds>;

                    using constants_data_type = typename base_type::constants_data_type;

                    poseidon1_optimized_constants() {
                        // The optimized round constants are obtained by repeatedly moving constants
                        // backwards through an MDS layer. That requires MDS^{-1}; the sparse matrices are
                        // derived independently from the MDS itself. Use crypto3's matrix inverse here
                        // rather than carrying another Gauss-Jordan implementation in the hash code.
                        const matrix_type &mds_matrix = this->get_mds_matrix();
                        const matrix_type mds_inverse = algebra::inverse(mds_matrix);
                        const auto partial_round_constants = load_partial_round_constants();

                        compute_equivalent_round_constants(partial_round_constants, mds_inverse);
                        compute_equivalent_matrices(mds_matrix);
                    }

                    inline const element_type &get_partial_round_constant(std::size_t round) const {
                        BOOST_ASSERT(round < part_rounds - 1);
                        // compute_equivalent_round_constants() rewrites the partial-round constants by
                        // pushing later full-state constants backward through MDS^{-1}. After that rewrite,
                        // all partial rounds except the last carry one scalar constant on state[0];
                        // constants for the other state elements were folded into earlier rounds.
                        return optimized_partial_round_constants[round];
                    }

                    inline void add_first_partial_round_constants(state_type &state) const {
                        // The first optimized partial step is still a full-state constant addition. It is
                        // the accumulation point for constants pushed back from later partial rounds.
                        for (std::size_t i = 0; i < state_words; ++i) {
                            state[i] += first_partial_round_constants[i];
                        }
                    }

                    inline void product_with_pre_sparse_matrix(state_type &state) const {
                        // This one dense product bridges the full-state constant addition at the start of
                        // the partial-round section into the sparse schedule used for the remaining
                        // partial rounds.
                        round_functions_type::product_with_matrix(state, pre_sparse_matrix);
                    }

                    inline void product_with_sparse_matrix(state_type &state, std::size_t round) const {
                        BOOST_ASSERT(round < part_rounds);

                        // Each optimized partial-round matrix has the form:
                        //
                        //   [ first_row[0]  first_row[1]  first_row[2]  ... ]
                        //   [ v[0]          1             0             ... ]
                        //   [ v[1]          0             1             ... ]
                        //   [ ...           ...           ...           I   ]
                        //
                        // Multiplying by it only needs one dense dot product for state[0]. Every other
                        // coordinate keeps its old value and receives old_state[0] * v[i - 1].
                        const element_type old_first_element = state[0];
                        element_type new_first_element = zero();
                        for (std::size_t i = 0; i < state_words; ++i) {
                            new_first_element += sparse_first_rows[round][i] * state[i];
                        }

                        state[0] = new_first_element;
                        for (std::size_t i = 1; i < state_words; ++i) {
                            state[i] += old_first_element * sparse_columns[round][i - 1];
                        }
                    }

                private:
                    using partial_round_constants_matrix_type = std::array<state_type, part_rounds>;

                    static constexpr element_type zero() {
                        return element_type(0u);
                    }

                    static constexpr element_type one() {
                        return element_type(1u);
                    }

                    static matrix_type zero_matrix() {
                        matrix_type result;
                        for (std::size_t i = 0; i < state_words; ++i) {
                            for (std::size_t j = 0; j < state_words; ++j) {
                                result[i][j] = zero();
                            }
                        }
                        return result;
                    }

                    static partial_round_constants_matrix_type load_partial_round_constants() {
                        // The checked-in round constants are laid out as full rounds, then partial rounds,
                        // then final full rounds. Only the middle slice is rewritten by the optimization.
                        partial_round_constants_matrix_type result;
                        for (std::size_t round = 0; round < part_rounds; ++round) {
                            for (std::size_t i = 0; i < state_words; ++i) {
                                result[round][i] = constants_data_type::round_constants[half_full_rounds + round][i];
                            }
                        }
                        return result;
                    }

                    static partial_matrix_type submatrix_inverse(const matrix_type &matrix) {
                        // The sparse factorization needs the inverse of the lower-right (t-1) x (t-1)
                        // block of the current dense matrix. That block describes how non-S-box state
                        // elements mix with each other.
                        partial_matrix_type submatrix;

                        for (std::size_t i = 0; i < state_words - 1; ++i) {
                            for (std::size_t j = 0; j < state_words - 1; ++j) {
                                submatrix[i][j] = matrix[i + 1][j + 1];
                            }
                        }

                        // Delegate the inversion itself to crypto3's algebra matrix helper; this keeps the
                        // optimized Poseidon code focused on the Poseidon-specific sparse factorization.
                        return algebra::inverse(submatrix);
                    }

                    static partial_vector_type partial_matrix_vector_mul(const partial_matrix_type &matrix,
                                                                         const partial_vector_type &vector) {
                        // Product on the non-first coordinates. It is used to solve for w_hat, the values
                        // that go into the sparse first row.
                        partial_vector_type result;
                        result.fill(zero());
                        for (std::size_t i = 0; i < state_words - 1; ++i) {
                            for (std::size_t j = 0; j < state_words - 1; ++j) {
                                result[i] += matrix[i][j] * vector[j];
                            }
                        }
                        return result;
                    }

                    void compute_equivalent_round_constants(const partial_round_constants_matrix_type &round_constants,
                                                            const matrix_type &mds_inverse) {
                        // Dense partial rounds add a full vector of constants before applying the S-box to
                        // state[0] and multiplying by MDS. Since coordinates 1..t-1 are linear in partial
                        // rounds, their constants can be moved backwards through the previous MDS layer.
                        //
                        // Working from the last partial round to the first:
                        //   1. multiply the constants to be moved by MDS^{-1};
                        //   2. keep coordinate 0 as the scalar constant for the optimized round;
                        //   3. add coordinates 1..t-1 into the preceding round's full-state constants.
                        //
                        // At the end, first_partial_round_constants is the only full-state constant vector
                        // left in the optimized partial-round section.
                        state_type current_round_constants = round_constants[part_rounds - 1];
                        std::array<element_type, part_rounds> optimized_constants;
                        optimized_constants.fill(zero());

                        for (std::size_t round = part_rounds - 1; round > 0; --round) {
                            const std::size_t source_round = round - 1;
                            const state_type previous_round_constants =
                                round_functions_type::matrix_vector_mul(mds_inverse, current_round_constants);

                            optimized_constants[source_round + 1] = previous_round_constants[0];
                            current_round_constants = round_constants[source_round];
                            for (std::size_t i = 1; i < state_words; ++i) {
                                current_round_constants[i] += previous_round_constants[i];
                            }
                        }

                        first_partial_round_constants = current_round_constants;
                        for (std::size_t i = 0; i < part_rounds - 1; ++i) {
                            optimized_partial_round_constants[i] = optimized_constants[i + 1];
                        }
                    }

                    void compute_equivalent_matrices(const matrix_type &mds) {
                        // The dense partial section would apply MDS after every partial S-box. The
                        // optimization replaces those dense products with:
                        //
                        //   pre_sparse_matrix, then part_rounds sparse matrices.
                        //
                        // The derivation follows the standard Poseidon1 sparse-MDS decomposition. At each
                        // step, multiplied_matrix is split so that the next sparse matrix has an arbitrary
                        // first row, an identity lower-right block, and only one non-zero column below the
                        // first element.
                        const matrix_type mds_transposed = algebra::transpose(mds);
                        matrix_type multiplied_matrix = mds_transposed;
                        matrix_type current_matrix = zero_matrix();

                        sparse_vectors_type equivalent_sparse_columns;
                        sparse_vectors_type equivalent_sparse_rows;

                        for (std::size_t round = 0; round < part_rounds; ++round) {
                            // Values from the first row become the below-diagonal sparse column after the
                            // final reversal. They are used as v in:
                            //   state[i] += old_state[0] * v[i - 1].
                            for (std::size_t i = 0; i < state_words - 1; ++i) {
                                equivalent_sparse_columns[round][i] = multiplied_matrix[0][i + 1];
                            }

                            // Values from the first column describe how state[0] flows into the other
                            // coordinates. Solving against the lower-right block gives w_hat, which is
                            // embedded into the sparse first row.
                            partial_vector_type column;
                            for (std::size_t i = 0; i < state_words - 1; ++i) {
                                column[i] = multiplied_matrix[i + 1][0];
                            }

                            const partial_matrix_type submatrix_inverse_value = submatrix_inverse(multiplied_matrix);
                            equivalent_sparse_rows[round] = partial_matrix_vector_mul(submatrix_inverse_value, column);

                            // current_matrix is the dense transition left after extracting one sparse
                            // matrix. Its first row and first column are cleared except for [0][0], so
                            // the non-first coordinates are carried by the identity part of the sparse
                            // matrix.
                            current_matrix = multiplied_matrix;
                            current_matrix[0][0] = one();
                            for (std::size_t i = 1; i < state_words; ++i) {
                                current_matrix[i][0] = zero();
                            }
                            for (std::size_t j = 1; j < state_words; ++j) {
                                current_matrix[0][j] = zero();
                            }

                            multiplied_matrix = algebra::matmul(mds_transposed, current_matrix);
                        }

                        // The factorization is produced in reverse order. Store the bridge matrix and then
                        // reverse the sparse vectors so the permutation can apply them in forward order.
                        pre_sparse_matrix = algebra::transpose(current_matrix);
                        for (std::size_t round = 0; round < part_rounds; ++round) {
                            const std::size_t reversed_round = part_rounds - 1 - round;
                            sparse_columns[round] = equivalent_sparse_columns[reversed_round];

                            // The sparse first row is [MDS[0][0], w_hat...]. product_with_sparse_matrix()
                            // uses it for the one dense dot product that updates state[0].
                            sparse_first_rows[round][0] = mds[0][0];
                            for (std::size_t i = 1; i < state_words; ++i) {
                                sparse_first_rows[round][i] = equivalent_sparse_rows[reversed_round][i - 1];
                            }
                        }
                    }

                    matrix_type pre_sparse_matrix;
                    state_type first_partial_round_constants;
                    partial_round_constants_type optimized_partial_round_constants;
                    sparse_vectors_type sparse_columns;
                    sparse_first_rows_type sparse_first_rows;
                };

            }    // namespace detail
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON1_OPTIMIZED_CONSTANTS_HPP
