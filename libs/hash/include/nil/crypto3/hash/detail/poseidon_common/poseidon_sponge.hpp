//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_DETAIL_POSEIDON_COMMON_SPONGE_HPP
#define CRYPTO3_HASH_DETAIL_POSEIDON_COMMON_SPONGE_HPP

#include <algorithm>
#include <cstddef>
#include <type_traits>

#include <boost/assert.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {

                enum class poseidon_sponge_absorb_mode { overwrite, add };

                enum class poseidon_sponge_padding_mode { pad10, padding_free };

                /*!
                 * @brief Mode-parameterized sponge for Poseidon permutations.
                 *
                 * Absorb modes:
                 * - overwrite: rate cells are replaced by input cells.
                 * - add: input cells are added into the existing rate cells. This is the usual field
                 *   analogue of XOR-mode binary sponges.
                 *
                 * Padding modes:
                 * - pad10: variable-length-safe two-case padding. A partial final block gets a rate
                 *   sentinel and zero suffix in overwrite mode. A full final block increments the first
                 *   capacity cell before the final permutation.
                 * - padding_free: fixed-length-only mode. Full blocks are permuted as-is; a partial final
                 *   block overwrites/adds only its present cells and leaves the rest of the rate unchanged.
                 *
                 * A streaming Pad10 implementation must delay the most recent full block. Otherwise a
                 * message whose length is exactly a multiple of the rate would be permuted before the
                 * capacity-domain padding marker is applied.
                 *
                 * The digest output is intentionally limited to OUT <= RATE: digest() copies one squeezed
                 * rate block. Supporting larger outputs requires a multi-squeeze digest path.
                 */
                template<typename PolicyType,
                         typename PermutationType,
                         poseidon_sponge_absorb_mode AbsorbMode = poseidon_sponge_absorb_mode::overwrite,
                         poseidon_sponge_padding_mode PaddingMode = poseidon_sponge_padding_mode::pad10>
                class poseidon_sponge_construction {
                    using policy_type = PolicyType;
                    using permutation_type = PermutationType;

                public:
                    // one field element
                    using word_type = typename policy_type::word_type;
                    // full Poseidon state: rate + capacity
                    using state_type = typename policy_type::state_type;
                    // one absorb block, exactly rate field elements
                    using block_type = typename policy_type::block_type;
                    using digest_type = typename policy_type::digest_type;

                    constexpr static const std::size_t word_bits = policy_type::word_bits;
                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t block_words = policy_type::block_words;
                    constexpr static const std::size_t digest_words = policy_type::digest_words;
                    constexpr static const std::size_t state_bits = policy_type::state_bits;
                    constexpr static const std::size_t block_bits = policy_type::block_bits;

                    static_assert(block_words > 0, "Poseidon sponge requires a positive rate.");
                    static_assert(block_words < state_words, "Poseidon sponge requires at least one capacity element.");
                    static_assert(digest_words <= block_words,
                                  "Poseidon sponge intentionally supports only OUT <= RATE; "
                                  "extend digest() to squeeze multiple blocks before raising this limit.");

                    poseidon_sponge_construction() {
                        reset();
                    }

                    void absorb(const block_type &block) {
                        BOOST_ASSERT_MSG(!finalized_, "Cannot absorb after Poseidon sponge finalization.");

                        flush_pending_full_block();
                        pending_full_block_ = block;
                        has_pending_full_block_ = true;
                    }

                    void absorb_with_padding(const block_type &block = block_type(),
                                             const std::size_t last_block_words_filled = 0) {
                        BOOST_ASSERT_MSG(!finalized_, "Poseidon sponge was already finalized.");
                        BOOST_ASSERT_MSG(last_block_words_filled <= block_words,
                                         "Final Poseidon sponge block exceeds the rate.");

                        if constexpr (PaddingMode == poseidon_sponge_padding_mode::pad10) {
                            finalize_pad10(block, last_block_words_filled);
                        } else {
                            finalize_padding_free(block, last_block_words_filled);
                        }
                    }

                    block_type squeeze() {
                        finalize_if_needed();

                        block_type block;
                        std::copy(state_.begin(), state_.begin() + block_words, block.begin());
                        permutation_type::permute(state_);
                        return block;
                    }

                    digest_type digest() {
                        const block_type block = squeeze();
                        return make_digest(block);
                    }

                    void reset() {
                        state_ = policy_type::iv_generator::generate();
                        pending_full_block_ = block_type();
                        has_pending_full_block_ = false;
                        finalized_ = false;
                    }

                    const state_type &state() const {
                        return state_;
                    }

                private:
                    static constexpr word_type zero() {
                        return word_type(0u);
                    }

                    static constexpr word_type one() {
                        return word_type(1u);
                    }

                    static digest_type make_digest(const block_type &block) {
                        if constexpr (digest_words == 1 && std::is_same<digest_type, word_type>::value) {
                            return block[0];
                        } else {
                            digest_type digest;
                            for (std::size_t i = 0; i < digest_words; ++i) {
                                digest[i] = block[i];
                            }
                            return digest;
                        }
                    }

                    void flush_pending_full_block() {
                        if (!has_pending_full_block_) {
                            return;
                        }

                        absorb_rate_block(pending_full_block_);
                        permutation_type::permute(state_);
                        has_pending_full_block_ = false;
                    }

                    void finalize_pad10(const block_type &block, std::size_t words_filled) {
                        if (words_filled == block_words) {
                            flush_pending_full_block();
                            absorb_pad10_final_full_block(block);
                            return;
                        }

                        if (words_filled == 0 && has_pending_full_block_) {
                            const block_type final_block = pending_full_block_;
                            has_pending_full_block_ = false;
                            absorb_pad10_final_full_block(final_block);
                            return;
                        }

                        flush_pending_full_block();
                        absorb_pad10_final_partial_block(block, words_filled);
                    }

                    void finalize_padding_free(const block_type &block, std::size_t words_filled) {
                        if (words_filled == block_words) {
                            flush_pending_full_block();
                            absorb_padding_free_final_full_block(block);
                            return;
                        }

                        if (words_filled == 0) {
                            flush_pending_full_block();
                            finalized_ = true;
                            return;
                        }

                        flush_pending_full_block();
                        absorb_padding_free_final_partial_block(block, words_filled);
                    }

                    void absorb_pad10_final_full_block(const block_type &block) {
                        absorb_rate_block(block);
                        state_[block_words] += one();
                        permutation_type::permute(state_);
                        finalized_ = true;
                    }

                    void absorb_pad10_final_partial_block(const block_type &block, std::size_t words_filled) {
                        for (std::size_t i = 0; i < words_filled; ++i) {
                            absorb_rate_word(i, block[i]);
                        }

                        if constexpr (AbsorbMode == poseidon_sponge_absorb_mode::overwrite) {
                            state_[words_filled] = one();
                            for (std::size_t i = words_filled + 1; i < block_words; ++i) {
                                state_[i] = zero();
                            }
                        } else {
                            state_[words_filled] += one();
                        }

                        permutation_type::permute(state_);
                        finalized_ = true;
                    }

                    void absorb_padding_free_final_full_block(const block_type &block) {
                        absorb_rate_block(block);
                        permutation_type::permute(state_);
                        finalized_ = true;
                    }

                    void absorb_padding_free_final_partial_block(const block_type &block, std::size_t words_filled) {
                        for (std::size_t i = 0; i < words_filled; ++i) {
                            absorb_rate_word(i, block[i]);
                        }
                        permutation_type::permute(state_);
                        finalized_ = true;
                    }

                    void absorb_rate_block(const block_type &block) {
                        for (std::size_t i = 0; i < block_words; ++i) {
                            absorb_rate_word(i, block[i]);
                        }
                    }

                    void absorb_rate_word(std::size_t index, const word_type &word) {
                        if constexpr (AbsorbMode == poseidon_sponge_absorb_mode::overwrite) {
                            state_[index] = word;
                        } else {
                            state_[index] += word;
                        }
                    }

                    void finalize_if_needed() {
                        // finalize as if the final block is empty
                        if (!finalized_) {
                            absorb_with_padding(block_type(), 0);
                        }
                    }

                    state_type state_;
                    block_type pending_full_block_;
                    bool has_pending_full_block_;
                    bool finalized_;
                };

                template<typename PolicyType, typename PermutationType>
                using poseidon_pad10_sponge_construction =
                    poseidon_sponge_construction<PolicyType,
                                                 PermutationType,
                                                 poseidon_sponge_absorb_mode::overwrite,
                                                 poseidon_sponge_padding_mode::pad10>;

                template<typename PolicyType, typename PermutationType>
                using poseidon_padding_free_sponge_construction =
                    poseidon_sponge_construction<PolicyType,
                                                 PermutationType,
                                                 poseidon_sponge_absorb_mode::overwrite,
                                                 poseidon_sponge_padding_mode::padding_free>;

            }    // namespace detail
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_DETAIL_POSEIDON_COMMON_SPONGE_HPP
