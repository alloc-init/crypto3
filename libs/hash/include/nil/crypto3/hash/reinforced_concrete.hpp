//---------------------------------------------------------------------------//
// Copyright (c) 2023 Mikhail Komarov <nemo@nil.foundation>
//
// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_REINFORCED_COCNRETE_HPP
#define CRYPTO3_HASH_REINFORCED_COCNRETE_HPP

#include <nil/crypto3/hash/detail/reinforced_concrete/reinforced_concrete_policy.hpp>
#include <nil/crypto3/hash/detail/reinforced_concrete/reinforced_concrete_functions.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            template<typename FieldType>
            struct reinforced_concrete_compressor {
                typedef detail::reinforced_concrete_functions<FieldType> policy_type;
                typedef typename policy_type::element_type element_type;

                constexpr static const std::size_t word_bits = policy_type::word_bits;
                typedef typename policy_type::word_type word_type;

                constexpr static const std::size_t digest_bits = policy_type::digest_bits;
                typedef typename policy_type::digest_type digest_type;

                constexpr static const std::size_t block_words = policy_type::block_words;
                constexpr static const std::size_t block_bits = policy_type::block_bits;
                typedef typename policy_type::block_type block_type;

                constexpr static const std::size_t state_words = policy_type::state_words;
                constexpr static const std::size_t state_bits = policy_type::state_bits;
                typedef typename policy_type::state_type state_type;

                static inline void process_block(state_type &state, block_type &block) {
                    for (int i = 0; i < block_words; ++i) {
                        state[i] ^= block[i];
                    }

                    policy_type::permute(state);
                }
            };
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_REINFORCED_COCNRETE_HPP