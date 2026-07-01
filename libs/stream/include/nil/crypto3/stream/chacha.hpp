//---------------------------------------------------------------------------//
// Copyright (c) 2019 Mikhail Komarov <nemo@nil.foundation>
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

#ifndef CRYPTO3_STREAM_CHACHA_HPP
#define CRYPTO3_STREAM_CHACHA_HPP

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include <boost/static_assert.hpp>

#include <nil/crypto3/stream/detail/chacha/chacha_functions.hpp>

namespace nil {
    namespace crypto3 {
        namespace stream {
            template<std::size_t IVBits, std::size_t KeyBits, std::size_t Rounds>
            struct chacha_finalizer {
                typedef detail::chacha_functions<Rounds, IVBits, KeyBits> policy_type;

            public:
                constexpr static const std::size_t rounds = policy_type::rounds;

                constexpr static const std::size_t block_bits = policy_type::block_bits;
                constexpr static const std::size_t block_size = policy_type::block_size;
                typedef typename policy_type::block_type block_type;

                constexpr static const std::size_t min_key_schedule_bits = policy_type::min_key_schedule_bits;
                constexpr static const std::size_t min_key_schedule_size = policy_type::min_key_schedule_size;
                typedef typename policy_type::key_schedule_type key_schedule_type;
                typedef typename key_schedule_type::value_type word_type;

                constexpr static const std::size_t iv_bits = policy_type::iv_bits;
                typedef typename policy_type::iv_type iv_type;

                constexpr static const std::size_t min_key_bits = policy_type::min_key_bits;
                constexpr static const std::size_t max_key_bits = policy_type::max_key_bits;
                constexpr static const std::size_t key_bits = policy_type::key_bits;
                typedef typename policy_type::key_type key_type;

                template<typename InputRange, typename OutputRange>
                void process(InputRange &in, OutputRange &out, key_schedule_type &schedule, block_type &block) {
                    for (std::size_t i = 0; i != block_size; ++i) {
                        out[i] = in[i] ^ block[i];
                    }
                }
            };
            /*!
             * @brief DJB's ChaCha (https://cr.yp.to/chacha.html)
             * @tparam Rounds Amount of rounds
             * @ingroup stream
             * @note Currently only 8, 12 or 20 rounds are supported, all others
             * will throw an exception
             */
            template<std::size_t IVBits = 64, std::size_t KeyBits = 128, std::size_t Rounds = 20>
            class chacha {
                typedef detail::chacha_functions<Rounds, IVBits, KeyBits> policy_type;

            public:
                constexpr static const std::size_t rounds = policy_type::rounds;

                constexpr static const std::size_t block_bits = policy_type::block_bits;
                constexpr static const std::size_t block_size = policy_type::block_size;
                typedef typename policy_type::block_type block_type;

                constexpr static const std::size_t min_key_schedule_bits = policy_type::min_key_schedule_bits;
                constexpr static const std::size_t min_key_schedule_size = policy_type::min_key_schedule_size;
                typedef typename policy_type::key_schedule_type key_schedule_type;
                typedef typename key_schedule_type::value_type word_type;

                constexpr static const std::size_t iv_bits = policy_type::iv_bits;
                typedef typename policy_type::iv_type iv_type;

                constexpr static const std::size_t min_key_bits = policy_type::min_key_bits;
                constexpr static const std::size_t max_key_bits = policy_type::max_key_bits;
                constexpr static const std::size_t key_bits = policy_type::key_bits;
                typedef typename policy_type::key_type key_type;

                chacha(block_type &block, key_schedule_type &schedule, const key_type &key,
                       const iv_type &iv = iv_type(), std::uint64_t initial_counter = 0) :
                    block_offset(0), ietf_counter_exhausted(false) {
                    policy_type::schedule_key(schedule, key);
                    policy_type::schedule_iv(schedule, iv, initial_counter);
                    generate_next_block(block, schedule);
                }

                template<typename InputIterator, typename OutputIterator>
                OutputIterator process(InputIterator first, InputIterator last, OutputIterator out,
                                       key_schedule_type &schedule, block_type &block) {
                    while (first != last) {
                        if (block_offset == block_size) {
                            generate_next_block(block, schedule);
                            block_offset = 0;
                        }

                        *out = *first ^ block[block_offset];

                        ++first;
                        ++out;
                        ++block_offset;
                    }

                    return out;
                }

                template<typename InputIterator, typename OutputIterator>
                OutputIterator process_n(InputIterator first, std::size_t length, OutputIterator out,
                                         key_schedule_type &schedule, block_type &block) {
                    while (length != 0) {
                        if (block_offset == block_size) {
                            generate_next_block(block, schedule);
                            block_offset = 0;
                        }

                        *out = *first ^ block[block_offset];

                        ++first;
                        ++out;
                        --length;
                        ++block_offset;
                    }

                    return out;
                }

                template<typename InputRange, typename OutputRange>
                void process(InputRange &in, OutputRange &out, key_schedule_type &schedule, block_type &block) {
                    process_n(in, block_size, out, schedule, block);
                }

                void seek(block_type &block, key_schedule_type &schedule, std::uint64_t offset) {
                    BOOST_STATIC_ASSERT(IVBits == 64 || IVBits == 96);

                    const std::uint64_t counter = offset / block_size;
                    block_offset = offset % block_size;
                    if (IVBits == 96) {
                        if (counter > ietf_counter_limit()) {
                            throw_ietf_counter_exhausted();
                        }
                        schedule[12] = static_cast<word_type>(counter);
                    } else {
                        schedule[12] = static_cast<word_type>(counter);
                        schedule[13] = static_cast<word_type>(counter >> 32);
                    }
                    ietf_counter_exhausted = false;
                    generate_next_block(block, schedule);
                }

            private:
                static constexpr std::uint64_t ietf_counter_limit() {
                    return std::numeric_limits<std::uint32_t>::max();
                }

                static void throw_ietf_counter_exhausted() {
                    throw std::out_of_range("ChaCha20 IETF counter exhausted");
                }

                void generate_next_block(block_type &block, key_schedule_type &schedule) {
                    if (IVBits == 96) {
                        if (ietf_counter_exhausted) {
                            throw_ietf_counter_exhausted();
                        }
                        if (schedule[12] == ietf_counter_limit()) {
                            policy_type::generate_block_without_counter_increment(block, schedule);
                            ietf_counter_exhausted = true;
                            return;
                        }
                    }

                    policy_type::generate_block(block, schedule);
                    ietf_counter_exhausted = false;
                }

                std::size_t block_offset;
                bool ietf_counter_exhausted;
            };

            template<std::size_t KeyBits = 256, std::size_t Rounds = 20>
            using ietf_chacha = chacha<96, KeyBits, Rounds>;

            template<std::size_t KeyBits = 256, std::size_t Rounds = 20>
            using original_chacha = chacha<64, KeyBits, Rounds>;

            using chacha20 = ietf_chacha<256, 20>;
            using original_chacha20 = original_chacha<256, 20>;

            class chacha20_cipher {
            public:
                typedef chacha20 cipher_type;

                constexpr static const std::size_t block_size = cipher_type::block_size;
                typedef typename cipher_type::key_type key_type;
                typedef typename cipher_type::iv_type nonce_type;
                typedef typename cipher_type::iv_type iv_type;

                chacha20_cipher(const key_type &key, const nonce_type &nonce, std::uint32_t initial_counter_value) :
                    initial_counter(initial_counter_value), block(), schedule(),
                    cipher(block, schedule, key, nonce, initial_counter_value) {
                }

                chacha20_cipher(const chacha20_cipher &) = delete;
                chacha20_cipher &operator=(const chacha20_cipher &) = delete;

                template<typename InputIterator, typename OutputIterator>
                OutputIterator process(InputIterator first, InputIterator last, OutputIterator out) {
                    return cipher.process(first, last, out, schedule, block);
                }

                void seek(std::uint64_t byte_offset) {
                    const std::uint64_t initial_offset =
                        static_cast<std::uint64_t>(initial_counter) * block_size;
                    const std::uint64_t max_offset =
                        static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max()) * block_size +
                        (block_size - 1);

                    if (byte_offset > max_offset - initial_offset) {
                        throw std::out_of_range("ChaCha20 IETF counter exhausted");
                    }

                    cipher.seek(block, schedule, initial_offset + byte_offset);
                }

            private:
                std::uint32_t initial_counter;
                typename cipher_type::block_type block;
                typename cipher_type::key_schedule_type schedule;
                cipher_type cipher;
            };
        }    // namespace stream
    }        // namespace crypto3
}    // namespace nil

#endif
