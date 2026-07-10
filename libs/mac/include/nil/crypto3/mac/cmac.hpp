//---------------------------------------------------------------------------//
// Copyright (c) 2026 Alloc Init
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

#ifndef CRYPTO3_MAC_CMAC_HPP
#define CRYPTO3_MAC_CMAC_HPP

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <type_traits>

#include <nil/crypto3/mac/mac_key.hpp>

namespace nil {
    namespace crypto3 {
        namespace mac {
            namespace detail {
                template<std::size_t BlockBits>
                struct cmac_polynomial;

                template<>
                struct cmac_polynomial<64> {
                    constexpr static const std::uint8_t value = 0x1b;
                };

                template<>
                struct cmac_polynomial<128> {
                    constexpr static const std::uint8_t value = 0x87;
                };

                template<typename BlockCipher>
                class cmac_state {
                    typedef typename BlockCipher::block_type block_type;
                    constexpr static const std::size_t block_octets = BlockCipher::block_bits / CHAR_BIT;

                public:
                    cmac_state() : position(0) {
                        chaining_value.fill(0);
                        buffer.fill(0);
                    }

                    template<typename InputIterator>
                    void update(InputIterator first, InputIterator last, const BlockCipher &cipher) {
                        while (first != last) {
                            if (position == block_octets) {
                                process_buffer(cipher);
                            }

                            buffer[position++] = static_cast<std::uint8_t>(*first++);
                        }
                    }

                    block_type finalize(const BlockCipher &cipher, const block_type &k1, const block_type &k2) const {
                        block_type final_block = buffer;

                        if (position == block_octets) {
                            xor_block(final_block, k1);
                        } else {
                            final_block[position] = 0x80;
                            std::fill(final_block.begin() + position + 1, final_block.end(), 0);
                            xor_block(final_block, k2);
                        }

                        xor_block(final_block, chaining_value);
                        return cipher.encrypt(final_block);
                    }

                private:
                    void process_buffer(const BlockCipher &cipher) {
                        xor_block(buffer, chaining_value);
                        chaining_value = cipher.encrypt(buffer);
                        buffer.fill(0);
                        position = 0;
                    }

                    static void xor_block(block_type &out, const block_type &in) {
                        for (std::size_t i = 0; i < out.size(); ++i) {
                            out[i] ^= in[i];
                        }
                    }

                    block_type chaining_value;
                    block_type buffer;
                    std::size_t position;
                };
            }    // namespace detail

            /*!
             * @brief CMAC, also known as OMAC1.
             * @tparam BlockCipher block cipher with byte-oriented block and key types.
             * @ingroup mac
             */
            template<typename BlockCipher>
            struct cmac {
                typedef BlockCipher cipher_type;

                constexpr static const std::size_t block_bits = cipher_type::block_bits;
                constexpr static const std::size_t block_octets = block_bits / CHAR_BIT;
                constexpr static const std::size_t block_words = cipher_type::block_words;
                typedef typename cipher_type::block_type block_type;

                constexpr static const std::size_t key_bits = cipher_type::key_bits;
                constexpr static const std::size_t key_octets = key_bits / CHAR_BIT;
                typedef typename cipher_type::key_type key_type;

                constexpr static const std::size_t digest_bits = block_bits;
                constexpr static const std::size_t digest_octets = block_octets;
                typedef block_type digest_type;
            };

            template<typename BlockCipher>
            using omac1 = cmac<BlockCipher>;

            template<typename BlockCipher>
            struct mac_key<cmac<BlockCipher>> {
                typedef cmac<BlockCipher> policy_type;
                typedef typename policy_type::cipher_type cipher_type;
                typedef typename policy_type::key_type key_type;
                typedef typename policy_type::block_type block_type;
                typedef typename policy_type::digest_type digest_type;
                typedef detail::cmac_state<cipher_type> accumulator_type;

                template<typename KeyRange>
                explicit mac_key(const KeyRange &key) :
                    cipher(process_key(key.cbegin(), key.cend())), k1(generate_subkey(cipher.encrypt(zero_block()))),
                    k2(generate_subkey(k1)) {
                }

                inline void init_accumulator(accumulator_type &acc) const {
                    acc = accumulator_type();
                }

                template<typename InputRange>
                inline void update(accumulator_type &acc, const InputRange &range) const {
                    update(acc, range.cbegin(), range.cend());
                }

                template<typename InputIterator>
                inline void update(accumulator_type &acc, InputIterator first, InputIterator last) const {
                    acc.update(first, last, cipher);
                }

                inline digest_type compute(const accumulator_type &acc) const {
                    return acc.finalize(cipher, k1, k2);
                }

            private:
                static_assert(policy_type::block_bits == 64 || policy_type::block_bits == 128,
                              "CMAC supports 64-bit and 128-bit block ciphers");
                static_assert(std::is_same<typename block_type::value_type, std::uint8_t>::value,
                              "CMAC currently requires byte-oriented block cipher blocks");
                static_assert(std::is_same<typename key_type::value_type, std::uint8_t>::value,
                              "CMAC currently requires byte-oriented block cipher keys");

                template<typename InputIterator>
                static key_type process_key(InputIterator first, InputIterator last) {
                    if (std::distance(first, last) !=
                        static_cast<typename std::iterator_traits<InputIterator>::difference_type>(
                            policy_type::key_octets)) {
                        throw std::invalid_argument("CMAC key size does not match the block cipher key size");
                    }

                    key_type key = {0};
                    std::copy(first, last, key.begin());
                    return key;
                }

                static block_type zero_block() {
                    block_type block = {0};
                    return block;
                }

                static block_type generate_subkey(const block_type &input) {
                    block_type out = {0};
                    std::uint8_t carry = 0;

                    for (std::size_t i = policy_type::block_octets; i-- > 0;) {
                        const std::uint8_t next_carry = input[i] >> 7;
                        out[i] = static_cast<std::uint8_t>((input[i] << 1) | carry);
                        carry = next_carry;
                    }

                    if (carry != 0) {
                        out[policy_type::block_octets - 1] ^= detail::cmac_polynomial<policy_type::block_bits>::value;
                    }

                    return out;
                }

                cipher_type cipher;
                block_type k1;
                block_type k2;
            };
        }    // namespace mac
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MAC_CMAC_HPP
