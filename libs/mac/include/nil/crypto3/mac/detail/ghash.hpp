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

#ifndef CRYPTO3_MAC_DETAIL_GHASH_HPP
#define CRYPTO3_MAC_DETAIL_GHASH_HPP

#include <array>
#include <cstddef>
#include <cstdint>

namespace nil {
    namespace crypto3 {
        namespace mac {
            namespace detail {
                namespace ghash {
                    constexpr static const std::size_t block_octets = 16;
                    typedef std::array<std::uint8_t, block_octets> block_type;

                    inline block_type zero_block() {
                        block_type block = {0};
                        return block;
                    }

                    inline void xor_block(block_type &out, const block_type &in) {
                        for (std::size_t i = 0; i < out.size(); ++i) {
                            out[i] ^= in[i];
                        }
                    }

                    // GHASH scans bits most-significant first. Return 0xff when the selected bit is set and 0x00
                    // otherwise, so callers can apply conditional XORs without a data-dependent branch.
                    inline std::uint8_t bit_mask(const block_type &block, std::size_t bit_index) {
                        const std::uint8_t bit =
                            static_cast<std::uint8_t>((block[bit_index / 8] >> (7 - (bit_index % 8))) & 1U);
                        return static_cast<std::uint8_t>(0U - bit);
                    }

                    inline void shift_right_one(block_type &block) {
                        std::uint8_t carry = 0;
                        for (std::size_t i = 0; i < block.size(); ++i) {
                            const std::uint8_t next_carry = static_cast<std::uint8_t>(block[i] & 1U);
                            block[i] = static_cast<std::uint8_t>((block[i] >> 1) | (carry << 7));
                            carry = next_carry;
                        }
                    }

                    inline block_type multiply(const block_type &x, const block_type &y) {
                        block_type z = zero_block();
                        block_type v = y;

                        for (std::size_t i = 0; i < 128; ++i) {
                            const std::uint8_t x_mask = bit_mask(x, i);
                            for (std::size_t j = 0; j < z.size(); ++j) {
                                z[j] ^= static_cast<std::uint8_t>(v[j] & x_mask);
                            }

                            const std::uint8_t reduction_mask = static_cast<std::uint8_t>(0U - (v[15] & 1U));
                            shift_right_one(v);
                            v[0] ^= static_cast<std::uint8_t>(0xe1U & reduction_mask);
                        }

                        return z;
                    }

                    inline void store_big_endian_64(block_type &out, std::size_t offset, std::uint64_t value) {
                        for (std::size_t i = 0; i < 8; ++i) {
                            out[offset + 7 - i] = static_cast<std::uint8_t>(value >> (8 * i));
                        }
                    }

                    // GHASH always finishes with one block containing the bit lengths of the authenticated data and
                    // ciphertext. These lengths are mixed into the hash.
                    inline block_type length_block(std::uint64_t associated_data_bits, std::uint64_t ciphertext_bits) {
                        block_type block = zero_block();
                        store_big_endian_64(block, 0, associated_data_bits);
                        store_big_endian_64(block, 8, ciphertext_bits);
                        return block;
                    }

                    // GHASH evaluates a polynomial over GF(2^128) using the hash subkey H. Each 16-byte input block
                    // is folded into the running value as:
                    //
                    //     Y_0 = 0
                    //     Y_i = (Y_{i-1} xor X_i) * H
                    //
                    // Callers feed AAD blocks first, then ciphertext blocks, padding each stream with zeros to a
                    // 16-byte boundary. The final block contains the bit lengths of those two streams.
                    class state {
                    public:
                        explicit state(const block_type &hash_subkey) :
                            hash_subkey(hash_subkey), value(zero_block()), buffer(zero_block()), buffer_size(0) {
                        }

                        template<typename InputIterator>
                        void update(InputIterator first, InputIterator last) {
                            while (first != last) {
                                buffer[buffer_size++] = static_cast<std::uint8_t>(*first++);
                                if (buffer_size == buffer.size()) {
                                    process_buffer();
                                }
                            }
                        }

                        void update_octet(std::uint8_t octet) {
                            buffer[buffer_size++] = octet;
                            if (buffer_size == buffer.size()) {
                                process_buffer();
                            }
                        }

                        void pad() {
                            if (buffer_size != 0) {
                                process_buffer();
                            }
                        }

                        void update_block(const block_type &block) {
                            xor_block(value, block);
                            value = multiply(value, hash_subkey);
                        }

                        // Add the mandatory final length block:
                        //     uint64_be(len(AAD) in bits) || uint64_be(len(ciphertext) in bits)
                        void update_lengths(std::uint64_t associated_data_bits, std::uint64_t ciphertext_bits) {
                            update_block(length_block(associated_data_bits, ciphertext_bits));
                        }

                        block_type digest() const {
                            return value;
                        }

                    private:
                        void process_buffer() {
                            update_block(buffer);
                            buffer = zero_block();
                            buffer_size = 0;
                        }

                        block_type hash_subkey;
                        // running hash value
                        block_type value;
                        // Staging area for ordinary input bytes until they form a 16-byte GHASH block. The mandatory
                        // final length block is created separately by length_block(...) and passed to
                        // update_block(...).
                        block_type buffer;
                        std::size_t buffer_size;
                    };
                }    // namespace ghash
            }    // namespace detail
        }    // namespace mac
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MAC_DETAIL_GHASH_HPP
