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

#ifndef CRYPTO3_STREAM_CHACHA_IMPL_HPP
#define CRYPTO3_STREAM_CHACHA_IMPL_HPP

#include <nil/crypto3/stream/detail/chacha/chacha_policy.hpp>

#define CHACHA_QUARTER_ROUND(a, b, c, d) \
    do {                                 \
        a += b;                          \
        d ^= a;                          \
        d = policy_type::template rotl<16>(d);    \
        c += d;                          \
        b ^= c;                          \
        b = policy_type::template rotl<12>(b);    \
        a += b;                          \
        d ^= a;                          \
        d = policy_type::template rotl<8>(d);     \
        c += d;                          \
        b ^= c;                          \
        b = policy_type::template rotl<7>(b);     \
    } while (0)

namespace nil {
    namespace crypto3 {
        namespace stream {
            namespace detail {
                template<std::size_t Round, std::size_t IVSize, std::size_t KeyBits>
                struct chacha_impl {
                    typedef chacha_policy<Round, IVSize, KeyBits> policy_type;

                    constexpr static const std::size_t word_bits = policy_type::word_bits;
                    typedef typename policy_type::word_type word_type;

                    constexpr static const std::size_t rounds = policy_type::rounds;

                    constexpr static const std::size_t min_key_schedule_bits = policy_type::key_schedule_bits;
                    constexpr static const std::size_t min_key_schedule_size = policy_type::key_schedule_size;
                    typedef typename policy_type::key_schedule_type key_schedule_type;

                    constexpr static const std::size_t block_bits = policy_type::block_bits;
                    constexpr static const std::size_t block_size = policy_type::block_size;
                    typedef typename policy_type::block_type block_type;

                    static void chacha_x8_original(std::array<std::uint8_t, block_size * 8> &block,
                                                   key_schedule_type &input) {
                        chacha_blocks(block.data(), 8, input, counter_mode::original);
                    }

                    static void chacha_x8_ietf(std::array<std::uint8_t, block_size * 8> &block,
                                               key_schedule_type &input) {
                        chacha_blocks(block.data(), 8, input, counter_mode::ietf);
                    }

                    static void chacha_x8(std::array<std::uint8_t, block_size * 8> &block,
                                          key_schedule_type &input) {
                        BOOST_STATIC_ASSERT(IVSize == 64 || IVSize == 96);
                        if (IVSize == 96) {
                            chacha_x8_ietf(block, input);
                        } else {
                            chacha_x8_original(block, input);
                        }
                    }

                    static void chacha_x4_original(std::array<std::uint8_t, block_size * 4> &block,
                                                   key_schedule_type &input) {
                        chacha_blocks(block.data(), 4, input, counter_mode::original);
                    }

                    static void chacha_x4_ietf(std::array<std::uint8_t, block_size * 4> &block,
                                               key_schedule_type &input) {
                        chacha_blocks(block.data(), 4, input, counter_mode::ietf);
                    }

                    static void chacha_x4(std::array<std::uint8_t, block_size * 4> &block,
                                          key_schedule_type &input) {
                        BOOST_STATIC_ASSERT(IVSize == 64 || IVSize == 96);
                        if (IVSize == 96) {
                            chacha_x4_ietf(block, input);
                        } else {
                            chacha_x4_original(block, input);
                        }
                    }

                    static void chacha_block_original(block_type &block, key_schedule_type &input) {
                        chacha_block(block.data(), input);
                        increment_counter(input, counter_mode::original);
                    }

                    static void chacha_block_ietf(block_type &block, key_schedule_type &input) {
                        chacha_block(block.data(), input);
                        increment_counter(input, counter_mode::ietf);
                    }

                    static void chacha_block(block_type &block, key_schedule_type &input) {
                        BOOST_STATIC_ASSERT(IVSize == 64 || IVSize == 96);
                        if (IVSize == 96) {
                            chacha_block_ietf(block, input);
                        } else {
                            chacha_block_original(block, input);
                        }
                    }

                private:
                    enum class counter_mode { original, ietf };

                    static void chacha_blocks(std::uint8_t *out, std::size_t blocks, key_schedule_type &input,
                                              counter_mode mode) {
                        for (std::size_t i = 0; i != blocks; ++i) {
                            chacha_block(out + block_size * i, input);
                            increment_counter(input, mode);
                        }
                    }

                    static void increment_counter(key_schedule_type &input, counter_mode mode) {
                        ++input[12];
                        if (mode == counter_mode::original && input[12] == 0) {
                            ++input[13];
                        }
                    }

                    static void chacha_block(std::uint8_t *out, const key_schedule_type &input) {
                        word_type x00 = input[0], x01 = input[1], x02 = input[2], x03 = input[3], x04 = input[4],
                                  x05 = input[5], x06 = input[6], x07 = input[7], x08 = input[8], x09 = input[9],
                                  x10 = input[10], x11 = input[11], x12 = input[12], x13 = input[13],
                                  x14 = input[14], x15 = input[15];

                        for (std::size_t r = 0; r != rounds / 2; ++r) {
                            CHACHA_QUARTER_ROUND(x00, x04, x08, x12);
                            CHACHA_QUARTER_ROUND(x01, x05, x09, x13);
                            CHACHA_QUARTER_ROUND(x02, x06, x10, x14);
                            CHACHA_QUARTER_ROUND(x03, x07, x11, x15);

                            CHACHA_QUARTER_ROUND(x00, x05, x10, x15);
                            CHACHA_QUARTER_ROUND(x01, x06, x11, x12);
                            CHACHA_QUARTER_ROUND(x02, x07, x08, x13);
                            CHACHA_QUARTER_ROUND(x03, x04, x09, x14);
                        }

                        boost::endian::store_little_u32(out + 4 * 0, x00 + input[0]);
                        boost::endian::store_little_u32(out + 4 * 1, x01 + input[1]);
                        boost::endian::store_little_u32(out + 4 * 2, x02 + input[2]);
                        boost::endian::store_little_u32(out + 4 * 3, x03 + input[3]);
                        boost::endian::store_little_u32(out + 4 * 4, x04 + input[4]);
                        boost::endian::store_little_u32(out + 4 * 5, x05 + input[5]);
                        boost::endian::store_little_u32(out + 4 * 6, x06 + input[6]);
                        boost::endian::store_little_u32(out + 4 * 7, x07 + input[7]);
                        boost::endian::store_little_u32(out + 4 * 8, x08 + input[8]);
                        boost::endian::store_little_u32(out + 4 * 9, x09 + input[9]);
                        boost::endian::store_little_u32(out + 4 * 10, x10 + input[10]);
                        boost::endian::store_little_u32(out + 4 * 11, x11 + input[11]);
                        boost::endian::store_little_u32(out + 4 * 12, x12 + input[12]);
                        boost::endian::store_little_u32(out + 4 * 13, x13 + input[13]);
                        boost::endian::store_little_u32(out + 4 * 14, x14 + input[14]);
                        boost::endian::store_little_u32(out + 4 * 15, x15 + input[15]);
                    }
                };
            }    // namespace detail
        }        // namespace stream
    }            // namespace crypto3
}    // namespace nil

#undef CHACHA_QUARTER_ROUND
#endif    // CRYPTO3_CHACHA_IMPL_HPP
