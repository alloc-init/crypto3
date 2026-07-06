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

#ifndef CRYPTO3_MAC_POLY1305_HPP
#define CRYPTO3_MAC_POLY1305_HPP

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <iterator>
#include <stdexcept>

#include <boost/multiprecision/cpp_int.hpp>

#include <nil/crypto3/mac/mac_key.hpp>

namespace nil {
    namespace crypto3 {
        namespace mac {
            namespace detail {
                class poly1305_state {
                    typedef boost::multiprecision::cpp_int integer_type;

                public:
                    poly1305_state() : accumulator(0), buffer_size(0) {
                        buffer.fill(0);
                    }

                    template<typename InputIterator>
                    void update(InputIterator first, InputIterator last, const integer_type &r) {
                        while (first != last) {
                            buffer[buffer_size++] = static_cast<std::uint8_t>(*first++);
                            if (buffer_size == buffer.size()) {
                                process_buffer(r);
                            }
                        }
                    }

                    integer_type finalized_accumulator(const integer_type &r) const {
                        poly1305_state finalized = *this;
                        if (finalized.buffer_size != 0) {
                            finalized.process_buffer(r);
                        }
                        return finalized.accumulator;
                    }

                private:
                    void process_buffer(const integer_type &r) {
                        integer_type block_value = 0;
                        for (std::size_t i = 0; i < buffer_size; ++i) {
                            block_value += integer_type(buffer[i]) << (8 * i);
                        }
                        // RFC 8439 processes each limb as little_endian(block || 0x01).
                        // Add that implicit 1 byte numerically instead of materializing it.
                        block_value += integer_type(1) << (8 * buffer_size);

                        accumulator = ((accumulator + block_value) * r) % modulus();
                        buffer_size = 0;
                    }

                    static integer_type modulus() {
                        return (integer_type(1) << 130) - 5;
                    }

                    integer_type accumulator;
                    std::array<std::uint8_t, 16> buffer;
                    std::size_t buffer_size;
                };

                template<typename InputIterator>
                boost::multiprecision::cpp_int load_little_endian_integer(InputIterator first, std::size_t size) {
                    boost::multiprecision::cpp_int result = 0;
                    for (std::size_t i = 0; i < size; ++i, ++first) {
                        result += boost::multiprecision::cpp_int(static_cast<std::uint8_t>(*first)) << (8 * i);
                    }
                    return result;
                }

                inline std::array<std::uint8_t, 16> store_little_endian_128(boost::multiprecision::cpp_int value) {
                    std::array<std::uint8_t, 16> out = {0};
                    for (std::size_t i = 0; i < out.size(); ++i) {
                        out[i] = static_cast<std::uint8_t>(value & 0xff);
                        value >>= 8;
                    }
                    return out;
                }

                inline void clamp_poly1305_r(std::array<std::uint8_t, 16> &r) {
                    r[3] &= 15;
                    r[7] &= 15;
                    r[11] &= 15;
                    r[15] &= 15;

                    r[4] &= 252;
                    r[8] &= 252;
                    r[12] &= 252;
                }
            }    // namespace detail

            struct poly1305 {
                constexpr static const std::size_t key_bits = 256;
                constexpr static const std::size_t key_octets = key_bits / CHAR_BIT;
                typedef std::array<std::uint8_t, key_octets> key_type;

                constexpr static const std::size_t block_bits = 128;
                constexpr static const std::size_t block_octets = block_bits / CHAR_BIT;
                constexpr static const std::size_t block_words = block_octets;
                typedef std::array<std::uint8_t, block_octets> block_type;

                constexpr static const std::size_t digest_bits = 128;
                constexpr static const std::size_t digest_octets = digest_bits / CHAR_BIT;
                typedef std::array<std::uint8_t, digest_octets> digest_type;
            };

            template<>
            struct mac_key<poly1305> {
                typedef poly1305 policy_type;
                typedef typename policy_type::digest_type digest_type;
                typedef detail::poly1305_state accumulator_type;

                template<typename KeyRange>
                explicit mac_key(const KeyRange &key) : r(0), s(0) {
                    process_key(key.cbegin(), key.cend());
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
                    acc.update(first, last, r);
                }

                inline digest_type compute(accumulator_type &acc) const {
                    boost::multiprecision::cpp_int tag = acc.finalized_accumulator(r) + s;
                    tag &= (boost::multiprecision::cpp_int(1) << policy_type::digest_bits) - 1;
                    return detail::store_little_endian_128(tag);
                }

            private:
                template<typename InputIterator>
                void process_key(InputIterator first, InputIterator last) {
                    if (std::distance(first, last) !=
                        static_cast<typename std::iterator_traits<InputIterator>::difference_type>(
                            policy_type::key_octets)) {
                        throw std::invalid_argument("Poly1305 key must be exactly 32 bytes");
                    }

                    std::array<std::uint8_t, 16> r_bytes = {0};
                    std::array<std::uint8_t, 16> s_bytes = {0};

                    InputIterator current = first;
                    for (std::size_t i = 0; i < r_bytes.size(); ++i, ++current) {
                        r_bytes[i] = static_cast<std::uint8_t>(*current);
                    }
                    for (std::size_t i = 0; i < s_bytes.size(); ++i, ++current) {
                        s_bytes[i] = static_cast<std::uint8_t>(*current);
                    }

                    detail::clamp_poly1305_r(r_bytes);
                    r = detail::load_little_endian_integer(r_bytes.cbegin(), r_bytes.size());
                    s = detail::load_little_endian_integer(s_bytes.cbegin(), s_bytes.size());
                }

                boost::multiprecision::cpp_int r;
                boost::multiprecision::cpp_int s;
            };
        }    // namespace mac
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MAC_POLY1305_HPP
