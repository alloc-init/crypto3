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
                inline std::uint64_t load_little_endian_64(const std::uint8_t *in) {
                    std::uint64_t result = 0;
                    for (std::size_t i = 0; i < 8; ++i) {
                        result |= std::uint64_t(in[i]) << (8 * i);
                    }
                    return result;
                }

                inline void store_little_endian_64(std::uint8_t *out, std::uint64_t value) {
                    for (std::size_t i = 0; i < 8; ++i) {
                        out[i] = static_cast<std::uint8_t>(value >> (8 * i));
                    }
                }

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

                struct poly1305_reference_key_schedule {
                    boost::multiprecision::cpp_int r;
                    boost::multiprecision::cpp_int s;
                };

                class poly1305_reference_state {
                    typedef boost::multiprecision::cpp_int integer_type;

                public:
                    poly1305_reference_state() : accumulator(0), buffer_size(0) {
                        buffer.fill(0);
                    }

                    template<typename InputIterator>
                    void update(InputIterator first, InputIterator last, const poly1305_reference_key_schedule &key) {
                        while (first != last) {
                            buffer[buffer_size++] = static_cast<std::uint8_t>(*first++);
                            if (buffer_size == buffer.size()) {
                                process_buffer(key.r);
                            }
                        }
                    }

                    integer_type finalized_accumulator(const poly1305_reference_key_schedule &key) const {
                        poly1305_reference_state finalized = *this;
                        if (finalized.buffer_size != 0) {
                            finalized.process_buffer(key.r);
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

                struct poly1305_reference_backend {
                    typedef poly1305_reference_key_schedule key_schedule_type;
                    typedef poly1305_reference_state accumulator_type;

                    static key_schedule_type process_key(const std::array<std::uint8_t, 32> &key) {
                        std::array<std::uint8_t, 16> r_bytes = {0};
                        std::array<std::uint8_t, 16> s_bytes = {0};

                        std::copy(key.begin(), key.begin() + r_bytes.size(), r_bytes.begin());
                        std::copy(key.begin() + r_bytes.size(), key.end(), s_bytes.begin());

                        clamp_poly1305_r(r_bytes);
                        return {load_little_endian_integer(r_bytes.cbegin(), r_bytes.size()),
                                load_little_endian_integer(s_bytes.cbegin(), s_bytes.size())};
                    }

                    template<typename InputIterator>
                    static void update(accumulator_type &acc,
                                       InputIterator first,
                                       InputIterator last,
                                       const key_schedule_type &key) {
                        acc.update(first, last, key);
                    }

                    static std::array<std::uint8_t, 16> compute(const accumulator_type &acc,
                                                                const key_schedule_type &key) {
                        boost::multiprecision::cpp_int tag = acc.finalized_accumulator(key) + key.s;
                        tag &= (boost::multiprecision::cpp_int(1) << 128) - 1;
                        return store_little_endian_128(tag);
                    }
                };

#if defined(__SIZEOF_INT128__)
                typedef unsigned __int128 poly1305_uint128_type;

                struct poly1305_optimized_key_schedule {
                    std::uint64_t r0;
                    std::uint64_t r1;
                    std::uint64_t r2;
                    std::uint64_t s1;
                    std::uint64_t s2;
                    std::uint64_t pad0;
                    std::uint64_t pad1;
                };

                class poly1305_optimized_state {
                public:
                    poly1305_optimized_state() : h0(0), h1(0), h2(0), buffer_size(0) {
                        buffer.fill(0);
                    }

                    template<typename InputIterator>
                    void update(InputIterator first, InputIterator last, const poly1305_optimized_key_schedule &key) {
                        while (first != last) {
                            buffer[buffer_size++] = static_cast<std::uint8_t>(*first++);
                            if (buffer_size == buffer.size()) {
                                process_block(buffer.data(), key, false);
                                buffer_size = 0;
                            }
                        }
                    }

                    std::array<std::uint8_t, 16> finalize(const poly1305_optimized_key_schedule &key) const {
                        poly1305_optimized_state finalized = *this;
                        if (finalized.buffer_size != 0) {
                            finalized.buffer[finalized.buffer_size] = 1;
                            std::fill(finalized.buffer.begin() + finalized.buffer_size + 1, finalized.buffer.end(), 0);
                            finalized.process_block(finalized.buffer.data(), key, true);
                        }

                        return finalized.finish(key);
                    }

                private:
                    void process_block(const std::uint8_t *block,
                                       const poly1305_optimized_key_schedule &key,
                                       bool partial) {
                        constexpr std::uint64_t mask44 = 0xfffffffffff;
                        constexpr std::uint64_t mask42 = 0x3ffffffffff;

                        const std::uint64_t t0 = load_little_endian_64(block);
                        const std::uint64_t t1 = load_little_endian_64(block + 8);
                        const std::uint64_t hibit = partial ? 0 : (std::uint64_t(1) << 40);

                        h0 += t0 & mask44;
                        h1 += ((t0 >> 44) | (t1 << 20)) & mask44;
                        h2 += ((t1 >> 24) & mask42) | hibit;

                        poly1305_uint128_type d0 = poly1305_uint128_type(h0) * key.r0 +
                                                   poly1305_uint128_type(h1) * key.s2 +
                                                   poly1305_uint128_type(h2) * key.s1;
                        poly1305_uint128_type d1 = poly1305_uint128_type(h0) * key.r1 +
                                                   poly1305_uint128_type(h1) * key.r0 +
                                                   poly1305_uint128_type(h2) * key.s2;
                        poly1305_uint128_type d2 = poly1305_uint128_type(h0) * key.r2 +
                                                   poly1305_uint128_type(h1) * key.r1 +
                                                   poly1305_uint128_type(h2) * key.r0;

                        std::uint64_t carry = static_cast<std::uint64_t>(d0 >> 44);
                        h0 = static_cast<std::uint64_t>(d0) & mask44;
                        d1 += carry;

                        carry = static_cast<std::uint64_t>(d1 >> 44);
                        h1 = static_cast<std::uint64_t>(d1) & mask44;
                        d2 += carry;

                        carry = static_cast<std::uint64_t>(d2 >> 42);
                        h2 = static_cast<std::uint64_t>(d2) & mask42;
                        h0 += carry * 5;

                        carry = h0 >> 44;
                        h0 &= mask44;
                        h1 += carry;
                    }

                    std::array<std::uint8_t, 16> finish(const poly1305_optimized_key_schedule &key) {
                        constexpr std::uint64_t mask44 = 0xfffffffffff;
                        constexpr std::uint64_t mask42 = 0x3ffffffffff;

                        std::uint64_t carry = h1 >> 44;
                        h1 &= mask44;
                        h2 += carry;

                        carry = h2 >> 42;
                        h2 &= mask42;
                        h0 += carry * 5;

                        carry = h0 >> 44;
                        h0 &= mask44;
                        h1 += carry;

                        carry = h1 >> 44;
                        h1 &= mask44;
                        h2 += carry;

                        carry = h2 >> 42;
                        h2 &= mask42;
                        h0 += carry * 5;

                        carry = h0 >> 44;
                        h0 &= mask44;
                        h1 += carry;

                        std::uint64_t g0 = h0 + 5;
                        carry = g0 >> 44;
                        g0 &= mask44;

                        std::uint64_t g1 = h1 + carry;
                        carry = g1 >> 44;
                        g1 &= mask44;

                        std::uint64_t g2 = h2 + carry - (std::uint64_t(1) << 42);

                        std::uint64_t mask = (g2 >> 63) - 1;
                        g0 &= mask;
                        g1 &= mask;
                        g2 &= mask;

                        mask = ~mask;
                        h0 = (h0 & mask) | g0;
                        h1 = (h1 & mask) | g1;
                        h2 = (h2 & mask) | g2;

                        h0 += key.pad0 & mask44;
                        carry = h0 >> 44;
                        h0 &= mask44;

                        h1 += (((key.pad0 >> 44) | (key.pad1 << 20)) & mask44) + carry;
                        carry = h1 >> 44;
                        h1 &= mask44;

                        h2 += ((key.pad1 >> 24) & mask42) + carry;
                        h2 &= mask42;

                        const std::uint64_t out0 = h0 | (h1 << 44);
                        const std::uint64_t out1 = (h1 >> 20) | (h2 << 24);

                        std::array<std::uint8_t, 16> out = {0};
                        store_little_endian_64(out.data(), out0);
                        store_little_endian_64(out.data() + 8, out1);

                        return out;
                    }

                    std::uint64_t h0;
                    std::uint64_t h1;
                    std::uint64_t h2;
                    std::array<std::uint8_t, 16> buffer;
                    std::size_t buffer_size;
                };

                struct poly1305_optimized_backend {
                    typedef poly1305_optimized_key_schedule key_schedule_type;
                    typedef poly1305_optimized_state accumulator_type;

                    static key_schedule_type process_key(const std::array<std::uint8_t, 32> &key) {
                        const std::uint64_t t0 = load_little_endian_64(key.data());
                        const std::uint64_t t1 = load_little_endian_64(key.data() + 8);

                        key_schedule_type schedule = {t0 & 0xffc0fffffff,
                                                      ((t0 >> 44) | (t1 << 20)) & 0xfffffc0ffff,
                                                      (t1 >> 24) & 0x00ffffffc0f,
                                                      0,
                                                      0,
                                                      load_little_endian_64(key.data() + 16),
                                                      load_little_endian_64(key.data() + 24)};

                        schedule.s1 = schedule.r1 * (5 << 2);
                        schedule.s2 = schedule.r2 * (5 << 2);
                        return schedule;
                    }

                    template<typename InputIterator>
                    static void update(accumulator_type &acc,
                                       InputIterator first,
                                       InputIterator last,
                                       const key_schedule_type &key) {
                        acc.update(first, last, key);
                    }

                    static std::array<std::uint8_t, 16> compute(const accumulator_type &acc,
                                                                const key_schedule_type &key) {
                        return acc.finalize(key);
                    }
                };

                typedef poly1305_optimized_backend poly1305_default_backend;
#else
                typedef poly1305_reference_backend poly1305_default_backend;
#endif
            }    // namespace detail

            template<typename Backend = detail::poly1305_default_backend>
            struct basic_poly1305 {
                typedef Backend backend_type;

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

            typedef basic_poly1305<> poly1305;
            typedef basic_poly1305<detail::poly1305_reference_backend> poly1305_reference;

#if defined(__SIZEOF_INT128__)
            typedef basic_poly1305<detail::poly1305_optimized_backend> poly1305_optimized;
#endif

            template<typename Backend>
            struct mac_key<basic_poly1305<Backend>> {
                typedef basic_poly1305<Backend> policy_type;
                typedef typename policy_type::backend_type backend_type;
                typedef typename policy_type::digest_type digest_type;
                typedef typename backend_type::accumulator_type accumulator_type;
                typedef typename backend_type::key_schedule_type key_schedule_type;

                template<typename KeyRange>
                explicit mac_key(const KeyRange &key) : schedule(process_key(key.cbegin(), key.cend())) {
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
                    backend_type::update(acc, first, last, schedule);
                }

                inline digest_type compute(accumulator_type &acc) const {
                    return backend_type::compute(acc, schedule);
                }

            private:
                template<typename InputIterator>
                static key_schedule_type process_key(InputIterator first, InputIterator last) {
                    if (std::distance(first, last) !=
                        static_cast<typename std::iterator_traits<InputIterator>::difference_type>(
                            policy_type::key_octets)) {
                        throw std::invalid_argument("Poly1305 key must be exactly 32 bytes");
                    }

                    typename policy_type::key_type key_bytes = {0};
                    for (std::size_t i = 0; i < key_bytes.size(); ++i, ++first) {
                        key_bytes[i] = static_cast<std::uint8_t>(*first);
                    }

                    return backend_type::process_key(key_bytes);
                }

                key_schedule_type schedule;
            };
        }    // namespace mac
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MAC_POLY1305_HPP
