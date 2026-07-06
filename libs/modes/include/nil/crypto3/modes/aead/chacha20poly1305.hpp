//---------------------------------------------------------------------------//
// Copyright (c) 2019 Mikhail Komarov <nemo@nil.foundation>
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

#ifndef CRYPTO3_MODE_AEAD_CHACHA20_POLY1305_HPP
#define CRYPTO3_MODE_AEAD_CHACHA20_POLY1305_HPP

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <vector>

#include <nil/crypto3/mac/algorithm/compute.hpp>
#include <nil/crypto3/mac/poly1305.hpp>
#include <nil/crypto3/stream/chacha.hpp>

namespace nil {
    namespace crypto3 {
        namespace modes {
            namespace aead {
                class chacha20poly1305 {
                public:
                    typedef stream::chacha20 stream_type;
                    typedef stream_type::key_type key_type;
                    typedef stream_type::iv_type nonce_type;
                    typedef mac::poly1305::key_type poly1305_key_type;
                    typedef mac::poly1305::digest_type tag_type;

                    constexpr static const std::size_t key_size = stream_type::key_bits / CHAR_BIT;
                    constexpr static const std::size_t nonce_size = stream_type::iv_bits / CHAR_BIT;
                    constexpr static const std::size_t tag_size = mac::poly1305::digest_octets;

                    static poly1305_key_type poly1305_key_gen(const key_type &key, const nonce_type &nonce) {
                        std::array<std::uint8_t, stream_type::block_size> zeros = {0};
                        std::array<std::uint8_t, stream_type::block_size> block = {0};

                        stream::chacha20_cipher cipher(key, nonce, 0);
                        cipher.process(zeros.begin(), zeros.end(), block.begin());

                        poly1305_key_type one_time_key = {0};
                        std::copy(block.begin(), block.begin() + one_time_key.size(), one_time_key.begin());
                        return one_time_key;
                    }

                    template<typename PlaintextIterator, typename AadIterator, typename OutputIterator>
                    static tag_type encrypt(PlaintextIterator plaintext_first, PlaintextIterator plaintext_last,
                                            AadIterator aad_first, AadIterator aad_last, const key_type &key,
                                            const nonce_type &nonce, OutputIterator ciphertext_out) {
                        const std::vector<std::uint8_t> plaintext = to_byte_vector(plaintext_first, plaintext_last);
                        const std::vector<std::uint8_t> aad = to_byte_vector(aad_first, aad_last);
                        ensure_plaintext_size(plaintext.size());

                        std::vector<std::uint8_t> ciphertext(plaintext.size());
                        stream::chacha20_cipher cipher(key, nonce, 1);
                        cipher.process(plaintext.begin(), plaintext.end(), ciphertext.begin());

                        std::copy(ciphertext.begin(), ciphertext.end(), ciphertext_out);
                        return compute_tag(aad, ciphertext, key, nonce);
                    }

                    template<typename PlaintextRange, typename AadRange, typename OutputIterator>
                    static tag_type encrypt(const PlaintextRange &plaintext, const AadRange &aad, const key_type &key,
                                            const nonce_type &nonce, OutputIterator ciphertext_out) {
                        return encrypt(std::begin(plaintext), std::end(plaintext), std::begin(aad), std::end(aad), key,
                                       nonce, ciphertext_out);
                    }

                    template<typename CiphertextIterator, typename AadIterator, typename OutputIterator>
                    static bool decrypt(CiphertextIterator ciphertext_first, CiphertextIterator ciphertext_last,
                                        AadIterator aad_first, AadIterator aad_last, const tag_type &tag,
                                        const key_type &key, const nonce_type &nonce, OutputIterator plaintext_out) {
                        const std::vector<std::uint8_t> ciphertext = to_byte_vector(ciphertext_first, ciphertext_last);
                        const std::vector<std::uint8_t> aad = to_byte_vector(aad_first, aad_last);
                        ensure_plaintext_size(ciphertext.size());

                        const tag_type expected_tag = compute_tag(aad, ciphertext, key, nonce);
                        if (!constant_time_equal(expected_tag, tag)) {
                            return false;
                        }

                        std::vector<std::uint8_t> plaintext(ciphertext.size());
                        stream::chacha20_cipher cipher(key, nonce, 1);
                        cipher.process(ciphertext.begin(), ciphertext.end(), plaintext.begin());
                        std::copy(plaintext.begin(), plaintext.end(), plaintext_out);
                        return true;
                    }

                    template<typename CiphertextRange, typename AadRange, typename OutputIterator>
                    static bool decrypt(const CiphertextRange &ciphertext, const AadRange &aad, const tag_type &tag,
                                        const key_type &key, const nonce_type &nonce, OutputIterator plaintext_out) {
                        return decrypt(std::begin(ciphertext), std::end(ciphertext), std::begin(aad), std::end(aad),
                                       tag, key, nonce, plaintext_out);
                    }

                    template<typename ByteRange>
                    static bool constant_time_equal(const ByteRange &lhs, const ByteRange &rhs) {
                        if (lhs.size() != rhs.size()) {
                            return false;
                        }

                        std::uint8_t diff = 0;
                        for (std::size_t i = 0; i != lhs.size(); ++i) {
                            diff |= static_cast<std::uint8_t>(lhs[i] ^ rhs[i]);
                        }
                        return diff == 0;
                    }

                private:
                    static constexpr std::uint64_t max_plaintext_size() {
                        return static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max()) *
                               stream_type::block_size;
                    }

                    static void ensure_plaintext_size(std::size_t size) {
                        if (size > max_plaintext_size()) {
                            throw std::out_of_range("ChaCha20-Poly1305 plaintext is too large");
                        }
                    }

                    template<typename InputIterator>
                    static std::vector<std::uint8_t> to_byte_vector(InputIterator first, InputIterator last) {
                        std::vector<std::uint8_t> bytes;
                        for (; first != last; ++first) {
                            bytes.push_back(static_cast<std::uint8_t>(*first));
                        }
                        return bytes;
                    }

                    static void append_pad16(std::vector<std::uint8_t> &out, std::size_t size) {
                        const std::size_t padding = size % 16 == 0 ? 0 : 16 - size % 16;
                        out.insert(out.end(), padding, 0);
                    }

                    static void append_u64_le(std::vector<std::uint8_t> &out, std::uint64_t value) {
                        for (std::size_t i = 0; i != 8; ++i) {
                            out.push_back(static_cast<std::uint8_t>((value >> (8 * i)) & 0xff));
                        }
                    }

                    static std::vector<std::uint8_t> mac_data(const std::vector<std::uint8_t> &aad,
                                                              const std::vector<std::uint8_t> &ciphertext) {
                        std::vector<std::uint8_t> data;
                        data.reserve(aad.size() + ciphertext.size() + 32);

                        data.insert(data.end(), aad.begin(), aad.end());
                        append_pad16(data, aad.size());
                        data.insert(data.end(), ciphertext.begin(), ciphertext.end());
                        append_pad16(data, ciphertext.size());
                        append_u64_le(data, static_cast<std::uint64_t>(aad.size()));
                        append_u64_le(data, static_cast<std::uint64_t>(ciphertext.size()));
                        return data;
                    }

                    static tag_type compute_tag(const std::vector<std::uint8_t> &aad,
                                                const std::vector<std::uint8_t> &ciphertext, const key_type &key,
                                                const nonce_type &nonce) {
                        const poly1305_key_type one_time_key = poly1305_key_gen(key, nonce);
                        const mac::mac_key<mac::poly1305> poly1305_key(one_time_key);
                        const std::vector<std::uint8_t> data = mac_data(aad, ciphertext);
                        return compute<mac::poly1305>(data, poly1305_key);
                    }
                };
            }    // namespace aead
        }    // namespace modes

        namespace stream {
            namespace modes {
                typedef ::nil::crypto3::modes::aead::chacha20poly1305 chacha20poly1305;
            }    // namespace modes
        }    // namespace stream
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MODE_AEAD_CHACHA20_POLY1305_HPP
