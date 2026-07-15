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

#ifndef CRYPTO3_MAC_GMAC_HPP
#define CRYPTO3_MAC_GMAC_HPP

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <nil/crypto3/mac/detail/ghash.hpp>
#include <nil/crypto3/mac/mac_key.hpp>

namespace nil {
    namespace crypto3 {
        namespace mac {
            namespace detail {
                constexpr inline bool valid_gmac_tag_bits(std::size_t tag_bits) {
                    return tag_bits == 128 || tag_bits == 120 || tag_bits == 112 || tag_bits == 104 || tag_bits == 96 ||
                           tag_bits == 64 || tag_bits == 32;
                }

                class ghash_state {
                public:
                    ghash_state() : ghash(ghash::zero_block()), associated_data_octets(0) {
                    }

                    explicit ghash_state(const ghash::block_type &hash_subkey) :
                        ghash(hash_subkey), associated_data_octets(0) {
                    }

                    template<typename InputIterator>
                    void update(InputIterator first, InputIterator last) {
                        while (first != last) {
                            increment_associated_data_size();
                            ghash.update_octet(static_cast<std::uint8_t>(*first++));
                        }
                    }

                    template<typename Cipher, std::size_t TagBits>
                    std::array<std::uint8_t, TagBits / CHAR_BIT> finalize(const Cipher &cipher,
                                                                          const ghash::block_type &j0) const {
                        ghash_state finalized = *this;
                        finalized.ghash.pad();
                        // GMAC authenticates data without encrypting a payload, so the GCM ciphertext length is zero.
                        finalized.ghash.update_lengths(finalized.associated_data_octets * CHAR_BIT, 0);

                        ghash::block_type tag = cipher.encrypt(j0);
                        ghash::xor_block(tag, finalized.ghash.digest());

                        std::array<std::uint8_t, TagBits / CHAR_BIT> truncated_tag = {0};
                        std::copy(tag.begin(), tag.begin() + truncated_tag.size(), truncated_tag.begin());
                        return truncated_tag;
                    }

                private:
                    void increment_associated_data_size() {
                        if (associated_data_octets == std::numeric_limits<std::uint64_t>::max() / CHAR_BIT) {
                            throw std::length_error("GMAC associated data is too long");
                        }

                        ++associated_data_octets;
                    }

                    ghash::state ghash;
                    std::uint64_t associated_data_octets;
                };
            }    // namespace detail

            /*!
             * @brief GMAC, the authentication-only specialization of GCM.
             *
             * GMAC uses a 128-bit block cipher, normally AES, in two roles. First it derives the GHASH subkey
             * H = AES_K(0^128). Then it derives the pre-counter block J0 from the IV: 96-bit IVs use
             * J0 = IV || 0x00000001, while other IV lengths are GHASHed with H as specified by GCM.
             *
             * The input passed to compute<gmac>(...) is authenticated data, not plaintext or ciphertext. GMAC hashes
             * that data as AAD and uses an empty ciphertext:
             *
             *     S = GHASH_H(AAD || pad || uint64_be(len(AAD)) || uint64_be(0))
             *     T = MSB_t(AES_K(J0) xor S)
             *
             * The caller is responsible for never reusing the same IV with the same key.
             *
             * @tparam BlockCipher 128-bit block cipher with byte-oriented block and key types.
             * @tparam TagBits output tag size in bits.
             * @ingroup mac
             */
            template<typename BlockCipher, std::size_t TagBits = 128>
            struct gmac {
                typedef BlockCipher cipher_type;

                constexpr static const std::size_t block_bits = cipher_type::block_bits;
                constexpr static const std::size_t block_octets = block_bits / CHAR_BIT;
                typedef typename cipher_type::block_type block_type;

                constexpr static const std::size_t key_bits = cipher_type::key_bits;
                constexpr static const std::size_t key_octets = key_bits / CHAR_BIT;
                typedef typename cipher_type::key_type key_type;

                constexpr static const std::size_t digest_bits = TagBits;
                constexpr static const std::size_t digest_octets = digest_bits / CHAR_BIT;
                typedef std::array<std::uint8_t, digest_octets> digest_type;

                static_assert(block_bits == 128, "GMAC requires a 128-bit block cipher");
                static_assert(TagBits % CHAR_BIT == 0, "GMAC tag size must be byte-aligned");
                static_assert(detail::valid_gmac_tag_bits(TagBits),
                              "GMAC tag size must be one of 128, 120, 112, 104, 96, 64, or 32 bits");
                static_assert(std::is_same<typename block_type::value_type, std::uint8_t>::value,
                              "GMAC currently requires byte-oriented block cipher blocks");
                static_assert(std::is_same<typename key_type::value_type, std::uint8_t>::value,
                              "GMAC currently requires byte-oriented block cipher keys");
            };

            template<typename BlockCipher, std::size_t TagBits>
            struct mac_key<gmac<BlockCipher, TagBits>> {
                typedef gmac<BlockCipher, TagBits> policy_type;
                typedef typename policy_type::cipher_type cipher_type;
                typedef typename policy_type::key_type key_type;
                typedef typename policy_type::digest_type digest_type;
                typedef detail::ghash_state accumulator_type;

                template<typename KeyRange, typename IvRange>
                mac_key(const KeyRange &key, const IvRange &iv) :
                    cipher(process_key(key.cbegin(), key.cend())),
                    hash_subkey(cipher.encrypt(detail::ghash::zero_block())),
                    j0(process_iv(iv.cbegin(), iv.cend(), hash_subkey)) {
                }

                inline void init_accumulator(accumulator_type &acc) const {
                    acc = accumulator_type(hash_subkey);
                }

                template<typename InputRange>
                inline void update(accumulator_type &acc, const InputRange &range) const {
                    update(acc, range.cbegin(), range.cend());
                }

                template<typename InputIterator>
                inline void update(accumulator_type &acc, InputIterator first, InputIterator last) const {
                    acc.update(first, last);
                }

                inline digest_type compute(const accumulator_type &acc) const {
                    return acc.template finalize<cipher_type, policy_type::digest_bits>(cipher, j0);
                }

            private:
                template<typename InputIterator>
                static key_type process_key(InputIterator first, InputIterator last) {
                    key_type key = {0};
                    typename key_type::iterator key_out = key.begin();

                    while (first != last) {
                        if (key_out == key.end()) {
                            throw std::invalid_argument("GMAC key size does not match the block cipher key size");
                        }

                        *key_out++ = static_cast<std::uint8_t>(*first++);
                    }

                    if (key_out != key.end()) {
                        throw std::invalid_argument("GMAC key size does not match the block cipher key size");
                    }

                    return key;
                }

                template<typename InputIterator>
                static detail::ghash::block_type
                    process_iv(InputIterator first, InputIterator last, const detail::ghash::block_type &hash_subkey) {
                    // Derive the GCM/GMAC pre-counter block J0 from the IV. The hash subkey H is AES_K(0^128).
                    // For the common 96-bit IV case, J0 is IV || 0x00000001. Other IV lengths are GHASHed using H.
                    std::vector<std::uint8_t> iv;
                    while (first != last) {
                        if (iv.size() == std::numeric_limits<std::uint64_t>::max() / CHAR_BIT) {
                            throw std::length_error("GMAC IV is too long");
                        }

                        iv.push_back(static_cast<std::uint8_t>(*first++));
                    }

                    if (iv.empty()) {
                        throw std::invalid_argument("GMAC IV must not be empty");
                    }

                    if (iv.size() == 12) {
                        detail::ghash::block_type block = detail::ghash::zero_block();
                        std::copy(iv.begin(), iv.end(), block.begin());
                        block[15] = 1;
                        return block;
                    }

                    detail::ghash::state ghash(hash_subkey);
                    ghash.update(iv.begin(), iv.end());
                    ghash.pad();
                    ghash.update_lengths(0, static_cast<std::uint64_t>(iv.size() * CHAR_BIT));
                    return ghash.digest();
                }

                cipher_type cipher;
                detail::ghash::block_type hash_subkey;
                // pre-counter block
                detail::ghash::block_type j0;
            };
        }    // namespace mac
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MAC_GMAC_HPP
