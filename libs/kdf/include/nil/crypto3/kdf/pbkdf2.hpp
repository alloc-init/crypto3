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

#ifndef CRYPTO3_KDF_PBKDF2_HPP
#define CRYPTO3_KDF_PBKDF2_HPP

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <vector>

#include <nil/crypto3/mac/algorithm/compute.hpp>
#include <nil/crypto3/mac/mac_key.hpp>

namespace nil {
    namespace crypto3 {
        namespace kdf {
            namespace detail {
                template<typename InputIterator>
                std::vector<std::uint8_t> copy_octets(InputIterator first, InputIterator last) {
                    std::vector<std::uint8_t> out;
                    while (first != last) {
                        out.push_back(static_cast<std::uint8_t>(*first++));
                    }
                    return out;
                }

                inline void append_big_endian_32(std::vector<std::uint8_t> &out, std::uint32_t value) {
                    out.push_back(static_cast<std::uint8_t>(value >> 24));
                    out.push_back(static_cast<std::uint8_t>(value >> 16));
                    out.push_back(static_cast<std::uint8_t>(value >> 8));
                    out.push_back(static_cast<std::uint8_t>(value));
                }
            }    // namespace detail

            /*!
             * @brief PBKDF2 from RFC 8018 section 5.2.
             *
             * MessageAuthenticationCode is the PBKDF2 PRF. In normal use this is HMAC over a hash function, for
             * example mac::hmac<hashes::sha1> for the RFC 6070 test vectors or mac::hmac<hashes::sha2<256>> for a
             * modern SHA-256 based instance.
             *
             * @tparam MessageAuthenticationCode PRF used by PBKDF2.
             * @ingroup kdf
             */
            template<typename MessageAuthenticationCode>
            struct pbkdf2 {
                typedef MessageAuthenticationCode mac_type;
                typedef mac::mac_key<mac_type> mac_key_type;
                typedef typename mac_type::digest_type prf_result_type;

                constexpr static const std::size_t digest_bits = mac_type::digest_bits;
                constexpr static const std::size_t digest_octets = digest_bits / CHAR_BIT;

                static_assert(digest_bits % CHAR_BIT == 0, "PBKDF2 requires a byte-aligned PRF output");

                template<typename PasswordIterator, typename SaltIterator, typename OutputIterator>
                static OutputIterator derive(PasswordIterator password_first,
                                             PasswordIterator password_last,
                                             SaltIterator salt_first,
                                             SaltIterator salt_last,
                                             std::size_t iterations,
                                             std::size_t derived_key_octets,
                                             OutputIterator out) {
                    if (iterations == 0) {
                        throw std::invalid_argument("PBKDF2 iteration count must be greater than zero");
                    }
                    if (derived_key_octets == 0) {
                        throw std::invalid_argument("PBKDF2 derived key length must be greater than zero");
                    }

                    const std::size_t full_blocks = derived_key_octets / digest_octets;
                    const bool has_partial_block = (derived_key_octets % digest_octets) != 0;
                    if (full_blocks > std::numeric_limits<std::uint32_t>::max() ||
                        (full_blocks == std::numeric_limits<std::uint32_t>::max() && has_partial_block)) {
                        throw std::length_error("PBKDF2 derived key length is too large");
                    }

                    const std::vector<std::uint8_t> password = detail::copy_octets(password_first, password_last);
                    const std::vector<std::uint8_t> salt = detail::copy_octets(salt_first, salt_last);
                    const mac_key_type key(password);

                    const std::size_t block_count = full_blocks + (has_partial_block ? 1 : 0);
                    std::vector<std::uint8_t> salt_with_counter;
                    salt_with_counter.reserve(salt.size() + 4);

                    // RFC 8018 defines each derived-key block as:
                    //
                    //     T_i = U_1 xor U_2 xor ... xor U_c
                    //     U_1 = PRF(P, S || INT(i))
                    //     U_j = PRF(P, U_{j-1})
                    //
                    // where P is the password, S is the salt, c is the iteration count, and INT(i) is the
                    // one-based block index encoded as a 32-bit big-endian integer. The final T_i block is
                    // truncated when the requested derived key length is not a multiple of the PRF output size.
                    for (std::size_t block_index = 1; block_index <= block_count; ++block_index) {
                        salt_with_counter.assign(salt.begin(), salt.end());
                        detail::append_big_endian_32(salt_with_counter, static_cast<std::uint32_t>(block_index));

                        prf_result_type u = compute<mac_type>(salt_with_counter, key);
                        prf_result_type block = u;

                        for (std::size_t i = 1; i < iterations; ++i) {
                            u = compute<mac_type>(u, key);
                            for (std::size_t j = 0; j < block.size(); ++j) {
                                block[j] ^= u[j];
                            }
                        }

                        const std::size_t output_offset = (block_index - 1) * digest_octets;
                        const std::size_t output_remaining = derived_key_octets - output_offset;
                        const std::size_t output_this_block = std::min(block.size(), output_remaining);
                        out = std::copy_n(block.begin(), output_this_block, out);
                    }

                    return out;
                }

                template<typename PasswordRange, typename SaltRange, typename OutputIterator>
                static OutputIterator derive(const PasswordRange &password,
                                             const SaltRange &salt,
                                             std::size_t iterations,
                                             std::size_t derived_key_octets,
                                             OutputIterator out) {
                    return derive(password.cbegin(), password.cend(), salt.cbegin(), salt.cend(), iterations,
                                  derived_key_octets, out);
                }

                template<typename PasswordRange, typename SaltRange>
                static std::vector<std::uint8_t> derive(const PasswordRange &password,
                                                        const SaltRange &salt,
                                                        std::size_t iterations,
                                                        std::size_t derived_key_octets) {
                    std::vector<std::uint8_t> derived_key;
                    derived_key.reserve(derived_key_octets);
                    derive(password, salt, iterations, derived_key_octets, std::back_inserter(derived_key));
                    return derived_key;
                }
            };
        }    // namespace kdf
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_KDF_PBKDF2_HPP
