//---------------------------------------------------------------------------//
// Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
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

#ifndef CRYPTO3_RANDOM_CHACHA_HPP
#define CRYPTO3_RANDOM_CHACHA_HPP

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <tuple>
#include <vector>

#include <boost/config.hpp>
#include <boost/noncopyable.hpp>

#include <nil/crypto3/hash/sha2.hpp>

#include <nil/crypto3/mac/algorithm/compute.hpp>
#include <nil/crypto3/mac/hmac.hpp>

#include <nil/crypto3/stream/chacha.hpp>

namespace nil {
    namespace crypto3 {
        namespace random {

            /*!
             * @brief
             * @tparam MessageAuthenticationCode
             *
             * ChaCha_RNG is a very fast but completely ad-hoc RNG created by
             * creating a 256-bit random value and using it as a key for ChaCha20.
             *
             * The RNG maintains two 256-bit keys, one for HMAC_SHA256 (HK) and the
             * other for ChaCha20 (CK). To compute a new key in response to
             * reseeding request or add_entropy calls, ChaCha_RNG computes
             *   CK' = HMAC_SHA256(HK, input_material)
             * Then a new HK' is computed by running ChaCha20 with the new key to
             * output 32 bytes:
             *   HK' = ChaCha20(CK')
             *
             * Output is produced with ChaCha20 under CK'. Counter 0 is reserved for
             * refreshing HK, and caller-visible output starts at counter 1.
             *
             * The first HK (before seeding occurs) is taken as the all zero value.
             *
             * @warning This RNG construction is probably fine but is non-standard.
             * The primary reason to use it is in cases where the other RNGs are
             * not fast enough.
             */
            template<typename StreamCipher = stream::chacha20_cipher,
                     typename MessageAuthenticationCode = mac::hmac<hashes::sha2<256>>>
            struct chacha : private boost::noncopyable {
                typedef StreamCipher stream_cipher_type;
                typedef MessageAuthenticationCode mac_type;
                typedef mac::mac_key<mac_type> key_type;
                typedef typename stream_cipher_type::key_type cipher_key_type;
                typedef typename stream_cipher_type::nonce_type nonce_type;

                typedef std::vector<std::uint8_t> result_type;

                constexpr static const std::size_t result_size = stream_cipher_type::block_size;
                constexpr static const std::uint64_t max_output_bytes_per_seed =
                    static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max()) * result_size;
                constexpr static const bool has_fixed_range = false;

                static_assert(mac_type::digest_bits == std::tuple_size<cipher_key_type>::value * CHAR_BIT,
                              "ChaCha RNG requires MAC output to match the ChaCha20 key size");

                /** Returns the smallest value that the \random_device can produce. */
                static constexpr std::size_t min BOOST_PREVENT_MACRO_SUBSTITUTION() {
                    return 0;
                }

                /** Returns the largest value that the \random_device can produce. */
                static constexpr std::size_t max BOOST_PREVENT_MACRO_SUBSTITUTION() {
                    return std::numeric_limits<std::size_t>::max();
                }

                /** Constructs a @c random_device, optionally using the default device. */
                chacha() {
                    seed();
                }

                /**
                 * Constructs a @c random_device, optionally using the given token as an
                 * access specification (for example, a URL) to some implementation-defined
                 * service for monitoring a stochastic process.
                 */
                template<typename SeedSinglePassRange>
                explicit chacha(const SeedSinglePassRange &token) {
                    reset_state();
                    seed(token);
                }

                ~chacha() {
                    hk.fill(0);
                    ck.fill(0);
                }

                /** default seeds the underlying generator. */
                void seed() {
                    reset_state();
                    const std::array<std::uint8_t, 0> empty_seed = {};
                    update(empty_seed.begin(), empty_seed.end());
                }

                /** Seeds the underlying generator with first and last. */
                template<typename InputIterator>
                void seed(InputIterator first, InputIterator last) {
                    update(first, last);
                }

                template<typename SeedSinglePassRange>
                void seed(const SeedSinglePassRange &seed_material) {
                    seed(std::begin(seed_material), std::end(seed_material));
                }

                /**
                 * Returns: An entropy estimate for the random numbers returned by
                 * operator(), in the range min() to log2( max()+1). A deterministic
                 * random number generator (e.g. a pseudo-random number engine)
                 * has entropy 0.
                 *
                 * Throws: Nothing.
                 */
                double entropy() const {
                    return 0.0;
                }

                /** Returns a random value in the range [min, max]. */
                result_type operator()() {
                    result_type result(result_size);
                    generate(result.begin(), result.end());
                    return result;
                }

                /** Fills a range with random values. */
                template<class Iter>
                void generate(Iter begin, Iter end) {
                    if (begin == end) {
                        return;
                    }

                    validate_output_request(1);
                    stream_cipher_type cipher(ck, nonce, 1);
                    cipher.seek(output_offset);

                    std::array<std::uint8_t, result_size> zeros = {0};
                    std::array<std::uint8_t, result_size> block = {0};

                    while (begin != end) {
                        std::size_t chunk = 0;
                        Iter probe = begin;
                        while (probe != end && chunk != result_size) {
                            ++probe;
                            ++chunk;
                        }

                        validate_output_request(chunk);
                        cipher.process(zeros.begin(), zeros.begin() + chunk, block.begin());

                        for (std::size_t i = 0; i != chunk; ++i) {
                            *begin = block[i];
                            ++begin;
                        }

                        output_offset += chunk;
                    }
                }

            protected:
                void set_output_offset_for_testing(std::uint64_t offset) {
                    output_offset = offset;
                }

                void reset_state() {
                    hk.fill(0);
                    ck.fill(0);
                    nonce.fill(0);
                    output_offset = 0;
                }

                template<typename InputIterator>
                void update(InputIterator first, InputIterator last) {
                    key_type hmac_key(hk);
                    typename mac_type::digest_type digest =
                        static_cast<typename mac_type::digest_type>(compute<mac_type>(first, last, hmac_key));

                    std::copy(digest.begin(), digest.end(), ck.begin());
                    refresh_hmac_key();
                    output_offset = 0;
                }

                void refresh_hmac_key() {
                    std::array<std::uint8_t, std::tuple_size<cipher_key_type>::value> zeros = {0};
                    stream_cipher_type cipher(ck, nonce, 0);
                    cipher.process(zeros.begin(), zeros.end(), hk.begin());
                }

                void validate_output_request(std::size_t size) const {
                    if (output_offset > max_output_bytes_per_seed || size > max_output_bytes_per_seed - output_offset) {
                        throw std::out_of_range("ChaCha RNG output limit reached; reseed required");
                    }
                }

                cipher_key_type hk;
                cipher_key_type ck;
                nonce_type nonce;
                std::uint64_t output_offset;
            };
        }    // namespace random
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_RANDOM_CHACHA_HPP
