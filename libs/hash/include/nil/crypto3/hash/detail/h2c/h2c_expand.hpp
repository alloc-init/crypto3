//---------------------------------------------------------------------------//
// Copyright (c) 2020-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Ilias Khairullin <ilias@nil.foundation>
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

#ifndef CRYPTO3_HASH_H2C_EXPAND_HPP
#define CRYPTO3_HASH_H2C_EXPAND_HPP

#include <nil/crypto3/hash/detail/h2c/h2c_suites.hpp>
#include <nil/crypto3/hash/detail/h2c/h2c_sgn0.hpp>

#include <nil/crypto3/algebra/algorithms/strxor.hpp>

#include <nil/crypto3/algebra/curves/bls12.hpp>

#include <nil/crypto3/hash/algorithm/hash.hpp>
#include <nil/crypto3/hash/accumulators/hash.hpp>

#include <boost/assert.hpp>
#include <boost/static_assert.hpp>
#include <boost/concept/assert.hpp>

#include <array>
#include <type_traits>
#include <iterator>
#include <vector>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            using namespace nil::crypto3::detail;

            template<std::size_t k, typename HashType,
                    /// HashType::digest_type is required to be uint8_t[]
                    typename = typename std::enable_if<
                            std::is_same<std::uint8_t, typename HashType::digest_type::value_type>::value>::type>
            class expand_message_xmd {
                // https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-10#section-5.4.1
                BOOST_STATIC_ASSERT_MSG(HashType::block_bits % 8 == 0, "r_in_bytes is not a multiple of 8");
                BOOST_STATIC_ASSERT_MSG(HashType::digest_bits % 8 == 0, "b_in_bytes is not a multiple of 8");
                BOOST_STATIC_ASSERT_MSG(HashType::digest_bits >= 2 * k,
                                        "K-bit collision resistance is not fulfilled");

                constexpr static const std::size_t b_in_bytes = HashType::digest_bits / 8;
                constexpr static const std::size_t r_in_bytes = HashType::block_bits / 8;

                constexpr static const std::array<std::uint8_t, r_in_bytes> Z_pad{0};

                template<typename RangeType>
                static inline void append_range(std::vector<std::uint8_t> &out, const RangeType &range) {
                    for (std::uint8_t value : range) {
                        out.push_back(value);
                    }
                }

            public:
                template<typename InputMsgType, typename InputDstType, typename OutputType,
                        typename = typename std::enable_if<
                                std::is_same<std::uint8_t, typename InputMsgType::value_type>::value &&
                                std::is_same<std::uint8_t, typename InputDstType::value_type>::value &&
                                std::is_same<std::uint8_t, typename OutputType::value_type>::value>::type>
                static inline void process(const std::size_t len_in_bytes, const InputMsgType &msg,
                                           const InputDstType &dst, OutputType &uniform_bytes) {
                    BOOST_CONCEPT_ASSERT((boost::SinglePassRangeConcept<InputMsgType>));
                    BOOST_CONCEPT_ASSERT((boost::SinglePassRangeConcept<InputDstType>));
                    BOOST_CONCEPT_ASSERT((boost::SinglePassRangeConcept<OutputType>));
                    BOOST_CONCEPT_ASSERT((boost::WriteableRangeConcept<OutputType>));

                    // https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-10#section-5.4.1
                    BOOST_ASSERT(len_in_bytes < 0x10000);
                    const std::size_t dst_len = static_cast<std::size_t>(std::distance(dst.begin(), dst.end()));
                    BOOST_ASSERT(dst_len >= 16 && dst_len <= 255);
                    BOOST_ASSERT(
                            std::size_t(std::distance(uniform_bytes.begin(), uniform_bytes.end())) >= len_in_bytes);

                    const std::array<std::uint8_t, 2> l_i_b_str = {
                            static_cast<std::uint8_t>(len_in_bytes >> 8u),
                            static_cast<std::uint8_t>(len_in_bytes % 0x100)};
                    const std::size_t ell = static_cast<std::size_t>(len_in_bytes / b_in_bytes) +
                                            static_cast<std::size_t>(len_in_bytes % b_in_bytes != 0);

                    // https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-10#section-5.4.1
                    BOOST_ASSERT(ell <= 255);

                    // TODO: use accumulators when they will be fixed
                    // accumulator_set<HashType> b0_acc;
                    // hash<HashType>(Z_pad, b0_acc);
                    // hash<HashType>(msg, b0_acc);
                    // hash<HashType>(l_i_b_str, b0_acc);
                    // hash<HashType>(std::array<std::uint8_t, 1> {0}, b0_acc);
                    // hash<HashType>(dst, b0_acc);
                    // hash<HashType>(std::array<std::uint8_t, 1> {static_cast<std::uint8_t>(dst.size())},
                    // b0_acc); typename HashType::digest_type b0 =
                    // accumulators::extract::hash<HashType>(b0_acc);
                    std::vector<std::uint8_t> msg_prime;
                    msg_prime.reserve(Z_pad.size() + l_i_b_str.size() + 1 + dst_len + 1);
                    append_range(msg_prime, Z_pad);
                    append_range(msg_prime, msg);
                    append_range(msg_prime, l_i_b_str);
                    msg_prime.push_back(static_cast<std::uint8_t>(0));
                    append_range(msg_prime, dst);
                    msg_prime.push_back(static_cast<std::uint8_t>(dst_len));
                    typename HashType::digest_type b0 = hash<HashType>(msg_prime);

                    // TODO: use accumulators when they will be fixed
                    // accumulator_set<HashType> bi_acc;
                    // hash<HashType>(b0, bi_acc);
                    // hash<HashType>(std::array<std::uint8_t, 1> {1}, bi_acc);
                    // hash<HashType>(dst, bi_acc);
                    // hash<HashType>(std::array<std::uint8_t, 1> {static_cast<std::uint8_t>(dst.size())},
                    // bi_acc); typename HashType::digest_type bi =
                    // accumulators::extract::hash<HashType>(bi_acc); std::copy(bi.begin(), bi.end(),
                    // uniform_bytes.begin());
                    std::vector<std::uint8_t> b_i_str;
                    b_i_str.reserve(b0.size() + 1 + dst_len + 1);
                    append_range(b_i_str, b0);
                    b_i_str.push_back(static_cast<std::uint8_t>(1));
                    append_range(b_i_str, dst);
                    b_i_str.push_back(static_cast<std::uint8_t>(dst_len));
                    typename HashType::digest_type bi = hash<HashType>(b_i_str);
                    std::copy(bi.begin(), bi.end(), uniform_bytes.begin());

                    typename HashType::digest_type xored_b;
                    for (std::size_t i = 2; i <= ell; i++) {
                        // TODO: use accumulators when they will be fixed
                        // accumulator_set<HashType> bi_acc;
                        // strxor(b0, bi, xored_b);
                        // hash<HashType>(xored_b, bi_acc);
                        // hash<HashType>(std::array<std::uint8_t, 1> {static_cast<std::uint8_t>(i)}, bi_acc);
                        // hash<HashType>(dst, bi_acc);
                        // hash<HashType>(std::array<std::uint8_t, 1> {static_cast<std::uint8_t>(dst.size())},
                        //                bi_acc);
                        // bi = accumulators::extract::hash<HashType>(bi_acc);
                        // std::copy(bi.begin(), bi.end(), uniform_bytes.begin() + (i - 1) * b_in_bytes);
                        nil::crypto3::algebra::strxor(b0, bi, xored_b);
                        std::vector<std::uint8_t> b_i_str;
                        b_i_str.reserve(xored_b.size() + 1 + dst_len + 1);
                        append_range(b_i_str, xored_b);
                        b_i_str.push_back(static_cast<std::uint8_t>(i));
                        append_range(b_i_str, dst);
                        b_i_str.push_back(static_cast<std::uint8_t>(dst_len));
                        bi = hash<HashType>(b_i_str);
                        std::copy(bi.begin(), bi.end(), uniform_bytes.begin() + (i - 1) * b_in_bytes);
                    }
                }
            };
        } // namespace haseh
    } // namespace crypto3
} // namespace nil

#endif    // CRYPTO3_HASH_H@C_EXPAND_HPP
