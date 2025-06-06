//---------------------------------------------------------------------------//
// Copyright (c) 2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2021 Ilias Khairullin <ilias@nil.foundation>
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

#ifndef CRYPTO3_HASH_H2C_HPP
#define CRYPTO3_HASH_H2C_HPP

#include <string>
#include <vector>

#include <nil/crypto3/hash/accumulators/hash.hpp>
#include <nil/crypto3/hash/detail/h2f/h2f_suites.hpp>
#include <nil/crypto3/hash/detail/h2c/h2c_suites.hpp>
#include <nil/crypto3/hash/detail/h2c/h2c_functions.hpp>
#include <nil/crypto3/hash/detail/stream_processors/stream_processors_enum.hpp>
#include <nil/crypto3/hash/h2f.hpp>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            template<typename GroupType, typename HashType, std::size_t _k = 128,
                     uniformity_count_t _uniformity_count = uniformity_count_t::uniform_count,
                     expand_msg_variant_t _expand_msg_variant = expand_msg_variant_t::rfc_xmd>
            struct h2c_default_params {
                constexpr static uniformity_count_t uniformity_count = _uniformity_count;
                constexpr static expand_msg_variant_t expand_msg_variant = _expand_msg_variant;
                constexpr static std::size_t k = _k;

                typedef std::vector<std::uint8_t> dst_type;
                static inline dst_type dst = []() {
                    using internal_suite_type = h2f_suite<typename GroupType::field_type, HashType, k>;
                    std::string default_tag_str = "QUUX-V01-CS02-with-";
                    dst_type dst(default_tag_str.begin(), default_tag_str.end());
                    dst.insert(dst.end(), internal_suite_type::suite_id.begin(), internal_suite_type::suite_id.end());
                    return dst;
                }();
            };

            /*!
             * @brief Hashing to Elliptic Curves
             * https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11
             *
             * @tparam GroupType
             * @tparam Params
             */
            template<typename GroupType, typename HashType = sha2<256>,
                     typename ParamsType = h2c_default_params<GroupType, HashType>>
            struct h2c {
                typedef h2c_suite<GroupType> suite_type;

                typedef typename suite_type::group_type group_type;
                typedef typename suite_type::group_value_type group_value_type;
                typedef typename suite_type::field_type field_type;
                typedef typename suite_type::field_value_type field_value_type;
                typedef typename suite_type::modular_type modular_type;
                typedef typename suite_type::integral_type integral_type;

                typedef h2f<field_type, HashType, ParamsType> hash_type;

                typedef typename hash_type::accumulator_type accumulator_type;
                typedef group_value_type result_type;
                typedef result_type digest_type;

                constexpr static std::size_t digest_bits = group_type::value_bits;

                struct construction {
                    struct params_type {
                        typedef nil::marshalling::option::big_endian digest_endian;
                    };

                    typedef void type;
                };

                constexpr static detail::stream_processor_type stream_processor =
                    detail::stream_processor_type::raw_delegating;
                using accumulator_tag = accumulators::tag::forwarding_hash<h2c<GroupType, HashType, ParamsType>>;

                static inline void init_accumulator(accumulator_type &acc) {
                    hash_type::init_accumulator(acc);
                }

                template<typename InputRange>
                static inline void update(accumulator_type &acc, const InputRange &range) {
                    hash_type::update(acc, range);
                }

                template<typename InputIterator>
                static inline void update(accumulator_type &acc, InputIterator first, InputIterator last) {
                    hash_type::update(acc, first, last);
                }

                static inline result_type process(accumulator_type &acc) {
                    return detail::ep_map<group_type, ParamsType::uniformity_count>(hash_type::process(acc));
                }
            };
        }    // namespace hashes
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_HASH_H2C_HPP
