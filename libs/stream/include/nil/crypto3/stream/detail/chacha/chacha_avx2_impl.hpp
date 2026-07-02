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

#ifndef CRYPTO3_STREAM_CHACHA_AVX2_IMPL_HPP
#define CRYPTO3_STREAM_CHACHA_AVX2_IMPL_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include <nil/crypto3/detail/config.hpp>

#include <nil/crypto3/stream/detail/chacha/chacha_policy.hpp>

#include <immintrin.h>

namespace nil {
    namespace crypto3 {
        namespace stream {
            namespace detail {
                template<std::size_t Round, std::size_t IVSize, std::size_t KeyBits>
                struct chacha_avx2_impl {
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
                                                   key_schedule_type &schedule) {
                        chacha_x8_impl(block.data(), schedule, counter_mode::original, 8);
                    }

                    static void chacha_x8_ietf(std::array<std::uint8_t, block_size * 8> &block,
                                               key_schedule_type &schedule) {
                        chacha_x8_impl(block.data(), schedule, counter_mode::ietf, 8);
                    }

                    static void chacha_x8(std::array<std::uint8_t, block_size * 8> &block,
                                          key_schedule_type &schedule) {
                        static_assert(IVSize == 64 || IVSize == 96, "ChaCha supports only 64-bit or 96-bit IVs");
                        if (IVSize == 96) {
                            chacha_x8_ietf(block, schedule);
                        } else {
                            chacha_x8_original(block, schedule);
                        }
                    }

                    static void chacha_x4_original(std::array<std::uint8_t, block_size * 4> &block,
                                                   key_schedule_type &schedule) {
                        std::array<std::uint8_t, block_size * 8> tmp = {0};
                        key_schedule_type working_schedule = schedule;

                        chacha_x8_impl(tmp.data(), working_schedule, counter_mode::original, 0);
                        std::copy(tmp.begin(), tmp.begin() + block.size(), block.begin());
                        advance_counter(schedule, counter_mode::original, 4);
                    }

                    static void chacha_x4_ietf(std::array<std::uint8_t, block_size * 4> &block,
                                               key_schedule_type &schedule) {
                        validate_can_advance_counter(schedule, counter_mode::ietf, 4);

                        std::array<std::uint8_t, block_size * 8> tmp = {0};
                        key_schedule_type working_schedule = schedule;

                        chacha_x8_impl(tmp.data(), working_schedule, counter_mode::ietf, 0);
                        std::copy(tmp.begin(), tmp.begin() + block.size(), block.begin());
                        advance_counter(schedule, counter_mode::ietf, 4);
                    }

                    static void chacha_x4(std::array<std::uint8_t, block_size * 4> &block,
                                          key_schedule_type &schedule) {
                        static_assert(IVSize == 64 || IVSize == 96, "ChaCha supports only 64-bit or 96-bit IVs");
                        if (IVSize == 96) {
                            chacha_x4_ietf(block, schedule);
                        } else {
                            chacha_x4_original(block, schedule);
                        }
                    }

                private:
                    enum class counter_mode { original, ietf };

                    static void validate_can_advance_counter(const key_schedule_type &schedule, counter_mode mode,
                                                             std::size_t blocks) {
                        const std::uint64_t max_counter = std::numeric_limits<std::uint32_t>::max();

                        if (mode == counter_mode::ietf && blocks != 0 &&
                            (blocks > max_counter || schedule[12] > max_counter - blocks)) {
                            throw std::out_of_range("ChaCha20 IETF counter exhausted");
                        }
                    }

                    static void advance_counter(key_schedule_type &schedule, counter_mode mode, std::size_t blocks) {
                        const word_type old_counter = schedule[12];
                        schedule[12] += static_cast<word_type>(blocks);
                        if (mode == counter_mode::original && schedule[12] < old_counter) {
                            ++schedule[13];
                        }
                    }

                    static BOOST_ATTRIBUTE_TARGET("avx2") void chacha_x8_impl(std::uint8_t *out,
                                                                              key_schedule_type &schedule,
                                                                              counter_mode mode,
                                                                              std::size_t advance_blocks) {
                        validate_can_advance_counter(schedule, mode, advance_blocks);

                        _mm256_zeroupper();

                        const __m256i CTR0 = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

                        const word_type C = 0xFFFFFFFF - schedule[12];
                        const __m256i CTR1 = _mm256_set_epi32(C < 7, C < 6, C < 5, C < 4, C < 3, C < 2, C < 1, 0);

                        const __m256i input0 = _mm256_set1_epi32(schedule[0]);
                        const __m256i input1 = _mm256_set1_epi32(schedule[1]);
                        const __m256i input2 = _mm256_set1_epi32(schedule[2]);
                        const __m256i input3 = _mm256_set1_epi32(schedule[3]);
                        const __m256i input4 = _mm256_set1_epi32(schedule[4]);
                        const __m256i input5 = _mm256_set1_epi32(schedule[5]);
                        const __m256i input6 = _mm256_set1_epi32(schedule[6]);
                        const __m256i input7 = _mm256_set1_epi32(schedule[7]);
                        const __m256i input8 = _mm256_set1_epi32(schedule[8]);
                        const __m256i input9 = _mm256_set1_epi32(schedule[9]);
                        const __m256i input10 = _mm256_set1_epi32(schedule[10]);
                        const __m256i input11 = _mm256_set1_epi32(schedule[11]);
                        const __m256i input12 = _mm256_add_epi32(_mm256_set1_epi32(schedule[12]), CTR0);
                        const __m256i input13 = _mm256_add_epi32(_mm256_set1_epi32(schedule[13]), CTR1);
                        const __m256i input14 = _mm256_set1_epi32(schedule[14]);
                        const __m256i input15 = _mm256_set1_epi32(schedule[15]);

                        __m256i R00 = input0;
                        __m256i R01 = input1;
                        __m256i R02 = input2;
                        __m256i R03 = input3;
                        __m256i R04 = input4;
                        __m256i R05 = input5;
                        __m256i R06 = input6;
                        __m256i R07 = input7;
                        __m256i R08 = input8;
                        __m256i R09 = input9;
                        __m256i R10 = input10;
                        __m256i R11 = input11;
                        __m256i R12 = input12;
                        __m256i R13 = input13;
                        __m256i R14 = input14;
                        __m256i R15 = input15;

#define CRYPTO3_CHACHA_AVX2_ADD(dst, src) dst = _mm256_add_epi32(dst, src)
#define CRYPTO3_CHACHA_AVX2_XOR(dst, src) dst = _mm256_xor_si256(dst, src)

                        for (std::size_t r = 0; r != rounds / 2; ++r) {
                            CRYPTO3_CHACHA_AVX2_ADD(R00, R04);
                            CRYPTO3_CHACHA_AVX2_ADD(R01, R05);
                            CRYPTO3_CHACHA_AVX2_ADD(R02, R06);
                            CRYPTO3_CHACHA_AVX2_ADD(R03, R07);

                            CRYPTO3_CHACHA_AVX2_XOR(R12, R00);
                            CRYPTO3_CHACHA_AVX2_XOR(R13, R01);
                            CRYPTO3_CHACHA_AVX2_XOR(R14, R02);
                            CRYPTO3_CHACHA_AVX2_XOR(R15, R03);

                            const __m256i shuf_rotl_16 =
                                _mm256_set_epi8(13, 12, 15, 14, 9, 8, 11, 10, 5, 4, 7, 6, 1, 0, 3, 2, 13, 12, 15, 14, 9,
                                                8, 11, 10, 5, 4, 7, 6, 1, 0, 3, 2);

                            R12 = _mm256_shuffle_epi8(R12, shuf_rotl_16);
                            R13 = _mm256_shuffle_epi8(R13, shuf_rotl_16);
                            R14 = _mm256_shuffle_epi8(R14, shuf_rotl_16);
                            R15 = _mm256_shuffle_epi8(R15, shuf_rotl_16);

                            CRYPTO3_CHACHA_AVX2_ADD(R08, R12);
                            CRYPTO3_CHACHA_AVX2_ADD(R09, R13);
                            CRYPTO3_CHACHA_AVX2_ADD(R10, R14);
                            CRYPTO3_CHACHA_AVX2_ADD(R11, R15);

                            CRYPTO3_CHACHA_AVX2_XOR(R04, R08);
                            CRYPTO3_CHACHA_AVX2_XOR(R05, R09);
                            CRYPTO3_CHACHA_AVX2_XOR(R06, R10);
                            CRYPTO3_CHACHA_AVX2_XOR(R07, R11);

                            R04 = _mm256_or_si256(_mm256_slli_epi32(R04, 12), _mm256_srli_epi32(R04, 32 - 12));
                            R05 = _mm256_or_si256(_mm256_slli_epi32(R05, 12), _mm256_srli_epi32(R05, 32 - 12));
                            R06 = _mm256_or_si256(_mm256_slli_epi32(R06, 12), _mm256_srli_epi32(R06, 32 - 12));
                            R07 = _mm256_or_si256(_mm256_slli_epi32(R07, 12), _mm256_srli_epi32(R07, 32 - 12));

                            CRYPTO3_CHACHA_AVX2_ADD(R00, R04);
                            CRYPTO3_CHACHA_AVX2_ADD(R01, R05);
                            CRYPTO3_CHACHA_AVX2_ADD(R02, R06);
                            CRYPTO3_CHACHA_AVX2_ADD(R03, R07);

                            CRYPTO3_CHACHA_AVX2_XOR(R12, R00);
                            CRYPTO3_CHACHA_AVX2_XOR(R13, R01);
                            CRYPTO3_CHACHA_AVX2_XOR(R14, R02);
                            CRYPTO3_CHACHA_AVX2_XOR(R15, R03);

                            const __m256i shuf_rotl_8 =
                                _mm256_set_epi8(14, 13, 12, 15, 10, 9, 8, 11, 6, 5, 4, 7, 2, 1, 0, 3, 14, 13, 12, 15,
                                                10, 9, 8, 11, 6, 5, 4, 7, 2, 1, 0, 3);

                            R12 = _mm256_shuffle_epi8(R12, shuf_rotl_8);
                            R13 = _mm256_shuffle_epi8(R13, shuf_rotl_8);
                            R14 = _mm256_shuffle_epi8(R14, shuf_rotl_8);
                            R15 = _mm256_shuffle_epi8(R15, shuf_rotl_8);

                            CRYPTO3_CHACHA_AVX2_ADD(R08, R12);
                            CRYPTO3_CHACHA_AVX2_ADD(R09, R13);
                            CRYPTO3_CHACHA_AVX2_ADD(R10, R14);
                            CRYPTO3_CHACHA_AVX2_ADD(R11, R15);

                            CRYPTO3_CHACHA_AVX2_XOR(R04, R08);
                            CRYPTO3_CHACHA_AVX2_XOR(R05, R09);
                            CRYPTO3_CHACHA_AVX2_XOR(R06, R10);
                            CRYPTO3_CHACHA_AVX2_XOR(R07, R11);

                            R04 = _mm256_or_si256(_mm256_slli_epi32(R04, 7), _mm256_srli_epi32(R04, 32 - 7));
                            R05 = _mm256_or_si256(_mm256_slli_epi32(R05, 7), _mm256_srli_epi32(R05, 32 - 7));
                            R06 = _mm256_or_si256(_mm256_slli_epi32(R06, 7), _mm256_srli_epi32(R06, 32 - 7));
                            R07 = _mm256_or_si256(_mm256_slli_epi32(R07, 7), _mm256_srli_epi32(R07, 32 - 7));

                            CRYPTO3_CHACHA_AVX2_ADD(R00, R05);
                            CRYPTO3_CHACHA_AVX2_ADD(R01, R06);
                            CRYPTO3_CHACHA_AVX2_ADD(R02, R07);
                            CRYPTO3_CHACHA_AVX2_ADD(R03, R04);

                            CRYPTO3_CHACHA_AVX2_XOR(R15, R00);
                            CRYPTO3_CHACHA_AVX2_XOR(R12, R01);
                            CRYPTO3_CHACHA_AVX2_XOR(R13, R02);
                            CRYPTO3_CHACHA_AVX2_XOR(R14, R03);

                            R15 = _mm256_shuffle_epi8(R15, shuf_rotl_16);
                            R12 = _mm256_shuffle_epi8(R12, shuf_rotl_16);
                            R13 = _mm256_shuffle_epi8(R13, shuf_rotl_16);
                            R14 = _mm256_shuffle_epi8(R14, shuf_rotl_16);

                            CRYPTO3_CHACHA_AVX2_ADD(R10, R15);
                            CRYPTO3_CHACHA_AVX2_ADD(R11, R12);
                            CRYPTO3_CHACHA_AVX2_ADD(R08, R13);
                            CRYPTO3_CHACHA_AVX2_ADD(R09, R14);

                            CRYPTO3_CHACHA_AVX2_XOR(R05, R10);
                            CRYPTO3_CHACHA_AVX2_XOR(R06, R11);
                            CRYPTO3_CHACHA_AVX2_XOR(R07, R08);
                            CRYPTO3_CHACHA_AVX2_XOR(R04, R09);

                            R05 = _mm256_or_si256(_mm256_slli_epi32(R05, 12), _mm256_srli_epi32(R05, 32 - 12));
                            R06 = _mm256_or_si256(_mm256_slli_epi32(R06, 12), _mm256_srli_epi32(R06, 32 - 12));
                            R07 = _mm256_or_si256(_mm256_slli_epi32(R07, 12), _mm256_srli_epi32(R07, 32 - 12));
                            R04 = _mm256_or_si256(_mm256_slli_epi32(R04, 12), _mm256_srli_epi32(R04, 32 - 12));

                            CRYPTO3_CHACHA_AVX2_ADD(R00, R05);
                            CRYPTO3_CHACHA_AVX2_ADD(R01, R06);
                            CRYPTO3_CHACHA_AVX2_ADD(R02, R07);
                            CRYPTO3_CHACHA_AVX2_ADD(R03, R04);

                            CRYPTO3_CHACHA_AVX2_XOR(R15, R00);
                            CRYPTO3_CHACHA_AVX2_XOR(R12, R01);
                            CRYPTO3_CHACHA_AVX2_XOR(R13, R02);
                            CRYPTO3_CHACHA_AVX2_XOR(R14, R03);

                            R15 = _mm256_shuffle_epi8(R15, shuf_rotl_8);
                            R12 = _mm256_shuffle_epi8(R12, shuf_rotl_8);
                            R13 = _mm256_shuffle_epi8(R13, shuf_rotl_8);
                            R14 = _mm256_shuffle_epi8(R14, shuf_rotl_8);

                            CRYPTO3_CHACHA_AVX2_ADD(R10, R15);
                            CRYPTO3_CHACHA_AVX2_ADD(R11, R12);
                            CRYPTO3_CHACHA_AVX2_ADD(R08, R13);
                            CRYPTO3_CHACHA_AVX2_ADD(R09, R14);

                            CRYPTO3_CHACHA_AVX2_XOR(R05, R10);
                            CRYPTO3_CHACHA_AVX2_XOR(R06, R11);
                            CRYPTO3_CHACHA_AVX2_XOR(R07, R08);
                            CRYPTO3_CHACHA_AVX2_XOR(R04, R09);

                            R05 = _mm256_or_si256(_mm256_slli_epi32(R05, 7), _mm256_srli_epi32(R05, 32 - 7));
                            R06 = _mm256_or_si256(_mm256_slli_epi32(R06, 7), _mm256_srli_epi32(R06, 32 - 7));
                            R07 = _mm256_or_si256(_mm256_slli_epi32(R07, 7), _mm256_srli_epi32(R07, 32 - 7));
                            R04 = _mm256_or_si256(_mm256_slli_epi32(R04, 7), _mm256_srli_epi32(R04, 32 - 7));
                        }

                        CRYPTO3_CHACHA_AVX2_ADD(R00, input0);
                        CRYPTO3_CHACHA_AVX2_ADD(R01, input1);
                        CRYPTO3_CHACHA_AVX2_ADD(R02, input2);
                        CRYPTO3_CHACHA_AVX2_ADD(R03, input3);
                        CRYPTO3_CHACHA_AVX2_ADD(R04, input4);
                        CRYPTO3_CHACHA_AVX2_ADD(R05, input5);
                        CRYPTO3_CHACHA_AVX2_ADD(R06, input6);
                        CRYPTO3_CHACHA_AVX2_ADD(R07, input7);
                        CRYPTO3_CHACHA_AVX2_ADD(R08, input8);
                        CRYPTO3_CHACHA_AVX2_ADD(R09, input9);
                        CRYPTO3_CHACHA_AVX2_ADD(R10, input10);
                        CRYPTO3_CHACHA_AVX2_ADD(R11, input11);
                        CRYPTO3_CHACHA_AVX2_ADD(R12, input12);
                        CRYPTO3_CHACHA_AVX2_ADD(R13, input13);
                        CRYPTO3_CHACHA_AVX2_ADD(R14, input14);
                        CRYPTO3_CHACHA_AVX2_ADD(R15, input15);

                        __m256i T0 = _mm256_unpacklo_epi32(R00, R01);
                        __m256i T1 = _mm256_unpacklo_epi32(R02, R03);
                        __m256i T2 = _mm256_unpackhi_epi32(R00, R01);
                        __m256i T3 = _mm256_unpackhi_epi32(R02, R03);

                        R00 = _mm256_unpacklo_epi64(T0, T1);
                        R01 = _mm256_unpackhi_epi64(T0, T1);
                        R02 = _mm256_unpacklo_epi64(T2, T3);
                        R03 = _mm256_unpackhi_epi64(T2, T3);

                        T0 = _mm256_unpacklo_epi32(R04, R05);
                        T1 = _mm256_unpacklo_epi32(R06, R07);
                        T2 = _mm256_unpackhi_epi32(R04, R05);
                        T3 = _mm256_unpackhi_epi32(R06, R07);

                        R04 = _mm256_unpacklo_epi64(T0, T1);
                        R05 = _mm256_unpackhi_epi64(T0, T1);
                        R06 = _mm256_unpacklo_epi64(T2, T3);
                        R07 = _mm256_unpackhi_epi64(T2, T3);

                        T0 = _mm256_unpacklo_epi32(R08, R09);
                        T1 = _mm256_unpacklo_epi32(R10, R11);
                        T2 = _mm256_unpackhi_epi32(R08, R09);
                        T3 = _mm256_unpackhi_epi32(R10, R11);

                        R08 = _mm256_unpacklo_epi64(T0, T1);
                        R09 = _mm256_unpackhi_epi64(T0, T1);
                        R10 = _mm256_unpacklo_epi64(T2, T3);
                        R11 = _mm256_unpackhi_epi64(T2, T3);

                        T0 = _mm256_unpacklo_epi32(R12, R13);
                        T1 = _mm256_unpacklo_epi32(R14, R15);
                        T2 = _mm256_unpackhi_epi32(R12, R13);
                        T3 = _mm256_unpackhi_epi32(R14, R15);

                        R12 = _mm256_unpacklo_epi64(T0, T1);
                        R13 = _mm256_unpackhi_epi64(T0, T1);
                        R14 = _mm256_unpacklo_epi64(T2, T3);
                        R15 = _mm256_unpackhi_epi64(T2, T3);

                        __m256i *output_mm = reinterpret_cast<__m256i *>(out);

                        _mm256_storeu_si256(output_mm, _mm256_permute2x128_si256(R00, R04, 0 + (2 << 4)));
                        _mm256_storeu_si256(output_mm + 1, _mm256_permute2x128_si256(R08, R12, 0 + (2 << 4)));
                        _mm256_storeu_si256(output_mm + 2, _mm256_permute2x128_si256(R01, R05, 0 + (2 << 4)));
                        _mm256_storeu_si256(output_mm + 3, _mm256_permute2x128_si256(R09, R13, 0 + (2 << 4)));
                        _mm256_storeu_si256(output_mm + 4, _mm256_permute2x128_si256(R02, R06, 0 + (2 << 4)));
                        _mm256_storeu_si256(output_mm + 5, _mm256_permute2x128_si256(R10, R14, 0 + (2 << 4)));
                        _mm256_storeu_si256(output_mm + 6, _mm256_permute2x128_si256(R03, R07, 0 + (2 << 4)));
                        _mm256_storeu_si256(output_mm + 7, _mm256_permute2x128_si256(R11, R15, 0 + (2 << 4)));

                        _mm256_storeu_si256(output_mm + 8, _mm256_permute2x128_si256(R00, R04, 1 + (3 << 4)));
                        _mm256_storeu_si256(output_mm + 9, _mm256_permute2x128_si256(R08, R12, 1 + (3 << 4)));
                        _mm256_storeu_si256(output_mm + 10, _mm256_permute2x128_si256(R01, R05, 1 + (3 << 4)));
                        _mm256_storeu_si256(output_mm + 11, _mm256_permute2x128_si256(R09, R13, 1 + (3 << 4)));
                        _mm256_storeu_si256(output_mm + 12, _mm256_permute2x128_si256(R02, R06, 1 + (3 << 4)));
                        _mm256_storeu_si256(output_mm + 13, _mm256_permute2x128_si256(R10, R14, 1 + (3 << 4)));
                        _mm256_storeu_si256(output_mm + 14, _mm256_permute2x128_si256(R03, R07, 1 + (3 << 4)));
                        _mm256_storeu_si256(output_mm + 15, _mm256_permute2x128_si256(R11, R15, 1 + (3 << 4)));

                        _mm256_zeroupper();

                        advance_counter(schedule, mode, advance_blocks);

#undef CRYPTO3_CHACHA_AVX2_ADD
#undef CRYPTO3_CHACHA_AVX2_XOR
                    }
                };
            }    // namespace detail
        }    // namespace stream
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_CHACHA_AVX2_IMPL_HPP
