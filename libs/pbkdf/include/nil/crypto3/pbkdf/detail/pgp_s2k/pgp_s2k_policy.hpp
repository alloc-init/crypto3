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

#ifndef CRYPTO3_PBKDF_PGP_S2K_POLICY_HPP
#define CRYPTO3_PBKDF_PGP_S2K_POLICY_HPP

#include <boost/integer.hpp>

namespace nil {
    namespace crypto3 {
        namespace pbkdf {
            namespace detail {
                template<typename HashType>
                struct pgp_s2k_policy {
                    typedef HashType hash_type;

                    constexpr static const std::size_t rounds = 256;
                    typedef typename boost::uint_t<CHAR_BIT * 4>::exact round_constant_type;
                    typedef typename std::array<round_constant_type, rounds> round_constants_type;
                    /*
                     * PGP stores the iteration count as a single byte
                     * Thus it can only actually take on one of 256 values, based on the
                     * formula in RFC 4880 section 3.6.1.3
                     */
                    constexpr static const round_constants_type constants = {
                        1024,     1088,     1152,     1216,     1280,     1344,     1408,     1472,     1536,
                        1600,     1664,     1728,     1792,     1856,     1920,     1984,     2048,     2176,
                        2304,     2432,     2560,     2688,     2816,     2944,     3072,     3200,     3328,
                        3456,     3584,     3712,     3840,     3968,     4096,     4352,     4608,     4864,
                        5120,     5376,     5632,     5888,     6144,     6400,     6656,     6912,     7168,
                        7424,     7680,     7936,     8192,     8704,     9216,     9728,     10240,    10752,
                        11264,    11776,    12288,    12800,    13312,    13824,    14336,    14848,    15360,
                        15872,    16384,    17408,    18432,    19456,    20480,    21504,    22528,    23552,
                        24576,    25600,    26624,    27648,    28672,    29696,    30720,    31744,    32768,
                        34816,    36864,    38912,    40960,    43008,    45056,    47104,    49152,    51200,
                        53248,    55296,    57344,    59392,    61440,    63488,    65536,    69632,    73728,
                        77824,    81920,    86016,    90112,    94208,    98304,    102400,   106496,   110592,
                        114688,   118784,   122880,   126976,   131072,   139264,   147456,   155648,   163840,
                        172032,   180224,   188416,   196608,   204800,   212992,   221184,   229376,   237568,
                        245760,   253952,   262144,   278528,   294912,   311296,   327680,   344064,   360448,
                        376832,   393216,   409600,   425984,   442368,   458752,   475136,   491520,   507904,
                        524288,   557056,   589824,   622592,   655360,   688128,   720896,   753664,   786432,
                        819200,   851968,   884736,   917504,   950272,   983040,   1015808,  1048576,  1114112,
                        1179648,  1245184,  1310720,  1376256,  1441792,  1507328,  1572864,  1638400,  1703936,
                        1769472,  1835008,  1900544,  1966080,  2031616,  2097152,  2228224,  2359296,  2490368,
                        2621440,  2752512,  2883584,  3014656,  3145728,  3276800,  3407872,  3538944,  3670016,
                        3801088,  3932160,  4063232,  4194304,  4456448,  4718592,  4980736,  5242880,  5505024,
                        5767168,  6029312,  6291456,  6553600,  6815744,  7077888,  7340032,  7602176,  7864320,
                        8126464,  8388608,  8912896,  9437184,  9961472,  10485760, 11010048, 11534336, 12058624,
                        12582912, 13107200, 13631488, 14155776, 14680064, 15204352, 15728640, 16252928, 16777216,
                        17825792, 18874368, 19922944, 20971520, 22020096, 23068672, 24117248, 25165824, 26214400,
                        27262976, 28311552, 29360128, 30408704, 31457280, 32505856, 33554432, 35651584, 37748736,
                        39845888, 41943040, 44040192, 46137344, 48234496, 50331648, 52428800, 54525952, 56623104,
                        58720256, 60817408, 62914560, 65011712};
                };
            }    // namespace detail
        }        // namespace pbkdf
    }            // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_PGP_S2K_POLICY_HPP
