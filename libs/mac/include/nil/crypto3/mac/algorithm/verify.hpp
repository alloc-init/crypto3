//---------------------------------------------------------------------------//

// Copyright (c) 2019-2020 Mikhail Komarov <nemo@nil.foundation>
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

#ifndef CRYPTO3_MAC_VERIFY_HPP
#define CRYPTO3_MAC_VERIFY_HPP

#include <cstdint>

#include <nil/crypto3/mac/algorithm/compute.hpp>
#include <nil/crypto3/mac/algorithm/mac.hpp>

#include <nil/crypto3/mac/mac_key.hpp>

namespace nil {
    namespace crypto3 {
        namespace mac {
            namespace detail {
                template<typename ActualRange, typename ExpectedIterator>
                bool constant_time_equal(ActualRange const &actual,
                                         ExpectedIterator expected_first,
                                         ExpectedIterator expected_last) {
                    // Tag length is public. Only the equal-length comparison below is constant-time with respect to
                    // tag contents.
                    std::uint8_t difference = 0;
                    for (std::uint8_t actual_octet : actual) {
                        if (expected_first == expected_last) {
                            return false;
                        }

                        difference |=
                            static_cast<std::uint8_t>(actual_octet ^ static_cast<std::uint8_t>(*expected_first++));
                    }

                    return expected_first == expected_last && difference == 0;
                }
            }    // namespace detail
        }    // namespace mac

        template<typename Mac, typename InputIterator, typename ExpectedIterator>
        bool verify(InputIterator first,
                    InputIterator last,
                    const mac::mac_key<Mac> &key,
                    ExpectedIterator expected_first,
                    ExpectedIterator expected_last) {
            const typename Mac::digest_type actual = compute<Mac>(first, last, key);
            return mac::detail::constant_time_equal(actual, expected_first, expected_last);
        }

        template<typename Mac, typename SinglePassRange, typename ExpectedRange>
        bool verify(const SinglePassRange &range, const mac::mac_key<Mac> &key, const ExpectedRange &expected) {
            return verify<Mac>(range.cbegin(), range.cend(), key, expected.cbegin(), expected.cend());
        }
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MAC_VERIFY_HPP
