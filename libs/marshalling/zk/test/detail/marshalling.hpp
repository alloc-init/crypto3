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

#ifndef CRYPTO3_MARSHALLING_ZK_TEST_DETAIL_MARSHALLING_HPP
#define CRYPTO3_MARSHALLING_ZK_TEST_DETAIL_MARSHALLING_HPP

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/status_type.hpp>
#include <nil/marshalling/types/integral.hpp>

#include <nil/crypto3/marshalling/multiprecision/types/integral.hpp>

namespace nil::crypto3::marshalling::test_tools {

    template<typename Marshalled>
    std::vector<std::uint8_t> encode(const Marshalled &filled) {
        std::vector<std::uint8_t> bytes(filled.length());
        auto output = bytes.begin();
        BOOST_REQUIRE(filled.write(output, bytes.size()) == nil::marshalling::status_type::success);
        BOOST_CHECK(output == bytes.end());
        return bytes;
    }

    template<typename Marshalled>
    Marshalled decode(const std::vector<std::uint8_t> &bytes) {
        Marshalled filled;
        auto input = bytes.begin();
        BOOST_REQUIRE(filled.read(input, bytes.size()) == nil::marshalling::status_type::success);
        // The test buffer holds one complete object, so no trailing bytes are expected.
        BOOST_CHECK(input == bytes.end());
        return filled;
    }

    // Write raw integers so malformed field encodings are not reduced before reaching the decoder.
    template<typename Endianness, typename Integral>
    void write_integer(std::vector<std::uint8_t> &bytes, std::size_t offset, const Integral &value) {
        using type_base = nil::marshalling::field_type<Endianness>;
        using integral_type =
            std::conditional_t<std::is_integral_v<Integral>, nil::marshalling::types::integral<type_base, Integral>,
                               types::integral<type_base, Integral>>;
        const integral_type raw(value);
        auto output = bytes.begin() + offset;
        BOOST_REQUIRE(raw.write(output, bytes.size() - offset) == nil::marshalling::status_type::success);
    }

}    // namespace nil::crypto3::marshalling::test_tools

#endif    // CRYPTO3_MARSHALLING_ZK_TEST_DETAIL_MARSHALLING_HPP
