//---------------------------------------------------------------------------//
// Copyright (c) 2026
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

#define BOOST_TEST_MODULE integer_sqrt_test

#include <cstdint>
#include <limits>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/math/detail/integer_sqrt.hpp>

namespace {
    namespace detail = nil::crypto3::math::detail;

    static_assert(detail::ceil_sqrt(std::uint32_t(0)) == 0);
    static_assert(detail::ceil_sqrt(std::uint32_t(1)) == 1);
    static_assert(detail::ceil_sqrt(std::uint32_t(2)) == 2);
    static_assert(detail::ceil_sqrt(std::uint32_t(16)) == 4);
    static_assert(detail::ceil_sqrt(std::uint32_t(17)) == 5);
}    // namespace

BOOST_AUTO_TEST_SUITE(integer_sqrt_test_suite)

BOOST_AUTO_TEST_CASE(ceil_sqrt_handles_exact_squares_and_their_neighbors) {
    BOOST_CHECK_EQUAL(detail::ceil_sqrt(std::uint64_t(3)), 2);
    BOOST_CHECK_EQUAL(detail::ceil_sqrt(std::uint64_t(4)), 2);
    BOOST_CHECK_EQUAL(detail::ceil_sqrt(std::uint64_t(5)), 3);
    BOOST_CHECK_EQUAL(detail::ceil_sqrt(std::uint64_t(15)), 4);
    BOOST_CHECK_EQUAL(detail::ceil_sqrt(std::uint64_t(16)), 4);
    BOOST_CHECK_EQUAL(detail::ceil_sqrt(std::uint64_t(17)), 5);
}

BOOST_AUTO_TEST_CASE(ceil_sqrt_handles_the_largest_unsigned_value_without_overflow) {
    const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
    BOOST_CHECK_EQUAL(detail::ceil_sqrt(maximum), std::uint64_t(1) << 32);
}

BOOST_AUTO_TEST_SUITE_END()
