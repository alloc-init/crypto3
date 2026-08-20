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

#ifndef CRYPTO3_MATH_DETAIL_INTEGER_SQRT_HPP
#define CRYPTO3_MATH_DETAIL_INTEGER_SQRT_HPP

#include <bit>
#include <concepts>

namespace nil::crypto3::math::detail {

    /** Return the smallest integer whose square is at least value. */
    template<std::unsigned_integral Integer>
    constexpr Integer ceil_sqrt(Integer value) {
        if (value < 2) {
            return value;
        }

        const Integer original_value = value;
        Integer root = 0;
        const unsigned shift = (std::bit_width(value) - 1) & ~1U;
        Integer bit = Integer {1} << shift;

        // Process one base-four digit per iteration. This constructs floor(sqrt(value)) using only shifts,
        // additions, and subtractions, and never forms an overflowing trial square.
        while (bit != 0) {
            if (value >= root + bit) {
                value = value - (root + bit);
                root = (root >> 1) + bit;
            } else {
                root >>= 1;
            }
            bit >>= 2;
        }

        // floor(sqrt(value))^2 cannot overflow because it is at most the original value.
        return root + static_cast<Integer>(root * root != original_value);
    }

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_DETAIL_INTEGER_SQRT_HPP
