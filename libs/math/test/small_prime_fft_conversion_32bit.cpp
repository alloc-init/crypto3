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

#include <boost/config.hpp>
// Exercise the fallback with 32-bit Boost limbs even on a 64-bit host. Configure
// Boost first so subsequent includes cannot re-enable its 128-bit integer support.
#undef BOOST_HAS_INT128

#include <nil/crypto3/math/detail/bn254_small_prime_fft.hpp>

#include <array>

using converter_type = nil::crypto3::math::detail::bn254_small_prime_fft_converter;
using base_field_type = nil::crypto3::algebra::fields::alt_bn128_base_field<254>;
static_assert(base_field_type::integral_type::backend_type::limb_bits == 32);

void instantiate_32bit_conversion() {
    // This call must discard the 64-bit Montgomery helper rather than instantiate it.
    std::array<converter_type::value_type, 1> input = {converter_type::value_type::zero()};
    std::array<double, 12> output {};
    converter_type(1).convert(output, input, 1);
}
