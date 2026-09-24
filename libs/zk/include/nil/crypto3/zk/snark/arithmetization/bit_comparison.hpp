//---------------------------------------------------------------------------//
// Copyright (c) 2026 Alloc Init Labs Inc.
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

#ifndef CRYPTO3_ZK_BIT_COMPARISON_HPP
#define CRYPTO3_ZK_BIT_COMPARISON_HPP

#include <cstddef>
#include <stdexcept>

#include <boost/multiprecision/integer.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {

                struct bit_comparison_counts {
                    std::size_t marker_count = 0;
                    // Marker updates and zero-product checks; the caller determines their constraint cost.
                    std::size_t product_count = 0;
                };

                /**
                 * Visit the operations enforcing sum_j 2^j * b[j] < bound for Boolean little-endian bits.
                 * The positive integer bound is independent of the arithmetic field's characteristic.
                 * The caller constrains Booleanity and supplies the field's prime characteristic (>= 2).
                 *
                 * The caller maintains a prefix marker e, initially 1. After the preceding checks, e is
                 * nonzero exactly when the processed prefix equals the bound's prefix. In descending order:
                 * - reuse_marker(j) sets e = b[j], without allocating a marker or adding a product relation.
                 * - update_marker(j) allocates a fresh marker, zero exactly when e*b[j] is zero, and uses it
                 *   as the new e. For example, next = e*b[j] or next^2 = e*b[j], as the caller requires.
                 * - check_zero(begin, end) enforces e * sum_{begin <= j < end} b[j] = 0.
                 *
                 * Bit positions address the supplied bit sequence, not witness indices. Callers own variable
                 * allocation, constraint representation, witness generation and any additional certification.
                 * Zero runs are split into groups shorter than the characteristic to prevent modular wraparound.
                 * The final check is at the bound's lowest set bit; lower bits cannot affect the strict bound.
                 *
                 * Return the numbers of fresh markers and product relations. No callbacks run if the input
                 * width is smaller than the bound's bit length, including an empty bit sequence.
                 * Invalid bounds or characteristics throw std::invalid_argument before any callback runs.
                 */
                template<typename Integer, typename Characteristic, typename ReuseMarker, typename UpdateMarker,
                         typename CheckZero>
                constexpr bit_comparison_counts
                    for_each_bit_comparison(const Integer &bound, std::size_t bit_count,
                                            const Characteristic &characteristic, ReuseMarker &&reuse_marker,
                                            UpdateMarker &&update_marker, CheckZero &&check_zero) {
                    if (bound <= 0 || characteristic < 2) {
                        throw std::invalid_argument("bit comparison: require a positive bound and characteristic >= 2");
                    }
                    bit_comparison_counts counts;
                    const std::size_t bound_bits = boost::multiprecision::msb(bound) + 1;
                    if (bit_count < bound_bits) {
                        return counts;
                    }

                    const std::size_t last_one = boost::multiprecision::lsb(bound);
                    std::size_t position = bit_count;
                    while (position > last_one) {
                        if (position > bound_bits || !boost::multiprecision::bit_test(bound, position - 1)) {
                            const auto end = position;
                            do {
                                --position;
                            } while (position > last_one && end - position < characteristic - 1 &&
                                     (position > bound_bits || !boost::multiprecision::bit_test(bound, position - 1)));
                            check_zero(position, end);
                            ++counts.product_count;
                            continue;
                        }

                        --position;
                        if (position == last_one) {
                            check_zero(position, position + 1);
                            ++counts.product_count;
                        } else if (position == bound_bits - 1) {
                            reuse_marker(position);
                        } else {
                            update_marker(position);
                            ++counts.marker_count;
                            ++counts.product_count;
                        }
                    }
                    return counts;
                }

            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_BIT_COMPARISON_HPP
