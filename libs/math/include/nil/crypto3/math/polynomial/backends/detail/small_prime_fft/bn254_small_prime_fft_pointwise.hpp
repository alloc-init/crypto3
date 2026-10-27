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

#ifndef CRYPTO3_MATH_DETAIL_BN254_SMALL_PRIME_FFT_POINTWISE_HPP
#define CRYPTO3_MATH_DETAIL_BN254_SMALL_PRIME_FFT_POINTWISE_HPP

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>

#include <boost/config.hpp>

#include <nil/crypto3/math/polynomial/backends/detail/small_prime_fft/small_prime_fft_arithmetic.hpp>

namespace nil::crypto3::math::detail {

    /**
     * Pointwise products in the BN254 tower u^2=-1, v^3=9+u, w^2=v, over one small modulus q.
     * This is a quotient ring; irreducibility of the tower polynomials modulo q is not required.
     * Coordinate channels have the same order as the converter, with one prime per call.
     *
     * Reuse the field library's Karatsuba tower identities and quadratic-extension squaring
     * identity with scalar, four-lane or eight-lane residue kernels. The existing field elements
     * bind their base arithmetic to the target modulus and cannot represent these small-prime intermediates.
     *
     * @pre Input cells are integer-valued doubles of magnitude < 2^53. scale is an integer with
     * |scale| < q/2, normally the FFT plan's centered 1/N. Buffers are disjoint or identical.
     * The thread's rounding mode stays unchanged during the call. Output cells are in (-q,q).
     * No buffers are allocated and no constant-time guarantee is provided.
     */
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
// The packed kernels use vector coordinates in the existing tower arrays.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
    class bn254_small_prime_fft_pointwise {

        template<bool UseFma, typename Value, typename Arithmetic>
        struct kernel {
            using quadratic_type = std::array<Value, 2>;
            using cubic_type = std::array<quadratic_type, 3>;
            using tower_type = std::array<cubic_type, 2>;

            const Arithmetic &arithmetic;

            Value reduce(Value value) const {
                if constexpr (UseFma) {
                    return arithmetic.reduce_fma(value);
                } else {
                    return arithmetic.reduce_integer(value);
                }
            }

            Value multiply(Value a, Value b) const {
                if constexpr (UseFma) {
                    return arithmetic.multiply_fma(a, b);
                } else {
                    return arithmetic.multiply_integer(a, b);
                }
            }

            BOOST_FORCEINLINE void multiply_pair(Value &a, Value &b, Value wa, Value wb) const {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                if constexpr (UseFma && requires { arithmetic.multiply_fma_pair(a, b, wa, wb); }) {
                    arithmetic.multiply_fma_pair(a, b, wa, wb);
                } else
#endif
                {
                    a = multiply(a, wa);
                    b = multiply(b, wb);
                }
            }

            template<bool Reduce = true>
            Value add(Value a, Value b) const {
                if constexpr (Reduce) {
                    return reduce(a + b);
                } else {
                    return a + b;
                }
            }

            // Three reduced coordinates give an exact integer of magnitude < 3*q.
            // Callers can defer this reduction to the next bounded sum or final scaling.
            template<bool Reduce = true>
            Value subtract(Value a, Value b, Value c) const {
                if constexpr (Reduce) {
                    return reduce(a - b - c);
                } else {
                    return a - b - c;
                }
            }

            template<bool Reduce = true, typename T, std::size_t N>
            std::array<T, N> add(const std::array<T, N> &a, const std::array<T, N> &b) const {
                std::array<T, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                    result[i] = add<Reduce>(a[i], b[i]);
                }
                return result;
            }

            template<bool Reduce = true, typename T, std::size_t N>
            std::array<T, N> subtract(const std::array<T, N> &a, const std::array<T, N> &b,
                                      const std::array<T, N> &c) const {
                std::array<T, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                    result[i] = subtract<Reduce>(a[i], b[i], c[i]);
                }
                return result;
            }

            // Quadratic inputs and returned coordinates are in (-q,q). Pair sums are < 2*q;
            // the Karatsuba product is < 4*q^2 and its remainder < 1.5*q. Correcting by the
            // two diagonal products gives a sum < 3.5*q, hence exact in binary64 before reduction.
            BOOST_FORCEINLINE quadratic_type multiply(const quadratic_type &a, const quadratic_type &b) const {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
                // Keep the two independent products in registers and interleave their FMA chains.
                Value diagonal0 = a[0];
                Value diagonal1 = a[1];
                multiply_pair(diagonal0, diagonal1, b[0], b[1]);
                const Value cross = multiply(a[0] + a[1], b[0] + b[1]);
                return {reduce(diagonal0 - diagonal1), reduce(cross - diagonal0 - diagonal1)};
            }

            quadratic_type multiply_xi(const quadratic_type &a) const {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
                // (9+u)*(a0+a1*u) = (9*a0-a1) + (a0+9*a1)*u. Reduce 8*a first:
                // |8*a| < 8*q < 2^53, whereas forming 9*a directly could lose an integer bit.
                return {reduce(reduce(8 * a[0]) + a[0] - a[1]), reduce(reduce(8 * a[1]) + a[1] + a[0])};
            }

            quadratic_type add_multiply_xi(const quadratic_type &a, const quadratic_type &b) const {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
                // a + (9+u)*b, with |a| < 3*q and |b| < q. After reducing 8*b,
                // each four-term sum is < 6*q < 2^53, so only its final reduction is needed.
                return {reduce(reduce(8 * b[0]) + b[0] - b[1] + a[0]), reduce(reduce(8 * b[1]) + b[1] + b[0] + a[1])};
            }

            cubic_type multiply_v(const cubic_type &a) const {
                return {multiply_xi(a[2]), a[0], a[1]};
            }

            cubic_type multiply(const cubic_type &a, const cubic_type &b) const {
                const auto diagonal0 = multiply(a[0], b[0]);
                const auto diagonal1 = multiply(a[1], b[1]);
                const auto diagonal2 = multiply(a[2], b[2]);
                // cross12 is reduced before multiplication by 8 inside multiply_xi.
                // The other cross terms stay < 3*q until these final sums; the cubic result is < q.
                // Finish each output pair before forming the next cross term to limit live temporaries.
                cubic_type result;
                {
                    const auto cross12 = subtract(multiply(add(a[1], a[2]), add(b[1], b[2])), diagonal1, diagonal2);
                    result[0] = add_multiply_xi(diagonal0, cross12);
                }
                {
                    const auto cross01 =
                        subtract<false>(multiply(add(a[0], a[1]), add(b[0], b[1])), diagonal0, diagonal1);
                    result[1] = add_multiply_xi(cross01, diagonal2);
                }
                {
                    const auto cross02 =
                        subtract<false>(multiply(add(a[0], a[2]), add(b[0], b[2])), diagonal0, diagonal2);
                    result[2] = add(cross02, diagonal1);
                }
                return result;
            }

            // Both tower operations return coordinates of magnitude < 3*q. The final
            // |scale| < q/2 gives a product < 1.5*q^2, whose FMA remainder is already in (-q,q).
            // Keeping these sums unreduced is also safe when two tower products are added below.
            tower_type multiply(const tower_type &a, const tower_type &b) const {
                const auto diagonal0 = multiply(a[0], b[0]);
                const auto diagonal1 = multiply(a[1], b[1]);
                const auto cross = subtract<false>(multiply(add(a[0], a[1]), add(b[0], b[1])), diagonal0, diagonal1);
                return {add<false>(diagonal0, multiply_v(diagonal1)), cross};
            }

            tower_type square(const tower_type &a) const {
                // The quadratic-extension square needs two cubic products instead of three:
                // c0 = (a0+a1)*(a0+v*a1) - a0*a1 - v*a0*a1; c1 = 2*a0*a1.
                const auto cross = multiply(a[0], a[1]);
                const auto combined = multiply(add(a[0], a[1]), add(a[0], multiply_v(a[1])));
                return {subtract<false>(combined, cross, multiply_v(cross)), add<false>(cross, cross)};
            }

            template<bool Square>
            BOOST_FORCEINLINE void load(tower_type &a, tower_type &b, std::span<const double> left,
                                        std::span<const double> right, std::size_t stride, std::size_t index) const {
                constexpr std::size_t width = sizeof(Value) / sizeof(double);
                for (std::size_t w = 0; w < 2; ++w) {
                    for (std::size_t v = 0; v < 3; ++v) {
                        for (std::size_t u = 0; u < 2; ++u) {
                            const std::size_t offset = (6 * w + 2 * v + u) * stride + index;
                            if constexpr (width == 1) {
                                a[w][v][u] = reduce(left[offset]);
                                if constexpr (!Square) {
                                    b[w][v][u] = reduce(right[offset]);
                                }
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                            } else {
                                a[w][v][u] = reduce(Arithmetic::load(left.data() + offset));
                                if constexpr (!Square) {
                                    b[w][v][u] = reduce(Arithmetic::load(right.data() + offset));
                                }
#endif
                            }
                        }
                    }
                }
            }
        };

    public:
        explicit bn254_small_prime_fft_pointwise(std::uint64_t modulus) : arithmetic_(modulus) {
        }

        /** Overwrite left with scale*left*right. Both buffers must contain exactly 12*stride cells. */
        void multiply(std::span<double> left, std::span<const double> right, std::size_t stride, double scale) const {
            check_buffers(left, right, stride, scale);
            if (arithmetic_.fma_supported()) {
                apply<true, false>(left, left, right, scale);
            } else {
                apply<false, false>(left, left, right, scale);
            }
        }

        /** Overwrite values with scale*values^2, loading each tower element once. */
        void square(std::span<double> values, std::size_t stride, double scale) const {
            check_buffers(values, values, stride, scale);
            if (arithmetic_.fma_supported()) {
                apply<true, true>(values, values, values, scale);
            } else {
                apply<false, true>(values, values, values, scale);
            }
        }

        /** Compute scale*(a*b+c*d). Every buffer has 12*stride cells; output may equal any input. */
        void multiply_sum(std::span<double> output, std::span<const double> a, std::span<const double> b,
                          std::span<const double> c, std::span<const double> d, std::size_t stride,
                          double scale) const {
            check_buffers(output, a, stride, scale);
            check_buffers(output, b, stride, scale);
            check_buffers(output, c, stride, scale);
            check_buffers(output, d, stride, scale);
            if (arithmetic_.fma_supported()) {
                apply<true, false, true>(output, a, b, scale, c, d);
            } else {
                apply<false, false, true>(output, a, b, scale, c, d);
            }
        }

    private:
        void check_buffers(std::span<const double> left, std::span<const double> right, std::size_t stride,
                           double scale) const {
            if (left.size() % 12 != 0 || left.size() / 12 != stride || left.size() != right.size()) {
                throw std::invalid_argument("bn254_small_prime_fft_pointwise: invalid buffer dimensions");
            }
            assert(std::isfinite(scale) && std::trunc(scale) == scale &&
                   std::abs(scale) < static_cast<double>(arithmetic_.modulus()) / 2);
        }

        template<bool UseFma, bool Square, bool Sum = false>
        void apply(std::span<double> output, std::span<const double> left, std::span<const double> right, double scale,
                   std::span<const double> extra_left = {}, std::span<const double> extra_right = {}) const {
            // Derive the indexing bound from the buffer extent checked by the public entry points.
            const std::size_t stride = output.size() / 12;
            std::size_t offset = 0;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
                using wide_kernel = small_prime_fft_arithmetic::wide_fma_kernel;
                offset = stride - stride % wide_kernel::lanes;
                if (offset != 0) {
                    const wide_kernel wide(arithmetic_);
                    apply_range<true, Square, Sum>(wide, output, left, right, extra_left, extra_right, 0, offset,
                                                   wide_kernel::broadcast(scale));
                }
#endif
                using packed_kernel = small_prime_fft_arithmetic::packed_fma_kernel;
                const std::size_t packed_end = stride - stride % packed_kernel::lanes;
                if (offset < packed_end) {
                    const packed_kernel packed(arithmetic_);
                    apply_range<true, Square, Sum>(packed, output, left, right, extra_left, extra_right, offset,
                                                   packed_end, packed_kernel::broadcast(scale));
                    offset = packed_end;
                }
            }
#endif
            apply_range<UseFma, Square, Sum>(arithmetic_, output, left, right, extra_left, extra_right, offset, stride,
                                             scale);
        }

        template<bool UseFma, bool Square, bool Sum, typename Arithmetic, typename Value>
        static void apply_range(const Arithmetic &arithmetic, std::span<double> output, std::span<const double> left,
                                std::span<const double> right, std::span<const double> extra_left,
                                std::span<const double> extra_right, std::size_t begin, std::size_t end, Value scale) {
            constexpr std::size_t width = sizeof(Value) / sizeof(double);
            const std::size_t stride = output.size() / 12;
            const kernel<UseFma, Value, Arithmetic> operations {arithmetic};
            using tower_type = typename kernel<UseFma, Value, Arithmetic>::tower_type;
            for (std::size_t i = begin; i < end; i += width) {
                tower_type a, b;
                operations.template load<Square>(a, b, left, right, stride, i);
                tower_type result;
                if constexpr (Square) {
                    result = operations.square(a);
                } else {
                    result = operations.multiply(a, b);
                    if constexpr (Sum) {
                        // Reuse the input slots after the first product. Load both remaining
                        // inputs before storing anything so output may still alias any input.
                        operations.template load<false>(a, b, extra_left, extra_right, stride, i);
                        // Each product coordinate has magnitude < 3*q; their sum is exact below 6*q.
                        // Reduce that sum before the usual inverse-transform normalization.
                        result = operations.add(result, operations.multiply(a, b));
                    }
                }
                for (std::size_t w = 0; w < 2; ++w) {
                    for (std::size_t v = 0; v < 3; ++v) {
                        // Fold inverse-transform normalization into this existing write pass,
                        // pairing the independent scales of the two quadratic coordinates.
                        const std::size_t offset = (6 * w + 2 * v) * stride + i;
                        Value first = result[w][v][0];
                        Value second = result[w][v][1];
                        operations.multiply_pair(first, second, scale, scale);
                        if constexpr (width == 1) {
                            output[offset] = first;
                            output[offset + stride] = second;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                        } else {
                            Arithmetic::store(output.data() + offset, first);
                            Arithmetic::store(output.data() + offset + stride, second);
#endif
                        }
                    }
                }
            }
        }

        small_prime_fft_arithmetic arithmetic_;
    };
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
#pragma GCC diagnostic pop
#endif

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_DETAIL_BN254_SMALL_PRIME_FFT_POINTWISE_HPP
