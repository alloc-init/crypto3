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

#ifndef CRYPTO3_MATH_DETAIL_SMALL_PRIME_FFT_PLAN_HPP
#define CRYPTO3_MATH_DETAIL_SMALL_PRIME_FFT_PLAN_HPP

#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include <boost/config.hpp>

#include <nil/crypto3/math/polynomial/backends/detail/small_prime_fft/small_prime_fft_arithmetic.hpp>
#include <nil/crypto3/math/polynomial/backends/detail/small_prime_fft/small_prime_fft_params.hpp>
#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>
#include <nil/crypto3/multiprecision/inverse.hpp>

#ifndef CRYPTO3_SMALL_PRIME_FFT_RADIX
#define CRYPTO3_SMALL_PRIME_FFT_RADIX 4
#endif
static_assert(CRYPTO3_SMALL_PRIME_FFT_RADIX == 4 || CRYPTO3_SMALL_PRIME_FFT_RADIX == 8 ||
                  CRYPTO3_SMALL_PRIME_FFT_RADIX == 16,
              "small-prime FFT stage radix must be 4, 8 or 16");

namespace nil::crypto3::math::detail {

    /**
     * In-place, power-of-two transforms for one prime from the catalogue. The forward transform
     * maps natural-order coefficients to evaluations at omega^k in bit-reversed index order.
     * Both directions use this order, so convolution needs no permutation passes.
     *
     * A reusable table of max(1,N/2) centered twiddles serves both directions. Radix-4 recursion
     * fuses pairs of layers, with radix-2 leaves for odd depths and special cases for unit twiddles.
     * AVX2/FMA builds fuse 8-, 16- and 32-point leaves; AVX-512 adds paired 32-/64-point
     * leaves and eight-lane stages. CRYPTO3_SMALL_PRIME_FFT_RADIX selects 4 (default), 8
     * or 16 for large SIMD stages; the scalar/integer fallback keeps the radix-4 schedule.
     * Transform calls allocate no buffers; recursion uses O(log N) stack space.
     *
     * @pre Inputs are integer-valued doubles in (-3*q,3*q). The calling thread's rounding mode
     * must remain unchanged during a call. Each call selects FMA or integer arithmetic once.
     * The forward output is in (-3*q,3*q); the normalized inverse output is in [0,q).
     * This operation has no constant-time guarantee.
     */
    class small_prime_fft_plan {
    public:
        small_prime_fft_plan(std::size_t prime_index, std::size_t size) :
            arithmetic_(checked_prime(prime_index, size).modulus), size_(size),
            twiddles_(std::max(std::size_t(1), size / 2), 1.0) {
            using wide_type = boost::multiprecision::uint128_modular_t;
            const auto &prime = small_prime_fft_primes[prime_index];
            wide_type root = prime.root;
            const unsigned depth = std::countr_zero(size_);
            for (unsigned level = prime.two_adicity; level > depth; --level) {
                root = root * root % prime.modulus;
            }
            omega_ = root.convert_to<std::uint64_t>();
            size_inverse_ = centered(
                static_cast<double>(boost::multiprecision::inverse_mod(wide_type(size_), wide_type(prime.modulus))
                                        .convert_to<std::uint64_t>()));

            // At a layer with 2^d blocks, block j uses omega_(2^(d+1))^bitreverse(j,d).
            // Each new table segment is the preceding prefix times a primitive root of twice
            // the previous order: [1], [i], [w8,i*w8], [w16,i*w16,w8*w16,i*w8*w16], ... .
            for (std::size_t half = size_ / 4; half != 0; half >>= 1) {
                twiddles_[half] = centered(static_cast<double>(root.convert_to<std::uint64_t>()));
                root = root * root % prime.modulus;
            }
            for (std::size_t half = 1; half < size_ / 2; half <<= 1) {
                const double step = twiddles_[half];
                for (std::size_t i = 1; i < half; ++i) {
                    twiddles_[half + i] = centered(arithmetic_.multiply_integer(twiddles_[i], step));
                }
            }
        }

        std::size_t size() const {
            return size_;
        }

        std::uint64_t modulus() const {
            return arithmetic_.modulus();
        }

        std::uint64_t omega() const {
            return omega_;
        }

        /** Exact centered representative of 1/N modulo q, also usable to fuse scaling into pointwise work. */
        double size_inverse() const {
            return size_inverse_;
        }

        /** Require exactly N entries; invalid dimensions throw before any writes. N=1 is the identity. */
        void fft(std::span<double> values) const {
            check_buffer(values);
            if (arithmetic_.fma_supported()) {
                forward_recursive<true>(values.data(), size_, 0);
            } else {
                forward_recursive<false>(values.data(), size_, 0);
            }
        }

        /**
         * Require exactly N entries in the forward spectrum order. By default, scale by 1/N and
         * return canonical residues. With normalize=false, return N times the inverse in (-2*q,2*q),
         * allowing the caller to fuse the scaling elsewhere. Invalid dimensions throw before writes.
         */
        void inverse_fft(std::span<double> values, bool normalize = true) const {
            check_buffer(values);
            if (arithmetic_.fma_supported()) {
                inverse_impl<true>(values, normalize);
            } else {
                inverse_impl<false>(values, normalize);
            }
        }

        /**
         * Unscaled inverse for pointwise products whose integer representatives are in (-q,q).
         * The caller has already folded 1/N into its pointwise scale. This stricter input bound
         * permits omitting the first leaf reductions; the output still lies in (-2*q,2*q).
         * Dimension errors throw before writes, as in inverse_fft().
         */
        void inverse_fft_pointwise(std::span<double> values) const {
            check_buffer<1>(values);
            if (arithmetic_.fma_supported()) {
                inverse_impl<true, true>(values, false);
            } else {
                inverse_impl<false, true>(values, false);
            }
        }

    private:
        static const small_prime_fft_prime &checked_prime(std::size_t prime_index, std::size_t size) {
            if (prime_index >= small_prime_fft_primes.size() || !std::has_single_bit(size) ||
                std::countr_zero(size) > small_prime_fft_primes[prime_index].two_adicity) {
                throw std::invalid_argument("small_prime_fft_plan: invalid prime index or transform size");
            }
            return small_prime_fft_primes[prime_index];
        }

        template<unsigned Bound = 3>
        void check_buffer(std::span<const double> values) const {
            if (values.size() != size_ || twiddles_.empty()) {
                throw std::invalid_argument("small_prime_fft_plan: invalid buffer dimensions or moved-from plan");
            }
#ifndef NDEBUG
            const double limit = Bound * static_cast<double>(modulus());
            for (double value : values) {
                assert(std::isfinite(value) && std::trunc(value) == value && std::abs(value) < limit);
            }
#endif
        }

        double centered(double value) const {
            const double q = static_cast<double>(modulus());
            return value > q / 2 ? value - q : (value < -q / 2 ? value + q : value);
        }

        BOOST_FORCEINLINE double inverse_twiddle(std::size_t index) const {
            if (index == 0) {
                return 1.0;
            }
            // Within [b,2*b), reflected exponents add to half the root order, hence w^-1=-w_reflected.
            const std::size_t base = std::bit_floor(index);
            return -twiddles_[base + (2 * base - 1 - index)];
        }

        // Keep arithmetic selection outside the butterfly loops, with compile-time dispatch inside them.
        template<bool UseFma, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static Value reduce(const Arithmetic &arithmetic, Value value) {
            if constexpr (UseFma) {
                return arithmetic.reduce_fma(value);
            } else {
                return arithmetic.reduce_integer(value);
            }
        }

        template<bool UseFma, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static Value multiply(const Arithmetic &arithmetic, Value left, Value right) {
            if constexpr (UseFma) {
                return arithmetic.multiply_fma(left, right);
            } else {
                return arithmetic.multiply_integer(left, right);
            }
        }

        template<bool UseFma, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void multiply_pair(const Arithmetic &arithmetic, Value &a, Value &b, Value wa,
                                                    Value wb) {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma && requires { arithmetic.multiply_fma_pair(a, b, wa, wb); }) {
                arithmetic.multiply_fma_pair(a, b, wa, wb);
            } else
#endif
            {
                a = multiply<UseFma>(arithmetic, a, wa);
                b = multiply<UseFma>(arithmetic, b, wb);
            }
        }

        template<bool UseFma, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void forward_butterfly(const Arithmetic &arithmetic, Value &v0, Value &v1, Value &v2,
                                                        Value &v3, Value w2, Value w, Value iw) {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
            const Value a0 = reduce<UseFma>(arithmetic, v0);
            const Value a1 = v1;
            Value a2 = v2;
            Value a3 = v3;
            if constexpr (UnitTwiddle) {
                a2 = reduce<UseFma>(arithmetic, a2);
                a3 = reduce<UseFma>(arithmetic, a3);
            } else {
                multiply_pair<UseFma>(arithmetic, a2, a3, w2, w2);
            }
            const Value even_sum = a0 + a2;
            const Value even_difference = a0 - a2;
            Value odd_sum = a1 + a3;
            Value odd_difference = a1 - a3;
            if constexpr (UnitTwiddle) {
                odd_sum = reduce<UseFma>(arithmetic, odd_sum);
                odd_difference = multiply<UseFma>(arithmetic, odd_difference, iw);
            } else {
                multiply_pair<UseFma>(arithmetic, odd_sum, odd_difference, w, iw);
            }
            // Input magnitudes < 3*q and weights < q/2 imply |a2|,|a3| < q and
            // |a1 +/- a3| < 4*q. Every product is < 2*q^2, so the odd terms are < q.
            // Thus all outputs are < 3*q, with every integer sum exact below 2^53.
            v0 = even_sum + odd_sum;
            v1 = even_sum - odd_sum;
            v2 = even_difference + odd_difference;
            v3 = even_difference - odd_difference;
        }

#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
        BOOST_FORCEINLINE static void transpose4(__m256d &v0, __m256d &v1, __m256d &v2, __m256d &v3) {
            const __m256d a = _mm256_unpacklo_pd(v0, v1);
            const __m256d b = _mm256_unpackhi_pd(v0, v1);
            const __m256d c = _mm256_unpacklo_pd(v2, v3);
            const __m256d d = _mm256_unpackhi_pd(v2, v3);
            v0 = _mm256_permute2f128_pd(a, c, 0x20);
            v1 = _mm256_permute2f128_pd(b, d, 0x20);
            v2 = _mm256_permute2f128_pd(a, c, 0x31);
            v3 = _mm256_permute2f128_pd(b, d, 0x31);
        }

        // Aligned table indices, with unaligned loads: the owning vector needs no special allocator.
        template<bool Inverse>
        BOOST_FORCEINLINE __m256d twiddle4(std::size_t index) const {
            if constexpr (!Inverse) {
                return _mm256_loadu_pd(twiddles_.data() + index);
            } else {
                if (index == 0) {
                    return _mm256_setr_pd(1.0, inverse_twiddle(1), inverse_twiddle(2), inverse_twiddle(3));
                }
                // Four consecutive entries stay in one [b,2*b) segment. Reflect, reverse, negate.
                const std::size_t base = std::bit_floor(index);
                const __m256d reflected = _mm256_loadu_pd(twiddles_.data() + 3 * base - index - 4);
                return _mm256_sub_pd(_mm256_setzero_pd(), _mm256_permute4x64_pd(reflected, 0x1b));
            }
        }

        template<bool Inverse>
        BOOST_FORCEINLINE void leaf_twiddles(std::size_t node, __m256d &w2, __m256d &w, __m256d &iw) const {
            w2 = twiddle4<Inverse>(4 * node);
            const __m256d low = twiddle4<Inverse>(8 * node);
            const __m256d high = twiddle4<Inverse>(8 * node + 4);
            w = _mm256_permute4x64_pd(_mm256_unpacklo_pd(low, high), 0xd8);
            iw = _mm256_permute4x64_pd(_mm256_unpackhi_pd(low, high), 0xd8);
        }

        template<bool UnitTwiddle>
        BOOST_FORCEINLINE void forward16(const small_prime_fft_arithmetic::packed_fma_kernel &packed, __m256d &v0,
                                         __m256d &v1, __m256d &v2, __m256d &v3, std::size_t node) const {
            forward_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, _mm256_set1_pd(twiddles_[node]),
                                                 _mm256_set1_pd(twiddles_[2 * node]),
                                                 _mm256_set1_pd(twiddles_[2 * node + 1]));
            // The four child transforms now run across SIMD lanes, without scalar leaves or reloads.
            transpose4(v0, v1, v2, v3);
            __m256d w2, w, iw;
            leaf_twiddles<false>(node, w2, w, iw);
            forward_butterfly<true, false>(packed, v0, v1, v2, v3, w2, w, iw);
            transpose4(v0, v1, v2, v3);
        }

        template<bool UnitTwiddle, bool BoundedInput>
        BOOST_FORCEINLINE void inverse16(const small_prime_fft_arithmetic::packed_fma_kernel &packed, __m256d &v0,
                                         __m256d &v1, __m256d &v2, __m256d &v3, std::size_t node) const {
            transpose4(v0, v1, v2, v3);
            __m256d w2, w, iw;
            leaf_twiddles<true>(node, w2, w, iw);
            inverse_butterfly<true, false>(packed, v0, v1, v2, v3, w2, w, iw, !BoundedInput);
            transpose4(v0, v1, v2, v3);
            inverse_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, _mm256_set1_pd(inverse_twiddle(node)),
                                                 _mm256_set1_pd(inverse_twiddle(2 * node)),
                                                 _mm256_set1_pd(inverse_twiddle(2 * node + 1)), false);
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput = false>
        void base16(double *values, std::size_t node) const {
            const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
            __m256d v0 = _mm256_loadu_pd(values);
            __m256d v1 = _mm256_loadu_pd(values + 4);
            __m256d v2 = _mm256_loadu_pd(values + 8);
            __m256d v3 = _mm256_loadu_pd(values + 12);
            if constexpr (Inverse) {
                inverse16<UnitTwiddle, BoundedInput>(packed, v0, v1, v2, v3, node);
            } else {
                forward16<UnitTwiddle>(packed, v0, v1, v2, v3, node);
            }
            _mm256_storeu_pd(values, v0);
            _mm256_storeu_pd(values + 4, v1);
            _mm256_storeu_pd(values + 8, v2);
            _mm256_storeu_pd(values + 12, v3);
        }

        template<bool Inverse, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void radix2_packed(const Arithmetic &packed, Value &left, Value &right, Value weight) {
            if constexpr (Inverse) {
                // Both subtransforms already return magnitudes < 2*q. No input reductions are needed.
                const Value sum = left + right;
                const Value difference = left - right;
                left = packed.reduce_fma(sum);
                right = UnitTwiddle ? packed.reduce_fma(difference) : packed.multiply_fma(difference, weight);
            } else {
                const Value a = packed.reduce_fma(left);
                const Value b = UnitTwiddle ? packed.reduce_fma(right) : packed.multiply_fma(right, weight);
                left = a + b;
                right = a - b;
            }
        }

        template<bool Inverse, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void stage4(const Arithmetic &arithmetic, Value &v0, Value &v1, Value &v2, Value &v3,
                                             const double *weights) {
            const Value w2 = Arithmetic::broadcast(weights[0]);
            const Value w = Arithmetic::broadcast(weights[1]);
            const Value iw = Arithmetic::broadcast(weights[2]);
            if constexpr (Inverse) {
                inverse_butterfly<true, UnitTwiddle>(arithmetic, v0, v1, v2, v3, w2, w, iw, false);
            } else {
                forward_butterfly<true, UnitTwiddle>(arithmetic, v0, v1, v2, v3, w2, w, iw);
            }
        }

        template<bool Inverse>
        BOOST_FORCEINLINE void stage_weights(double *weights, std::size_t node) const {
            weights[0] = Inverse ? inverse_twiddle(node) : twiddles_[node];
            weights[1] = Inverse ? inverse_twiddle(2 * node) : twiddles_[2 * node];
            weights[2] = Inverse ? inverse_twiddle(2 * node + 1) : twiddles_[2 * node + 1];
        }

        template<bool Inverse, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE void stage4(const Arithmetic &arithmetic, Value &v0, Value &v1, Value &v2, Value &v3,
                                      std::size_t node) const {
            double weights[3];
            stage_weights<Inverse>(weights, node);
            stage4<Inverse, UnitTwiddle>(arithmetic, v0, v1, v2, v3, weights);
        }

        template<unsigned Radix, bool Inverse, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void fused_outer(const Arithmetic &arithmetic, Value *v, const double *weights) {
            if constexpr (Radix == 8) {
                const Value weight = Arithmetic::broadcast(weights[0]);
                radix2_packed<Inverse, UnitTwiddle>(arithmetic, v[0], v[4], weight);
                radix2_packed<Inverse, UnitTwiddle>(arithmetic, v[1], v[5], weight);
                radix2_packed<Inverse, UnitTwiddle>(arithmetic, v[2], v[6], weight);
                radix2_packed<Inverse, UnitTwiddle>(arithmetic, v[3], v[7], weight);
            } else {
                stage4<Inverse, UnitTwiddle>(arithmetic, v[0], v[4], v[8], v[12], weights);
                stage4<Inverse, UnitTwiddle>(arithmetic, v[1], v[5], v[9], v[13], weights);
                stage4<Inverse, UnitTwiddle>(arithmetic, v[2], v[6], v[10], v[14], weights);
                stage4<Inverse, UnitTwiddle>(arithmetic, v[3], v[7], v[11], v[15], weights);
            }
        }

        // Fuse three/four binary layers without intermediate array passes. Every constituent
        // butterfly retains the forward <3*q and inverse <2*q bounds; no wider lazy range is used.
        // UnitTwiddle is valid only at node zero, hence only the first child inherits it.
        template<unsigned Radix, bool Inverse, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void fused_stage(const Arithmetic &arithmetic, Value *v, const double *weights) {
            static_assert(Radix == 8 || Radix == 16);
            if constexpr (!Inverse) {
                fused_outer<Radix, Inverse, UnitTwiddle>(arithmetic, v, weights);
            }
            constexpr std::size_t outer_weights = Radix == 8 ? 1 : 3;
            stage4<Inverse, UnitTwiddle>(arithmetic, v[0], v[1], v[2], v[3], weights + outer_weights);
            stage4<Inverse, false>(arithmetic, v[4], v[5], v[6], v[7], weights + outer_weights + 3);
            if constexpr (Radix == 16) {
                stage4<Inverse, false>(arithmetic, v[8], v[9], v[10], v[11], weights + outer_weights + 6);
                stage4<Inverse, false>(arithmetic, v[12], v[13], v[14], v[15], weights + outer_weights + 9);
            }
            if constexpr (Inverse) {
                fused_outer<Radix, Inverse, UnitTwiddle>(arithmetic, v, weights);
            }
        }

        template<unsigned Radix, bool Inverse, bool UnitTwiddle>
        void large_radix(double *values, std::size_t stride, std::size_t node) const {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
            using arithmetic_type = small_prime_fft_arithmetic::wide_fma_kernel;
#else
            using arithmetic_type = small_prime_fft_arithmetic::packed_fma_kernel;
#endif
            const arithmetic_type arithmetic(arithmetic_);
            // Dispatch requires stride >=16, a power of two divisible by the SIMD width.
            assert(stride >= 16 && stride % arithmetic_type::lanes == 0);
            // The node is fixed throughout this coefficient loop. Resolve reflections and
            // table lookups once, retaining scalar weights to avoid forcing Radix-1 live vectors.
            double weights[Radix - 1];
            if constexpr (Radix == 8) {
                weights[0] = Inverse ? inverse_twiddle(node) : twiddles_[node];
            } else {
                stage_weights<Inverse>(weights, node);
            }
            constexpr std::size_t outer_weights = Radix == 8 ? 1 : 3;
#pragma GCC unroll 4
            for (std::size_t child = 0; child < Radix / 4; ++child) {
                stage_weights<Inverse>(weights + outer_weights + 3 * child, Radix / 4 * node + child);
            }
            for (std::size_t i = 0; i < stride; i += arithmetic_type::lanes) {
                arithmetic_type::value_type v[Radix];
#pragma GCC unroll 16
                for (std::size_t part = 0; part < Radix; ++part) {
                    v[part] = arithmetic_type::load(values + part * stride + i);
                }
                fused_stage<Radix, Inverse, UnitTwiddle>(arithmetic, v, weights);
#pragma GCC unroll 16
                for (std::size_t part = 0; part < Radix; ++part) {
                    arithmetic_type::store(values + part * stride + i, v[part]);
                }
            }
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput = false>
        void base32(double *values, std::size_t node) const {
            const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
            __m256d v0 = _mm256_loadu_pd(values);
            __m256d v1 = _mm256_loadu_pd(values + 4);
            __m256d v2 = _mm256_loadu_pd(values + 8);
            __m256d v3 = _mm256_loadu_pd(values + 12);
            __m256d v4 = _mm256_loadu_pd(values + 16);
            __m256d v5 = _mm256_loadu_pd(values + 20);
            __m256d v6 = _mm256_loadu_pd(values + 24);
            __m256d v7 = _mm256_loadu_pd(values + 28);
            const __m256d weight = _mm256_set1_pd(Inverse ? inverse_twiddle(node) : twiddles_[node]);
            // One radix-2 layer and two register-resident 16-point transforms cover all five layers.
            if constexpr (Inverse) {
                inverse16<UnitTwiddle, BoundedInput>(packed, v0, v1, v2, v3, 2 * node);
                inverse16<false, BoundedInput>(packed, v4, v5, v6, v7, 2 * node + 1);
            }
            radix2_packed<Inverse, UnitTwiddle>(packed, v0, v4, weight);
            radix2_packed<Inverse, UnitTwiddle>(packed, v1, v5, weight);
            radix2_packed<Inverse, UnitTwiddle>(packed, v2, v6, weight);
            radix2_packed<Inverse, UnitTwiddle>(packed, v3, v7, weight);
            if constexpr (!Inverse) {
                forward16<UnitTwiddle>(packed, v0, v1, v2, v3, 2 * node);
                forward16<false>(packed, v4, v5, v6, v7, 2 * node + 1);
            }
            _mm256_storeu_pd(values, v0);
            _mm256_storeu_pd(values + 4, v1);
            _mm256_storeu_pd(values + 8, v2);
            _mm256_storeu_pd(values + 12, v3);
            _mm256_storeu_pd(values + 16, v4);
            _mm256_storeu_pd(values + 20, v5);
            _mm256_storeu_pd(values + 24, v6);
            _mm256_storeu_pd(values + 28, v7);
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput = false>
        void base8(double *values, std::size_t node) const {
            const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
            const __m256d low = _mm256_loadu_pd(values);
            const __m256d high = _mm256_loadu_pd(values + 4);
            __m256d v0, v1, v2, v3;
            if constexpr (Inverse) {
                __m256d a = _mm256_permute4x64_pd(_mm256_unpacklo_pd(low, high), 0xd8);
                __m256d b = _mm256_permute4x64_pd(_mm256_unpackhi_pd(low, high), 0xd8);
                if constexpr (!BoundedInput) {
                    a = packed.reduce_fma(a);
                    b = packed.reduce_fma(b);
                }
                const __m256d sum = packed.reduce_fma(a + b);
                const __m256d difference = packed.multiply_fma(a - b, twiddle4<true>(4 * node));
                const __m256d even = _mm256_unpacklo_pd(sum, difference);
                const __m256d odd = _mm256_unpackhi_pd(sum, difference);
                v0 = _mm256_permute2f128_pd(even, even, 0x00);
                v1 = _mm256_permute2f128_pd(odd, odd, 0x00);
                v2 = _mm256_permute2f128_pd(even, even, 0x11);
                v3 = _mm256_permute2f128_pd(odd, odd, 0x11);
                inverse_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, _mm256_set1_pd(inverse_twiddle(node)),
                                                     _mm256_set1_pd(inverse_twiddle(2 * node)),
                                                     _mm256_set1_pd(inverse_twiddle(2 * node + 1)), false);
                _mm256_storeu_pd(values, _mm256_permute2f128_pd(v0, v1, 0x20));
                _mm256_storeu_pd(values + 4, _mm256_permute2f128_pd(v2, v3, 0x20));
            } else {
                // Two independent columns occupy the lower lanes; the duplicated upper lanes are discarded.
                v0 = _mm256_permute2f128_pd(low, low, 0x00);
                v1 = _mm256_permute2f128_pd(low, low, 0x11);
                v2 = _mm256_permute2f128_pd(high, high, 0x00);
                v3 = _mm256_permute2f128_pd(high, high, 0x11);
                forward_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, _mm256_set1_pd(twiddles_[node]),
                                                     _mm256_set1_pd(twiddles_[2 * node]),
                                                     _mm256_set1_pd(twiddles_[2 * node + 1]));
                const __m256d a = packed.reduce_fma(
                    _mm256_permute2f128_pd(_mm256_unpacklo_pd(v0, v1), _mm256_unpacklo_pd(v2, v3), 0x20));
                const __m256d b = packed.multiply_fma(
                    _mm256_permute2f128_pd(_mm256_unpackhi_pd(v0, v1), _mm256_unpackhi_pd(v2, v3), 0x20),
                    twiddle4<false>(4 * node));
                const __m256d even = _mm256_unpacklo_pd(a + b, a - b);
                const __m256d odd = _mm256_unpackhi_pd(a + b, a - b);
                _mm256_storeu_pd(values, _mm256_permute2f128_pd(even, odd, 0x20));
                _mm256_storeu_pd(values + 4, _mm256_permute2f128_pd(even, odd, 0x31));
            }
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput = false>
        void basecase(double *values, std::size_t length, std::size_t node) const {
            if (length == 16) {
                base16<Inverse, UnitTwiddle, BoundedInput>(values, node);
            } else if (length == 32) {
                base32<Inverse, UnitTwiddle, BoundedInput>(values, node);
            } else {
                base8<Inverse, UnitTwiddle, BoundedInput>(values, node);
            }
        }

#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
        BOOST_FORCEINLINE static __m512d join4(__m256d low, __m256d high) {
            return _mm512_insertf64x4(_mm512_castpd256_pd512(low), high, 1);
        }

        // Two independent 4x4 transposes, one in each 256-bit half. The stored spectrum
        // remains bit-reversed; these permutations only change the lanes during a leaf.
        BOOST_FORCEINLINE static void transpose4(__m512d &v0, __m512d &v1, __m512d &v2, __m512d &v3) {
            const __m512d a = _mm512_unpacklo_pd(v0, v1);
            const __m512d b = _mm512_unpackhi_pd(v0, v1);
            const __m512d c = _mm512_unpacklo_pd(v2, v3);
            const __m512d d = _mm512_unpackhi_pd(v2, v3);
            const __m512i low = _mm512_setr_epi64(0, 1, 8, 9, 4, 5, 12, 13);
            const __m512i high = _mm512_setr_epi64(2, 3, 10, 11, 6, 7, 14, 15);
            v0 = _mm512_permutex2var_pd(a, low, c);
            v1 = _mm512_permutex2var_pd(b, low, d);
            v2 = _mm512_permutex2var_pd(a, high, c);
            v3 = _mm512_permutex2var_pd(b, high, d);
        }

        template<bool Inverse>
        BOOST_FORCEINLINE __m512d paired_twiddle(std::size_t low, std::size_t high) const {
            return join4(_mm256_set1_pd(Inverse ? inverse_twiddle(low) : twiddles_[low]),
                         _mm256_set1_pd(Inverse ? inverse_twiddle(high) : twiddles_[high]));
        }

        // Keep the complete leaf in one function: outlining this step forces live vectors
        // through memory across the call. Recursive transform boundaries remain uninlined.
        template<bool Inverse>
        BOOST_FORCEINLINE void paired_outer(const small_prime_fft_arithmetic::wide_fma_kernel &arithmetic, __m512d &v0,
                                            __m512d &v1, __m512d &v2, __m512d &v3, std::size_t node) const {
            const __m512d w2 = paired_twiddle<Inverse>(node, node + 1);
            const __m512d w = paired_twiddle<Inverse>(2 * node, 2 * node + 2);
            const __m512d iw = paired_twiddle<Inverse>(2 * node + 1, 2 * node + 3);
            if constexpr (Inverse) {
                inverse_butterfly<true, false>(arithmetic, v0, v1, v2, v3, w2, w, iw, false);
            } else {
                forward_butterfly<true, false>(arithmetic, v0, v1, v2, v3, w2, w, iw);
            }
        }

        template<bool Inverse, bool BoundedInput>
        BOOST_FORCEINLINE void base16_pair(const small_prime_fft_arithmetic::wide_fma_kernel &arithmetic, __m512d &v0,
                                           __m512d &v1, __m512d &v2, __m512d &v3, std::size_t node) const {
            // Group matching rows of two adjacent 16-point transforms in the two vector halves.
            const __m512d a = _mm512_shuffle_f64x2(v0, v2, 0x44);
            const __m512d b = _mm512_shuffle_f64x2(v0, v2, 0xee);
            const __m512d c = _mm512_shuffle_f64x2(v1, v3, 0x44);
            const __m512d d = _mm512_shuffle_f64x2(v1, v3, 0xee);
            v0 = a;
            v1 = b;
            v2 = c;
            v3 = d;
            if constexpr (!Inverse) {
                paired_outer<false>(arithmetic, v0, v1, v2, v3, node);
            }
            transpose4(v0, v1, v2, v3);
            __m256d lo_w2, lo_w, lo_iw, hi_w2, hi_w, hi_iw;
            leaf_twiddles<Inverse>(node, lo_w2, lo_w, lo_iw);
            leaf_twiddles<Inverse>(node + 1, hi_w2, hi_w, hi_iw);
            if constexpr (Inverse) {
                inverse_butterfly<true, false>(arithmetic, v0, v1, v2, v3, join4(lo_w2, hi_w2), join4(lo_w, hi_w),
                                               join4(lo_iw, hi_iw), !BoundedInput);
            } else {
                forward_butterfly<true, false>(arithmetic, v0, v1, v2, v3, join4(lo_w2, hi_w2), join4(lo_w, hi_w),
                                               join4(lo_iw, hi_iw));
            }
            transpose4(v0, v1, v2, v3);
            if constexpr (Inverse) {
                paired_outer<true>(arithmetic, v0, v1, v2, v3, node);
            }
            // Restore contiguous bit-reversed order in the same four input/output vectors.
            const __m512d out0 = _mm512_shuffle_f64x2(v0, v1, 0x44);
            const __m512d out1 = _mm512_shuffle_f64x2(v2, v3, 0x44);
            const __m512d out2 = _mm512_shuffle_f64x2(v0, v1, 0xee);
            const __m512d out3 = _mm512_shuffle_f64x2(v2, v3, 0xee);
            v0 = out0;
            v1 = out1;
            v2 = out2;
            v3 = out3;
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput>
        void base32_wide(double *values, std::size_t node) const {
            const small_prime_fft_arithmetic::wide_fma_kernel arithmetic(arithmetic_);
            __m512d v0 = _mm512_loadu_pd(values);
            __m512d v1 = _mm512_loadu_pd(values + 8);
            __m512d v2 = _mm512_loadu_pd(values + 16);
            __m512d v3 = _mm512_loadu_pd(values + 24);
            const __m512d weight = _mm512_set1_pd(Inverse ? inverse_twiddle(node) : twiddles_[node]);
            if constexpr (Inverse) {
                base16_pair<true, BoundedInput>(arithmetic, v0, v1, v2, v3, 2 * node);
            }
            radix2_packed<Inverse, UnitTwiddle>(arithmetic, v0, v2, weight);
            radix2_packed<Inverse, UnitTwiddle>(arithmetic, v1, v3, weight);
            if constexpr (!Inverse) {
                base16_pair<false, false>(arithmetic, v0, v1, v2, v3, 2 * node);
            }
            _mm512_storeu_pd(values, v0);
            _mm512_storeu_pd(values + 8, v1);
            _mm512_storeu_pd(values + 16, v2);
            _mm512_storeu_pd(values + 24, v3);
        }

        template<bool Inverse, bool UnitTwiddle, bool BoundedInput>
        void base64_wide(double *values, std::size_t node) const {
            const small_prime_fft_arithmetic::wide_fma_kernel arithmetic(arithmetic_);
            __m512d v[8];
#pragma GCC unroll 8
            for (std::size_t i = 0; i < 8; ++i) {
                v[i] = _mm512_loadu_pd(values + 8 * i);
            }
            if constexpr (!Inverse) {
                stage4<false, UnitTwiddle>(arithmetic, v[0], v[2], v[4], v[6], node);
                stage4<false, UnitTwiddle>(arithmetic, v[1], v[3], v[5], v[7], node);
            }
            base16_pair<Inverse, BoundedInput>(arithmetic, v[0], v[1], v[2], v[3], 4 * node);
            base16_pair<Inverse, BoundedInput>(arithmetic, v[4], v[5], v[6], v[7], 4 * node + 2);
            if constexpr (Inverse) {
                stage4<true, UnitTwiddle>(arithmetic, v[0], v[2], v[4], v[6], node);
                stage4<true, UnitTwiddle>(arithmetic, v[1], v[3], v[5], v[7], node);
            }
#pragma GCC unroll 8
            for (std::size_t i = 0; i < 8; ++i) {
                _mm512_storeu_pd(values + 8 * i, v[i]);
            }
        }

        template<bool Inverse, bool BoundedInput>
        void wide_basecase(double *values, std::size_t length, std::size_t node) const {
            if (length == 32) {
                if (node == 0) {
                    base32_wide<Inverse, true, BoundedInput>(values, node);
                } else {
                    base32_wide<Inverse, false, BoundedInput>(values, node);
                }
            } else {
                if (node == 0) {
                    base64_wide<Inverse, true, BoundedInput>(values, node);
                } else {
                    base64_wide<Inverse, false, BoundedInput>(values, node);
                }
            }
        }
#endif
#endif

        template<bool UseFma, bool UnitTwiddle>
        void forward_radix4(double *values, std::size_t quarter, std::size_t node, std::size_t first,
                            std::size_t last) const {
            const double w2 = twiddles_[node];
            const double w = twiddles_[2 * node];
            const double iw = twiddles_[2 * node + 1];
            std::size_t i = first;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
            if constexpr (UseFma) {
                if (quarter >= 8) {
                    const small_prime_fft_arithmetic::wide_fma_kernel packed(arithmetic_);
                    const __m512d packed_w2 = _mm512_set1_pd(w2);
                    const __m512d packed_w = _mm512_set1_pd(w);
                    const __m512d packed_iw = _mm512_set1_pd(iw);
#pragma GCC unroll 2
                    for (; last - i >= 8; i += 8) {
                        __m512d v0 = _mm512_loadu_pd(values + i);
                        __m512d v1 = _mm512_loadu_pd(values + i + quarter);
                        __m512d v2 = _mm512_loadu_pd(values + i + 2 * quarter);
                        __m512d v3 = _mm512_loadu_pd(values + i + 3 * quarter);
                        forward_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, packed_w2, packed_w, packed_iw);
                        _mm512_storeu_pd(values + i, v0);
                        _mm512_storeu_pd(values + i + quarter, v1);
                        _mm512_storeu_pd(values + i + 2 * quarter, v2);
                        _mm512_storeu_pd(values + i + 3 * quarter, v3);
                    }
                }
            }
#endif
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
                if (quarter >= 4) {
                    const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
                    const __m256d packed_w2 = _mm256_set1_pd(w2);
                    const __m256d packed_w = _mm256_set1_pd(w);
                    const __m256d packed_iw = _mm256_set1_pd(iw);
#pragma GCC unroll 2
                    for (; last - i >= 4; i += 4) {
                        __m256d v0 = _mm256_loadu_pd(values + i);
                        __m256d v1 = _mm256_loadu_pd(values + i + quarter);
                        __m256d v2 = _mm256_loadu_pd(values + i + 2 * quarter);
                        __m256d v3 = _mm256_loadu_pd(values + i + 3 * quarter);
                        forward_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, packed_w2, packed_w, packed_iw);
                        _mm256_storeu_pd(values + i, v0);
                        _mm256_storeu_pd(values + i + quarter, v1);
                        _mm256_storeu_pd(values + i + 2 * quarter, v2);
                        _mm256_storeu_pd(values + i + 3 * quarter, v3);
                    }
                }
            }
#endif
            for (; i < last; ++i) {
                forward_butterfly<UseFma, UnitTwiddle>(arithmetic_, values[i], values[i + quarter],
                                                       values[i + 2 * quarter], values[i + 3 * quarter], w2, w, iw);
            }
        }

#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
        template<bool Inverse, bool BoundedInput>
        void four_layer_block(double *values, std::size_t length, std::size_t node) const {
            const std::size_t quarter = length / 4;
            const std::size_t sixteenth = length / 16;
            // Sixteen strips of 64 doubles occupy 8 KiB. Complete both radix-4 passes on
            // each strip group while it is hot, before recursing into its sixteen children.
            // This only changes visitation order; every butterfly retains the same inputs and bounds.
            constexpr std::size_t strip_size = 64;
            const auto pass = [&](double *base, std::size_t stride, std::size_t j, std::size_t first,
                                  std::size_t last) {
                if constexpr (Inverse) {
                    if (j == 0) {
                        inverse_radix4<true, true, BoundedInput>(base, stride, j, first, last);
                    } else {
                        inverse_radix4<true, false, BoundedInput>(base, stride, j, first, last);
                    }
                } else {
                    if (j == 0) {
                        forward_radix4<true, true>(base, stride, j, first, last);
                    } else {
                        forward_radix4<true, false>(base, stride, j, first, last);
                    }
                }
            };
            for (std::size_t first = 0; first < sixteenth; first += strip_size) {
                const std::size_t last = std::min(first + strip_size, sixteenth);
                if constexpr (Inverse) {
                    for (std::size_t child = 0; child < 4; ++child) {
                        pass(values + child * quarter, sixteenth, 4 * node + child, first, last);
                    }
                }
                for (std::size_t part = 0; part < 4; ++part) {
                    pass(values, quarter, node, first + part * sixteenth, last + part * sixteenth);
                }
                if constexpr (!Inverse) {
                    for (std::size_t child = 0; child < 4; ++child) {
                        pass(values + child * quarter, sixteenth, 4 * node + child, first, last);
                    }
                }
            }
        }
#endif

        template<bool UseFma>
        void forward_recursive(double *values, std::size_t length, std::size_t node) const {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
                if (length == 32 || length == 64) {
                    wide_basecase<false, false>(values, length, node);
                    return;
                }
#endif
                if (length >= 8 && length <= 32) {
                    if (node == 0) {
                        basecase<false, true>(values, length, node);
                    } else {
                        basecase<false, false>(values, length, node);
                    }
                    return;
                }
#if CRYPTO3_SMALL_PRIME_FFT_RADIX != 4
                constexpr std::size_t radix = CRYPTO3_SMALL_PRIME_FFT_RADIX;
                if (length >= 16 * radix) {
                    const auto stride = length / radix;
                    if (node == 0) {
                        large_radix<radix, false, true>(values, stride, node);
                    } else {
                        large_radix<radix, false, false>(values, stride, node);
                    }
                    for (std::size_t child = 0; child < radix; ++child) {
                        forward_recursive<true>(values + child * stride, stride, radix * node + child);
                    }
                    return;
                }
#endif
                if (length >= 4096) {
                    four_layer_block<false, false>(values, length, node);
                    for (std::size_t child = 0; child < 16; ++child) {
                        forward_recursive<true>(values + child * (length / 16), length / 16, 16 * node + child);
                    }
                    return;
                }
            }
#endif
            if (length == 1) {
                return;
            }
            if (length == 2) {
                const double left = reduce<UseFma>(arithmetic_, values[0]);
                const double right = node == 0 ? reduce<UseFma>(arithmetic_, values[1]) :
                                                 multiply<UseFma>(arithmetic_, values[1], twiddles_[node]);
                values[0] = left + right;
                values[1] = left - right;
                return;
            }
            const std::size_t quarter = length / 4;
            if (node == 0) {
                forward_radix4<UseFma, true>(values, quarter, node, 0, quarter);
            } else {
                forward_radix4<UseFma, false>(values, quarter, node, 0, quarter);
            }
            if (quarter > 1) {
                for (std::size_t child = 0; child < 4; ++child) {
                    forward_recursive<UseFma>(values + child * quarter, quarter, 4 * node + child);
                }
            }
        }

        template<bool UseFma, bool UnitTwiddle, typename Arithmetic, typename Value>
        BOOST_FORCEINLINE static void inverse_butterfly(const Arithmetic &arithmetic, Value &v0, Value &v1, Value &v2,
                                                        Value &v3, Value w2, Value w, Value iw, bool leaf) {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
            // Leaf inputs may be forward residues of magnitude < 3*q. Reduce them before
            // adding four entries; non-leaf inputs already satisfy the inverse bound < 2*q.
            const Value b0 = leaf ? reduce<UseFma>(arithmetic, v0) : v0;
            const Value b1 = leaf ? reduce<UseFma>(arithmetic, v1) : v1;
            const Value b2 = leaf ? reduce<UseFma>(arithmetic, v2) : v2;
            const Value b3 = leaf ? reduce<UseFma>(arithmetic, v3) : v3;
            const Value even_sum = b0 + b1;
            const Value odd_sum = b2 + b3;
            Value left = b0 - b1;
            Value right = b2 - b3;
            if constexpr (UnitTwiddle) {
                left = reduce<UseFma>(arithmetic, left);
                right = multiply<UseFma>(arithmetic, right, iw);
            } else {
                multiply_pair<UseFma>(arithmetic, left, right, w, iw);
            }
            // Four inputs < 2*q give sums < 8*q < 2^53. Times a centered weight this is
            // < 4*q^2, so the a2 remainder is < 1.5*q. Left and right are each < q;
            // their sum is < 2*q and their weighted difference is < q. The invariant closes.
            v0 = reduce<UseFma>(arithmetic, even_sum + odd_sum);
            v1 = left + right;
            v2 = even_sum - odd_sum;
            v3 = left - right;
            if constexpr (UnitTwiddle) {
                v2 = reduce<UseFma>(arithmetic, v2);
                v3 = reduce<UseFma>(arithmetic, v3);
            } else {
                multiply_pair<UseFma>(arithmetic, v2, v3, w2, w2);
            }
        }

        template<bool UseFma, bool UnitTwiddle, bool BoundedInput>
        void inverse_radix4(double *values, std::size_t quarter, std::size_t node, std::size_t first,
                            std::size_t last) const {
            const double w2 = inverse_twiddle(node);
            const double w = inverse_twiddle(2 * node);
            const double iw = inverse_twiddle(2 * node + 1);
            std::size_t i = first;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
            if constexpr (UseFma) {
                if (quarter >= 8) {
                    const small_prime_fft_arithmetic::wide_fma_kernel packed(arithmetic_);
                    const __m512d packed_w2 = _mm512_set1_pd(w2);
                    const __m512d packed_w = _mm512_set1_pd(w);
                    const __m512d packed_iw = _mm512_set1_pd(iw);
#pragma GCC unroll 2
                    for (; last - i >= 8; i += 8) {
                        __m512d v0 = _mm512_loadu_pd(values + i);
                        __m512d v1 = _mm512_loadu_pd(values + i + quarter);
                        __m512d v2 = _mm512_loadu_pd(values + i + 2 * quarter);
                        __m512d v3 = _mm512_loadu_pd(values + i + 3 * quarter);
                        inverse_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, packed_w2, packed_w, packed_iw,
                                                             false);
                        _mm512_storeu_pd(values + i, v0);
                        _mm512_storeu_pd(values + i + quarter, v1);
                        _mm512_storeu_pd(values + i + 2 * quarter, v2);
                        _mm512_storeu_pd(values + i + 3 * quarter, v3);
                    }
                }
            }
#endif
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
                if (quarter >= 4) {
                    const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
                    const __m256d packed_w2 = _mm256_set1_pd(w2);
                    const __m256d packed_w = _mm256_set1_pd(w);
                    const __m256d packed_iw = _mm256_set1_pd(iw);
#pragma GCC unroll 2
                    for (; last - i >= 4; i += 4) {
                        __m256d v0 = _mm256_loadu_pd(values + i);
                        __m256d v1 = _mm256_loadu_pd(values + i + quarter);
                        __m256d v2 = _mm256_loadu_pd(values + i + 2 * quarter);
                        __m256d v3 = _mm256_loadu_pd(values + i + 3 * quarter);
                        inverse_butterfly<true, UnitTwiddle>(packed, v0, v1, v2, v3, packed_w2, packed_w, packed_iw,
                                                             false);
                        _mm256_storeu_pd(values + i, v0);
                        _mm256_storeu_pd(values + i + quarter, v1);
                        _mm256_storeu_pd(values + i + 2 * quarter, v2);
                        _mm256_storeu_pd(values + i + 3 * quarter, v3);
                    }
                }
            }
#endif
            for (; i < last; ++i) {
                inverse_butterfly<UseFma, UnitTwiddle>(arithmetic_, values[i], values[i + quarter],
                                                       values[i + 2 * quarter], values[i + 3 * quarter], w2, w, iw,
                                                       quarter == 1 && !BoundedInput);
            }
        }

        template<bool UseFma, bool BoundedInput = false>
        void inverse_recursive(double *values, std::size_t length, std::size_t node) const {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
            if constexpr (UseFma) {
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
                if (length == 32 || length == 64) {
                    wide_basecase<true, BoundedInput>(values, length, node);
                    return;
                }
#endif
                if (length >= 8 && length <= 32) {
                    if (node == 0) {
                        basecase<true, true, BoundedInput>(values, length, node);
                    } else {
                        basecase<true, false, BoundedInput>(values, length, node);
                    }
                    return;
                }
#if CRYPTO3_SMALL_PRIME_FFT_RADIX != 4
                constexpr std::size_t radix = CRYPTO3_SMALL_PRIME_FFT_RADIX;
                if (length >= 16 * radix) {
                    const auto stride = length / radix;
                    for (std::size_t child = 0; child < radix; ++child) {
                        inverse_recursive<true, BoundedInput>(values + child * stride, stride, radix * node + child);
                    }
                    if (node == 0) {
                        large_radix<radix, true, true>(values, stride, node);
                    } else {
                        large_radix<radix, true, false>(values, stride, node);
                    }
                    return;
                }
#endif
                if (length >= 4096) {
                    for (std::size_t child = 0; child < 16; ++child) {
                        inverse_recursive<true, BoundedInput>(values + child * (length / 16), length / 16,
                                                              16 * node + child);
                    }
                    four_layer_block<true, BoundedInput>(values, length, node);
                    return;
                }
            }
#endif
            if (length == 1) {
                if constexpr (!BoundedInput) {
                    values[0] = reduce<UseFma>(arithmetic_, values[0]);
                }
                return;
            }
            if (length == 2) {
                const double left = BoundedInput ? values[0] : reduce<UseFma>(arithmetic_, values[0]);
                const double right = BoundedInput ? values[1] : reduce<UseFma>(arithmetic_, values[1]);
                values[0] = reduce<UseFma>(arithmetic_, left + right);
                values[1] = node == 0 ? reduce<UseFma>(arithmetic_, left - right) :
                                        multiply<UseFma>(arithmetic_, left - right, inverse_twiddle(node));
                return;
            }
            const std::size_t quarter = length / 4;
            if (quarter > 1) {
                for (std::size_t child = 0; child < 4; ++child) {
                    inverse_recursive<UseFma, BoundedInput>(values + child * quarter, quarter, 4 * node + child);
                }
            }
            if (node == 0) {
                inverse_radix4<UseFma, true, BoundedInput>(values, quarter, node, 0, quarter);
            } else {
                inverse_radix4<UseFma, false, BoundedInput>(values, quarter, node, 0, quarter);
            }
        }

        template<bool UseFma, bool BoundedInput = false>
        void inverse_impl(std::span<double> values, bool normalize) const {
            inverse_recursive<UseFma, BoundedInput>(values.data(), size_, 0);
            if (normalize) {
                const double q = static_cast<double>(modulus());
                std::size_t i = 0;
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED)
                if constexpr (UseFma) {
                    const small_prime_fft_arithmetic::wide_fma_kernel packed(arithmetic_);
                    const __m512d scale = _mm512_set1_pd(size_inverse_);
                    const __m512d packed_q = _mm512_set1_pd(q);
                    for (; values.size() - i >= 8; i += 8) {
                        __m512d value = packed.multiply_fma(_mm512_loadu_pd(values.data() + i), scale);
                        const auto negative = _mm512_cmp_pd_mask(value, _mm512_setzero_pd(), _CMP_LT_OQ);
                        value = _mm512_mask_add_pd(value, negative, value, packed_q);
                        _mm512_storeu_pd(values.data() + i, value);
                    }
                }
#endif
#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
                if constexpr (UseFma) {
                    if (values.size() >= 4) {
                        const small_prime_fft_arithmetic::packed_fma_kernel packed(arithmetic_);
                        const __m256d scale = _mm256_set1_pd(size_inverse_);
                        const __m256d packed_q = _mm256_set1_pd(q);
                        for (; values.size() - i >= 4; i += 4) {
                            __m256d value = packed.multiply_fma(_mm256_loadu_pd(values.data() + i), scale);
                            const __m256d negative = _mm256_cmp_pd(value, _mm256_setzero_pd(), _CMP_LT_OQ);
                            value = _mm256_add_pd(value, _mm256_and_pd(negative, packed_q));
                            _mm256_storeu_pd(values.data() + i, value);
                        }
                    }
                }
#endif
                for (; i < values.size(); ++i) {
                    double &value = values[i];
                    value = multiply<UseFma>(arithmetic_, value, size_inverse_);
                    value = value < 0 ? value + q : value;
                }
            }
        }

        small_prime_fft_arithmetic arithmetic_;
        std::size_t size_;
        std::uint64_t omega_;
        double size_inverse_;
        std::vector<double> twiddles_;
    };

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_DETAIL_SMALL_PRIME_FFT_PLAN_HPP
