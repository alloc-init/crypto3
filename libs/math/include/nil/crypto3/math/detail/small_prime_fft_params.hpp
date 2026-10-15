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

#ifndef CRYPTO3_MATH_DETAIL_SMALL_PRIME_FFT_PARAMS_HPP
#define CRYPTO3_MATH_DETAIL_SMALL_PRIME_FFT_PARAMS_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>

#include <boost/multiprecision/cpp_int.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/field_order.hpp>

namespace nil::crypto3::math::detail {

    struct small_prime_fft_prime {
        std::uint64_t modulus;
        std::uint64_t root;    // Exact order 2^two_adicity; shorter roots are obtained by squaring.
        unsigned two_adicity;
    };

    // These distinct primes q = k*2^s + 1 (k odd) balance residue arithmetic, FFT length and CRT capacity:
    //  - 2^49 < q < 2^50 matches the FMA kernel's certified range and leaves three headroom bits in
    //    binary64's 53-bit significand for exact integer-valued lazy sums with magnitude below 8*q.
    //    Smaller primes contribute fewer CRT bits, generally requiring more primes and transforms
    //    to reconstruct the same full-width BN254 Fp12 product.
    //  - s ranges from 41 to 44, providing roots for every power-of-two transform length up to 2^s.
    //    Entries are ordered by decreasing s, then decreasing q, favoring root capacity first and
    //    more CRT bits among primes with the same root capacity.
    //  - Twelve entries provide spare CRT capacity for larger convolutions; the 32768-point
    //    full-width BN254 Fp12 workload needs eleven. Selection takes only the shortest prefix whose
    //    product M exceeds twice the signed coefficient bound, so unused entries incur no transforms.
    // Other prime/root choices meeting these conditions can also work. Tests certify this catalogue's
    // primality, exact root orders and FMA reciprocal error bounds; arithmetic checks remain separate
    // from these integer reconstruction parameters, which use no floating-point arithmetic.
    inline constexpr std::array<small_prime_fft_prime, 12> small_prime_fft_primes = {{
        {1108307720798209ULL, 194751219211145ULL, 44},
        {659706976665601ULL, 422061946989193ULL, 43},
        {1086317488242689ULL, 813873581740013ULL, 42},
        {910395627798529ULL, 212465101045679ULL, 42},
        {699289395265537ULL, 643876812462936ULL, 42},
        {1022545813831681ULL, 136965991847555ULL, 41},
        {1013749720809473ULL, 549289857321204ULL, 41},
        {868614185943041ULL, 322646709537195ULL, 41},
        {864216139431937ULL, 46940759627233ULL, 41},
        {837827860365313ULL, 128509036413366ULL, 41},
        {815837627809793ULL, 584069832576029ULL, 41},
        {600333348765697ULL, 295748699958414ULL, 41},
    }};

    struct small_prime_fft_crt_params {
        std::size_t prime_count;
        boost::multiprecision::cpp_int modulus;
    };

    /**
     * Select the shortest nonempty catalogue prefix with M > 2*coefficient_bound, where every
     * prime supports transform_size. Then each integer in [-coefficient_bound, coefficient_bound]
     * has a unique centered CRT representative. The bound must cover the integer expression being
     * reconstructed, before reduction modulo the target field characteristic.
     *
     * Return nullopt for a negative bound, a non-power-of-two size (including zero), insufficient
     * root order, or insufficient CRT capacity. Arbitrary-precision setup arithmetic avoids overflow.
     */
    inline std::optional<small_prime_fft_crt_params>
        select_small_prime_fft_crt_params(std::size_t transform_size,
                                          const boost::multiprecision::cpp_int &coefficient_bound) {
        if (!std::has_single_bit(transform_size) || coefficient_bound < 0) {
            return std::nullopt;
        }
        const unsigned depth = std::countr_zero(transform_size);
        const boost::multiprecision::cpp_int twice_bound = 2 * coefficient_bound;
        boost::multiprecision::cpp_int modulus = 1;
        for (std::size_t i = 0; i < small_prime_fft_primes.size(); ++i) {
            if (depth > small_prime_fft_primes[i].two_adicity) {
                return std::nullopt;
            }
            modulus *= small_prime_fft_primes[i].modulus;
            if (modulus > twice_bound) {
                return small_prime_fft_crt_params {i + 1, std::move(modulus)};
            }
        }
        return std::nullopt;
    }

    /**
     * Bound one coordinate of a sum of term_count BN254 Fp12 products, with input coordinates
     * lifted canonically into [0, p-1]. In the tower u^2=-1, v^3=9+u, w^2=v, use the equivalent
     * basis w^j*u^k (0<=j<6, 0<=k<2). At output w^j there are j+1 unwrapped pairs of w powers,
     * with absolute coordinate weight 2 each, and 5-j wrapped pairs, with weight 2*(9+1) each.
     * Thus each output coordinate has weight 2*(j+1)+20*(5-j) <= 102. The triangle inequality
     * gives C = 102*term_count*(p-1)^2, including the signs introduced by u^2=-1.
     *
     * This is a bound on the final integer tower product, not on lazy arithmetic intermediates.
     * Cyclic input folding, digit splitting and other input lifts require their own bounds.
     */
    inline boost::multiprecision::cpp_int bn254_fp12_convolution_bound(std::size_t term_count) {
        const boost::multiprecision::cpp_int largest_coordinate =
            algebra::fields::field_characteristic<algebra::fields::alt_bn128_base_field<254>>() - 1;
        return 102 * boost::multiprecision::cpp_int(term_count) * largest_coordinate * largest_coordinate;
    }

    struct small_prime_fft_params {
        std::size_t left_size;
        std::size_t right_size;
        std::size_t result_size;
        std::size_t transform_size;
        std::size_t spectrum_bytes;    // One spectrum: 12 coordinates * prime_count * transform_size 64-bit cells.
        boost::multiprecision::cpp_int coefficient_bound;
        small_prime_fft_crt_params crt;
    };

    /**
     * Plan a full BN254 Fp12 product, or its first coefficient_count coefficients. Only input prefixes
     * that can affect that output are used, but the transform still contains their COMPLETE linear
     * product: shortening it to the output count would permit cyclic aliasing into low coefficients.
     * Square uses the same plan with equal operand lengths. No transform buffers are allocated here.
     *
     * Return nullopt when either input length or coefficient_count is zero (the caller handles a zero
     * result without a transform), when size arithmetic would overflow, when the catalogue cannot
     * supply the roots/CRT capacity, or when one spectrum exceeds max_spectrum_bytes. This byte limit
     * covers only one spectrum, not all future operation scratch, twiddle tables or prepared operands.
     */
    inline std::optional<small_prime_fft_params>
        make_bn254_small_prime_fft_params(std::size_t left_size, std::size_t right_size,
                                          std::size_t coefficient_count = std::numeric_limits<std::size_t>::max(),
                                          std::size_t max_spectrum_bytes = std::numeric_limits<std::size_t>::max()) {
        left_size = std::min(left_size, coefficient_count);
        right_size = std::min(right_size, coefficient_count);
        if (left_size == 0 || right_size == 0 ||
            left_size > std::numeric_limits<std::size_t>::max() - (right_size - 1)) {
            return std::nullopt;
        }
        const std::size_t product_size = left_size + (right_size - 1);
        // std::bit_ceil requires the result to be representable in size_t.
        if (product_size > (std::size_t(1) << (std::numeric_limits<std::size_t>::digits - 1))) {
            return std::nullopt;
        }
        const std::size_t transform_size = std::bit_ceil(product_size);
        boost::multiprecision::cpp_int bound = bn254_fp12_convolution_bound(std::min(left_size, right_size));
        auto crt = select_small_prime_fft_crt_params(transform_size, bound);
        if (!crt) {
            return std::nullopt;
        }
        const std::size_t bytes_per_index = 12 * crt->prime_count * sizeof(std::uint64_t);
        if (transform_size > max_spectrum_bytes / bytes_per_index) {
            return std::nullopt;
        }
        return small_prime_fft_params {left_size,
                                       right_size,
                                       std::min(product_size, coefficient_count),
                                       transform_size,
                                       transform_size * bytes_per_index,
                                       std::move(bound),
                                       std::move(*crt)};
    }

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_DETAIL_SMALL_PRIME_FFT_PARAMS_HPP
