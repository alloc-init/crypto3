// Copyright (c) 2026. Distributed under the MIT License.

#include <nil/crypto3/math/polynomial/backends/detail/small_prime_fft/bn254_small_prime_fft_crt.hpp>

#ifdef CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED
constexpr bool fft_avx2 = true;
#else
constexpr bool fft_avx2 = false;
#endif
#ifdef CRYPTO3_SMALL_PRIME_FFT_AVX512_SELECTED
constexpr bool fft_avx512 = true;
#else
constexpr bool fft_avx512 = false;
#endif
#ifdef CRYPTO3_SMALL_PRIME_CRT_IFMA_SELECTED
constexpr bool crt_ifma = true;
#else
constexpr bool crt_ifma = false;
#endif

static_assert(fft_avx2 == bool(EXPECT_FFT_AVX2), "AVX2/FMA FFT selection does not match configuration");
static_assert(fft_avx512 == bool(EXPECT_FFT_AVX512), "AVX512F/DQ FFT selection does not match configuration");
static_assert(crt_ifma == bool(EXPECT_CRT_IFMA), "IFMA CRT selection does not match configuration");

int main() {
    return 0;
}
