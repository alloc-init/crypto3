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

#define BOOST_TEST_MODULE small_prime_fft_params_test

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/math/detail/small_prime_fft_params.hpp>

namespace {
    namespace detail = nil::crypto3::math::detail;
    namespace fields = nil::crypto3::algebra::fields;
    using boost::multiprecision::cpp_int;
    using fq = fields::alt_bn128_base_field<254>;
    using fq12 = fields::fp12_2over3over2<fq>;

    // Coordinates are [1,u,v,uv,v^2,uv^2,w,uw,vw,uvw,v^2w,uv^2w].
    std::size_t coordinate_index(unsigned w_power, unsigned u_power) {
        return 6 * (w_power % 2) + 2 * (w_power / 2) + u_power;
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(small_prime_fft_params_test_suite)

BOOST_AUTO_TEST_CASE(catalogue_has_distinct_certified_primes_and_exact_roots) {
    for (std::size_t i = 0; i < detail::small_prime_fft_primes.size(); ++i) {
        const auto &entry = detail::small_prime_fft_primes[i];
        BOOST_TEST_CONTEXT("prime index " << i) {
            BOOST_CHECK_GE(entry.modulus, std::uint64_t(1) << 49);
            BOOST_CHECK_LT(entry.modulus, std::uint64_t(1) << 50);
            BOOST_CHECK_EQUAL(std::countr_zero(entry.modulus - 1), entry.two_adicity);
            BOOST_CHECK_GT(entry.root, 0);
            BOOST_CHECK_LT(entry.root, entry.modulus);
            for (std::size_t j = 0; j < i; ++j) {
                BOOST_CHECK_NE(entry.modulus, detail::small_prime_fft_primes[j].modulus);
            }

            const cpp_int modulus = entry.modulus;
            const cpp_int order = cpp_int(1) << entry.two_adicity;
            const cpp_int root = entry.root;
            const cpp_int half_power = boost::multiprecision::powm(root, cpp_int(order / 2), modulus);
            BOOST_REQUIRE_EQUAL(half_power, modulus - 1);
            BOOST_REQUIRE_EQUAL(boost::multiprecision::powm(root, order, modulus), 1);

            // Every prime divisor r of modulus must give root order exactly 2^s, hence r >= 2^s+1.
            // Since (2^s)^2 > modulus, a composite modulus is impossible. This certifies primality,
            // rather than relying on a probabilistic primality test or assuming the table is prime.
            BOOST_REQUIRE_GT(order * order, modulus);
            for (unsigned depth = 0; depth <= entry.two_adicity; ++depth) {
                const cpp_int short_order = cpp_int(1) << depth;
                const cpp_int short_root = boost::multiprecision::powm(root, cpp_int(order / short_order), modulus);
                BOOST_CHECK_EQUAL(boost::multiprecision::powm(short_root, short_order, modulus), 1);
                if (depth != 0) {
                    BOOST_CHECK_EQUAL(boost::multiprecision::powm(short_root, cpp_int(short_order / 2), modulus),
                                      modulus - 1);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(crt_selection_checks_exact_capacity_boundaries) {
    cpp_int product = 1;
    for (std::size_t i = 0; i < detail::small_prime_fft_primes.size(); ++i) {
        product *= detail::small_prime_fft_primes[i].modulus;
        // The odd product M distinguishes all integers from -(M-1)/2 through +(M-1)/2.
        const auto fits = detail::select_small_prime_fft_crt_params(1, (product - 1) / 2);
        BOOST_REQUIRE(fits);
        BOOST_CHECK_EQUAL(fits->prime_count, i + 1);
        BOOST_CHECK_EQUAL(fits->modulus, product);
        const auto next = detail::select_small_prime_fft_crt_params(1, (product + 1) / 2);
        if (i + 1 < detail::small_prime_fft_primes.size()) {
            BOOST_REQUIRE(next);
            BOOST_CHECK_EQUAL(next->prime_count, i + 2);
        } else {
            BOOST_CHECK(!next);
        }
    }
    BOOST_CHECK(!detail::select_small_prime_fft_crt_params(0, 1));
    BOOST_CHECK(!detail::select_small_prime_fft_crt_params(3, 1));
    BOOST_CHECK(!detail::select_small_prime_fft_crt_params(1, -1));
    const auto zero_bound = detail::select_small_prime_fft_crt_params(1, 0);
    BOOST_REQUIRE(zero_bound);
    BOOST_CHECK_EQUAL(zero_bound->prime_count, 1);
}

BOOST_AUTO_TEST_CASE(tower_bound_matches_all_basis_products_in_crypto3_coordinate_order) {
    const cpp_int modulus = fields::field_characteristic<fq>();
    std::array<unsigned, 12> weights {};
    for (unsigned a = 0; a < 12; ++a) {
        for (unsigned b = 0; b < 12; ++b) {
            fq12::value_type left = fq12::value_type::zero();
            fq12::value_type right = fq12::value_type::zero();
            left.coordinate(a) = fq::value_type::one();
            right.coordinate(b) = fq::value_type::one();
            const auto actual = left * right;

            const unsigned w_power = a / 6 + (a % 6 / 2) * 2 + b / 6 + (b % 6 / 2) * 2;
            const unsigned u_power = a % 2 + b % 2;
            const bool wraps = w_power >= 6;
            std::array<int, 12> expected {};
            // Reduce w^6 to 9+u, then reduce each u^2 to -1 using signed integer arithmetic.
            for (unsigned term = 0; term < (wraps ? 2U : 1U); ++term) {
                const unsigned exponent = u_power + term;
                const int magnitude = wraps && term == 0 ? 9 : 1;
                expected[coordinate_index(w_power % 6, exponent % 2)] +=
                    (exponent / 2 % 2 == 0 ? magnitude : -magnitude);
            }
            for (unsigned c = 0; c < 12; ++c) {
                cpp_int coordinate(actual.coordinate(c).to_integral().backend().to_cpp_int());
                if (coordinate > modulus / 2) {
                    coordinate -= modulus;
                }
                BOOST_CHECK_EQUAL(coordinate, expected[c]);
                weights[c] += expected[c] < 0 ? -expected[c] : expected[c];
            }
        }
    }
    const std::array<unsigned, 12> expected_weights = {102, 102, 66, 66, 30, 30, 84, 84, 48, 48, 12, 12};
    BOOST_CHECK_EQUAL_COLLECTIONS(weights.begin(), weights.end(), expected_weights.begin(), expected_weights.end());
    const cpp_int maximum = std::numeric_limits<std::size_t>::max();
    const cpp_int expected =
        *std::max_element(weights.begin(), weights.end()) * maximum * (modulus - 1) * (modulus - 1);
    BOOST_CHECK_EQUAL(detail::bn254_fp12_convolution_bound(std::numeric_limits<std::size_t>::max()), expected);
    BOOST_CHECK_EQUAL(detail::bn254_fp12_convolution_bound(0), 0);
}

BOOST_AUTO_TEST_CASE(workload_plan_covers_signed_fp12_products) {
    const auto plan = detail::make_bn254_small_prime_fft_params(14939, 14939);
    BOOST_REQUIRE(plan);
    BOOST_CHECK_EQUAL(plan->left_size, 14939);
    BOOST_CHECK_EQUAL(plan->right_size, 14939);
    BOOST_CHECK_EQUAL(plan->result_size, 29877);
    BOOST_CHECK_EQUAL(plan->transform_size, 32768);
    BOOST_CHECK_EQUAL(plan->crt.prime_count, 11);
    BOOST_CHECK_EQUAL(plan->spectrum_bytes, 33 * 1024 * 1024);
    const cpp_int modulus = fields::field_characteristic<fq>();
    const cpp_int twice_bound = cpp_int(204) * 14939 * (modulus - 1) * (modulus - 1);
    BOOST_CHECK_EQUAL(2 * plan->coefficient_bound, twice_bound);
    BOOST_CHECK_EQUAL(boost::multiprecision::msb(twice_bound) + 1, 529);
    BOOST_CHECK_GT(plan->crt.modulus, twice_bound);
    BOOST_CHECK_LE(plan->crt.modulus / detail::small_prime_fft_primes[10].modulus, twice_bound);
}

BOOST_AUTO_TEST_CASE(linear_and_low_shapes_do_not_alias_high_coefficients) {
    for (std::size_t left : {1, 2, 3, 8, 9, 16, 17}) {
        for (std::size_t right : {1, 2, 7, 8, 9, 16, 33}) {
            for (std::size_t count : {1, 2, 7, 8, 9, 16, 17, 64}) {
                const auto plan = detail::make_bn254_small_prime_fft_params(left, right, count);
                BOOST_REQUIRE(plan);
                const std::size_t used_left = std::min(left, count);
                const std::size_t used_right = std::min(right, count);
                const std::size_t full_size = used_left + used_right - 1;
                BOOST_CHECK_EQUAL(plan->left_size, used_left);
                BOOST_CHECK_EQUAL(plan->right_size, used_right);
                BOOST_CHECK_EQUAL(plan->result_size, std::min(count, full_size));
                BOOST_CHECK(std::has_single_bit(plan->transform_size));
                BOOST_CHECK_GE(plan->transform_size, full_size);
                BOOST_CHECK_LT(plan->transform_size / 2, full_size);
                BOOST_CHECK_EQUAL(plan->coefficient_bound,
                                  detail::bn254_fp12_convolution_bound(std::min(used_left, used_right)));
            }
        }
    }
    const auto low = detail::make_bn254_small_prime_fft_params(100, 100, 8);
    BOOST_REQUIRE(low);
    BOOST_CHECK_EQUAL(low->result_size, 8);
    BOOST_CHECK_EQUAL(low->transform_size, 16);
}

BOOST_AUTO_TEST_CASE(empty_overflow_and_storage_limits_are_checked_before_allocation) {
    constexpr auto maximum = std::numeric_limits<std::size_t>::max();
    constexpr auto top_bit = std::size_t(1) << (std::numeric_limits<std::size_t>::digits - 1);
    BOOST_CHECK(!detail::make_bn254_small_prime_fft_params(0, 1));
    BOOST_CHECK(!detail::make_bn254_small_prime_fft_params(1, 0));
    BOOST_CHECK(!detail::make_bn254_small_prime_fft_params(1, 1, 0));
    BOOST_CHECK(!detail::make_bn254_small_prime_fft_params(maximum, 2));
    BOOST_CHECK(!detail::make_bn254_small_prime_fft_params(2, maximum));
    BOOST_CHECK(!detail::make_bn254_small_prime_fft_params(top_bit, 2));
    BOOST_CHECK(!detail::make_bn254_small_prime_fft_params(1, maximum));
    // Trimming before size arithmetic makes this small low product valid, even with huge input lengths.
    const auto trimmed = detail::make_bn254_small_prime_fft_params(maximum, maximum, 1);
    BOOST_REQUIRE(trimmed);
    BOOST_CHECK_EQUAL(trimmed->transform_size, 1);
    const auto plan = detail::make_bn254_small_prime_fft_params(14939, 14939);
    BOOST_REQUIRE(plan);
    BOOST_CHECK(detail::make_bn254_small_prime_fft_params(14939, 14939, maximum, plan->spectrum_bytes));
    BOOST_CHECK(!detail::make_bn254_small_prime_fft_params(14939, 14939, maximum, plan->spectrum_bytes - 1));
    BOOST_CHECK(!detail::make_bn254_small_prime_fft_params(1, 1, 1, 0));
}

BOOST_AUTO_TEST_CASE(root_capacity_and_larger_prime_prefix_are_checked_without_buffers) {
    // These lengths exercise metadata only; the planner never allocates the associated spectra.
    if constexpr (std::numeric_limits<std::size_t>::digits > 44) {
        const std::size_t largest = std::uint64_t(1) << 41;
        const auto plan = detail::make_bn254_small_prime_fft_params(largest / 2, largest / 2);
        BOOST_REQUIRE(plan);
        BOOST_CHECK_EQUAL(plan->transform_size, largest);
        BOOST_CHECK_EQUAL(plan->crt.prime_count, 12);
        BOOST_CHECK(!detail::make_bn254_small_prime_fft_params(largest, 2));
        // At depth 42 the first five primes support the root, but cannot cover a single Fp12 product.
        BOOST_CHECK(!detail::select_small_prime_fft_crt_params(largest * 2, detail::bn254_fp12_convolution_bound(1)));
        BOOST_CHECK(detail::select_small_prime_fft_crt_params(std::uint64_t(1) << 44, 1));
        BOOST_CHECK(!detail::select_small_prime_fft_crt_params(std::uint64_t(1) << 45, 1));
    }
}

BOOST_AUTO_TEST_SUITE_END()
