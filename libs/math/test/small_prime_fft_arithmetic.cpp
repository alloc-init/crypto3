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

#define BOOST_TEST_MODULE small_prime_fft_arithmetic_test

#include <nil/crypto3/math/detail/small_prime_fft_arithmetic.hpp>
#include <nil/crypto3/math/detail/small_prime_fft_params.hpp>

#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/test/unit_test.hpp>

namespace {
    namespace detail = nil::crypto3::math::detail;
    using boost::multiprecision::cpp_int;
    using boost::multiprecision::cpp_rational;
    using arithmetic_type = detail::small_prime_fft_arithmetic;

    struct restore_rounding {
        int saved = std::fegetround();
        ~restore_rounding() {
            std::fesetround(saved);
        }
    };

    cpp_int positive_residue(cpp_int value, std::uint64_t modulus) {
        value %= modulus;
        return value < 0 ? cpp_int(value + modulus) : value;
    }

    void check_residue(double value, const cpp_int &expected, std::uint64_t modulus, double bound) {
        BOOST_REQUIRE(std::isfinite(value));
        BOOST_REQUIRE_EQUAL(value, std::trunc(value));
        BOOST_REQUIRE_LT(std::abs(value), bound);
        BOOST_CHECK_EQUAL(positive_residue(static_cast<std::int64_t>(value), modulus), expected);
    }

    void check_product(const arithmetic_type &arithmetic, std::int64_t left, std::int64_t right) {
        const auto modulus = arithmetic.modulus();
        const cpp_int product = cpp_int(left) * right;
        const cpp_int expected = positive_residue(product, modulus);
        const cpp_int magnitude = boost::multiprecision::abs(product);
        check_residue(arithmetic.multiply_integer(static_cast<double>(left), static_cast<double>(right)), expected,
                      modulus, static_cast<double>(modulus));
        if (arithmetic.fma_supported() && magnitude < 4 * cpp_int(modulus) * modulus) {
            const double bound = magnitude < 2 * cpp_int(modulus) * modulus ? double(modulus) : 1.5 * modulus;
            check_residue(arithmetic.multiply_fma(static_cast<double>(left), static_cast<double>(right)), expected,
                          modulus, bound);
        }
    }

    void check_reduction(const arithmetic_type &arithmetic, std::int64_t value) {
        const auto modulus = arithmetic.modulus();
        const cpp_int expected = positive_residue(value, modulus);
        check_residue(arithmetic.reduce_integer(static_cast<double>(value)), expected, modulus, double(modulus));
        if (arithmetic.fma_supported()) {
            check_residue(arithmetic.reduce_fma(static_cast<double>(value)), expected, modulus, double(modulus));
        }
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(small_prime_fft_arithmetic_test_suite)

BOOST_AUTO_TEST_CASE(reciprocals_have_exact_error_certificates) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    for (const auto &prime : detail::small_prime_fft_primes) {
        BOOST_TEST_CONTEXT("modulus " << prime.modulus) {
            const arithmetic_type arithmetic(prime.modulus);
            BOOST_CHECK_EQUAL(arithmetic.modulus(), prime.modulus);
            // Decode the stored binary64 value and verify its error using independent rational arithmetic.
            const auto bits = std::bit_cast<std::uint64_t>(arithmetic.reciprocal());
            BOOST_CHECK_EQUAL(bits >> 63, 0);
            const int exponent = int((bits >> 52) & 0x7ff) - 1023 - 52;
            BOOST_REQUIRE_EQUAL(exponent, -102);
            const cpp_int significand = (bits & ((std::uint64_t(1) << 52) - 1)) | (std::uint64_t(1) << 52);
            const cpp_rational reciprocal = cpp_rational(significand) / (cpp_int(1) << -exponent);
            const cpp_rational error = boost::multiprecision::abs(reciprocal - cpp_rational(1) / prime.modulus);
            BOOST_CHECK_LT(error, cpp_rational(1) / (cpp_int(1) << 103));

            for (unsigned factor : {2, 4}) {
                const cpp_int product_limit = factor * cpp_int(prime.modulus) * prime.modulus;
                const unsigned product_exponent = boost::multiprecision::msb(product_limit) + 1;
                const cpp_rational rounding_error = factor == 2 ? cpp_rational(5) / 8 : cpp_rational(3) / 4;
                const cpp_rational quotient_error =
                    product_limit * error + (cpp_int(1) << (product_exponent - 54)) * reciprocal + rounding_error;
                BOOST_CHECK_LT(quotient_error, factor == 2 ? cpp_rational(1) : cpp_rational(3) / 2);
            }
#if defined(__FAST_MATH__) || defined(__ASSOCIATIVE_MATH__) || defined(__RECIPROCAL_MATH__) || defined(_M_FP_FAST) || \
    (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0) || FLT_EVAL_METHOD != 0
            BOOST_CHECK(!arithmetic.fma_supported());
#else
            BOOST_REQUIRE(arithmetic.fma_supported());
#endif
        }
    }
}

BOOST_AUTO_TEST_CASE(products_cover_signs_zeros_and_lazy_range_boundaries) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    for (const auto &prime : detail::small_prime_fft_primes) {
        const arithmetic_type arithmetic(prime.modulus);
        const auto q = static_cast<std::int64_t>(prime.modulus);
        std::vector<std::int64_t> values {0,
                                          1,
                                          q / 2,
                                          q / 2 + 1,
                                          q - 1,
                                          q,
                                          q + 1,
                                          3 * q / 2,
                                          3 * q / 2 + 1,
                                          2 * q - 1,
                                          2 * q,
                                          2 * q + 1,
                                          3 * q - 1,
                                          3 * q,
                                          3 * q + 1,
                                          4 * q - 1,
                                          (std::int64_t(1) << 53) - 1};
        const auto positive_count = values.size();
        for (std::size_t i = 1; i < positive_count; ++i) {
            values.push_back(-values[i]);
        }
        for (auto left : values) {
            for (auto right : values) {
                BOOST_TEST_CONTEXT("q=" << q << " left=" << left << " right=" << right) {
                    // The integer path covers all pairs. The FMA path is checked exactly when its product
                    // precondition holds, including (2*q-1)*(2*q+1) = 4*q^2-1 and both signs.
                    check_product(arithmetic, left, right);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(random_products_match_the_integer_oracle) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    std::mt19937_64 random(0xc0ffee);
    for (const auto &prime : detail::small_prime_fft_primes) {
        const arithmetic_type arithmetic(prime.modulus);
        const auto q = static_cast<std::int64_t>(prime.modulus);
        for (unsigned i = 0; i < 2000; ++i) {
            const auto left = static_cast<std::int64_t>(random() % (8 * q - 1)) - (4 * q - 1);
            const auto right = static_cast<std::int64_t>(random() % (2 * q - 1)) - (q - 1);
            BOOST_TEST_CONTEXT("q=" << q << " left=" << left << " right=" << right) {
                check_product(arithmetic, left, right);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(reductions_cover_rounding_boundaries_and_large_integer_representatives) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    std::mt19937_64 random(0x123456);
    constexpr std::int64_t limit = std::int64_t(1) << 53;
    for (const auto &prime : detail::small_prime_fft_primes) {
        const arithmetic_type arithmetic(prime.modulus);
        const auto q = static_cast<std::int64_t>(prime.modulus);
        for (std::int64_t multiple = -15; multiple <= 15; ++multiple) {
            for (std::int64_t offset :
                 {-q / 2 - 1, -q / 2, std::int64_t(-1), std::int64_t(0), std::int64_t(1), q / 2, q / 2 + 1}) {
                const auto value = multiple * q + offset;
                if (value > -limit && value < limit) {
                    BOOST_TEST_CONTEXT("q=" << q << " value=" << value) {
                        check_reduction(arithmetic, value);
                    }
                }
            }
        }
        check_reduction(arithmetic, limit - 1);
        check_reduction(arithmetic, 1 - limit);
        for (unsigned i = 0; i < 2000; ++i) {
            check_reduction(arithmetic, static_cast<std::int64_t>(random() % (2 * limit - 1)) - (limit - 1));
        }
    }
}

BOOST_AUTO_TEST_CASE(integer_fallback_and_setup_work_in_all_rounding_modes) {
    restore_rounding restore;
    for (const auto &prime : detail::small_prime_fft_primes) {
        BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
        const arithmetic_type original(prime.modulus);
        const auto reciprocal_bits = std::bit_cast<std::uint64_t>(original.reciprocal());
        const bool supported_at_nearest = original.fma_supported();
        const auto q = static_cast<std::int64_t>(prime.modulus);
        for (int mode : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
            BOOST_REQUIRE_EQUAL(std::fesetround(mode), 0);
            const arithmetic_type current(prime.modulus);
            BOOST_CHECK_EQUAL(std::fegetround(), mode);
            BOOST_CHECK_EQUAL(std::bit_cast<std::uint64_t>(current.reciprocal()), reciprocal_bits);
            BOOST_CHECK_EQUAL(current.fma_supported(), mode == FE_TONEAREST && supported_at_nearest);
            BOOST_CHECK_EQUAL(original.fma_supported(), mode == FE_TONEAREST && supported_at_nearest);
            check_product(current, -2 * q + 1, 2 * q + 1);
            check_product(original, (std::int64_t(1) << 53) - 1, -q);
            check_reduction(current, -3 * q - 1);
            check_reduction(original, (std::int64_t(1) << 53) - 1);
            BOOST_CHECK_EQUAL(std::fegetround(), mode);
        }
    }
}

BOOST_AUTO_TEST_CASE(invalid_moduli_are_rejected_and_valid_endpoints_work) {
    for (auto modulus :
         {std::uint64_t(0), std::uint64_t(1), (std::uint64_t(1) << 49) - 1, std::uint64_t(1) << 49,
          std::uint64_t(1) << 50, (std::uint64_t(1) << 50) + 1, std::numeric_limits<std::uint64_t>::max()}) {
        BOOST_CHECK_THROW(arithmetic_type {modulus}, std::invalid_argument);
    }
    for (auto modulus : {(std::uint64_t(1) << 49) + 1, (std::uint64_t(1) << 50) - 1}) {
        const arithmetic_type arithmetic(modulus);
        const auto q = static_cast<std::int64_t>(modulus);
        check_product(arithmetic, q - 1, -q + 1);
        check_reduction(arithmetic, -4 * q + 1);
    }
}

#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX2_SELECTED)
template<typename Packed>
void check_packed_products() {
    constexpr auto lanes = Packed::lanes;
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    std::mt19937_64 random(0x41565832);
    for (const auto &prime : detail::small_prime_fft_primes) {
        const arithmetic_type arithmetic(prime.modulus);
        if (!arithmetic.fma_supported()) {
            continue;
        }
        const Packed packed(arithmetic);
        const auto q = static_cast<std::int64_t>(prime.modulus);
        const std::array<std::array<std::int64_t, 2>, 12> boundaries {{{0, 4 * q - 1},
                                                                       {q - 1, 1 - q},
                                                                       {2 * q - 1, 2 * q + 1},
                                                                       {1 - 2 * q, 2 * q + 1},
                                                                       {4 * q - 1, q - 1},
                                                                       {1 - 4 * q, 1 - q},
                                                                       {3 * q - 1, q / 2},
                                                                       {q / 2 + 1, 1 - 3 * q},
                                                                       {(std::int64_t(1) << 53) - 1, 1},
                                                                       {1 - (std::int64_t(1) << 53), -1},
                                                                       {1, -1},
                                                                       {-q, q}}};
        for (std::size_t batch = 0; batch < 1003; ++batch) {
            alignas(64) std::array<double, lanes + 2> left {}, right {}, output;
            output.fill(-19.0);
            for (std::size_t lane = 0; lane < lanes; ++lane) {
                const auto a = batch < 3 ? boundaries[(lanes * batch + lane) % boundaries.size()][0] :
                                           static_cast<std::int64_t>(random() % (8 * q - 1)) - (4 * q - 1);
                const auto b = batch < 3 ? boundaries[(lanes * batch + lane) % boundaries.size()][1] :
                                           static_cast<std::int64_t>(random() % (2 * q - 1)) - (q - 1);
                left[lane + 1] = static_cast<double>(a);
                right[lane + 1] = static_cast<double>(b);
            }
            Packed::store(output.data() + 1,
                          packed.multiply_fma(Packed::load(left.data() + 1), Packed::load(right.data() + 1)));
            std::array<double, lanes> paired, reversed, reversed_left, reversed_right;
            for (std::size_t lane = 0; lane < lanes; ++lane) {
                reversed_left[lane] = left[lanes - lane];
                reversed_right[lane] = right[lanes - lane];
            }
            auto a = Packed::load(left.data() + 1);
            auto b = Packed::load(reversed_left.data());
            packed.multiply_fma_pair(a, b, Packed::load(right.data() + 1), Packed::load(reversed_right.data()));
            Packed::store(paired.data(), a);
            Packed::store(reversed.data(), b);
            for (std::size_t lane = 1; lane <= lanes; ++lane) {
                const cpp_int product =
                    cpp_int(static_cast<std::int64_t>(left[lane])) * static_cast<std::int64_t>(right[lane]);
                const double bound = boost::multiprecision::abs(product) < 2 * cpp_int(q) * q ? double(q) : 1.5 * q;
                BOOST_TEST_CONTEXT("q=" << q << " batch=" << batch << " lane=" << lane) {
                    check_residue(output[lane], positive_residue(product, prime.modulus), prime.modulus, bound);
                    BOOST_CHECK_EQUAL(output[lane], arithmetic.multiply_fma(left[lane], right[lane]));
                    BOOST_CHECK_EQUAL(paired[lane - 1], output[lane]);
                    BOOST_CHECK_EQUAL(reversed[lanes - lane], output[lane]);
                }
            }
            BOOST_CHECK_EQUAL(output.front(), -19.0);
            BOOST_CHECK_EQUAL(output.back(), -19.0);
        }
    }
    BOOST_CHECK_EQUAL(std::fegetround(), FE_TONEAREST);
}

template<typename Packed>
void check_packed_reductions() {
    constexpr auto lanes = Packed::lanes;
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    std::mt19937_64 random(0x41565833);
    constexpr std::int64_t limit = std::int64_t(1) << 53;
    for (const auto &prime : detail::small_prime_fft_primes) {
        const arithmetic_type arithmetic(prime.modulus);
        if (!arithmetic.fma_supported()) {
            continue;
        }
        const Packed packed(arithmetic);
        const auto q = static_cast<std::int64_t>(prime.modulus);
        for (std::size_t batch = 0; batch < 1001; ++batch) {
            alignas(64) std::array<double, lanes + 2> input {}, output;
            output.fill(-23.0);
            for (std::size_t lane = 1; lane <= lanes; ++lane) {
                input[lane] = static_cast<double>(static_cast<std::int64_t>(random() % (2 * limit - 1)) - (limit - 1));
            }
            if (batch == 0) {
                input[1] = double(limit - 1);
                input[2] = double(1 - limit);
                input[3] = double(q / 2);
                input[4] = double(-q / 2 - 1);
            }
            Packed::store(output.data() + 1, packed.reduce_fma(Packed::load(input.data() + 1)));
            for (std::size_t lane = 1; lane <= lanes; ++lane) {
                BOOST_TEST_CONTEXT("q=" << q << " batch=" << batch << " lane=" << lane) {
                    check_residue(output[lane], positive_residue(static_cast<std::int64_t>(input[lane]), prime.modulus),
                                  prime.modulus, double(q));
                    BOOST_CHECK_EQUAL(output[lane], arithmetic.reduce_fma(input[lane]));
                }
            }
            BOOST_CHECK_EQUAL(output.front(), -23.0);
            BOOST_CHECK_EQUAL(output.back(), -23.0);
        }
    }
    BOOST_CHECK_EQUAL(std::fegetround(), FE_TONEAREST);
}

BOOST_AUTO_TEST_CASE(packed_products_match_scalar_and_exact_integer_oracles) {
    check_packed_products<arithmetic_type::packed_fma_kernel>();
}

BOOST_AUTO_TEST_CASE(packed_reductions_preserve_lane_order_and_large_integer_bounds) {
    check_packed_reductions<arithmetic_type::packed_fma_kernel>();
}

#if defined(CRYPTO3_SMALL_PRIME_FFT_AVX512_AVAILABLE)
BOOST_AUTO_TEST_CASE(wide_products_match_scalar_and_exact_integer_oracles) {
    check_packed_products<arithmetic_type::wide_fma_kernel>();
}

BOOST_AUTO_TEST_CASE(wide_reductions_preserve_lane_order_and_large_integer_bounds) {
    check_packed_reductions<arithmetic_type::wide_fma_kernel>();
}
#endif

#endif

BOOST_AUTO_TEST_SUITE_END()
