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

#define BOOST_TEST_MODULE small_prime_fft_transform_test

#include <nil/crypto3/math/detail/bn254_small_prime_fft.hpp>
#include <nil/crypto3/math/detail/field_utils.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/test/unit_test.hpp>

namespace {
    namespace detail = nil::crypto3::math::detail;
    using boost::multiprecision::cpp_int;
    using plan_type = detail::small_prime_fft_plan;

    struct restore_rounding {
        int saved = std::fegetround();
        ~restore_rounding() {
            std::fesetround(saved);
        }
    };

    cpp_int canonical(const cpp_int &value, const cpp_int &modulus) {
        cpp_int result = value % modulus;
        return result < 0 ? result + modulus : result;
    }

    cpp_int as_integer(double value) {
        BOOST_REQUIRE(std::isfinite(value));
        BOOST_REQUIRE(std::abs(value) < 0x1p53);
        BOOST_REQUIRE_EQUAL(std::trunc(value), value);
        return static_cast<std::int64_t>(value);
    }

    cpp_int root_for(std::size_t prime, std::size_t size) {
        const auto &entry = detail::small_prime_fft_primes[prime];
        const cpp_int exponent = cpp_int(1) << (entry.two_adicity - std::countr_zero(size));
        return boost::multiprecision::powm(cpp_int(entry.root), exponent, cpp_int(entry.modulus));
    }

    // Direct O(N^2) integer DFT, independent of the butterfly schedule and floating-point kernels.
    std::vector<cpp_int> dft(std::span<const double> input, std::size_t prime, bool inverse = false) {
        const cpp_int q = detail::small_prime_fft_primes[prime].modulus;
        cpp_int root = root_for(prime, input.size());
        cpp_int scale = 1;
        if (inverse) {
            root = boost::multiprecision::powm(root, q - 2, q);
            scale = boost::multiprecision::powm(cpp_int(input.size()), q - 2, q);
        }
        std::vector<cpp_int> result(input.size());
        cpp_int step = 1;
        for (std::size_t k = 0; k < input.size(); ++k) {
            cpp_int sum = 0;
            cpp_int power = 1;
            for (double value : input) {
                sum += static_cast<std::int64_t>(value) * power;
                power = power * step % q;
            }
            result[k] = canonical(sum * scale, q);
            step = step * root % q;
        }
        return result;
    }

    std::vector<double> lazy_input(std::size_t size, std::uint64_t modulus, std::mt19937_64 &random) {
        const auto q = static_cast<std::int64_t>(modulus);
        const std::array<std::int64_t, 13> edges {3 * q - 1, 1 - 3 * q, 0,  1,      -1,    q - 1, q,
                                                  q + 1,     -q - 1,    -q, -q + 1, 2 * q, -2 * q};
        std::vector<double> result(size);
        for (std::size_t i = 0; i < size; ++i) {
            const auto value = i % 2 == 0 ? edges[(i / 2) % edges.size()] :
                                            static_cast<std::int64_t>(random() % (6 * q - 1)) - (3 * q - 1);
            result[i] = static_cast<double>(value);
        }
        return result;
    }

    std::vector<double> bit_reversed(std::span<const double> input) {
        std::vector<double> result(input.size());
        const unsigned depth = std::countr_zero(input.size());
        for (std::size_t i = 0; i < input.size(); ++i) {
            result[detail::bitreverse(i, depth)] = input[i];
        }
        return result;
    }

    void check_residue(double actual, const cpp_int &expected, std::uint64_t modulus, unsigned bound,
                       bool normalized = false) {
        const cpp_int value = as_integer(actual);
        BOOST_CHECK(std::abs(actual) < bound * static_cast<double>(modulus));
        if (normalized) {
            BOOST_CHECK(actual >= 0);
        }
        BOOST_CHECK_EQUAL(canonical(value, cpp_int(modulus)), expected);
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(small_prime_fft_transform_tests)

BOOST_AUTO_TEST_CASE(forward_matches_integer_dft_for_every_prime) {
    std::mt19937_64 random(0x46465401);
    for (std::size_t prime = 0; prime < detail::small_prime_fft_primes.size(); ++prime) {
        for (std::size_t size = 1; size <= 64; size *= 2) {
            BOOST_TEST_CONTEXT("prime=" << prime << " N=" << size) {
                const plan_type plan(prime, size);
                const auto q = detail::small_prime_fft_primes[prime].modulus;
                BOOST_CHECK_EQUAL(plan.size(), size);
                BOOST_CHECK_EQUAL(plan.modulus(), q);
                BOOST_CHECK_EQUAL(plan.omega(), root_for(prime, size));
                BOOST_CHECK(std::abs(plan.size_inverse()) < static_cast<double>(q) / 2);
                BOOST_CHECK_EQUAL(canonical(as_integer(plan.size_inverse()) * size, cpp_int(q)), 1);

                auto values = lazy_input(size, q, random);
                const auto expected = dft(values, prime);
                plan.fft(values);
                for (std::size_t k = 0; k < size; ++k) {
                    check_residue(values[detail::bitreverse(k, std::countr_zero(size))], expected[k], q, 3);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(inverse_matches_independent_spectra) {
    std::mt19937_64 random(0x46465402);
    for (std::size_t prime = 0; prime < detail::small_prime_fft_primes.size(); ++prime) {
        for (std::size_t size = 1; size <= 64; size *= 2) {
            BOOST_TEST_CONTEXT("prime=" << prime << " N=" << size) {
                const plan_type plan(prime, size);
                const auto q = plan.modulus();
                const auto spectrum = lazy_input(size, q, random);
                const auto expected = dft(spectrum, prime, true);
                auto normalized = bit_reversed(spectrum);
                auto unnormalized = normalized;
                plan.inverse_fft(normalized);
                plan.inverse_fft(unnormalized, false);
                for (std::size_t i = 0; i < size; ++i) {
                    check_residue(normalized[i], expected[i], q, 1, true);
                    check_residue(unnormalized[i], canonical(expected[i] * size, cpp_int(q)), q, 2);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(basis_vectors_pin_spectrum_order) {
    constexpr std::size_t size = 16;
    for (std::size_t prime = 0; prime < detail::small_prime_fft_primes.size(); ++prime) {
        const plan_type plan(prime, size);
        const cpp_int q = plan.modulus();
        const cpp_int root = root_for(prime, size);
        for (std::size_t coefficient = 0; coefficient < size; ++coefficient) {
            std::vector<double> values(size, 0.0);
            values[coefficient] = 1.0;
            plan.fft(values);
            for (std::size_t k = 0; k < size; ++k) {
                const cpp_int expected = boost::multiprecision::powm(root, cpp_int(coefficient * k), q);
                check_residue(values[detail::bitreverse(k, 4)], expected, plan.modulus(), 3);
            }
            plan.inverse_fft(values);
            for (std::size_t i = 0; i < size; ++i) {
                BOOST_CHECK_EQUAL(values[i], i == coefficient ? 1.0 : 0.0);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(constant_and_alternating_extremes) {
    for (std::size_t prime = 0; prime < detail::small_prime_fft_primes.size(); ++prime) {
        const auto q = detail::small_prime_fft_primes[prime].modulus;
        const auto signed_q = static_cast<std::int64_t>(q);
        for (std::size_t size = 1; size <= 128; size *= 2) {
            const plan_type plan(prime, size);
            for (const std::int64_t value :
                 {std::int64_t(0), std::int64_t(1), signed_q, -signed_q, 3 * signed_q - 1, 1 - 3 * signed_q}) {
                for (bool alternating : {false, true}) {
                    std::vector<double> values(size);
                    for (std::size_t i = 0; i < size; ++i) {
                        values[i] = static_cast<double>(alternating && i % 2 ? -value : value);
                    }
                    plan.fft(values);
                    for (std::size_t k = 0; k < size; ++k) {
                        const std::size_t nonzero_index = alternating ? size / 2 : 0;
                        const cpp_int expected = k == nonzero_index ? canonical(cpp_int(value) * size, cpp_int(q)) : 0;
                        check_residue(values[detail::bitreverse(k, std::countr_zero(size))], expected, q, 3);
                    }
                    plan.inverse_fft(values);
                    for (std::size_t i = 0; i < size; ++i) {
                        const cpp_int expected = canonical(cpp_int(alternating && i % 2 ? -value : value), cpp_int(q));
                        check_residue(values[i], expected, q, 1, true);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(round_trips_at_larger_even_and_odd_depths) {
    std::mt19937_64 random(0x46465403);
    for (std::size_t prime = 0; prime < detail::small_prime_fft_primes.size(); ++prime) {
        for (const std::size_t size : {128, 256, 512, 1024, 16384, 32768}) {
            BOOST_TEST_CONTEXT("prime=" << prime << " N=" << size) {
                const plan_type plan(prime, size);
                const auto input = lazy_input(size, plan.modulus(), random);
                auto values = input;
                plan.fft(values);
                for (double value : values) {
                    as_integer(value);
                    BOOST_CHECK(std::abs(value) < 3 * static_cast<double>(plan.modulus()));
                }
                plan.inverse_fft(values);
                for (std::size_t i = 0; i < size; ++i) {
                    check_residue(values[i], canonical(as_integer(input[i]), cpp_int(plan.modulus())), plan.modulus(),
                                  1, true);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(convolution_matches_schoolbook_with_both_scaling_paths) {
    std::mt19937_64 random(0x46465404);
    const std::array<std::pair<std::size_t, std::size_t>, 7> shapes {
        {{1, 1}, {1, 17}, {2, 3}, {7, 9}, {17, 33}, {32, 32}, {65, 31}}};
    for (std::size_t prime = 0; prime < detail::small_prime_fft_primes.size(); ++prime) {
        for (const auto [left_size, right_size] : shapes) {
            BOOST_TEST_CONTEXT("prime=" << prime << " shape=" << left_size << "x" << right_size) {
                const auto q = detail::small_prime_fft_primes[prime].modulus;
                const auto left = lazy_input(left_size, q, random);
                const auto right = lazy_input(right_size, q, random);
                const std::size_t size = std::bit_ceil(left_size + right_size - 1);
                std::vector<cpp_int> expected(size);
                for (std::size_t i = 0; i < left_size; ++i) {
                    for (std::size_t j = 0; j < right_size; ++j) {
                        expected[i + j] += as_integer(left[i]) * as_integer(right[j]);
                    }
                }
                for (auto &value : expected) {
                    value = canonical(value, cpp_int(q));
                }

                const plan_type plan(prime, size);
                const detail::small_prime_fft_arithmetic arithmetic(q);
                std::vector<double> a(size, 0.0), b(size, 0.0);
                std::copy(left.begin(), left.end(), a.begin());
                std::copy(right.begin(), right.end(), b.begin());
                plan.fft(a);
                plan.fft(b);
                for (std::size_t i = 0; i < size; ++i) {
                    a[i] = arithmetic.multiply_integer(a[i], b[i]);
                    b[i] = arithmetic.multiply_integer(a[i], plan.size_inverse());
                }
                auto unnormalized = a;
                plan.inverse_fft(a);
                plan.inverse_fft(b, false);
                plan.inverse_fft(unnormalized, false);
                for (std::size_t i = 0; i < size; ++i) {
                    check_residue(a[i], expected[i], q, 1, true);
                    check_residue(b[i], expected[i], q, 2);
                    check_residue(unnormalized[i], canonical(expected[i] * size, cpp_int(q)), q, 2);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(rounding_modes_preserve_exactness_and_environment) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    std::mt19937_64 random(0x46465405);
    for (std::size_t prime = 0; prime < detail::small_prime_fft_primes.size(); ++prime) {
        for (const std::size_t size : {32, 64}) {
            const plan_type original(prime, size);
            const auto q = original.modulus();
            const auto input = lazy_input(size, q, random);
            const auto expected = dft(input, prime);
            for (const int mode : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
                BOOST_TEST_CONTEXT("prime=" << prime << " N=" << size << " rounding=" << mode) {
                    BOOST_REQUIRE_EQUAL(std::fesetround(mode), 0);
                    const plan_type constructed_in_mode(prime, size);
                    for (const plan_type *plan : {&original, &constructed_in_mode}) {
                        auto values = input;
                        plan->fft(values);
                        for (std::size_t k = 0; k < size; ++k) {
                            check_residue(values[detail::bitreverse(k, std::countr_zero(size))], expected[k], q, 3);
                        }
                        plan->inverse_fft(values);
                        for (std::size_t i = 0; i < size; ++i) {
                            check_residue(values[i], canonical(as_integer(input[i]), cpp_int(q)), q, 1, true);
                        }
                        BOOST_CHECK_EQUAL(std::fegetround(), mode);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(plan_reuse_and_in_place_buffer_boundaries) {
    constexpr std::size_t size = 32;
    const plan_type original(0, size);
    auto copy = original;
    const plan_type moved = std::move(copy);
    std::mt19937_64 random(0x46465406);
    for (const plan_type *plan : {&original, &moved}) {
        for (unsigned trial = 0; trial < 3; ++trial) {
            const auto input = lazy_input(size, plan->modulus(), random);
            std::vector<double> storage(size + 2, -17.0);
            auto values = std::span<double>(storage).subspan(1, size);
            std::copy(input.begin(), input.end(), values.begin());
            const double *address = values.data();
            plan->fft(values);
            plan->inverse_fft(values);
            BOOST_CHECK(values.data() == address);
            BOOST_CHECK_EQUAL(storage.front(), -17.0);
            BOOST_CHECK_EQUAL(storage.back(), -17.0);
            for (std::size_t i = 0; i < size; ++i) {
                check_residue(values[i], canonical(as_integer(input[i]), cpp_int(plan->modulus())), plan->modulus(), 1,
                              true);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(invalid_parameters_and_dimensions_throw_before_writes) {
    BOOST_CHECK_THROW(plan_type(detail::small_prime_fft_primes.size(), 1), std::invalid_argument);
    BOOST_CHECK_THROW(plan_type(std::numeric_limits<std::size_t>::max(), 1), std::invalid_argument);
    for (const std::size_t size :
         {std::size_t(0), std::size_t(3), std::size_t(6), std::numeric_limits<std::size_t>::max()}) {
        BOOST_CHECK_THROW(plan_type(0, size), std::invalid_argument);
    }
    for (std::size_t prime = 0; prime < detail::small_prime_fft_primes.size(); ++prime) {
        const unsigned depth = detail::small_prime_fft_primes[prime].two_adicity + 1;
        if (depth < std::numeric_limits<std::size_t>::digits) {
            BOOST_CHECK_THROW(plan_type(prime, std::size_t(1) << depth), std::invalid_argument);
        }
    }
    const plan_type plan(0, 8);
    for (const std::size_t size : {0, 7, 9}) {
        std::vector<double> values(size, 42.0);
        const auto before = values;
        BOOST_CHECK_THROW(plan.fft(values), std::invalid_argument);
        BOOST_CHECK_THROW(plan.inverse_fft(values), std::invalid_argument);
        BOOST_CHECK_THROW(plan.inverse_fft(values, false), std::invalid_argument);
        BOOST_CHECK_EQUAL_COLLECTIONS(values.begin(), values.end(), before.begin(), before.end());
    }
}

BOOST_AUTO_TEST_SUITE_END()
