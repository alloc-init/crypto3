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

#define BOOST_TEST_MODULE small_prime_fft_conversion_test

#include <nil/crypto3/math/polynomial/backends/detail/small_prime_fft/bn254_small_prime_fft_conversion.hpp>

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
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/test/unit_test.hpp>

namespace {
    namespace detail = nil::crypto3::math::detail;
    using boost::multiprecision::cpp_int;
    using fq = nil::crypto3::algebra::fields::alt_bn128_base_field<254>;
    using converter_type = detail::bn254_small_prime_fft_converter;
    using value_type = converter_type::value_type;
    using coordinates_type = std::array<cpp_int, 12>;
    constexpr double sentinel = -12345.0;

    struct restore_rounding {
        int saved = std::fegetround();
        ~restore_rounding() {
            std::fesetround(saved);
        }
    };

    std::vector<value_type> make_coefficients(std::span<const coordinates_type> coordinates) {
        std::vector<value_type> result(coordinates.size(), value_type::zero());
        for (std::size_t i = 0; i < result.size(); ++i) {
            // Assign the tower components directly, independently of coordinate() in the converter.
            for (std::size_t w = 0; w < 2; ++w) {
                for (std::size_t v = 0; v < 3; ++v) {
                    for (std::size_t u = 0; u < 2; ++u) {
                        result[i].data[w].data[v].data[u] =
                            fq::value_type(fq::integral_type(coordinates[i][6 * w + 2 * v + u].str()));
                    }
                }
            }
        }
        return result;
    }

    void check_output(std::span<const double> output, std::size_t prime_count,
                      std::span<const coordinates_type> coordinates, std::size_t stride) {
        const cpp_int characteristic(fq::modulus.backend().to_cpp_int());
        BOOST_REQUIRE_EQUAL(output.size(), prime_count * 12 * stride);
        for (std::size_t prime = 0; prime < prime_count; ++prime) {
            const auto modulus = detail::small_prime_fft_primes[prime].modulus;
            for (std::size_t coordinate = 0; coordinate < 12; ++coordinate) {
                for (std::size_t i = 0; i < stride; ++i) {
                    BOOST_TEST_CONTEXT("prime=" << prime << " coordinate=" << coordinate << " coefficient=" << i) {
                        const double actual = output[(prime * 12 + coordinate) * stride + i];
                        BOOST_REQUIRE(std::isfinite(actual));
                        BOOST_REQUIRE_EQUAL(actual, std::trunc(actual));
                        BOOST_REQUIRE_GE(actual, 0.0);
                        BOOST_REQUIRE_LT(actual, static_cast<double>(modulus));
                        cpp_int expected = 0;
                        if (i < coordinates.size()) {
                            expected = (coordinates[i][coordinate] % characteristic) % modulus;
                        }
                        BOOST_CHECK_EQUAL(static_cast<std::uint64_t>(actual), expected);
                    }
                }
            }
        }
    }

    std::vector<double> check_conversion(const converter_type &converter, std::span<const value_type> coefficients,
                                         std::span<const coordinates_type> coordinates, std::size_t stride) {
        const std::size_t cells = converter.prime_count() * 12 * stride;
        std::vector<double> storage(cells + 2, sentinel);
        auto output = std::span<double>(storage).subspan(1, cells);
        converter.convert(output, coefficients, stride);
        BOOST_CHECK_EQUAL(storage.front(), sentinel);
        BOOST_CHECK_EQUAL(storage.back(), sentinel);
        check_output(output, converter.prime_count(), coordinates, stride);
        return storage;
    }

    std::vector<coordinates_type> random_coordinates(std::size_t count) {
        std::mt19937_64 random(0x51a11);
        std::vector<coordinates_type> result(count);
        for (auto &coefficient : result) {
            for (auto &coordinate : coefficient) {
                for (unsigned limb = 0; limb < 4; ++limb) {
                    coordinate <<= 64;
                    coordinate += random();
                }
                coordinate &= (cpp_int(1) << 254) - 1;
            }
        }
        return result;
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(small_prime_fft_conversion_test_suite)

BOOST_AUTO_TEST_CASE(basis_and_coordinate_layout_are_preserved_for_every_catalogue_prefix) {
    std::vector<coordinates_type> coordinates(14);
    for (std::size_t c = 0; c < 12; ++c) {
        coordinates[c][c] = 1;
        coordinates[12][c] = 100 + c;
    }
    const auto coefficients = make_coefficients(coordinates);
    for (std::size_t count = 1; count <= detail::small_prime_fft_primes.size(); ++count) {
        const converter_type converter(count);
        BOOST_CHECK_EQUAL(converter.prime_count(), count);
        check_conversion(converter, coefficients, coordinates, 19);
    }
}

BOOST_AUTO_TEST_CASE(canonical_values_and_chunk_boundaries_match_integer_oracle) {
    const cpp_int characteristic(fq::modulus.backend().to_cpp_int());
    std::vector<cpp_int> values = {
        0, 1, characteristic - 1, characteristic, characteristic + 1, (cpp_int(1) << 254) - 1};
    // Every bit position and its neighbors cover chunk boundaries, carries, and sparse coordinates.
    for (unsigned bit = 0; bit < 254; ++bit) {
        const cpp_int power = cpp_int(1) << bit;
        values.push_back(power - 1);
        values.push_back(power);
        values.push_back(power + 1);
    }
    // Also cross limb boundaries in the stored Montgomery value a*R mod p.
    // Multiplying by R^-1 constructs canonical inputs independently of the decoder.
    using base_backend_type = fq::integral_type::backend_type;
    const cpp_int radix = cpp_int(1) << (base_backend_type::limb_bits * base_backend_type::internal_limb_count);
    const cpp_int inverse_radix = boost::multiprecision::powm(radix, cpp_int(characteristic - 2), characteristic);
    values.push_back(inverse_radix);
    values.push_back(((characteristic - 1) * inverse_radix) % characteristic);
    for (unsigned bit = base_backend_type::limb_bits; bit < 254; bit += base_backend_type::limb_bits) {
        for (int delta : {-1, 0, 1}) {
            values.push_back((((cpp_int(1) << bit) + delta) * inverse_radix) % characteristic);
        }
    }
    for (const auto &prime : detail::small_prime_fft_primes) {
        for (int delta : {-1, 0, 1}) {
            values.push_back(cpp_int(prime.modulus) + delta);
            values.push_back((characteristic / prime.modulus) * prime.modulus + delta);
        }
    }
    std::vector<coordinates_type> coordinates((values.size() + 11) / 12);
    for (std::size_t i = 0; i < values.size(); ++i) {
        coordinates[i / 12][i % 12] = values[i];
    }
    auto coefficients = make_coefficients(coordinates);
    // Exercise a field negation as well as construction from canonical and noncanonical integers.
    coefficients[0].data[0].data[1].data[0] = -fq::value_type::one();
    BOOST_REQUIRE_EQUAL(coordinates[0][2], characteristic - 1);
    check_conversion(converter_type(12), coefficients, coordinates, std::bit_ceil(coordinates.size()));
}

BOOST_AUTO_TEST_CASE(random_full_width_coordinates_match_integer_oracle) {
    const auto coordinates = random_coordinates(129);
    const auto coefficients = make_coefficients(coordinates);
    const converter_type converter(12);
    check_conversion(converter, coefficients, coordinates, 129);
    check_conversion(converter, coefficients, coordinates, 256);
}

BOOST_AUTO_TEST_CASE(rounding_modes_preserve_exact_conversion_with_reused_context) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    const auto coordinates = random_coordinates(33);
    const auto coefficients = make_coefficients(coordinates);
    const converter_type original(12);
    const auto baseline = check_conversion(original, coefficients, coordinates, 64);
    for (int mode : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
        BOOST_REQUIRE_EQUAL(std::fesetround(mode), 0);
        const converter_type current(12);
        for (const auto *converter : {&original, &current}) {
            const auto actual = check_conversion(*converter, coefficients, coordinates, 64);
            BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), baseline.begin(), baseline.end());
            BOOST_CHECK_EQUAL(std::fegetround(), mode);
        }
    }
}

BOOST_AUTO_TEST_CASE(reused_buffers_are_zero_padded_and_empty_inputs_work) {
    const auto coordinates = random_coordinates(33);
    const auto coefficients = make_coefficients(coordinates);
    const converter_type converter(3);
    constexpr std::size_t stride = 64;
    std::vector<double> output(3 * 12 * stride, sentinel);
    for (std::size_t count : {std::size_t(33), std::size_t(1), std::size_t(0)}) {
        converter.convert(output, std::span<const value_type>(coefficients).first(count), stride);
        check_output(output, 3, std::span<const coordinates_type>(coordinates).first(count), stride);
    }
    // A zero coefficient also overwrites every active cell of a previous nonzero conversion.
    converter.convert(output, coefficients, stride);
    const std::array<coordinates_type, 1> zero {};
    converter.convert(output, make_coefficients(zero), stride);
    check_output(output, 3, zero, stride);
    converter.convert({}, {}, 0);
}

BOOST_AUTO_TEST_CASE(invalid_dimensions_and_prime_counts_are_rejected_before_writes) {
    for (std::size_t count : {std::size_t(0), std::size_t(13), std::numeric_limits<std::size_t>::max()}) {
        BOOST_CHECK_THROW(converter_type {count}, std::invalid_argument);
    }
    const converter_type converter(1);
    const std::array<value_type, 2> coefficients = {value_type::one(), value_type::one()};
    std::vector<double> output(24, sentinel);
    BOOST_CHECK_THROW(converter.convert(std::span<double>(output).first(23), coefficients, 2), std::invalid_argument);
    BOOST_CHECK_THROW(converter.convert(output, coefficients, 1), std::invalid_argument);
    BOOST_CHECK_THROW(converter.convert(output, {}, 1), std::invalid_argument);
    BOOST_CHECK_THROW(converter.convert(std::span<double>(output).first(12), coefficients, 1), std::invalid_argument);
    BOOST_CHECK_THROW(converter.convert(output, {}, std::numeric_limits<std::size_t>::max()), std::invalid_argument);
    BOOST_CHECK_THROW(converter.convert(output, {}, 0), std::invalid_argument);
    BOOST_CHECK_THROW(converter.convert({}, coefficients, 2), std::invalid_argument);
    BOOST_CHECK(std::all_of(output.begin(), output.end(), [](double value) { return value == sentinel; }));
}

BOOST_AUTO_TEST_CASE(conversion_uses_existing_plan_sizes_and_input_prefixes) {
    const auto plan = detail::make_bn254_small_prime_fft_params(35, 17, 19);
    BOOST_REQUIRE(plan);
    const converter_type converter(plan->crt.prime_count);
    const auto coordinates = random_coordinates(35);
    const auto coefficients = make_coefficients(coordinates);
    std::vector<double> output(plan->spectrum_bytes / sizeof(double), sentinel);
    for (std::size_t count : {plan->left_size, plan->right_size}) {
        converter.convert(output, std::span<const value_type>(coefficients).first(count), plan->transform_size);
        check_output(output, plan->crt.prime_count, std::span<const coordinates_type>(coordinates).first(count),
                     plan->transform_size);
    }
}

BOOST_AUTO_TEST_SUITE_END()
