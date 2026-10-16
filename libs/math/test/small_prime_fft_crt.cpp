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

#define BOOST_TEST_MODULE small_prime_fft_crt_test

#include <nil/crypto3/math/detail/bn254_small_prime_fft.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
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
    using reconstructor_type = detail::bn254_small_prime_fft_reconstructor;
    using value_type = reconstructor_type::value_type;
    using coordinates_type = std::array<cpp_int, 12>;

    struct restore_rounding {
        int saved = std::fegetround();
        ~restore_rounding() {
            std::fesetround(saved);
        }
    };

    cpp_int canonical(const cpp_int &value, const cpp_int &modulus) {
        cpp_int result = value % modulus;
        if (result < 0) {
            result += modulus;
        }
        return result;
    }

    detail::small_prime_fft_crt_params make_params(std::size_t prime_count) {
        cpp_int product = 1;
        for (std::size_t i = 0; i < prime_count; ++i) {
            product *= detail::small_prime_fft_primes[i].modulus;
        }
        return {prime_count, product};
    }

    std::vector<double> make_residues(std::span<const coordinates_type> coordinates, std::size_t prime_count,
                                      std::size_t stride, bool lazy = false) {
        // Padding is intentionally not a valid residue: reconstruction must read only the requested prefix.
        std::vector<double> result(prime_count * 12 * stride, std::numeric_limits<double>::quiet_NaN());
        for (std::size_t prime = 0; prime < prime_count; ++prime) {
            const cpp_int modulus = detail::small_prime_fft_primes[prime].modulus;
            for (std::size_t coordinate = 0; coordinate < 12; ++coordinate) {
                for (std::size_t i = 0; i < coordinates.size(); ++i) {
                    cpp_int residue = canonical(coordinates[i][coordinate], modulus);
                    if (lazy) {
                        // Reach large positive and negative representatives while staying below 2^53.
                        residue += ((prime + coordinate + i) % 2 == 0 ? 7 : -8) * modulus;
                    }
                    result[(prime * 12 + coordinate) * stride + i] =
                        static_cast<double>(residue.convert_to<std::int64_t>());
                }
            }
        }
        return result;
    }

    void check_output(std::span<const value_type> output, std::span<const coordinates_type> coordinates) {
        BOOST_REQUIRE_EQUAL(output.size(), coordinates.size());
        const cpp_int characteristic(fq::modulus.backend().to_cpp_int());
        for (std::size_t i = 0; i < output.size(); ++i) {
            for (std::size_t w = 0; w < 2; ++w) {
                for (std::size_t v = 0; v < 3; ++v) {
                    for (std::size_t u = 0; u < 2; ++u) {
                        const std::size_t coordinate = 6 * w + 2 * v + u;
                        BOOST_TEST_CONTEXT("coefficient=" << i << " coordinate=" << coordinate) {
                            // Inspect the tower directly, independently of coordinate() in reconstruction.
                            const cpp_int actual(
                                output[i].data[w].data[v].data[u].to_integral().backend().to_cpp_int());
                            BOOST_CHECK_EQUAL(actual, canonical(coordinates[i][coordinate], characteristic));
                        }
                    }
                }
            }
        }
    }

    void check_reconstruction(const reconstructor_type &reconstructor, std::span<const coordinates_type> coordinates,
                              std::size_t stride, bool lazy = false) {
        const auto residues = make_residues(coordinates, reconstructor.prime_count(), stride, lazy);
        const auto sentinel = value_type::one();
        std::vector<value_type> storage(coordinates.size() + 2, sentinel);
        auto output = std::span<value_type>(storage).subspan(1, coordinates.size());
        reconstructor.reconstruct(output, residues, stride);
        BOOST_CHECK(storage.front() == sentinel);
        BOOST_CHECK(storage.back() == sentinel);
        check_output(output, coordinates);
    }

    std::vector<coordinates_type> random_coordinates(std::size_t count, const cpp_int &bound) {
        std::mt19937_64 random(0xc47);
        std::vector<coordinates_type> result(count);
        for (auto &coefficient : result) {
            for (auto &coordinate : coefficient) {
                for (unsigned limb = 0; limb < 10; ++limb) {
                    coordinate <<= 64;
                    coordinate += random();
                }
                coordinate = coordinate % (2 * bound + 1) - bound;
            }
        }
        return result;
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(small_prime_fft_crt_test_suite)

BOOST_AUTO_TEST_CASE(every_catalogue_prefix_recovers_signed_boundaries_and_carries) {
    const cpp_int characteristic(fq::modulus.backend().to_cpp_int());
    for (std::size_t count = 1; count <= detail::small_prime_fft_primes.size(); ++count) {
        BOOST_TEST_CONTEXT("prime count=" << count) {
            const auto params = make_params(count);
            const cpp_int bound = params.modulus / 2;
            const reconstructor_type reconstructor(params, bound);
            BOOST_CHECK_EQUAL(reconstructor.prime_count(), count);
            std::vector<cpp_int> values = {0, 1, -1, bound, bound - 1, -bound, 1 - bound};
            for (cpp_int value : {cpp_int(characteristic - 1), characteristic, cpp_int(characteristic + 1)}) {
                if (value <= bound) {
                    values.push_back(value);
                    values.push_back(-value);
                }
            }
            for (unsigned bit = 1; bit <= boost::multiprecision::msb(bound); ++bit) {
                if (bit % 64 != 0 && bit % 50 != 0) {
                    continue;
                }
                for (int delta : {-1, 0, 1}) {
                    const cpp_int value = (cpp_int(1) << bit) + delta;
                    if (value <= bound) {
                        values.push_back(value);
                        values.push_back(-value);
                    }
                }
            }
            // All CRT digits q_i-1 produce the largest cofactor sum, just below count*M.
            // The associated signed integer is -sum(M/q_i), which exercises every correction bit.
            cpp_int largest_sum_value = 0;
            for (std::size_t prime = 0; prime < count; ++prime) {
                largest_sum_value -= params.modulus / detail::small_prime_fft_primes[prime].modulus;
            }
            values.push_back(largest_sum_value);
            std::vector<coordinates_type> coordinates(values.size());
            for (std::size_t i = 0; i < coordinates.size(); ++i) {
                for (std::size_t c = 0; c < 12; ++c) {
                    coordinates[i][c] = values[(i + c) % values.size()];
                }
            }
            check_reconstruction(reconstructor, coordinates, coordinates.size() + 3);
            check_reconstruction(reconstructor, coordinates, coordinates.size(), true);
        }
    }
}

BOOST_AUTO_TEST_CASE(random_signed_integers_cover_all_prime_counts_and_partial_blocks) {
    for (std::size_t count = 1; count <= detail::small_prime_fft_primes.size(); ++count) {
        const auto params = make_params(count);
        const cpp_int bound = params.modulus / 2;
        const reconstructor_type reconstructor(params, bound);
        const auto coordinates = random_coordinates(65, bound);
        check_reconstruction(reconstructor, coordinates, 65);
        check_reconstruction(reconstructor, coordinates, 96, true);
    }
}

BOOST_AUTO_TEST_CASE(wide_export_preserves_split_boundaries_and_field_multiples) {
    const cpp_int characteristic(fq::modulus.backend().to_cpp_int());
    for (std::size_t count = 8; count <= detail::small_prime_fft_primes.size(); ++count) {
        const auto params = make_params(count);
        const cpp_int bound = params.modulus / 2;
        const reconstructor_type reconstructor(params, bound);
        std::vector<cpp_int> values {0, bound, -bound};
        for (unsigned bit : {254, 256, 384, 507, 508, 512}) {
            const cpp_int power = cpp_int(1) << bit;
            const cpp_int multiple = (power / characteristic) * characteristic;
            for (const cpp_int &anchor : {power, multiple}) {
                for (int delta : {-1, 0, 1}) {
                    const cpp_int value = anchor + delta;
                    if (value >= 0 && value <= bound) {
                        values.push_back(value);
                        values.push_back(-value);
                    }
                }
            }
        }
        // Exercise nonzero high and low parts together, including the largest high part.
        const cpp_int max_high = bound >> 384;
        const std::vector<cpp_int> lows {
            0, 1, characteristic - 1, characteristic, characteristic + 1, (cpp_int(1) << 384) - 1};
        for (const cpp_int &high : {cpp_int(1), max_high}) {
            for (const cpp_int &low : lows) {
                const cpp_int value = (high << 384) + low;
                if (value <= bound) {
                    values.push_back(value);
                    values.push_back(-value);
                }
            }
        }
        std::vector<coordinates_type> coordinates(values.size());
        for (std::size_t i = 0; i < coordinates.size(); ++i) {
            for (std::size_t c = 0; c < 12; ++c) {
                coordinates[i][c] = values[(i + c) % values.size()];
            }
        }
        check_reconstruction(reconstructor, coordinates, coordinates.size() + 3, true);
    }
}

BOOST_AUTO_TEST_CASE(cofactor_accumulation_preserves_carries_for_every_active_width) {
    for (std::size_t count = 1; count <= detail::small_prime_fft_primes.size(); ++count) {
        const auto params = make_params(count);
        const cpp_int bound = params.modulus / 2;
        const reconstructor_type reconstructor(params, bound);
        std::vector<cpp_int> values;
        for (unsigned bit : {0, 1, 31, 32, 48, 49}) {
            for (std::size_t parity = 0; parity < 2; ++parity) {
                cpp_int sum = 0;
                for (std::size_t prime = 0; prime < count; ++prime) {
                    const auto q = detail::small_prime_fft_primes[prime].modulus;
                    const auto digit = prime % 2 == parity ? q - 1 : (std::uint64_t(1) << bit) - 1;
                    sum += (params.modulus / q) * digit;
                }
                cpp_int value = canonical(sum, params.modulus);
                if (value > bound) {
                    value -= params.modulus;
                }
                values.push_back(value);
            }
        }
        std::vector<coordinates_type> coordinates(values.size());
        for (std::size_t i = 0; i < coordinates.size(); ++i) {
            for (std::size_t c = 0; c < 12; ++c) {
                coordinates[i][c] = values[(i + c) % values.size()];
            }
        }
        check_reconstruction(reconstructor, coordinates, 32, true);
    }
}

BOOST_AUTO_TEST_CASE(packed_lanes_and_block_tails_preserve_prefixes_and_rounding_modes) {
    restore_rounding restore;
    for (std::size_t count = 1; count <= detail::small_prime_fft_primes.size(); ++count) {
        const auto params = make_params(count);
        const cpp_int bound = params.modulus / 2;
        const reconstructor_type reconstructor(params, bound);
        const auto coordinates = random_coordinates(67, bound);
        for (int mode : {FE_TONEAREST, FE_TOWARDZERO}) {
            BOOST_REQUIRE_EQUAL(std::fesetround(mode), 0);
            for (const std::size_t size : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 31, 32, 33, 34, 35, 63, 64, 65, 66, 67}) {
                check_reconstruction(reconstructor, std::span<const coordinates_type>(coordinates).first(size), 71,
                                     true);
                BOOST_CHECK_EQUAL(std::fegetround(), mode);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(wide_bound_violations_are_rejected_before_field_reduction) {
    const auto params = make_params(12);
    const cpp_int characteristic(fq::modulus.backend().to_cpp_int());
    const cpp_int bound = (cpp_int(1) << 384) - 1;
    const reconstructor_type reconstructor(params, bound);
    std::array<coordinates_type, 1> coordinates {};
    for (int sign : {-1, 1}) {
        coordinates[0][0] = sign * bound;
        check_reconstruction(reconstructor, coordinates, 1, true);
        const std::vector<cpp_int> invalid {bound + 1, (bound / characteristic + 1) * characteristic};
        for (const cpp_int &value : invalid) {
            coordinates[0][0] = sign * value;
            const auto residues = make_residues(coordinates, 12, 1, true);
            std::array<value_type, 1> output {value_type::one()};
            BOOST_CHECK_THROW(reconstructor.reconstruct(output, residues, 1), std::overflow_error);
        }
    }
}

BOOST_AUTO_TEST_CASE(workload_plan_bound_is_recovered_with_both_signs) {
    const auto plan = detail::make_bn254_small_prime_fft_params(14939, 14939);
    BOOST_REQUIRE(plan);
    const reconstructor_type reconstructor(plan->crt, plan->coefficient_bound);
    auto coordinates = random_coordinates(33, plan->coefficient_bound);
    coordinates[0][0] = plan->coefficient_bound;
    coordinates[0][1] = -plan->coefficient_bound;
    check_reconstruction(reconstructor, coordinates, 64, true);
}

BOOST_AUTO_TEST_CASE(conversion_round_trip_preserves_fp12_coordinates_and_requested_prefixes) {
    const auto plan = detail::make_bn254_small_prime_fft_params(65, 33);
    BOOST_REQUIRE(plan);
    const cpp_int characteristic(fq::modulus.backend().to_cpp_int());
    auto coordinates = random_coordinates(65, characteristic - 1);
    std::vector<value_type> input(coordinates.size(), value_type::zero());
    for (std::size_t i = 0; i < input.size(); ++i) {
        for (std::size_t w = 0; w < 2; ++w) {
            for (std::size_t v = 0; v < 3; ++v) {
                for (std::size_t u = 0; u < 2; ++u) {
                    auto &coordinate = coordinates[i][6 * w + 2 * v + u];
                    coordinate = canonical(coordinate, characteristic);
                    input[i].data[w].data[v].data[u] = fq::value_type(fq::integral_type(coordinate.str()));
                }
            }
        }
    }
    const detail::bn254_small_prime_fft_converter converter(plan->crt.prime_count);
    const reconstructor_type reconstructor(plan->crt, plan->coefficient_bound);
    std::vector<double> residues(plan->spectrum_bytes / sizeof(double));
    converter.convert(residues, input, plan->transform_size);
    std::vector<value_type> output(input.size(), value_type::one());
    for (std::size_t count : {std::size_t(65), std::size_t(33), std::size_t(1), std::size_t(0)}) {
        reconstructor.reconstruct(std::span<value_type>(output).first(count), residues, plan->transform_size);
        check_output(std::span<const value_type>(output).first(count),
                     std::span<const coordinates_type>(coordinates).first(count));
    }
    converter.convert(residues, {}, plan->transform_size);
    reconstructor.reconstruct(output, residues, plan->transform_size);
    BOOST_CHECK(std::all_of(output.begin(), output.end(), [](const auto &value) { return value.is_zero(); }));
    reconstructor.reconstruct({}, {}, 0);
}

BOOST_AUTO_TEST_CASE(rounding_modes_preserve_reconstruction_with_reused_and_fresh_contexts) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    const auto params = make_params(12);
    const cpp_int bound = params.modulus / 2;
    const auto coordinates = random_coordinates(33, bound);
    const reconstructor_type original(params, bound);
    for (int mode : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
        BOOST_REQUIRE_EQUAL(std::fesetround(mode), 0);
        const reconstructor_type current(params, bound);
        for (const auto *reconstructor : {&original, &current}) {
            check_reconstruction(*reconstructor, coordinates, 64, true);
            BOOST_CHECK_EQUAL(std::fegetround(), mode);
        }
    }
}

BOOST_AUTO_TEST_CASE(inconsistent_parameters_and_insufficient_capacity_are_rejected) {
    const auto valid = make_params(1);
    for (std::size_t count : {std::size_t(0), std::size_t(13), std::numeric_limits<std::size_t>::max()}) {
        const detail::small_prime_fft_crt_params invalid {count, valid.modulus};
        BOOST_CHECK_THROW(reconstructor_type(invalid, 0), std::invalid_argument);
    }
    for (cpp_int modulus : {cpp_int(0), cpp_int(-1), cpp_int(valid.modulus - 1), cpp_int(valid.modulus + 1)}) {
        const detail::small_prime_fft_crt_params invalid {1, modulus};
        BOOST_CHECK_THROW(reconstructor_type(invalid, 0), std::invalid_argument);
    }
    BOOST_CHECK_THROW(reconstructor_type(valid, -1), std::invalid_argument);
    BOOST_CHECK_THROW(reconstructor_type(valid, (valid.modulus + 1) / 2), std::invalid_argument);
    BOOST_CHECK_THROW(reconstructor_type(valid, cpp_int(1) << 10000), std::invalid_argument);
    const reconstructor_type zero_bound(valid, 0);
    const std::array<coordinates_type, 1> zero {};
    check_reconstruction(zero_bound, zero, 1);
}

BOOST_AUTO_TEST_CASE(centered_coordinates_outside_the_declared_bound_are_rejected) {
    const reconstructor_type reconstructor(make_params(1), 10);
    std::array<coordinates_type, 1> coordinates {};
    for (int value : {-10, 10}) {
        coordinates[0][0] = value;
        check_reconstruction(reconstructor, coordinates, 1);
    }
    for (int value : {-11, 11}) {
        coordinates[0][0] = value;
        const auto residues = make_residues(coordinates, 1, 1);
        std::array<value_type, 1> output {value_type::one()};
        BOOST_CHECK_THROW(reconstructor.reconstruct(output, residues, 1), std::overflow_error);
    }
}

BOOST_AUTO_TEST_CASE(invalid_dimensions_are_rejected_before_output_writes) {
    const reconstructor_type reconstructor(make_params(1), 10);
    std::vector<double> residues(24, 0);
    std::vector<value_type> output(2, value_type::one());
    BOOST_CHECK_THROW(reconstructor.reconstruct(output, std::span<const double>(residues).first(23), 2),
                      std::invalid_argument);
    BOOST_CHECK_THROW(reconstructor.reconstruct(output, residues, 1), std::invalid_argument);
    BOOST_CHECK_THROW(reconstructor.reconstruct({}, residues, 1), std::invalid_argument);
    BOOST_CHECK_THROW(reconstructor.reconstruct(output, std::span<const double>(residues).first(12), 1),
                      std::invalid_argument);
    BOOST_CHECK_THROW(reconstructor.reconstruct({}, residues, std::numeric_limits<std::size_t>::max()),
                      std::invalid_argument);
    BOOST_CHECK_THROW(reconstructor.reconstruct({}, residues, 0), std::invalid_argument);
    BOOST_CHECK_THROW(reconstructor.reconstruct(output, {}, 2), std::invalid_argument);
    BOOST_CHECK(
        std::all_of(output.begin(), output.end(), [](const auto &value) { return value == value_type::one(); }));
}

BOOST_AUTO_TEST_SUITE_END()
