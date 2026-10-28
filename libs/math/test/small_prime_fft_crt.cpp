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

#include <nil/crypto3/math/polynomial/backends/detail/small_prime_fft/bn254_small_prime_fft_crt.hpp>
#include <nil/crypto3/math/polynomial/backends/detail/small_prime_fft/small_prime_fft_plan.hpp>

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

    template<std::size_t Lanes>
    struct scalar_ifma52_operations {
        using limb = nil::crypto3::algebra::fields::detail::fp12_fast::limb;
        using wide = nil::crypto3::algebra::fields::detail::fp12_fast::wide_limb;
        static constexpr std::size_t lanes = Lanes;
        static constexpr limb mask = (limb(1) << 52) - 1;
        using value_type = std::array<limb, lanes>;
        using mask_type = std::array<bool, lanes>;
        static value_type broadcast(limb x) {
            value_type r;
            r.fill(x);
            return r;
        }
        static value_type load(const limb *x) {
            value_type r;
            std::copy_n(x, lanes, r.begin());
            return r;
        }
        static void store(limb *x, const value_type &a) {
            std::copy(a.begin(), a.end(), x);
        }
        static value_type add(value_type a, const value_type &b) {
            for (std::size_t i = 0; i < lanes; ++i) {
                if (a[i] > std::numeric_limits<limb>::max() - b[i]) {
                    throw std::overflow_error("IFMA lane overflow");
                }
                a[i] += b[i];
            }
            return a;
        }
        static value_type subtract(value_type a, const value_type &b) {
            for (std::size_t i = 0; i < lanes; ++i)
                a[i] -= b[i];
            return a;
        }
        static value_type bit_and(value_type a, const value_type &b) {
            for (std::size_t i = 0; i < lanes; ++i)
                a[i] &= b[i];
            return a;
        }
        static value_type shift52(value_type a) {
            for (auto &x : a)
                x >>= 52;
            return a;
        }
        template<unsigned Bits>
        static value_type shift_right(value_type a) {
            for (auto &x : a)
                x >>= Bits;
            return a;
        }
        template<unsigned Bits>
        static value_type shift_left(value_type a) {
            for (auto &x : a) {
                if (x > (std::numeric_limits<limb>::max() >> Bits))
                    throw std::overflow_error("IFMA shift overflow");
                x <<= Bits;
            }
            return a;
        }
        static mask_type less(const value_type &a, const value_type &b) {
            mask_type r;
            for (std::size_t i = 0; i < lanes; ++i)
                r[i] = a[i] < b[i];
            return r;
        }
        static value_type select(const mask_type &mask, value_type yes, const value_type &no) {
            for (std::size_t i = 0; i < lanes; ++i)
                if (!mask[i])
                    yes[i] = no[i];
            return yes;
        }
        static value_type multiply_add_low(value_type sum, const value_type &a, const value_type &b) {
            value_type product;
            for (std::size_t i = 0; i < lanes; ++i)
                product[i] = limb(wide(a[i] & mask) * (b[i] & mask)) & mask;
            return add(sum, product);
        }
        static value_type multiply_add_high(value_type sum, const value_type &a, const value_type &b) {
            value_type product;
            for (std::size_t i = 0; i < lanes; ++i)
                product[i] = limb((wide(a[i] & mask) * (b[i] & mask)) >> 52);
            return add(sum, product);
        }
    };

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
                                      std::size_t stride, bool lazy = false, bool weighted = false) {
        // Padding is intentionally not a valid residue: reconstruction must read only the requested prefix.
        std::vector<double> result(prime_count * 12 * stride, std::numeric_limits<double>::quiet_NaN());
        for (std::size_t prime = 0; prime < prime_count; ++prime) {
            const cpp_int modulus = detail::small_prime_fft_primes[prime].modulus;
            cpp_int inverse = 1;
            if (weighted) {
                // Compute the weight independently of the reconstructor's cached inverses.
                const cpp_int cofactor = make_params(prime_count).modulus / modulus;
                inverse = boost::multiprecision::powm(cpp_int(cofactor % modulus), cpp_int(modulus - 2), modulus);
            }
            for (std::size_t coordinate = 0; coordinate < 12; ++coordinate) {
                for (std::size_t i = 0; i < coordinates.size(); ++i) {
                    cpp_int residue = canonical(coordinates[i][coordinate] * inverse, modulus);
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
                            const cpp_int stored(
                                output[i].data[w].data[v].data[u].data.backend().base_data().to_cpp_int());
                            BOOST_CHECK_LT(stored, characteristic);
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

        const auto weighted = make_residues(coordinates, reconstructor.prime_count(), stride, lazy, true);
        std::fill(output.begin(), output.end(), sentinel);
        reconstructor.reconstruct_weighted(output, weighted, stride);
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
    template<typename Operations>
    void check_crt_accumulation() {
        using params = reconstructor_type::field_type::extension_policy::fast_params;
        using accumulator_type = detail::small_prime_fft_crt_accumulator<params, Operations>;
        using radix = typename accumulator_type::radix;
        const cpp_int p(fq::modulus.backend().to_cpp_int());
        std::mt19937_64 random(0x52c47);
        for (std::size_t count = 1; count <= detail::small_prime_fft_primes.size(); ++count) {
            accumulator_type accumulator;
            std::array<cpp_int, Operations::lanes> expected {};
            for (std::size_t prime = 0; prime < count; ++prime) {
                const auto q = detail::small_prime_fft_primes[prime].modulus;
                cpp_int weight = p - 1;
                if (prime % 3 != 0) {
                    weight = 0;
                    for (unsigned i = 0; i < 4; ++i)
                        weight = (weight << 64) + random();
                    weight %= p;
                }
                std::array<std::uint64_t, 4> words {};
                boost::multiprecision::export_bits(weight, words.begin(), 64, false);
                const auto packed = radix::pack(words.data());
                std::array<std::uint64_t, Operations::lanes> digits;
                for (std::size_t lane = 0; lane < Operations::lanes; ++lane) {
                    digits[lane] = lane == 0 ? 0 : lane == 1 ? q - 1 : random() % q;
                    expected[lane] += weight * digits[lane];
                }
                accumulator.add(packed, Operations::load(digits.data()));
            }
            std::array<std::array<std::uint64_t, Operations::lanes>, 6> columns;
            for (std::size_t digit = 0; digit < columns.size(); ++digit)
                Operations::store(columns[digit].data(), accumulator.columns[digit]);
            std::array<cpp_int, Operations::lanes> correction;
            std::array<std::array<std::uint64_t, Operations::lanes>, 5> correction_digits;
            for (std::size_t lane = 0; lane < Operations::lanes; ++lane) {
                cpp_int actual = 0;
                for (std::size_t i = columns.size(); i-- > 0;)
                    actual = (actual << 52) + columns[i][lane];
                BOOST_CHECK_EQUAL(actual, expected[lane]);
                correction[lane] = lane % 2 ? p - 1 : cpp_int(lane);
                cpp_int c = correction[lane];
                for (std::size_t i = 0; i < 5; ++i) {
                    correction_digits[i][lane] = (c & radix::mask).template convert_to<std::uint64_t>();
                    c >>= 52;
                }
            }
            std::array<typename Operations::value_type, 5> vector_correction;
            for (std::size_t i = 0; i < 5; ++i)
                vector_correction[i] = Operations::load(correction_digits[i].data());
            const auto reduced = accumulator.reduce(vector_correction);
            std::array<std::array<std::uint64_t, Operations::lanes>, 5> digits;
            for (std::size_t i = 0; i < 5; ++i)
                Operations::store(digits[i].data(), reduced[i]);
            const cpp_int s_inverse = boost::multiprecision::powm(cpp_int((cpp_int(1) << 260) % p), cpp_int(p - 2), p);
            for (std::size_t lane = 0; lane < Operations::lanes; ++lane) {
                cpp_int actual = 0;
                for (std::size_t i = 5; i-- > 0;) {
                    BOOST_CHECK_LE(digits[i][lane], radix::mask);
                    actual = (actual << 52) + digits[i][lane];
                }
                BOOST_CHECK_LT(actual, p);
                BOOST_CHECK_EQUAL(actual, ((expected[lane] + correction[lane]) * s_inverse) % p);
            }
        }
    }
    template<typename Operations>
    void check_crt_certificates() {
        using params = reconstructor_type::field_type::extension_policy::fast_params;
        using accumulator_type = detail::small_prime_fft_crt_accumulator<params, Operations>;
        constexpr std::uint64_t scale = std::uint64_t(1) << 60;
        std::mt19937_64 random(0xc471f1);
        for (std::size_t count = 1; count <= detail::small_prime_fft_primes.size(); ++count) {
            const auto plan = make_params(count);
            for (unsigned trial = 0; trial < 32; ++trial) {
                auto lower = Operations::broadcast(0);
                std::array<cpp_int, Operations::lanes> sums {}, expected_lower {};
                for (std::size_t prime = 0; prime < count; ++prime) {
                    const auto q = detail::small_prime_fft_primes[prime].modulus;
                    const auto reciprocal = ((cpp_int(1) << 110) / q).convert_to<std::uint64_t>();
                    std::array<std::uint64_t, Operations::lanes> input;
                    for (std::size_t lane = 0; lane < Operations::lanes; ++lane) {
                        input[lane] = trial == 0 ? 0 : trial == 1 ? q - 1 : random() % q;
                        sums[lane] += cpp_int(input[lane]) * (plan.modulus / q);
                        expected_lower[lane] += (cpp_int(input[lane]) * reciprocal) >> 50;
                    }
                    lower = Operations::add(
                        lower, accumulator_type::certificate_term(Operations::load(input.data()), reciprocal));
                }
                std::array<std::uint64_t, Operations::lanes> actual_lower;
                Operations::store(actual_lower.data(), lower);
                for (std::size_t lane = 0; lane < Operations::lanes; ++lane)
                    BOOST_CHECK_EQUAL(actual_lower[lane], expected_lower[lane]);
                for (std::uint64_t bound : {std::uint64_t(0), std::uint64_t(1), scale / 8, scale / 2 - 1}) {
                    const auto certificate = accumulator_type::certify(lower, bound, count);
                    std::array<std::uint64_t, Operations::lanes> quotients, accepted;
                    Operations::store(quotients.data(), certificate.quotients);
                    Operations::store(accepted.data(), certificate.accepted);
                    for (std::size_t lane = 0; lane < Operations::lanes; ++lane) {
                        const cpp_int quotient = (expected_lower[lane] + scale / 2) / scale;
                        const cpp_int center = quotient * scale;
                        const bool expected = expected_lower[lane] + bound >= center &&
                                              expected_lower[lane] + 2 * count <= center + bound;
                        BOOST_CHECK_EQUAL(quotients[lane], quotient);
                        BOOST_CHECK_EQUAL(bool(accepted[lane]), expected);
                        if (accepted[lane]) {
                            const cpp_int difference = sums[lane] - quotient * plan.modulus;
                            BOOST_CHECK_LE(boost::multiprecision::abs(difference) * scale,
                                           cpp_int(bound) * plan.modulus);
                        }
                    }
                }
            }
        }
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(small_prime_fft_crt_test_suite)

BOOST_AUTO_TEST_CASE(radix52_accumulation_matches_full_integer_sums) {
    check_crt_accumulation<scalar_ifma52_operations<8>>();
    check_crt_certificates<scalar_ifma52_operations<8>>();
#if defined(CRYPTO3_SMALL_PRIME_CRT_IFMA_SELECTED)
    check_crt_accumulation<nil::crypto3::algebra::fields::detail::fp12_fast::ifma52_operations<8>>();
    check_crt_certificates<nil::crypto3::algebra::fields::detail::fp12_fast::ifma52_operations<8>>();
#endif
}

BOOST_AUTO_TEST_CASE(combined_normalization_is_centered_and_exact_in_every_rounding_mode) {
    restore_rounding restore;
    for (int mode : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
        BOOST_REQUIRE_EQUAL(std::fesetround(mode), 0);
        for (std::size_t count = 1; count <= detail::small_prime_fft_primes.size(); ++count) {
            const auto params = make_params(count);
            const reconstructor_type reconstructor(params, params.modulus / 2);
            for (std::size_t prime = 0; prime < count; ++prime) {
                const cpp_int q = detail::small_prime_fft_primes[prime].modulus;
                for (std::size_t size : {1, 2, 4, 8, 32, 256}) {
                    const detail::small_prime_fft_plan plan(prime, size);
                    const double weight = reconstructor.weighted_normalization(prime, plan.size_inverse());
                    BOOST_CHECK_LT(std::abs(weight), q.convert_to<double>() / 2);
                    const cpp_int exact_weight = static_cast<std::int64_t>(weight);
                    BOOST_CHECK_EQUAL(canonical(exact_weight * size * (params.modulus / q), q), 1);
                }
            }
            BOOST_CHECK_THROW(reconstructor.weighted_normalization(count, 1), std::out_of_range);
        }
        BOOST_CHECK_EQUAL(std::fegetround(), mode);
    }
}

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

BOOST_AUTO_TEST_CASE(direct_export_covers_every_centering_correction_and_maximum_digits) {
    for (std::size_t count = 1; count <= detail::small_prime_fft_primes.size(); ++count) {
        const auto params = make_params(count);
        const cpp_int bound = params.modulus / 2;
        const reconstructor_type reconstructor(params, bound);
        for (std::size_t quotient = 0; quotient <= count; ++quotient) {
            cpp_int sum = 0;
            for (std::size_t prime = 0; prime < count; ++prime) {
                const cpp_int q = detail::small_prime_fft_primes[prime].modulus;
                const cpp_int digit = (q - 1) * quotient / count;
                sum += digit * (params.modulus / q);
            }
            BOOST_REQUIRE_EQUAL((sum + bound) / params.modulus, quotient);
            std::array<coordinates_type, 1> coordinates;
            coordinates[0].fill(sum - quotient * params.modulus);
            check_reconstruction(reconstructor, coordinates, 3, true);
        }
    }
}

BOOST_AUTO_TEST_CASE(exact_fallback_preserves_tiny_bounds_and_values_adjacent_to_bounds) {
    for (std::size_t count = 1; count <= detail::small_prime_fft_primes.size(); ++count) {
        const auto params = make_params(count);
        // Include bounds smaller than the certificate's resolution and bounds near M/2.
        const std::vector<cpp_int> bounds {
            0, 1, params.modulus >> 60, (params.modulus >> 60) * 32 + 7, params.modulus / 4, params.modulus / 2};
        for (const auto &bound : bounds) {
            const reconstructor_type reconstructor(params, bound);
            auto coordinates = random_coordinates(33, bound);
            coordinates[0][0] = bound;
            coordinates[0][1] = -bound;
            check_reconstruction(reconstructor, coordinates, 35, true);
            if (bound == params.modulus / 2) {
                continue;
            }
            for (int sign : {-1, 1}) {
                coordinates[0][0] = sign * (bound + 1);
                for (bool weighted : {false, true}) {
                    const auto residues = make_residues(coordinates, count, 35, true, weighted);
                    std::vector<value_type> output(coordinates.size(), value_type::one());
                    if (weighted) {
                        BOOST_CHECK_THROW(reconstructor.reconstruct_weighted(output, residues, 35),
                                          std::overflow_error);
                    } else {
                        BOOST_CHECK_THROW(reconstructor.reconstruct(output, residues, 35), std::overflow_error);
                    }
                    BOOST_CHECK(std::all_of(output.begin(), output.end(),
                                            [](const auto &value) { return value == value_type::one(); }));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(bound_failure_preserves_completed_prefixes_across_block_boundaries) {
    restore_rounding restore;
    const auto plan = detail::make_bn254_small_prime_fft_params(14939, 14939);
    BOOST_REQUIRE(plan);
    for (const cpp_int &bound : {cpp_int(0), plan->coefficient_bound}) {
        const reconstructor_type reconstructor(plan->crt, bound);
        for (int mode : {FE_TONEAREST, FE_UPWARD}) {
            BOOST_REQUIRE_EQUAL(std::fesetround(mode), 0);
            for (std::size_t bad : {0, 1, 7, 8, 9, 15, 16, 17, 31, 32, 33}) {
                for (int sign : {-1, 1}) {
                    std::vector<coordinates_type> coordinates(35);
                    coordinates[bad][0] = sign * (bound + 1);
                    for (bool weighted : {false, true}) {
                        const auto residues = make_residues(coordinates, plan->crt.prime_count, 37, true, weighted);
                        std::vector<value_type> output(coordinates.size(), value_type::one());
                        if (weighted) {
                            BOOST_CHECK_THROW(reconstructor.reconstruct_weighted(output, residues, 37),
                                              std::overflow_error);
                        } else {
                            BOOST_CHECK_THROW(reconstructor.reconstruct(output, residues, 37), std::overflow_error);
                        }
                        for (std::size_t i = 0; i < output.size(); ++i) {
                            BOOST_CHECK(output[i] == (i < bad ? value_type::zero() : value_type::one()));
                        }
                    }
                }
            }
        }
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
    for (auto reconstruct : {&reconstructor_type::reconstruct, &reconstructor_type::reconstruct_weighted}) {
        BOOST_CHECK_THROW((reconstructor.*reconstruct)(output, std::span<const double>(residues).first(23), 2),
                          std::invalid_argument);
        BOOST_CHECK_THROW((reconstructor.*reconstruct)(output, residues, 1), std::invalid_argument);
        BOOST_CHECK_THROW((reconstructor.*reconstruct)({}, residues, 1), std::invalid_argument);
        BOOST_CHECK_THROW((reconstructor.*reconstruct)(output, std::span<const double>(residues).first(12), 1),
                          std::invalid_argument);
        BOOST_CHECK_THROW((reconstructor.*reconstruct)({}, residues, std::numeric_limits<std::size_t>::max()),
                          std::invalid_argument);
        BOOST_CHECK_THROW((reconstructor.*reconstruct)({}, residues, 0), std::invalid_argument);
        BOOST_CHECK_THROW((reconstructor.*reconstruct)(output, {}, 2), std::invalid_argument);
        BOOST_CHECK(
            std::all_of(output.begin(), output.end(), [](const auto &value) { return value == value_type::one(); }));
        BOOST_CHECK_NO_THROW((reconstructor.*reconstruct)({}, {}, 0));
    }
}

BOOST_AUTO_TEST_SUITE_END()
