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

#define BOOST_TEST_MODULE small_prime_fft_backend_test

#include <nil/crypto3/math/polynomial/backends/bn254_small_prime_fft_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_modular_arithmetic.hpp>

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
    namespace math = nil::crypto3::math;
    namespace detail = math::detail;
    namespace arithmetic = math::polynomial_arithmetic;
    using boost::multiprecision::cpp_int;
    using backend_type = arithmetic::bn254_small_prime_fft_backend;
    using polynomial_type = backend_type::polynomial_type;
    using value_type = backend_type::value_type;
    using fq = nil::crypto3::algebra::fields::alt_bn128_base_field<254>;
    using context_type = arithmetic::polynomial_context<backend_type>;
    using reference_backend = arithmetic::schoolbook_backend<value_type>;
    using coordinates_type = std::array<cpp_int, 12>;

    static_assert(arithmetic::PolynomialBackend<backend_type>);
    static_assert(math::detail::SupportsDivrem<backend_type>);
    static_assert(math::detail::SupportsPreparedLowProduct<backend_type>);
    static_assert(math::detail::SupportsPreparedCyclicRemainder<backend_type>);

    template<typename Context, std::size_t Columns>
    concept SupportsMatrixProduct =
        requires(Context &context, std::array<polynomial_type *, 2 * Columns> output,
                 std::array<const polynomial_type *, 4> left, std::array<const polynomial_type *, 2 * Columns> right) {
            { context.template try_multiply_matrix<Columns>(output, left, right) } -> std::same_as<bool>;
        };
    static_assert(SupportsMatrixProduct<context_type, 1> && SupportsMatrixProduct<context_type, 2>);
    static_assert(!SupportsMatrixProduct<context_type, 3>);
    static_assert(!SupportsMatrixProduct<arithmetic::polynomial_context<reference_backend>, 2>);

    struct prepared_product_counts {
        std::size_t preparations = 0;
        std::size_t prepared_products = 0;
        std::size_t ordinary_products = 0;
        std::size_t cyclic_preparations = 0;
        std::size_t cyclic_remainders = 0;
    };

    // Observe dispatch at the backend boundary while keeping all arithmetic in the real backend,
    // following the existing division tests' counting-backend convention.
    struct counting_prepared_backend : backend_type {
        counting_prepared_backend(std::size_t size, prepared_product_counts &counts) :
            backend_type(size), counts(counts) {
        }

        auto prepare_low_product(const polynomial_type &fixed, std::size_t variable_count, std::size_t count) {
            ++counts.preparations;
            return backend_type::prepare_low_product(fixed, variable_count, count);
        }

        bool try_multiply_low_prepared(polynomial_type &output, const polynomial_type &input,
                                       const prepared_low_product_type &prepared, std::size_t count) {
            const bool accepted = backend_type::try_multiply_low_prepared(output, input, prepared, count);
            counts.prepared_products += accepted;
            return accepted;
        }

        void multiply_low(polynomial_type &output, const polynomial_type &left, const polynomial_type &right,
                          std::size_t count) {
            ++counts.ordinary_products;
            backend_type::multiply_low(output, left, right, count);
        }

        auto prepare_cyclic_remainder(const polynomial_type &divisor, std::size_t quotient_count) {
            ++counts.cyclic_preparations;
            return backend_type::prepare_cyclic_remainder(divisor, quotient_count);
        }

        bool try_cyclic_remainder(polynomial_type &output, const polynomial_type &dividend,
                                  const polynomial_type &quotient, const prepared_cyclic_remainder_type &prepared) {
            const bool accepted = backend_type::try_cyclic_remainder(output, dividend, quotient, prepared);
            counts.cyclic_remainders += accepted;
            return accepted;
        }

        prepared_product_counts &counts;
    };

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

    // Expand all 144 coordinate products and reduce u^2, then w^2, then v^3.
    // This oracle uses neither the Karatsuba identities nor the production field multiplication.
    coordinates_type tower_product(const coordinates_type &a, const coordinates_type &b) {
        coordinates_type result {};
        for (std::size_t i = 0; i < 12; ++i) {
            for (std::size_t j = 0; j < 12; ++j) {
                const std::size_t u = i % 2 + j % 2;
                const std::size_t w = i / 6 + j / 6;
                const std::size_t v = (i % 6) / 2 + (j % 6) / 2 + w / 2;
                cpp_int product = a[i] * b[j];
                if (u == 2) {
                    product = -product;
                }
                const std::size_t index = 6 * (w % 2) + 2 * (v % 3) + u % 2;
                if (v < 3) {
                    result[index] += product;
                } else {
                    result[index] += 9 * product;
                    result[index ^ 1] += u % 2 ? -product : product;
                }
            }
        }
        return result;
    }

    value_type field_value(const coordinates_type &coordinates) {
        const cpp_int p(fq::modulus.backend().to_cpp_int());
        value_type result = value_type::zero();
        for (std::size_t c = 0; c < 12; ++c) {
            result.data[c / 6].data[(c % 6) / 2].data[c % 2] =
                fq::value_type(fq::integral_type(canonical(coordinates[c], p).str()));
        }
        return result;
    }

    polynomial_type random_polynomial(std::size_t size, std::mt19937_64 &random) {
        const cpp_int p(fq::modulus.backend().to_cpp_int());
        polynomial_type result(size, value_type::zero());
        for (std::size_t i = 0; i < size; ++i) {
            coordinates_type coordinates {};
            for (std::size_t c = 0; c < 12; ++c) {
                for (unsigned limb = 0; limb < 4; ++limb) {
                    coordinates[c] = (coordinates[c] << 64) + random();
                }
                switch ((i + c) % 7) {
                    case 0:
                        coordinates[c] = p - 1;
                        break;
                    case 1:
                        coordinates[c] = 0;
                        break;
                    case 2:
                        coordinates[c] = 1;
                        break;
                    case 3:
                        coordinates[c] = cpp_int(1) << 253;
                        break;
                    case 4:
                        coordinates[c] = p - 2;
                        break;
                    default:
                        break;
                }
            }
            result[i] = field_value(coordinates);
        }
        return result;
    }

    polynomial_type expected_product(const polynomial_type &left, const polynomial_type &right) {
        polynomial_type result;
        reference_backend {}.multiply(result, left, right);
        return result;
    }

    template<std::size_t Columns>
    std::array<polynomial_type, 2 * Columns>
        expected_matrix_product(const std::array<polynomial_type, 4> &left,
                                const std::array<polynomial_type, 2 * Columns> &right) {
        std::array<polynomial_type, 2 * Columns> result;
        for (std::size_t row = 0; row < 2; ++row) {
            for (std::size_t column = 0; column < Columns; ++column) {
                math::addition(result[row * Columns + column], expected_product(left[2 * row], right[column]),
                               expected_product(left[2 * row + 1], right[Columns + column]));
            }
        }
        return result;
    }

    template<std::size_t Columns>
    bool try_matrix_product(backend_type &backend, std::array<polynomial_type, 2 * Columns> &output,
                            const std::array<polynomial_type, 4> &left,
                            const std::array<polynomial_type, 2 * Columns> &right) {
        std::array<polynomial_type *, 2 * Columns> outputs;
        std::array<const polynomial_type *, 4> left_entries;
        std::array<const polynomial_type *, 2 * Columns> right_entries;
        for (std::size_t i = 0; i < left.size(); ++i) {
            left_entries[i] = &left[i];
        }
        for (std::size_t i = 0; i < output.size(); ++i) {
            outputs[i] = &output[i];
            right_entries[i] = &right[i];
        }
        return backend.try_multiply_matrix<Columns>(outputs, left_entries, right_entries);
    }

    polynomial_type expected_low(polynomial_type product, std::size_t count) {
        if (count == 0) {
            return {value_type::zero()};
        }
        if (product.size() > count) {
            product.resize(count);
        }
        math::condense(product);
        return product;
    }

    void check_polynomial(const polynomial_type &actual, const polynomial_type &expected) {
        BOOST_REQUIRE(!actual.empty());
        BOOST_CHECK(actual == expected);
        BOOST_CHECK(actual.size() == 1 || !actual[actual.size() - 1].is_zero());
    }

    void check_pointwise(std::span<const double> actual, std::span<const coordinates_type> expected,
                         std::uint64_t modulus, double scale) {
        for (std::size_t i = 0; i < expected.size(); ++i) {
            for (std::size_t c = 0; c < 12; ++c) {
                const double value = actual[c * expected.size() + i];
                BOOST_REQUIRE(std::isfinite(value));
                BOOST_REQUIRE(std::abs(value) < static_cast<double>(modulus));
                BOOST_REQUIRE_EQUAL(value, std::trunc(value));
                BOOST_CHECK_EQUAL(canonical(cpp_int(static_cast<std::int64_t>(value)), cpp_int(modulus)),
                                  canonical(expected[i][c] * static_cast<std::int64_t>(scale), cpp_int(modulus)));
            }
        }
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(small_prime_fft_backend_tests)

BOOST_AUTO_TEST_CASE(pointwise_products_and_squares_match_expanded_integer_tower) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    std::mt19937_64 random(0x46465410);
    for (const auto &prime : detail::small_prime_fft_primes) {
        const auto q = static_cast<std::int64_t>(prime.modulus);
        const detail::bn254_small_prime_fft_pointwise pointwise(prime.modulus);
        // Cover both vector widths, every eight-lane remainder, and multiple wide blocks.
        for (const std::size_t stride : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 31, 32, 33}) {
            std::vector<double> left(12 * stride), right(12 * stride);
            std::vector<coordinates_type> products(stride), squares(stride);
            for (std::size_t i = 0; i < stride; ++i) {
                coordinates_type a {}, b {};
                for (std::size_t c = 0; c < 12; ++c) {
                    std::int64_t x;
                    if (i % 4 == 0 && c % 2 == 0) {
                        x = (std::int64_t(1) << 53) - 1;
                    } else {
                        x = c % 3 == 0 ? 3 * q - 1 : static_cast<std::int64_t>(random() % (6 * q - 1)) - (3 * q - 1);
                    }
                    std::int64_t y;
                    if (i % 4 == 1 && c % 2 == 1) {
                        y = 1 - (std::int64_t(1) << 53);
                    } else {
                        y = c % 3 == 1 ? 1 - 3 * q : static_cast<std::int64_t>(random() % (6 * q - 1)) - (3 * q - 1);
                    }
                    a[c] = x;
                    b[c] = y;
                    left[c * stride + i] = static_cast<double>(x);
                    right[c * stride + i] = static_cast<double>(y);
                }
                products[i] = tower_product(a, b);
                squares[i] = tower_product(a, a);
            }
            for (double scale : {0.0, 1.0, -1.0, static_cast<double>((q - 1) / 2), -static_cast<double>((q - 1) / 2)}) {
                std::vector<double> storage(left.size() + 2, -19.0);
                auto actual = std::span<double>(storage).subspan(1, left.size());
                std::copy(left.begin(), left.end(), actual.begin());
                pointwise.multiply(actual, right, stride, scale);
                check_pointwise(actual, products, prime.modulus, scale);
                std::copy(left.begin(), left.end(), actual.begin());
                pointwise.square(actual, stride, scale);
                check_pointwise(actual, squares, prime.modulus, scale);
                std::copy(left.begin(), left.end(), actual.begin());
                pointwise.multiply(actual, actual, stride, scale);
                check_pointwise(actual, squares, prime.modulus, scale);
                BOOST_CHECK_EQUAL(storage.front(), -19.0);
                BOOST_CHECK_EQUAL(storage.back(), -19.0);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(pointwise_invalid_dimensions_throw_before_writes) {
    const detail::bn254_small_prime_fft_pointwise pointwise(detail::small_prime_fft_primes[0].modulus);
    std::vector<double> left(24, 7.0), right(24, 8.0);
    const auto before = left;
    BOOST_CHECK_THROW(pointwise.multiply(left, right, 1, 1.0), std::invalid_argument);
    BOOST_CHECK_THROW(pointwise.multiply(left, std::span<const double>(right).first(23), 2, 1.0),
                      std::invalid_argument);
    BOOST_CHECK_THROW(pointwise.square(std::span<double>(left).first(23), 2, 1.0), std::invalid_argument);
    BOOST_CHECK_THROW(pointwise.multiply(left, right, std::numeric_limits<std::size_t>::max(), 1.0),
                      std::invalid_argument);
    BOOST_CHECK_THROW(pointwise.square(left, std::numeric_limits<std::size_t>::max(), 1.0), std::invalid_argument);
    BOOST_CHECK_THROW(pointwise.multiply_sum(left, right, right, right, right, 1, 1.0), std::invalid_argument);
    BOOST_CHECK_THROW(
        pointwise.multiply_sum(left, right, right, right, std::span<const double>(right).first(23), 2, 1.0),
        std::invalid_argument);
    BOOST_CHECK_EQUAL_COLLECTIONS(left.begin(), left.end(), before.begin(), before.end());
    BOOST_CHECK_NO_THROW(pointwise.multiply({}, {}, 0, 1.0));
    BOOST_CHECK_NO_THROW(pointwise.square({}, 0, 1.0));
    BOOST_CHECK_NO_THROW(pointwise.multiply_sum({}, {}, {}, {}, {}, 0, 1.0));
}

BOOST_AUTO_TEST_CASE(pointwise_sums_match_expanded_integer_tower_and_preserve_aliases) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    std::mt19937_64 random(0x46465421);
    for (const auto &prime : detail::small_prime_fft_primes) {
        const auto q = static_cast<std::int64_t>(prime.modulus);
        const detail::bn254_small_prime_fft_pointwise pointwise(prime.modulus);
        for (const std::size_t stride : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 31, 32, 33}) {
            std::array<std::vector<double>, 4> inputs;
            for (auto &input : inputs) {
                input.resize(12 * stride);
            }
            std::vector<coordinates_type> expected(stride);
            for (std::size_t i = 0; i < stride; ++i) {
                std::array<coordinates_type, 4> values;
                for (std::size_t operand = 0; operand < inputs.size(); ++operand) {
                    for (std::size_t c = 0; c < 12; ++c) {
                        std::int64_t value = static_cast<std::int64_t>(random() % (6 * q - 1)) - (3 * q - 1);
                        if ((c + i + operand) % 7 == 0) {
                            value = (std::int64_t(1) << 53) - 1;
                            if (operand % 2) {
                                value = -value;
                            }
                        }
                        values[operand][c] = value;
                        inputs[operand][c * stride + i] = static_cast<double>(value);
                    }
                }
                expected[i] = tower_product(values[0], values[1]);
                const auto second = tower_product(values[2], values[3]);
                for (std::size_t c = 0; c < 12; ++c) {
                    expected[i][c] += second[c];
                }
            }
            for (const double scale :
                 {0.0, 1.0, -1.0, static_cast<double>((q - 1) / 2), -static_cast<double>((q - 1) / 2)}) {
                std::vector<double> storage(12 * stride + 2, -19.0);
                auto output = std::span<double>(storage).subspan(1, 12 * stride);
                pointwise.multiply_sum(output, inputs[0], inputs[1], inputs[2], inputs[3], stride, scale);
                check_pointwise(output, expected, prime.modulus, scale);
                BOOST_CHECK_EQUAL(storage.front(), -19.0);
                BOOST_CHECK_EQUAL(storage.back(), -19.0);
                for (std::size_t alias = 0; alias < inputs.size(); ++alias) {
                    auto aliased = inputs;
                    pointwise.multiply_sum(aliased[alias], aliased[0], aliased[1], aliased[2], aliased[3], stride,
                                           scale);
                    check_pointwise(aliased[alias], expected, prime.modulus, scale);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(pointwise_boundary_products_and_cancellation_preserve_rounding_modes) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    // Multiple wide blocks and a scalar tail, with all tower coordinates active.
    constexpr std::size_t stride = 65;
    for (const auto &prime : detail::small_prime_fft_primes) {
        const auto q = static_cast<std::int64_t>(prime.modulus);
        const detail::bn254_small_prime_fft_pointwise pointwise(prime.modulus);
        const std::array<std::int64_t, 11> edges {0,
                                                  1,
                                                  -1,
                                                  q - 1,
                                                  1 - q,
                                                  (q - 1) / 2,
                                                  -(q - 1) / 2,
                                                  3 * q - 1,
                                                  1 - 3 * q,
                                                  (std::int64_t(1) << 53) - 1,
                                                  1 - (std::int64_t(1) << 53)};
        std::vector<double> left(12 * stride), right(12 * stride), negative_left(12 * stride);
        std::vector<coordinates_type> products(stride), squares(stride), zero(stride);
        for (std::size_t i = 0; i < stride; ++i) {
            coordinates_type a {}, b {};
            for (std::size_t c = 0; c < 12; ++c) {
                const auto x = edges[(i + 3 * c) % edges.size()];
                const auto y = edges[(i / edges.size() + 7 * c) % edges.size()];
                a[c] = x;
                b[c] = y;
                left[c * stride + i] = static_cast<double>(x);
                right[c * stride + i] = static_cast<double>(y);
                negative_left[c * stride + i] = static_cast<double>(-x);
            }
            products[i] = tower_product(a, b);
            squares[i] = tower_product(a, a);
        }
        for (const int mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
            BOOST_REQUIRE_EQUAL(std::fesetround(mode), 0);
            // The largest allowed scales exercise the bound on unreduced tower outputs.
            for (const double scale : {static_cast<double>((q - 1) / 2), -static_cast<double>((q - 1) / 2)}) {
                auto actual = left;
                pointwise.multiply(actual, right, stride, scale);
                check_pointwise(actual, products, prime.modulus, scale);
                actual = left;
                pointwise.square(actual, stride, scale);
                check_pointwise(actual, squares, prime.modulus, scale);
                pointwise.multiply_sum(actual, left, right, negative_left, right, stride, scale);
                check_pointwise(actual, zero, prime.modulus, scale);
                BOOST_CHECK_EQUAL(std::fegetround(), mode);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(fused_matrix_products_match_schoolbook_and_preserve_aliases) {
    backend_type backend(256);
    std::mt19937_64 random(0x46465422);
    for (const std::size_t size : {3, 17, 64, 8}) {
        std::array<polynomial_type, 4> left, right, actual;
        for (std::size_t i = 0; i < left.size(); ++i) {
            left[i] = random_polynomial(size, random);
            right[i] = random_polynomial(size, random);
        }
        const auto expected = expected_matrix_product<2>(left, right);
        BOOST_REQUIRE(try_matrix_product<2>(backend, actual, left, right));
        for (std::size_t i = 0; i < actual.size(); ++i) {
            check_polynomial(actual[i], expected[i]);
        }
        auto left_alias = left;
        BOOST_REQUIRE(try_matrix_product<2>(backend, left_alias, left_alias, right));
        auto right_alias = right;
        BOOST_REQUIRE(try_matrix_product<2>(backend, right_alias, left, right_alias));
        for (std::size_t i = 0; i < actual.size(); ++i) {
            check_polynomial(left_alias[i], expected[i]);
            check_polynomial(right_alias[i], expected[i]);
        }
        std::array<polynomial_type, 2> pair {right[0], right[1]};
        const auto pair_expected = expected_matrix_product<1>(left, pair);
        BOOST_REQUIRE(try_matrix_product<1>(backend, pair, left, pair));
        check_polynomial(pair[0], pair_expected[0]);
        check_polynomial(pair[1], pair_expected[1]);
    }
    // Entries differing by one coefficient can straddle a power-of-two transform boundary.
    std::array<polynomial_type, 4> left, right, actual;
    for (std::size_t i = 0; i < left.size(); ++i) {
        left[i] = random_polynomial(16 + i % 2, random);
        right[i] = random_polynomial(16 + (i / 2), random);
    }
    const auto expected = expected_matrix_product<2>(left, right);
    BOOST_REQUIRE(try_matrix_product<2>(backend, actual, left, right));
    for (std::size_t i = 0; i < actual.size(); ++i) {
        check_polynomial(actual[i], expected[i]);
    }
}

BOOST_AUTO_TEST_CASE(fused_matrix_products_check_the_sum_bound_and_workspace_before_writes) {
    std::mt19937_64 random(0x46465423);
    backend_type backend(64);
    std::array<polynomial_type, 4> left, right, output;
    for (std::size_t i = 0; i < left.size(); ++i) {
        left[i] = random_polynomial(16, random);
        right[i] = random_polynomial(16, random);
    }
    const auto expected = expected_matrix_product<2>(left, right);
    BOOST_REQUIRE(try_matrix_product<2>(backend, output, left, right));
    for (std::size_t i = 0; i < output.size(); ++i) {
        check_polynomial(output[i], expected[i]);
    }
    const auto before = output;
    const auto check_unchanged = [&]() {
        for (std::size_t i = 0; i < output.size(); ++i) {
            check_polynomial(output[i], before[i]);
        }
    };
    const std::size_t bytes = 8 * 12 * backend.prime_count() * 32 * sizeof(double);
    backend_type exact_budget(64, std::numeric_limits<std::size_t>::max(), bytes);
    BOOST_REQUIRE(try_matrix_product<2>(exact_budget, output, left, right));
    check_unchanged();
    for (const std::size_t budget : {std::size_t(0), bytes - 1}) {
        backend_type limited(64, std::numeric_limits<std::size_t>::max(), budget);
        BOOST_CHECK(!try_matrix_product<2>(limited, output, left, right));
        check_unchanged();
    }
    // Both individual products fit, but 17+17 terms exceed the plan's 32-term CRT bound.
    for (std::size_t i = 0; i < left.size(); ++i) {
        left[i] = random_polynomial(17, random);
        right[i] = random_polynomial(17, random);
    }
    BOOST_CHECK(!try_matrix_product<2>(backend, output, left, right));
    check_unchanged();
    for (const std::size_t size : {1, 2, 65}) {
        left[0] = random_polynomial(size, random);
        BOOST_CHECK(!try_matrix_product<2>(backend, output, left, right));
        check_unchanged();
    }
    left[0] = polynomial_type(3, value_type::zero());
    BOOST_CHECK(!try_matrix_product<2>(backend, output, left, right));
    check_unchanged();

    std::array<polynomial_type *, 4> outputs {&output[0], &output[1], &output[2], &output[3]};
    std::array<const polynomial_type *, 4> inputs {&right[0], &right[1], &right[2], &right[3]};
    outputs[3] = outputs[0];
    BOOST_CHECK(!backend.try_multiply_matrix<2>(outputs, inputs, inputs));
    check_unchanged();
    outputs[3] = nullptr;
    BOOST_CHECK(!backend.try_multiply_matrix<2>(outputs, inputs, inputs));
    check_unchanged();
    outputs[3] = &output[3];
    inputs[3] = nullptr;
    BOOST_CHECK(!backend.try_multiply_matrix<2>(outputs, inputs, inputs));
    check_unchanged();

    backend_type large_backend(512);
    for (std::size_t i = 0; i < left.size(); ++i) {
        left[i] = random_polynomial(i == 0 ? 3 : 40, random);
        right[i] = random_polynomial(i == 0 ? 3 : 40, random);
    }
    BOOST_CHECK(!try_matrix_product<2>(large_backend, output, left, right));
    check_unchanged();
}

BOOST_AUTO_TEST_CASE(fused_matrix_products_share_identical_inputs_and_survive_backend_moves) {
    std::mt19937_64 random(0x46465424);
    polynomial_type input = random_polynomial(16, random);
    const auto product = expected_product(input, input);
    polynomial_type expected;
    math::addition(expected, product, product);
    backend_type ordinary(64);
    const std::size_t one_spectrum_bytes = 12 * ordinary.prime_count() * 32 * sizeof(double);
    backend_type backend(64, std::numeric_limits<std::size_t>::max(), one_spectrum_bytes);
    std::array<polynomial_type, 4> output;
    const std::array<polynomial_type *, 4> outputs {&input, &output[1], &output[2], &output[3]};
    const std::array<const polynomial_type *, 4> inputs {&input, &input, &input, &input};
    BOOST_REQUIRE(backend.try_multiply_matrix<2>(outputs, inputs, inputs));
    check_polynomial(input, expected);
    for (std::size_t i = 1; i < output.size(); ++i) {
        check_polynomial(output[i], expected);
    }
    input = random_polynomial(16, random);
    auto moved = std::move(backend);
    BOOST_CHECK(!backend.try_multiply_matrix<2>(outputs, inputs, inputs));
    const auto square = expected_product(input, input);
    math::addition(expected, square, square);
    BOOST_REQUIRE(moved.try_multiply_matrix<2>(outputs, inputs, inputs));
    check_polynomial(input, expected);
    for (std::size_t i = 1; i < output.size(); ++i) {
        check_polynomial(output[i], expected);
    }
}

BOOST_AUTO_TEST_CASE(fused_matrix_products_preserve_cancellation_and_rounding_modes) {
    restore_rounding restore;
    std::mt19937_64 random(0x46465425);
    backend_type backend(64);
    std::array<polynomial_type, 4> left, right, output;
    left[0] = left[1] = random_polynomial(15, random);
    left[2] = left[3] = random_polynomial(16, random);
    right[0] = random_polynomial(17, random);
    right[1] = random_polynomial(16, random);
    const polynomial_type zero {value_type::zero()};
    math::subtraction(right[2], zero, right[0]);
    math::subtraction(right[3], zero, right[1]);
    for (const int mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
        BOOST_REQUIRE_EQUAL(std::fesetround(mode), 0);
        BOOST_REQUIRE(try_matrix_product<2>(backend, output, left, right));
        BOOST_CHECK_EQUAL(std::fegetround(), mode);
        for (const auto &result : output) {
            check_polynomial(result, zero);
        }
    }
}

BOOST_AUTO_TEST_CASE(all_tower_basis_products_cross_the_conversion_and_crt_pipeline) {
    backend_type backend(8);
    polynomial_type actual;
    // Degree-two monomials keep these checks on the transform path when scalar and linear
    // operands use direct arithmetic. Their product has its only nonzero coefficient at X^4.
    for (std::size_t i = 0; i < 12; ++i) {
        coordinates_type a {};
        a[i] = 1;
        const polynomial_type left {value_type::zero(), value_type::zero(), field_value(a)};
        for (std::size_t j = 0; j < 12; ++j) {
            coordinates_type b {};
            b[j] = 1;
            const polynomial_type right {value_type::zero(), value_type::zero(), field_value(b)};
            backend.multiply(actual, left, right);
            polynomial_type expected(5, value_type::zero());
            expected[4] = field_value(tower_product(a, b));
            check_polynomial(actual, expected);
        }
        backend.square(actual, left);
        polynomial_type expected(5, value_type::zero());
        expected[4] = field_value(tower_product(a, a));
        check_polynomial(actual, expected);
    }
}

BOOST_AUTO_TEST_CASE(short_operands_preserve_low_products_aliases_and_capacity) {
    backend_type backend(512);
    std::mt19937_64 random(0x46465440);
    for (const std::size_t short_size : {1, 2, 3}) {
        const auto short_operand = random_polynomial(short_size, random);
        for (const std::size_t long_size : {17, 64, 255}) {
            const auto long_operand = random_polynomial(long_size, random);
            const auto full = expected_product(short_operand, long_operand);
            for (const std::size_t count :
                 {std::size_t(0), std::size_t(1), std::size_t(2), std::size_t(3), std::size_t(17), full.size()}) {
                const auto expected = expected_low(full, count);
                auto left_alias = short_operand;
                auto right_alias = long_operand;
                backend.multiply_low(left_alias, left_alias, long_operand, count);
                backend.multiply_low(right_alias, short_operand, right_alias, count);
                check_polynomial(left_alias, expected);
                check_polynomial(right_alias, expected);
                backend.multiply_low(right_alias, long_operand, short_operand, count);
                check_polynomial(right_alias, expected);
            }
        }
    }

    backend_type limited(32);
    const auto linear = random_polynomial(2, random);
    const auto at_limit = random_polynomial(31, random);
    polynomial_type result;
    limited.multiply(result, linear, at_limit);
    check_polynomial(result, expected_product(linear, at_limit));
    const auto before = result;
    BOOST_CHECK_THROW(limited.multiply(result, linear, random_polynomial(32, random)), std::invalid_argument);
    BOOST_CHECK_THROW(
        limited.multiply_low(result, random_polynomial(33, random), polynomial_type {value_type::one()}, 33),
        std::invalid_argument);
    check_polynomial(result, before);
    backend_type moved(std::move(limited));
    BOOST_CHECK_THROW(limited.multiply(result, linear, at_limit), std::invalid_argument);
    check_polynomial(result, before);
}

BOOST_AUTO_TEST_CASE(full_width_products_match_native_fp12_schoolbook) {
    backend_type backend(512);
    std::mt19937_64 random(0x46465411);
    const std::array<std::pair<std::size_t, std::size_t>, 9> shapes {
        {{1, 1}, {1, 17}, {2, 3}, {7, 9}, {8, 9}, {17, 33}, {32, 32}, {65, 31}, {257, 129}}};
    polynomial_type result;
    for (const auto [a, b] : shapes) {
        BOOST_TEST_CONTEXT("shape=" << a << "x" << b) {
            const auto left = random_polynomial(a, random);
            const auto right = random_polynomial(b, random);
            backend.multiply(result, left, right);
            check_polynomial(result, expected_product(left, right));
        }
    }
    // Canonical p-1 coordinates force large integer lifts and signed CRT reconstruction.
    coordinates_type coordinates;
    coordinates.fill(cpp_int(fq::modulus.backend().to_cpp_int()) - 1);
    const polynomial_type dense(17, field_value(coordinates));
    const polynomial_type separate_copy = dense;
    backend.multiply(result, dense, separate_copy);
    check_polynomial(result, expected_product(dense, separate_copy));
}

BOOST_AUTO_TEST_CASE(dedicated_squares_match_multiplication_and_schoolbook) {
    backend_type backend(256);
    std::mt19937_64 random(0x46465412);
    for (const std::size_t size : {1, 2, 7, 16, 33, 65}) {
        const auto input = random_polynomial(size, random);
        const polynomial_type separate_copy = input;
        const auto expected = expected_product(input, input);
        polynomial_type squared, multiplied;
        backend.square(squared, input);
        backend.multiply(multiplied, input, separate_copy);
        check_polynomial(squared, expected);
        check_polynomial(multiplied, expected);
        squared = input;
        backend.square(squared, squared);
        check_polynomial(squared, expected);
    }
}

BOOST_AUTO_TEST_CASE(low_products_trim_inputs_and_do_not_wrap_high_coefficients) {
    backend_type backend(128);
    std::mt19937_64 random(0x46465413);
    const auto left = random_polynomial(17, random);
    const auto right = random_polynomial(33, random);
    const auto full = expected_product(left, right);
    polynomial_type actual;
    for (const std::size_t count : {0, 1, 2, 7, 8, 9, 16, 17, 32, 49, 64}) {
        backend.multiply_low(actual, left, right, count);
        check_polynomial(actual, expected_low(full, count));
        backend.multiply_low(actual, left, left, count);
        check_polynomial(actual, expected_low(expected_product(left, left), count));
    }
    // (1+X^8)^2 mod X^9 = 1+2*X^8. A length-16 cyclic convolution would corrupt its constant term.
    polynomial_type sparse(9, value_type::zero());
    sparse[0] = sparse[8] = value_type::one();
    backend.multiply_low(actual, sparse, sparse, 9);
    check_polynomial(actual, expected_low(expected_product(sparse, sparse), 9));
    backend_type small_backend(4);
    small_backend.multiply_low(actual, left, right, 2);
    check_polynomial(actual, expected_low(full, 2));
}

BOOST_AUTO_TEST_CASE(context_interface_and_all_output_aliases) {
    context_type context {backend_type(32)};
    std::mt19937_64 random(0x46465414);
    const auto left = random_polynomial(7, random);
    const auto right = random_polynomial(9, random);
    const auto expected = expected_product(left, right);
    polynomial_type result = left;
    math::multiplication(result, result, right, context);
    check_polynomial(result, expected);
    result = right;
    math::multiplication(result, left, result, context);
    check_polynomial(result, expected);
    result = left;
    math::multiplication(result, result, result, context);
    check_polynomial(result, expected_product(left, left));
    result = left;
    math::multiply_low(result, result, right, 5, context);
    check_polynomial(result, expected_low(expected, 5));
    result = right;
    math::multiply_low(result, left, result, 5, context);
    check_polynomial(result, expected_low(expected, 5));
    result = left;
    math::multiply_low(result, result, result, 5, context);
    check_polynomial(result, expected_low(expected_product(left, left), 5));
    result = left;
    math::square(result, result, context);
    check_polynomial(result, expected_product(left, left));
}

BOOST_AUTO_TEST_CASE(zeros_constants_and_cancellation_produce_canonical_results) {
    backend_type backend(32);
    const polynomial_type zero {value_type::zero()}, one {value_type::one()};
    const polynomial_type a {value_type::one(), -value_type::one()};
    const polynomial_type b {value_type::one(), value_type::one()};
    polynomial_type result = a;
    backend.multiply(result, zero, b);
    check_polynomial(result, zero);
    backend.multiply(result, a, zero);
    check_polynomial(result, zero);
    backend.square(result, zero);
    check_polynomial(result, zero);
    backend.multiply_low(result, a, b, 0);
    check_polynomial(result, zero);
    backend.multiply(result, a, one);
    check_polynomial(result, a);
    backend.multiply_low(result, a, b, 2);
    check_polynomial(result, one);
    polynomial_type padded = a;
    padded.resize(12, value_type::zero());
    backend.multiply(result, padded, b);
    check_polynomial(result, expected_product(a, b));
    polynomial_type zero_prefix(17, value_type::zero());
    zero_prefix[16] = value_type::one();
    backend.multiply_low(result, zero_prefix, b, 4);
    check_polynomial(result, zero);
}

BOOST_AUTO_TEST_CASE(repeated_size_changes_and_owned_backend_state) {
    backend_type original(256);
    std::mt19937_64 random(0x46465415);
    polynomial_type actual;
    const auto initial = random_polynomial(65, random);
    original.square(actual, initial);
    auto copy = original;
    backend_type moved = std::move(copy);
    for (backend_type *backend : {&original, &moved}) {
        for (const std::size_t size : {33, 1, 17, 2, 65, 7, 33, 1}) {
            const auto left = random_polynomial(size, random);
            const auto right = random_polynomial(size + 1, random);
            const auto expected = expected_product(left, right);
            backend->multiply(actual, left, right);
            check_polynomial(actual, expected);
            backend->square(actual, right);
            check_polynomial(actual, expected_product(right, right));
            backend->multiply_low(actual, left, right, 3);
            check_polynomial(actual, expected_low(expected, 3));
        }
    }
}

BOOST_AUTO_TEST_CASE(rounding_modes_preserve_pipeline_results_and_environment) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    backend_type original(16);
    std::mt19937_64 random(0x46465416);
    const auto left = random_polynomial(3, random);
    const auto right = random_polynomial(5, random);
    const auto product = expected_product(left, right);
    const auto squared = expected_product(right, right);
    for (const int mode : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
        BOOST_REQUIRE_EQUAL(std::fesetround(mode), 0);
        backend_type constructed_in_mode(16);
        for (backend_type *backend : {&original, &constructed_in_mode}) {
            polynomial_type actual;
            backend->multiply(actual, left, right);
            check_polynomial(actual, product);
            backend->square(actual, right);
            check_polynomial(actual, squared);
            backend->multiply_low(actual, left, right, 3);
            check_polynomial(actual, expected_low(product, 3));
            BOOST_CHECK_EQUAL(std::fegetround(), mode);
        }
    }
}

BOOST_AUTO_TEST_CASE(configured_capacity_and_byte_limits_are_checked) {
    for (const std::size_t size :
         {std::size_t(0), std::size_t(3), std::size_t(6), std::numeric_limits<std::size_t>::max()}) {
        BOOST_CHECK_THROW(backend_type {size}, std::invalid_argument);
    }
    BOOST_CHECK_THROW(backend_type(std::size_t(1) << 20, 0), std::invalid_argument);
    const auto params = detail::make_bn254_small_prime_fft_params(16, 17);
    BOOST_REQUIRE(params.has_value());
    BOOST_CHECK_THROW(backend_type(32, params->spectrum_bytes - 1), std::invalid_argument);
    backend_type backend(32, params->spectrum_bytes);
    BOOST_CHECK_EQUAL(backend.max_transform_size(), 32);
    BOOST_CHECK_EQUAL(backend.prime_count(), params->crt.prime_count);
    std::mt19937_64 random(0x46465417);
    const auto large = random_polynomial(17, random);
    polynomial_type result {value_type::one()};
    const auto before = result;
    BOOST_CHECK_THROW(backend.multiply(result, large, large), std::invalid_argument);
    BOOST_CHECK_THROW(backend.square(result, large), std::invalid_argument);
    BOOST_CHECK_THROW(backend.multiply_low(result, large, large, 17), std::invalid_argument);
    check_polynomial(result, before);
    polynomial_type alias = large;
    BOOST_CHECK_THROW(backend.square(alias, alias), std::invalid_argument);
    check_polynomial(alias, large);
    backend.multiply_low(result, large, large, 16);
    check_polynomial(result, expected_low(expected_product(large, large), 16));
}

BOOST_AUTO_TEST_CASE(existing_modular_multiplication_uses_the_new_backend) {
    arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    context_type context {backend_type(64), options};
    arithmetic::polynomial_context<reference_backend> reference {reference_backend {}, options};
    std::mt19937_64 random(0x46465418);
    auto divisor = random_polynomial(9, random);
    divisor[8] = value_type::one();
    math::polynomial_divisor_context<backend_type> fixed(divisor, 8, context);
    math::polynomial_divisor_context<reference_backend> reference_fixed(divisor, 8, reference);
    for (unsigned trial = 0; trial < 3; ++trial) {
        const auto left = random_polynomial(8, random);
        const auto right = random_polynomial(8, random);
        polynomial_type actual, expected;
        math::mulmod(expected, left, right, reference_fixed, reference);
        math::mulmod(actual, left, right, fixed, context);
        check_polynomial(actual, expected);
        actual = left;
        math::mulmod(actual, actual, right, fixed, context);
        check_polynomial(actual, expected);
        math::squaremod(expected, left, reference_fixed, reference);
        actual = left;
        math::squaremod(actual, actual, fixed, context);
        check_polynomial(actual, expected);
    }
}

BOOST_AUTO_TEST_CASE(prepared_products_match_native_and_ordinary_low_products) {
    context_type context {backend_type(256)};
    std::mt19937_64 random(0x46465420);
    const std::array<std::pair<std::size_t, std::size_t>, 5> shapes {{{1, 1}, {3, 5}, {7, 9}, {17, 33}, {65, 31}}};
    for (const auto [a, b] : shapes) {
        const auto input = random_polynomial(a, random);
        const auto fixed = random_polynomial(b, random);
        for (const std::size_t count : {0, 1, 2, 7, 8, 9, 16, 17, 32, 49, 96}) {
            BOOST_TEST_CONTEXT("shape=" << a << "x" << b << ", count=" << count) {
                const auto prepared = context.prepare_low_product(fixed, input.size(), count);
                BOOST_REQUIRE(prepared.has_value());
                BOOST_CHECK_EQUAL(prepared->coefficient_count(), count);
                const std::size_t size = count == 0 ? 0 : std::bit_ceil(std::min(a, count) + std::min(b, count) - 1);
                BOOST_CHECK_EQUAL(prepared->transform_size(), size);
                polynomial_type variable = input;
                for (unsigned repeat = 0; repeat < 3; ++repeat) {
                    const auto expected = expected_low(expected_product(variable, fixed), count);
                    polynomial_type ordinary, actual;
                    context.multiply_low(ordinary, variable, fixed, count);
                    check_polynomial(ordinary, expected);
                    BOOST_REQUIRE(context.try_multiply_low_prepared(actual, variable, *prepared, count));
                    check_polynomial(actual, expected);
                    actual = variable;
                    BOOST_REQUIRE(context.try_multiply_low_prepared(actual, actual, *prepared, count));
                    check_polynomial(actual, expected);
                    variable[0] += value_type::one();
                }
            }
        }
    }
    // Large canonical lifts exercise signed CRT through the prepared path as well.
    coordinates_type coordinates;
    coordinates.fill(cpp_int(fq::modulus.backend().to_cpp_int()) - 1);
    const polynomial_type dense(17, field_value(coordinates));
    const auto prepared = context.prepare_low_product(dense, dense.size(), dense.size());
    BOOST_REQUIRE(prepared.has_value());
    polynomial_type actual;
    BOOST_REQUIRE(context.try_multiply_low_prepared(actual, dense, *prepared, dense.size()));
    check_polynomial(actual, expected_low(expected_product(dense, dense), dense.size()));
}

BOOST_AUTO_TEST_CASE(prepared_products_handle_zero_and_truncated_prefixes) {
    context_type context {backend_type(4)};
    const polynomial_type zero {value_type::zero()}, one {value_type::one()};
    polynomial_type input(17, value_type::one()), actual;
    const auto prepared = context.prepare_low_product(input, input.size(), 2);
    BOOST_REQUIRE(prepared.has_value());
    BOOST_CHECK_EQUAL(prepared->transform_size(), 4);
    BOOST_REQUIRE(context.try_multiply_low_prepared(actual, input, *prepared, 2));
    check_polynomial(actual, expected_low(expected_product(input, input), 2));

    context_type tiny {backend_type(1)};
    // Zero variable prefixes need no compatible transform, even in a smaller backend.
    BOOST_REQUIRE(tiny.try_multiply_low_prepared(actual, zero, *prepared, 2));
    check_polynomial(actual, zero);
    polynomial_type zero_prefix(17, value_type::zero());
    zero_prefix[16] = value_type::one();
    BOOST_REQUIRE(tiny.try_multiply_low_prepared(actual, zero_prefix, *prepared, 2));
    check_polynomial(actual, zero);
    for (const auto &fixed : {zero, zero_prefix}) {
        const auto empty = context.prepare_low_product(fixed, input.size(), 2);
        BOOST_REQUIRE(empty.has_value());
        BOOST_CHECK_EQUAL(empty->transform_size(), 0);
        BOOST_CHECK_EQUAL(empty->storage_bytes(), 0);
        BOOST_REQUIRE(tiny.try_multiply_low_prepared(actual, input, *empty, 2));
        check_polynomial(actual, zero);
        actual = one;
        BOOST_CHECK(!tiny.try_multiply_low_prepared(actual, input, *empty, 1));
        check_polynomial(actual, one);
    }
    const auto empty = context.prepare_low_product(input, 0, 0);
    BOOST_REQUIRE(empty.has_value());
    BOOST_CHECK_EQUAL(empty->transform_size(), 0);
    BOOST_CHECK_EQUAL(empty->storage_bytes(), 0);
    BOOST_REQUIRE(tiny.try_multiply_low_prepared(actual, input, *empty, 0));
    check_polynomial(actual, zero);
    BOOST_CHECK(!context.prepare_low_product(zero, 0, 2).has_value());

    // The retained prefixes can end in zeros, although the original polynomials are canonical.
    const polynomial_type fixed {value_type(2), value_type::zero(), value_type(7)};
    const polynomial_type variable {value_type(3), value_type::zero(), value_type::one()};
    const auto trailing = context.prepare_low_product(fixed, variable.size(), 2);
    BOOST_REQUIRE(trailing.has_value());
    BOOST_REQUIRE(context.try_multiply_low_prepared(actual, variable, *trailing, 2));
    check_polynomial(actual, polynomial_type {value_type(6)});

    context_type larger {backend_type(32)};
    polynomial_type sparse(9, value_type::zero());
    sparse[0] = sparse[8] = value_type::one();
    const auto no_wrap = larger.prepare_low_product(sparse, sparse.size(), 9);
    BOOST_REQUIRE(no_wrap.has_value());
    BOOST_CHECK_EQUAL(no_wrap->transform_size(), 32);
    BOOST_REQUIRE(larger.try_multiply_low_prepared(actual, sparse, *no_wrap, 9));
    check_polynomial(actual, expected_low(expected_product(sparse, sparse), 9));
}

BOOST_AUTO_TEST_CASE(prepared_products_require_matching_precision_and_transform) {
    context_type context {backend_type(32)};
    std::mt19937_64 random(0x46465421);
    const auto fixed = random_polynomial(3, random);
    const auto prepared = context.prepare_low_product(fixed, 3, 10);
    BOOST_REQUIRE(prepared.has_value());
    BOOST_CHECK_EQUAL(prepared->transform_size(), 8);
    for (const std::size_t size : {3, 4, 5, 6}) {
        const auto variable = random_polynomial(size, random);
        polynomial_type actual;
        BOOST_REQUIRE(context.try_multiply_low_prepared(actual, variable, *prepared, 10));
        check_polynomial(actual, expected_product(variable, fixed));
    }
    const auto variable = random_polynomial(3, random);
    polynomial_type actual {value_type::one()};
    const auto before = actual;
    for (const std::size_t count : {0, 9, 11}) {
        BOOST_CHECK(!context.try_multiply_low_prepared(actual, variable, *prepared, count));
        check_polynomial(actual, before);
    }
    for (const std::size_t size : {1, 2, 7, 33}) {
        const auto input = random_polynomial(size, random);
        actual = input;
        BOOST_CHECK(!context.try_multiply_low_prepared(actual, actual, *prepared, 10));
        check_polynomial(actual, input);
    }
    context_type too_small {backend_type(4)};
    actual = before;
    BOOST_CHECK(!too_small.try_multiply_low_prepared(actual, variable, *prepared, 10));
    check_polynomial(actual, before);
    BOOST_CHECK(!too_small.prepare_low_product(fixed, 3, 3).has_value());
    BOOST_CHECK(!context.prepare_low_product(fixed, 0, 3).has_value());
    BOOST_CHECK(!context
                     .prepare_low_product(fixed, std::numeric_limits<std::size_t>::max(),
                                          std::numeric_limits<std::size_t>::max())
                     .has_value());
}

BOOST_AUTO_TEST_CASE(prepared_products_own_spectra_across_contexts_and_moves) {
    std::mt19937_64 random(0x46465422);
    const auto fixed = random_polynomial(5, random);
    const auto variable = random_polynomial(7, random);
    const auto prepared = [&] {
        context_type temporary {backend_type(64)};
        polynomial_type source = fixed;
        auto result = temporary.prepare_low_product(source, variable.size(), 17);
        source[0] += value_type::one();
        return result;
    }();
    BOOST_REQUIRE(prepared.has_value());
    auto copy = *prepared;
    const auto moved = std::move(copy);
    const auto expected = expected_product(variable, fixed);
    for (const std::size_t size : {16, 128}) {
        context_type context {backend_type(size)};
        for (const auto *snapshot : {&*prepared, &moved}) {
            polynomial_type actual;
            BOOST_REQUIRE(context.try_multiply_low_prepared(actual, variable, *snapshot, 17));
            check_polynomial(actual, expected);
            auto changed_fixed = fixed;
            changed_fixed[0] += value_type::one();
            const auto changed = context.prepare_low_product(changed_fixed, variable.size(), 17);
            BOOST_REQUIRE(changed.has_value());
            BOOST_REQUIRE(context.try_multiply_low_prepared(actual, variable, *changed, 17));
            check_polynomial(actual, expected_product(variable, changed_fixed));
            BOOST_REQUIRE(context.try_multiply_low_prepared(actual, variable, *snapshot, 17));
            check_polynomial(actual, expected);
        }
        polynomial_type actual = variable;
        BOOST_CHECK(!context.try_multiply_low_prepared(actual, actual, copy, 17));
        check_polynomial(actual, variable);
    }
    backend_type source(32);
    backend_type destination = std::move(source);
    polynomial_type actual = variable;
    BOOST_CHECK(!source.prepare_low_product(fixed, variable.size(), 17).has_value());
    BOOST_CHECK(!source.try_multiply_low_prepared(actual, actual, *prepared, 17));
    check_polynomial(actual, variable);
    BOOST_REQUIRE(destination.try_multiply_low_prepared(actual, variable, *prepared, 17));
    check_polynomial(actual, expected);
}

BOOST_AUTO_TEST_CASE(small_preparations_do_not_retain_large_workspaces) {
    backend_type backend(256);
    polynomial_type actual, large(65, value_type::one());
    backend.multiply(actual, large, large);
    const polynomial_type small(9, value_type::one());
    const auto prepared = backend.prepare_low_product(small, small.size(), small.size());
    BOOST_REQUIRE(prepared.has_value());
    BOOST_CHECK_EQUAL(prepared->transform_size(), 32);
    const std::size_t bytes_per_position = 12 * backend.prime_count() * sizeof(double);
    BOOST_CHECK_GE(prepared->storage_bytes(), 32 * bytes_per_position);
    BOOST_CHECK_LT(prepared->storage_bytes(), 256 * bytes_per_position);
    BOOST_REQUIRE(backend.try_multiply_low_prepared(actual, small, *prepared, small.size()));
    polynomial_type expected(small.size(), value_type::zero());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        expected[i] = value_type(i + 1);
    }
    check_polynomial(actual, expected);
}

BOOST_AUTO_TEST_CASE(prepared_products_survive_rounding_mode_changes) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    context_type context {backend_type(32)};
    std::mt19937_64 random(0x46465423);
    const auto fixed = random_polynomial(5, random);
    const auto variable = random_polynomial(9, random);
    const auto expected = expected_low(expected_product(variable, fixed), 9);
    for (const int prepare_mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
        BOOST_REQUIRE_EQUAL(std::fesetround(prepare_mode), 0);
        context_type preparing {backend_type(32)};
        const auto prepared = preparing.prepare_low_product(fixed, variable.size(), 9);
        BOOST_REQUIRE(prepared.has_value());
        BOOST_CHECK_EQUAL(std::fegetround(), prepare_mode);
        for (const int use_mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
            BOOST_REQUIRE_EQUAL(std::fesetround(use_mode), 0);
            polynomial_type actual;
            BOOST_REQUIRE(context.try_multiply_low_prepared(actual, variable, *prepared, 9));
            check_polynomial(actual, expected);
            BOOST_CHECK_EQUAL(std::fegetround(), use_mode);
        }
    }
}

BOOST_AUTO_TEST_CASE(divisor_context_reuses_prepared_products_and_falls_back) {
    prepared_product_counts counts;
    arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    using counted_context = arithmetic::polynomial_context<counting_prepared_backend>;
    counted_context preparing {counting_prepared_backend(64, counts), options};
    std::mt19937_64 random(0x46465424);
    auto divisor = random_polynomial(9, random);
    divisor[8] = value_type::one();
    const math::polynomial_divisor_context<counting_prepared_backend> fixed(divisor, 7, preparing);
    BOOST_CHECK_EQUAL(counts.preparations, 1);
    BOOST_CHECK_EQUAL(counts.cyclic_preparations, 1);
    auto copy = fixed;
    const auto moved = std::move(copy);
    counted_context context {counting_prepared_backend(32, counts), options};
    const auto remainder = random_polynomial(6, random);
    for (const auto *snapshot : {&fixed, &moved}) {
        for (const std::size_t size : {7, 6, 2, 1, 7}) {
            auto quotient = random_polynomial(size, random);
            quotient[size - 1] = value_type::one();
            auto dividend = expected_product(quotient, divisor);
            math::addition(dividend, dividend, remainder);
            counts = {};
            polynomial_type actual_quotient, actual_remainder;
            math::divrem(actual_quotient, actual_remainder, dividend, *snapshot, context);
            check_polynomial(actual_quotient, quotient);
            check_polynomial(actual_remainder, remainder);
            const std::size_t inverse_accepted = size == 7 ? 1 : 0;
            const std::size_t cyclic_accepted = size > 1 ? 1 : 0;
            BOOST_CHECK_EQUAL(counts.prepared_products, inverse_accepted);
            BOOST_CHECK_EQUAL(counts.cyclic_remainders, cyclic_accepted);
            BOOST_CHECK_EQUAL(counts.ordinary_products, 2 - inverse_accepted - cyclic_accepted);
            BOOST_CHECK_EQUAL(counts.preparations, 0);
            BOOST_CHECK_EQUAL(counts.cyclic_preparations, 0);
            actual_quotient = dividend;
            math::divrem(actual_quotient, actual_remainder, actual_quotient, *snapshot, context);
            check_polynomial(actual_quotient, quotient);
            check_polynomial(actual_remainder, remainder);
            actual_remainder = dividend;
            math::divrem(actual_quotient, actual_remainder, actual_remainder, *snapshot, context);
            check_polynomial(actual_quotient, quotient);
            check_polynomial(actual_remainder, remainder);
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_divisor_context_handles_sparse_fixed_prefixes) {
    prepared_product_counts counts;
    arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    arithmetic::polynomial_context<counting_prepared_backend> context {counting_prepared_backend(32, counts), options};
    std::mt19937_64 random(0x46465425);
    const auto quotient = random_polynomial(7, random);
    const auto remainder = random_polynomial(3, random);
    for (const auto &constant : {value_type::zero(), value_type::one()}) {
        polynomial_type divisor(9, value_type::zero());
        divisor[0] = constant;
        divisor[8] = value_type::one();
        const math::polynomial_divisor_context<counting_prepared_backend> fixed(divisor, 7, context);
        auto dividend = expected_product(quotient, divisor);
        math::addition(dividend, dividend, remainder);
        counts = {};
        polynomial_type actual_quotient, actual_remainder;
        math::divrem(actual_quotient, actual_remainder, dividend, fixed, context);
        check_polynomial(actual_quotient, quotient);
        check_polynomial(actual_remainder, remainder);
        BOOST_CHECK_EQUAL(counts.prepared_products, constant.is_zero() ? 2 : 1);
        BOOST_CHECK_EQUAL(counts.cyclic_remainders, constant.is_zero() ? 0 : 1);
        BOOST_CHECK_EQUAL(counts.ordinary_products, 0);
        BOOST_CHECK_EQUAL(counts.preparations, 0);
        BOOST_CHECK_EQUAL(counts.cyclic_preparations, 0);
    }
}

BOOST_AUTO_TEST_CASE(cyclic_remainders_match_native_division_and_preserve_aliases) {
    context_type context {backend_type(256)};
    std::mt19937_64 random(0x46465430);
    for (const std::size_t degree : {4, 6, 7, 8, 12, 15, 16, 24, 31, 32, 48, 63, 64}) {
        const auto divisor = random_polynomial(degree + 1, random);
        auto quotient = random_polynomial(degree - 1, random);
        const auto remainder = random_polynomial(degree, random);
        const auto prepared = context.prepare_cyclic_remainder(divisor, quotient.size());
        BOOST_REQUIRE(prepared.has_value());
        BOOST_CHECK_EQUAL(prepared->divisor_degree(), degree);
        BOOST_CHECK_EQUAL(prepared->transform_size(), std::bit_ceil(degree));
        for (unsigned repeat = 0; repeat < 3; ++repeat) {
            BOOST_TEST_CONTEXT("degree=" << degree << ", repeat=" << repeat) {
                auto dividend = expected_product(quotient, divisor);
                math::addition(dividend, dividend, remainder);
                polynomial_type expected_quotient, expected_remainder, actual;
                math::division(expected_quotient, expected_remainder, dividend, divisor);
                check_polynomial(expected_quotient, quotient);
                check_polynomial(expected_remainder, remainder);
                BOOST_REQUIRE(context.try_cyclic_remainder(actual, dividend, quotient, *prepared));
                check_polynomial(actual, expected_remainder);
                // All N cells are reconstructed and condensed, so this also checks cancellation above degree(H).
                BOOST_CHECK_LT(actual.size(), divisor.size());
                actual = dividend;
                BOOST_REQUIRE(context.try_cyclic_remainder(actual, actual, quotient, *prepared));
                check_polynomial(actual, remainder);
                actual = quotient;
                BOOST_REQUIRE(context.try_cyclic_remainder(actual, dividend, actual, *prepared));
                check_polynomial(actual, remainder);
                quotient[0] += value_type::one();
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(cyclic_remainders_fold_complete_divisors_with_canonical_coordinates) {
    context_type context {backend_type(16)};
    coordinates_type coordinates;
    coordinates.fill(cpp_int(fq::modulus.backend().to_cpp_int()) - 1);
    const polynomial_type divisor(9, field_value(coordinates)), quotient(7, field_value(coordinates));
    const polynomial_type remainder {value_type(19), field_value(coordinates)};
    const auto prepared = context.prepare_cyclic_remainder(divisor, quotient.size());
    BOOST_REQUIRE(prepared.has_value());
    BOOST_CHECK_EQUAL(prepared->transform_size(), 8);
    auto dividend = expected_product(quotient, divisor);
    math::addition(dividend, dividend, remainder);
    polynomial_type actual;
    BOOST_REQUIRE(context.try_cyclic_remainder(actual, dividend, quotient, *prepared));
    check_polynomial(actual, remainder);

    // H=X^8-1 folds to zero, while V=Q*H+R has nonzero high coefficients that must still be folded.
    polynomial_type vanishing(9, value_type::zero()), sparse_quotient(7, value_type::zero());
    vanishing[0] = -value_type::one();
    vanishing[8] = sparse_quotient[6] = value_type::one();
    const auto zero_spectrum = context.prepare_cyclic_remainder(vanishing, sparse_quotient.size());
    BOOST_REQUIRE(zero_spectrum.has_value());
    BOOST_CHECK_EQUAL(zero_spectrum->transform_size(), 8);
    dividend = expected_product(sparse_quotient, vanishing);
    math::addition(dividend, dividend, remainder);
    BOOST_REQUIRE(context.try_cyclic_remainder(actual, dividend, sparse_quotient, *zero_spectrum));
    check_polynomial(actual, remainder);
    dividend = expected_product(sparse_quotient, vanishing);
    BOOST_REQUIRE(context.try_cyclic_remainder(actual, dividend, sparse_quotient, *zero_spectrum));
    check_polynomial(actual, polynomial_type {value_type::zero()});
}

BOOST_AUTO_TEST_CASE(cyclic_remainders_reject_unsupported_shapes_without_writes) {
    context_type context {backend_type(16)};
    std::mt19937_64 random(0x46465431);
    const auto divisor = random_polynomial(9, random);
    const auto prepared = context.prepare_cyclic_remainder(divisor, 7);
    BOOST_REQUIRE(prepared.has_value());
    BOOST_CHECK(!context.prepare_cyclic_remainder(divisor, 0).has_value());
    BOOST_CHECK(!context.prepare_cyclic_remainder(divisor, 1).has_value());
    BOOST_CHECK(!context.prepare_cyclic_remainder(polynomial_type {value_type(7)}, 7).has_value());
    BOOST_CHECK(!context.prepare_cyclic_remainder(polynomial_type {value_type(7), value_type(11)}, 7).has_value());
    polynomial_type monomial(9, value_type::zero());
    monomial[8] = value_type::one();
    BOOST_CHECK(!context.prepare_cyclic_remainder(monomial, 7).has_value());
    // A degree-five divisor and four-coefficient quotient already need only an eight-point low product.
    BOOST_CHECK(!context.prepare_cyclic_remainder(random_polynomial(6, random), 4).has_value());
    const auto capped = context.prepare_cyclic_remainder(divisor, std::numeric_limits<std::size_t>::max());
    BOOST_REQUIRE(capped.has_value());
    BOOST_CHECK_EQUAL(capped->transform_size(), 8);

    const auto reject = [&](const polynomial_type &dividend, const polynomial_type &quotient) {
        polynomial_type actual {value_type(99)};
        BOOST_CHECK(!context.try_cyclic_remainder(actual, dividend, quotient, *prepared));
        check_polynomial(actual, polynomial_type {value_type(99)});
        actual = dividend;
        BOOST_CHECK(!context.try_cyclic_remainder(actual, actual, quotient, *prepared));
        check_polynomial(actual, dividend);
        actual = quotient;
        BOOST_CHECK(!context.try_cyclic_remainder(actual, dividend, actual, *prepared));
        check_polynomial(actual, quotient);
    };
    for (const std::size_t count : {1, 8, 9}) {
        const auto quotient = random_polynomial(count, random);
        reject(expected_product(quotient, divisor), quotient);
    }
    reject(random_polynomial(8, random), polynomial_type {value_type::zero()});
    const auto quotient = random_polynomial(7, random);
    const auto dividend = expected_product(quotient, divisor);
    reject(dividend, random_polynomial(6, random));
    for (const std::size_t size : {4, 8}) {
        context_type too_small {backend_type(size)};
        BOOST_CHECK(!too_small.prepare_cyclic_remainder(divisor, quotient.size()).has_value());
        polynomial_type actual = dividend;
        BOOST_CHECK(!too_small.try_cyclic_remainder(actual, actual, quotient, *prepared));
        check_polynomial(actual, dividend);
    }
    polynomial_type actual;
    BOOST_REQUIRE(context.try_cyclic_remainder(actual, dividend, quotient, *capped));
    check_polynomial(actual, polynomial_type {value_type::zero()});
}

BOOST_AUTO_TEST_CASE(cyclic_preparations_own_spectra_across_contexts_and_moves) {
    std::mt19937_64 random(0x46465432);
    const auto divisor = random_polynomial(9, random);
    const auto quotient = random_polynomial(7, random);
    const auto remainder = random_polynomial(8, random);
    auto dividend = expected_product(quotient, divisor);
    math::addition(dividend, dividend, remainder);
    const auto prepared = [&] {
        context_type temporary {backend_type(64)};
        auto source = divisor;
        auto result = temporary.prepare_cyclic_remainder(source, quotient.size());
        source[0] += value_type::one();
        return result;
    }();
    BOOST_REQUIRE(prepared.has_value());
    auto copy = *prepared;
    const auto moved = std::move(copy);
    for (const std::size_t size : {16, 128}) {
        context_type context {backend_type(size)};
        for (const auto *snapshot : {&*prepared, &moved}) {
            polynomial_type actual;
            BOOST_REQUIRE(context.try_cyclic_remainder(actual, dividend, quotient, *snapshot));
            check_polynomial(actual, remainder);
            auto changed_divisor = divisor;
            changed_divisor[0] += value_type::one();
            const auto changed = context.prepare_cyclic_remainder(changed_divisor, quotient.size());
            BOOST_REQUIRE(changed.has_value());
            const auto changed_dividend = expected_product(quotient, changed_divisor);
            BOOST_REQUIRE(context.try_cyclic_remainder(actual, changed_dividend, quotient, *changed));
            check_polynomial(actual, polynomial_type {value_type::zero()});
            BOOST_REQUIRE(context.try_cyclic_remainder(actual, dividend, quotient, *snapshot));
            check_polynomial(actual, remainder);
        }
        polynomial_type actual = dividend;
        BOOST_CHECK(!context.try_cyclic_remainder(actual, actual, quotient, copy));
        check_polynomial(actual, dividend);
    }
    backend_type source(32);
    backend_type destination = std::move(source);
    polynomial_type actual = dividend;
    BOOST_CHECK(!source.prepare_cyclic_remainder(divisor, quotient.size()).has_value());
    BOOST_CHECK(!source.try_cyclic_remainder(actual, actual, quotient, *prepared));
    check_polynomial(actual, dividend);
    BOOST_REQUIRE(destination.try_cyclic_remainder(actual, dividend, quotient, *prepared));
    check_polynomial(actual, remainder);
}

BOOST_AUTO_TEST_CASE(small_cyclic_preparations_do_not_retain_large_workspaces) {
    backend_type backend(256);
    polynomial_type actual, large(65, value_type::one());
    backend.square(actual, large);
    const polynomial_type divisor(9, value_type::one()), quotient(7, value_type::one());
    const auto prepared = backend.prepare_cyclic_remainder(divisor, quotient.size());
    BOOST_REQUIRE(prepared.has_value());
    BOOST_CHECK_EQUAL(prepared->transform_size(), 8);
    const std::size_t bytes_per_position = 12 * backend.prime_count() * sizeof(double);
    BOOST_CHECK_GE(prepared->storage_bytes(), 8 * bytes_per_position);
    BOOST_CHECK_LT(prepared->storage_bytes(), 256 * bytes_per_position);
    auto dividend = expected_product(quotient, divisor);
    dividend[0] += value_type(19);
    BOOST_REQUIRE(backend.try_cyclic_remainder(actual, dividend, quotient, *prepared));
    check_polynomial(actual, polynomial_type {value_type(19)});
}

BOOST_AUTO_TEST_CASE(cyclic_remainders_survive_rounding_mode_changes) {
    restore_rounding restore;
    BOOST_REQUIRE_EQUAL(std::fesetround(FE_TONEAREST), 0);
    context_type context {backend_type(16)};
    std::mt19937_64 random(0x46465433);
    const auto divisor = random_polynomial(9, random);
    const auto quotient = random_polynomial(7, random);
    const auto remainder = random_polynomial(8, random);
    auto dividend = expected_product(quotient, divisor);
    math::addition(dividend, dividend, remainder);
    for (const int prepare_mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
        BOOST_REQUIRE_EQUAL(std::fesetround(prepare_mode), 0);
        context_type preparing {backend_type(16)};
        const auto prepared = preparing.prepare_cyclic_remainder(divisor, quotient.size());
        BOOST_REQUIRE(prepared.has_value());
        BOOST_CHECK_EQUAL(std::fegetround(), prepare_mode);
        for (const int use_mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
            BOOST_REQUIRE_EQUAL(std::fesetround(use_mode), 0);
            polynomial_type actual;
            BOOST_REQUIRE(context.try_cyclic_remainder(actual, dividend, quotient, *prepared));
            check_polynomial(actual, remainder);
            BOOST_CHECK_EQUAL(std::fegetround(), use_mode);
        }
    }
}

BOOST_AUTO_TEST_CASE(cyclic_divisor_context_falls_back_outside_reduced_product_range) {
    prepared_product_counts counts;
    arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    arithmetic::polynomial_context<counting_prepared_backend> context {counting_prepared_backend(64, counts), options};
    std::mt19937_64 random(0x46465434);
    const auto divisor = random_polynomial(9, random);
    const math::polynomial_divisor_context<counting_prepared_backend> fixed(divisor, 9, context);
    BOOST_CHECK_EQUAL(counts.preparations, 1);
    BOOST_CHECK_EQUAL(counts.cyclic_preparations, 1);
    const auto remainder = random_polynomial(8, random);
    for (const std::size_t count : {7, 8, 9, 1, 7}) {
        const auto quotient = random_polynomial(count, random);
        auto dividend = expected_product(quotient, divisor);
        math::addition(dividend, dividend, remainder);
        counts = {};
        polynomial_type actual_quotient, actual_remainder;
        math::divrem(actual_quotient, actual_remainder, dividend, fixed, context);
        check_polynomial(actual_quotient, quotient);
        check_polynomial(actual_remainder, remainder);
        const std::size_t inverse_accepted = count == 9 ? 1 : 0;
        const std::size_t cyclic_accepted = count == 7 ? 1 : 0;
        BOOST_CHECK_EQUAL(counts.prepared_products, inverse_accepted);
        BOOST_CHECK_EQUAL(counts.cyclic_remainders, cyclic_accepted);
        BOOST_CHECK_EQUAL(counts.ordinary_products, 2 - inverse_accepted - cyclic_accepted);
        BOOST_CHECK_EQUAL(counts.preparations + counts.cyclic_preparations, 0);
    }
}

BOOST_AUTO_TEST_CASE(divisor_context_keeps_prepared_low_products_when_cyclic_is_not_shorter) {
    prepared_product_counts counts;
    arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    arithmetic::polynomial_context<counting_prepared_backend> context {counting_prepared_backend(16, counts), options};
    std::mt19937_64 random(0x46465435);
    const auto divisor = random_polynomial(6, random);
    const math::polynomial_divisor_context<counting_prepared_backend> fixed(divisor, 4, context);
    BOOST_CHECK_EQUAL(counts.preparations, 2);
    BOOST_CHECK_EQUAL(counts.cyclic_preparations, 1);
    const auto remainder = random_polynomial(5, random);
    for (unsigned repeat = 0; repeat < 3; ++repeat) {
        const auto quotient = random_polynomial(4, random);
        auto dividend = expected_product(quotient, divisor);
        math::addition(dividend, dividend, remainder);
        counts = {};
        polynomial_type actual_quotient, actual_remainder;
        math::divrem(actual_quotient, actual_remainder, dividend, fixed, context);
        check_polynomial(actual_quotient, quotient);
        check_polynomial(actual_remainder, remainder);
        BOOST_CHECK_EQUAL(counts.prepared_products, 2);
        BOOST_CHECK_EQUAL(counts.cyclic_remainders + counts.ordinary_products, 0);
        BOOST_CHECK_EQUAL(counts.preparations + counts.cyclic_preparations, 0);
    }
}

BOOST_AUTO_TEST_SUITE_END()
