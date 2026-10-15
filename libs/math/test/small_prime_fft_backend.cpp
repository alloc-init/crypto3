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
    constexpr std::size_t stride = 9;
    std::mt19937_64 random(0x46465410);
    for (const auto &prime : detail::small_prime_fft_primes) {
        const auto q = static_cast<std::int64_t>(prime.modulus);
        const detail::bn254_small_prime_fft_pointwise pointwise(prime.modulus);
        std::vector<double> left(12 * stride), right(12 * stride);
        std::vector<coordinates_type> products(stride), squares(stride);
        for (std::size_t i = 0; i < stride; ++i) {
            coordinates_type a {}, b {};
            for (std::size_t c = 0; c < 12; ++c) {
                const std::int64_t x =
                    c % 3 == 0 ? 3 * q - 1 : static_cast<std::int64_t>(random() % (6 * q - 1)) - (3 * q - 1);
                const std::int64_t y =
                    c % 3 == 1 ? 1 - 3 * q : static_cast<std::int64_t>(random() % (6 * q - 1)) - (3 * q - 1);
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

BOOST_AUTO_TEST_CASE(pointwise_invalid_dimensions_throw_before_writes) {
    const detail::bn254_small_prime_fft_pointwise pointwise(detail::small_prime_fft_primes[0].modulus);
    std::vector<double> left(24, 7.0), right(24, 8.0);
    const auto before = left;
    BOOST_CHECK_THROW(pointwise.multiply(left, right, 1, 1.0), std::invalid_argument);
    BOOST_CHECK_THROW(pointwise.multiply(left, std::span<const double>(right).first(23), 2, 1.0),
                      std::invalid_argument);
    BOOST_CHECK_THROW(pointwise.square(std::span<double>(left).first(23), 2, 1.0), std::invalid_argument);
    BOOST_CHECK_THROW(pointwise.square(left, std::numeric_limits<std::size_t>::max(), 1.0), std::invalid_argument);
    BOOST_CHECK_EQUAL_COLLECTIONS(left.begin(), left.end(), before.begin(), before.end());
    BOOST_CHECK_NO_THROW(pointwise.multiply({}, {}, 0, 1.0));
    BOOST_CHECK_NO_THROW(pointwise.square({}, 0, 1.0));
}

BOOST_AUTO_TEST_CASE(all_tower_basis_products_cross_the_conversion_and_crt_pipeline) {
    backend_type backend(1);
    polynomial_type actual;
    for (std::size_t i = 0; i < 12; ++i) {
        coordinates_type a {};
        a[i] = 1;
        const polynomial_type left {field_value(a)};
        for (std::size_t j = 0; j < 12; ++j) {
            coordinates_type b {};
            b[j] = 1;
            const polynomial_type right {field_value(b)};
            backend.multiply(actual, left, right);
            check_polynomial(actual, polynomial_type {field_value(tower_product(a, b))});
        }
        backend.square(actual, left);
        check_polynomial(actual, polynomial_type {field_value(tower_product(a, a))});
    }
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

BOOST_AUTO_TEST_SUITE_END()
