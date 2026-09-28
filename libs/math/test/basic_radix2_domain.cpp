//---------------------------------------------------------------------------//
// Copyright (c) 2024 Dmitrii Tabalin <d.tabalin@nil.foundation>
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

#define BOOST_TEST_MODULE basic_radix2_domain_test

#include <vector>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/scalar_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/bls12.hpp>
#include <nil/crypto3/math/algorithms/make_evaluation_domain.hpp>
#include <nil/crypto3/math/domains/basic_radix2_domain.hpp>
#include <nil/crypto3/math/polynomial/polynomial.hpp>
#include <nil/crypto3/math/polynomial/polynomial_dfs.hpp>
#include <nil/crypto3/math/polynomial/operations/evaluate.hpp>
#include <nil/crypto3/math/polynomial/operations/shift.hpp>
#include <nil/crypto3/math/algorithms/unity_root.hpp>
#include <nil/crypto3/math/domains/detail/basic_radix2_domain_aux.hpp>

#include <nil/crypto3/algebra/random_element.hpp>

using namespace nil::crypto3::algebra;
using namespace nil::crypto3::math;

typedef fields::bls12_fr<381> FieldType;
using lagrange_fields = boost::mpl::list<fields::alt_bn128_scalar_field<254>, fields::bls12_fr<381>>;
namespace math_detail = nil::crypto3::math::detail;

BOOST_AUTO_TEST_SUITE(basic_radix2_domain_test_suit)

BOOST_AUTO_TEST_CASE_TEMPLATE(lagrange_matches_independent_reference, Field, lagrange_fields) {
    using value_type = typename Field::value_type;

    for (const std::size_t m : {2, 4, 8, 32}) {
        std::vector<value_type> points;
        math_detail::create_fft_cache<Field>(m, unity_root<Field>(m), points);
        for (const value_type &t : {value_type::zero(), value_type(2), value_type(17)}) {
            BOOST_TEST_CONTEXT("domain size " << m << ", evaluation point " << t) {
                const auto values = math_detail::basic_radix2_evaluate_all_lagrange_polynomials<Field>(m, t);
                value_type vanishing_at_t;
                const auto cached_values =
                    math_detail::basic_radix2_evaluate_all_lagrange_polynomials<Field>(points, t, vanishing_at_t);
                BOOST_REQUIRE_EQUAL(values.size(), m);
                BOOST_CHECK(cached_values == values);
                BOOST_CHECK_EQUAL(vanishing_at_t, t.pow(m) - value_type::one());
                for (std::size_t i = 0; i < m; ++i) {
                    BOOST_CHECK_EQUAL(values[i], evaluate_lagrange_polynomial(points, t, m, i));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(lagrange_domain_points_select_one_basis_polynomial, Field, lagrange_fields) {
    using value_type = typename Field::value_type;

    for (const std::size_t m : {2, 8, 16}) {
        std::vector<value_type> points;
        math_detail::create_fft_cache<Field>(m, unity_root<Field>(m), points);
        for (std::size_t point_index = 0; point_index < m; ++point_index) {
            BOOST_TEST_CONTEXT("domain size " << m << ", point index " << point_index) {
                const auto values =
                    math_detail::basic_radix2_evaluate_all_lagrange_polynomials<Field>(m, points[point_index]);
                value_type vanishing_at_t = value_type::one();
                const auto cached_values = math_detail::basic_radix2_evaluate_all_lagrange_polynomials<Field>(
                    points, points[point_index], vanishing_at_t);
                BOOST_REQUIRE_EQUAL(values.size(), m);
                BOOST_CHECK(cached_values == values);
                BOOST_CHECK(vanishing_at_t.is_zero());
                for (std::size_t i = 0; i < m; ++i) {
                    BOOST_CHECK_EQUAL(values[i], i == point_index ? value_type::one() : value_type::zero());
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(lagrange_singleton_domain, Field, lagrange_fields) {
    using value_type = typename Field::value_type;

    const std::vector<value_type> points {value_type::one()};
    for (const value_type &t : {value_type::zero(), value_type::one(), value_type(17)}) {
        const auto values = math_detail::basic_radix2_evaluate_all_lagrange_polynomials<Field>(1, t);
        value_type vanishing_at_t;
        const auto cached_values =
            math_detail::basic_radix2_evaluate_all_lagrange_polynomials<Field>(points, t, vanishing_at_t);
        BOOST_CHECK(values == points);
        BOOST_CHECK(cached_values == points);
        BOOST_CHECK_EQUAL(vanishing_at_t, t - value_type::one());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(lagrange_reconstructs_a_polynomial_at_size_256, Field, lagrange_fields) {
    using value_type = typename Field::value_type;

    constexpr std::size_t m = 256;
    std::vector<value_type> points;
    math_detail::create_fft_cache<Field>(m, unity_root<Field>(m), points);
    std::vector<value_type> coefficients(m);
    for (std::size_t i = 0; i < m; ++i) {
        coefficients[i] = value_type(i + 1);
    }
    // Evaluate a degree-(m-1) polynomial directly, independently of the domain's FFT and Lagrange routines.
    std::vector<value_type> evaluations(m);
    for (std::size_t i = 0; i < m; ++i) {
        evaluations[i] = evaluate_polynomial(coefficients, points[i], m);
    }
    for (const value_type &t : {value_type::zero(), value_type(17)}) {
        const auto values = math_detail::basic_radix2_evaluate_all_lagrange_polynomials<Field>(m, t);
        BOOST_REQUIRE_EQUAL(values.size(), m);
        value_type reconstructed = value_type::zero();
        for (std::size_t i = 0; i < m; ++i) {
            reconstructed += values[i] * evaluations[i];
        }
        BOOST_CHECK_EQUAL(reconstructed, evaluate_polynomial(coefficients, t, m));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(lagrange_rejects_invalid_domain_sizes, Field, lagrange_fields) {
    using value_type = typename Field::value_type;

    for (const std::size_t m : {0, 3, 6}) {
        BOOST_CHECK_THROW(math_detail::basic_radix2_evaluate_all_lagrange_polynomials<Field>(m, value_type(2)),
                          std::invalid_argument);
        const std::vector<value_type> points(m, value_type::one());
        value_type vanishing_at_t;
        BOOST_CHECK_THROW(
            math_detail::basic_radix2_evaluate_all_lagrange_polynomials<Field>(points, value_type(2), vanishing_at_t),
            std::invalid_argument);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(lagrange_public_interfaces_match_independent_reference, Field, lagrange_fields) {
    using value_type = typename Field::value_type;

    for (const std::size_t m : {2, 8, 32}) {
        const basic_radix2_domain<Field> domain(m);
        const evaluation_domain<Field> &base_domain = domain;
        std::vector<value_type> points(m);
        for (std::size_t i = 0; i < m; ++i) {
            points[i] = domain.get_domain_element(i);
        }
        for (const value_type &t : {value_type::zero(), value_type::one(), value_type(17), points.back()}) {
            BOOST_TEST_CONTEXT("domain size " << m << ", evaluation point " << t) {
                std::vector<value_type> expected(m);
                for (std::size_t i = 0; i < m; ++i) {
                    expected[i] = evaluate_lagrange_polynomial(points, t, m, i);
                }
                BOOST_CHECK(domain.evaluate_all_lagrange_polynomials(t) == expected);
                BOOST_CHECK(base_domain.evaluate_all_lagrange_polynomials(t) == expected);

                value_type vanishing_at_t;
                BOOST_CHECK(domain.evaluate_all_lagrange_polynomials(t, vanishing_at_t) == expected);
                BOOST_CHECK_EQUAL(vanishing_at_t, t.pow(m) - value_type::one());
                value_type base_vanishing_at_t;
                BOOST_CHECK(base_domain.evaluate_all_lagrange_polynomials(t, base_vanishing_at_t) == expected);
                BOOST_CHECK_EQUAL(base_vanishing_at_t, vanishing_at_t);

                // The powers overload still uses the inverse FFT and must keep the same coefficient order.
                std::vector<value_type> t_powers(m, value_type::one());
                for (std::size_t i = 1; i < m; ++i) {
                    t_powers[i] = t_powers[i - 1] * t;
                }
                BOOST_CHECK(domain.evaluate_all_lagrange_polynomials(t_powers.cbegin(), t_powers.cend()) == expected);
                BOOST_CHECK(base_domain.evaluate_all_lagrange_polynomials(t_powers.cbegin(), t_powers.cend()) ==
                            expected);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(lagrange_combined_output_can_alias_input, Field, lagrange_fields) {
    using value_type = typename Field::value_type;

    const basic_radix2_domain<Field> domain(8);
    const evaluation_domain<Field> &base_domain = domain;
    for (const value_type &t : {value_type::zero(), value_type::one(), value_type(17), domain.get_domain_element(7)}) {
        BOOST_TEST_CONTEXT("evaluation point " << t) {
            const auto expected = domain.evaluate_all_lagrange_polynomials(t);
            value_type t_and_vanishing = t;
            BOOST_CHECK(base_domain.evaluate_all_lagrange_polynomials(t_and_vanishing, t_and_vanishing) == expected);
            BOOST_CHECK_EQUAL(t_and_vanishing, t.pow(8) - value_type::one());
        }
    }
}

BOOST_AUTO_TEST_CASE(basic_radix2_domain_benchmark, *boost::unit_test::disabled()) {
    using value_type = FieldType::value_type;
    const std::size_t fft_count = 5;
    const std::array<std::size_t, fft_count> fft_sizes = {1 << 16, 1 << 17, 1 << 18, 1 << 19, 1 << 20};
    std::array<std::vector<value_type>, fft_count> test_data;
    std::chrono::time_point<std::chrono::high_resolution_clock> gen_start(std::chrono::high_resolution_clock::now());
    for (std::size_t i = 0; i < fft_count; ++i) {
        test_data[i].resize(fft_sizes[i]);
        for (std::size_t j = 0; j < fft_sizes[i]; ++j) {
            test_data[i][j] = nil::crypto3::algebra::random_element<FieldType>();
        }
    }
    std::cout << "Generation: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                       gen_start)
                     .count()
              << " ms" << std::endl;

    // manually calculate the power, saving all the intermediate powers
    std::chrono::time_point<std::chrono::high_resolution_clock> cache_start(std::chrono::high_resolution_clock::now());
    std::vector<std::shared_ptr<std::vector<value_type>>> omega_powers(fft_count);
    for (std::size_t i = 0; i < fft_count; i++) {
        omega_powers[i].reset(new std::vector<value_type>);
        omega_powers[i]->resize(fft_sizes[i]);
        (*omega_powers[i])[0] = unity_root<FieldType>(fft_sizes[i]);
        for (std::size_t j = 1; j < fft_sizes[i]; j++) {
            (*omega_powers[i])[j] = (*omega_powers[i])[j - 1] * (*omega_powers[i])[0];
        }
    }
    std::cout << "Cache: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                       cache_start)
                     .count()
              << " ms" << std::endl;

    std::chrono::time_point<std::chrono::high_resolution_clock> start_fft(std::chrono::high_resolution_clock::now());
    for (std::size_t i = 0; i < fft_count; ++i) {
        nil::crypto3::math::detail::basic_radix2_fft<FieldType>(
            test_data[i],
            unity_root<FieldType>(fft_sizes[i]));    // omega_powers[i]);
    }

    std::cout << "Uncached FFT: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                       start_fft)
                     .count()
              << " ms" << std::endl;

    std::chrono::time_point<std::chrono::high_resolution_clock> start_cached(std::chrono::high_resolution_clock::now());
    for (std::size_t i = 0; i < fft_count; ++i) {
        nil::crypto3::math::detail::basic_radix2_fft<FieldType>(
            test_data[i], unity_root<FieldType>(fft_sizes[i]), omega_powers[i]);
    }
    std::cout << "Cached FFT: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                       start_cached)
                     .count()
              << " ms" << std::endl;
}

BOOST_AUTO_TEST_CASE(fft_vs_multiplication_benchmark) {
    using value_type = FieldType::value_type;
    const std::size_t fft_size = 1 << 16;
    std::vector<value_type> test_data(fft_size);
    std::chrono::time_point<std::chrono::high_resolution_clock> gen_start(std::chrono::high_resolution_clock::now());
    for (std::size_t i = 0; i < fft_size; ++i) {
        test_data[i] = nil::crypto3::algebra::random_element<FieldType>();
    }
    std::vector<value_type> duped_data(test_data);
    value_type random_mult = nil::crypto3::algebra::random_element<FieldType>();
    std::cout << "Generation: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                       gen_start)
                     .count()
              << " ms" << std::endl;

    std::chrono::time_point<std::chrono::high_resolution_clock> start_fft(std::chrono::high_resolution_clock::now());
    nil::crypto3::math::detail::basic_radix2_fft<FieldType>(test_data, unity_root<FieldType>(fft_size));
    std::cout << "FFT: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                       start_fft)
                     .count()
              << " ms" << std::endl;

    std::chrono::time_point<std::chrono::high_resolution_clock> start_mult(std::chrono::high_resolution_clock::now());

    for (std::size_t i = 0; i < fft_size; ++i) {
        duped_data[i] *= random_mult;
    }
    std::cout << "Multiplication: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                       start_mult)
                     .count()
              << " ms" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
