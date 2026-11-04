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

#define BOOST_TEST_MODULE distinct_degree_factorization_test

#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/babybear.hpp>
#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/factorization/distinct_degree_factorization.hpp>
#include <nil/crypto3/math/polynomial/factorization/kaltofen_shoup_distinct_degree_factorization.hpp>
#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/random/algebraic_engine.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq_value_type = fq_field_type::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;
    using small_value = fields::babybear::value_type;
    using small_backend = polynomial_arithmetic::schoolbook_backend<small_value>;
    using small_polynomial = small_backend::polynomial_type;

    // The primitive generator makes X^n - a irreducible for the degrees used below: their prime divisors
    // divide |F|-1, and |F| = 1 mod 4. Coprime powers of the generator give distinct such factors.
    small_polynomial certificate_binomial(std::size_t degree, std::size_t generator_power = 1) {
        small_polynomial result(degree + 1, small_value::zero());
        result.front() =
            -small_value(fields::arithmetic_params<fields::babybear>::multiplicative_generator).pow(generator_power);
        result.back() = small_value::one();
        return result;
    }

    template<typename ValueType>
    struct counting_backend_impl : polynomial_arithmetic::schoolbook_backend<ValueType> {
        using base_type = polynomial_arithmetic::schoolbook_backend<ValueType>;
        using polynomial_type = typename base_type::polynomial_type;
        explicit counting_backend_impl(std::size_t &calls) : calls(calls) {
        }
        void multiply(polynomial_type &out, const polynomial_type &a, const polynomial_type &b) const {
            ++calls;
            base_type::multiply(out, a, b);
        }
        void square(polynomial_type &out, const polynomial_type &a) const {
            ++calls;
            base_type::square(out, a);
        }
        void multiply_low(polynomial_type &out, const polynomial_type &a, const polynomial_type &b,
                          std::size_t count) const {
            ++calls;
            base_type::multiply_low(out, a, b, count);
        }
        std::size_t &calls;
    };
    using counting_backend = counting_backend_impl<fq_value_type>;
    using small_counting_backend = counting_backend_impl<small_value>;

    // Select a deterministic dense irreducible using the independent reference DDF. A binomial would always have
    // a gapped derivative-GCD chain and would not exercise the public certificate-to-discovery handoff.
    math::polynomial<fq_value_type> dense_irreducible() {
        using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
        polynomial_arithmetic::polynomial_context<backend_type> context;
        nil::crypto3::random::algebraic_engine<fq_field_type> generator(4819);
        for (std::size_t attempt = 0; attempt < 64; ++attempt) {
            backend_type::polynomial_type input(10, fq_value_type::one());
            for (std::size_t i = 0; i < 9; ++i) {
                input[i] = generator();
            }
            const auto factors = math::distinct_degree_factorization_reference<backend_type>(input, context);
            if (factors.factors.size() == 1 && factors.factors.front().irreducible_factor_degree == 9) {
                return input;
            }
        }
        throw std::runtime_error("deterministic irreducible fixture search exhausted");
    }

    fq12_value_type fq12_value(std::size_t first_coordinate) {
        fq12_value_type value = fq12_value_type::zero();
        for (std::size_t i = 0; i < fq12_field_type::arity; ++i) {
            value.coordinate(i) = fq_value_type(first_coordinate + i);
        }
        return value;
    }

    template<polynomial_arithmetic::PolynomialBackend Backend>
    typename Backend::polynomial_type multiply(Backend &backend, const typename Backend::polynomial_type &left,
                                               const typename Backend::polynomial_type &right) {
        typename Backend::polynomial_type result;
        backend.multiply(result, left, right);
        return result;
    }

    template<polynomial_arithmetic::PolynomialBackend Backend>
    typename Backend::polynomial_type
        reconstruct(Backend &backend,
                    const math::distinct_degree_factorization_result<typename Backend::polynomial_type> &result) {
        typename Backend::polynomial_type reconstructed = {result.leading_coefficient};
        for (const auto &factor : result.factors) {
            reconstructed = multiply(backend, reconstructed, factor.polynomial);
        }
        return reconstructed;
    }

    template<polynomial_arithmetic::PolynomialBackend Backend>
    void check_reference_and_kaltofen_shoup_agree(
        const typename Backend::polynomial_type &input,
        polynomial_arithmetic::polynomial_context<Backend> &reference_context,
        polynomial_arithmetic::polynomial_context<Backend> &kaltofen_shoup_context) {
        const auto reference = math::distinct_degree_factorization_reference<Backend>(input, reference_context);
        const auto kaltofen_shoup =
            math::distinct_degree_factorization_kaltofen_shoup<Backend>(input, kaltofen_shoup_context);
        BOOST_CHECK(kaltofen_shoup == reference);
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(distinct_degree_factorization_test_suite)

BOOST_AUTO_TEST_CASE(residual_certificate_plans_reuse_completed_blocks_and_skip_excluded_targets) {
    struct example {
        std::size_t degree, excluded, jump, construction, remainders, applications;
        std::vector<std::size_t> targets;
    };
    for (const auto &test : std::vector<example> {{13351, 783, 9, 0, 2, 18, {13351, 1027}},
                                                  {11614, 1479, 22, 5, 1, 9, {11614, 5807}},
                                                  {14028, 696, 10, 2, 7, 31, {14028, 7014, 4676, 2004}},
                                                  {10690, 2436, 30, 2, 3, 6, {10690, 5345}},
                                                  {14429, 348, 11, 7, 0, 15, {14429}}}) {
        const auto plan = math::detail::kaltofen_shoup_plan_residual_certificate(test.degree, 87, test.excluded);
        BOOST_CHECK(plan.targets == test.targets);
        BOOST_CHECK_EQUAL(plan.jump_blocks, test.jump);
        BOOST_CHECK_EQUAL(plan.construction_compositions, test.construction);
        BOOST_CHECK_EQUAL(plan.remainder_compositions, test.remainders);
        BOOST_CHECK_EQUAL(plan.jump_compositions, test.applications);
    }
    BOOST_CHECK(math::detail::kaltofen_shoup_plan_residual_certificate(9, 2, 4).targets.empty());
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_plan_residual_certificate(0, 1, 1), std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_plan_residual_certificate(9, 0, 1), std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_plan_residual_certificate(9, 2, 1), std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_plan_residual_certificate(9, 2, 3), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(residual_certificate_proves_prime_prime_power_and_composite_degrees) {
    small_backend backend;
    for (std::size_t degree : {5, 8, 9, 15}) {
        const auto remaining = certificate_binomial(degree);
        const auto original = multiply(backend, remaining, small_polynomial {-small_value::one(), small_value::one()});
        polynomial_arithmetic::polynomial_context<small_backend> context;
        const auto expected = math::distinct_degree_factorization_reference<small_backend>(remaining, context);
        BOOST_REQUIRE_EQUAL(expected.factors.size(), 1);
        BOOST_REQUIRE_EQUAL(expected.factors.front().irreducible_factor_degree, degree);
        const math::polynomial_frobenius_context<small_backend> frobenius(original, context);
        const std::size_t block = degree == 5 ? 1 : 2;
        const math::detail::kaltofen_shoup_frobenius_precomputation<small_backend> table(block, frobenius, context);
        const auto giant = table.baby_step(block);
        const auto saved_giant = giant;
        const auto saved_remaining = remaining;
        std::vector<small_polynomial> saved_babies;
        for (std::size_t i = 0; i <= block; ++i)
            saved_babies.push_back(table.baby_step(i));
        BOOST_CHECK(math::detail::kaltofen_shoup_certify_residual(remaining, giant, block, table,
                                                                  frobenius.divisor_context(), context));
        BOOST_CHECK(giant == saved_giant);
        BOOST_CHECK(remaining == saved_remaining);
        for (std::size_t i = 0; i <= block; ++i)
            BOOST_CHECK(table.baby_step(i) == saved_babies[i]);
    }
}

BOOST_AUTO_TEST_CASE(residual_certificate_detects_both_failure_paths_including_odd_factor_counts) {
    small_backend backend;
    // Mixed degrees fail h[n] = X. Equal-degree pairs and triples pass it and require the proper-degree GCD.
    for (const auto &degrees : std::vector<std::vector<std::size_t>> {{2, 3}, {3, 3}, {3, 3, 3}, {4, 4}}) {
        small_polynomial remaining {small_value::one()};
        for (std::size_t i = 0; i < degrees.size(); ++i)
            remaining = multiply(backend, remaining, certificate_binomial(degrees[i], 1 + 6 * i));
        const auto original = multiply(backend, remaining, small_polynomial {-small_value::one(), small_value::one()});
        polynomial_arithmetic::polynomial_context<small_backend> context;
        const math::polynomial_frobenius_context<small_backend> frobenius(original, context);
        const math::detail::kaltofen_shoup_frobenius_precomputation<small_backend> table(1, frobenius, context);
        const auto giant = table.baby_step(1);
        const auto saved_remaining = remaining;
        BOOST_CHECK(!math::detail::kaltofen_shoup_certify_residual(remaining, giant, 1, table,
                                                                   frobenius.divisor_context(), context));
        BOOST_CHECK(remaining == saved_remaining);
        BOOST_CHECK(giant == table.baby_step(1));
        const math::polynomial_frobenius_context<small_backend> residual_frobenius(remaining, context);
        small_polynomial power {small_value::zero(), small_value::one()};
        for (std::size_t i = 0; i < remaining.degree(); ++i)
            math::frobenius_map(power, power, residual_frobenius, context);
        const small_polynomial x {small_value::zero(), small_value::one()};
        BOOST_CHECK_EQUAL(power == x, degrees.front() == degrees.back());
    }
}

BOOST_AUTO_TEST_CASE(residual_certificate_skips_even_residuals_and_existing_parity_proofs_without_extra_work) {
    small_backend backend;
    struct example {
        small_polynomial remaining;
        bool odd;
    };
    const std::vector<example> examples {{multiply(backend, certificate_binomial(3), certificate_binomial(5)), false},
                                         {certificate_binomial(5), true}};
    for (const auto &test : examples) {
        polynomial_arithmetic::polynomial_context<small_backend> reference;
        const auto expected = math::distinct_degree_factorization_reference<small_backend>(test.remaining, reference);
        BOOST_REQUIRE_EQUAL(expected.factors.size(), test.odd ? 1 : 2);
        const auto original =
            multiply(backend, test.remaining, small_polynomial {-small_value::one(), small_value::one()});
        std::size_t setup_calls = 0;
        polynomial_arithmetic::polynomial_context<small_counting_backend> preparation {
            small_counting_backend(setup_calls)};
        const math::polynomial_frobenius_context<small_counting_backend> frobenius(original, preparation);
        const math::detail::kaltofen_shoup_frobenius_precomputation<small_counting_backend> table(1, frobenius,
                                                                                                  preparation);
        // At this boundary, odd parity proves the degree-five residual irreducible. Even parity rules out a
        // successful certificate on the pair. Unknown parity must retain the ordinary bounded attempt policy.
        for (bool known : {false, true}) {
            std::size_t disabled_calls = 0, enabled_calls = 0;
            for (bool enable : {false, true}) {
                polynomial_arithmetic::polynomial_context_options options;
                options.distinct_degree_certificate_minimum_excluded_degree = enable ? 1 : 0;
                auto &calls = enable ? enabled_calls : disabled_calls;
                polynomial_arithmetic::polynomial_context<small_counting_backend> context(small_counting_backend(calls),
                                                                                          options);
                math::detail::kaltofen_shoup_residual_certificate_state state;
                std::vector<math::distinct_degree_factor<small_polynomial>> groups;
                const auto status = math::detail::kaltofen_shoup_continue_discovery(
                    &groups, test.remaining, 1, table.baby_step(1), table, frobenius.divisor_context(), context,
                    [](const auto &) { return math::factorization_control::continue_factorization; }, state,
                    known ? std::optional<bool>(test.odd) : std::nullopt);
                BOOST_CHECK(status == math::factorization_control::continue_factorization);
                BOOST_CHECK(groups == expected.factors);
                BOOST_CHECK_EQUAL(state.attempts, enable && !known ? 1 : 0);
                BOOST_CHECK_EQUAL(state.last_attempted_degree, enable && !known ? test.remaining.degree() : 0);
            }
            if (known) {
                BOOST_CHECK_EQUAL(enabled_calls, disabled_calls);
            } else {
                BOOST_CHECK_GT(enabled_calls, disabled_calls);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(residual_certificate_budget_preserves_discovery_and_stop_prefixes) {
    small_backend backend;
    auto remaining = multiply(backend, certificate_binomial(3), certificate_binomial(5));
    remaining = multiply(backend, remaining, certificate_binomial(24));
    const auto original = multiply(backend, remaining, small_polynomial {-small_value::one(), small_value::one()});
    polynomial_arithmetic::polynomial_context<small_backend> preparation;
    const auto expected = math::distinct_degree_factorization_reference<small_backend>(remaining, preparation);
    BOOST_REQUIRE_EQUAL(expected.factors.size(), 3);
    const math::polynomial_frobenius_context<small_backend> frobenius(original, preparation);
    const math::detail::kaltofen_shoup_frobenius_precomputation<small_backend> table(1, frobenius, preparation);
    for (const auto parity : {std::optional<bool> {}, std::optional<bool>(true)}) {
        for (std::size_t threshold : {0, 1, 4, 100}) {
            for (std::size_t budget : {0, 1, 2, 3, 100}) {
                for (std::size_t stop_after = 0; stop_after <= 3; ++stop_after) {
                    for (bool collect : {false, true}) {
                        polynomial_arithmetic::polynomial_context_options options;
                        options.distinct_degree_certificate_minimum_excluded_degree = threshold;
                        options.distinct_degree_certificate_maximum_attempts = budget;
                        polynomial_arithmetic::polynomial_context<small_backend> context(small_backend {}, options);
                        math::detail::kaltofen_shoup_residual_certificate_state state;
                        std::vector<math::distinct_degree_factor<small_polynomial>> reported, collected;
                        const auto control = math::detail::kaltofen_shoup_continue_discovery<small_backend>(
                            collect ? &collected : nullptr, remaining, 1, table.baby_step(1), table,
                            frobenius.divisor_context(), context,
                            [&](const auto &group) {
                                reported.push_back(group);
                                return reported.size() == stop_after ?
                                           math::factorization_control::stop_factorization :
                                           math::factorization_control::continue_factorization;
                            },
                            state, parity);
                        BOOST_CHECK_EQUAL(control == math::factorization_control::continue_factorization,
                                          stop_after == 0);
                        BOOST_REQUIRE_EQUAL(reported.size(), stop_after == 0 ? 3 : stop_after);
                        for (std::size_t i = 0; i < reported.size(); ++i)
                            BOOST_CHECK(reported[i] == expected.factors[i]);
                        if (collect)
                            BOOST_CHECK(collected == reported);
                        BOOST_CHECK_LE(state.attempts, budget);
                        if (stop_after == 0) {
                            const auto attempts =
                                threshold == 1 ? std::min(budget, parity ? std::size_t(2) : std::size_t(3)) :
                                threshold == 4 ? std::min(budget, parity ? std::size_t(1) : std::size_t(2)) :
                                                 0;
                            BOOST_CHECK_EQUAL(state.attempts, attempts);
                            // Many empty blocks on each composite residual must not spend another attempt.
                            const auto last_degree =
                                attempts == 0                   ? 0 :
                                threshold == 1 && attempts == 1 ? 32 :
                                !parity && ((threshold == 1 && attempts == 2) || (threshold == 4 && attempts == 1)) ?
                                                                  29 :
                                                                  24;
                            BOOST_CHECK_EQUAL(state.last_attempted_degree, last_degree);
                        }
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(residual_certificate_preserves_its_budget_until_even_parity_changes_to_odd) {
    small_backend backend;
    auto remaining = multiply(backend, certificate_binomial(3), certificate_binomial(3, 7));
    remaining = multiply(backend, remaining, certificate_binomial(5));
    remaining = multiply(backend, remaining, certificate_binomial(24));
    polynomial_arithmetic::polynomial_context<small_backend> preparation;
    const auto expected = math::distinct_degree_factorization_reference<small_backend>(remaining, preparation);
    BOOST_REQUIRE_EQUAL(expected.factors.size(), 3);
    BOOST_REQUIRE_EQUAL(expected.factors.front().polynomial.degree(), 6);
    const math::polynomial_frobenius_context<small_backend> frobenius(remaining, preparation);
    for (std::size_t block : {1, 2}) {
        const math::detail::kaltofen_shoup_frobenius_precomputation<small_backend> table(block, frobenius, preparation);
        for (const auto parity : {std::optional<bool> {}, std::optional<bool>(false)}) {
            for (std::size_t stop_after = 0; stop_after <= 3; ++stop_after) {
                polynomial_arithmetic::polynomial_context_options options;
                options.distinct_degree_certificate_minimum_excluded_degree = 1;
                polynomial_arithmetic::polynomial_context<small_backend> context(small_backend {}, options);
                math::detail::kaltofen_shoup_residual_certificate_state state;
                std::vector<math::distinct_degree_factor<small_polynomial>> reported;
                const auto status = math::detail::kaltofen_shoup_continue_discovery<small_backend>(
                    nullptr, remaining, block, table.baby_step(block), table, frobenius.divisor_context(), context,
                    [&](const auto &group) {
                        reported.push_back(group);
                        // The first group removes two cubics, so parity stays even. Removing the quintic flips it
                        // to odd, allowing the unused attempt to certify the final degree-24 factor.
                        BOOST_CHECK_EQUAL(state.attempts, !parity || reported.size() == 3 ? 1 : 0);
                        return reported.size() == stop_after ? math::factorization_control::stop_factorization :
                                                               math::factorization_control::continue_factorization;
                    },
                    state, parity);
                BOOST_CHECK_EQUAL(status == math::factorization_control::continue_factorization, stop_after == 0);
                BOOST_REQUIRE_EQUAL(reported.size(), stop_after == 0 ? 3 : stop_after);
                BOOST_CHECK(std::equal(reported.begin(), reported.end(), expected.factors.begin()));
                BOOST_CHECK_EQUAL(state.last_attempted_degree,
                                  !parity                            ? 35 :
                                  stop_after == 0 || stop_after == 3 ? 24 :
                                                                       0);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(reference_and_kaltofen_shoup_factorizations_separate_known_irreducible_degrees) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type linear_one = {fq_value_type(1), fq_value_type::one()};
    const polynomial_type linear_two = {fq_value_type(2), fq_value_type::one()};
    const polynomial_type degree_one_factors = multiply(backend, linear_one, linear_two);

    // Three is a multiplicative generator of BN254 Fq, hence it is neither a square nor a cube. Therefore X^2 - 3
    // and X^3 - 3 are irreducible over this field.
    const polynomial_type degree_two_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type degree_three_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                                 fq_value_type::one()};
    const fq_value_type leading_coefficient(11);

    polynomial_type input = multiply(backend, degree_one_factors, degree_two_factor);
    input = multiply(backend, input, degree_three_factor);
    for (fq_value_type &coefficient : input) {
        coefficient *= leading_coefficient;
    }

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const auto result = math::distinct_degree_factorization_reference<backend_type>(input, arithmetic_context);
    const auto kaltofen_shoup_result =
        math::distinct_degree_factorization_kaltofen_shoup<backend_type>(input, arithmetic_context);

    BOOST_CHECK(result.complete);
    BOOST_CHECK(result.leading_coefficient == leading_coefficient);
    BOOST_REQUIRE_EQUAL(result.factors.size(), 3);
    BOOST_CHECK(result.factors[0] == math::distinct_degree_factor<polynomial_type>({degree_one_factors, 1}));
    BOOST_CHECK(result.factors[1] == math::distinct_degree_factor<polynomial_type>({degree_two_factor, 2}));
    BOOST_CHECK(result.factors[2] == math::distinct_degree_factor<polynomial_type>({degree_three_factor, 3}));
    BOOST_CHECK(reconstruct(backend, result) == input);
    BOOST_CHECK(kaltofen_shoup_result == result);
}

BOOST_AUTO_TEST_CASE(staged_reference_and_kaltofen_shoup_factorizations_stop_after_the_reported_degree) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type degree_one_factor = {fq_value_type(1), fq_value_type::one()};
    const polynomial_type degree_two_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type degree_three_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                                 fq_value_type::one()};
    polynomial_type input = multiply(backend, degree_one_factor, degree_two_factor);
    input = multiply(backend, input, degree_three_factor);

    auto make_callback = [](std::size_t &callback_count) {
        return [&callback_count](const math::distinct_degree_factor<polynomial_type> &factor) {
            ++callback_count;
            return factor.irreducible_factor_degree == 2 ? math::factorization_control::stop_factorization :
                                                           math::factorization_control::continue_factorization;
        };
    };

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::size_t reference_callback_count = 0;
    auto reference_callback = make_callback(reference_callback_count);
    const auto reference_result =
        math::distinct_degree_factorization_reference<backend_type>(input, arithmetic_context, reference_callback);

    std::size_t kaltofen_shoup_callback_count = 0;
    auto kaltofen_shoup_callback = make_callback(kaltofen_shoup_callback_count);
    const auto kaltofen_shoup_result = math::distinct_degree_factorization_kaltofen_shoup<backend_type>(
        input, arithmetic_context, kaltofen_shoup_callback);

    BOOST_CHECK(!reference_result.complete);
    BOOST_CHECK_EQUAL(reference_callback_count, 2);
    BOOST_CHECK_EQUAL(kaltofen_shoup_callback_count, 2);
    BOOST_CHECK(kaltofen_shoup_result == reference_result);
    BOOST_REQUIRE_EQUAL(reference_result.factors.size(), 2);
    BOOST_CHECK(reference_result.factors[0] == math::distinct_degree_factor<polynomial_type>({degree_one_factor, 1}));
    BOOST_CHECK(reference_result.factors[1] == math::distinct_degree_factor<polynomial_type>({degree_two_factor, 2}));
}

BOOST_AUTO_TEST_CASE(reference_and_kaltofen_shoup_factorizations_reject_a_repeated_factor) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type factor = {fq_value_type(1), fq_value_type::one()};
    const polynomial_type input = multiply(backend, factor, factor);
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;

    BOOST_CHECK_THROW(math::distinct_degree_factorization_reference<backend_type>(input, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::distinct_degree_factorization_kaltofen_shoup<backend_type>(input, arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(zero_and_constant_polynomials_have_no_distinct_degree_factors) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const auto zero = math::distinct_degree_factorization_reference<backend_type>(
        polynomial_type {fq_value_type::zero()}, arithmetic_context);
    const auto kaltofen_shoup_zero = math::distinct_degree_factorization_kaltofen_shoup<backend_type>(
        polynomial_type {fq_value_type::zero()}, arithmetic_context);
    BOOST_CHECK(zero.complete);
    BOOST_CHECK(zero.leading_coefficient == fq_value_type::zero());
    BOOST_CHECK(zero.factors.empty());
    BOOST_CHECK(kaltofen_shoup_zero == zero);

    const auto constant = math::distinct_degree_factorization_reference<backend_type>(
        polynomial_type {fq_value_type(13)}, arithmetic_context);
    const auto kaltofen_shoup_constant = math::distinct_degree_factorization_kaltofen_shoup<backend_type>(
        polynomial_type {fq_value_type(13)}, arithmetic_context);
    BOOST_CHECK(constant.complete);
    BOOST_CHECK(constant.leading_coefficient == fq_value_type(13));
    BOOST_CHECK(constant.factors.empty());
    BOOST_CHECK(kaltofen_shoup_constant == constant);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_matches_the_reference_across_square_free_degrees) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> reference_context;
    polynomial_arithmetic::polynomial_context<backend_type> kaltofen_shoup_context;
    // Each polynomial is c*X^n - 3 with nonzero c and n below the field characteristic. Any common root of the
    // polynomial and its derivative would have to be zero, but zero is not a root, so every input is square-free.
    for (const std::size_t degree : {1U, 2U, 3U, 4U, 5U, 8U, 9U, 12U, 20U}) {
        BOOST_TEST_CONTEXT("degree = " << degree) {
            polynomial_type input(degree + 1, fq_value_type::zero());
            input.front() = -fq_value_type(3);
            input.back() = fq_value_type(degree + 1);
            check_reference_and_kaltofen_shoup_agree(input, reference_context, kaltofen_shoup_context);
        }
    }
}

BOOST_AUTO_TEST_CASE(reference_and_kaltofen_shoup_factorizations_support_native_extension_field_coefficients) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    schoolbook_backend backend;
    const polynomial_type first_linear_factor = {fq12_value(1), fq12_value_type::one()};
    const polynomial_type second_linear_factor = {fq12_value(20), fq12_value_type::one()};
    const polynomial_type input = multiply(backend, first_linear_factor, second_linear_factor);

    polynomial_arithmetic::polynomial_context<schoolbook_backend> reference_context;
    polynomial_arithmetic::polynomial_context<schoolbook_backend> kaltofen_shoup_context;
    const auto reference = math::distinct_degree_factorization_reference<schoolbook_backend>(input, reference_context);
    const auto kaltofen_shoup =
        math::distinct_degree_factorization_kaltofen_shoup<schoolbook_backend>(input, kaltofen_shoup_context);

    // Products of representatives modulo this degree-two polynomial contain at most three coefficients.
    polynomial_arithmetic::polynomial_context<mixed_radix_backend> mixed_radix_context {mixed_radix_backend(3)};
    const auto mixed_radix_kaltofen_shoup =
        math::distinct_degree_factorization_kaltofen_shoup<mixed_radix_backend>(input, mixed_radix_context);

    BOOST_CHECK(reference.complete);
    BOOST_REQUIRE_EQUAL(reference.factors.size(), 1);
    BOOST_CHECK_EQUAL(reference.factors.front().irreducible_factor_degree, 1);
    BOOST_CHECK(reference.factors.front().polynomial == input);
    BOOST_CHECK(kaltofen_shoup == reference);
    BOOST_CHECK(mixed_radix_kaltofen_shoup == reference);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_block_size_balances_baby_and_giant_steps) {
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(0), 1);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(1), 1);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(2), 1);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(8), 2);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(9), 3);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(18), 3);
    BOOST_CHECK_EQUAL(math::detail::kaltofen_shoup_block_size(19), 4);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_precomputation_builds_baby_and_giant_frobenius_steps) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    constexpr std::size_t block_size = 3;
    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type x = {fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(block_size, frobenius_context,
                                                                                       arithmetic_context);

    BOOST_CHECK_EQUAL(precomputation.block_size(), block_size);
    polynomial_type expected_baby_step(x);
    for (std::size_t index = 0; index <= block_size; ++index) {
        BOOST_CHECK(precomputation.baby_step(index) == expected_baby_step);
        if (index != block_size) {
            math::frobenius_map(expected_baby_step, expected_baby_step, frobenius_context, arithmetic_context);
        }
    }

    polynomial_type expected_giant_step(x);
    polynomial_type giant_step = precomputation.baby_step(block_size);
    for (std::size_t block = 1; block <= 3; ++block) {
        math::frobenius_map(expected_giant_step, expected_giant_step, block_size, frobenius_context,
                            arithmetic_context);
        BOOST_CHECK(giant_step == expected_giant_step);
        if (block != 3) {
            precomputation.apply_block_frobenius(giant_step, giant_step, frobenius_context.divisor_context(),
                                                 arithmetic_context);
        }
    }
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_giant_steps_outlive_the_initial_prepared_cache) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using precomputation_type = math::detail::kaltofen_shoup_frobenius_precomputation<backend_type>;
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_prepared_byte_limit = 1024 * 1024;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
    const polynomial_type modulus = {fq12_value(2), fq12_value(3), fq12_value(4), fq12_value_type::zero(),
                                     fq12_value_type::one()};
    std::optional<math::polynomial_divisor_context<backend_type>> divisor;
    std::vector<polynomial_type> babies;
    polynomial_type expected;
    {
        math::polynomial_frobenius_context<backend_type> initial(modulus, context);
        BOOST_REQUIRE_GT(initial.composition_precomputation().prepared_coefficient_bytes(), 0);
        babies = precomputation_type::make_baby_steps(3, initial, context, {});
        math::frobenius_map(expected, babies.back(), 3, initial, context);
        divisor.emplace(initial.divisor_context());
    }
    const precomputation_type precomputation(std::move(babies), *divisor, context);
    polynomial_type giant = precomputation.baby_step(3);
    precomputation.apply_block_frobenius(giant, giant, *divisor, context);
    BOOST_CHECK(giant == expected);
    BOOST_CHECK_THROW(precomputation_type(std::vector<polynomial_type> {}, *divisor, context), std::invalid_argument);
    BOOST_CHECK_THROW(precomputation_type(std::vector<polynomial_type> {{fq12_value_type::one()}}, *divisor, context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_precomputation_rejects_a_zero_block_size) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);

    BOOST_CHECK_THROW(
        math::detail::kaltofen_shoup_frobenius_precomputation<backend_type>(0, frobenius_context, arithmetic_context),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_precomputation_reduces_the_initial_baby_step) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type linear_divisor = {fq_value_type(5), fq_value_type::one()};
    math::polynomial_frobenius_context<backend_type> linear_frobenius_context(linear_divisor, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> linear_precomputation(
        1, linear_frobenius_context, arithmetic_context);
    BOOST_CHECK(linear_precomputation.baby_step(0) == polynomial_type({-fq_value_type(5)}));

    const polynomial_type constant_divisor = {fq_value_type(7)};
    math::polynomial_frobenius_context<backend_type> constant_frobenius_context(constant_divisor, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> constant_precomputation(
        1, constant_frobenius_context, arithmetic_context);
    BOOST_CHECK(constant_precomputation.baby_step(0) == polynomial_type({fq_value_type::zero()}));
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_coarse_blocks_extract_degree_intervals) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    // Three generates the multiplicative group of BN254 Fq, so 3 is not a square and neither 3 nor 9 = 3^2 is a
    // cube. Consequently X^2 - 3, X^3 - 3, and X^3 - 9 are irreducible over this field.
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type first_cubic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                                fq_value_type::one()};
    const polynomial_type second_cubic_factor = {-fq_value_type(9), fq_value_type::zero(), fq_value_type::zero(),
                                                 fq_value_type::one()};
    const polynomial_type first_block = multiply(backend, linear_factor, quadratic_factor);
    const polynomial_type second_block = multiply(backend, first_cubic_factor, second_cubic_factor);
    const polynomial_type input = multiply(backend, first_block, second_block);

    constexpr std::size_t block_size = 2;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(input, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(block_size, frobenius_context,
                                                                                       arithmetic_context);

    polynomial_type giant_step = precomputation.baby_step(block_size);
    polynomial_type factor;
    math::detail::kaltofen_shoup_coarse_block_factor(factor, input, giant_step, block_size, precomputation,
                                                     frobenius_context.divisor_context(), arithmetic_context);
    BOOST_CHECK(factor == first_block);

    precomputation.apply_block_frobenius(giant_step, giant_step, frobenius_context.divisor_context(),
                                         arithmetic_context);
    math::detail::kaltofen_shoup_coarse_block_factor(factor, second_block, giant_step, 1, precomputation,
                                                     frobenius_context.divisor_context(), arithmetic_context);
    BOOST_CHECK(factor == second_block);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_coarse_block_rejects_an_invalid_degree_count) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(divisor, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(2, frobenius_context,
                                                                                       arithmetic_context);
    const polynomial_type giant_step = precomputation.baby_step(2);
    polynomial_type output;

    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_coarse_block_factor(output, divisor, giant_step, 0, precomputation,
                                                                       frobenius_context.divisor_context(),
                                                                       arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_coarse_block_factor(output, divisor, giant_step, 3, precomputation,
                                                                       frobenius_context.divisor_context(),
                                                                       arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_fine_splitting_separates_exact_degrees) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    // Three generates the multiplicative group of BN254 Fq, so 3 is not a square and neither 3 nor 9 = 3^2 is a
    // cube. Consequently X^2 - 3, X^3 - 3, and X^3 - 9 are irreducible over this field.
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type first_cubic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                                fq_value_type::one()};
    const polynomial_type second_cubic_factor = {-fq_value_type(9), fq_value_type::zero(), fq_value_type::zero(),
                                                 fq_value_type::one()};
    const polynomial_type first_block = multiply(backend, linear_factor, quadratic_factor);
    const polynomial_type second_block = multiply(backend, first_cubic_factor, second_cubic_factor);
    const polynomial_type input = multiply(backend, first_block, second_block);

    constexpr std::size_t block_size = 2;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(input, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(block_size, frobenius_context,
                                                                                       arithmetic_context);

    std::vector<factor_type> factors;
    std::size_t callback_count = 0;
    auto record_factor = [&](const factor_type &) {
        ++callback_count;
        return math::factorization_control::continue_factorization;
    };

    polynomial_type giant_step = precomputation.baby_step(block_size);
    BOOST_CHECK(math::detail::kaltofen_shoup_split_coarse_block(&factors, first_block, giant_step, 1, block_size,
                                                                precomputation, arithmetic_context, record_factor) ==
                math::factorization_control::continue_factorization);

    precomputation.apply_block_frobenius(giant_step, giant_step, frobenius_context.divisor_context(),
                                         arithmetic_context);
    BOOST_CHECK(math::detail::kaltofen_shoup_split_coarse_block(&factors, second_block, giant_step, 3, 1,
                                                                precomputation, arithmetic_context, record_factor) ==
                math::factorization_control::continue_factorization);

    BOOST_CHECK_EQUAL(callback_count, 3);
    BOOST_REQUIRE_EQUAL(factors.size(), 3);
    BOOST_CHECK(factors[0] == factor_type({linear_factor, 1}));
    BOOST_CHECK(factors[1] == factor_type({quadratic_factor, 2}));
    BOOST_CHECK(factors[2] == factor_type({second_block, 3}));
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_fine_splitting_honors_an_early_stop) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type input = multiply(backend, linear_factor, quadratic_factor);

    constexpr std::size_t block_size = 2;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(input, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(block_size, frobenius_context,
                                                                                       arithmetic_context);
    const polynomial_type giant_step = precomputation.baby_step(block_size);
    std::vector<factor_type> factors;

    const auto control = math::detail::kaltofen_shoup_split_coarse_block(
        &factors, input, giant_step, 1, block_size, precomputation, arithmetic_context,
        [](const factor_type &) { return math::factorization_control::stop_factorization; });

    BOOST_CHECK(control == math::factorization_control::stop_factorization);
    BOOST_REQUIRE_EQUAL(factors.size(), 1);
    BOOST_CHECK(factors.front() == factor_type({linear_factor, 1}));
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_fine_splitting_rejects_invalid_degree_intervals) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type input = multiply(backend, linear_factor, quadratic_factor);

    constexpr std::size_t block_size = 2;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_frobenius_context<backend_type> frobenius_context(input, arithmetic_context);
    math::detail::kaltofen_shoup_frobenius_precomputation<backend_type> precomputation(block_size, frobenius_context,
                                                                                       arithmetic_context);
    const polynomial_type giant_step = precomputation.baby_step(block_size);
    std::vector<factor_type> factors;
    auto continue_factorization = [](const factor_type &) {
        return math::factorization_control::continue_factorization;
    };

    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_split_coarse_block(&factors, input, giant_step, 0, 1, precomputation,
                                                                      arithmetic_context, continue_factorization),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_split_coarse_block(&factors, input, giant_step, 1, 0, precomputation,
                                                                      arithmetic_context, continue_factorization),
                      std::invalid_argument);
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_split_coarse_block(&factors, input, giant_step, 1, 3, precomputation,
                                                                      arithmetic_context, continue_factorization),
                      std::invalid_argument);

    // The stated interval contains only degree one, so the unclassified quadratic must be rejected.
    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_split_coarse_block(&factors, input, giant_step, 1, 1, precomputation,
                                                                      arithmetic_context, continue_factorization),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_internal_driver_processes_multiple_blocks) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type first_cubic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                                fq_value_type::one()};
    const polynomial_type second_cubic_factor = {-fq_value_type(9), fq_value_type::zero(), fq_value_type::zero(),
                                                 fq_value_type::one()};
    const polynomial_type degree_three_factors = multiply(backend, first_cubic_factor, second_cubic_factor);
    polynomial_type input = multiply(backend, linear_factor, quadratic_factor);
    input = multiply(backend, input, degree_three_factors);

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::vector<factor_type> factors;
    std::size_t callback_count = 0;
    // A block size of two forces a full degree-1-to-2 block followed by a partial block containing degree three.
    const auto control = math::detail::kaltofen_shoup_factor_monic_square_free(
        &factors, input, 2, arithmetic_context, [&](const factor_type &) {
            ++callback_count;
            return math::factorization_control::continue_factorization;
        });

    BOOST_CHECK(control == math::factorization_control::continue_factorization);
    BOOST_CHECK_EQUAL(callback_count, 3);
    BOOST_REQUIRE_EQUAL(factors.size(), 3);
    BOOST_CHECK(factors[0] == factor_type({linear_factor, 1}));
    BOOST_CHECK(factors[1] == factor_type({quadratic_factor, 2}));
    BOOST_CHECK(factors[2] == factor_type({degree_three_factors, 3}));
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_internal_driver_emits_the_final_irreducible_factor) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    const polynomial_type cubic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                          fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::vector<factor_type> factors;
    const auto control = math::detail::kaltofen_shoup_factor_monic_square_free(
        &factors, cubic_factor, 2, arithmetic_context,
        [](const factor_type &) { return math::factorization_control::continue_factorization; });

    BOOST_CHECK(control == math::factorization_control::continue_factorization);
    BOOST_REQUIRE_EQUAL(factors.size(), 1);
    BOOST_CHECK(factors.front() == factor_type({cubic_factor, 3}));
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_internal_driver_rejects_a_zero_block_size) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    const polynomial_type input = {fq_value_type::one(), fq_value_type::one()};
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::vector<factor_type> factors;

    BOOST_CHECK_THROW(math::detail::kaltofen_shoup_factor_monic_square_free(
                          &factors, input, 0, arithmetic_context,
                          [](const factor_type &) { return math::factorization_control::continue_factorization; }),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_internal_driver_propagates_an_early_stop) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    backend_type backend;
    const polynomial_type linear_factor = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type quadratic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type cubic_factor = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                          fq_value_type::one()};
    polynomial_type input = multiply(backend, linear_factor, quadratic_factor);
    input = multiply(backend, input, cubic_factor);

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    std::vector<factor_type> factors;
    const auto control = math::detail::kaltofen_shoup_factor_monic_square_free(
        &factors, input, 2, arithmetic_context, [](const factor_type &factor) {
            return factor.irreducible_factor_degree == 2 ? math::factorization_control::stop_factorization :
                                                           math::factorization_control::continue_factorization;
        });

    BOOST_CHECK(control == math::factorization_control::stop_factorization);
    BOOST_REQUIRE_EQUAL(factors.size(), 2);
    BOOST_CHECK(factors[0] == factor_type({linear_factor, 1}));
    BOOST_CHECK(factors[1] == factor_type({quadratic_factor, 2}));
}

BOOST_AUTO_TEST_CASE(kaltofen_shoup_streaming_matches_collected_groups_and_stop_prefixes) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using factor_type = math::distinct_degree_factor<polynomial_type>;

    backend_type backend;
    const polynomial_type linear = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type quadratic = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type cubic = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                   fq_value_type::one()};
    const polynomial_type second_cubic = {-fq_value_type(9), fq_value_type::zero(), fq_value_type::zero(),
                                          fq_value_type::one()};
    const auto mixed = multiply(backend, multiply(backend, linear, quadratic), multiply(backend, cubic, second_cubic));
    for (const auto &input : {linear, cubic, mixed}) {
        for (std::size_t screening_limit : {0, 1, 100}) {
            for (std::size_t block_size : {1, 2}) {
                polynomial_arithmetic::polynomial_context_options options;
                options.distinct_degree_screening_limit = screening_limit;
                polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
                std::vector<factor_type> expected;
                BOOST_REQUIRE(math::detail::kaltofen_shoup_factor_monic_square_free(
                                  &expected, input, block_size, context, [&](const factor_type &factor) {
                                      // The public collecting path must append before reporting the group.
                                      BOOST_REQUIRE(!expected.empty());
                                      BOOST_CHECK(&factor == &expected.back());
                                      return math::factorization_control::continue_factorization;
                                  }) == math::factorization_control::continue_factorization);

                for (std::size_t stop_after = 0; stop_after <= expected.size(); ++stop_after) {
                    polynomial_arithmetic::polynomial_context<backend_type> streaming_context(backend_type {}, options);
                    std::size_t reported = 0;
                    const auto control = math::detail::kaltofen_shoup_factor_monic_square_free(
                        nullptr, input, block_size, streaming_context, [&](const factor_type &factor) {
                            BOOST_REQUIRE_LT(reported, expected.size());
                            BOOST_CHECK(factor == expected[reported++]);
                            return reported == stop_after ? math::factorization_control::stop_factorization :
                                                            math::factorization_control::continue_factorization;
                        });
                    BOOST_CHECK_EQUAL(reported, stop_after ? stop_after : expected.size());
                    BOOST_CHECK(control == (stop_after ? math::factorization_control::stop_factorization :
                                                         math::factorization_control::continue_factorization));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(screening_continuation_matches_reference_and_unscreened_blocks) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using mixed_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    backend_type backend;
    const polynomial_type linear = {fq_value_type::one(), fq_value_type::one()};
    const polynomial_type quadratic = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type cubic = {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(),
                                   fq_value_type::one()};
    const polynomial_type second_cubic = {-fq_value_type(9), fq_value_type::zero(), fq_value_type::zero(),
                                          fq_value_type::one()};
    const auto higher_degrees = multiply(backend, quadratic, multiply(backend, cubic, second_cubic));
    const auto mixed_degrees = multiply(backend, linear, higher_degrees);
    auto nonmonic = mixed_degrees;
    math::scalar_multiplication(nonmonic, nonmonic, fq_value_type(7));
    const std::vector<polynomial_type> inputs = {
        {fq_value_type::zero()},
        {fq_value_type(7)},
        linear,
        quadratic,
        cubic,
        multiply(backend, linear, polynomial_type {fq_value_type::zero(), fq_value_type::one()}),
        higher_degrees,
        mixed_degrees,
        nonmonic};
    for (const auto &input : inputs) {
        polynomial_arithmetic::polynomial_context<backend_type> reference_context;
        const auto reference = math::distinct_degree_factorization_reference<backend_type>(input, reference_context);
        for (std::size_t limit : {0, 1, 2, 100}) {
            polynomial_arithmetic::polynomial_context_options options;
            options.distinct_degree_screening_limit = limit;
            options.basecase_divisor_coefficient_cutoff = 0;
            options.basecase_quotient_coefficient_cutoff = 0;
            polynomial_arithmetic::polynomial_context<backend_type> school(backend_type {}, options);
            polynomial_arithmetic::polynomial_context<mixed_backend> mixed(mixed_backend(26), options);
            BOOST_CHECK(math::distinct_degree_factorization_kaltofen_shoup<backend_type>(input, school) == reference);
            BOOST_CHECK(math::distinct_degree_factorization_kaltofen_shoup<mixed_backend>(input, mixed) == reference);
        }
    }
    // Smaller explicit block sizes force continuation across several degree blocks after the linear screen.
    for (std::size_t block_size : {1, 2}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.distinct_degree_screening_limit = 1;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        const auto reference = math::distinct_degree_factorization_reference<backend_type>(mixed_degrees, context);
        std::vector<math::distinct_degree_factor<polynomial_type>> groups;
        const auto control = math::detail::kaltofen_shoup_factor_monic_square_free(
            &groups, mixed_degrees, block_size, context,
            [](const auto &) { return math::factorization_control::continue_factorization; });
        BOOST_CHECK(control == math::factorization_control::continue_factorization);
        BOOST_CHECK(groups == reference.factors);
    }
}

BOOST_AUTO_TEST_CASE(screening_preserves_full_fq12_groups_and_early_stops) {
    using school_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename school_backend::polynomial_type;
    school_backend backend;
    const polynomial_type first = {-fq12_value(1), fq12_value_type::one()};
    const polynomial_type second = {-fq12_value(20), fq12_value_type::one()};
    const auto input = multiply(backend, first, second);
    polynomial_arithmetic::polynomial_context_options options;
    options.distinct_degree_screening_limit = 1;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<school_backend> school(school_backend {}, options);
    polynomial_arithmetic::polynomial_context<mixed_backend> mixed(mixed_backend(18), options);
    std::size_t callbacks = 0;
    const auto callback = [&](const auto &group) {
        ++callbacks;
        // Screening emits a reducible degree-one group; irreducible splitting remains the next stage's job.
        BOOST_CHECK_EQUAL(group.irreducible_factor_degree, 1);
        BOOST_CHECK(group.polynomial == input);
        return math::factorization_control::stop_factorization;
    };
    const auto reference = math::distinct_degree_factorization_kaltofen_shoup<school_backend>(input, school, callback);
    const auto result = math::distinct_degree_factorization_kaltofen_shoup<mixed_backend>(input, mixed, callback);
    BOOST_CHECK(!result.complete);
    BOOST_CHECK(result == reference);
    BOOST_CHECK_EQUAL(callbacks, 2);
}

BOOST_AUTO_TEST_CASE(parity_bound_is_strict_and_does_not_overflow) {
    for (std::size_t next : {1, 2, 87, 4873}) {
        BOOST_CHECK(math::detail::kaltofen_shoup_requires_degree_test(3 * next, next, true));
        BOOST_CHECK(!math::detail::kaltofen_shoup_requires_degree_test(3 * next - 1, next, true));
        BOOST_CHECK(math::detail::kaltofen_shoup_requires_degree_test(3 * next - 1, next, false));
        BOOST_CHECK(math::detail::kaltofen_shoup_requires_degree_test(3 * next - 1, next, std::nullopt));
        BOOST_CHECK(!math::detail::kaltofen_shoup_requires_degree_test(2 * next - 1, next, false));
    }
    const auto maximum = std::numeric_limits<std::size_t>::max();
    BOOST_CHECK(!math::detail::kaltofen_shoup_requires_degree_test(maximum, maximum / 3 + 1, true));
    BOOST_CHECK(math::detail::kaltofen_shoup_requires_degree_test(maximum, maximum / 3, true));
    BOOST_CHECK(!math::detail::kaltofen_shoup_requires_degree_test(0, 1, true));
}

BOOST_AUTO_TEST_CASE(parity_tracks_factor_counts_and_preserves_blocks_and_stop_prefixes) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = backend_type::polynomial_type;
    using group_type = math::distinct_degree_factor<polynomial_type>;
    backend_type backend;
    const polynomial_type linear {fq_value_type(1), fq_value_type::one()};
    const polynomial_type second_linear {fq_value_type(2), fq_value_type::one()};
    const polynomial_type quadratic {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type cubic {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(), fq_value_type::one()};
    const auto large = dense_irreducible();
    const std::vector<polynomial_type> inputs {
        large, multiply(backend, linear, large), multiply(backend, multiply(backend, linear, second_linear), large),
        multiply(backend, quadratic, cubic), multiply(backend, multiply(backend, linear, quadratic), cubic)};
    const std::vector<std::size_t> counts {1, 2, 3, 2, 3};
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        polynomial_arithmetic::polynomial_context<backend_type> reference_context;
        const auto expected =
            math::distinct_degree_factorization_reference<backend_type>(inputs[index], reference_context);
        for (std::size_t block : {1, 2, 3, 4}) {
            for (std::size_t screening : {0, 1, 4}) {
                for (std::size_t certificate_threshold : {0, 1}) {
                    polynomial_arithmetic::polynomial_context_options options;
                    options.distinct_degree_screening_limit = screening;
                    options.distinct_degree_certificate_minimum_excluded_degree = certificate_threshold;
                    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
                    for (std::size_t stop = 0; stop <= expected.factors.size(); ++stop) {
                        std::vector<group_type> collected;
                        const auto result = math::detail::kaltofen_shoup_factor_monic_square_free(
                            &collected, inputs[index], block, context,
                            [&](const group_type &group) {
                                BOOST_CHECK(&group == &collected.back());
                                return collected.size() == stop ? math::factorization_control::stop_factorization :
                                                                  math::factorization_control::continue_factorization;
                            },
                            counts[index] % 2 != 0);
                        const auto size = stop ? stop : expected.factors.size();
                        BOOST_REQUIRE_EQUAL(collected.size(), size);
                        BOOST_CHECK(std::equal(collected.begin(), collected.end(), expected.factors.begin()));
                        BOOST_CHECK(result == (stop ? math::factorization_control::stop_factorization :
                                                      math::factorization_control::continue_factorization));
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(parity_removes_polynomial_work_and_even_composites_keep_the_same_work) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = backend_type::polynomial_type;
    backend_type backend;
    const polynomial_type linear {fq_value_type(1), fq_value_type::one()};
    const polynomial_type second_linear {fq_value_type(2), fq_value_type::one()};
    const polynomial_type quadratic {-fq_value_type(3), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type cubic {-fq_value_type(3), fq_value_type::zero(), fq_value_type::zero(), fq_value_type::one()};
    const auto large = dense_irreducible();
    const std::vector<polynomial_type> inputs {large, multiply(backend, linear, large),
                                               multiply(backend, multiply(backend, linear, second_linear), large),
                                               multiply(backend, quadratic, cubic)};
    const std::vector<bool> parity {true, false, true, false};
    for (std::size_t i = 0; i < inputs.size(); ++i) {
        std::size_t baseline_calls = 0, candidate_calls = 0;
        polynomial_arithmetic::polynomial_context_options options;
        options.distinct_degree_screening_limit = 1;
        polynomial_arithmetic::polynomial_context<counting_backend> baseline(counting_backend(baseline_calls), options);
        polynomial_arithmetic::polynomial_context<counting_backend> candidate(counting_backend(candidate_calls),
                                                                              options);
        std::vector<math::distinct_degree_factor<polynomial_type>> before, after;
        const auto proceed = [](const auto &) { return math::factorization_control::continue_factorization; };
        math::detail::kaltofen_shoup_factor_monic_square_free(&before, inputs[i], 1, baseline, proceed);
        math::detail::kaltofen_shoup_factor_monic_square_free(&after, inputs[i], 1, candidate, proceed, parity[i]);
        BOOST_CHECK(before == after);
        if (i < 3) {
            BOOST_CHECK_LT(candidate_calls, baseline_calls);
        } else {
            BOOST_CHECK_EQUAL(candidate_calls, baseline_calls);
            options.distinct_degree_certificate_minimum_excluded_degree = 1;
            std::size_t guarded_calls = 0;
            polynomial_arithmetic::polynomial_context<counting_backend> guarded(counting_backend(guarded_calls),
                                                                                options);
            std::vector<math::distinct_degree_factor<polynomial_type>> with_certificate;
            math::detail::kaltofen_shoup_factor_monic_square_free(&with_certificate, inputs[i], 1, guarded, proceed,
                                                                  parity[i]);
            BOOST_CHECK(with_certificate == before);
            BOOST_CHECK_EQUAL(guarded_calls, candidate_calls);
        }
    }
}

BOOST_AUTO_TEST_CASE(public_discovery_uses_the_gcd_certificate_and_preserves_unknown_fallback) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = backend_type::polynomial_type;
    const auto input = dense_irreducible();
    polynomial_arithmetic::polynomial_context<backend_type> context;
    polynomial_type monic;
    fq_value_type leading;
    std::optional<bool> parity;
    BOOST_REQUIRE(
        math::detail::prepare_square_free_factorization_input<backend_type>(monic, leading, input, context, &parity));
    BOOST_REQUIRE(parity.has_value());
    BOOST_CHECK(*parity);
    check_reference_and_kaltofen_shoup_agree(input, context, context);
    const polynomial_type irregular {fq_value_type::one(), fq_value_type::zero(), fq_value_type::zero(),
                                     fq_value_type::one()};
    BOOST_REQUIRE(math::detail::prepare_square_free_factorization_input<backend_type>(monic, leading, irregular,
                                                                                      context, &parity));
    BOOST_CHECK(!parity.has_value());
    check_reference_and_kaltofen_shoup_agree(irregular, context, context);
}

BOOST_AUTO_TEST_SUITE_END()
