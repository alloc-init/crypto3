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

#define BOOST_TEST_MODULE polynomial_composition_test

#include <atomic>
#include <cstddef>
#include <exception>
#include <future>
#include <limits>
#include <optional>
#include <thread>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_composition.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/bn254_small_prime_fft_backend.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq_value_type = fq_field_type::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;

    fq12_value_type fq12_value(std::size_t first_coordinate) {
        fq12_value_type value = fq12_value_type::zero();
        for (std::size_t i = 0; i < fq12_field_type::arity; ++i) {
            value.coordinate(i) = fq_value_type(first_coordinate + i);
        }
        return value;
    }

    struct preparation_execution_state {
        std::size_t threads = 4;
        std::size_t calls = 0;
        bool fail = false;
        std::atomic<std::size_t> active_workers = 0;
    };

    template<typename Preparation>
    void run_preparation_panels(const Preparation &preparation, preparation_execution_state &state) {
        ++state.calls;
        const std::size_t count = preparation.panel_count();
        BOOST_REQUIRE_GT(count, 0);
        BOOST_CHECK_THROW(preparation.prepare_panels(1, 0), std::out_of_range);
        BOOST_CHECK_THROW(preparation.prepare_panels(0, count + 1), std::out_of_range);
        preparation.prepare_panels(count, count);
        const std::size_t workers = std::min(state.threads, count);
        std::vector<std::atomic<std::size_t>> visits(count);
        std::vector<std::exception_ptr> errors(workers);
        std::vector<std::jthread> threads;
        for (std::size_t rank = 0; rank < workers; ++rank) {
            threads.emplace_back([&, rank] {
                ++state.active_workers;
                try {
                    // Interleave complete panels in reverse order; no two workers write the same panel.
                    for (std::size_t offset = rank; offset < count; offset += workers) {
                        const std::size_t panel = count - 1 - offset;
                        preparation.prepare_panels(panel, panel + 1);
                        ++visits[panel];
                        if (state.fail && rank == 0) {
                            throw std::runtime_error("application preparation failure");
                        }
                    }
                } catch (...) {
                    errors[rank] = std::current_exception();
                }
                --state.active_workers;
            });
        }
        for (auto &thread : threads) {
            thread.join();
        }
        for (const auto &error : errors) {
            if (error) {
                std::rethrow_exception(error);
            }
        }
        for (const auto &visited : visits) {
            BOOST_CHECK_EQUAL(visited.load(), 1);
        }
    }

    template<typename Value>
    struct preparation_backend : polynomial_arithmetic::schoolbook_backend<Value> {
        explicit preparation_backend(preparation_execution_state &state) : state(state) {
        }

        template<typename Preparation>
        void prepare_composition_coefficients(const Preparation &preparation) {
            run_preparation_panels(preparation, state);
        }

        preparation_execution_state &state;
    };

    using grouped_backend_type = polynomial_arithmetic::bn254_small_prime_fft_backend;
    using grouped_polynomial_type = grouped_backend_type::polynomial_type;
    static_assert(math::detail::SupportsPreparedProductSum<grouped_backend_type>);
    static_assert(
        !math::detail::SupportsPreparedProductSum<polynomial_arithmetic::schoolbook_backend<fq12_value_type>>);
    static_assert(!math::detail::SupportsPreparedProductSum<
                  polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>>);

    grouped_polynomial_type grouped_polynomial(std::size_t size, std::size_t seed) {
        grouped_polynomial_type result(size, fq12_value_type::zero());
        for (std::size_t i = 0; i < size; ++i) {
            const auto value = fq12_value(seed + 13 * i);
            result[i] = i % 3 == 1 ? -value : value;
        }
        return result;
    }

    struct grouped_operation_counts {
        std::size_t append_attempts = 0;
        std::size_t products = 0;
        std::size_t finishes = 0;
        std::size_t remainders = 0;
        std::size_t decline_on_append = 0;
        bool decline_finish = false;
        bool decline_preparation = false;
    };

    struct grouped_counting_backend : grouped_backend_type {
        grouped_counting_backend(std::size_t size, grouped_operation_counts &counts) :
            grouped_backend_type(size), counts(counts) {
        }

        std::optional<product_accumulator_type>
            make_product_accumulator(std::size_t count, std::size_t products,
                                     std::size_t bytes = std::numeric_limits<std::size_t>::max()) const {
            if (counts.decline_preparation) {
                return std::nullopt;
            }
            return grouped_backend_type::make_product_accumulator(count, products, bytes);
        }

        bool try_accumulate_prepared_product(product_accumulator_type &accumulator, const polynomial_type &input,
                                             const prepared_low_product_type &prepared) {
            if (++counts.append_attempts == counts.decline_on_append) {
                return false;
            }
            const bool used = grouped_backend_type::try_accumulate_prepared_product(accumulator, input, prepared);
            counts.products += used;
            return used;
        }

        bool try_finish_product_sum(polynomial_type &output, product_accumulator_type &accumulator) {
            if (counts.decline_finish) {
                return false;
            }
            const bool used = grouped_backend_type::try_finish_product_sum(output, accumulator);
            counts.finishes += used;
            return used;
        }

        bool try_cyclic_remainder(polynomial_type &output, const polynomial_type &dividend,
                                  const polynomial_type &quotient, const prepared_cyclic_remainder_type &prepared) {
            const bool used = grouped_backend_type::try_cyclic_remainder(output, dividend, quotient, prepared);
            counts.remainders += used;
            return used;
        }

        grouped_operation_counts &counts;
    };

    struct composition_execution_counts {
        std::size_t batches = 0;
        bool fail = false;
    };

    struct partitioned_composition_backend : grouped_backend_type {
        partitioned_composition_backend(std::size_t size, composition_execution_counts &counts,
                                        std::size_t threads = 1) :
            grouped_backend_type(size), counts(counts), threads(threads) {
        }

        template<typename Evaluation>
        void evaluate_composition_coefficients(const Evaluation &evaluation) {
            ++counts.batches;
            // Application scheduling stays in this test. Each worker visits a disjoint reverse stride;
            // the calling thread participates, and every helper completes before return or rethrow.
            const auto visit = [&](std::size_t rank) {
                for (std::size_t offset = rank; offset < evaluation.panel_count(); offset += threads) {
                    const auto panel = evaluation.panel_count() - 1 - offset;
                    evaluation.evaluate_panels(panel, panel + 1);
                    if (counts.fail) {
                        throw std::runtime_error("application coefficient evaluation failed");
                    }
                }
            };
            std::vector<std::future<void>> helpers;
            helpers.reserve(threads - 1);
            for (std::size_t rank = 1; rank < threads; ++rank) {
                helpers.push_back(std::async(std::launch::async, [&, rank] { visit(rank); }));
            }
            std::exception_ptr failure;
            try {
                visit(0);
            } catch (...) {
                failure = std::current_exception();
            }
            for (auto &helper : helpers) {
                try {
                    helper.get();
                } catch (...) {
                    if (!failure) {
                        failure = std::current_exception();
                    }
                }
            }
            if (failure) {
                std::rethrow_exception(failure);
            }
        }

        composition_execution_counts &counts;
        std::size_t threads;
    };
}    // namespace

BOOST_AUTO_TEST_SUITE(polynomial_composition_test_suite)

BOOST_AUTO_TEST_CASE(reference_composition_matches_a_known_result_and_may_alias_inputs) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    // Modulo X^2 + 1, (1 + 2Y + 3Y^2) evaluated at Y = X is -2 + 2X.
    const polynomial_type outer = {fq_value_type(1), fq_value_type(2), fq_value_type(3)};
    const polynomial_type inner = {fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type expected = {-fq_value_type(2), fq_value_type(2)};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);

    polynomial_type result;
    math::compose_mod_reference(result, outer, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(result == expected);

    polynomial_type outer_alias = outer;
    math::compose_mod_reference(outer_alias, outer_alias, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(outer_alias == expected);

    polynomial_type inner_alias = inner;
    math::compose_mod_reference(inner_alias, outer, inner_alias, divisor_context, arithmetic_context);
    BOOST_CHECK(inner_alias == expected);
}

BOOST_AUTO_TEST_CASE(reference_composition_reduces_the_inner_polynomial) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type outer = {fq_value_type(2), fq_value_type(3), fq_value_type(4)};
    const polynomial_type inner = {fq_value_type(5), fq_value_type(6), fq_value_type(7), fq_value_type(8)};
    const polynomial_type divisor = {fq_value_type(1), fq_value_type(1), fq_value_type(1)};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 2, arithmetic_context);

    polynomial_type reduced_inner;
    math::remainder(reduced_inner, inner, divisor_context, arithmetic_context);

    polynomial_type from_unreduced;
    polynomial_type from_reduced;
    math::compose_mod_reference(from_unreduced, outer, inner, divisor_context, arithmetic_context);
    math::compose_mod_reference(from_reduced, outer, reduced_inner, divisor_context, arithmetic_context);
    BOOST_CHECK(from_unreduced == from_reduced);

    polynomial_type inner_alias = inner;
    math::compose_mod_reference(inner_alias, outer, inner_alias, divisor_context, arithmetic_context);
    BOOST_CHECK(inner_alias == from_unreduced);
}

BOOST_AUTO_TEST_CASE(composition_handles_zero_constant_and_constant_modulus_cases) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type(1)};
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);
    const polynomial_type inner = {fq_value_type(4), fq_value_type(5), fq_value_type(6), fq_value_type(7)};

    polynomial_type result;
    math::compose_mod_reference(result, polynomial_type {fq_value_type::zero()}, inner, divisor_context,
                                arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));
    math::compose_mod(result, polynomial_type {fq_value_type::zero()}, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));

    math::compose_mod_reference(result, polynomial_type {fq_value_type(7)}, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type(7)}));
    math::compose_mod(result, polynomial_type {fq_value_type(7)}, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type(7)}));

    const polynomial_type constant_divisor = {fq_value_type(11)};
    math::polynomial_divisor_context<backend_type> constant_context(constant_divisor, 1, arithmetic_context);
    math::compose_mod_reference(result, polynomial_type {fq_value_type(1), fq_value_type(2)}, inner, constant_context,
                                arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));
    math::compose_mod(result, polynomial_type {fq_value_type(1), fq_value_type(2)}, inner, constant_context,
                      arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));
}

BOOST_AUTO_TEST_CASE(composition_rejects_insufficient_inverse_precision_without_changing_output) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type {}, options);

    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);
    const polynomial_type outer = {fq_value_type(2), fq_value_type(3)};
    const polynomial_type inner = {fq_value_type(4), fq_value_type(5), fq_value_type(6), fq_value_type(7)};

    polynomial_type output = {fq_value_type(17)};
    BOOST_CHECK_THROW(math::compose_mod_reference(output, outer, inner, divisor_context, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK(output == polynomial_type({fq_value_type(17)}));

    BOOST_CHECK_THROW(math::compose_mod(output, outer, inner, divisor_context, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK(output == polynomial_type({fq_value_type(17)}));
}

BOOST_AUTO_TEST_CASE(brent_kung_composition_matches_the_reference_at_different_memory_limits) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type outer = {fq_value_type(1), fq_value_type(2), fq_value_type(3), fq_value_type(4),
                                   fq_value_type(5), fq_value_type(6), fq_value_type(7), fq_value_type(8),
                                   fq_value_type(9), fq_value_type(10)};
    const polynomial_type inner = {fq_value_type(11), fq_value_type(12), fq_value_type(13),
                                   fq_value_type(14), fq_value_type(15), fq_value_type(16)};
    const polynomial_type divisor = {fq_value_type(17), fq_value_type(18), fq_value_type(19),   fq_value_type(20),
                                     fq_value_type(21), fq_value_type(22), fq_value_type::one()};

    polynomial_arithmetic::polynomial_context_options reference_options;
    reference_options.basecase_divisor_coefficient_cutoff = 0;
    reference_options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> reference_context(backend_type {}, reference_options);
    math::polynomial_divisor_context<backend_type> reference_divisor(divisor, 5, reference_context);
    polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, reference_divisor, reference_context);

    for (const std::size_t cached_power_limit : {1, 2, 3, 4, 100}) {
        polynomial_arithmetic::polynomial_context_options options = reference_options;
        options.modular_composition_cached_power_limit = cached_power_limit;
        polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> divisor_context(divisor, 5, arithmetic_context);

        polynomial_type result;
        math::compose_mod(result, outer, inner, divisor_context, arithmetic_context);
        BOOST_CHECK(result == expected);
    }
}

BOOST_AUTO_TEST_CASE(brent_kung_precomputation_is_reusable_and_enforces_its_outer_size) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type first_outer = {fq_value_type(1), fq_value_type(2), fq_value_type(3), fq_value_type(4),
                                         fq_value_type(5), fq_value_type(6), fq_value_type(7), fq_value_type(8),
                                         fq_value_type(9), fq_value_type(10)};
    const polynomial_type second_outer = {fq_value_type(11), fq_value_type(12), fq_value_type(13),
                                          fq_value_type(14), fq_value_type(15), fq_value_type(16)};
    const polynomial_type inner = {fq_value_type(17), fq_value_type(18), fq_value_type(19),
                                   fq_value_type(20), fq_value_type(21), fq_value_type(22)};
    const polynomial_type divisor = {fq_value_type(23), fq_value_type(24), fq_value_type(25),   fq_value_type(26),
                                     fq_value_type(27), fq_value_type(28), fq_value_type::one()};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 5, arithmetic_context);
    math::polynomial_composition_precomputation<backend_type> precomputation(inner, first_outer.size(), divisor_context,
                                                                             arithmetic_context);

    for (const polynomial_type &outer : {first_outer, second_outer}) {
        polynomial_type expected;
        math::compose_mod_reference(expected, outer, inner, divisor_context, arithmetic_context);

        polynomial_type result;
        math::compose_mod(result, outer, precomputation, divisor_context, arithmetic_context);
        BOOST_CHECK(result == expected);

        polynomial_type aliased = outer;
        math::compose_mod(aliased, aliased, precomputation, divisor_context, arithmetic_context);
        BOOST_CHECK(aliased == expected);
    }

    polynomial_type oversized_outer = first_outer;
    oversized_outer.push_back(fq_value_type(11));
    polynomial_type unchanged = {fq_value_type(17)};
    BOOST_CHECK_THROW(
        math::compose_mod(unchanged, oversized_outer, precomputation, divisor_context, arithmetic_context),
        std::invalid_argument);
    BOOST_CHECK(unchanged == polynomial_type({fq_value_type(17)}));

    const polynomial_type constant_divisor = {fq_value_type(29)};
    math::polynomial_divisor_context<backend_type> constant_divisor_context(constant_divisor, 1, arithmetic_context);
    math::polynomial_composition_precomputation<backend_type> constant_precomputation(
        inner, first_outer.size(), constant_divisor_context, arithmetic_context);
    BOOST_CHECK_THROW(math::compose_mod(unchanged, oversized_outer, constant_precomputation, constant_divisor_context,
                                        arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK(unchanged == polynomial_type({fq_value_type(17)}));

    BOOST_CHECK_THROW(
        math::polynomial_composition_precomputation<backend_type>(inner, 0, divisor_context, arithmetic_context),
        std::invalid_argument);

    polynomial_arithmetic::polynomial_context_options invalid_options;
    invalid_options.modular_composition_cached_power_limit = 0;
    polynomial_arithmetic::polynomial_context<backend_type> invalid_arithmetic_context(backend_type {},
                                                                                       invalid_options);
    math::polynomial_divisor_context<backend_type> invalid_divisor_context(divisor, 5, invalid_arithmetic_context);
    BOOST_CHECK_THROW(math::polynomial_composition_precomputation<backend_type>(
                          inner, first_outer.size(), invalid_divisor_context, invalid_arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(brent_kung_composition_may_alias_inputs_and_rejects_a_zero_cache_limit) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type outer = {fq_value_type(1), fq_value_type(2), fq_value_type(3), fq_value_type(4)};
    const polynomial_type inner = {fq_value_type(5), fq_value_type(6), fq_value_type(7), fq_value_type(8)};
    const polynomial_type divisor = {fq_value_type(9), fq_value_type(10), fq_value_type(11), fq_value_type::one()};

    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_cached_power_limit = 2;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type {}, options);
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 2, arithmetic_context);
    polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, divisor_context, arithmetic_context);

    polynomial_type outer_alias = outer;
    math::compose_mod(outer_alias, outer_alias, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(outer_alias == expected);

    polynomial_type inner_alias = inner;
    math::compose_mod(inner_alias, outer, inner_alias, divisor_context, arithmetic_context);
    BOOST_CHECK(inner_alias == expected);

    options.modular_composition_cached_power_limit = 0;
    polynomial_arithmetic::polynomial_context<backend_type> invalid_context(backend_type {}, options);
    math::polynomial_divisor_context<backend_type> invalid_divisor(divisor, 2, invalid_context);
    polynomial_type output = {fq_value_type(17)};
    BOOST_CHECK_THROW(math::compose_mod(output, outer, inner, invalid_divisor, invalid_context), std::invalid_argument);
    BOOST_CHECK(output == polynomial_type({fq_value_type(17)}));
}

BOOST_AUTO_TEST_CASE(composition_matches_between_schoolbook_and_mixed_radix_backends) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    const polynomial_type outer = {fq12_value(1), fq12_value(13), fq12_value(25), fq12_value(37)};
    const polynomial_type inner = {fq12_value(49), fq12_value(61), fq12_value(73)};
    const polynomial_type divisor = {fq12_value(85), fq12_value(97), fq12_value(109), fq12_value(121)};

    polynomial_arithmetic::polynomial_context<schoolbook_backend> schoolbook_context;
    math::polynomial_divisor_context<schoolbook_backend> schoolbook_divisor(divisor, 2, schoolbook_context);
    polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, schoolbook_divisor, schoolbook_context);

    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    options.modular_composition_cached_power_limit = 2;
    polynomial_arithmetic::polynomial_context<mixed_radix_backend> mixed_radix_context(mixed_radix_backend(18),
                                                                                       options);
    math::polynomial_divisor_context<mixed_radix_backend> mixed_radix_divisor(divisor, 2, mixed_radix_context);
    polynomial_type result;
    math::compose_mod_reference(result, outer, inner, mixed_radix_divisor, mixed_radix_context);
    BOOST_CHECK(result == expected);

    math::compose_mod(result, outer, inner, mixed_radix_divisor, mixed_radix_context);
    BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_CASE(explicit_composition_blocks_preserve_the_outer_contract_and_match_horner) {
    using backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial = backend::polynomial_type;
    polynomial modulus(12), inner(11), outer(17);
    for (std::size_t i = 0; i < modulus.size(); ++i)
        modulus[i] = fq12_value(13 * i + 2);
    modulus.back() = fq12_value_type::one();
    for (std::size_t i = 0; i < inner.size(); ++i)
        inner[i] = fq12_value(17 * i + 3);
    for (std::size_t i = 0; i < outer.size(); ++i)
        outer[i] = fq12_value(19 * i + 5);

    for (const std::size_t prepared_bytes : {0U, 1024U * 1024U}) {
        for (const std::size_t power_limit : {1U, 7U, 20U}) {
            polynomial_arithmetic::polynomial_context_options options;
            options.modular_composition_cached_power_limit = power_limit;
            options.modular_composition_prepared_byte_limit = prepared_bytes;
            options.modular_composition_block_size = 9;    // Larger than ceil(sqrt(17)) = 5.
            polynomial_arithmetic::polynomial_context<backend> context(backend {}, options);
            const math::polynomial_divisor_context<backend> divisor(modulus, 10, context);
            const math::polynomial_composition_precomputation<backend> cache(inner, outer.size(), divisor, context);
            BOOST_CHECK_EQUAL(cache.block_size(), std::min<std::size_t>(9, power_limit));
            BOOST_CHECK_EQUAL(cache.maximum_outer_coefficient_count(), outer.size());
            if (power_limit > 1 && prepared_bytes != 0)
                BOOST_CHECK_GT(cache.prepared_coefficient_bytes(), 0);
            else
                BOOST_CHECK_EQUAL(cache.prepared_coefficient_bytes(), 0);
            polynomial expected, result = outer;
            math::compose_mod_reference(expected, outer, inner, divisor, context);
            math::compose_mod(result, result, cache, divisor, context);
            BOOST_CHECK(result == expected);
            math::compose_mod_reference(expected, expected, inner, divisor, context);
            math::compose_mod(result, result, cache, divisor, context);
            BOOST_CHECK(result == expected);

            const math::polynomial_composition_precomputation<backend> short_cache(inner, 3, divisor, context);
            BOOST_CHECK_EQUAL(short_cache.block_size(), std::min<std::size_t>(3, power_limit));
            BOOST_CHECK_THROW(math::compose_mod(result, outer, short_cache, divisor, context), std::invalid_argument);
        }
    }
    polynomial_arithmetic::polynomial_context<backend> automatic;
    const math::polynomial_divisor_context<backend> divisor(modulus, 10, automatic);
    const math::polynomial_composition_precomputation<backend> automatic_cache(inner, outer.size(), divisor, automatic);
    BOOST_CHECK_EQUAL(automatic_cache.block_size(), 5);
    polynomial_arithmetic::polynomial_context_options invalid_options;
    invalid_options.modular_composition_cached_power_limit = 0;
    invalid_options.modular_composition_block_size = 9;
    polynomial_arithmetic::polynomial_context<backend> invalid(backend {}, invalid_options);
    polynomial result = outer;
    BOOST_CHECK_THROW(math::compose_mod(result, outer, inner, divisor, invalid), std::invalid_argument);
    BOOST_CHECK(result == outer);
}

BOOST_AUTO_TEST_CASE(per_cache_block_limits_preserve_outer_size_and_other_caches) {
    using backend = grouped_backend_type;
    using polynomial = backend::polynomial_type;
    using reference_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    auto modulus = grouped_polynomial(13, 31);
    modulus[0] = modulus[1] = fq12_value_type::zero();    // A repeated zero root and a nonmonic leading coefficient.
    const auto inner = grouped_polynomial(12, 7), outer = grouped_polynomial(37, 43);
    polynomial_arithmetic::polynomial_context<reference_backend> reference;
    const math::polynomial_divisor_context<reference_backend> reference_divisor(modulus, 11, reference);
    polynomial expected;
    math::compose_mod_reference(expected, outer, inner, reference_divisor, reference);

    for (const std::size_t requested : {0U, 3U, 9U}) {
        for (const std::size_t power_limit : {2U, 64U}) {
            for (const std::size_t prepared_bytes : {0U, 4U * 1024U * 1024U}) {
                polynomial_arithmetic::polynomial_context_options options;
                options.modular_composition_block_size = requested;
                options.modular_composition_cached_power_limit = power_limit;
                options.modular_composition_prepared_byte_limit = prepared_bytes;
                options.modular_composition_giant_step_group_size = 8;
                options.modular_composition_giant_step_prepared_byte_limit = 4 * 1024 * 1024;
                polynomial_arithmetic::polynomial_context<backend> context {backend(32), options};
                const math::polynomial_divisor_context<backend> divisor(modulus, 11, context);
                const math::polynomial_composition_precomputation<backend> reusable(inner, outer.size(), divisor,
                                                                                    context);
                const auto reusable_width = reusable.block_size();
                for (const std::size_t limit : {1U, 4U, 20U}) {
                    BOOST_TEST_CONTEXT("requested=" << requested << ", power_limit=" << power_limit
                                                    << ", prepared_bytes=" << prepared_bytes << ", limit=" << limit) {
                        const math::polynomial_composition_precomputation<backend> temporary(
                            inner, outer.size(), divisor, context, prepared_bytes, limit);
                        BOOST_CHECK_EQUAL(temporary.block_size(), std::min(reusable_width, limit));
                        BOOST_CHECK_EQUAL(temporary.maximum_outer_coefficient_count(), outer.size());
                        BOOST_CHECK_LE(temporary.prepared_coefficient_bytes(), prepared_bytes);
                        polynomial result = outer;
                        math::compose_mod(result, result, temporary, divisor, context);
                        BOOST_CHECK(result == expected);
                        math::compose_mod(result, outer, reusable, divisor, context);
                        BOOST_CHECK(result == expected);
                        BOOST_CHECK_EQUAL(reusable.block_size(), reusable_width);
                        BOOST_CHECK_EQUAL(context.options().modular_composition_block_size, requested);
                        BOOST_CHECK_EQUAL(context.options().modular_composition_cached_power_limit, power_limit);
                    }
                }
                BOOST_CHECK_THROW(math::polynomial_composition_precomputation<backend>(inner, outer.size(), divisor,
                                                                                       context, prepared_bytes, 0),
                                  std::invalid_argument);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_cache_obeys_its_budget_and_pads_short_powers) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using precomputation_type = math::polynomial_composition_precomputation<backend_type>;
    using kernel = fields::detail::fp12_fast::fp12_fast<fq12_field_type::extension_policy::fast_params>;
    // k=5, inner=X: the nonconstant baby steps have lengths 2,3,4,5. Padding is included in the budget.
    const std::size_t required = 4 * 5 * sizeof(kernel::prepared_product_operand);
    const polynomial_type inner = {fq12_value_type::zero(), fq12_value_type::one()};
    polynomial_type divisor(8, fq12_value_type::zero());
    divisor[0] = fq12_value(3);
    divisor.back() = fq12_value(17);    // Nonmonic divisor.
    polynomial_type outer(21);
    for (std::size_t i = 0; i < outer.size(); ++i) {
        outer[i] = fq12_value(23 + 12 * i);
    }
    for (std::size_t limit :
         {std::size_t(0), required - 1, required, required + 1, std::numeric_limits<std::size_t>::max()}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = 5;
        options.modular_composition_prepared_byte_limit = limit;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> division(divisor, 6, context);
        const precomputation_type cache(inner, outer.size(), division, context);
        BOOST_CHECK_EQUAL(cache.block_size(), 5);
        BOOST_CHECK_EQUAL(cache.prepared_coefficient_bytes(), limit >= required ? required : 0);
        polynomial_type expected, actual;
        math::compose_mod_reference(expected, outer, inner, division, context);
        math::compose_mod(actual, outer, cache, division, context);
        BOOST_CHECK(actual == expected);

        // A temporary cache may use less prepared storage without changing this context or any existing cache.
        // Its explicit cap cannot override the context's own limit.
        for (std::size_t local_limit : {std::size_t(0), required - 1, required, required + 1}) {
            const precomputation_type temporary(inner, outer.size(), division, context, local_limit);
            BOOST_CHECK_EQUAL(temporary.prepared_coefficient_bytes(),
                              std::min(limit, local_limit) >= required ? required : 0);
            math::compose_mod(actual, outer, temporary, division, context);
            BOOST_CHECK(actual == expected);
            BOOST_CHECK_EQUAL(context.options().modular_composition_prepared_byte_limit, limit);
            BOOST_CHECK_EQUAL(cache.prepared_coefficient_bytes(), limit >= required ? required : 0);
        }

        // Check the shared single-block entry point too, including an existing accumulator and cancellation.
        for (std::size_t block = 0; block < 5; ++block) {
            const std::size_t first = block * 5;
            const std::size_t end = std::min(first + 5, outer.size());
            polynomial_type slice(outer.begin() + first, outer.begin() + end), evaluated;
            math::compose_mod_reference(evaluated, slice, inner, division, context);
            actual = {fq12_value(77), fq12_value(91)};
            math::addition(expected, actual, evaluated);
            math::detail::add_composition_block(actual, outer, cache, block, context);
            BOOST_CHECK(actual == expected);
            math::subtraction(actual, polynomial_type {}, evaluated);
            math::detail::add_composition_block(actual, outer, cache, block, context);
            BOOST_CHECK(actual == polynomial_type {fq12_value_type::zero()});
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_cache_owns_its_data_and_survives_copy_and_move) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using precomputation_type = math::polynomial_composition_precomputation<backend_type>;
    const polynomial_type inner = {fq12_value(1), fq12_value(13), fq12_value(25)};
    const polynomial_type divisor = {fq12_value(37), fq12_value(49), fq12_value(61), fq12_value(73), fq12_value(85)};
    std::optional<precomputation_type> saved;
    {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_prepared_byte_limit = 1024 * 1024;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> division(divisor, 3, context);
        polynomial_type temporary_inner = inner;
        saved.emplace(temporary_inner, 25, division, context);
        temporary_inner.assign(1, fq12_value(999));
    }
    BOOST_REQUIRE_GT(saved->prepared_coefficient_bytes(), 0);
    auto copy = *saved;
    const std::size_t bytes = copy.prepared_coefficient_bytes();
    saved.reset();
    auto moved = std::move(copy);
    BOOST_CHECK_EQUAL(moved.prepared_coefficient_bytes(), bytes);
    polynomial_arithmetic::polynomial_context<backend_type> context;
    math::polynomial_divisor_context<backend_type> division(divisor, 3, context);
    for (std::size_t count : {1, 4, 5, 6, 19, 20, 21, 25}) {
        polynomial_type outer(count);
        for (std::size_t i = 0; i < count; ++i) {
            outer[i] = i % 3 == 0 ? fq12_value_type::zero() : fq12_value(i * 12 + count);
        }
        outer.back() = fq12_value(123);
        polynomial_type expected, actual;
        math::compose_mod_reference(expected, outer, inner, division, context);
        math::compose_mod(actual, outer, moved, division, context);
        BOOST_CHECK(actual == expected);
        math::compose_mod(outer, outer, moved, division, context);
        BOOST_CHECK(outer == expected);
    }
    polynomial_type oversized(26, fq12_value(1)), unchanged {fq12_value(91)};
    BOOST_CHECK_THROW(math::compose_mod(unchanged, oversized, moved, division, context), std::invalid_argument);
    BOOST_CHECK(unchanged == polynomial_type {fq12_value(91)});
}

BOOST_AUTO_TEST_CASE(prepared_composition_dispatch_counts_nonzero_nonconstant_terms) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    constexpr std::size_t block_size = 17;
    polynomial_type divisor(10), inner(9);
    for (std::size_t i = 0; i < divisor.size(); ++i)
        divisor[i] = fq12_value(12 * i + 1);
    for (std::size_t i = 0; i < inner.size(); ++i)
        inner[i] = fq12_value(12 * i + 121);

    for (std::size_t lanes : {0, 4, 8}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = block_size;
        options.modular_composition_prepared_byte_limit = 4 * 1024 * 1024;
        options.modular_composition_simd_lanes = lanes;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> division(divisor, 8, context);
        const math::polynomial_composition_precomputation<backend_type> cache(inner, block_size * block_size, division,
                                                                              context);
        BOOST_REQUIRE_EQUAL(cache.block_size(), block_size);
        BOOST_REQUIRE(cache.prepared_coefficients().available());
#if defined(CRYPTO3_HAS_AVX512IFMA)
        BOOST_CHECK_EQUAL(cache.prepared_coefficients().simd_lanes(), lanes);
#else
        BOOST_CHECK_EQUAL(cache.prepared_coefficients().simd_lanes(), 0);
#endif
        for (std::size_t active_terms : {0, 1, 7, 8, 9, 16}) {
            polynomial_type outer(2 * block_size, fq12_value_type::zero());
            outer[0] = fq12_value(13);    // A nonzero constant must not move a seven-term sum above the cutoff.
            outer.back() = fq12_value(25);
            // Place active terms at the end of the block, including a lone high-exponent term.
            for (std::size_t j = 1; j <= active_terms; ++j)
                outer[block_size - j] = fq12_value(12 * j + 37);
            BOOST_CHECK_EQUAL(cache.prepared_coefficients().use_for_block(outer, 0), active_terms >= 8);
            BOOST_CHECK(!cache.prepared_coefficients().use_for_block(outer, 1));

            polynomial_type slice(outer.begin(), outer.begin() + block_size), evaluated;
            math::condense(slice);
            math::compose_mod_reference(evaluated, slice, inner, division, context);
            polynomial_type actual(11, fq12_value(61)), expected;
            math::addition(expected, actual, evaluated);
            math::detail::add_composition_block(actual, outer, cache, 0, context);
            BOOST_CHECK(actual == expected);
            math::subtraction(actual, polynomial_type {}, evaluated);
            math::detail::add_composition_block(actual, outer, cache, 0, context);
            BOOST_CHECK(actual == polynomial_type {fq12_value_type::zero()});
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_mixes_scalar_blocks_with_dense_batches_and_partial_tails) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    constexpr std::size_t block_size = 17;
    // Five adjacent dense blocks cross a four-block batch boundary. Subsequent blocks alternate between paths,
    // including an all-zero block, a constant block, and seven-term sums immediately below the cutoff.
    const std::array<std::size_t, 14> active_terms = {8, 16, 8, 9, 16, 7, 8, 0, 1, 8, 16, 7, 0, 8};
    polynomial_type divisor(10), inner(9);
    for (std::size_t i = 0; i < divisor.size(); ++i)
        divisor[i] = fq12_value(12 * i + 1);
    for (std::size_t i = 0; i < inner.size(); ++i)
        inner[i] = fq12_value(12 * i + 121);

    for (std::size_t lanes : {0, 4, 8}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = block_size;
        options.modular_composition_prepared_byte_limit = 4 * 1024 * 1024;
        options.modular_composition_simd_lanes = lanes;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> division(divisor, 8, context);
        const math::polynomial_composition_precomputation<backend_type> cache(inner, block_size * block_size, division,
                                                                              context);
        BOOST_REQUIRE_EQUAL(cache.block_size(), block_size);
        BOOST_REQUIRE(cache.prepared_coefficients().available());
        for (std::size_t tail_size : {1, 2, 8, 9, 17}) {
            polynomial_type outer(active_terms.size() * block_size + tail_size, fq12_value_type::zero());
            for (std::size_t block = 0; block < active_terms.size(); ++block) {
                const std::size_t first = block * block_size;
                if (block != 7)
                    outer[first] = fq12_value(12 * first + 217);
                for (std::size_t j = 1; j <= active_terms[block]; ++j)
                    outer[first + block_size - j] = fq12_value(12 * (first + j) + 229);
                BOOST_CHECK_EQUAL(cache.prepared_coefficients().use_for_block(outer, block), active_terms[block] >= 8);
            }
            for (std::size_t i = active_terms.size() * block_size; i < outer.size(); ++i)
                outer[i] = fq12_value(12 * i + 241);
            BOOST_CHECK_EQUAL(cache.prepared_coefficients().use_for_block(outer, active_terms.size()), tail_size >= 9);
            polynomial_type expected, actual;
            math::compose_mod_reference(expected, outer, inner, division, context);
            math::compose_mod(actual, outer, cache, division, context);
            BOOST_CHECK(actual == expected);
            math::compose_mod(outer, outer, cache, division, context);
            BOOST_CHECK(outer == expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_tiles_match_horner_at_block_and_coefficient_boundaries) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    // Dense full blocks contain eight nonconstant terms, so they still exercise the prepared kernels.
    constexpr std::size_t block_size = 9;
    for (std::size_t degree : {1, 3, 4, 5, 9}) {
        polynomial_type divisor(degree + 1), inner(degree);
        for (std::size_t i = 0; i <= degree; ++i) {
            divisor[i] = fq12_value(i * 12 + 1);
            if (i < degree) {
                inner[i] = fq12_value(i * 12 + 73);
            }
        }
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = block_size;
        options.modular_composition_prepared_byte_limit = 1024 * 1024;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> division(divisor, std::max(std::size_t(1), degree - 1), context);
        const math::polynomial_composition_precomputation<backend_type> cache(inner, 82, division, context);
        BOOST_REQUIRE_EQUAL(cache.block_size(), block_size);
        BOOST_REQUIRE_GT(cache.prepared_coefficient_bytes(), 0);
        for (std::size_t count : {1, 8, 9, 10, 35, 36, 37, 44, 45, 46, 71, 72, 73, 81, 82}) {
            for (bool sparse : {false, true}) {
                polynomial_type outer(count, fq12_value_type::zero());
                for (std::size_t i = 0; i < count; ++i) {
                    if (!sparse || i % block_size == 0 || i + 1 == count) {
                        outer[i] = fq12_value(1 + 12 * i);
                    }
                }
                polynomial_type expected, actual;
                math::compose_mod_reference(expected, outer, inner, division, context);
                math::compose_mod(actual, outer, cache, division, context);
                BOOST_CHECK(actual == expected);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_supports_123_terms_and_the_small_prime_backend) {
    using backend_type = polynomial_arithmetic::bn254_small_prime_fft_backend;
    using polynomial_type = typename backend_type::polynomial_type;
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_cached_power_limit = 123;
    options.modular_composition_prepared_byte_limit = 4 * 1024 * 1024;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type(32), options);
    const polynomial_type divisor = {fq12_value(1),  fq12_value(13), fq12_value(25),
                                     fq12_value(37), fq12_value(49), fq12_value(61)};
    const polynomial_type inner = {fq12_value(73), fq12_value(85), fq12_value(97), fq12_value(109), fq12_value(121)};
    math::polynomial_divisor_context<backend_type> division(divisor, 4, context);
    const math::polynomial_composition_precomputation<backend_type> cache(inner, 122 * 122 + 1, division, context);
    BOOST_REQUIRE_EQUAL(cache.block_size(), 123);
    BOOST_REQUIRE_GT(cache.prepared_coefficient_bytes(), 0);
    using reference_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    polynomial_arithmetic::polynomial_context<reference_backend> reference_context;
    math::polynomial_divisor_context<reference_backend> reference_division(divisor, 4, reference_context);
    for (std::size_t count : {122, 123, 124, 491, 492, 493, 615}) {
        polynomial_type outer(count);
        for (std::size_t i = 0; i < count; ++i) {
            outer[i] = fq12_value(i * 12 + 1);
        }
        polynomial_type expected, actual;
        math::compose_mod_reference(expected, outer, inner, reference_division, reference_context);
        math::compose_mod(actual, outer, cache, division, context);
        BOOST_CHECK(actual == expected);
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_fallbacks_handle_constant_moduli_and_unsupported_fields) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_prepared_byte_limit = 1024 * 1024;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
    const polynomial_type inner = {fq12_value(1), fq12_value(13)};
    const polynomial_type outer = {fq12_value(25), fq12_value(37), fq12_value(49)};
    math::polynomial_divisor_context<backend_type> constant_division({fq12_value(7)}, 1, context);
    const math::polynomial_composition_precomputation<backend_type> constant_cache(inner, 3, constant_division,
                                                                                   context);
    BOOST_CHECK_EQUAL(constant_cache.prepared_coefficient_bytes(), 0);
    polynomial_type result = outer;
    math::compose_mod(result, result, constant_cache, constant_division, context);
    BOOST_CHECK(result == polynomial_type {fq12_value_type::zero()});

    options.modular_composition_cached_power_limit = 1;
    polynomial_arithmetic::polynomial_context<backend_type> scalar_context(backend_type {}, options);
    math::polynomial_divisor_context<backend_type> division({fq12_value(1), fq12_value(13), fq12_value(25)}, 1,
                                                            scalar_context);
    const math::polynomial_composition_precomputation<backend_type> scalar_cache(inner, 3, division, scalar_context);
    BOOST_CHECK_EQUAL(scalar_cache.prepared_coefficient_bytes(), 0);
    polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, division, scalar_context);
    math::compose_mod(result, outer, scalar_cache, division, scalar_context);
    BOOST_CHECK(result == expected);

    using base_backend = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    options.modular_composition_cached_power_limit = 3;
    polynomial_arithmetic::polynomial_context<base_backend> base_context(base_backend {}, options);
    math::polynomial_divisor_context<base_backend> base_division({fq_value_type(1), fq_value_type(1)}, 1, base_context);
    const math::polynomial_composition_precomputation<base_backend> base_cache({fq_value_type(2)}, 3, base_division,
                                                                               base_context);
    BOOST_CHECK_EQUAL(base_cache.prepared_coefficient_bytes(), 0);
    base_backend::polynomial_type base_result;
    math::compose_mod(base_result, {fq_value_type(1), fq_value_type(2), fq_value_type(3)}, base_cache, base_division,
                      base_context);
    BOOST_CHECK(base_result == base_backend::polynomial_type {fq_value_type(17)});
}

BOOST_AUTO_TEST_CASE(prepared_composition_handles_degenerate_inputs_and_one_off_aliases) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_prepared_byte_limit = 1024 * 1024;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
    const polynomial_type divisor = {fq12_value(1), fq12_value(13), fq12_value(25), fq12_value(37), fq12_value(49)};
    math::polynomial_divisor_context<backend_type> division(divisor, 3, context);
    for (const polynomial_type inner : {polynomial_type {fq12_value_type::zero()}, polynomial_type {fq12_value(61)},
                                        polynomial_type {fq12_value(73), fq12_value(85), fq12_value(97)}}) {
        const math::polynomial_composition_precomputation<backend_type> cache(inner, 25, division, context);
        BOOST_REQUIRE_GT(cache.prepared_coefficient_bytes(), 0);
        for (const polynomial_type outer :
             {polynomial_type {fq12_value_type::zero()}, polynomial_type {fq12_value(109)},
              polynomial_type {-inner[0], fq12_value_type::one()},
              polynomial_type {fq12_value(121), fq12_value_type::zero(), fq12_value(133)}}) {
            polynomial_type expected, result;
            math::compose_mod_reference(expected, outer, inner, division, context);
            math::compose_mod(result, outer, cache, division, context);
            BOOST_CHECK(result == expected);
            auto outer_alias = outer;
            math::compose_mod(outer_alias, outer_alias, inner, division, context);
            BOOST_CHECK(outer_alias == expected);
            auto inner_alias = inner;
            math::compose_mod(inner_alias, outer, inner_alias, division, context);
            BOOST_CHECK(inner_alias == expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_simd_panels_match_scalar_at_lane_and_term_boundaries) {
    using polynomial_type = math::polynomial<fq12_value_type>;
    using cache_type = math::detail::prepared_composition_coefficients<fq12_value_type>;
    constexpr std::size_t budget = 16 * 1024 * 1024;
    for (std::size_t width : {3, 4, 5, 7, 8, 9, 15, 17}) {
        std::vector<polynomial_type> powers(18);
        powers[0] = {fq12_value_type::one()};
        for (std::size_t j = 1; j < powers.size(); ++j) {
            powers[j].resize(j % 3 == 0 ? 1 : width);
            for (std::size_t c = 0; c < powers[j].size(); ++c)
                powers[j][c] = fq12_value(12 * (j * width + c) + 1);
        }
        cache_type scalar(powers, budget, 0);
        for (std::size_t lanes : {4, 8}) {
            cache_type candidate(powers, budget, lanes);
            BOOST_REQUIRE(candidate.available());
#if defined(CRYPTO3_HAS_AVX512IFMA)
            BOOST_CHECK_EQUAL(candidate.simd_lanes(), width >= 4 ? lanes : 0);
            if (width >= 4) {
                const std::size_t expected =
                    ((width + lanes - 1) / lanes) * lanes * 23 * 17 * 5 * sizeof(std::uint64_t);
                BOOST_CHECK_EQUAL(candidate.bytes(), expected);
            }
#else
            BOOST_CHECK_EQUAL(candidate.simd_lanes(), 0);
#endif
            BOOST_CHECK_LE(candidate.bytes(), budget);
            const auto bytes = candidate.bytes();
            cache_type exact(powers, bytes, lanes), small(powers, bytes - 1, lanes);
            BOOST_CHECK(exact.available());
            BOOST_CHECK_LE(small.bytes(), bytes - 1);
            auto copied = candidate;
            auto moved = std::move(copied);
            typename cache_type::workspace scratch, reference_scratch;
            std::array<polynomial_type, 4> actual, expected;
            for (bool sparse : {false, true}) {
                polynomial_type outer(5 * powers.size() - 7);
                for (std::size_t i = 0; i < outer.size(); ++i) {
                    outer[i] = sparse && i % 7 ? fq12_value_type::zero() : fq12_value(12 * i + 1);
                }
                // Include a wholly constant block amid active blocks, then a short final block.
                for (std::size_t j = 1; j < powers.size(); ++j)
                    outer[powers.size() + j] = fq12_value_type::zero();
                for (std::size_t first : {0, 1, 4}) {
                    const std::size_t count = std::min(std::size_t(4), std::size_t(5) - first);
                    moved.evaluate_blocks(actual, outer, first, count, scratch);
                    scalar.evaluate_blocks(expected, outer, first, count, reference_scratch);
                    for (std::size_t i = 0; i < count; ++i)
                        BOOST_CHECK(actual[i] == expected[i]);
                }
            }
        }
        BOOST_CHECK_THROW(cache_type(powers, budget, 3), std::invalid_argument);
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_panels_support_independent_ranges_and_workspace_reuse) {
    using polynomial_type = grouped_polynomial_type;
    using cache_type = math::detail::prepared_composition_coefficients<fq12_value_type>;
    constexpr std::size_t block_size = 10, blocks = 3;
    for (std::size_t width : {3, 9, 17}) {
        std::vector<polynomial_type> powers(block_size);
        powers[0] = {fq12_value_type::one()};
        for (std::size_t exponent = 1; exponent < powers.size(); ++exponent) {
            powers[exponent] = grouped_polynomial(exponent % 3 == 0 ? 1 : width, 13 * exponent);
        }
        for (std::size_t lanes : {0, 4, 8}) {
            const cache_type cache(powers, 4 * 1024 * 1024, lanes);
            BOOST_REQUIRE(cache.available());
            typename cache_type::workspace scratch;
            std::array<polynomial_type, cache_type::block_tile_size> output;
            output.back() = {fq12_value(77)};
            for (unsigned pattern = 0; pattern < 3; ++pattern) {
                polynomial_type outer(block_size * blocks - 2, fq12_value_type::zero());
                for (std::size_t i = 0; i < outer.size(); ++i) {
                    if (pattern == 0 || (pattern == 1 && i % 3 == 0) || i % block_size == 0) {
                        outer[i] = fq12_value(17 + 12 * i);
                    }
                }
                const auto evaluation = cache.prepare_blocks(output, outer, 0, blocks, scratch);
                const auto prepared_output = output;
                BOOST_CHECK_THROW(evaluation.evaluate_panels(1, 0), std::out_of_range);
                BOOST_CHECK_THROW(evaluation.evaluate_panels(0, evaluation.panel_count() + 1), std::out_of_range);
                BOOST_CHECK(output == prepared_output);
                evaluation.evaluate_panels(evaluation.panel_count(), evaluation.panel_count());
                for (std::size_t panel = evaluation.panel_count(); panel-- > 0;) {
                    evaluation.evaluate_panels(panel, panel + 1);
                }
                for (std::size_t block = 0; block < blocks; ++block) {
                    BOOST_CHECK_EQUAL(output[block].size(), prepared_output[block].size());
                }
                evaluation.finish();
                for (std::size_t block = 0; block < blocks; ++block) {
                    polynomial_type expected(width, fq12_value_type::zero());
                    for (std::size_t exponent = 0;
                         exponent < block_size && block * block_size + exponent < outer.size();
                         ++exponent) {
                        for (std::size_t coefficient = 0; coefficient < powers[exponent].size(); ++coefficient) {
                            expected[coefficient] +=
                                powers[exponent][coefficient] * outer[block * block_size + exponent];
                        }
                    }
                    math::condense(expected);
                    BOOST_CHECK(output[block] == expected);
                }
                BOOST_CHECK(output.back() == polynomial_type {fq12_value(77)});
                if (pattern == 2) {
                    BOOST_CHECK_EQUAL(evaluation.panel_count(), 0);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(application_composition_evaluation_partitions_panels_across_threads) {
    using backend_type = partitioned_composition_backend;
    using polynomial_type = backend_type::polynomial_type;
    using cache_type = math::detail::prepared_composition_coefficients<fq12_value_type>;
    constexpr std::size_t block_size = 10, blocks = 3;
    std::vector<polynomial_type> powers(block_size);
    powers[0] = {fq12_value_type::one()};
    for (std::size_t exponent = 1; exponent < powers.size(); ++exponent) {
        // Nine SIMD8 panels give every worker a nonempty range, including a partial final panel.
        powers[exponent] = grouped_polynomial(exponent % 3 == 0 ? 64 : 65, 13 * exponent);
    }
    for (std::size_t lanes : {0, 4, 8}) {
        const cache_type cache(powers, 4 * 1024 * 1024, lanes);
        BOOST_REQUIRE(cache.available());
        for (std::size_t threads : {1, 2, 4, 8}) {
            composition_execution_counts counts;
            polynomial_arithmetic::polynomial_context<backend_type> context(backend_type(64, counts, threads));
            typename cache_type::workspace scratch, reference_scratch;
            std::array<polynomial_type, cache_type::block_tile_size> output, expected;
            for (unsigned pattern = 0; pattern < 3; ++pattern) {
                polynomial_type outer(block_size * blocks - 2, fq12_value_type::zero());
                for (std::size_t i = 0; i < outer.size(); ++i) {
                    if (pattern == 0 || (pattern == 1 && i % 3 == 0) || i % block_size == 0) {
                        outer[i] = fq12_value(17 + 12 * i);
                    }
                }
                cache.evaluate_blocks(expected, outer, 0, blocks, reference_scratch);
                const auto evaluation = cache.prepare_blocks(output, outer, 0, blocks, scratch);
                context.evaluate_composition_coefficients(evaluation);
                evaluation.finish();
                BOOST_CHECK(output == expected);
            }
            BOOST_CHECK_EQUAL(counts.batches, 3);
        }
    }
}

BOOST_AUTO_TEST_CASE(application_composition_evaluation_preserves_grouping_aliases_and_failure_reuse) {
    using backend_type = partitioned_composition_backend;
    using polynomial_type = backend_type::polynomial_type;
    using reference_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    constexpr std::size_t degree = 19;
    const auto modulus = grouped_polynomial(degree + 1, 41);
    const auto inner = grouped_polynomial(degree, 73);
    const auto outer = grouped_polynomial(9 * 9 + 3, 101);
    polynomial_arithmetic::polynomial_context<reference_backend> reference;
    math::polynomial_divisor_context<reference_backend> reference_divisor(modulus, degree - 1, reference);
    polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, reference_divisor, reference);
    for (std::size_t group : {0, 8}) {
        for (std::size_t bytes : {0, 4 * 1024 * 1024}) {
            for (std::size_t threads : {1, 2, 4, 8}) {
                composition_execution_counts counts;
                polynomial_arithmetic::polynomial_context_options options;
                options.modular_composition_block_size = 9;
                options.modular_composition_prepared_byte_limit = bytes;
                options.modular_composition_giant_step_group_size = group;
                options.modular_composition_giant_step_prepared_byte_limit = 16 * 1024 * 1024;
                polynomial_arithmetic::polynomial_context<backend_type> context(backend_type(64, counts, threads),
                                                                                options);
                math::polynomial_divisor_context<backend_type> divisor(modulus, degree - 1, context);
                const math::polynomial_composition_precomputation<backend_type> cache(inner, outer.size(), divisor,
                                                                                      context);
                polynomial_type output = outer;
                math::compose_mod(output, output, cache, divisor, context);
                BOOST_CHECK(output == expected);
                BOOST_CHECK_EQUAL(counts.batches != 0, bytes != 0);
                if (bytes != 0) {
                    counts.fail = true;
                    output = outer;
                    BOOST_CHECK_THROW(math::compose_mod(output, output, cache, divisor, context), std::runtime_error);
                    BOOST_CHECK(output == outer);
                    counts.fail = false;
                    math::compose_mod(output, output, cache, divisor, context);
                    BOOST_CHECK(output == expected);

                    // Exercise the single-block helper independently of the batched composition paths.
                    const polynomial_type block(outer.begin(), outer.begin() + 9);
                    polynomial_type expected_block;
                    math::compose_mod_reference(expected_block, block, inner, reference_divisor, reference);
                    const auto before = counts.batches;
                    output = {fq12_value_type::zero()};
                    math::detail::add_composition_block(output, outer, cache, 0, context);
                    BOOST_CHECK_EQUAL(counts.batches, before + 1);
                    BOOST_CHECK(output == expected_block);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_simd_complete_composition_matches_horner_and_aliases) {
    using backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = backend::polynomial_type;
    polynomial_type divisor(10), inner(8), outer(83);
    for (std::size_t i = 0; i < divisor.size(); ++i)
        divisor[i] = fq12_value(12 * i + 1);
    for (std::size_t i = 0; i < inner.size(); ++i)
        inner[i] = fq12_value(12 * i + 121);
    for (std::size_t i = 0; i < outer.size(); ++i)
        outer[i] = i % 3 ? fq12_value(12 * i + 217) : fq12_value_type::zero();
    for (std::size_t lanes : {0, 4, 8}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = 17;
        options.modular_composition_prepared_byte_limit = 4 * 1024 * 1024;
        options.modular_composition_simd_lanes = lanes;
        polynomial_arithmetic::polynomial_context<backend> context(backend {}, options);
        math::polynomial_divisor_context<backend> division(divisor, 8, context);
        const math::polynomial_composition_precomputation<backend> cache(inner, 257, division, context);
        polynomial_type expected, actual = outer;
        math::compose_mod_reference(expected, outer, inner, division, context);
        math::compose_mod(actual, actual, cache, division, context);
        BOOST_CHECK(actual == expected);
    }
}

BOOST_AUTO_TEST_CASE(prepared_simd_cache_handles_maximal_residues_in_246_coefficient_blocks) {
    using polynomial_type = math::polynomial<fq12_value_type>;
    using cache_type = math::detail::prepared_composition_coefficients<fq12_value_type>;
    constexpr std::size_t block_size = 246, budget = 16 * 1024 * 1024;
    auto maximum = fq_value_type::modulus_params.get_mod_obj().get_mod();
    eval_subtract(maximum, boost::multiprecision::limb_type(1));
    auto maximal_flat_real = maximum;
    eval_subtract(maximal_flat_real, boost::multiprecision::limb_type(9));
    fq12_value_type maximal_tower, maximal_flat;
    for (std::size_t w = 0; w < 2; ++w) {
        for (std::size_t v = 0; v < 3; ++v) {
            maximal_flat.data[w].data[v].data[0].data.backend().base_data() = maximal_flat_real;
            maximal_flat.data[w].data[v].data[1].data.backend().base_data() = maximum;
            for (std::size_t u = 0; u < 2; ++u) {
                maximal_tower.data[w].data[v].data[u].data.backend().base_data() = maximum;
            }
        }
    }
    for (const std::size_t width : {5, 9}) {
        std::vector<polynomial_type> powers(block_size);
        powers[0] = {fq12_value_type::one()};
        for (std::size_t exponent = 1; exponent < block_size; ++exponent) {
            // Short powers and explicit zero coefficients both require zero padding in the cache.
            powers[exponent].resize(exponent % 11 == 0 ? 1 : width);
            for (std::size_t coefficient = 0; coefficient < powers[exponent].size(); ++coefficient) {
                const std::array<fq12_value_type, 4> values = {maximal_tower, maximal_flat, fq12_value_type::zero(),
                                                               fq12_value(12 * exponent + coefficient)};
                powers[exponent][coefficient] = values[(exponent + coefficient) % values.size()];
            }
        }
        polynomial_type outer(3 * block_size + 17);
        for (std::size_t i = 0; i < outer.size(); ++i) {
            outer[i] = i % 13 == 0 ? fq12_value_type::zero() : (i % 2 == 0 ? maximal_tower : maximal_flat);
        }
        // Reuse scratch across dense blocks, a constant block, and a partial final block.
        for (std::size_t i = block_size + 1; i < 2 * block_size; ++i) {
            outer[i] = fq12_value_type::zero();
        }
        std::array<polynomial_type, 4> expected;
        for (std::size_t block = 0; block < expected.size(); ++block) {
            const std::size_t first = block * block_size, count = std::min(block_size, outer.size() - first);
            expected[block].assign(width, fq12_value_type::zero());
            expected[block][0] = outer[first];
            // Direct coefficient sums avoid both prepared product decompositions.
            for (std::size_t exponent = 1; exponent < count; ++exponent) {
                for (std::size_t coefficient = 0; coefficient < powers[exponent].size(); ++coefficient) {
                    expected[block][coefficient] += outer[first + exponent] * powers[exponent][coefficient];
                }
            }
            math::condense(expected[block]);
        }
        for (const std::size_t lanes : {0, 4, 8}) {
            const cache_type cache(powers, budget, lanes);
            BOOST_REQUIRE(cache.available());
#if defined(CRYPTO3_HAS_AVX512IFMA)
            BOOST_CHECK_EQUAL(cache.simd_lanes(), lanes);
#else
            BOOST_CHECK_EQUAL(cache.simd_lanes(), 0);
#endif
            typename cache_type::workspace scratch;
            std::array<polynomial_type, 4> actual;
            cache.evaluate_blocks(actual, outer, 0, 4, scratch);
            BOOST_CHECK(actual == expected);
            cache.evaluate_blocks(actual, outer, 1, 3, scratch);
            for (std::size_t i = 0; i < 3; ++i) {
                BOOST_CHECK(actual[i] == expected[i + 1]);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_simd_configuration_keeps_six_limb_fields_on_scalar_kernel) {
    using value = fields::fp12_2over3over2<fields::bls12<381>>::value_type;
    using polynomial = math::polynomial<value>;
    using cache_type = math::detail::prepared_composition_coefficients<value>;
    std::vector<polynomial> powers(9, polynomial(4, value::one()));
    powers[0] = {value::one()};
    const cache_type cache(powers, 1 << 20, 8);
    BOOST_REQUIRE(cache.available());
    BOOST_CHECK_EQUAL(cache.simd_lanes(), 0);
    typename cache_type::workspace scratch;
    std::array<polynomial, 4> output;
    cache.evaluate_blocks(output, polynomial(9, value::one()), 0, 1, scratch);
    BOOST_CHECK(output[0] == polynomial({value(9), value(8), value(8), value(8)}));
}

BOOST_AUTO_TEST_CASE(grouped_composition_matches_reference_for_partial_groups_and_reducible_nonmonic_moduli) {
    using backend_type = grouped_backend_type;
    using polynomial_type = grouped_polynomial_type;
    using reference_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    polynomial_arithmetic::polynomial_context<reference_backend> reference;
    // X^2*D has a repeated zero root; D's extension-field leading coefficient is nonmonic.
    auto divisor = grouped_polynomial(18, 31);
    divisor.get_storage().insert(divisor.begin(), 2, fq12_value_type::zero());
    const auto inner = grouped_polynomial(25, 7);
    math::polynomial_divisor_context<reference_backend> reference_divisor(divisor, 18, reference);
    std::vector<polynomial_type> outers, expected;
    for (const std::size_t size : {1, 2, 9, 10, 35, 36, 37, 71, 72, 73, 89}) {
        outers.push_back(grouped_polynomial(size, 43));
        expected.emplace_back();
        math::compose_mod_reference(expected.back(), outers.back(), inner, reference_divisor, reference);
    }
    for (const std::size_t group : {1, 4, 8}) {
        for (const std::size_t coefficient_bytes : {0U, 4U * 1024U * 1024U}) {
            BOOST_TEST_CONTEXT("group=" << group << ", coefficient_bytes=" << coefficient_bytes) {
                polynomial_arithmetic::polynomial_context_options options;
                options.basecase_divisor_coefficient_cutoff = 0;
                options.basecase_quotient_coefficient_cutoff = 0;
                options.modular_composition_block_size = 9;
                options.modular_composition_prepared_byte_limit = coefficient_bytes;
                options.modular_composition_giant_step_group_size = group;
                options.modular_composition_giant_step_prepared_byte_limit = 4 * 1024 * 1024;
                polynomial_arithmetic::polynomial_context<backend_type> context {backend_type(64), options};
                math::polynomial_divisor_context<backend_type> fixed(divisor, 18, context);
                math::polynomial_composition_precomputation<backend_type> cache(inner, 89, fixed, context);
                BOOST_REQUIRE_EQUAL(cache.giant_step_group_size(), group);
                BOOST_CHECK_EQUAL(cache.prepared_coefficients().available(), coefficient_bytes != 0);
                for (std::size_t i = 0; i < outers.size(); ++i) {
                    polynomial_type actual;
                    math::compose_mod(actual, outers[i], cache, fixed, context);
                    BOOST_CHECK(actual == expected[i]);
                    actual = outers[i];
                    math::compose_mod(actual, actual, cache, fixed, context);
                    BOOST_CHECK(actual == expected[i]);
                    BOOST_CHECK_EQUAL(cache.prepared_giant_steps().accumulator->product_count(), 0);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(grouped_composition_reduces_once_per_group_without_adding_products) {
    using backend_type = grouped_counting_backend;
    using polynomial_type = grouped_polynomial_type;
    const auto outer = grouped_polynomial(540, 17);    // 60 blocks, hence 59 giant products.
    const auto inner = grouped_polynomial(19, 3);
    auto divisor = grouped_polynomial(20, 71);
    divisor.back() = fq12_value_type::one();
    polynomial_type expected;
    for (const std::size_t group : {0, 1, 4, 8}) {
        grouped_operation_counts counts;
        polynomial_arithmetic::polynomial_context_options options;
        options.basecase_divisor_coefficient_cutoff = 0;
        options.basecase_quotient_coefficient_cutoff = 0;
        options.modular_composition_block_size = 9;
        options.modular_composition_prepared_byte_limit = 4 * 1024 * 1024;
        options.modular_composition_giant_step_group_size = group;
        options.modular_composition_giant_step_prepared_byte_limit = 4 * 1024 * 1024;
        polynomial_arithmetic::polynomial_context<backend_type> context {backend_type(64, counts), options};
        math::polynomial_divisor_context<backend_type> fixed(divisor, 18, context);
        math::polynomial_composition_precomputation<backend_type> cache(inner, outer.size(), fixed, context);
        counts = {};
        polynomial_type actual;
        math::compose_mod(actual, outer, cache, fixed, context);
        if (group == 0) {
            expected = actual;
            BOOST_CHECK_EQUAL(counts.products, 0);
            BOOST_CHECK_EQUAL(counts.finishes, 0);
            BOOST_CHECK_EQUAL(counts.remainders, 59);
        } else {
            BOOST_CHECK(actual == expected);
            BOOST_CHECK_EQUAL(counts.products, 59);
            BOOST_CHECK_EQUAL(counts.finishes, (59 + group - 1) / group);
            BOOST_CHECK_EQUAL(counts.remainders, counts.finishes);
        }
    }
}

BOOST_AUTO_TEST_CASE(grouped_composition_respects_cache_budgets_and_backend_capabilities) {
    using backend_type = grouped_backend_type;
    const auto outer = grouped_polynomial(89, 11);
    const auto inner = grouped_polynomial(19, 3);
    const auto divisor = grouped_polynomial(20, 7);
    polynomial_arithmetic::polynomial_context<backend_type> ordinary {backend_type(64)};
    math::polynomial_divisor_context<backend_type> ordinary_fixed(divisor, 18, ordinary);
    grouped_polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, ordinary_fixed, ordinary);
    const auto shape = ordinary.make_product_accumulator(19, 4);
    BOOST_REQUIRE(shape.has_value());
    const std::size_t required = 4 * shape->required_spectrum_bytes();
    for (const std::size_t limit : {std::size_t(0), required - 1, required}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_block_size = 9;
        options.modular_composition_giant_step_group_size = 4;
        options.modular_composition_giant_step_prepared_byte_limit = limit;
        polynomial_arithmetic::polynomial_context<backend_type> context {backend_type(64), options};
        math::polynomial_divisor_context<backend_type> fixed(divisor, 18, context);
        math::polynomial_composition_precomputation<backend_type> cache(inner, outer.size(), fixed, context);
        BOOST_CHECK_EQUAL(cache.giant_step_group_size(), limit == required ? 4 : 0);
        BOOST_CHECK_EQUAL(cache.prepared_giant_step_bytes(), limit == required ? required : 0);
        grouped_polynomial_type actual;
        math::compose_mod(actual, outer, cache, fixed, context);
        BOOST_CHECK(actual == expected);
    }
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_block_size = 9;
    options.modular_composition_giant_step_group_size = 8;
    options.modular_composition_giant_step_prepared_byte_limit = 4 * 1024 * 1024;
    using reference_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    polynomial_arithmetic::polynomial_context<reference_backend> unsupported {reference_backend {}, options};
    math::polynomial_divisor_context<reference_backend> unsupported_fixed(divisor, 18, unsupported);
    math::polynomial_composition_precomputation<reference_backend> unsupported_cache(inner, outer.size(),
                                                                                     unsupported_fixed, unsupported);
    BOOST_CHECK_EQUAL(unsupported_cache.giant_step_group_size(), 0);
    BOOST_CHECK_EQUAL(unsupported_cache.prepared_giant_step_bytes(), 0);
    grouped_polynomial_type actual;
    math::compose_mod(actual, outer, unsupported_cache, unsupported_fixed, unsupported);
    BOOST_CHECK(actual == expected);

    polynomial_arithmetic::polynomial_context<backend_type> context {backend_type(64), options};
    math::polynomial_divisor_context<backend_type> fixed(divisor, 18, context);
    math::polynomial_composition_precomputation<backend_type> one_block(inner, 9, fixed, context);
    BOOST_CHECK_EQUAL(one_block.giant_step_group_size(), 0);
    math::polynomial_composition_precomputation<backend_type> two_blocks(inner, 18, fixed, context);
    BOOST_CHECK_EQUAL(two_blocks.giant_step_group_size(), 1);
}

BOOST_AUTO_TEST_CASE(grouped_composition_falls_back_after_partial_work_without_changing_aliased_input) {
    using backend_type = grouped_counting_backend;
    const auto outer = grouped_polynomial(89, 11);
    const auto inner = grouped_polynomial(19, 3);
    const auto divisor = grouped_polynomial(20, 7);
    grouped_operation_counts counts;
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_block_size = 9;
    options.modular_composition_prepared_byte_limit = 4 * 1024 * 1024;
    options.modular_composition_giant_step_group_size = 4;
    options.modular_composition_giant_step_prepared_byte_limit = 4 * 1024 * 1024;
    polynomial_arithmetic::polynomial_context<backend_type> context {backend_type(64, counts), options};
    math::polynomial_divisor_context<backend_type> fixed(divisor, 18, context);
    math::polynomial_composition_precomputation<backend_type> cache(inner, outer.size(), fixed, context);
    grouped_polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, fixed, context);
    for (const std::size_t failure : {1, 3, 5}) {
        counts = {};
        counts.decline_on_append = failure;
        auto actual = outer;
        math::compose_mod(actual, actual, cache, fixed, context);
        BOOST_CHECK(actual == expected);
        BOOST_CHECK_EQUAL(counts.append_attempts, failure);
    }
    counts = {};
    counts.decline_finish = true;
    auto actual = outer;
    math::compose_mod(actual, actual, cache, fixed, context);
    BOOST_CHECK(actual == expected);
    BOOST_CHECK_GT(counts.products, 0);
    BOOST_CHECK_EQUAL(counts.finishes, 0);
    counts = {};
    counts.decline_preparation = true;
    math::polynomial_composition_precomputation<backend_type> declined(inner, outer.size(), fixed, context);
    BOOST_CHECK_EQUAL(declined.giant_step_group_size(), 0);
    math::compose_mod(actual, outer, declined, fixed, context);
    BOOST_CHECK(actual == expected);
}

BOOST_AUTO_TEST_CASE(grouped_composition_mixes_sparse_blocks_and_owns_prepared_state) {
    using backend_type = grouped_backend_type;
    using cache_type = math::polynomial_composition_precomputation<backend_type>;
    const auto inner = grouped_polynomial(19, 3);
    const auto divisor = grouped_polynomial(20, 7);
    auto outer = grouped_polynomial(89, 11);
    // Mix dense tiles, scalar blocks and a partial tail across both grouping and tile boundaries.
    for (std::size_t block : {1, 4, 6}) {
        for (std::size_t i = 1; i < 9; ++i) {
            outer[block * 9 + i] = fq12_value_type::zero();
        }
    }
    auto cache = [&] {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_block_size = 9;
        options.modular_composition_prepared_byte_limit = 4 * 1024 * 1024;
        options.modular_composition_giant_step_group_size = 8;
        options.modular_composition_giant_step_prepared_byte_limit = 4 * 1024 * 1024;
        polynomial_arithmetic::polynomial_context<backend_type> context {backend_type(128), options};
        auto source_inner = inner;
        math::polynomial_divisor_context<backend_type> fixed(divisor, 18, context);
        cache_type snapshot(source_inner, outer.size(), fixed, context);
        source_inner[0] += fq12_value_type::one();
        return snapshot;
    }();
    auto copied = cache;
    auto moved = std::move(copied);
    polynomial_arithmetic::polynomial_context<backend_type> context {backend_type(64)};
    math::polynomial_divisor_context<backend_type> fixed(divisor, 18, context);
    grouped_polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, fixed, context);
    for (const auto *snapshot : {&cache, &moved}) {
        BOOST_REQUIRE_EQUAL(snapshot->giant_step_group_size(), 8);
        for (unsigned repeat = 0; repeat < 2; ++repeat) {
            auto actual = outer;
            math::compose_mod(actual, actual, *snapshot, fixed, context);
            BOOST_CHECK(actual == expected);
            BOOST_CHECK_EQUAL(snapshot->prepared_giant_steps().accumulator->product_count(), 0);
            BOOST_CHECK_EQUAL(snapshot->prepared_giant_steps().accumulator->storage_bytes(), 0);
        }
    }
}

BOOST_AUTO_TEST_CASE(grouped_composition_handles_zero_constant_and_nilpotent_giants) {
    using backend_type = grouped_backend_type;
    using polynomial_type = grouped_polynomial_type;
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_block_size = 3;
    options.modular_composition_giant_step_group_size = 8;
    options.modular_composition_giant_step_prepared_byte_limit = 4 * 1024 * 1024;
    polynomial_arithmetic::polynomial_context<backend_type> context {backend_type(32), options};
    polynomial_type divisor(8, fq12_value_type::zero());
    divisor.back() = fq12_value_type(3);    // 3*X^7: powers of X eventually vanish.
    math::polynomial_divisor_context<backend_type> fixed(divisor, 6, context);
    const auto outer = grouped_polynomial(60, 11);
    for (const polynomial_type inner :
         {polynomial_type {fq12_value_type::zero()}, polynomial_type {fq12_value_type::one()},
          polynomial_type {fq12_value(3)}, polynomial_type {fq12_value_type::zero(), fq12_value_type::one()}}) {
        polynomial_type expected;
        math::compose_mod_reference(expected, outer, inner, fixed, context);
        auto actual = inner;
        math::compose_mod(actual, outer, actual, fixed, context);
        BOOST_CHECK(actual == expected);
        actual = outer;
        math::compose_mod(actual, actual, inner, fixed, context);
        BOOST_CHECK(actual == expected);
    }
    math::polynomial_divisor_context<backend_type> constant(polynomial_type {fq12_value_type(3)}, 1, context);
    math::polynomial_composition_precomputation<backend_type> cache(outer, outer.size(), constant, context);
    BOOST_CHECK_EQUAL(cache.giant_step_group_size(), 0);
    polynomial_type actual;
    math::compose_mod(actual, outer, cache, constant, context);
    BOOST_CHECK(actual == polynomial_type {fq12_value_type::zero()});
}

BOOST_AUTO_TEST_CASE(coefficient_preparation_panels_match_direct_sums_with_padding_and_tails) {
    using polynomial_type = math::polynomial<fq12_value_type>;
    using cache_type = math::detail::prepared_composition_coefficients<fq12_value_type>;
    for (std::size_t width : {3u, 8u, 11u, 17u}) {
        for (std::size_t terms : {1u, 8u, 12u}) {
            std::vector<polynomial_type> powers {polynomial_type {fq12_value_type::one()}};
            polynomial_type outer(terms + 1, fq12_value_type::zero());
            outer[0] = fq12_value(19);
            polynomial_type expected(width, fq12_value_type::zero());
            expected[0] = outer[0];
            for (std::size_t exponent = 1; exponent <= terms; ++exponent) {
                const std::size_t length = exponent == 1 ? width : width - exponent % 3;
                powers.emplace_back(length, fq12_value_type::zero());
                outer[exponent] = fq12_value(3 * exponent + 1);
                for (std::size_t coefficient = 0; coefficient < length; ++coefficient) {
                    powers.back()[coefficient] = fq12_value(5 * coefficient + exponent);
                    expected[coefficient] += powers.back()[coefficient] * outer[exponent];
                }
            }
            math::condense(expected);
            for (std::size_t lanes : {0u, 4u, 8u}) {
                const cache_type serial(powers, 4 * 1024 * 1024, lanes);
                for (std::size_t threads : {1u, 3u, 8u}) {
                    preparation_execution_state state;
                    state.threads = threads;
                    const cache_type parallel(powers, 4 * 1024 * 1024, lanes, [&](const auto &preparation) {
                        BOOST_CHECK_EQUAL(preparation.coefficient_count(), width);
                        BOOST_CHECK_EQUAL(preparation.term_count(), terms);
                        run_preparation_panels(preparation, state);
                    });
                    BOOST_CHECK_EQUAL(state.calls, 1);
                    BOOST_CHECK_EQUAL(state.active_workers.load(), 0);
                    BOOST_CHECK_EQUAL(parallel.bytes(), serial.bytes());
                    BOOST_CHECK_EQUAL(parallel.simd_lanes(), serial.simd_lanes());
                    std::array<polynomial_type, cache_type::block_tile_size> output;
                    cache_type::workspace scratch;
                    parallel.evaluate_blocks(output, outer, 0, 1, scratch);
                    BOOST_CHECK(output[0] == expected);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(coefficient_preparation_hook_preserves_fallback_and_joins_on_failure) {
    using backend_type = preparation_backend<fq12_value_type>;
    using polynomial_type = backend_type::polynomial_type;
    using cache_type = math::polynomial_composition_precomputation<backend_type>;
    const auto inner = grouped_polynomial(5, 7);
    const auto modulus = grouped_polynomial(14, 3);
    const auto outer = grouped_polynomial(37, 5);
    const auto saved_inner = inner;
    const auto saved_modulus = modulus;
    for (std::size_t limit : {0u, 1u, 1024u * 1024u}) {
        preparation_execution_state state;
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_block_size = 17;
        options.modular_composition_prepared_byte_limit = limit;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type(state), options);
        math::polynomial_divisor_context<backend_type> divisor(modulus, modulus.size(), context);
        polynomial_type expected;
        math::compose_mod_reference(expected, outer, inner, divisor, context);
        const cache_type cache(inner, outer.size(), divisor, context);
        BOOST_CHECK_EQUAL(state.calls, limit > 1 ? 1 : 0);
        polynomial_type result;
        math::compose_mod(result, outer, cache, divisor, context);
        BOOST_CHECK(result == expected);

        if (limit > 1) {
            state.fail = true;
            BOOST_CHECK_THROW(cache_type(inner, outer.size(), divisor, context), std::runtime_error);
            BOOST_CHECK_EQUAL(state.active_workers.load(), 0);
            state.fail = false;
            const cache_type retry(inner, outer.size(), divisor, context);
            math::compose_mod(result, outer, retry, divisor, context);
            BOOST_CHECK(result == expected);
            // A failed construction cannot invalidate an earlier immutable cache.
            math::compose_mod(result, outer, cache, divisor, context);
            BOOST_CHECK(result == expected);
        }
        BOOST_CHECK(inner == saved_inner);
        BOOST_CHECK(modulus == saved_modulus);
    }
}

BOOST_AUTO_TEST_CASE(coefficient_preparation_skips_unavailable_and_unsupported_caches) {
    using polynomial_type = math::polynomial<fq12_value_type>;
    using cache_type = math::detail::prepared_composition_coefficients<fq12_value_type>;
    std::size_t calls = 0;
    const auto dispatch = [&](const auto &preparation) {
        ++calls;
        preparation.prepare_panels(0, preparation.panel_count());
    };
    const std::vector<polynomial_type> powers {polynomial_type {fq12_value_type::one()}, grouped_polynomial(11, 3),
                                               grouped_polynomial(4, 7)};
    for (std::size_t limit : {0u, 1u}) {
        const cache_type cache(powers, limit, 8, dispatch);
        BOOST_CHECK(!cache.available());
    }
    const cache_type empty(std::vector<polynomial_type> {}, 1024, 8, dispatch);
    const cache_type constant(std::vector<polynomial_type> {powers[0]}, 1024, 8, dispatch);
    BOOST_CHECK(!empty.available());
    BOOST_CHECK(!constant.available());
    BOOST_CHECK_EQUAL(calls, 0);
    BOOST_CHECK_THROW(cache_type(powers, 0, 2, dispatch), std::invalid_argument);

    using backend_type = preparation_backend<fq_value_type>;
    preparation_execution_state state;
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_prepared_byte_limit = 1024 * 1024;
    options.modular_composition_block_size = 9;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type(state), options);
    const backend_type::polynomial_type modulus {fq_value_type(1), fq_value_type(2), fq_value_type(1)};
    math::polynomial_divisor_context<backend_type> divisor(modulus, 3, context);
    const math::polynomial_composition_precomputation<backend_type> cache({fq_value_type(3), fq_value_type(1)}, 17,
                                                                          divisor, context);
    BOOST_CHECK_EQUAL(cache.prepared_coefficient_bytes(), 0);
    BOOST_CHECK_EQUAL(state.calls, 0);
}

BOOST_AUTO_TEST_SUITE_END()
