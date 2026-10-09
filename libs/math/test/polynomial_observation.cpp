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

#define BOOST_TEST_MODULE polynomial_observation_test

#include <algorithm>
#include <array>
#include <limits>
#include <vector>
#include <stdexcept>
#include <type_traits>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/factorization/complete_factorization.hpp>
#include <nil/crypto3/random/algebraic_engine.hpp>
#include <nil/crypto3/math/polynomial/reconstruction/polynomial_x_norm_reconstruction.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace pa = math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using stage = pa::polynomial_stage;
    using outcome = pa::polynomial_scope_outcome;
    using reason = pa::polynomial_rejection_reason;
    using metric = pa::polynomial_metric;
    using field_type = fields::alt_bn128_base_field<254>;
    using value_type = field_type::value_type;
    using backend_type = pa::schoolbook_backend<value_type>;
    using polynomial_type = backend_type::polynomial_type;

    enum class event_kind { begin, end, progress };

    struct recorded_event {
        event_kind kind;
        pa::polynomial_scope_token scope;
        stage operation;
        outcome result = outcome::completed;
        reason rejection = reason::unspecified;
        std::size_t completed = 0;
        std::size_t total = 0;
        pa::polynomial_metadata metadata = {};
    };

    // Bounded test sink: recording overflow does not interfere with nesting or balanced callbacks.
    // No Boost assertions, allocation, timing or stream operations occur inside noexcept callbacks.
    template<bool Algorithms = true, bool Kernels = true, std::size_t Capacity = 32, bool Detailed = true>
    struct recording_observer {
        static constexpr bool observe_algorithms = Algorithms;
        static constexpr bool observe_kernels = Kernels;

        std::array<recorded_event, Capacity> events {};
        std::array<std::uint64_t, 32> parents {};
        std::size_t event_count = 0;
        std::size_t depth = 0;
        std::uint64_t next_id = 1;
        bool balanced = true;

        recording_observer() = default;
        recording_observer(const recording_observer &) = delete;
        recording_observer &operator=(const recording_observer &) = delete;

        void record(recorded_event event) noexcept {
            if constexpr (!Detailed) {
                // Pipeline tests retain high-level records only. All callbacks still check nesting, including
                // dropped arithmetic scopes and the thousands of exponent steps for extension-field orders.
                if ((event.operation > stage::invocation && event.operation < stage::powmod) ||
                    event.operation == stage::multiply_by_x ||
                    (event.kind == event_kind::progress &&
                     (event.operation == stage::powmod || event.operation == stage::powmod_x))) {
                    return;
                }
            }
            if (event_count < events.size()) {
                events[event_count] = event;
            }
            ++event_count;
        }

        pa::polynomial_scope_token on_begin(const pa::polynomial_scope_begin &event) noexcept {
            pa::polynomial_scope_token token {next_id++, depth == 0 ? 0 : parents[depth - 1]};
            parents[depth++] = token.id;
            record(
                {event_kind::begin, token, event.stage, outcome::completed, reason::unspecified, 0, 0, event.metadata});
            return token;
        }

        void on_end(const pa::polynomial_scope_end &event) noexcept {
            balanced &= depth != 0 && parents[depth - 1] == event.scope.id;
            if (depth != 0) {
                --depth;
            }
            balanced &= event.scope.parent_id == (depth == 0 ? 0 : parents[depth - 1]);
            record({event_kind::end, event.scope, event.stage, event.outcome, event.rejection, 0, 0, event.metadata});
        }

        void on_progress(const pa::polynomial_scope_progress &event) noexcept {
            balanced &= depth != 0 && parents[depth - 1] == event.scope.id;
            record({event_kind::progress, event.scope, event.stage, outcome::completed, reason::unspecified,
                    event.completed, event.total, event.metadata});
        }
    };

    using observer_type = recording_observer<>;
    using observed_context = pa::polynomial_context<backend_type, observer_type>;
    using plain_context = pa::polynomial_context<backend_type>;

    struct disabled_observer {
        static constexpr bool observe_algorithms = false;
        static constexpr bool observe_kernels = false;
        // Deliberately no callbacks.
    };

    struct throwing_observer : recording_observer<> {
        void on_end(const pa::polynomial_scope_end &) {
        }
    };

    struct incomplete_observer {
        static constexpr bool observe_algorithms = true;
        static constexpr bool observe_kernels = false;
    };

    static_assert(pa::PolynomialObserver<observer_type>);
    static_assert(pa::PolynomialObserver<disabled_observer>);
    static_assert(!pa::PolynomialObserver<throwing_observer>);
    static_assert(!pa::PolynomialObserver<incomplete_observer>);
    static_assert(!pa::PolynomialObserver<int>);
    static_assert(std::is_default_constructible_v<plain_context>);
    static_assert(std::is_same_v<decltype(pa::polynomial_context(backend_type {})), plain_context>);
    static_assert(std::is_same_v<decltype(pa::polynomial_context(backend_type {}, pa::polynomial_context_options {})),
                                 plain_context>);
    static_assert(!std::is_default_constructible_v<observed_context>);
    static_assert(!std::is_constructible_v<observed_context, backend_type>);
    static_assert(!std::is_constructible_v<observed_context, backend_type, observer_type &&>);
    static_assert(std::is_copy_constructible_v<observed_context>);
    using disabled_scope = decltype(std::declval<plain_context &>().observe<stage::invocation>());
    using enabled_scope = decltype(std::declval<observed_context &>().observe<stage::invocation>());
    static_assert(std::is_empty_v<disabled_scope>);
    static_assert(std::is_trivially_destructible_v<disabled_scope>);
    static_assert(std::is_nothrow_destructible_v<enabled_scope>);
    static_assert(!std::is_move_constructible_v<enabled_scope>);
    static_assert(!std::is_copy_constructible_v<enabled_scope>);

    template<typename Context>
    concept HasPreparedOperations = requires(Context &context, const typename Context::polynomial_type &p) {
        context.prepare_low_product(p, 5, 6);
        context.prepare_cyclic_remainder(p, 5);
    };

    using mixed_backend = pa::mixed_radix_backend<field_type>;
    using observed_mixed_context = pa::polynomial_context<mixed_backend, observer_type>;
    static_assert(HasPreparedOperations<observed_mixed_context>);
    static_assert(!HasPreparedOperations<observed_context>);

    template<bool Algorithms, bool Kernels>
    void check_levels() {
        recording_observer<Algorithms, Kernels> observer;
        pa::polynomial_context<backend_type, decltype(observer)> context(backend_type {}, observer);
        std::size_t root_metadata = 0, algorithm_metadata = 0, kernel_metadata = 0;
        {
            auto root = context.template observe<stage::invocation>([&]() noexcept {
                ++root_metadata;
                return pa::polynomial_metadata {};
            });
            {
                auto algorithm = context.template observe<stage::complete_factorization>([&]() noexcept {
                    ++algorithm_metadata;
                    return pa::polynomial_metadata {};
                });
                auto kernel = context.template observe<stage::fft_forward>([&]() noexcept {
                    ++kernel_metadata;
                    return pa::polynomial_metadata {};
                });
            }
        }
        BOOST_CHECK_EQUAL(root_metadata, Algorithms || Kernels);
        BOOST_CHECK_EQUAL(algorithm_metadata, Algorithms);
        BOOST_CHECK_EQUAL(kernel_metadata, Kernels);
        BOOST_CHECK_EQUAL(observer.event_count, 2 * ((Algorithms || Kernels) + Algorithms + Kernels));
        if constexpr (Algorithms || Kernels) {
            BOOST_CHECK_EQUAL(observer.events[1].scope.parent_id, observer.events[0].scope.id);
        }
        BOOST_CHECK(observer.balanced);
        BOOST_CHECK_EQUAL(observer.depth, 0);
    }

    // A final custom backend with no observer API, and a distinguishable dedicated square entry point.
    struct backend_counts {
        std::size_t products = 0, squares = 0, low_products = 0;
        bool throw_on_square = false;
        std::vector<std::array<std::size_t, 4>> operations;
    };

    template<typename Value>
    struct counting_backend final {
        using polynomial_type = math::polynomial<Value>;
        backend_counts &counts;
        pa::schoolbook_backend<Value> backend;

        void multiply(polynomial_type &output, const polynomial_type &a, const polynomial_type &b) {
            ++counts.products;
            counts.operations.push_back({0, a.size(), b.size(), 0});
            backend.multiply(output, a, b);
        }
        void square(polynomial_type &output, const polynomial_type &a) {
            ++counts.squares;
            counts.operations.push_back({1, a.size(), 0, 0});
            if (counts.throw_on_square) {
                throw std::domain_error("backend square failure");
            }
            backend.square(output, a);
        }
        void multiply_low(polynomial_type &output, const polynomial_type &a, const polynomial_type &b,
                          std::size_t count) {
            ++counts.low_products;
            counts.operations.push_back({2, a.size(), b.size(), count});
            backend.multiply_low(output, a, b, count);
        }
    };

    using custom_backend = counting_backend<value_type>;

    template<typename Backend>
    void check_arithmetic(Backend plain_backend, Backend observed_backend, const typename Backend::polynomial_type &p) {
        recording_observer<true, false> observer;
        pa::polynomial_context<Backend> plain(std::move(plain_backend));
        pa::polynomial_context<Backend, decltype(observer)> observed(std::move(observed_backend), observer);
        typename Backend::polynomial_type expected, actual;
        {
            auto root = observed.template observe<stage::invocation>();
            plain.multiply(expected, p, p);
            observed.multiply(actual, p, p);
            BOOST_CHECK(actual == expected);
            plain.square(expected, p);
            actual = p;
            observed.square(actual, actual);
            BOOST_CHECK(actual == expected);
            plain.multiply_low(expected, p, p, 2);
            actual = p;
            observed.multiply_low(actual, actual, actual, 2);
            BOOST_CHECK(actual == expected);
        }
        BOOST_REQUIRE_EQUAL(observer.event_count, 2);
        BOOST_CHECK(observer.balanced);
    }
}    // namespace

BOOST_AUTO_TEST_CASE(disabled_observation_has_empty_storage_and_never_evaluates_metadata) {
    pa::polynomial_context<backend_type, disabled_observer> context;
    struct original_context_layout {
        backend_type backend;
        pa::polynomial_context_options options;
    };
    static_assert(sizeof(plain_context) == sizeof(original_context_layout));
    std::size_t evaluations = 0;
    auto metadata = [&]() noexcept {
        ++evaluations;
        return pa::polynomial_metadata {};
    };
    {
        auto scope = context.observe<stage::invocation>(metadata);
        scope.set_result_metadata(metadata);
        scope.progress(1, 2, metadata);
        scope.progress(2, 2);
        scope.reject(reason::not_square);
        scope.callback_stop();
        scope.decline();
        disabled_observer observer;
        auto kernel = pa::observe_polynomial<stage::fft_forward>(observer, metadata);
        kernel.set_result_metadata(metadata);
        kernel.decline();
    }
    BOOST_CHECK_EQUAL(evaluations, 0);
}

BOOST_AUTO_TEST_CASE(algorithm_and_kernel_selection_are_independent) {
    check_levels<false, false>();
    check_levels<true, false>();
    check_levels<false, true>();
    check_levels<true, true>();
}

BOOST_AUTO_TEST_CASE(borrowed_observer_preserves_hierarchy_across_context_copies_and_progress) {
    observer_type observer;
    observed_context context(backend_type {}, observer);
    auto copy = context;
    {
        auto root = context.observe<stage::invocation>(
            []() noexcept { return pa::polynomial_metadata {{{metric::invocation_label, 42}}}; });
        for (std::size_t i = 0; i < 2; ++i) {
            {
                auto child = copy.observe<stage::gcd>();
                child.set_result_metadata(
                    []() noexcept { return pa::polynomial_metadata {{{metric::output_coefficients, 1}}}; });
            }
            // Progress follows the completed child, so its work is attributed before the heartbeat.
            root.progress(i + 1, 2);
        }
    }
    BOOST_REQUIRE_EQUAL(observer.event_count, 8);
    const auto &events = observer.events;
    BOOST_CHECK(events[0].kind == event_kind::begin);
    BOOST_CHECK_EQUAL(events[0].metadata[0].value, 42);
    BOOST_CHECK_EQUAL(events[0].scope.parent_id, 0);
    BOOST_CHECK_EQUAL(events[1].scope.parent_id, events[0].scope.id);
    BOOST_CHECK_EQUAL(events[2].scope.id, events[1].scope.id);
    BOOST_CHECK(events[2].kind == event_kind::end);
    BOOST_CHECK(events[2].metadata[0].metric == metric::output_coefficients);
    BOOST_CHECK_EQUAL(events[2].metadata[0].value, 1);
    BOOST_CHECK(events[3].kind == event_kind::progress);
    BOOST_CHECK_EQUAL(events[3].scope.id, events[0].scope.id);
    BOOST_CHECK_EQUAL(events[3].completed, 1);
    BOOST_CHECK_EQUAL(events[3].total, 2);
    BOOST_CHECK_NE(events[4].scope.id, events[1].scope.id);
    BOOST_CHECK_EQUAL(events[4].scope.parent_id, events[0].scope.id);
    BOOST_CHECK_EQUAL(events[6].completed, 2);
    BOOST_CHECK_EQUAL(events[7].scope.id, events[0].scope.id);
    BOOST_CHECK(observer.balanced);
    BOOST_CHECK_EQUAL(observer.depth, 0);
}

BOOST_AUTO_TEST_CASE(rejection_and_callback_stop_are_distinct_from_successful_completion) {
    observer_type observer;
    observed_context context(backend_type {}, observer);
    auto reject = [&] {
        auto scope = context.observe<stage::x_norm_recovery>();
        scope.reject(reason::not_square);
        return false;
    };
    auto stop = [&] {
        auto scope = context.observe<stage::complete_factorization>();
        scope.callback_stop();
        return false;
    };
    BOOST_CHECK(!reject());
    BOOST_CHECK(!stop());
    {
        auto scope = context.observe<stage::invocation>();
    }
    BOOST_REQUIRE_EQUAL(observer.event_count, 6);
    BOOST_CHECK(observer.events[1].result == outcome::rejected);
    BOOST_CHECK(observer.events[1].rejection == reason::not_square);
    BOOST_CHECK(observer.events[3].result == outcome::callback_stopped);
    BOOST_CHECK(observer.events[3].rejection == reason::unspecified);
    BOOST_CHECK(observer.events[5].result == outcome::completed);
    BOOST_CHECK(observer.balanced);
    BOOST_CHECK_EQUAL(observer.depth, 0);
}

BOOST_AUTO_TEST_CASE(unwinding_closes_every_scope_and_preserves_the_backend_exception) {
    observer_type observer;
    backend_counts counts;
    counts.throw_on_square = true;
    pa::polynomial_context<custom_backend, observer_type> context(custom_backend {counts, {}}, observer);
    polynomial_type output;
    try {
        auto outer = context.observe<stage::invocation>();
        auto inner = context.observe<stage::square>();
        inner.reject(reason::not_square);
        context.square(output, polynomial_type {1, 2});
        BOOST_FAIL("expected the backend exception");
    } catch (const std::domain_error &error) {
        BOOST_CHECK_EQUAL(error.what(), "backend square failure");
    }
    BOOST_REQUIRE_EQUAL(observer.event_count, 6);
    for (std::size_t i = 3; i < 6; ++i) {
        BOOST_CHECK(observer.events[i].result == outcome::exception);
        BOOST_CHECK(observer.events[i].rejection == reason::unspecified);
    }
    BOOST_CHECK(observer.balanced);
    BOOST_CHECK_EQUAL(observer.depth, 0);
}

BOOST_AUTO_TEST_CASE(caught_exceptions_and_scopes_opened_during_unwinding_do_not_report_false_failures) {
    observer_type observer;
    observed_context context(backend_type {}, observer);
    struct observe_cleanup {
        observed_context &context;
        ~observe_cleanup() noexcept {
            auto scope = context.observe<stage::gcd>();
        }
    };
    {
        auto root = context.observe<stage::invocation>();
        try {
            observe_cleanup cleanup {context};
            throw std::domain_error("caught");
        } catch (const std::domain_error &) {
        }
    }
    BOOST_REQUIRE_EQUAL(observer.event_count, 4);
    BOOST_CHECK(observer.events[2].result == outcome::completed);
    BOOST_CHECK(observer.events[3].result == outcome::completed);
    BOOST_CHECK(observer.balanced);
    BOOST_CHECK_EQUAL(observer.depth, 0);
}

BOOST_AUTO_TEST_CASE(dropping_records_does_not_require_unbounded_event_storage) {
    observer_type observer;
    observed_context context(backend_type {}, observer);
    {
        auto root = context.observe<stage::invocation>();
        for (std::size_t i = 0; i < 1000; ++i) {
            auto child = context.observe<stage::multiply>();
        }
    }
    BOOST_CHECK_EQUAL(observer.event_count, 2002);
    BOOST_CHECK_EQUAL(observer.next_id, 1002);
    BOOST_CHECK(observer.balanced);
    BOOST_CHECK_EQUAL(observer.depth, 0);
}

BOOST_AUTO_TEST_CASE(observed_context_preserves_schoolbook_and_extension_field_arithmetic) {
    check_arithmetic(backend_type {}, backend_type {}, polynomial_type {1, 2, 3});
    using extension_field = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using extension_value = extension_field::value_type;
    using extension_backend = pa::mixed_radix_backend<field_type, extension_value>;
    auto a = extension_value::zero(), b = extension_value::zero();
    for (std::size_t i = 0; i < extension_field::arity; ++i) {
        a.coordinate(i) = value_type(i + 1);
        b.coordinate(i) = value_type(2 * i + 1);
    }
    check_arithmetic(extension_backend(9), extension_backend(9), extension_backend::polynomial_type {a, b, a});
}

BOOST_AUTO_TEST_CASE(custom_backend_needs_no_observer_methods_and_retains_dedicated_squaring) {
    observer_type observer;
    backend_counts counts;
    pa::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 7;
    pa::polynomial_context context(custom_backend {counts, {}}, observer, options);
    const polynomial_type input {1, 2};
    polynomial_type output;
    context.multiply(output, input, input);
    context.square(output, input);
    context.multiply_low(output, input, input, 2);
    BOOST_CHECK_EQUAL(counts.products, 1);
    BOOST_CHECK_EQUAL(counts.squares, 1);
    BOOST_CHECK_EQUAL(counts.low_products, 1);
    BOOST_CHECK_EQUAL(context.options().basecase_divisor_coefficient_cutoff, 7);
    BOOST_CHECK_EQUAL(context.options().gcd_half_gcd_cutoff, 0);
    // Context scopes are available even when a backend has no detailed instrumentation overloads.
    BOOST_REQUIRE_EQUAL(observer.event_count, 6);
    BOOST_CHECK(observer.events[0].operation == stage::multiply);
    BOOST_CHECK(observer.events[2].operation == stage::square);
    BOOST_CHECK(observer.events[4].operation == stage::multiply_low);
    for (std::size_t i : {1, 3, 5}) {
        BOOST_CHECK(observer.events[i].metadata[1].metric == metric::detailed_backend);
        BOOST_CHECK_EQUAL(observer.events[i].metadata[1].value, 0);
    }
    BOOST_CHECK(observer.balanced);
}

BOOST_AUTO_TEST_CASE(observed_context_keeps_prepared_fast_paths_and_declines_available) {
    recording_observer<true, true, 1024> observer;
    pa::polynomial_context context(mixed_backend(18), observer);
    backend_type reference;
    const polynomial_type divisor {2, 3, 5, 7, 11, 13, 17}, quotient {2, 3, 5, 7, 11};
    polynomial_type dividend, output, expected;
    const auto low = context.prepare_low_product(divisor, quotient.size(), 6);
    BOOST_REQUIRE(low.has_value());
    reference.multiply_low(expected, quotient, divisor, 6);
    BOOST_REQUIRE(context.try_multiply_low_prepared(output, quotient, *low, 6));
    BOOST_CHECK(output == expected);
    BOOST_CHECK(!context.try_multiply_low_prepared(output, quotient, *low, 7));
    BOOST_CHECK(output == expected);
    const auto cyclic = context.prepare_cyclic_remainder(divisor, quotient.size());
    BOOST_REQUIRE(cyclic.has_value());
    reference.multiply(dividend, divisor, quotient);
    dividend[0] += value_type(19);
    output = dividend;
    BOOST_REQUIRE(context.try_cyclic_remainder(output, output, quotient, *cyclic));
    BOOST_CHECK(output == polynomial_type({value_type(19)}));
    BOOST_CHECK(!context.prepare_cyclic_remainder(divisor, 1).has_value());
    BOOST_CHECK(!context.try_cyclic_remainder(output, dividend, polynomial_type {1}, *cyclic));
    BOOST_CHECK(output == polynomial_type({value_type(19)}));
}

namespace {
    // Larger bounded recording only for small pipeline fixtures; the production observer remains caller-owned.
    using arithmetic_observer = recording_observer<true, false, 2048>;
    using path = pa::polynomial_arithmetic_path;

    std::size_t measurement(const recorded_event &event, metric key) {
        for (const auto &entry : event.metadata) {
            if (entry.metric == key) {
                return entry.value;
            }
        }
        return std::numeric_limits<std::size_t>::max();
    }

    template<typename Observer>
    std::size_t count_events(const Observer &observer, stage operation, event_kind kind = event_kind::end) {
        BOOST_REQUIRE_LE(observer.event_count, observer.events.size());
        return std::count_if(observer.events.begin(), observer.events.begin() + observer.event_count,
                             [&](const auto &event) { return event.operation == operation && event.kind == kind; });
    }

    template<typename Observer>
    bool reported_path(const Observer &observer, stage operation, path selected) {
        BOOST_REQUIRE_LE(observer.event_count, observer.events.size());
        return std::any_of(observer.events.begin(), observer.events.begin() + observer.event_count,
                           [&](const auto &event) {
                               return event.operation == operation && event.kind == event_kind::end &&
                                      measurement(event, metric::arithmetic_path) == static_cast<std::size_t>(selected);
                           });
    }

    template<typename Backend>
    void check_division_observation(Backend backend, const typename Backend::polynomial_type &divisor,
                                    bool expect_cyclic) {
        using polynomial = typename Backend::polynomial_type;
        using value = typename polynomial::value_type;
        pa::polynomial_context_options options;
        options.basecase_divisor_coefficient_cutoff = 0;
        options.basecase_quotient_coefficient_cutoff = 0;
        pa::polynomial_context<Backend> plain(backend, options);
        arithmetic_observer observer;
        pa::polynomial_context<Backend, arithmetic_observer> observed(std::move(backend), observer, options);
        const math::polynomial_divisor_context<Backend> plain_divisor(divisor, 5, plain);
        const math::polynomial_divisor_context<Backend> observed_divisor(divisor, 5, observed);
        // Preparation does not borrow the observer: the same owned snapshot works with either context.
        pa::schoolbook_backend<value> reference;
        for (std::size_t quotient_size : {5, 4, 1}) {
            const polynomial known_quotient(quotient_size, value::one());
            polynomial dividend;
            reference.multiply(dividend, known_quotient, divisor);
            dividend[0] += value::one();
            polynomial expected_quotient, expected_remainder, actual_remainder;
            math::divrem(expected_quotient, expected_remainder, dividend, observed_divisor, plain);
            auto actual_quotient = dividend;
            math::divrem(actual_quotient, actual_remainder, actual_quotient, observed_divisor, observed);
            BOOST_CHECK(actual_quotient == known_quotient);
            BOOST_CHECK(actual_quotient == expected_quotient);
            BOOST_CHECK(actual_remainder == polynomial({value::one()}));
            BOOST_CHECK(actual_remainder == expected_remainder);
            // A snapshot prepared by an unobserved context must also work in an observed call.
            math::divrem(actual_quotient, actual_remainder, dividend, plain_divisor, observed);
            BOOST_CHECK(actual_quotient == expected_quotient);
            BOOST_CHECK(actual_remainder == expected_remainder);
        }
        BOOST_CHECK_EQUAL(count_events(observer, stage::divisor_preparation), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::inverse_series), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::division), 6);
        BOOST_CHECK_EQUAL(count_events(observer, stage::quotient_estimation), 6);
        BOOST_CHECK_EQUAL(count_events(observer, stage::remainder_reconstruction), 6);
        BOOST_CHECK_EQUAL(count_events(observer, stage::basecase_division), 0);
        BOOST_CHECK(reported_path(observer, stage::division, path::reciprocal));
        BOOST_CHECK(reported_path(observer, stage::remainder_reconstruction, path::low_product));
        BOOST_CHECK_EQUAL(reported_path(observer, stage::remainder_reconstruction, path::cyclic), expect_cyclic);
        BOOST_CHECK(observer.balanced);
        BOOST_CHECK_EQUAL(observer.depth, 0);
        std::size_t last_precision = 1;
        for (std::size_t i = 0; i < observer.event_count; ++i) {
            const auto &event = observer.events[i];
            if (event.kind == event_kind::progress && event.operation == stage::inverse_series) {
                BOOST_CHECK_GT(event.completed, last_precision);
                BOOST_CHECK_LE(event.completed, 5);
                BOOST_CHECK_EQUAL(event.total, 5);
                last_precision = event.completed;
            }
            if (event.kind == event_kind::end && event.operation == stage::division) {
                BOOST_CHECK(event.result == outcome::completed);
                BOOST_CHECK_EQUAL(measurement(event, metric::remainder_coefficients), 1);
            }
        }
        BOOST_CHECK_EQUAL(last_precision, 5);
    }

    template<typename Backend>
    void check_gcd_observation(Backend backend) {
        using polynomial = typename Backend::polynomial_type;
        const polynomial common {2, 1};
        polynomial left, right;
        backend.multiply(left, common, polynomial {2, 11, 7, 19, 3, 17, 23, 5, 13, 29, 31, 37, 41});
        backend.multiply(right, common, polynomial {7, 31, 3, 23, 5, 17, 43, 11, 19, 47, 29, 37});
        pa::polynomial_context<Backend> plain(backend);
        polynomial expected;
        math::gcd(expected, left, right, plain);
        BOOST_CHECK(expected == common);
        for (std::size_t half_gcd_cutoff : {0, 2}) {
            arithmetic_observer observer;
            pa::polynomial_context_options options;
            options.gcd_half_gcd_cutoff = half_gcd_cutoff;
            options.half_gcd_basecase_cutoff = 4;
            pa::polynomial_context<Backend, arithmetic_observer> context(backend, observer, options);
            polynomial result = left;
            math::gcd(result, result, right, context);
            BOOST_CHECK(result == expected);
            BOOST_CHECK_EQUAL(count_events(observer, stage::gcd), 1);
            BOOST_CHECK_GT(count_events(observer, stage::euclidean_step), 0);
            if (half_gcd_cutoff == 0) {
                BOOST_CHECK_EQUAL(count_events(observer, stage::half_gcd), 0);
                BOOST_CHECK_EQUAL(count_events(observer, stage::half_gcd_recursion), 0);
            } else {
                BOOST_CHECK_GT(count_events(observer, stage::half_gcd), 0);
                BOOST_CHECK_GT(count_events(observer, stage::half_gcd_recursion), 0);
                BOOST_CHECK_GT(count_events(observer, stage::half_gcd_basecase), 0);
            }
            for (std::size_t i = 0; i < observer.event_count; ++i) {
                const auto &event = observer.events[i];
                if (event.kind == event_kind::begin && event.operation == stage::euclidean_step) {
                    BOOST_CHECK_GE(measurement(event, metric::input_coefficients),
                                   measurement(event, metric::second_input_coefficients));
                }
                if (event.kind == event_kind::end && event.operation == stage::euclidean_step) {
                    BOOST_CHECK_NE(measurement(event, metric::quotient_coefficients),
                                   std::numeric_limits<std::size_t>::max());
                    BOOST_CHECK(event.result == outcome::completed);
                }
            }
            BOOST_CHECK(observer.balanced);
            BOOST_CHECK_EQUAL(observer.depth, 0);
        }
    }

    polynomial_type repeated_factor_input() {
        backend_type backend;
        polynomial_type input {value_type(7)};
        // 7 * (X + 1) * (X + 2)^2 * (X + 3)^3: three ordered, nonempty multiplicity groups.
        for (std::size_t multiplicity = 1; multiplicity <= 3; ++multiplicity) {
            for (std::size_t repeat = 0; repeat < multiplicity; ++repeat) {
                backend.multiply(input, input, polynomial_type {value_type(multiplicity), value_type::one()});
            }
        }
        return input;
    }
}    // namespace

BOOST_AUTO_TEST_CASE(explicit_phase_finish_is_idempotent_and_leaves_temporaries_alive) {
    observer_type observer;
    observed_context context(backend_type {}, observer);
    std::size_t ignored_metadata = 0;
    {
        auto root = context.observe<stage::invocation>();
        auto phase = context.observe<stage::gcd>();
        polynomial_type temporary {1, 2};
        phase.finish();
        phase.finish();
        phase.progress(1, 1, [&]() noexcept {
            ++ignored_metadata;
            return pa::polynomial_metadata {};
        });
        phase.set_result_metadata([&]() noexcept {
            ++ignored_metadata;
            return pa::polynomial_metadata {};
        });
        phase.reject(reason::not_square);
        BOOST_CHECK_EQUAL(temporary.size(), 2);
        auto next = context.observe<stage::division>();
    }
    BOOST_REQUIRE_EQUAL(observer.event_count, 6);
    BOOST_CHECK_EQUAL(observer.events[3].scope.parent_id, observer.events[0].scope.id);
    BOOST_CHECK(observer.events[2].result == outcome::completed);
    BOOST_CHECK_EQUAL(ignored_metadata, 0);
    BOOST_CHECK(observer.balanced);
}

BOOST_AUTO_TEST_CASE(observed_division_reuses_reciprocal_and_preserves_cyclic_and_fallback_paths) {
    check_division_observation(backend_type {}, polynomial_type {2, 3, 5, 7, 11, 13, 17}, false);
    using extension_field = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using extension_value = extension_field::value_type;
    using extension_backend = pa::mixed_radix_backend<field_type, extension_value>;
    extension_backend::polynomial_type divisor(7, extension_value::zero());
    for (std::size_t i = 0; i < divisor.size(); ++i) {
        for (std::size_t j = 0; j < extension_field::arity; ++j) {
            divisor[i].coordinate(j) = value_type(1 + 17 * i + j);
        }
    }
    check_division_observation(extension_backend(18), divisor, true);
}

BOOST_AUTO_TEST_CASE(observed_division_preserves_shortcuts_validation_and_exception_outcomes) {
    arithmetic_observer observer;
    pa::polynomial_context<backend_type, arithmetic_observer> context(backend_type {}, observer);
    const polynomial_type divisor {2, 3, 5, 7}, input {11, 13, 17, 19, 23, 29, 31};
    // Basecase selection still precedes precision validation.
    const math::polynomial_divisor_context<backend_type> prepared(divisor, 1, context);
    polynomial_type q, r, expected_q, expected_r;
    math::division(expected_q, expected_r, input, divisor);
    math::divrem(q, r, input, prepared, context);
    BOOST_CHECK(q == expected_q);
    BOOST_CHECK(r == expected_r);
    BOOST_CHECK(reported_path(observer, stage::division, path::basecase));
    math::divrem(q, r, polynomial_type {1, 2}, prepared, context);
    BOOST_CHECK(math::is_zero(q));
    BOOST_CHECK(r == polynomial_type({value_type(1), value_type(2)}));
    BOOST_CHECK(reported_path(observer, stage::division, path::trivial));
    BOOST_CHECK_THROW(math::divrem(q, q, input, prepared, context), std::invalid_argument);
    pa::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    pa::polynomial_context<backend_type, arithmetic_observer> forced(backend_type {}, observer, options);
    BOOST_CHECK_THROW(math::divrem(q, r, input, prepared, forced), std::invalid_argument);
    const math::polynomial_divisor_context<backend_type> constant(polynomial_type {value_type(2)}, 7, forced);
    math::divrem(q, r, input, constant, forced);
    BOOST_CHECK(math::is_zero(r));
    BOOST_CHECK(reported_path(observer, stage::remainder_reconstruction, path::trivial));
    const polynomial_type sentinel {value_type(99)};
    q = sentinel;
    BOOST_CHECK_THROW(math::exact_division(q, input, prepared, context), std::invalid_argument);
    BOOST_CHECK(q == sentinel);
    BOOST_CHECK_THROW(
        (math::polynomial_divisor_context<backend_type>(polynomial_type {value_type::zero()}, 1, context)),
        std::invalid_argument);
    BOOST_REQUIRE_LE(observer.event_count, observer.events.size());
    std::size_t exception_count = 0;
    for (std::size_t i = 0; i < observer.event_count; ++i) {
        const auto &event = observer.events[i];
        if (event.kind == event_kind::end) {
            BOOST_CHECK(event.result != outcome::rejected);
            exception_count += event.result == outcome::exception;
        }
    }
    BOOST_CHECK_EQUAL(exception_count, 4);
    BOOST_CHECK(observer.balanced);
    BOOST_CHECK_EQUAL(observer.depth, 0);
}

BOOST_AUTO_TEST_CASE(observed_gcd_preserves_classical_and_opt_in_half_gcd_results) {
    check_gcd_observation(backend_type {});
    check_gcd_observation(mixed_backend(18));
}

BOOST_AUTO_TEST_CASE(square_free_observation_preserves_callbacks_backend_calls_and_stop_boundaries) {
    const auto input = repeated_factor_input();
    using factor_type = math::polynomial_factor<polynomial_type>;
    for (int mode : {0, 1, 2}) {    // Complete, callback-requested stop, callback exception.
        backend_counts plain_counts, observed_counts;
        pa::polynomial_context<custom_backend> plain(custom_backend {plain_counts, {}});
        arithmetic_observer observer;
        pa::polynomial_context<custom_backend, arithmetic_observer> context(custom_backend {observed_counts, {}},
                                                                            observer);
        std::vector<factor_type> plain_factors, observed_factors;
        auto run = [&](auto &arithmetic, auto &factors) {
            return math::square_free_factorization(input, arithmetic, [&](const factor_type &factor) {
                factors.push_back(factor);
                // Callback arithmetic must be attributed beneath the callback, in both stopped and complete runs.
                polynomial_type common;
                math::gcd(common, factor.polynomial, factor.polynomial, arithmetic);
                if (factors.size() == 2) {
                    if (mode == 2) {
                        throw std::runtime_error("factor callback failure");
                    }
                    if (mode == 1) {
                        return math::factorization_control::stop_factorization;
                    }
                }
                return math::factorization_control::continue_factorization;
            });
        };
        if (mode == 2) {
            BOOST_CHECK_THROW(run(plain, plain_factors), std::runtime_error);
            BOOST_CHECK_THROW(run(context, observed_factors), std::runtime_error);
        } else {
            const auto expected = run(plain, plain_factors);
            const auto actual = run(context, observed_factors);
            BOOST_CHECK(actual == expected);
            BOOST_CHECK_EQUAL(actual.complete, mode == 0);
            BOOST_CHECK(actual.leading_coefficient == value_type(7));
        }
        BOOST_CHECK(observed_factors == plain_factors);
        BOOST_CHECK(observed_counts.operations == plain_counts.operations);
        BOOST_CHECK_EQUAL(count_events(observer, stage::square_free_factorization), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::square_free_normalization), 2);
        BOOST_CHECK_EQUAL(count_events(observer, stage::derivative), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::factorization_exact_quotient), mode == 0 ? 7 : 4);
        BOOST_CHECK_EQUAL(count_events(observer, stage::multiplicity, event_kind::progress), mode == 0 ? 3 : 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::factor_callback, event_kind::progress), mode == 0 ? 3 :
                                                                                                mode == 1 ? 2 :
                                                                                                            1);
        BOOST_REQUIRE_LE(observer.event_count, observer.events.size());
        const auto &end = observer.events[observer.event_count - 1];
        BOOST_CHECK(end.operation == stage::square_free_factorization);
        BOOST_CHECK(end.result == (mode == 0 ? outcome::completed :
                                   mode == 1 ? outcome::callback_stopped :
                                               outcome::exception));
        BOOST_CHECK_EQUAL(measurement(end, metric::factor_count), mode == 0 ? 3 : 2);
        for (std::size_t i = 0; i + 1 < observer.event_count; ++i) {
            const auto &event = observer.events[i];
            if (event.kind == event_kind::begin && event.operation == stage::factor_callback) {
                BOOST_CHECK(observer.events[i + 1].operation == stage::gcd);
                BOOST_CHECK_EQUAL(observer.events[i + 1].scope.parent_id, event.scope.id);
            }
        }
        BOOST_CHECK(observer.balanced);
        BOOST_CHECK_EQUAL(observer.depth, 0);
    }
}

BOOST_AUTO_TEST_CASE(square_free_constant_shortcuts_and_disabled_algorithm_hooks) {
    for (value_type constant : {value_type::zero(), value_type(7)}) {
        arithmetic_observer observer;
        pa::polynomial_context<backend_type, arithmetic_observer> context(backend_type {}, observer);
        const auto result = math::square_free_factorization(polynomial_type {constant}, context);
        BOOST_CHECK(result.leading_coefficient == constant);
        BOOST_CHECK(result.factors.empty());
        BOOST_CHECK_EQUAL(count_events(observer, stage::gcd), 0);
        BOOST_CHECK_EQUAL(count_events(observer, stage::divisor_preparation), 0);
        BOOST_CHECK_EQUAL(count_events(observer, stage::factor_callback), 0);
        BOOST_CHECK(observer.balanced);
    }
    const auto input = repeated_factor_input();
    plain_context plain;
    const auto expected = math::square_free_factorization(input, plain);
    pa::polynomial_context<backend_type, disabled_observer> disabled;
    BOOST_CHECK(math::square_free_factorization(input, disabled) == expected);
    recording_observer<false, true, 8192> kernels_only;
    pa::polynomial_context<backend_type, decltype(kernels_only)> context(backend_type {}, kernels_only);
    BOOST_CHECK(math::square_free_factorization(input, context) == expected);
    BOOST_REQUIRE_GT(kernels_only.event_count, 0);
    BOOST_REQUIRE_LE(kernels_only.event_count, kernels_only.events.size());
    for (std::size_t i = 0; i < kernels_only.event_count; ++i)
        BOOST_CHECK(kernels_only.events[i].operation >= stage::multiply);
    BOOST_CHECK(kernels_only.balanced);
    BOOST_CHECK_EQUAL(kernels_only.depth, 0);
}

namespace {
    using pipeline_observer = recording_observer<true, false, 2048, false>;

    template<typename Observer>
    void check_balanced(const Observer &observer) {
        BOOST_REQUIRE_LE(observer.event_count, observer.events.size());
        BOOST_CHECK(observer.balanced);
        BOOST_CHECK_EQUAL(observer.depth, 0);
    }

    template<typename Observer>
    std::size_t count_outcome(const Observer &observer, stage operation, outcome result,
                              reason rejection = reason::unspecified) {
        BOOST_REQUIRE_LE(observer.event_count, observer.events.size());
        return std::count_if(observer.events.begin(), observer.events.begin() + observer.event_count,
                             [&](const auto &event) {
                                 return event.kind == event_kind::end && event.operation == operation &&
                                        event.result == result && event.rejection == rejection;
                             });
    }

    template<typename Backend>
    void check_composition_observation(Backend backend, const typename Backend::polynomial_type &inner) {
        using polynomial = typename Backend::polynomial_type;
        using value = typename polynomial::value_type;
        pa::polynomial_context_options options;
        options.modular_composition_cached_power_limit = 2;
        options.basecase_divisor_coefficient_cutoff = 0;
        options.basecase_quotient_coefficient_cutoff = 0;
        pa::polynomial_context<Backend> plain(backend, options);
        arithmetic_observer observer;
        pa::polynomial_context<Backend, arithmetic_observer> observed(backend, observer, options);
        const polynomial h {value(1), value(2), value(0), value(1), value(1)};
        const polynomial outer {value(2), value(3), value(5), value(7), value(11), value(13), value(17)};
        const math::polynomial_divisor_context<Backend> divisor(h, 6, plain);
        const math::polynomial_composition_precomputation<Backend> cache(inner, outer.size(), divisor, observed);
        BOOST_CHECK_EQUAL(cache.block_size(), 2);
        BOOST_CHECK_EQUAL(count_events(observer, stage::composition_cache), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::composition_cached_powers), 1);
        // Prepared data owns no observer: it works across observed and unobserved contexts in both directions.
        polynomial expected, reference, actual = outer;
        math::compose_mod(expected, outer, cache, divisor, plain);
        math::compose_mod_reference(reference, outer, inner, divisor, plain);
        BOOST_CHECK(expected == reference);
        math::compose_mod(actual, actual, cache, divisor, observed);
        BOOST_CHECK(actual == expected);
        BOOST_CHECK_EQUAL(count_events(observer, stage::composition_linear_combination), 4);
        BOOST_CHECK_EQUAL(count_events(observer, stage::composition_giant_step), 3);
        BOOST_CHECK_EQUAL(count_events(observer, stage::compose_mod_cached, event_kind::progress), 4);
        BOOST_CHECK_EQUAL(count_events(observer, stage::composition_finalization), 1);
        const math::polynomial_composition_precomputation<Backend> plain_cache(inner, outer.size(), divisor, plain);
        math::compose_mod(actual, outer, plain_cache, divisor, observed);
        BOOST_CHECK(actual == expected);
        BOOST_CHECK_EQUAL(count_events(observer, stage::composition_cache), 1);
        actual = inner;
        math::compose_mod(actual, outer, actual, divisor, observed);
        BOOST_CHECK(actual == expected);
        BOOST_CHECK_EQUAL(count_events(observer, stage::composition_cache), 2);
        BOOST_CHECK_THROW(math::compose_mod(actual, polynomial(8, value::one()), cache, divisor, observed),
                          std::invalid_argument);
        BOOST_CHECK_EQUAL(count_outcome(observer, stage::compose_mod_cached, outcome::exception), 1);
        check_balanced(observer);
    }

    // Scripted zero/constant samples, one unsuccessful character trial, then a proper GCD factor of X^2 - 1.
    const std::array split_samples = {value_type(0), value_type(0), value_type(1),  value_type(0),
                                      value_type(3), value_type(1), -value_type(1), value_type(1)};
}    // namespace

BOOST_AUTO_TEST_CASE(completed_step_counters_are_disabled_or_stop_at_finish) {
    plain_context plain;
    auto disabled = plain.observe<stage::powmod>();
    disabled.advance();
    BOOST_CHECK_EQUAL(disabled.completed_steps(), 0);
    observer_type observer;
    observed_context context(backend_type {}, observer);
    {
        auto scope = context.observe<stage::powmod>();
        scope.advance(2);
        scope.advance(2);
        BOOST_CHECK_EQUAL(scope.completed_steps(), 2);
        scope.finish();
        scope.advance(2);
        BOOST_CHECK_EQUAL(scope.completed_steps(), 2);
    }
    BOOST_REQUIRE_EQUAL(observer.event_count, 4);
    BOOST_CHECK_EQUAL(observer.events[2].completed, 2);
    BOOST_CHECK_EQUAL(observer.events[2].total, 2);
    check_balanced(observer);
}

BOOST_AUTO_TEST_CASE(exponent_observation_preserves_operations_and_reports_completed_reductions) {
    for (bool specialized : {false, true}) {
        backend_counts expected_counts, actual_counts;
        pa::polynomial_context<custom_backend> plain(custom_backend {expected_counts, {}});
        arithmetic_observer observer;
        pa::polynomial_context<custom_backend, arithmetic_observer> observed(custom_backend {actual_counts, {}},
                                                                             observer);
        const polynomial_type h {2, 3, 0, 5};
        const math::polynomial_divisor_context<custom_backend> divisor(h, 2, plain);
        actual_counts = expected_counts;    // Compare operations after the same precomputation.
        polynomial_type expected {0, 1}, actual = expected;
        if (specialized) {
            math::powmod_x(expected, 13, divisor, plain);
            math::powmod_x(actual, 13, divisor, observed);
        } else {
            math::powmod(expected, expected, 13, divisor, plain);
            math::powmod(actual, actual, 13, divisor, observed);
        }
        BOOST_CHECK(actual == expected);
        BOOST_CHECK(actual_counts.operations == expected_counts.operations);
        BOOST_CHECK_EQUAL(count_events(observer, stage::squaremod), 3);
        BOOST_CHECK_EQUAL(count_events(observer, stage::multiply_by_x), specialized ? 2 : 0);
        BOOST_CHECK_EQUAL(count_events(observer, stage::mulmod), specialized ? 0 : 2);
        const auto operation = specialized ? stage::powmod_x : stage::powmod;
        BOOST_CHECK_EQUAL(count_events(observer, operation, event_kind::progress), 4);
        std::size_t progress = 0, completed_squares = 0;
        for (std::size_t i = 0; i < observer.event_count; ++i) {
            const auto &event = observer.events[i];
            if (event.kind == event_kind::end && event.operation == stage::squaremod)
                ++completed_squares;
            if (event.kind == event_kind::progress && event.operation == operation) {
                BOOST_CHECK_EQUAL(event.completed, ++progress);
                BOOST_CHECK_EQUAL(completed_squares, specialized ? progress - 1 : std::min(progress, std::size_t(3)));
                BOOST_CHECK_EQUAL(event.total, specialized ? 4 : 0);
            }
            if (event.kind == event_kind::end && event.operation == operation) {
                BOOST_CHECK_EQUAL(measurement(event, metric::exponent_bits), 4);
                BOOST_CHECK_EQUAL(measurement(event, metric::output_coefficients), actual.size());
            }
        }
        check_balanced(observer);
    }
}

BOOST_AUTO_TEST_CASE(exponent_shortcuts_large_integers_and_exceptions_remain_observable) {
    plain_context plain;
    arithmetic_observer observer;
    pa::polynomial_context<backend_type, arithmetic_observer> context(backend_type {}, observer);
    const math::polynomial_divisor_context<backend_type> linear(polynomial_type {2, 1}, 1, plain);
    const math::polynomial_divisor_context<backend_type> constant(polynomial_type {value_type(7)}, 1, plain);
    polynomial_type result, expected;
    math::powmod_x(result, 0, constant, context);
    BOOST_CHECK(result == polynomial_type({value_type(0)}));
    math::powmod_x(result, 0, linear, context);
    BOOST_CHECK(result == polynomial_type({value_type(1)}));
    math::powmod_x(result, 1, linear, context);
    BOOST_CHECK(result == polynomial_type({-value_type(2)}));
    BOOST_CHECK_THROW(math::powmod_x(result, -1, linear, context), std::invalid_argument);
    BOOST_CHECK_THROW(math::powmod(result, polynomial_type {0, 1}, -1, linear, context), std::invalid_argument);
    const boost::multiprecision::cpp_int exponent = (boost::multiprecision::cpp_int(1) << 70) + 3;
    math::powmod_x(result, exponent, linear, context);
    math::powmod(expected, polynomial_type {0, 1}, exponent, linear, plain);
    BOOST_CHECK(result == expected);
    BOOST_CHECK_EQUAL(count_events(observer, stage::powmod_x, event_kind::progress), 72);    // 1 bit + 71 bits.
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::powmod_x, outcome::exception), 1);
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::powmod, outcome::exception), 1);
    check_balanced(observer);

    backend_counts counts;
    arithmetic_observer failing_observer;
    pa::polynomial_context<custom_backend, arithmetic_observer> failing(custom_backend {counts, {}}, failing_observer);
    const math::polynomial_divisor_context<custom_backend> divisor(polynomial_type {1, 1, 1}, 1, failing);
    counts.throw_on_square = true;
    BOOST_CHECK_THROW(math::powmod_x(result, 3, divisor, failing), std::domain_error);
    BOOST_CHECK_EQUAL(count_events(failing_observer, stage::powmod_x, event_kind::progress), 1);
    BOOST_CHECK_EQUAL(count_outcome(failing_observer, stage::squaremod, outcome::exception), 1);
    BOOST_CHECK_EQUAL(count_outcome(failing_observer, stage::powmod_x, outcome::exception), 1);
    check_balanced(failing_observer);
}

BOOST_AUTO_TEST_CASE(composition_observation_covers_schoolbook_and_mixed_extension_backends) {
    check_composition_observation(backend_type {}, polynomial_type {2, 3, 5});
    using extension_field = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using extension_backend = pa::mixed_radix_backend<field_type, extension_field::value_type>;
    auto coefficient = extension_field::value_type::zero();
    for (std::size_t i = 0; i < extension_field::arity; ++i)
        coefficient.coordinate(i) = value_type(i + 1);
    check_composition_observation(
        extension_backend(18),
        extension_backend::polynomial_type {coefficient, coefficient + coefficient, coefficient});
}

BOOST_AUTO_TEST_CASE(frobenius_observation_preserves_cached_and_uncached_composition_paths) {
    backend_counts plain_counts, observed_counts;
    pa::polynomial_context<custom_backend> plain(custom_backend {plain_counts, {}});
    pipeline_observer observer;
    pa::polynomial_context<custom_backend, pipeline_observer> observed(custom_backend {observed_counts, {}}, observer);
    const polynomial_type h {3, 1, 0, 0, 1};
    const math::polynomial_frobenius_context<custom_backend> expected_context(h, plain);
    const math::polynomial_frobenius_context<custom_backend> context(h, observed);
    BOOST_CHECK(expected_context.x_to_field_order() == context.x_to_field_order());
    BOOST_CHECK_EQUAL(count_events(observer, stage::initial_frobenius_power), 1);
    BOOST_CHECK_EQUAL(count_events(observer, stage::powmod_x), 1);
    BOOST_CHECK_EQUAL(count_events(observer, stage::composition_cache), 1);
    polynomial_type expected, actual {2, 3, 5};
    math::frobenius_map(expected, actual, 2, expected_context, plain);
    math::frobenius_map(actual, actual, 2, context, observed);
    BOOST_CHECK(actual == expected);
    BOOST_CHECK_EQUAL(count_events(observer, stage::composition_cache), 1);
    BOOST_CHECK_EQUAL(count_events(observer, stage::frobenius_iterations, event_kind::progress), 2);
    const polynomial_type unreduced {1, 2, 3, 4, 5, 6};
    math::frobenius_map(expected, unreduced, expected_context, plain);
    math::frobenius_map(actual, unreduced, context, observed);
    BOOST_CHECK(actual == expected);
    BOOST_CHECK_EQUAL(count_events(observer, stage::composition_cache), 2);
    BOOST_CHECK_EQUAL(count_events(observer, stage::powmod_x), 1);
    BOOST_CHECK(reported_path(observer, stage::frobenius_map, path::cached_composition));
    BOOST_CHECK(reported_path(observer, stage::frobenius_map, path::uncached_composition));
    BOOST_CHECK(plain_counts.operations == observed_counts.operations);
    check_balanced(observer);
}

BOOST_AUTO_TEST_CASE(equal_degree_attempts_rejections_stops_and_exceptions_preserve_sampling) {
    const polynomial_type input {-value_type(1), value_type(0), value_type(1)};
    for (int mode : {0, 1, 2}) {    // Complete, stop after first factor, throw in the first callback.
        pipeline_observer observer;
        backend_counts plain_counts, observed_counts;
        pa::polynomial_context<custom_backend> plain(custom_backend {plain_counts, {}});
        pa::polynomial_context<custom_backend, pipeline_observer> context(custom_backend {observed_counts, {}},
                                                                          observer);
        std::size_t plain_samples = 0, observed_samples = 0;
        auto plain_generator = [&] { return split_samples.at(plain_samples++); };
        auto observed_generator = [&] { return split_samples.at(observed_samples++); };
        std::vector<math::polynomial_factor<polynomial_type>> expected_factors, actual_factors;
        auto callback = [&](auto &factors, const auto &factor) {
            factors.push_back(factor);
            if (mode == 2)
                throw std::domain_error("callback failure");
            return mode == 1 ? math::factorization_control::stop_factorization :
                               math::factorization_control::continue_factorization;
        };
        auto expected_callback = [&](const auto &factor) { return callback(expected_factors, factor); };
        auto actual_callback = [&](const auto &factor) { return callback(actual_factors, factor); };
        if (mode == 2) {
            BOOST_CHECK_THROW(math::equal_degree_factorization(input, 1, plain, plain_generator, expected_callback),
                              std::domain_error);
            BOOST_CHECK_THROW(math::equal_degree_factorization(input, 1, context, observed_generator, actual_callback),
                              std::domain_error);
        } else {
            const auto expected = math::equal_degree_factorization(input, 1, plain, plain_generator, expected_callback);
            const auto actual =
                math::equal_degree_factorization(input, 1, context, observed_generator, actual_callback);
            BOOST_CHECK(actual == expected);
        }
        BOOST_CHECK(expected_factors == actual_factors);
        BOOST_CHECK(plain_counts.operations == observed_counts.operations);
        BOOST_CHECK_EQUAL(plain_samples, split_samples.size());
        BOOST_CHECK_EQUAL(observed_samples, plain_samples);
        BOOST_CHECK_EQUAL(count_events(observer, stage::random_candidate), 4);
        BOOST_CHECK_EQUAL(
            count_outcome(observer, stage::random_candidate, outcome::rejected, reason::constant_candidate), 2);
        BOOST_CHECK_EQUAL(count_events(observer, stage::splitting_attempt), 2);
        BOOST_CHECK_EQUAL(
            count_outcome(observer, stage::splitting_attempt, outcome::rejected, reason::unsuccessful_split), 1);
        BOOST_CHECK(reported_path(observer, stage::splitting_attempt, path::gcd_split));
        BOOST_CHECK_EQUAL(count_events(observer, stage::powmod), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::irreducible_factor), mode == 0 ? 2 : 1);
        const auto result = mode == 0 ? outcome::completed : mode == 1 ? outcome::callback_stopped : outcome::exception;
        BOOST_CHECK_EQUAL(count_outcome(observer, stage::equal_degree_factorization, result), 1);
        BOOST_CHECK_EQUAL(count_outcome(observer, stage::equal_degree_splitting, result), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::factor_callback, event_kind::progress), mode == 0 ? 2 :
                                                                                                mode == 1 ? 1 :
                                                                                                            0);
        check_balanced(observer);
    }
    // X splits X^2 - 1 through the power path, without an initial GCD factor.
    pipeline_observer observer;
    pa::polynomial_context<backend_type, pipeline_observer> context(backend_type {}, observer);
    const math::polynomial_divisor_context<backend_type> divisor(input, 2, context);
    const math::detail::cantor_zassenhaus_context<backend_type> split_context(1);
    std::size_t sample = 0;
    auto generator = [&] { return value_type(sample++); };
    polynomial_type factor;
    BOOST_CHECK(math::detail::try_cantor_zassenhaus_split(factor, split_context, divisor, context, generator));
    BOOST_CHECK_EQUAL(sample, 2);
    BOOST_CHECK(reported_path(observer, stage::splitting_attempt, path::power_split));
    check_balanced(observer);
}

BOOST_AUTO_TEST_CASE(generator_exceptions_leave_no_completed_trial_or_sample) {
    pipeline_observer observer;
    pa::polynomial_context<backend_type, pipeline_observer> context(backend_type {}, observer);
    std::size_t calls = 0;
    auto generator = [&]() -> value_type {
        ++calls;
        throw std::domain_error("generator failure");
    };
    BOOST_CHECK_THROW(math::equal_degree_factorization(polynomial_type {-value_type(1), 0, 1}, 1, context, generator),
                      std::domain_error);
    BOOST_CHECK_EQUAL(calls, 1);
    for (const auto operation : {stage::random_candidate, stage::splitting_attempt, stage::equal_degree_splitting,
                                 stage::equal_degree_factorization}) {
        BOOST_CHECK_EQUAL(count_outcome(observer, operation, outcome::exception), 1);
        BOOST_CHECK_EQUAL(count_events(observer, operation, event_kind::progress), 0);
    }
    check_balanced(observer);
}

BOOST_AUTO_TEST_CASE(distinct_degree_observation_distinguishes_degree_blocks_and_composition_blocks) {
    // Three generates Fq*, so both X^2 - 3 and X^3 - 3 are irreducible over this field.
    const polynomial_type quadratic {-value_type(3), 0, 1}, cubic {-value_type(3), 0, 0, 1};
    polynomial_type input;
    backend_type {}.multiply(input, quadratic, cubic);
    plain_context plain;
    const auto expected = math::distinct_degree_factorization_reference(input, plain);
    for (bool stop : {false, true}) {
        pipeline_observer observer;
        pa::polynomial_context<backend_type, pipeline_observer> context(backend_type {}, observer);
        std::vector<math::distinct_degree_factor<polynomial_type>> actual;
        const auto control =
            math::detail::kaltofen_shoup_factor_monic_square_free(actual, input, 1, context, [&](const auto &) {
                return stop ? math::factorization_control::stop_factorization :
                              math::factorization_control::continue_factorization;
            });
        BOOST_REQUIRE_EQUAL(actual.size(), stop ? 1 : expected.factors.size());
        BOOST_CHECK(std::equal(actual.begin(), actual.end(), expected.factors.begin()));
        BOOST_CHECK((control == math::factorization_control::stop_factorization) == stop);
        BOOST_CHECK_EQUAL(count_events(observer, stage::frobenius_preparation), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::kaltofen_shoup_preparation), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::frobenius_baby_steps, event_kind::progress), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::frobenius_giant_step, event_kind::progress), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::degree_block), 2);
        BOOST_CHECK_EQUAL(count_events(observer, stage::degree_block, event_kind::progress), stop ? 1 : 2);
        BOOST_CHECK_EQUAL(count_events(observer, stage::interval_product), 2);
        BOOST_CHECK_EQUAL(count_events(observer, stage::coarse_gcd), 2);
        BOOST_CHECK_EQUAL(count_outcome(observer, stage::kaltofen_shoup_factor,
                                        stop ? outcome::callback_stopped : outcome::completed),
                          1);
        bool saw_composition_size = false;
        for (std::size_t i = 0; i < observer.event_count; ++i) {
            const auto &event = observer.events[i];
            if (event.kind == event_kind::begin && event.operation == stage::degree_block)
                BOOST_CHECK_EQUAL(measurement(event, metric::degree_block_size), 1);
            if (event.kind == event_kind::end && event.operation == stage::composition_cache) {
                BOOST_CHECK_EQUAL(measurement(event, metric::composition_block_size), 3);
                saw_composition_size = true;
            }
        }
        BOOST_CHECK(saw_composition_size);
        check_balanced(observer);
    }
    // Public blocked and classical interfaces retain their callback stop and group ordering as well.
    for (bool reference : {false, true}) {
        pipeline_observer observer;
        pa::polynomial_context<backend_type, pipeline_observer> context(backend_type {}, observer);
        auto callback = [](const auto &) { return math::factorization_control::stop_factorization; };
        const auto actual = reference ? math::distinct_degree_factorization_reference(input, context, callback) :
                                        math::distinct_degree_factorization_kaltofen_shoup(input, context, callback);
        BOOST_CHECK(!actual.complete);
        BOOST_REQUIRE_EQUAL(actual.factors.size(), 1);
        BOOST_CHECK(actual.factors.front() == expected.factors.front());
        BOOST_CHECK_EQUAL(
            count_outcome(observer, reference ? stage::distinct_degree_reference : stage::distinct_degree_factorization,
                          outcome::callback_stopped),
            1);
        check_balanced(observer);
    }
}

namespace {
    template<typename Backend>
    void check_complete_observation(Backend plain_backend, Backend observed_backend) {
        backend_type builder;
        polynomial_type input {value_type(7)};
        for (const polynomial_type &factor :
             {polynomial_type {1, 1}, polynomial_type {2, 1}, polynomial_type {-value_type(3), 0, 1},
              polynomial_type {4, 1}, polynomial_type {4, 1}})
            builder.multiply(input, input, factor);
        for (int mode : {0, 1, 2}) {
            pa::polynomial_context<Backend> plain(plain_backend);
            pipeline_observer observer;
            pa::polynomial_context<Backend, pipeline_observer> context(observed_backend, observer);
            nil::crypto3::random::algebraic_engine<field_type> first_engine(17), second_engine(17);
            std::vector<value_type> first_samples, second_samples;
            auto first_generator = [&] {
                first_samples.push_back(first_engine());
                return first_samples.back();
            };
            auto second_generator = [&] {
                second_samples.push_back(second_engine());
                return second_samples.back();
            };
            std::vector<math::polynomial_factor<polynomial_type>> first_factors, second_factors;
            auto callback = [&](auto &factors, const auto &factor) {
                factors.push_back(factor);
                if (mode == 2)
                    throw std::domain_error("factor callback failure");
                return mode == 1 ? math::factorization_control::stop_factorization :
                                   math::factorization_control::continue_factorization;
            };
            auto first_callback = [&](const auto &factor) { return callback(first_factors, factor); };
            auto second_callback = [&](const auto &factor) {
                // Caller work is a child of factor_callback and remains in the factorization's inclusive total.
                auto work = context.template observe<stage::invocation>(
                    []() noexcept { return pa::polynomial_metadata {{{metric::invocation_label, 42}}}; });
                return callback(second_factors, factor);
            };
            {
                auto invocation = context.template observe<stage::invocation>();
                if (mode == 2) {
                    BOOST_CHECK_THROW(math::complete_factorization(input, plain, first_generator, first_callback),
                                      std::domain_error);
                    BOOST_CHECK_THROW(math::complete_factorization(input, context, second_generator, second_callback),
                                      std::domain_error);
                } else {
                    const auto expected = math::complete_factorization(input, plain, first_generator, first_callback);
                    const auto actual = math::complete_factorization(input, context, second_generator, second_callback);
                    BOOST_CHECK(actual == expected);
                    BOOST_CHECK_EQUAL(actual.complete, mode == 0);
                }
            }
            BOOST_CHECK(first_samples == second_samples);
            BOOST_CHECK(first_factors == second_factors);
            BOOST_CHECK_EQUAL(count_events(observer, stage::irreducible_factor), mode == 0 ? 4 : 1);
            BOOST_CHECK_EQUAL(count_outcome(observer, stage::complete_factorization,
                                            mode == 0 ? outcome::completed :
                                            mode == 1 ? outcome::callback_stopped :
                                                        outcome::exception),
                              1);
            // The first irreducible callback must follow collection of *all* distinct-degree groups of its component.
            const auto end = observer.events.begin() + observer.event_count;
            const auto group_end = std::find_if(observer.events.begin(), end, [](const auto &event) {
                return event.operation == stage::degree_group_collection && event.kind == event_kind::end;
            });
            const auto first_factor = std::find_if(observer.events.begin(), end, [](const auto &event) {
                return event.operation == stage::irreducible_factor && event.kind == event_kind::begin;
            });
            BOOST_REQUIRE(group_end != end);
            BOOST_REQUIRE(first_factor != end);
            BOOST_CHECK(group_end < first_factor);
            for (auto it = observer.events.begin(); it != end; ++it) {
                if (it->kind == event_kind::begin && it->operation == stage::invocation &&
                    measurement(*it, metric::invocation_label) == 42) {
                    const auto parent = std::find_if(observer.events.begin(), it, [&](const auto &event) {
                        return event.kind == event_kind::begin && event.scope.id == it->scope.parent_id;
                    });
                    BOOST_REQUIRE(parent != it);
                    BOOST_CHECK(parent->operation == stage::factor_callback);
                }
            }
            check_balanced(observer);
        }
    }
}    // namespace

BOOST_AUTO_TEST_CASE(complete_factorization_preserves_arithmetic_rng_callback_order_and_nesting) {
    backend_counts plain_counts, observed_counts;
    check_complete_observation(custom_backend {plain_counts, {}}, custom_backend {observed_counts, {}});
    BOOST_CHECK(plain_counts.operations == observed_counts.operations);
    check_complete_observation(mixed_backend(18), mixed_backend(18));
}

BOOST_AUTO_TEST_CASE(factorization_shortcuts_and_disabled_observers_need_no_rng_or_progress) {
    pipeline_observer observer;
    pa::polynomial_context<backend_type, pipeline_observer> context(backend_type {}, observer);
    auto generator = []() -> value_type { throw std::domain_error("unexpected RNG call"); };
    plain_context plain;
    for (const polynomial_type &input :
         {polynomial_type {value_type(0)}, polynomial_type {value_type(7)}, polynomial_type {2, 1}}) {
        BOOST_CHECK(math::complete_factorization(input, context, generator) ==
                    math::complete_factorization(input, plain, generator));
    }
    BOOST_CHECK_EQUAL(count_events(observer, stage::powmod_x), 0);
    BOOST_CHECK_EQUAL(count_events(observer, stage::random_candidate), 0);
    BOOST_CHECK_EQUAL(count_events(observer, stage::degree_block), 0);
    check_balanced(observer);
    recording_observer<false, true, 8192> kernels;
    pa::polynomial_context<backend_type, decltype(kernels)> kernels_only(backend_type {}, kernels);
    pa::polynomial_context<backend_type, disabled_observer> disabled;
    const polynomial_type input = repeated_factor_input();
    BOOST_CHECK(math::complete_factorization(input, kernels_only, generator) ==
                math::complete_factorization(input, disabled, generator));
    BOOST_REQUIRE_GT(kernels.event_count, 0);
    BOOST_REQUIRE_LE(kernels.event_count, kernels.events.size());
    for (std::size_t i = 0; i < kernels.event_count; ++i)
        BOOST_CHECK(kernels.events[i].operation >= stage::multiply);
    BOOST_CHECK(kernels.balanced);
    BOOST_CHECK_EQUAL(kernels.depth, 0);
}

namespace {
    using norm_value = fields::babybear::value_type;
    using norm_backend = pa::schoolbook_backend<norm_value>;
    using norm_polynomial = norm_backend::polynomial_type;
    using norm_counting_backend = counting_backend<norm_value>;
    using norm_observer = recording_observer<true, false, 8192, false>;

    norm_value norm_nonresidue(norm_value start = norm_value(2)) {
        while (start.is_square())
            start += norm_value::one();
        return start;
    }

    template<typename Representation>
    void check_same_representation(const std::optional<Representation> &actual,
                                   const std::optional<Representation> &expected) {
        BOOST_REQUIRE_EQUAL(actual.has_value(), expected.has_value());
        if (actual) {
            BOOST_CHECK(actual->p == expected->p);
            BOOST_CHECK(actual->q == expected->q);
        }
    }

    template<typename Observer>
    const recorded_event &begin_event(const Observer &observer, std::uint64_t id) {
        BOOST_REQUIRE_LE(observer.event_count, observer.events.size());
        const auto end = observer.events.begin() + observer.event_count;
        const auto event = std::find_if(observer.events.begin(), end, [&](const auto &item) {
            return item.kind == event_kind::begin && item.scope.id == id;
        });
        BOOST_REQUIRE(event != end);
        return *event;
    }
}    // namespace

BOOST_AUTO_TEST_CASE(norm_shortcuts_and_coefficient_filters_do_not_factor_or_sample) {
    const auto nonresidue = norm_nonresidue();
    const std::array inputs = {norm_polynomial {norm_value(0)},
                               norm_polynomial {norm_value(25)},
                               norm_polynomial {nonresidue},
                               norm_polynomial {nonresidue, 0, 1},
                               norm_polynomial {norm_value(1), norm_value(0), nonresidue},
                               norm_polynomial {norm_value(1), norm_value(0), norm_value(0), nonresidue}};
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        norm_observer observer;
        pa::polynomial_context<norm_backend> plain;
        pa::polynomial_context<norm_backend, norm_observer> context(norm_backend {}, observer);
        auto generator = []() -> norm_value { throw std::domain_error("unexpected sampling"); };
        const auto expected = math::recover_polynomial_x_norm_representation(inputs[index], plain, generator);
        const auto actual = math::recover_polynomial_x_norm_representation(inputs[index], context, generator);
        check_same_representation(actual, expected);
        BOOST_CHECK_EQUAL(actual.has_value(), index < 2);
        BOOST_CHECK_EQUAL(count_events(observer, stage::complete_factorization), 0);
        BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_final_verification), index < 2 ? 1 : 0);
        if (index >= 2) {
            const auto rejection = index == 2 ? reason::constant_not_square :
                                   index == 3 ? reason::constant_coefficient_not_square :
                                                reason::signed_leading_coefficient_not_square;
            BOOST_CHECK_EQUAL(count_outcome(observer, stage::x_norm_recovery, outcome::rejected, rejection), 1);
        }
        if (index == 3)
            BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_leading_filter), 0);
        check_balanced(observer);
    }
}

BOOST_AUTO_TEST_CASE(square_tests_report_zero_norm_shortcut_and_euler_paths) {
    pa::polynomial_context<norm_backend> plain;
    const norm_polynomial x {0, 1};
    const norm_polynomial irreducible {1, -norm_value(79), 1};
    for (norm_value leading : {norm_value(1), norm_value(9)}) {
        norm_polynomial modulus;
        math::scalar_multiplication(modulus, irreducible, leading);
        const math::polynomial_divisor_context<norm_backend> divisor(modulus, 1, plain);
        norm_observer observer;
        pa::polynomial_context<norm_backend, norm_observer> context(norm_backend {}, observer);
        BOOST_CHECK(math::is_square_mod(x, divisor, context) == math::is_square_mod(x, divisor, plain));
        BOOST_CHECK(reported_path(observer, stage::square_test, path::indeterminate_norm));
        BOOST_CHECK_EQUAL(count_events(observer, stage::powmod), 0);
        BOOST_CHECK(math::is_square_mod(norm_polynomial {norm_value(0)}, divisor, context));
        BOOST_CHECK(reported_path(observer, stage::square_test, path::zero));
        const norm_polynomial other {3, 5};
        BOOST_CHECK(math::is_square_mod(other, divisor, context) == math::is_square_mod(other, divisor, plain));
        BOOST_CHECK(reported_path(observer, stage::square_test, path::euler_criterion));
        BOOST_CHECK_EQUAL(count_events(observer, stage::powmod), 1);
        check_balanced(observer);
    }
}

BOOST_AUTO_TEST_CASE(square_root_observation_preserves_sampling_caches_and_correction_steps) {
    backend_counts expected_counts, actual_counts;
    pa::polynomial_context<norm_counting_backend> plain(norm_counting_backend {expected_counts, {}});
    norm_observer observer;
    pa::polynomial_context<norm_counting_backend, norm_observer> context(norm_counting_backend {actual_counts, {}},
                                                                         observer);
    const norm_polynomial x {0, 1}, h {-norm_nonresidue(), 0, 1};
    const math::polynomial_divisor_context<norm_counting_backend> divisor(h, 1, plain);
    actual_counts = expected_counts;
    std::size_t expected_samples = 0, actual_samples = 0;
    auto first = [&] { return ++expected_samples == 1 ? norm_polynomial {norm_value(1)} : x; };
    auto second = [&] { return ++actual_samples == 1 ? norm_polynomial {norm_value(1)} : x; };
    const math::polynomial_square_root_context<norm_counting_backend> expected_cache(divisor, plain, first);
    const math::polynomial_square_root_context<norm_counting_backend> cache(divisor, context, second);
    BOOST_CHECK_EQUAL(expected_samples, 2);
    BOOST_CHECK_EQUAL(actual_samples, expected_samples);
    BOOST_CHECK(cache.non_residue_to_odd_order() == expected_cache.non_residue_to_odd_order());
    BOOST_CHECK(expected_counts.operations == actual_counts.operations);
    BOOST_CHECK_EQUAL(
        count_outcome(observer, stage::nonresidue_attempt, outcome::rejected, reason::nonresidue_candidate_square), 1);
    BOOST_CHECK_EQUAL(count_events(observer, stage::nonresidue_search, event_kind::progress), 2);
    BOOST_CHECK_EQUAL(count_events(observer, stage::nonresidue_power), 1);
    norm_polynomial input;
    norm_backend {}.square(input, norm_polynomial {3, 5});
    math::remainder(input, input, divisor, plain);
    actual_counts = expected_counts;
    norm_polynomial expected, actual = input;
    // Reuse each cache with the other context: the cache neither owns nor stores an observer.
    BOOST_REQUIRE(math::square_root_mod(expected, input, cache, plain));
    BOOST_REQUIRE(math::square_root_mod(actual, actual, expected_cache, context));
    BOOST_CHECK(actual == expected);
    BOOST_CHECK(expected_counts.operations == actual_counts.operations);
    BOOST_CHECK_GT(count_events(observer, stage::square_root_correction), 0);
    BOOST_CHECK_EQUAL(count_events(observer, stage::square_root, event_kind::progress),
                      count_events(observer, stage::square_root_correction));
    BOOST_CHECK_EQUAL(count_events(observer, stage::square_root_preparation), 1);
    BOOST_CHECK(!math::square_root_mod(actual, x, cache, context));
    BOOST_CHECK(actual == norm_polynomial({norm_value(0)}));
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::square_root, outcome::rejected, reason::not_square), 1);
    BOOST_CHECK(math::square_root_mod(actual, actual, cache, context));
    BOOST_CHECK(reported_path(observer, stage::square_root, path::zero));
    for (std::size_t i = 0; i < observer.event_count; ++i) {
        const auto &event = observer.events[i];
        if (event.kind == event_kind::progress && event.operation == stage::square_root) {
            BOOST_REQUIRE_GT(i, 0);
            BOOST_CHECK(observer.events[i - 1].kind == event_kind::end);
            BOOST_CHECK(observer.events[i - 1].operation == stage::square_root_correction);
        }
    }
    check_balanced(observer);
}

BOOST_AUTO_TEST_CASE(square_root_preparation_exceptions_close_attempts_without_false_progress) {
    pa::polynomial_context<norm_backend> plain;
    const math::polynomial_divisor_context<norm_backend> divisor(norm_polynomial {-norm_nonresidue(), 0, 1}, 1, plain);
    for (bool throws : {false, true}) {
        norm_observer observer;
        pa::polynomial_context<norm_backend, norm_observer> context(norm_backend {}, observer);
        std::size_t calls = 0;
        auto generator = [&]() -> norm_polynomial {
            if (++calls == 1)
                return norm_polynomial {norm_value(1)};
            if (throws)
                throw std::domain_error("nonresidue generator failure");
            return norm_polynomial {1, 0};    // The existing validator rejects this noncanonical candidate.
        };
        if (throws) {
            BOOST_CHECK_THROW((math::polynomial_square_root_context<norm_backend>(divisor, context, generator)),
                              std::domain_error);
        } else {
            BOOST_CHECK_THROW((math::polynomial_square_root_context<norm_backend>(divisor, context, generator)),
                              std::invalid_argument);
        }
        BOOST_CHECK_EQUAL(calls, 2);
        BOOST_CHECK_EQUAL(count_events(observer, stage::nonresidue_search, event_kind::progress), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::nonresidue_attempt, event_kind::progress), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::nonresidue_power), 0);
        for (const auto operation :
             {stage::square_root_preparation, stage::nonresidue_search, stage::nonresidue_attempt})
            BOOST_CHECK_EQUAL(count_outcome(observer, operation, outcome::exception), 1);
        check_balanced(observer);
    }
}

BOOST_AUTO_TEST_CASE(rational_reconstruction_observation_preserves_bounds_aliasing_and_backend_calls) {
    const norm_polynomial modulus {3, 3, 3, 2, 1}, residue {2, 1, 1, 1};
    for (std::size_t bound : {1, 2}) {
        norm_observer observer;
        backend_counts expected_counts, actual_counts;
        pa::polynomial_context<norm_counting_backend> plain(norm_counting_backend {expected_counts, {}});
        pa::polynomial_context<norm_counting_backend, norm_observer> context(norm_counting_backend {actual_counts, {}},
                                                                             observer);
        norm_polynomial p = residue, q = modulus, expected_p = p, expected_q = q;
        const bool expected =
            math::rational_reconstruct(expected_p, expected_q, expected_p, expected_q, 0, bound, plain);
        const bool actual = math::rational_reconstruct(p, q, p, q, 0, bound, context);
        BOOST_CHECK_EQUAL(actual, expected);
        BOOST_CHECK_EQUAL(actual, bound == 2);
        BOOST_CHECK(p == expected_p && q == expected_q);
        BOOST_CHECK(expected_counts.operations == actual_counts.operations);
        if (!actual) {
            BOOST_CHECK(p == residue && q == modulus);
            BOOST_CHECK_EQUAL(count_outcome(observer, stage::rational_reconstruction, outcome::rejected,
                                            reason::rational_degree_bound),
                              1);
        }
        BOOST_CHECK_GT(count_events(observer, stage::rational_reconstruction_step), 0);
        BOOST_CHECK_EQUAL(count_events(observer, stage::rational_reconstruction_step),
                          count_events(observer, stage::rational_reconstruction, event_kind::progress));
        BOOST_CHECK_EQUAL(count_events(observer, stage::rational_reconstruction_normalization), actual ? 1 : 0);
        BOOST_CHECK_THROW(math::rational_reconstruct(p, p, residue, modulus, 0, bound, context), std::invalid_argument);
        BOOST_CHECK_EQUAL(count_outcome(observer, stage::rational_reconstruction, outcome::exception), 1);
        check_balanced(observer);
    }
}

namespace {
    template<typename Backend>
    void check_norm_pipeline(Backend plain_backend, Backend observed_backend,
                             const typename Backend::polynomial_type &input, bool succeeds, std::size_t factor_count) {
        using value = typename Backend::polynomial_type::value_type;
        using field = typename value::field_type;
        norm_observer observer;
        pa::polynomial_context<Backend> plain(plain_backend);
        pa::polynomial_context<Backend, norm_observer> context(observed_backend, observer);
        nil::crypto3::random::algebraic_engine<field> first_engine(173), second_engine(173);
        std::vector<value> first_samples, second_samples;
        auto first_generator = [&] {
            first_samples.push_back(first_engine());
            return first_samples.back();
        };
        auto second_generator = [&] {
            second_samples.push_back(second_engine());
            return second_samples.back();
        };
        const auto expected = math::recover_polynomial_x_norm_representation(input, plain, first_generator);
        const auto actual = math::recover_polynomial_x_norm_representation(input, context, second_generator);
        check_same_representation(actual, expected);
        BOOST_CHECK_EQUAL(actual.has_value(), succeeds);
        BOOST_CHECK(first_samples == second_samples);
        BOOST_CHECK_EQUAL(count_events(observer, stage::complete_factorization), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_factor), factor_count);
        BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_factor, event_kind::progress), factor_count);
        BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_final_verification), succeeds ? 1 : 0);
        if (succeeds) {
            // Use an independent plain context so checking the identity does not alter either compared call trace.
            pa::polynomial_context<pa::schoolbook_backend<value>> verification;
            BOOST_CHECK(math::evaluate_polynomial_x_norm(*actual, verification) == input);
            BOOST_CHECK(reported_path(observer, stage::x_norm_multiplicity, path::even_multiplicity));
            BOOST_CHECK(reported_path(observer, stage::x_norm_multiplicity, path::odd_multiplicity));
            BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_combination), 1);
            BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_multiply), factor_count - 1);
            BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_combination_level), factor_count == 3 ? 2 : 1);
            BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_combination, event_kind::progress),
                              factor_count == 3 ? 2 : 1);
            BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_leading_scalar), 1);
            BOOST_CHECK_EQUAL(count_outcome(observer, stage::x_norm_recovery, outcome::completed), 1);
            BOOST_CHECK_EQUAL(count_outcome(observer, stage::complete_factorization, outcome::completed), 1);
            BOOST_CHECK_GT(count_events(observer, stage::rational_reconstruction), 0);
            BOOST_CHECK_GT(count_events(observer, stage::x_norm_normalization), 0);
        } else {
            BOOST_CHECK_EQUAL(
                count_outcome(observer, stage::x_norm_recovery, outcome::rejected, reason::odd_factor_recovery), 1);
            BOOST_CHECK_EQUAL(
                count_outcome(observer, stage::x_norm_factor, outcome::rejected, reason::odd_factor_recovery), 1);
            BOOST_CHECK_EQUAL(
                count_outcome(observer, stage::x_norm_irreducible, outcome::rejected, reason::x_not_square), 1);
            BOOST_CHECK_EQUAL(count_outcome(observer, stage::complete_factorization, outcome::callback_stopped), 1);
            BOOST_CHECK_EQUAL(count_events(observer, stage::square_root_preparation), 0);
            BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_combination), 0);
        }
        // Norm work belongs inside factor_callback, and ultimately beneath the single top-level recovery call.
        // A factor's progress is emitted after its recovery has returned (and after insertion when successful).
        std::vector<std::pair<std::size_t, std::size_t>> factors;
        for (std::size_t i = 0; i < observer.event_count; ++i) {
            const auto &event = observer.events[i];
            if (event.kind == event_kind::begin && event.operation == stage::x_norm_factor) {
                const auto &parent = begin_event(observer, event.scope.parent_id);
                BOOST_CHECK(parent.operation == stage::factor_callback);
                factors.emplace_back(measurement(event, metric::factor_degree),
                                     measurement(event, metric::multiplicity));
                const auto &multiplicity = begin_event(observer, observer.events[i + 1].scope.id);
                BOOST_CHECK(multiplicity.operation == stage::x_norm_multiplicity);
            }
            if (event.kind == event_kind::progress && event.operation == stage::x_norm_factor) {
                BOOST_REQUIRE_GT(i, 0);
                BOOST_CHECK(observer.events[i - 1].kind == event_kind::end);
                BOOST_CHECK(observer.events[i - 1].operation == stage::x_norm_multiplicity);
            }
        }
        BOOST_CHECK_EQUAL(factors.size(), factor_count);
        if (succeeds) {
            // Reuse the same seeded generator for a direct factorization; callbacks do not change factor order.
            // The original result/consumption comparison above also covers sharing RNG with norm callbacks.
            nil::crypto3::random::algebraic_engine<field> generator(173);
            pa::polynomial_context<pa::schoolbook_backend<value>> order_context;
            const auto factored = math::complete_factorization(input, order_context, generator);
            BOOST_REQUIRE_EQUAL(factored.factors.size(), factors.size());
            for (std::size_t i = 0; i < factors.size(); ++i) {
                BOOST_CHECK_EQUAL(factors[i].first, factored.factors[i].polynomial.size() - 1);
                BOOST_CHECK_EQUAL(factors[i].second, factored.factors[i].multiplicity);
            }
        }
        check_balanced(observer);
    }
}    // namespace

BOOST_AUTO_TEST_CASE(norm_pipeline_observation_preserves_even_odd_recovery_combination_and_rng) {
    norm_polynomial input {norm_value(25)};
    const norm_polynomial quadratic {1, -norm_value(79), 1};
    const norm_polynomial even_factor {-norm_nonresidue(), 1};
    norm_backend builder;
    // 25 * (X-4) * (X-c)^2 * (X^2-79X+1)^3, c nonsquare: three representations, including an odd carry in the tree.
    for (const auto &factor :
         {norm_polynomial {-norm_value(4), 1}, even_factor, even_factor, quadratic, quadratic, quadratic})
        builder.multiply(input, input, factor);
    backend_counts first_counts, second_counts;
    check_norm_pipeline(norm_counting_backend {first_counts, {}}, norm_counting_backend {second_counts, {}}, input,
                        true, 3);
    BOOST_CHECK(first_counts.operations == second_counts.operations);
}

BOOST_AUTO_TEST_CASE(norm_pipeline_rejection_preserves_first_failing_factor_and_stop_point) {
    const auto first = norm_nonresidue(), second = norm_nonresidue(first + norm_value::one());
    norm_polynomial input;
    norm_backend {}.multiply(input, norm_polynomial {-first, 1}, norm_polynomial {-second, 1});
    BOOST_REQUIRE(input.front().is_square() && input.back().is_square());
    backend_counts first_counts, second_counts;
    check_norm_pipeline(norm_counting_backend {first_counts, {}}, norm_counting_backend {second_counts, {}}, input,
                        false, 1);
    BOOST_CHECK(first_counts.operations == second_counts.operations);
}

BOOST_AUTO_TEST_CASE(norm_pipeline_supports_observed_schoolbook_and_mixed_fp12_backends) {
    using field = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using value = field::value_type;
    using schoolbook = pa::schoolbook_backend<value>;
    using mixed = pa::mixed_radix_backend<field_type, value>;
    using polynomial = schoolbook::polynomial_type;
    auto root = value::zero();
    for (std::size_t i = 0; i < field::arity; ++i)
        root.coordinate(i) = value_type(i + 1);
    const polynomial odd {-root.squared(), value::one()};
    const polynomial even {-(root + value::one()).squared(), value::one()};
    polynomial input {value(25)};
    schoolbook backend;
    for (const auto &factor : {even, even, odd, odd, odd})
        backend.multiply(input, input, factor);
    check_norm_pipeline(schoolbook {}, schoolbook {}, input, true, 2);
    check_norm_pipeline(mixed(18), mixed(18), input, true, 2);
}

BOOST_AUTO_TEST_CASE(irreducible_norm_scalar_rejection_is_distinct_from_square_test_failure) {
    const auto nonresidue = norm_nonresidue();
    norm_observer observer;
    pa::polynomial_context<norm_backend, norm_observer> context(norm_backend {}, observer);
    pa::polynomial_context<norm_backend> plain;
    auto generator = [&] { return nonresidue; };
    const norm_polynomial input {-nonresidue, nonresidue};
    const auto expected = math::recover_irreducible_polynomial_x_norm_representation(input, plain, generator);
    const auto actual = math::recover_irreducible_polynomial_x_norm_representation(input, context, generator);
    check_same_representation(actual, expected);
    BOOST_CHECK(!actual);
    BOOST_CHECK_EQUAL(
        count_outcome(observer, stage::x_norm_irreducible, outcome::rejected, reason::nonsquare_normalization_scalar),
        1);
    BOOST_CHECK_EQUAL(
        count_outcome(observer, stage::x_norm_scalar_check, outcome::rejected, reason::nonsquare_normalization_scalar),
        1);
    BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_normalization), 0);
    BOOST_CHECK_EQUAL(count_events(observer, stage::square_root), 1);
    check_balanced(observer);
}

BOOST_AUTO_TEST_CASE(norm_exceptions_preserve_scope_balance_and_do_not_report_rejection) {
    norm_observer observer;
    pa::polynomial_context<norm_backend, norm_observer> context(norm_backend {}, observer);
    auto generator = []() -> norm_value { throw std::domain_error("recovery generator failure"); };
    // One irreducible quadratic reaches the norm callback without randomized equal-degree splitting.
    BOOST_CHECK_THROW(
        math::recover_polynomial_x_norm_representation(norm_polynomial {1, -norm_value(79), 1}, context, generator),
        std::domain_error);
    for (const auto operation : {stage::x_norm_recovery, stage::complete_factorization, stage::factor_callback,
                                 stage::x_norm_factor, stage::x_norm_irreducible, stage::square_root_preparation})
        BOOST_CHECK_EQUAL(count_outcome(observer, operation, outcome::exception), 1);
    BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_factor, event_kind::progress), 0);
    BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_final_verification), 0);
    check_balanced(observer);

    backend_counts counts;
    norm_observer failing_observer;
    pa::polynomial_context<norm_counting_backend, norm_observer> failing(norm_counting_backend {counts, {}},
                                                                         failing_observer);
    counts.throw_on_square = true;
    BOOST_CHECK_THROW(
        math::recover_polynomial_x_norm_representation(norm_polynomial {norm_value(0)}, failing, generator),
        std::domain_error);
    BOOST_CHECK_EQUAL(count_outcome(failing_observer, stage::x_norm_final_verification, outcome::exception), 1);
    BOOST_CHECK_EQUAL(count_outcome(failing_observer, stage::x_norm_recovery, outcome::exception), 1);
    check_balanced(failing_observer);
}

BOOST_AUTO_TEST_CASE(norm_disabled_and_kernel_only_contexts_keep_public_cache_and_call_syntax) {
    pa::polynomial_context<norm_backend, disabled_observer> disabled;
    recording_observer<false, true, 8192> observer;
    pa::polynomial_context<norm_backend, decltype(observer)> context(norm_backend {}, observer);
    const norm_polynomial input {1, -norm_value(79), 1};
    nil::crypto3::random::algebraic_engine<fields::babybear> first(19), second(19);
    check_same_representation(math::recover_polynomial_x_norm_representation(input, context, second),
                              math::recover_polynomial_x_norm_representation(input, disabled, first));
    BOOST_REQUIRE_GT(observer.event_count, 0);
    BOOST_REQUIRE_LE(observer.event_count, observer.events.size());
    for (std::size_t i = 0; i < observer.event_count; ++i)
        BOOST_CHECK(observer.events[i].operation >= stage::multiply);
    BOOST_CHECK(observer.balanced);
    BOOST_CHECK_EQUAL(observer.depth, 0);
}

namespace {
    using kernel_observer = recording_observer<false, true, 2048>;
    using extension_field = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using extension_backend = pa::mixed_radix_backend<field_type, extension_field::value_type>;

    // Populate all extension coordinates; FFT coverage must not accidentally stay in the base subfield.
    template<typename Value>
    Value dense_coefficient(std::size_t seed) {
        if constexpr (std::is_same_v<Value, value_type>) {
            return Value(seed);
        } else {
            auto result = Value::zero();
            for (std::size_t i = 0; i < extension_field::arity; ++i)
                result.coordinate(i) = value_type(seed + i);
            return result;
        }
    }

    template<typename Observer>
    std::size_t count_measurement(const Observer &observer, stage operation, metric key, std::size_t value,
                                  event_kind kind = event_kind::begin) {
        BOOST_REQUIRE_LE(observer.event_count, observer.events.size());
        return std::count_if(
            observer.events.begin(), observer.events.begin() + observer.event_count, [&](const auto &event) {
                return event.kind == kind && event.operation == operation && measurement(event, key) == value;
            });
    }

    template<typename Backend>
    void check_kernel_operations() {
        using polynomial = typename Backend::polynomial_type;
        using value = typename polynomial::value_type;
        kernel_observer observer;
        pa::polynomial_context context(Backend(18), observer);
        pa::schoolbook_backend<value> reference;
        const polynomial a {dense_coefficient<value>(1), dense_coefficient<value>(3), dense_coefficient<value>(5),
                            dense_coefficient<value>(7)};
        const polynomial b {dense_coefficient<value>(2), dense_coefficient<value>(4), dense_coefficient<value>(6),
                            dense_coefficient<value>(8)};
        polynomial actual = a, expected;
        reference.multiply(expected, a, b);
        context.multiply(actual, actual, b);
        BOOST_CHECK(actual == expected);
        reference.square(expected, a);
        actual = a;
        context.square(actual, actual);
        BOOST_CHECK(actual == expected);
        reference.multiply_low(expected, a, b, 3);
        actual = b;
        context.multiply_low(actual, a, actual, 3);
        BOOST_CHECK(actual == expected);
        BOOST_CHECK_EQUAL(count_events(observer, stage::multiply), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::square), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::multiply_low), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::fft_forward), 5);
        BOOST_CHECK_EQUAL(count_events(observer, stage::fft_inverse), 3);
        BOOST_CHECK_EQUAL(count_measurement(observer, stage::fft_forward, metric::transform_length, 9), 3);
        // A low product still transforms the full product of its prefixes, not a length-three cyclic product.
        BOOST_CHECK_EQUAL(count_measurement(observer, stage::fft_forward, metric::transform_length, 6), 2);
        BOOST_CHECK_EQUAL(count_events(observer, stage::pointwise_multiply), 2);
        BOOST_CHECK_EQUAL(count_events(observer, stage::pointwise_square), 1);
        BOOST_CHECK_EQUAL(count_events(observer, stage::output_finalization), 3);
        for (const auto operation : {stage::multiply, stage::square, stage::multiply_low})
            BOOST_CHECK_EQUAL(count_measurement(observer, operation, metric::detailed_backend, 1, event_kind::end), 1);
        for (std::size_t i = 0; i < observer.event_count; ++i) {
            const auto &event = observer.events[i];
            if (event.kind != event_kind::begin)
                continue;
            if (event.operation == stage::fft_normalization)
                BOOST_CHECK(begin_event(observer, event.scope.parent_id).operation == stage::fft_inverse);
            if (event.operation == stage::fft_forward || event.operation == stage::fft_inverse) {
                BOOST_REQUIRE_LT(i + 1, observer.event_count);
                BOOST_CHECK(observer.events[i + 1].operation == stage::buffer_preparation);
                BOOST_CHECK_EQUAL(observer.events[i + 1].scope.parent_id, event.scope.id);
            }
        }
        context.square(actual, polynomial {value::zero()});
        context.multiply(actual, a, polynomial {value::zero()});
        context.multiply_low(actual, a, b, 0);
        BOOST_CHECK(math::is_zero(actual));
        BOOST_CHECK_EQUAL(count_events(observer, stage::fft_forward), 5);
        BOOST_CHECK_EQUAL(count_events(observer, stage::fft_inverse), 3);
        check_balanced(observer);
    }

    template<typename Backend>
    void check_cached_squaremod_kernels(std::size_t n) {
        using polynomial = typename Backend::polynomial_type;
        using value = typename polynomial::value_type;
        pa::polynomial_context_options options;
        options.basecase_divisor_coefficient_cutoff = 0;
        options.basecase_quotient_coefficient_cutoff = 0;
        kernel_observer preparation, reuse;
        pa::polynomial_context setup(Backend(18), preparation, options);
        pa::polynomial_context warm(Backend(18), reuse, options);
        pa::polynomial_context<pa::schoolbook_backend<value>> reference;
        polynomial h(n + 1, value::zero()), input(n, value::zero());
        for (std::size_t i = 0; i <= n; ++i)
            h[i] = dense_coefficient<value>(i + 1);
        for (std::size_t i = 0; i < n; ++i)
            input[i] = dense_coefficient<value>(2 * i + 1);
        const math::polynomial_divisor_context<Backend> divisor(h, n - 1, setup);
        const math::polynomial_divisor_context<pa::schoolbook_backend<value>> reference_divisor(h, n - 1, reference);
        BOOST_CHECK_EQUAL(count_events(preparation, stage::prepare_low_product), 1);
        BOOST_CHECK_EQUAL(count_events(preparation, stage::prepare_cyclic_remainder), 1);
        BOOST_CHECK_EQUAL(
            count_measurement(preparation, stage::prepare_low_product, metric::transform_length, 18, event_kind::end),
            1);
        BOOST_CHECK_EQUAL(count_measurement(preparation, stage::prepare_cyclic_remainder, metric::transform_length, 9,
                                            event_kind::end),
                          1);
        BOOST_CHECK_EQUAL(count_measurement(preparation, stage::prepare_cyclic_remainder,
                                            metric::prepared_storage_bytes, 9 * sizeof(value), event_kind::end),
                          1);
        // Reuse across contexts and changing residues must leave both fixed spectra intact.
        for (std::size_t trial = 0; trial < 2; ++trial) {
            polynomial actual = input, expected;
            math::squaremod(expected, input, reference_divisor, reference);
            math::squaremod(actual, actual, divisor, warm);
            BOOST_CHECK(actual == expected);
            input[0] += value::one();
        }
        BOOST_CHECK_EQUAL(count_events(reuse, stage::prepare_low_product), 0);
        BOOST_CHECK_EQUAL(count_events(reuse, stage::prepare_cyclic_remainder), 0);
        BOOST_CHECK_EQUAL(count_outcome(reuse, stage::try_multiply_low_prepared, outcome::completed), 2);
        BOOST_CHECK_EQUAL(count_outcome(reuse, stage::try_cyclic_remainder, outcome::completed), 2);
        BOOST_CHECK_EQUAL(count_events(reuse, stage::multiply_low), 0);
        BOOST_CHECK_EQUAL(count_events(reuse, stage::low_product_fallback), 0);
        // Per square: two full transforms for the square, two for Q, two short transforms for R.
        for (const auto operation : {stage::fft_forward, stage::fft_inverse}) {
            BOOST_CHECK_EQUAL(count_events(reuse, operation), 6);
            BOOST_CHECK_EQUAL(count_measurement(reuse, operation, metric::transform_length, 18), 4);
            BOOST_CHECK_EQUAL(count_measurement(reuse, operation, metric::transform_length, 9), 2);
        }
        BOOST_CHECK_EQUAL(count_events(reuse, stage::pointwise_square), 2);
        BOOST_CHECK_EQUAL(count_events(reuse, stage::pointwise_multiply), 4);
        BOOST_CHECK_EQUAL(count_events(reuse, stage::cyclic_folding), 4);
        BOOST_CHECK_EQUAL(count_events(reuse, stage::coefficient_subtraction), 2);
        check_balanced(preparation);
        check_balanced(reuse);
    }
}    // namespace

BOOST_AUTO_TEST_CASE(mixed_radix_kernels_preserve_products_squares_low_products_and_aliasing) {
    check_kernel_operations<mixed_backend>();
    check_kernel_operations<extension_backend>();
}

BOOST_AUTO_TEST_CASE(cached_squaremod_reports_actual_full_and_short_transforms) {
    check_cached_squaremod_kernels<mixed_backend>(8);
    check_cached_squaremod_kernels<extension_backend>(9);    // M == degree(H) folds the leading coefficient too.
}

BOOST_AUTO_TEST_CASE(changing_quotient_lengths_reports_decline_and_existing_fallbacks) {
    pa::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    pa::polynomial_context plain(mixed_backend(18), options);
    const polynomial_type h {1, 2, 3, 4, 5, 6, 7, 8, 2}, remainder {3, 5};
    const math::polynomial_divisor_context<mixed_backend> divisor(h, 7, plain);
    backend_type reference;
    for (const polynomial_type &q : {polynomial_type {2, 3, 4}, polynomial_type {2}}) {
        kernel_observer observer;
        pa::polynomial_context context(mixed_backend(18), observer, options);
        polynomial_type dividend, actual_q, actual_r;
        reference.multiply(dividend, h, q);
        math::addition(dividend, dividend, remainder);
        math::divrem(actual_q, actual_r, dividend, divisor, context);
        BOOST_CHECK(actual_q == q);
        BOOST_CHECK(actual_r == remainder);
        BOOST_CHECK_EQUAL(count_outcome(observer, stage::try_multiply_low_prepared, outcome::declined), 1);
        BOOST_CHECK_EQUAL(count_measurement(observer, stage::low_product_fallback, metric::prepared_available, 1), 1);
        if (q.size() == 1) {
            BOOST_CHECK_EQUAL(count_outcome(observer, stage::try_cyclic_remainder, outcome::declined), 1);
            BOOST_CHECK_EQUAL(count_measurement(observer, stage::low_product_fallback, metric::prepared_available, 0),
                              1);
        } else {
            BOOST_CHECK_EQUAL(count_outcome(observer, stage::try_cyclic_remainder, outcome::completed), 1);
            BOOST_CHECK_EQUAL(count_events(observer, stage::low_product_fallback), 1);
        }
        const auto count = observer.event_count;
        math::divrem(actual_q, actual_r, remainder, divisor, context);
        BOOST_CHECK(math::is_zero(actual_q));
        BOOST_CHECK(actual_r == remainder);
        BOOST_CHECK_EQUAL(observer.event_count, count);    // No backend work below the modulus degree.
        check_balanced(observer);
    }
}

BOOST_AUTO_TEST_CASE(prepared_attempts_distinguish_zero_success_decline_and_owned_reuse) {
    kernel_observer observer;
    pa::polynomial_context context(mixed_backend(18), observer);
    const polynomial_type input {1, 2, 3}, fixed {4, 5, 6};
    polynomial_type output {7}, expected;
    const auto zero = context.prepare_low_product(polynomial_type {0}, 3, 3);
    BOOST_REQUIRE(zero);
    BOOST_REQUIRE(context.try_multiply_low_prepared(output, input, *zero, 3));
    BOOST_CHECK(math::is_zero(output));
    BOOST_CHECK(!context.try_multiply_low_prepared(output, input, *zero, 2));
    BOOST_CHECK(math::is_zero(output));
    BOOST_CHECK_EQUAL(count_events(observer, stage::fft_forward), 0);
    BOOST_CHECK_EQUAL(
        count_measurement(observer, stage::prepare_low_product, metric::prepared_storage_bytes, 0, event_kind::end), 1);
    auto prepared = context.prepare_low_product(fixed, 3, 3);
    BOOST_REQUIRE(prepared);
    auto owned = std::move(*prepared);
    BOOST_CHECK(!context.try_multiply_low_prepared(output, input, *prepared, 3));
    BOOST_CHECK(math::is_zero(output));
    backend_type reference;
    reference.multiply_low(expected, input, fixed, 3);
    for (std::size_t i = 0; i < 2; ++i) {
        BOOST_REQUIRE(context.try_multiply_low_prepared(output, input, owned, 3));
        BOOST_CHECK(output == expected);
    }
    BOOST_CHECK(!context.prepare_cyclic_remainder(polynomial_type {1, 2}, 1));
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::prepare_cyclic_remainder, outcome::declined), 1);
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::try_multiply_low_prepared, outcome::declined), 2);
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::try_multiply_low_prepared, outcome::completed), 3);
    BOOST_CHECK_EQUAL(count_events(observer, stage::fft_forward), 3);
    BOOST_CHECK_EQUAL(count_events(observer, stage::fft_inverse), 2);
    check_balanced(observer);
}

BOOST_AUTO_TEST_CASE(kernel_exceptions_balance_scopes_without_changing_outputs) {
    kernel_observer observer;
    pa::polynomial_context context(mixed_backend(3), observer);
    polynomial_type output {7}, input {1, 2, 3};
    BOOST_CHECK_THROW(context.square(output, input), std::invalid_argument);
    BOOST_CHECK(output == polynomial_type({value_type(7)}));
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::square, outcome::exception), 1);
    BOOST_CHECK_EQUAL(count_events(observer, stage::fft_forward), 0);
    math::mixed_radix_fft_plan<field_type> plan(3);
    std::vector<value_type> oversized(4, value_type(1)), workspace;
    BOOST_CHECK_THROW(plan.fft(oversized, workspace, observer), std::invalid_argument);
    BOOST_CHECK_EQUAL(oversized.size(), 4);
    BOOST_CHECK(workspace.empty());
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::fft_forward, outcome::exception), 1);
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::buffer_preparation, outcome::exception), 1);
    BOOST_CHECK_EQUAL(count_events(observer, stage::fft_inverse), 0);
    check_balanced(observer);
}

BOOST_AUTO_TEST_CASE(norm_recovery_attributes_nested_kernels_to_one_caller_scope) {
    recording_observer<true, true, 8192> observer;
    pa::polynomial_context context(mixed_backend(18), observer);
    pa::polynomial_context plain(mixed_backend(18));
    auto generator = []() -> value_type { throw std::domain_error("unexpected RNG call"); };
    const polynomial_type input {4, 4, 1};    // (X + 2)^2 takes the even-multiplicity path.
    const auto expected = math::recover_polynomial_x_norm_representation(input, plain, generator);
    const auto actual = math::recover_polynomial_x_norm_representation(input, context, generator);
    check_same_representation(actual, expected);
    BOOST_REQUIRE(actual);
    BOOST_CHECK_EQUAL(count_events(observer, stage::x_norm_recovery), 1);
    BOOST_CHECK_EQUAL(count_events(observer, stage::complete_factorization), 1);
    BOOST_CHECK_GT(count_events(observer, stage::fft_forward), 0);
    for (std::size_t i = 0; i < observer.event_count; ++i) {
        const auto &event = observer.events[i];
        if (event.kind != event_kind::begin)
            continue;
        auto id = event.scope.id;
        while (begin_event(observer, id).scope.parent_id != 0)
            id = begin_event(observer, id).scope.parent_id;
        BOOST_CHECK(begin_event(observer, id).operation == stage::x_norm_recovery);
    }
    check_balanced(observer);
}

namespace {
    // A legacy optional capability with no observer overloads or FFT/storage introspection.
    struct custom_prepared_backend : backend_type {
        struct prepared_low_product_type {
            polynomial_type fixed;
            std::size_t precision;
        };
        std::optional<prepared_low_product_type> prepare_low_product(const polynomial_type &fixed, std::size_t,
                                                                     std::size_t precision) {
            return prepared_low_product_type {fixed, precision};
        }
        bool try_multiply_low_prepared(polynomial_type &output, const polynomial_type &left,
                                       const prepared_low_product_type &prepared, std::size_t precision) {
            if (precision != prepared.precision)
                return false;
            multiply_low(output, left, prepared.fixed, precision);
            return true;
        }
    };
}    // namespace

BOOST_AUTO_TEST_CASE(custom_prepared_backend_keeps_optional_capability_without_observer_methods) {
    kernel_observer observer;
    pa::polynomial_context context(custom_prepared_backend {}, observer);
    const polynomial_type fixed {1, 2, 3}, left {2, 3, 4};
    const auto prepared = context.prepare_low_product(fixed, 3, 3);
    BOOST_REQUIRE(prepared);
    polynomial_type expected, actual;
    backend_type reference;
    reference.multiply_low(expected, left, fixed, 3);
    BOOST_REQUIRE(context.try_multiply_low_prepared(actual, left, *prepared, 3));
    BOOST_CHECK(actual == expected);
    BOOST_CHECK(!context.try_multiply_low_prepared(actual, left, *prepared, 2));
    BOOST_CHECK(actual == expected);
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::prepare_low_product, outcome::completed), 1);
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::try_multiply_low_prepared, outcome::completed), 1);
    BOOST_CHECK_EQUAL(count_outcome(observer, stage::try_multiply_low_prepared, outcome::declined), 1);
    BOOST_CHECK_EQUAL(observer.event_count, 6);
    for (std::size_t i = 0; i < observer.event_count; ++i) {
        BOOST_CHECK_EQUAL(measurement(observer.events[i], metric::transform_length),
                          std::numeric_limits<std::size_t>::max());
        BOOST_CHECK_EQUAL(measurement(observer.events[i], metric::prepared_storage_bytes),
                          std::numeric_limits<std::size_t>::max());
    }
    check_balanced(observer);
}
