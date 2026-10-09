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

#define BOOST_TEST_MODULE polynomial_observation_collector_test
#include <numeric>
#include <sstream>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/arithmetic_params/babybear.hpp>
#include <nil/crypto3/math/polynomial/polynomial_observation_collector.hpp>
#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/reconstruction/polynomial_x_norm_reconstruction.hpp>
#include <nil/crypto3/random/algebraic_engine.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace pa = math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;
    using stage = pa::polynomial_stage;
    using metric = pa::polynomial_metric;
    using outcome = pa::polynomial_scope_outcome;
    using reason = pa::polynomial_rejection_reason;
    using ns = std::chrono::nanoseconds;

    // Caller-owned clock state: timing tests use no mutable global state and no timing tolerances.
    struct manual_clock {
        using duration = ns;
        using time_point = std::chrono::time_point<manual_clock>;
        static constexpr bool is_steady = true;
        ns *time;
        time_point now() const noexcept {
            return time_point(*time);
        }
    };
    using collector_type = pa::polynomial_observation_collector<true, true, manual_clock>;
    static_assert(pa::PolynomialObserver<collector_type>);
    static_assert(!std::is_move_constructible_v<collector_type>);
    static_assert(!std::is_copy_constructible_v<collector_type>);

    auto empty_metadata = []() noexcept { return pa::polynomial_metadata {}; };

    template<typename Collector>
    std::size_t calls(const Collector &collector, stage operation) {
        std::size_t total = 0;
        for (std::size_t i = 0; i < collector.size(); ++i)
            if (collector[i].stage == operation)
                total += collector[i].calls;
        return total;
    }

    template<typename Collector>
    void check_totals(const Collector &collector) {
        ns roots {}, exclusive {};
        for (std::size_t i = 0; i < collector.size(); ++i) {
            const auto &t = collector[i];
            BOOST_CHECK_LE(t.parent_id, i);
            BOOST_CHECK_GE(t.inclusive.count(), t.exclusive.count());
            BOOST_CHECK_GE(t.exclusive.count(), 0);
            if (t.parent_id == 0)
                roots += t.inclusive;
            exclusive += t.exclusive;
            BOOST_CHECK_EQUAL(std::accumulate(t.outcomes.begin(), t.outcomes.end(), std::size_t(0)), t.calls);
        }
        BOOST_CHECK_EQUAL(roots.count(), exclusive.count());
        BOOST_CHECK(collector.balanced());
        BOOST_CHECK_EQUAL(collector.dropped_scopes(), 0);
        BOOST_CHECK_EQUAL(collector.dropped_metadata(), 0);
    }
}    // namespace

BOOST_AUTO_TEST_CASE(hierarchical_times_repeated_calls_and_completed_progress_are_exact) {
    ns time {};
    collector_type collector(8, 8, manual_clock {&time});
    {
        auto root = pa::observe_polynomial<stage::invocation>(collector, empty_metadata);
        time = ns(10);
        {
            auto child = pa::observe_polynomial<stage::powmod_x>(collector, empty_metadata);
            time = ns(20);
            child.progress(1, 3);
            child.progress(1, 3);    // A repeated heartbeat must not count the same work twice.
            time = ns(40);
            child.progress(3, 3);
        }
        time = ns(50);
        {
            auto child = pa::observe_polynomial<stage::powmod_x>(collector, empty_metadata);
            child.progress(2, 2);
            time = ns(60);
        }
        time = ns(100);
    }
    BOOST_REQUIRE_EQUAL(collector.size(), 2);
    BOOST_CHECK_EQUAL(collector[0].inclusive.count(), 100);
    BOOST_CHECK_EQUAL(collector[0].exclusive.count(), 60);
    BOOST_CHECK_EQUAL(collector[1].inclusive.count(), 40);
    BOOST_CHECK_EQUAL(collector[1].exclusive.count(), 40);
    BOOST_CHECK_EQUAL(collector[1].calls, 2);
    BOOST_CHECK_EQUAL(collector[1].progress_events, 4);
    BOOST_CHECK_EQUAL(collector[1].completed_work, 5);
    BOOST_CHECK_EQUAL(collector[1].first_begin.count(), 10);
    BOOST_CHECK_EQUAL(collector[1].first_end.count(), 40);
    BOOST_CHECK_EQUAL(collector[1].last_end.count(), 60);
    check_totals(collector);
}

BOOST_AUTO_TEST_CASE(labels_transform_lengths_metadata_and_outcomes_remain_distinct) {
    ns time {};
    collector_type collector(16, 8, manual_clock {&time});
    for (std::size_t label : {7, 8}) {
        auto root = pa::observe_polynomial<stage::invocation>(
            collector, [=]() noexcept { return pa::polynomial_metadata {{{metric::invocation_label, label}}}; });
        for (std::size_t length : {9, 18, 9}) {
            auto child = pa::observe_polynomial<stage::fft_forward>(collector, [=]() noexcept {
                return pa::polynomial_metadata {
                    {{metric::transform_length, length}, {metric::input_coefficients, label}}};
            });
            time += ns(10);
        }
        {
            auto child = pa::observe_polynomial<stage::x_norm_factor>(collector, empty_metadata);
            child.reject(reason::x_not_square);
            time += ns(5);
        }
        {
            auto child = pa::observe_polynomial<stage::try_cyclic_remainder>(collector, empty_metadata);
            child.decline();
        }
        root.callback_stop();
    }
    BOOST_REQUIRE_EQUAL(collector.size(), 10);
    BOOST_CHECK_EQUAL(collector[1].group.value, 9);
    BOOST_CHECK_EQUAL(collector[1].calls, 2);
    BOOST_CHECK_EQUAL(collector[2].group.value, 18);
    BOOST_CHECK_EQUAL(collector[3].first_rejection.count(), 35);
    BOOST_CHECK_EQUAL(collector[3].rejections[static_cast<std::size_t>(reason::x_not_square)], 1);
    BOOST_CHECK_EQUAL(collector[4].outcomes[static_cast<std::size_t>(outcome::declined)], 1);
    BOOST_CHECK_EQUAL(collector[0].outcomes[static_cast<std::size_t>(outcome::callback_stopped)], 1);
    BOOST_CHECK_EQUAL(collector[5].group.value, 8);
    check_totals(collector);
}

BOOST_AUTO_TEST_CASE(capacity_exhaustion_preserves_balance_and_rolls_time_into_the_parent) {
    for (const auto limits : {std::pair {1, 8}, std::pair {8, 1}, std::pair {0, 0}}) {
        ns time {};
        collector_type collector(limits.first, limits.second, manual_clock {&time});
        {
            auto root = pa::observe_polynomial<stage::invocation>(collector, empty_metadata);
            time = ns(10);
            {
                auto child = pa::observe_polynomial<stage::gcd>(collector, empty_metadata);
                auto grandchild = pa::observe_polynomial<stage::division>(collector, empty_metadata);
                grandchild.progress(1, 1);
                time = ns(30);
            }
            time = ns(50);
        }
        BOOST_CHECK(collector.balanced());
        BOOST_CHECK_EQUAL(collector.dropped_scopes(), limits.first == 0 ? 3 : 2);
        if (collector.size()) {
            BOOST_CHECK_EQUAL(collector[0].inclusive.count(), 50);
            BOOST_CHECK_EQUAL(collector[0].exclusive.count(), 50);
        }
        collector.reset();
        BOOST_CHECK_EQUAL(collector.size(), 0);
        BOOST_CHECK_EQUAL(collector.dropped_scopes(), 0);
        BOOST_CHECK_EQUAL(collector.elapsed().count(), 0);
    }
}

BOOST_AUTO_TEST_CASE(exceptions_balance_collection_and_reset_cannot_invalidate_an_active_scope) {
    ns time {};
    collector_type collector(8, 8, manual_clock {&time});
    BOOST_CHECK_THROW(([&] {
                          auto root = pa::observe_polynomial<stage::invocation>(collector, empty_metadata);
                          BOOST_CHECK_THROW(collector.reset(), std::logic_error);
                          auto child = pa::observe_polynomial<stage::square>(collector, empty_metadata);
                          time = ns(25);
                          throw std::domain_error("arithmetic failure");
                      }()),
                      std::domain_error);
    BOOST_REQUIRE_EQUAL(collector.size(), 2);
    for (std::size_t i = 0; i < collector.size(); ++i)
        BOOST_CHECK_EQUAL(collector[i].outcomes[static_cast<std::size_t>(outcome::exception)], 1);
    check_totals(collector);
}

BOOST_AUTO_TEST_CASE(csv_is_caller_owned_rectangular_and_stream_errors_stay_outside_callbacks) {
    ns time {};
    collector_type collector(8, 8, manual_clock {&time});
    {
        auto scope = pa::observe_polynomial<stage::x_norm_factor>(
            collector, []() noexcept { return pa::polynomial_metadata {{{metric::factor_degree, 3}}}; });
        scope.reject(reason::x_not_square);
    }
    std::ostringstream stream;
    collector.write_csv(stream);
    std::istringstream lines(stream.str());
    std::string line;
    std::size_t records = 0;
    while (std::getline(lines, line)) {
        BOOST_CHECK_EQUAL(std::count(line.begin(), line.end(), ','), 19);
        ++records;
    }
    BOOST_CHECK_EQUAL(records, 6);    // Header, scope, outcome, rejection, metric, summary.
    struct failing_buffer : std::streambuf {
        int_type overflow(int_type) override {
            return traits_type::eof();
        }
    } buffer;
    std::ostream failed(&buffer);
    failed.exceptions(std::ios::badbit | std::ios::failbit);
    BOOST_CHECK_THROW(collector.write_csv(failed), std::ios_base::failure);
    BOOST_CHECK(collector.balanced());
}

namespace {
    using field = fields::babybear;
    using value = field::value_type;
    using schoolbook = pa::schoolbook_backend<value>;
    using mixed = pa::mixed_radix_backend<field>;
    using polynomial = schoolbook::polynomial_type;

    value nonresidue(value start = value(2)) {
        while (start.is_square())
            start += value::one();
        return start;
    }

    template<typename Backend, bool Kernels>
    void check_integration(Backend backend) {
        pa::polynomial_observation_collector<true, Kernels> collector(2048);
        pa::polynomial_context observed(backend, collector);
        pa::polynomial_context plain(backend);
        // Repeated even and odd factors, including an irreducible quadratic needing modular root recovery.
        polynomial input {value(25)};
        for (const auto &factor : {polynomial {-value(4), 1}, polynomial {-value(9), 1}, polynomial {-value(9), 1},
                                   polynomial {1, -value(79), 1}})
            schoolbook {}.multiply(input, input, factor);
        for (bool stop : {false, true}) {
            nil::crypto3::random::algebraic_engine<field> first(173), second(173);
            std::vector<math::polynomial_factor<polynomial>> a, b;
            auto callback = [&](auto &factors, const auto &factor) {
                factors.push_back(factor);
                return stop ? math::factorization_control::stop_factorization :
                              math::factorization_control::continue_factorization;
            };
            const auto expected = math::complete_factorization(input, plain, first,
                                                               [&](const auto &factor) { return callback(a, factor); });
            const auto actual = math::complete_factorization(input, observed, second,
                                                             [&](const auto &factor) { return callback(b, factor); });
            BOOST_CHECK(actual == expected);
            BOOST_CHECK(a == b);
            BOOST_CHECK(first() == second());
            check_totals(collector);
            collector.reset();
        }
        const auto expected = [&] {
            nil::crypto3::random::algebraic_engine<field> generator(173);
            return math::recover_polynomial_x_norm_representation(input, plain, generator);
        }();
        nil::crypto3::random::algebraic_engine<field> generator(173);
        const auto actual = math::recover_polynomial_x_norm_representation(input, observed, generator);
        BOOST_REQUIRE(actual && expected);
        BOOST_CHECK(actual->p == expected->p && actual->q == expected->q);
        BOOST_CHECK_EQUAL(calls(collector, stage::x_norm_recovery), 1);
        BOOST_CHECK_GT(calls(collector, stage::rational_reconstruction), 0);
        check_totals(collector);
        collector.reset();
        // A necessary coefficient filter rejects before any factorization work starts.
        BOOST_CHECK(
            !math::recover_polynomial_x_norm_representation(polynomial {nonresidue(), value(1)}, observed, generator));
        BOOST_CHECK_EQUAL(calls(collector, stage::complete_factorization), 0);
        check_totals(collector);
        collector.reset();
        polynomial rejected;
        const auto c = nonresidue(), d = nonresidue(c + value::one());
        schoolbook {}.multiply(rejected, polynomial {-c, 1}, polynomial {-d, 1});
        BOOST_CHECK(!math::recover_polynomial_x_norm_representation(rejected, observed, generator));
        BOOST_CHECK_EQUAL(calls(collector, stage::x_norm_factor), 1);
        BOOST_CHECK_EQUAL(calls(collector, stage::square_root_preparation), 0);
        check_totals(collector);
    }
}    // namespace

BOOST_AUTO_TEST_CASE(collector_integrates_with_schoolbook_and_mixed_pipelines_at_both_levels) {
    check_integration<schoolbook, false>(schoolbook {});
    check_integration<schoolbook, true>(schoolbook {});
    check_integration<mixed, false>(mixed(32));
    check_integration<mixed, true>(mixed(32));
}

BOOST_AUTO_TEST_CASE(independent_workers_have_independent_collectors_and_scope_ids) {
    ns first_time {}, second_time {};
    collector_type first(8, 8, manual_clock {&first_time}), second(8, 8, manual_clock {&second_time});
    {
        auto a = pa::observe_polynomial<stage::invocation>(first, empty_metadata);
        {
            auto b = pa::observe_polynomial<stage::invocation>(second, empty_metadata);
            first_time = ns(25);
            second_time = ns(10);
        }
        first_time = ns(50);
    }
    BOOST_CHECK_EQUAL(first[0].inclusive.count(), 50);
    BOOST_CHECK_EQUAL(second[0].inclusive.count(), 10);
    BOOST_CHECK_EQUAL(first[0].parent_id, 0);
    BOOST_CHECK_EQUAL(second[0].parent_id, 0);
    check_totals(first);
    check_totals(second);
}

BOOST_AUTO_TEST_CASE(repeated_events_keep_bounded_rows_and_metadata_overflow_is_reported) {
    ns time {};
    collector_type collector(1, 1, manual_clock {&time});
    for (std::size_t i = 1; i <= 13; ++i) {
        auto scope = pa::observe_polynomial<stage::gcd>(
            collector, [=]() noexcept { return pa::polynomial_metadata {{{static_cast<metric>(i), i}}}; });
        time += ns(1);
    }
    BOOST_CHECK_EQUAL(collector.size(), 1);
    BOOST_CHECK_EQUAL(collector.dropped_metadata(), 1);
    for (std::size_t i = 0; i < 1000; ++i) {
        auto scope = pa::observe_polynomial<stage::gcd>(collector, empty_metadata);
        time += ns(1);
    }
    BOOST_CHECK_EQUAL(collector.size(), 1);
    BOOST_CHECK_EQUAL(collector[0].calls, 1013);
    BOOST_CHECK_EQUAL(collector.dropped_scopes(), 0);
    BOOST_CHECK_EQUAL(collector[0].inclusive.count(), 1013);
    BOOST_CHECK(collector.balanced());
}
