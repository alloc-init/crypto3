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

#ifndef CRYPTO3_MATH_POLYNOMIAL_OBSERVATION_COLLECTOR_HPP
#define CRYPTO3_MATH_POLYNOMIAL_OBSERVATION_COLLECTOR_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <ostream>
#include <stdexcept>
#include <vector>

#include <nil/crypto3/math/polynomial/polynomial_observation.hpp>

namespace nil::crypto3::math::polynomial_arithmetic {

    enum class polynomial_metadata_phase { begin, end, progress };

    struct polynomial_metric_summary {
        polynomial_metadata_phase phase = polynomial_metadata_phase::begin;
        polynomial_metric metric = polynomial_metric::none;
        std::size_t samples = 0;
        std::size_t minimum = 0;
        std::size_t maximum = 0;
    };

    struct polynomial_observation_totals {
        using duration = std::chrono::nanoseconds;
        std::size_t id = 0, parent_id = 0;
        polynomial_stage stage = polynomial_stage::invocation;
        // Separate caller labels and transform lengths without creating a row per loop index or operand size.
        polynomial_measurement group {};
        std::size_t calls = 0;
        std::array<std::size_t, static_cast<std::size_t>(polynomial_scope_outcome::declined) + 1> outcomes {};
        std::array<std::size_t, static_cast<std::size_t>(polynomial_rejection_reason::rational_recovery) + 1>
            rejections {};
        duration inclusive {}, exclusive {};
        duration first_begin {}, first_end {-1}, last_end {-1}, first_rejection {-1};
        std::size_t progress_events = 0, completed_work = 0;
        // Four fields per event phase; varying field sets share these bounded slots. Overflow is reported.
        std::array<polynomial_metric_summary, 12> metrics {};
    };

    /**
     * Optional, caller-owned wall-clock collector. Include this header explicitly; algorithms do not depend on it.
     * Storage is allocated at construction. Callbacks do not allocate, format, lock, or perform I/O.
     * Rows aggregate by parent row, stage, and invocation label or transform length. Repeated calls reuse rows.
     * Capacity exhaustion drops aggregation, not scope balancing; dropped work stays in its retained ancestor's
     * exclusive time. Inspect dropped_scopes()/dropped_metadata() before treating a report as complete.
     *
     * Clock must be monotonic and now() must not throw. Clock injection supports deterministic timing tests.
     * Keep this object alive and at the same address until all borrowing contexts and scopes are gone.
     * Use one collector and context per worker; callbacks, inspection, reset, and output are not thread-safe.
     * Use no_polynomial_observer for the disabled path, which needs no collector, clock, or allocation.
     */
    template<bool Algorithms = true, bool Kernels = false, typename Clock = std::chrono::steady_clock>
        requires(Algorithms || Kernels)
    class polynomial_observation_collector {
    public:
        static constexpr bool observe_algorithms = Algorithms;
        static constexpr bool observe_kernels = Kernels;
        using duration = polynomial_observation_totals::duration;
        static_assert(Clock::is_steady);
        static_assert(noexcept(std::declval<Clock &>().now()));

        explicit polynomial_observation_collector(std::size_t max_rows = 1024, std::size_t max_depth = 128,
                                                  Clock clock = {}) :
            rows_(max_rows), stack_(max_depth), clock_(std::move(clock)), epoch_(clock_.now()) {
        }
        polynomial_observation_collector(const polynomial_observation_collector &) = delete;
        polynomial_observation_collector &operator=(const polynomial_observation_collector &) = delete;
        polynomial_observation_collector(polynomial_observation_collector &&) = delete;
        polynomial_observation_collector &operator=(polynomial_observation_collector &&) = delete;

        polynomial_scope_token on_begin(const polynomial_scope_begin &event) noexcept {
            const auto now = clock_.now();
            const polynomial_scope_token token {next_id_++, active_id_};
            active_id_ = token.id;
            std::size_t row = 0;
            if (depth_ < stack_.size()) {
                // A missing parent means its entire subtree is unaggregated, rather than attached to a false root.
                const auto parent = depth_ == 0 ? 0 : stack_[depth_ - 1].row;
                if (depth_ == 0 || parent != 0)
                    row = find_row(parent, event);
                stack_[depth_] = {row, now, {}, 0};
            }
            ++depth_;
            if (row == 0) {
                ++dropped_scopes_;
            } else {
                auto &totals = rows_[row - 1].totals;
                if (totals.calls++ == 0)
                    totals.first_begin = since_epoch(now);
                add_metadata(totals, polynomial_metadata_phase::begin, event.metadata);
            }
            return token;
        }

        void on_end(const polynomial_scope_end &event) noexcept {
            const auto now = clock_.now();
            // Valid scopes close in stack order, including when aggregation capacity is exhausted.
            if (depth_ == 0 || event.scope.id != active_id_) {
                ++unbalanced_events_;
                return;
            }
            active_id_ = event.scope.parent_id;
            --depth_;
            if (depth_ >= stack_.size())
                return;
            const auto &frame = stack_[depth_];
            if (frame.row == 0)
                return;
            auto &totals = rows_[frame.row - 1].totals;
            const auto elapsed = std::chrono::duration_cast<duration>(now - frame.start);
            totals.inclusive += elapsed;
            totals.exclusive += elapsed - frame.children;
            ++totals.outcomes[static_cast<std::size_t>(event.outcome)];
            totals.last_end = since_epoch(now);
            if (totals.first_end.count() < 0)
                totals.first_end = totals.last_end;
            if (event.outcome == polynomial_scope_outcome::rejected) {
                ++totals.rejections[static_cast<std::size_t>(event.rejection)];
                if (totals.first_rejection.count() < 0)
                    totals.first_rejection = totals.last_end;
            }
            add_metadata(totals, polynomial_metadata_phase::end, event.metadata);
            if (depth_ != 0)
                stack_[depth_ - 1].children += elapsed;
        }

        void on_progress(const polynomial_scope_progress &event) noexcept {
            if (depth_ == 0 || event.scope.id != active_id_) {
                ++unbalanced_events_;
                return;
            }
            if (depth_ > stack_.size())
                return;
            auto &frame = stack_[depth_ - 1];
            if (frame.row == 0)
                return;
            auto &totals = rows_[frame.row - 1].totals;
            ++totals.progress_events;
            // Progress is cumulative within each invocation. Count only newly completed work, not repeated reports.
            if (event.completed >= frame.completed)
                totals.completed_work += event.completed - frame.completed;
            frame.completed = event.completed;
            add_metadata(totals, polynomial_metadata_phase::progress, event.metadata);
        }

        std::size_t size() const noexcept {
            return row_count_;
        }
        const polynomial_observation_totals &operator[](std::size_t index) const noexcept {
            return rows_[index].totals;
        }
        std::size_t dropped_scopes() const noexcept {
            return dropped_scopes_;
        }
        std::size_t dropped_metadata() const noexcept {
            return dropped_metadata_;
        }
        bool balanced() const noexcept {
            return depth_ == 0 && active_id_ == 0 && unbalanced_events_ == 0;
        }
        duration elapsed() const noexcept {
            return since_epoch(clock_.now());
        }

        // Call between invocations after writing any report that needs to be retained. Allocations are reused.
        void reset() {
            if (depth_ != 0)
                throw std::logic_error("cannot reset a polynomial collector with active scopes");
            std::fill_n(rows_.begin(), row_count_, row_type {});
            row_count_ = first_root_ = dropped_scopes_ = dropped_metadata_ = unbalanced_events_ = 0;
            next_id_ = 1;
            epoch_ = clock_.now();
        }

        // Presentation is caller-triggered. Stream failures propagate here, never from an algorithm callback.
        // One rectangular schema: scope rows hold totals; outcome/rejection/metric rows use detail/count/min/max.
        void write_csv(std::ostream &output) const {
            output << "record,id,parent_id,stage,group_metric,group_value,calls,inclusive_ns,exclusive_ns,"
                      "first_begin_ns,first_end_ns,last_end_ns,first_rejection_ns,progress_events,completed_work,"
                      "phase,detail,count,min,max\n";
            auto prefix = [&](const char *kind, const polynomial_observation_totals &t) {
                output << kind << ',' << t.id << ',' << t.parent_id << ',' << static_cast<unsigned>(t.stage) << ','
                       << static_cast<unsigned>(t.group.metric) << ',' << t.group.value << ',' << t.calls << ','
                       << t.inclusive.count() << ',' << t.exclusive.count() << ',' << t.first_begin.count() << ','
                       << t.first_end.count() << ',' << t.last_end.count() << ',' << t.first_rejection.count() << ','
                       << t.progress_events << ',' << t.completed_work;
            };
            for (std::size_t i = 0; i < row_count_; ++i) {
                const auto &t = rows_[i].totals;
                prefix("scope", t);
                output << ",,,,,\n";
                for (std::size_t j = 0; j < t.outcomes.size(); ++j)
                    if (t.outcomes[j]) {
                        prefix("outcome", t);
                        output << ",," << j << ',' << t.outcomes[j] << ",,\n";
                    }
                for (std::size_t j = 0; j < t.rejections.size(); ++j)
                    if (t.rejections[j]) {
                        prefix("rejection", t);
                        output << ",," << j << ',' << t.rejections[j] << ",,\n";
                    }
                for (const auto &m : t.metrics)
                    if (m.samples) {
                        prefix("metric", t);
                        output << ',' << static_cast<unsigned>(m.phase) << ',' << static_cast<unsigned>(m.metric) << ','
                               << m.samples << ',' << m.minimum << ',' << m.maximum << '\n';
                    }
            }
            // Summary uses the same columns: count=dropped scopes, min=dropped metadata, max=unbalanced events.
            prefix("summary", {});
            output << ",,," << dropped_scopes_ << ',' << dropped_metadata_ << ',' << unbalanced_events_ << '\n';
        }

    private:
        struct row_type {
            polynomial_observation_totals totals;
            std::size_t first_child = 0, next_sibling = 0;
        };
        struct frame_type {
            std::size_t row;
            typename Clock::time_point start;
            duration children;
            std::size_t completed;
        };

        duration since_epoch(typename Clock::time_point now) const noexcept {
            return std::chrono::duration_cast<duration>(now - epoch_);
        }

        std::size_t find_row(std::size_t parent, const polynomial_scope_begin &event) noexcept {
            polynomial_measurement group {};
            for (const auto &m : event.metadata)
                if (m.metric == polynomial_metric::transform_length ||
                    (event.stage == polynomial_stage::invocation && m.metric == polynomial_metric::invocation_label)) {
                    group = m;
                    break;
                }
            auto &first = parent == 0 ? first_root_ : rows_[parent - 1].first_child;
            for (auto id = first; id != 0; id = rows_[id - 1].next_sibling) {
                const auto &t = rows_[id - 1].totals;
                if (t.stage == event.stage && t.group.metric == group.metric && t.group.value == group.value)
                    return id;
            }
            if (row_count_ == rows_.size())
                return 0;
            const auto id = ++row_count_;
            auto &row = rows_[id - 1];
            row.totals.id = id;
            row.totals.parent_id = parent;
            row.totals.stage = event.stage;
            row.totals.group = group;
            row.next_sibling = first;
            first = id;
            return id;
        }

        void add_metadata(polynomial_observation_totals &totals, polynomial_metadata_phase phase,
                          const polynomial_metadata &metadata) noexcept {
            for (const auto &value : metadata) {
                if (value.metric == polynomial_metric::none)
                    continue;
                auto slot = std::find_if(totals.metrics.begin(), totals.metrics.end(), [&](const auto &m) {
                    return m.samples == 0 || (m.phase == phase && m.metric == value.metric);
                });
                if (slot == totals.metrics.end()) {
                    ++dropped_metadata_;
                    continue;
                }
                if (slot->samples++ == 0) {
                    slot->phase = phase;
                    slot->metric = value.metric;
                    slot->minimum = slot->maximum = value.value;
                } else {
                    slot->minimum = std::min(slot->minimum, value.value);
                    slot->maximum = std::max(slot->maximum, value.value);
                }
            }
        }

        std::vector<row_type> rows_;
        std::vector<frame_type> stack_;
        Clock clock_;
        typename Clock::time_point epoch_;
        std::size_t row_count_ = 0, first_root_ = 0, depth_ = 0;
        std::size_t dropped_scopes_ = 0, dropped_metadata_ = 0, unbalanced_events_ = 0;
        std::uint64_t next_id_ = 1, active_id_ = 0;
    };

}    // namespace nil::crypto3::math::polynomial_arithmetic
#endif    // CRYPTO3_MATH_POLYNOMIAL_OBSERVATION_COLLECTOR_HPP
