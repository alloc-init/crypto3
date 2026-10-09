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

#include <iostream>

#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/polynomial_observation_collector.hpp>
#include <nil/crypto3/math/polynomial/reconstruction/polynomial_x_norm_reconstruction.hpp>
#include <nil/crypto3/random/algebraic_engine.hpp>

namespace math = nil::crypto3::math;
namespace pa = math::polynomial_arithmetic;

// Sampling and output belong to the application. Only completed-work progress is sampled, so a long single
// operation can postpone the next heartbeat. Every arithmetic event is still aggregated without being printed.
class progress_observer : public pa::polynomial_observation_collector<true, false> {
    using base = pa::polynomial_observation_collector<true, false>;

public:
    explicit progress_observer(std::ostream &progress_output) : output_(progress_output) {
    }

    pa::polynomial_scope_token on_begin(const pa::polynomial_scope_begin &event) noexcept {
        const auto token = base::on_begin(event);
        if (event.stage == pa::polynomial_stage::invocation)
            first_factor_ = true;
        if (event.stage == pa::polynomial_stage::factor_callback && first_factor_) {
            first_factor_ = false;
            write([&] {
                output_ << "first_factor," << token.id << ',' << token.parent_id << ',' << elapsed().count() << '\n';
            });
        }
        return token;
    }

    void on_end(const pa::polynomial_scope_end &event) noexcept {
        base::on_end(event);
        if (event.stage == pa::polynomial_stage::x_norm_recovery ||
            event.stage == pa::polynomial_stage::complete_factorization)
            write([&] {
                output_ << "end," << event.scope.id << ',' << static_cast<unsigned>(event.stage) << ','
                        << elapsed().count() << ',' << static_cast<unsigned>(event.outcome) << ','
                        << static_cast<unsigned>(event.rejection) << '\n';
            });
    }

    void on_progress(const pa::polynomial_scope_progress &event) noexcept {
        base::on_progress(event);
        const auto now = elapsed();
        if (now < next_report_)
            return;
        next_report_ = now + std::chrono::seconds(30);
        write([&] {
            output_ << "progress," << event.scope.id << ',' << static_cast<unsigned>(event.stage) << ',' << now.count()
                    << ',' << event.completed << ',' << event.total << '\n';
            output_.flush();
        });
    }

    bool output_failed() const noexcept {
        return output_failed_;
    }

private:
    template<typename Write>
    void write(Write &&write_record) noexcept {
        if (output_failed_)
            return;
        try {
            write_record();
            output_failed_ = !output_;
        } catch (...) {
            // Observer failures must not change arithmetic results or mask an exception being unwound.
            output_failed_ = true;
        }
    }
    std::ostream &output_;
    duration next_report_ = std::chrono::seconds(30);
    bool first_factor_ = true, output_failed_ = false;
};

int main() {
    using field = nil::crypto3::algebra::fields::babybear;
    using value = field::value_type;
    using backend = pa::schoolbook_backend<value>;
    using polynomial = backend::polynomial_type;

    // Each concurrent worker should construct its own observer, arithmetic context, generator, and output sink.
    progress_observer observer(std::clog);
    pa::polynomial_context context(backend {}, observer);
    nil::crypto3::random::algebraic_engine<field> generator(173);
    const polynomial h {1, -value(79), 1};    // X^2 - 79X + 1 is an irreducible, representable norm over this field.

    // Caller labels identify invocations without adding application-specific concepts to the algorithms.
    {
        auto invocation = context.observe<pa::polynomial_stage::invocation>(
            []() noexcept { return pa::polynomial_metadata {{{pa::polynomial_metric::invocation_label, 1}}}; });
        const auto factors = math::complete_factorization(h, context, generator);
        if (!factors.complete)
            return 1;
    }
    {
        auto invocation = context.observe<pa::polynomial_stage::invocation>(
            []() noexcept { return pa::polynomial_metadata {{{pa::polynomial_metric::invocation_label, 2}}}; });
        const auto representation = math::recover_polynomial_x_norm_representation(h, context, generator);
        if (!representation)
            return 1;
    }
    // Write summaries after both scopes close. CSV contains sizes, outcomes and timings, never coefficients.
    observer.write_csv(std::cout);
    return observer.balanced() && !observer.output_failed() && observer.dropped_scopes() == 0 &&
                   observer.dropped_metadata() == 0 ?
               0 :
               1;
}
