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

#ifndef CRYPTO3_MATH_POLYNOMIAL_OBSERVATION_HPP
#define CRYPTO3_MATH_POLYNOMIAL_OBSERVATION_HPP

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <type_traits>
#include <utility>

namespace nil::crypto3::math::polynomial_arithmetic {

    // Stable identifiers: algorithm stages occupy 1..0x0fff; backend/kernel stages start at 0x1000.
    // An invocation is a caller-defined root, enabled when either observation level is enabled.
    enum class polynomial_stage : std::uint16_t {
        invocation = 0,
        gcd = 0x0100,
        euclidean_step = 0x0101,
        half_gcd = 0x0102,
        half_gcd_recursion = 0x0103,
        half_gcd_basecase = 0x0104,
        division = 0x0120,
        basecase_division = 0x0121,
        divisor_preparation = 0x0122,
        inverse_series = 0x0123,
        quotient_estimation = 0x0124,
        remainder_reconstruction = 0x0125,
        exact_division = 0x0126,
        factorization_exact_quotient = 0x0127,
        mulmod = 0x0140,
        squaremod = 0x0141,
        powmod = 0x0142,
        powmod_x = 0x0143,
        multiply_by_x = 0x0144,
        frobenius_preparation = 0x0150,
        initial_frobenius_power = 0x0151,
        frobenius_map = 0x0152,
        frobenius_iterations = 0x0153,
        composition_cache = 0x0160,
        composition_cached_powers = 0x0161,
        compose_mod_cached = 0x0162,
        composition_linear_combination = 0x0163,
        composition_giant_step = 0x0164,
        composition_finalization = 0x0165,
        compose_mod = 0x0166,
        compose_mod_reference = 0x0167,
        complete_factorization = 0x0200,
        degree_group_collection = 0x0201,
        square_free_factorization = 0x0220,
        square_free_normalization = 0x0221,
        derivative = 0x0222,
        multiplicity = 0x0223,
        factor_callback = 0x0224,
        distinct_degree_factorization = 0x0240,
        kaltofen_shoup_factor = 0x0241,
        kaltofen_shoup_preparation = 0x0242,
        frobenius_baby_steps = 0x0243,
        degree_block = 0x0244,
        frobenius_giant_step = 0x0245,
        interval_product = 0x0246,
        coarse_gcd = 0x0247,
        degree_refinement = 0x0248,
        classified_factor_removal = 0x0249,
        distinct_degree_reference = 0x024a,
        degree_group_callback = 0x024b,
        equal_degree_factorization = 0x0260,
        equal_degree_group = 0x0261,
        equal_degree_splitting = 0x0262,
        splitting_subproblem = 0x0263,
        splitting_attempt = 0x0264,
        random_candidate = 0x0265,
        subgroup_split = 0x0266,
        irreducible_factor = 0x0267,
        x_norm_recovery = 0x0300,
        x_norm_shortcut = 0x0301,
        x_norm_constant_filter = 0x0302,
        x_norm_leading_filter = 0x0303,
        x_norm_factor = 0x0304,
        x_norm_multiplicity = 0x0305,
        x_norm_irreducible = 0x0306,
        x_norm_square_test = 0x0307,
        x_norm_scalar_check = 0x0308,
        x_norm_normalization = 0x0309,
        x_norm_exact_check = 0x030a,
        x_norm_combination = 0x030b,
        x_norm_combination_level = 0x030c,
        x_norm_multiply = 0x030d,
        x_norm_leading_scalar = 0x030e,
        x_norm_final_verification = 0x030f,
        x_norm_evaluation = 0x0310,
        x_norm_factor_power = 0x0311,
        square_test = 0x0340,
        square_root_preparation = 0x0341,
        square_root_order = 0x0342,
        nonresidue_search = 0x0343,
        nonresidue_attempt = 0x0344,
        nonresidue_power = 0x0345,
        square_root = 0x0346,
        square_root_initial_powers = 0x0347,
        square_root_order_search = 0x0348,
        square_root_correction = 0x0349,
        rational_reconstruction = 0x0360,
        rational_reconstruction_step = 0x0361,
        rational_reconstruction_normalization = 0x0362,
        multiply = 0x1000,
        square = 0x1001,
        multiply_low = 0x1002,
        prepare_low_product = 0x1003,
        try_multiply_low_prepared = 0x1004,
        prepare_cyclic_remainder = 0x1005,
        try_cyclic_remainder = 0x1006,
        low_product_fallback = 0x1007,
        buffer_preparation = 0x1010,
        output_finalization = 0x1011,
        pointwise_multiply = 0x1012,
        pointwise_square = 0x1013,
        cyclic_folding = 0x1014,
        coefficient_subtraction = 0x1015,
        spectrum_compaction = 0x1016,
        fft_forward = 0x1100,
        fft_inverse = 0x1101,
        fft_normalization = 0x1102
    };

    enum class polynomial_metric : std::uint8_t {
        none = 0,
        input_coefficients = 1,
        second_input_coefficients = 2,
        output_coefficients = 3,
        modulus_degree = 4,
        iteration = 5,
        invocation_label = 6,
        quotient_coefficients = 7,
        remainder_coefficients = 8,
        inverse_precision = 9,
        multiplicity = 10,
        factor_degree = 11,
        factor_count = 12,
        arithmetic_path = 13,
        half_gcd_cutoff = 14,
        exponent_bits = 15,
        composition_block_size = 16,
        degree_block_size = 17,
        block_index = 18,
        degree_count = 19,
        pending_subproblems = 20,
        numerator_degree_bound = 21,
        denominator_degree_bound = 22,
        numerator_coefficients = 23,
        denominator_coefficients = 24,
        two_adicity = 25,
        cached_order = 26,
        p_coefficients = 27,
        q_coefficients = 28,
        second_p_coefficients = 29,
        second_q_coefficients = 30,
        transform_length = 31,
        prepared_storage_bytes = 32,
        accepted = 33,
        prepared_available = 34,
        low_coefficient_count = 35,
        detailed_backend = 36
    };

    // Diagnostic values for arithmetic_path metadata; these report decisions made by existing dispatch.
    enum class polynomial_arithmetic_path : std::size_t {
        trivial = 0,
        basecase = 1,
        reciprocal = 2,
        low_product = 3,
        cyclic = 4,
        cached_composition = 5,
        uncached_composition = 6,
        gcd_split = 7,
        power_split = 8,
        zero = 9,
        constant = 10,
        characteristic_two = 11,
        indeterminate_norm = 12,
        euler_criterion = 13,
        even_multiplicity = 14,
        odd_multiplicity = 15,
        supplied_nonresidue = 16,
        sampled_nonresidue = 17
    };

    struct polynomial_measurement {
        polynomial_metric metric = polynomial_metric::none;
        std::size_t value = 0;
    };

    // Fixed-size, allocation-free metadata. Unused entries have metric == none, distinct from a measured zero.
    // Values describe sizes or control flow, never coefficients, witnesses, or generator state.
    using polynomial_metadata = std::array<polynomial_measurement, 4>;

    enum class polynomial_scope_outcome : std::uint8_t {
        completed = 0,
        rejected = 1,
        callback_stopped = 2,
        exception = 3,
        // An optional acceleration was unavailable. This is distinct from mathematical rejection or an error.
        declined = 4
    };

    enum class polynomial_rejection_reason : std::uint8_t {
        unspecified = 0,
        not_square = 1,
        constant_candidate = 2,
        unsuccessful_split = 3,
        constant_not_square = 4,
        constant_coefficient_not_square = 5,
        signed_leading_coefficient_not_square = 6,
        x_not_square = 7,
        rational_degree_bound = 8,
        zero_normalization_scalar = 9,
        nonsquare_normalization_scalar = 10,
        odd_factor_recovery = 11,
        leading_scalar_not_square = 12,
        nonresidue_candidate_square = 13,
        rational_recovery = 14
    };

    struct polynomial_scope_token {
        std::uint64_t id = 0;
        std::uint64_t parent_id = 0;
    };

    struct polynomial_scope_begin {
        polynomial_stage stage;
        polynomial_metadata metadata = {};
    };

    struct polynomial_scope_end {
        polynomial_scope_token scope;
        polynomial_stage stage;
        polynomial_scope_outcome outcome;
        polynomial_rejection_reason rejection;
        polynomial_metadata metadata = {};
    };

    struct polynomial_scope_progress {
        polynomial_scope_token scope;
        polynomial_stage stage;
        // Cumulative completed work within this scope; total == 0 means the total is not known.
        // Emit only after the corresponding operation (including its reductions) has completed.
        std::size_t completed;
        std::size_t total;
        polynomial_metadata metadata = {};
    };

    struct no_polynomial_observer {
        static constexpr bool observe_algorithms = false;
        static constexpr bool observe_kernels = false;
    };

    /**
     * A caller-owned observer supplies compile-time observation levels and three noexcept callbacks.
     * on_begin returns a fresh nonzero ID and the nearest active observed parent's ID (zero for a root).
     * The observer maintains this nesting across contexts borrowing it, even if it drops stored records.
     * on_end closes the corresponding scope; on_progress neither opens nor closes a scope.
     *
     * Scopes nest in ordinary stack order. Callback work performed inside an algorithm belongs below that
     * algorithm's active scope. A collector's inclusive time spans begin to end, including observed children;
     * exclusive time subtracts direct children's inclusive times. Do not sum inclusive totals across levels.
     * Disabled children leave their work in the nearest observed ancestor's exclusive total. FFT scopes include
     * buffer preparation, and inverse FFT scopes also include normalization; those phases are nested children.
     *
     * Prepared-operation begin events count attempts. A completed end with accepted == 1 counts acceptance;
     * declined records an unavailable acceleration, not an arithmetic failure. Transform events record actual
     * calls, so an accepted zero-product shortcut need not contain any FFT. Preparation storage metadata counts
     * the retained spectrum allocation, not temporary workspace or process memory.
     *
     * Callbacks must not throw, including on allocation or output failure, and must not reenter the observed
     * context. They may aggregate instead of retaining events. Clocks, sampling, formatting and output are the
     * observer's responsibility; the hooks perform none of them. Use separate contexts and observers per worker.
     * The observer must outlive all borrowing contexts and scopes. Event references are valid only during the
     * callback. A type with both levels disabled need not implement callbacks at all.
     */
    template<typename Observer>
    concept PolynomialObserver =
        requires {
            typename std::bool_constant<Observer::observe_algorithms>;
            typename std::bool_constant<Observer::observe_kernels>;
        } && (!(Observer::observe_algorithms || Observer::observe_kernels) ||
              requires(Observer &observer, const polynomial_scope_begin &begin, const polynomial_scope_end &end,
                       const polynomial_scope_progress &progress) {
                  { observer.on_begin(begin) } noexcept -> std::same_as<polynomial_scope_token>;
                  { observer.on_end(end) } noexcept -> std::same_as<void>;
                  { observer.on_progress(progress) } noexcept -> std::same_as<void>;
              });

    namespace detail {
        template<typename Factory>
        concept PolynomialMetadataFactory = requires(Factory &&factory) {
            { std::forward<Factory>(factory)() } noexcept -> std::same_as<polynomial_metadata>;
        };

        template<PolynomialObserver Observer, polynomial_stage Stage>
        inline constexpr bool observe_polynomial_stage =
            Stage == polynomial_stage::invocation       ? Observer::observe_algorithms || Observer::observe_kernels :
            static_cast<std::uint16_t>(Stage) >= 0x1000 ? Observer::observe_kernels :
                                                          Observer::observe_algorithms;

        // Empty storage for disabled observation; enabled contexts borrow, rather than copy, the observer.
        template<PolynomialObserver Observer, bool Enabled = Observer::observe_algorithms || Observer::observe_kernels>
        struct polynomial_observer_reference {
            polynomial_observer_reference() = default;
            explicit polynomial_observer_reference(Observer &) noexcept {
            }
        };

        template<PolynomialObserver Observer>
        struct polynomial_observer_reference<Observer, true> {
            explicit polynomial_observer_reference(Observer &observer) noexcept : observer(&observer) {
            }

            Observer *observer;
        };
    }    // namespace detail

    /**
     * Stack-bound event scope. finish() or destruction emits exactly one end event, using exception outcome when
     * a new exception is unwinding, even if rejection or a callback stop had previously been recorded. finish()
     * permits a measured phase to end without changing the lifetime of arithmetic temporaries. Close scopes in
     * nesting order; after finish(), further updates (including lazy metadata) are ignored.
     * Metadata factories are lazy and noexcept: disabled scopes never invoke them or query exception state.
     */
    template<PolynomialObserver Observer, polynomial_stage Stage,
             bool Enabled = detail::observe_polynomial_stage<Observer, Stage>>
    class polynomial_observation_scope {
    public:
        static constexpr bool enabled = true;

        explicit polynomial_observation_scope(Observer &observer, const polynomial_metadata &metadata) noexcept :
            observer_(observer), scope_(observer.on_begin({Stage, metadata})),
            uncaught_exceptions_(std::uncaught_exceptions()) {
        }

        polynomial_observation_scope(const polynomial_observation_scope &) = delete;
        polynomial_observation_scope &operator=(const polynomial_observation_scope &) = delete;
        polynomial_observation_scope(polynomial_observation_scope &&) = delete;
        polynomial_observation_scope &operator=(polynomial_observation_scope &&) = delete;

        ~polynomial_observation_scope() noexcept {
            finish();
        }

        void finish() noexcept {
            if (!active_) {
                return;
            }
            if (std::uncaught_exceptions() > uncaught_exceptions_) {
                outcome_ = polynomial_scope_outcome::exception;
                rejection_ = polynomial_rejection_reason::unspecified;
            }
            observer_.on_end({scope_, Stage, outcome_, rejection_, result_metadata_});
            active_ = false;
        }

        void reject(polynomial_rejection_reason reason) noexcept {
            if (active_) {
                outcome_ = polynomial_scope_outcome::rejected;
                rejection_ = reason;
            }
        }

        void decline() noexcept {
            if (active_) {
                outcome_ = polynomial_scope_outcome::declined;
                rejection_ = polynomial_rejection_reason::unspecified;
            }
        }

        void callback_stop() noexcept {
            if (active_) {
                outcome_ = polynomial_scope_outcome::callback_stopped;
                rejection_ = polynomial_rejection_reason::unspecified;
            }
        }

        template<detail::PolynomialMetadataFactory Factory>
        void set_result_metadata(Factory &&metadata) noexcept {
            if (active_) {
                result_metadata_ = std::forward<Factory>(metadata)();
            }
        }

        template<detail::PolynomialMetadataFactory Factory>
        void progress(std::size_t completed, std::size_t total, Factory &&metadata) noexcept {
            if (active_) {
                observer_.on_progress({scope_, Stage, completed, total, std::forward<Factory>(metadata)()});
            }
        }

        void progress(std::size_t completed, std::size_t total = 0) noexcept {
            progress(completed, total, []() noexcept { return polynomial_metadata {}; });
        }

        // Count completed steps only when observation is enabled. No counter exists in the disabled specialization.
        // Unlike progress(completed, total), advance maintains this scope's counter; explicit progress values
        // do not update it. finish freezes the counter together with the other scope updates.
        void advance(std::size_t total = 0) noexcept {
            if (active_) {
                progress(++completed_steps_, total);
            }
        }

        std::size_t completed_steps() const noexcept {
            return completed_steps_;
        }

    private:
        Observer &observer_;
        polynomial_scope_token scope_;
        int uncaught_exceptions_;
        polynomial_scope_outcome outcome_ = polynomial_scope_outcome::completed;
        polynomial_rejection_reason rejection_ = polynomial_rejection_reason::unspecified;
        polynomial_metadata result_metadata_ = {};
        bool active_ = true;
        std::size_t completed_steps_ = 0;
    };

    template<PolynomialObserver Observer, polynomial_stage Stage>
    class polynomial_observation_scope<Observer, Stage, false> {
    public:
        static constexpr bool enabled = false;

        polynomial_observation_scope() = default;
        polynomial_observation_scope(const polynomial_observation_scope &) = delete;
        polynomial_observation_scope &operator=(const polynomial_observation_scope &) = delete;
        polynomial_observation_scope(polynomial_observation_scope &&) = delete;
        polynomial_observation_scope &operator=(polynomial_observation_scope &&) = delete;

        void finish() noexcept {
        }
        void reject(polynomial_rejection_reason) noexcept {
        }
        void decline() noexcept {
        }
        void callback_stop() noexcept {
        }
        template<detail::PolynomialMetadataFactory Factory>
        void set_result_metadata(Factory &&) noexcept {
        }
        template<detail::PolynomialMetadataFactory Factory>
        void progress(std::size_t, std::size_t, Factory &&) noexcept {
        }
        void progress(std::size_t, std::size_t = 0) noexcept {
        }
        void advance(std::size_t = 0) noexcept {
        }
        std::size_t completed_steps() const noexcept {
            return 0;
        }
    };

    // Low-level backends and FFT plans borrow the same observer for a call, without storing it or needing a context.
    // This is also the context's scope factory; disabled stages never evaluate metadata or read exception state.
    template<polynomial_stage Stage, PolynomialObserver Observer, detail::PolynomialMetadataFactory Factory>
    [[nodiscard]] auto observe_polynomial(Observer &observer, Factory &&metadata) noexcept {
        using scope_type = polynomial_observation_scope<Observer, Stage>;
        if constexpr (scope_type::enabled) {
            return scope_type(observer, std::forward<Factory>(metadata)());
        } else {
            return scope_type();
        }
    }

}    // namespace nil::crypto3::math::polynomial_arithmetic

#endif    // CRYPTO3_MATH_POLYNOMIAL_OBSERVATION_HPP
