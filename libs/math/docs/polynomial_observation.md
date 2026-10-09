# Polynomial observation {#math_polynomial_observation}

@tableofcontents

Polynomial observation measures the existing arithmetic, factorization, and `X`-norm recovery algorithms. It does
not select algorithms, change cutoffs, alter generator consumption or callback ordering, or disable prepared
operands. Existing `polynomial_context<Backend>` callers use the default no-op observer.

## Attaching an observer

The event contract is in `<nil/crypto3/math/polynomial/polynomial_observation.hpp>`. The optional wall-clock
collector is in `<nil/crypto3/math/polynomial/polynomial_observation_collector.hpp>`; algorithms do not include it.

```cpp
namespace pa = nil::crypto3::math::polynomial_arithmetic;
namespace math = nil::crypto3::math;

// The booleans enable algorithm scopes and backend/kernel scopes, respectively.
using observer_type = pa::polynomial_observation_collector<true, true>;
observer_type observer(2048, 128); // Maximum aggregation rows and active scope depth.
pa::polynomial_context<backend_type, observer_type> context(backend, observer);

{
    auto invocation = context.observe<pa::polynomial_stage::invocation>([]() noexcept {
        return pa::polynomial_metadata {{{pa::polynomial_metric::invocation_label, 42}}};
    });
    auto factors = math::complete_factorization(H, context, coefficient_generator);
}
observer.write_csv(output_stream);
observer.reset(); // Between invocations; retains allocated storage and starts a new clock epoch.

{
    auto invocation = context.observe<pa::polynomial_stage::invocation>([]() noexcept {
        return pa::polynomial_metadata {{{pa::polynomial_metric::invocation_label, 43}}};
    });
    auto representation = math::recover_polynomial_x_norm_representation(H, context, coefficient_generator);
}
observer.write_csv(output_stream);
```

`<true, false>` selects algorithm scopes, `<false, true>` selects kernel scopes, and `<true, true>` selects both.
These are compile-time choices. To disable observation, use the existing `polynomial_context<Backend>`; no collector
is constructed. Its empty scopes do not evaluate metadata factories, read clocks or exception state, allocate,
count events, lock, format, or perform I/O. Schoolbook and custom backends need no observer methods. Context-level
operation scopes still cover custom backends; detailed kernel phases require their optional observer overloads.

The [buildable example](../example/polynomial_observation.cpp) observes both public operations. It uses caller labels,
a caller-owned CSV stream, and a separate progress stream. Enable examples with `-DBUILD_EXAMPLES=ON` and build
`algebra_polynomial_observation_example`. The collector tests are registered as
`crypto3_math_polynomial_observation_collector_test`.

## Lifetime, concurrency, and exceptions

The context owns its backend and borrows the observer. Context copies borrow the same observer. Prepared divisor,
composition, and square-root data do not own an observer: later operations use the current arithmetic context.

An observer must outlive every borrowing context and active scope. The supplied collector cannot be copied or moved.
Use a separate observer, context, generator, and output sink inside each concurrent worker; neither callbacks nor
collector inspection/reset/output are synchronized. An application can collect into worker-local string streams and
combine them after joining its workers. Do not inspect a collector from another thread while it is receiving events.

A custom observer supplies two static booleans, `observe_algorithms` and `observe_kernels`, and:

* `on_begin(const polynomial_scope_begin&) noexcept` returns a fresh nonzero `polynomial_scope_token` with the
  nearest active observed parent's ID (zero for a root).
* `on_end(const polynomial_scope_end&) noexcept` closes that token with its outcome and optional rejection reason.
* `on_progress(const polynomial_scope_progress&) noexcept` reports cumulative completed work within that scope;
  `total == 0` means unknown. Event references are valid only during the callback.

Callbacks must not throw or reenter the observed context. Handle allocation/output failures within a custom observer.
Scope destructors report exception unwinding and cannot mask an arithmetic exception. Rejection, callback-requested
stopping, and declined optional acceleration are separate outcomes. A type with both flags false needs no callbacks.
The event hooks contain no clocks or presentation code. The optional collector uses `std::chrono::steady_clock`;
a custom monotonic clock with nonthrowing `now()` can be supplied as its third template argument for testing.

## Timing and bounded aggregation

Inclusive time covers a scope and its children. Exclusive time subtracts the inclusive times of its direct recorded
children. Sum exclusive totals, or compare scopes at a single level; summing inclusive totals across levels counts
nested work repeatedly. Disabled children remain part of their nearest observed ancestor's exclusive time.

Factorization includes work performed by its norm-recovery callbacks. In each square-free component, the complete
factorization algorithm collects degree groups before invoking its irreducible-factor callbacks. Observation preserves
that ordering. Norm work is visibly nested under `factor_callback`, so a report can distinguish factor discovery from
subsequent recovery. FFT scopes include buffer preparation; inverse FFT scopes also include normalization as a child.

The collector preallocates its row and stack storage at construction. Callbacks do not allocate or perform I/O. Rows
aggregate by parent row and stage, with separate rows for caller invocation labels and transform lengths. Operand
sizes and other metadata are summarized by sample count, minimum, and maximum for each event phase. Loop indices
therefore do not create an unbounded trace. Progress counters sum newly completed work within each invocation,
without counting a repeated cumulative progress report twice.

`size()` and `operator[]` expose the retained `polynomial_observation_totals`. Inspect `balanced()`,
`dropped_scopes()`, and `dropped_metadata()` when accepting a report. Each row has twelve metadata-summary slots;
excess distinct metric/phase pairs are counted as dropped metadata. Exhausting row or depth capacity drops aggregation
for that subtree while still balancing scope tokens. Its time remains in the nearest retained ancestor's exclusive
time. Increase capacities and rerun when complete detail is needed. `reset()` is permitted only without active scopes.
Use distinct labels for a bounded set of invocations, or write/reset between candidates to keep storage bounded.

Timestamps are integer nanoseconds relative to collector construction or reset, not UTC. `first_begin`, `first_end`,
`last_end`, and `first_rejection` describe retained rows; `-1` means an end/rejection has not occurred. To measure time
to the first delivered factor, subtract the invocation's `first_begin` from its first `factor_callback` begin. For
candidate rejection, use the enclosing `x_norm_recovery` rejection time. Rejected random splitting/nonresidue attempts
are normal search outcomes and do not by themselves mean that norm recovery rejected the input.

The collector does not retain every invocation's timeline. A custom observer can stream selected high-level begin/end
records while aggregating frequent arithmetic events. The example demonstrates first-factor and final-outcome records
and a roughly thirty-second heartbeat from completed-work progress. Sampling is an application policy; a long single
operation without a completed-work boundary can postpone a heartbeat. Stream errors are caught inside its callbacks.
Output work also contributes to the observed run's cost.

## CSV records

`write_csv(std::ostream&)` is explicitly invoked by the caller. Stream errors propagate from this method according to
the stream's exception mask, outside arithmetic callbacks. The method does not open files or choose an output path.

The rectangular CSV schema has these record types:

| Record | Meaning |
|---|---|
| `scope` | Hierarchical totals; `id`/`parent_id` identify aggregation rows, not individual event scope tokens. |
| `outcome` | `detail` is an outcome ID and `count` its number of closed scopes. |
| `rejection` | `detail` is a rejection-reason ID and `count` its occurrences. |
| `metric` | `phase` is 0 (begin), 1 (end), or 2 (progress); `detail` is the metric ID; `count`, `min`, `max` summarize it. |
| `summary` | `count` is dropped scopes, `min` dropped metadata, and `max` unbalanced events. |

Only `scope` records should be summed for timing/call totals. Detail rows repeat their scope's columns for context.
Group metadata distinguishes transform lengths and invocation labels. Metric IDs, outcomes, and stage IDs are the
stable enum values catalogued below. A zero measurement differs from a missing metric.

Prepared-operation begins count attempts. A completed end with `accepted == 1` counts acceptance; `declined` reports
unavailable acceleration. FFT events count actual transforms: an accepted zero-product shortcut needs none. Prepared
storage bytes describe the retained spectrum allocation, not peak workspace or process RSS. A fallback scope reports
whether an existing preparation was available when ordinary low-product arithmetic was used.

## Coverage and interpretation limits

| Area | Observed work |
|---|---|
| Square-free/GCD/division | Normalization, derivative, multiplicity, classical and optional half-GCD paths, Euclidean steps, reciprocal preparation/reuse, exact quotients, quotient estimation and remainder reconstruction. |
| Powers/Frobenius | Generic powers and specialized `powmod_x`, squares and multiplication by X, completed exponent steps, initial context setup and subsequent cached composition. |
| Composition | Cache construction and cached powers, coefficientwise block combinations, giant-step products, and finalization. Composition block size and factorization degree-block size are distinct metrics. |
| Distinct-degree factorization | Initial Frobenius setup, baby steps, block composition cache, giant steps, interval products, coarse GCDs, degree refinement, classified-factor removal, and termination. |
| Equal-degree/complete factorization | Random candidate attempts, powers/GCDs, pending subproblems, splits, irreducible outputs, factor callbacks and existing early stops. |
| Norm recovery | Shortcuts and filters, even/odd multiplicities, per-factor tests and recovery, square-root contexts/nonresidue searches/corrections, rational reconstruction, scalar normalization, exact checks, balanced combination, leading scalar, and final norm verification. |
| Backend/FFT | Multiply, dedicated square, low products, preparation/reuse/decline/fallback, buffers, forward/inverse transforms, pointwise operations, normalization, folding/subtraction, compaction, truncation and output finalization. |

Scope entry/exit adds observer overhead to an instrumented run. The collector reports wall time, including scheduling
and callback overhead, not CPU time. There are no per-coefficient or per-field-multiplication events. Scalar work,
allocation/deallocation, dispatch, and copies outside explicit phases remain in an enclosing scope's exclusive time.
FFT-plan construction before attaching the context, and member initialization before a preparation scope opens, need
an enclosing caller scope to be timed and have no finer breakdown. Custom backends without detailed overloads expose
only context operation totals. Peak RSS, process CPU time, and UTC timestamps are left to the application and must be
labelled as process-wide where applicable.

Instrumentation is not a performance optimization. Compare disabled, algorithm-only, and algorithm-plus-kernel builds
on identical inputs, with the same compiler flags and worker count. Separate collector/context construction and cache
preparation from warm operation timing; validate results outside timed regions, report individual trials and medians,
and avoid timing assertions in correctness tests. Small prime-field operations can have substantial relative observer
overhead; expensive extension-field operations need separate measurements. Small fixtures do not establish overhead
for complete large-degree factorization or norm recovery.

## Stage catalog

| ID | Name |
|---|---|
| 0 | `invocation` |
| 256 | `gcd` |
| 257 | `euclidean_step` |
| 258 | `half_gcd` |
| 259 | `half_gcd_recursion` |
| 260 | `half_gcd_basecase` |
| 288 | `division` |
| 289 | `basecase_division` |
| 290 | `divisor_preparation` |
| 291 | `inverse_series` |
| 292 | `quotient_estimation` |
| 293 | `remainder_reconstruction` |
| 294 | `exact_division` |
| 295 | `factorization_exact_quotient` |
| 320 | `mulmod` |
| 321 | `squaremod` |
| 322 | `powmod` |
| 323 | `powmod_x` |
| 324 | `multiply_by_x` |
| 336 | `frobenius_preparation` |
| 337 | `initial_frobenius_power` |
| 338 | `frobenius_map` |
| 339 | `frobenius_iterations` |
| 352 | `composition_cache` |
| 353 | `composition_cached_powers` |
| 354 | `compose_mod_cached` |
| 355 | `composition_linear_combination` |
| 356 | `composition_giant_step` |
| 357 | `composition_finalization` |
| 358 | `compose_mod` |
| 359 | `compose_mod_reference` |
| 512 | `complete_factorization` |
| 513 | `degree_group_collection` |
| 544 | `square_free_factorization` |
| 545 | `square_free_normalization` |
| 546 | `derivative` |
| 547 | `multiplicity` |
| 548 | `factor_callback` |
| 576 | `distinct_degree_factorization` |
| 577 | `kaltofen_shoup_factor` |
| 578 | `kaltofen_shoup_preparation` |
| 579 | `frobenius_baby_steps` |
| 580 | `degree_block` |
| 581 | `frobenius_giant_step` |
| 582 | `interval_product` |
| 583 | `coarse_gcd` |
| 584 | `degree_refinement` |
| 585 | `classified_factor_removal` |
| 586 | `distinct_degree_reference` |
| 587 | `degree_group_callback` |
| 608 | `equal_degree_factorization` |
| 609 | `equal_degree_group` |
| 610 | `equal_degree_splitting` |
| 611 | `splitting_subproblem` |
| 612 | `splitting_attempt` |
| 613 | `random_candidate` |
| 614 | `subgroup_split` |
| 615 | `irreducible_factor` |
| 768 | `x_norm_recovery` |
| 769 | `x_norm_shortcut` |
| 770 | `x_norm_constant_filter` |
| 771 | `x_norm_leading_filter` |
| 772 | `x_norm_factor` |
| 773 | `x_norm_multiplicity` |
| 774 | `x_norm_irreducible` |
| 775 | `x_norm_square_test` |
| 776 | `x_norm_scalar_check` |
| 777 | `x_norm_normalization` |
| 778 | `x_norm_exact_check` |
| 779 | `x_norm_combination` |
| 780 | `x_norm_combination_level` |
| 781 | `x_norm_multiply` |
| 782 | `x_norm_leading_scalar` |
| 783 | `x_norm_final_verification` |
| 784 | `x_norm_evaluation` |
| 785 | `x_norm_factor_power` |
| 832 | `square_test` |
| 833 | `square_root_preparation` |
| 834 | `square_root_order` |
| 835 | `nonresidue_search` |
| 836 | `nonresidue_attempt` |
| 837 | `nonresidue_power` |
| 838 | `square_root` |
| 839 | `square_root_initial_powers` |
| 840 | `square_root_order_search` |
| 841 | `square_root_correction` |
| 864 | `rational_reconstruction` |
| 865 | `rational_reconstruction_step` |
| 866 | `rational_reconstruction_normalization` |
| 4096 | `multiply` |
| 4097 | `square` |
| 4098 | `multiply_low` |
| 4099 | `prepare_low_product` |
| 4100 | `try_multiply_low_prepared` |
| 4101 | `prepare_cyclic_remainder` |
| 4102 | `try_cyclic_remainder` |
| 4103 | `low_product_fallback` |
| 4112 | `buffer_preparation` |
| 4113 | `output_finalization` |
| 4114 | `pointwise_multiply` |
| 4115 | `pointwise_square` |
| 4116 | `cyclic_folding` |
| 4117 | `coefficient_subtraction` |
| 4118 | `spectrum_compaction` |
| 4352 | `fft_forward` |
| 4353 | `fft_inverse` |
| 4354 | `fft_normalization` |

## Metric catalog

| ID | Name |
|---|---|
| 0 | `none` |
| 1 | `input_coefficients` |
| 2 | `second_input_coefficients` |
| 3 | `output_coefficients` |
| 4 | `modulus_degree` |
| 5 | `iteration` |
| 6 | `invocation_label` |
| 7 | `quotient_coefficients` |
| 8 | `remainder_coefficients` |
| 9 | `inverse_precision` |
| 10 | `multiplicity` |
| 11 | `factor_degree` |
| 12 | `factor_count` |
| 13 | `arithmetic_path` |
| 14 | `half_gcd_cutoff` |
| 15 | `exponent_bits` |
| 16 | `composition_block_size` |
| 17 | `degree_block_size` |
| 18 | `block_index` |
| 19 | `degree_count` |
| 20 | `pending_subproblems` |
| 21 | `numerator_degree_bound` |
| 22 | `denominator_degree_bound` |
| 23 | `numerator_coefficients` |
| 24 | `denominator_coefficients` |
| 25 | `two_adicity` |
| 26 | `cached_order` |
| 27 | `p_coefficients` |
| 28 | `q_coefficients` |
| 29 | `second_p_coefficients` |
| 30 | `second_q_coefficients` |
| 31 | `transform_length` |
| 32 | `prepared_storage_bytes` |
| 33 | `accepted` |
| 34 | `prepared_available` |
| 35 | `low_coefficient_count` |
| 36 | `detailed_backend` |

## Outcome catalog

| ID | Name |
|---|---|
| 0 | `completed` |
| 1 | `rejected` |
| 2 | `callback_stopped` |
| 3 | `exception` |
| 4 | `declined` |

## Rejection catalog

| ID | Name |
|---|---|
| 0 | `unspecified` |
| 1 | `not_square` |
| 2 | `constant_candidate` |
| 3 | `unsuccessful_split` |
| 4 | `constant_not_square` |
| 5 | `constant_coefficient_not_square` |
| 6 | `signed_leading_coefficient_not_square` |
| 7 | `x_not_square` |
| 8 | `rational_degree_bound` |
| 9 | `zero_normalization_scalar` |
| 10 | `nonsquare_normalization_scalar` |
| 11 | `odd_factor_recovery` |
| 12 | `leading_scalar_not_square` |
| 13 | `nonresidue_candidate_square` |
| 14 | `rational_recovery` |
