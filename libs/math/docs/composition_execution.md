# Application-controlled composition evaluation

Prepared modular composition can delegate coefficient evaluation to an optional
polynomial-backend operation. Crypto3 prepares operands, allocates outputs and
finalizes the result on the calling thread. The default backend evaluates every
coefficient panel serially. Thread ownership, scheduling and resource limits
belong to an application that supplies this optional operation.

The extension applies to the existing prepared coefficient path, including
grouped composition, batched Horner evaluation and single-block evaluation.
Sparse blocks, unsupported fields and unavailable coefficient caches retain
their ordinary arithmetic paths. Product preparation and grouped accumulation
remain separate backend capabilities.

## Backend contract

A backend may define:

```cpp
template<typename Evaluation>
void evaluate_composition_coefficients(const Evaluation &evaluation);
```

`polynomial_context` detects this operation at compile time. A backend without
it continues to satisfy `PolynomialBackend` and uses the serial evaluation.
There is no runtime executor pointer, thread count or scheduler in the library.

The supplied evaluation has two operations needed by the backend:

- `panel_count()` gives the number of independent coefficient panels.
- `evaluate_panels(first, end)` evaluates a half-open range of panel indices.
  It executes serial arithmetic and does not resize the outputs or workspace.

The backend must visit every panel exactly once, using disjoint ranges for
concurrent calls. An empty range is valid. An invalid range throws
`std::out_of_range` before writing any coefficients. Panel boundaries preserve
the selected scalar, SIMD4 or SIMD8 kernel, including partial final panels.

The evaluation borrows immutable coefficient-cache and prepared-operand data,
and writes separate output coefficients. It is valid only for the duration of
the backend operation. The backend must complete every access before returning
or throwing; all asynchronous work must be joined even when a worker fails.
It must not retain the evaluation or call its `finish()` method. The composition
caller condenses outputs after successful completion. On failure, temporary
outputs are discarded and the exception propagates.

This does not make a whole polynomial context concurrently reusable. Other
operations can mutate backend scratch storage. Use distinct contexts for
concurrent polynomial operations; within one prepared evaluation only the
documented panel ranges may run concurrently.

## Direct prepared-block evaluation

Callers already using the prepared coefficient cache can separate its serial
phases explicitly:

```cpp
const auto evaluation = cache.prepare_blocks(outputs, outer, first_block,
                                             block_count, scratch);
evaluation.evaluate_panels(0, evaluation.panel_count());
evaluation.finish();
```

The cache must be available, the requested blocks must exist, and `outer` must
not alias an output. The cache, workspace and output array must remain alive
and in place throughout evaluation. Finalization, workspace reuse and storage
destruction must wait for all range evaluations to complete. The original
`evaluate_blocks` entry point performs these three phases serially.

## Application integration

An application can extend its existing backend while inheriting its arithmetic
and prepared-product capabilities. For example, with an application-owned pool
whose `run` method joins every job before returning or rethrowing an exception:

```cpp
template<typename Backend, typename Pool>
class application_backend : public Backend {
public:
    application_backend(Backend backend, Pool &pool) :
        Backend(std::move(backend)), pool_(pool) {}

    template<typename Evaluation>
    void evaluate_composition_coefficients(const Evaluation &evaluation) {
        pool_.run(evaluation.panel_count(), [&](std::size_t first, std::size_t end) {
            evaluation.evaluate_panels(first, end);
        });
    }

private:
    Pool &pool_;
};
```

The pool must outlive every context that borrows it. Keep it alive across
compositions and recoveries, and use one application worker budget across
operations. A backend may choose serial execution for small batches without
changing the mathematical result. Each panel retains the original arithmetic
order; only independent output coefficients are scheduled differently.

## Correctness coverage

The regular `polynomial_composition` test uses a small application backend to
exercise disjoint panel evaluation with 1/2/4/8 threads. It compares results
with serial evaluation and reference arithmetic, including sparse and constant
blocks, partial SIMD panels, grouped composition, aliasing, exceptions and
workspace reuse. These checks run with benchmark tests disabled. Thread support
is linked only to the test executable; applications provide their own executor.
