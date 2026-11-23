# Application-controlled execution of prepared modular composition

An application can schedule two parts of prepared modular composition on its
own executor: building the coefficient cache and evaluating outer polynomial
blocks with that cache. Each hook is optional and defaults to serial execution
on the calling thread. Crypto3 allocates shared storage before dispatch and
finalizes the result after all panels complete.

The extension applies to the existing prepared coefficient path, including
grouped composition, batched Horner evaluation and single-block evaluation.
Sparse blocks, unsupported fields and unavailable coefficient caches retain
their ordinary arithmetic paths. Product preparation and grouped accumulation
remain separate backend capabilities.

## Backend contract

A backend can implement either or both of these operations.

```cpp
template<typename Preparation>
void prepare_composition_coefficients(const Preparation &preparation);

template<typename Evaluation>
void evaluate_composition_coefficients(const Evaluation &evaluation);
```

`polynomial_context` detects each operation at compile time. During cache
construction, `prepare_panels(first, end)` writes complete coefficient rows or
SIMD panels from immutable powers into preallocated storage. The descriptor
reports `panel_count()`, `coefficient_count()` and `term_count()` so the
application can choose a worker count. Crypto3 skips this hook when the cache
is disabled, unsupported or exceeds its byte budget.

The supplied evaluation has two operations needed by the backend:

- `panel_count()` gives the number of independent coefficient panels.
- `evaluate_panels(first, end)` evaluates a half-open range of panel indices.
  It executes serial arithmetic and does not resize the outputs or workspace.

For either hook, a successful dispatch must visit every panel exactly once.
Concurrent workers must use disjoint ranges and finish before the hook returns
or throws. Both descriptors borrow their storage for the duration of the hook.
An empty range is valid. An invalid range throws
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

This example adds both hooks to the BN254 small-prime FFT backend. The
application supplies a pool whose `run` method joins every job before
returning or propagating an exception.

```cpp
using fft_backend = nil::crypto3::math::polynomial_arithmetic::bn254_small_prime_fft_backend;

template<typename Pool>
class application_backend : public fft_backend {
public:
    application_backend(std::size_t max_transform_size, Pool &pool) :
        fft_backend(max_transform_size), pool_(pool) {}

    bool supports_product_size(std::size_t count) const {
        return fft_backend::supports_product_size(count);
    }

    template<typename Preparation>
    void prepare_composition_coefficients(const Preparation &preparation) {
        pool_.run(preparation.panel_count(), [&](std::size_t first, std::size_t end) {
            preparation.prepare_panels(first, end);
        });
    }

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
operations. A backend may choose serial execution for small batches. Each
panel retains the original arithmetic order.

The adapter forwards `supports_product_size()` because it uses the base
backend's multiplication plan. Half-GCD rational reconstruction checks this
query on the concrete backend. An adapter with tighter limits should report
its own supported sizes or omit the query to use classical reconstruction.

## Correctness coverage

The regular `polynomial_composition` test uses a small application backend to
exercise disjoint panel evaluation with 1/2/4/8 threads. It compares results
with serial evaluation and reference arithmetic, including sparse and constant
blocks, partial SIMD panels, grouped composition, aliasing, exceptions and
workspace reuse. Preparation tests cover padding, cache budgets, disjoint
writes, and retry after a worker throws. These checks run with benchmark tests
disabled. The test application supplies its own executor.
