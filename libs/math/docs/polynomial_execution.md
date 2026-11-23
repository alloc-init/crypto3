# Application-controlled small-prime polynomial products

`bn254_small_prime_fft_backend` continues to execute serially. Applications can
use `basic_bn254_small_prime_fft_backend<ExecutionPolicy>` to schedule independent
parts of its products on their own executor. Crypto3 creates no threads, owns no
thread pool and has no runtime thread-count setting. The application controls
worker lifetime, scheduling, concurrency limits and when small jobs run serially.

The existing name is an alias for the template with
`serial_small_prime_fft_execution`. Its constructors, arithmetic bounds, optional
capabilities and arithmetic loop order are retained. Higher-level division, GCD,
modular composition and recovery algorithms use the chosen backend through the
existing `polynomial_context` interface.

## Execution contract

The policy is stored by value and supplies:

```cpp
template<typename Task>
void execute(const Task &task);
```

The borrowed task exposes `stage()`, `task_count()`, `transform_size()` and
`execute_range(first, end)`. The last operation runs serial arithmetic for a
half-open range of **task indices**. An empty range is valid; reversed or
out-of-bounds ranges throw `std::out_of_range` before writes.

The policy must visit every index exactly once on success, using disjoint ranges
for concurrent calls. It must finish every access before returning or throwing,
including joining outstanding workers when another worker fails. It must not
retain the task, call another operation on the same backend while a task is
active, or allow the caller to change the borrowed inputs or output storage.

Tasks use immutable plans and prepared operands, disjoint writes, and local
kernel scratch. Shared spectrum/output allocation and final output swaps remain
on the calling thread; kernel-local scratch may allocate on workers.
These guarantees apply to ranges within one operation; a whole backend/context
still requires sequential use. Concurrent polynomial operations need distinct
backends and contexts.

| `small_prime_fft_execution_stage` | One task index represents |
|---|---|
| `conversion` | Up to 32 coefficient positions, across all primes and coordinates, including zero padding |
| `forward_transform` | One prime/coordinate channel of an operand being prepared |
| `prime_product` | One prime's forward transforms, Fp12 product or square, and inverse transforms |
| `prime_sum` | One prime's fused matrix product-sum and inverse transforms |
| `pointwise_product` | One prime's product for a grouped accumulator append |
| `inverse_transform` | One prime/coordinate channel during grouped-sum finalization |
| `reconstruction` | Up to 32 output coefficients, including CRT reconstruction and field export |

`transform_size()` is the current transform length, not the configured maximum.
The final coefficient task may be shorter than 32 positions. Grouping adjacent
indices in one range preserves the existing conversion and CRT tiling. Ordinary,
squared, prepared and cyclic products share the prime product pipeline. Keeping
each prime's transforms and product together preserves its locality. Matrix
outputs are processed sequentially, with independent prime jobs for each output.

Scalar/linear products, zero paths, polynomial condensation and cyclic field
subtraction retain their serial implementation. Threading within a single FFT
is not introduced. Available parallelism is bounded by the current task count;
prime product stages have one task per selected auxiliary prime.

## Application integration

This example assumes an application pool whose `run` method waits for all work,
propagates worker exceptions after joining, and passes disjoint half-open ranges:

```cpp
namespace arithmetic = nil::crypto3::math::polynomial_arithmetic;

template<typename Pool>
struct product_execution {
    Pool &pool;

    template<typename Task>
    void execute(const Task &task) {
        pool.run(task.task_count(), [&](std::size_t first, std::size_t end) {
            task.execute_range(first, end);
        });
    }
};

template<typename Pool>
auto make_context(Pool &pool, std::size_t max_transform_size) {
    using backend = arithmetic::basic_bn254_small_prime_fft_backend<product_execution<Pool>>;
    return arithmetic::polynomial_context<backend>(
        backend(max_transform_size, product_execution<Pool> {pool}));
}
```

The pool must outlive the contexts that borrow it. Reuse it across products and
recoveries, and coordinate its worker budget with concurrent grinding candidates.
A policy can use the stage, transform length and task count to select serial
execution or a subset of its workers. Prepared objects can cross compatible
backend instances with the same policy type; the application can change its
worker budget without changing that type. Moving the backend preserves its
policy under the policy's usual move semantics.

## Rounding, failure and storage

Floating-point fast-path selection occurs on the executing thread. Workers may
have different rounding modes; each must keep its mode fixed throughout a range
call. Kernels use their existing exact fallback when that thread cannot use the
fast path. Crypto3 does not change a worker's floating-point environment.

Output polynomials are committed after all scheduled work succeeds. If a policy
throws after joining its workers, outputs remain unchanged and backend scratch
can be reused. Prepared operands remain immutable. With a custom execution policy,
a grouped append also keeps its previous logical sum and product count; failed
finalization can be retried.

Those grouped guarantees have costs on the custom-policy path: an append waits
for pointwise products, then adds them to the accumulator serially; finalization
copies the accumulated spectrum into the backend's ordinary product buffer
before inverse transforms. The default serial policy retains the original fused
append and in-place finalization. The opted-in copy and serial addition must be
included in performance measurements, even if an application policy uses one
worker. There is no full spectrum allocation per worker; existing spectrum and
matrix budgets retain their meanings, and worker stacks/pool storage belong to
the application.

## Correctness coverage

The regular `small_prime_fft_backend` tests compare application execution with
schoolbook arithmetic using 1/2/4/8 workers, uneven ranges, coefficient tails,
output aliases, moved backends, repeated preparations, cancellation and grouped
reuse. They exercise differing worker rounding modes, failures in both the
calling thread and helpers, worker completion before exception propagation,
division, half-GCD and grouped composition. Conversion and CRT tests also check
their disjoint range boundaries directly. These tests run with benchmark tests
disabled; only the test application links thread support.

These checks establish arithmetic and lifecycle behavior. Performance estimates
require measurements with an application executor, including preparation,
scheduling, grouped copies, memory and concurrent candidate throughput.
