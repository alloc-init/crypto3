# Performance {#zk_perf}

The repository currently provides runtime test targets for LPC and Pedersen
commitments. It does not contain reproducible benchmark targets for Placeholder
or the R1CS proof systems, so this page does not publish historical timing or
proof-size figures.

@tableofcontents

## Building benchmarks

Configure a Release build with benchmark tests enabled:

```sh
cmake -S . -B build \
  -DBUILD_TESTS=ON \
  -DBUILD_BENCH_TESTS=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build --target zk_runtime_bench_tests --parallel
```

The aggregate target builds:

- `zk_lpc_bench_test`
- `zk_pedersen_bench_test`

Run them through CTest:

```sh
ctest --test-dir build \
  -R '^zk_(lpc|pedersen)_bench_test$' \
  --output-on-failure
```

## LPC workload

`test/bench_test/lpc.cpp` exercises polynomial commitment, proof generation,
and verification with several FRI step-list configurations. The workload uses
large polynomials and is intended for runtime profiling rather than a quick
unit-test pass.

## Interpreting results

Record at least the following when comparing measurements:

- Crypto3 revision and local changes.
- Compiler, compiler version, and optimization flags.
- Build type and enabled profiling/debug options.
- CPU, operating system, and available memory.
- Commitment parameters, polynomial size, and FRI step list.

Do not compare results from Debug and Release builds. Template instantiation,
large compile-time constants, FFTs, multiexponentiation, hashing, and Merkle
tree construction can contribute differently depending on the protocol and
parameter set.
