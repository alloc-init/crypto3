# Coefficient-combination baseline

This mode of the existing benchmark measures the scalar coefficient-block kernel
used by `compose_mod`. Run it from the Crypto3 checkout on the x86 machine:

```bash
python3 libs/math/test/bench_test/run_small_prime_fft_benchmark.py \
  --coefficient-combination --degree 14938 \
  --warmups 1 --trials 3 --timeout 7200
```

It uses one thread and pins Linux execution to the highest permitted CPU by
default; `--cpu N` selects another permitted CPU. The compiler defaults to `g++`
with `-O3 -DNDEBUG -std=c++23 -march=native`. Use `--cxx` and `--boost-include` to
select the compiler and Boost installation. On macOS, for example, use
`--cxx /usr/bin/clang++ --boost-include /opt/homebrew/include`.

The runner first checks 20 small fixtures, then measures the requested size.
`--check-only` runs those fixtures without the large measurement. All builds and
results go into a new directory under `build/coefficient-combination-benchmarks`.
The final `RESULT_ARCHIVE` contains the logs, CSV samples, summary, input hashes,
source hashes and diff, exact compiler command, CPU information, and thread
settings. Copy that archive back for analysis. It excludes the executable.

To transfer the current source without committing or publishing the branch:

```bash
python3 libs/math/test/bench_test/run_small_prime_fft_benchmark.py \
  --coefficient-combination --package /tmp/coefficient-combination-unit1.tar.gz
```

Extract the archive on x86 and run the same command from its
`crypto3-coefficient-combination-benchmark` directory. The runner verifies the
exported source hashes before building. The destination archive must be new.

## Workload and timing

The fixture is deterministic and uses all twelve coordinates of BN254 Fp12.
The exported `h`, `a`, and `b` files are the modulus, outer polynomial, and inner
polynomial, respectively. The current small-prime backend constructs the real
composition cache once. Subsequent trials reuse it and the same outer polynomial.

At degree 14,938 the block size is 123, there are 122 blocks, and the dense kernel
performs 221,321,408 Fp12 products. The constant power is added directly; the final
block is shorter. Blocks run in descending order. The first accumulator is zero;
later blocks start from a fixed dense polynomial representing an existing Horner
accumulator. This isolates the kernel without computing the giant-step products.

- `setup` rows separate fixture generation/export, backend preparation, divisor
  preparation, and composition-cache construction.
- Each `warmup` or `trial` row sums the elapsed wall and process CPU times of all
  block-kernel calls. Accumulator resets and checks happen outside these timers.
- Every block is checked coefficient by coefficient against a reverse-order
  dot-product reference. Small fixtures also use independent schoolbook/Horner
  composition. `verification` rows report this extra work separately; checking
  the full-size matrix adds substantial runtime.
- Results and references use reusable single-polynomial buffers. No complete
  matrix of output blocks is retained.

`summary.json` reports warm medians and ranges. `process_peak_rss_bytes` is the
process-wide lifetime high-water mark, including setup and verification; zero
means unavailable. `cache_coefficient_bytes` counts stored baby-step and giant-step
coefficients, excluding allocator overhead and backend buffers.
`output_and_reference_capacity_bytes` counts the two reusable result buffers.

This is a component baseline on synthetic inputs. It does not measure a whole
composition or recovery workflow, and the historical 323–325 second observations
came from a different concurrent run. Future candidates must use these same
inputs and timing boundaries and account for any preparation of each new outer
polynomial when reporting their speedups. No observer infrastructure is required
by this benchmark; detailed observation belongs on the instrumentation branch.
