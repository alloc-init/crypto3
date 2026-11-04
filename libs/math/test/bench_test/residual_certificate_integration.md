# Combined parity and residual-certificate benchmark

This package measures the actual library continuation with **parity alone** and **parity plus certificates**. Both
variants use the same pinned production source snapshot, without observer instrumentation. This is a feature comparison
within one build, not a comparison of two separately built commits. The default comparison covers **all remaining
distinct-degree groups** from each saved checkpoint, with two repetitions in reversed trial order.

The algorithm and caller controls are documented in [polynomial factorization](../../docs/polynomial_factorization.md).
The production default remains disabled. The package allows one attempt starting at the restored boundary; the
production parity gate skips an even-factor residual, and the natural stopping bound takes precedence. These states
were selected for the earlier experiment and do not calibrate a trigger for fresh grinding. An eligible composite
residual pays for its failed certificate and its complete fallback scan.

Before either timed variant, the harness restores parity from the original input's production square-free GCD and
the number of factors already emitted at the checkpoint. The completed factor prefix is checked against the pinned
capture: contiguous indices, completed outcomes, multiplicity one and the exact removed degree. Future residual
factors supply only validation outside the measured interval; they do not determine the production parity or trigger.
All five full-size inputs require available parity, so an unavailable result fails instead of silently changing the
comparison. Small checks also exercise the conservative unknown-parity path.

## What this measures

| Variant | Measured work |
|---|---|
| `disabled` | Parity-aware production continuation, with no certificate attempts |
| `prepared` | The same continuation, allowing one eligible certificate attempt with an additional prepared-cache allowance |
| `scalar` | The same continuation with zero prepared payload for the extra cache; ordinary powers still allocated |
| `certificate` | Isolated production certificate only; a failed result does not include a scan and is never reported as a continuation speedup |

The default variants are `disabled prepared`. Width is 246 on the five full-size states and 4 on small fixtures.
The existing composition budget is 8 GiB for full-size inputs. `prepared` gives the temporary cache the same cap;
this is a per-cache payload limit, not an RSS limit. Scalar mode preserves preparation on the existing giant-step cache.
Use `--scope next` to stop after the next group, or the default `--scope all` to return every remaining group.
`--minimum-excluded-degree` and `--maximum-attempts` configure the production policy; defaults here are 1 and 1.
The restored checkpoint is the first possible attempt boundary. The parity gate, natural bound or a larger threshold
can suppress the attempt.

The setup timer includes reading and validating the state, restoring parity and rebuilding the existing giant-step
cache; each stage also has a separate timer. Production obtains parity from its existing initial GCD, so this replay's
extra restoration GCD must not be added to a full-run projection as additional production work. The measured
interval includes all certificate planning, jump construction, temporary cache construction, arithmetic, fallback
scanning, and output-group copies. Serialized output validation happens afterward. Peak RSS is sampled immediately
after continuation and includes setup. There is no observer instrumentation in either variant. `plan.json` contains
planned composition counts for an isolated complete certificate, not observed counts when it exits early or discovery
later changes R.

The results retain source and binary hashes, compiler commands/features, thread settings, CPU affinity, input hashes,
plans, complete factor-group bytes, wall/CPU time, process peak RSS, restored parity, completed-prefix counts, attempt
counts, per-trial logs, host information and periodic resource/process-tree samples. The runner verifies state hashes
before each trial and afterward, and compares factor bytes and degree order across continuation variants. Trial order
alternates by case and repetition.
`comparison.json` keeps individual results and signed savings; a regression remains a negative saving.

These are saved-residual component timings. They exclude earlier discovery, initial Frobenius construction,
individual-factor norm recovery and fresh-counter search. They do not establish a total grinding speedup. Broad
algorithm/kernel attribution remains in the earlier experiment; the permanent instrumentation branch must be adapted
separately before profiling this integrated implementation. Do not apply a universal profiling-overhead correction.

## Prior experiment and current status

The earlier production-integration comparison, before parity integration, saved about 26.7–45.7 minutes on three
irreducible residuals of degrees 11,614, 14,028 and 14,429. Attempts on two composite residuals added about 4.4 and 7.2
minutes before ordinary discovery could resume. Those composite cases also belong to inputs with successful norm
recovery. Peak RSS increased from roughly 6.4 GiB to 11.6–11.8 GiB with the extra prepared cache. That was one paired
trial per input; the returned archive was `residual-integration-i_rhkzlk.tar.gz`, SHA-256
`669d6a90b68c0ce29d910d1cab2244b0e52170a4cfebfb5faca2b7d85ebad78b`.

The separate parity measurements already shorten the baseline on the three irreducible residuals. Do not add their
savings to the earlier certificate savings. The completed combined run measured the remaining certificate benefit
with parity enabled in both variants. Median continuation times from two paired rounds were:

| Residual factor degrees | Parity alone (min) | Parity + certificate (min) | Time reduction |
|---|---:|---:|---:|
| 5,890 + 7,461 | 38.72 | 38.83 | -0.28% |
| 11,614 | 18.75 | 6.83 | 63.58% |
| 14,028 | 31.78 | 15.50 | 51.23% |
| 5,288 + 5,402 | 20.90 | 20.88 | 0.05% |
| 14,429 | 35.55 | 8.97 | 74.77% |

Both composite controls made zero certificate attempts in both rounds and stayed near 6.4 GiB peak RSS. The larger
control was 0.20 seconds faster in one pair and 13.13 seconds slower in the other: its median 6.47-second slowdown
does not establish exact zero overhead. The smaller control differed by less than one second at the median. The three
irreducible residuals saved an additional 11.92, 16.28 and 26.58 minutes beyond parity, with peak RSS of 11.76–11.80 GiB.
These component percentages do not apply to a whole norm recovery. The benefit of an online grinding policy remains
unmeasured, and odd composite residuals can still incur failed-attempt costs.

All 20 full-size trials returned matching factor bytes and degree order. The run passed 114 C++ tests, nine runner
checks and 52 small-state trials. It used five single-threaded workers on an AMD EPYC 9454, GCC 15.2.0,
`-std=c++23 -O3 -DNDEBUG -march=native`, and BMI2/ADX/AVX512IFMA. No observation-overhead correction was applied.
The returned archive was `residual-combined-ysyp4jjg.tar.gz`, SHA-256
`37486eb24f4daf9e2913b2ea68b26e875ce4ad68ee7e8f9b45ae2f137b5b5e92`. It records the frozen uncommitted merge snapshot;
its source manifest, source patch, inputs and output hashes were checked against the distributed package.

Reference archive: `saved-residual-vo3g34sm.tar.gz`, SHA-256
`52ac775575974ea73e07cd93369e86905f3bbe84200811d0c5580734529fc717`.
The pinned reference package manifest has SHA-256
`de0ae157b1b4ca4b7b774b38775b4ba45f0b3abe2711454b73a6b5aa50eb8181`.
The package reuses that experiment's tested detachment, resource sampling and archive support. It includes pinned
small states and hashes of the full states, so it does not repeat their expensive construction.

## Build a package locally

From the repository root, a working tree, including an uncommitted merge, can be frozen without committing:

```sh
python3 libs/math/test/bench_test/run_residual_certificate_integration_benchmark.py package \
  --reference-package .test-logs/residual-certificate/packages/crypto3-saved-residual-x86-20261103 \
  --saved-run .test-logs/save-residual/analysis-vo3g34sm/raw-complete/saved-residual-vo3g34sm \
  --parity-run .test-logs/residual-parity/analysis-j4r19y2g/raw/residual-parity-j4r19y2g \
  --revision master --snapshot-working-tree \
  --output .test-logs/residual-certificate-integration/packages/crypto3-residual-combined-x86-20261104 \
  --checks .test-logs/residual-certificate-integration/checks
```

`package` exports the selected base revision, snapshots the complete production include/feature footprint, and verifies
the saved-state identities. With `--snapshot-working-tree` (also accepted as `--allow-dirty`), changed and new source
files are copied, deleted files are removed, and every frozen source is hashed. Provenance records the base revision,
HEAD, pending merge parent, status, source overlay and patch. This immutable snapshot is the measured source identity;
the base commit alone does not describe it. Without the flag, sources must match the selected revision. Packaging
performs no remote execution.

For a committed source revision, use `--revision HEAD` and omit `--snapshot-working-tree`. Choose an unused output
directory for a new distribution; existing pinned packages are never overwritten.

## Run on x86

The operator runs these commands on the benchmark machine after the tarball and checksum have been copied there:

```sh
cd ~
sha256sum -c crypto3-residual-combined-x86-20261104.tar.gz.sha256
tar -xzf crypto3-residual-combined-x86-20261104.tar.gz
cd crypto3-residual-combined-x86-20261104
./start.sh --workers 5 --repetitions 2 --timeout 18000 \
  --saved-run /home/riccardo/crypto3-saved-residual-x86-20261103/saved-residual-results/saved-residual-vo3g34sm
```

The launcher detaches with persistent logs and survives SSH disconnects; no tmux is needed. It prints `PID`, an exact
`FOLLOW_COMMAND`, a `CANCEL_COMMAND`, and the absolute `RESULTS_ARCHIVE` path. The supervisor reserves the last ten
minutes (or 20% for shorter timeouts) for cleanup and archiving. The final work, success, failure, cancellation or timeout
is archived automatically with a checksum. A hard timeout during archiving retains the source directory and prints
`PARTIAL_RESULTS_ARCHIVE`; it cannot promise a complete archive after SIGKILL or power loss.

By default, the runner finds the completed reference under
`~/crypto3-saved-residual-x86-20261103/saved-residual-results/`. If it was moved, add
`--saved-run /absolute/path/to/saved-residual-vo3g34sm`. Full states must match the pinned hashes. The run uses up to five
physical cores, one single-threaded process per case, with a 16 GiB admission allowance per worker plus a 16 GiB host
reserve. This admission estimate is not a hard memory limit. At least 15 GiB free storage is required for copied states
and archives.

The completed default full panel took 2 h 47 min 41 s through build, checks and arithmetic, excluding final archiving.
Allow roughly three hours on the same otherwise idle host; the timeout is five hours. There are
20 full trials: five cases, two variants and two repetitions. Including `scalar` or `certificate` adds trials; scalar
large-cache performance is not yet established. Select variants or cases explicitly to narrow the experiment:

```sh
./start.sh --variants disabled prepared scalar --scope next --workers 5 --timeout 21600
./start.sh --variants certificate --workers 5 --timeout 14400
```

`--cases` accepts case names from `panel.json`. Six arithmetic test suites, runner lifecycle checks and 52 small-state
trials always run before full trials. The small panel covers known odd and even parity, unknown parity, natural-bound
stopping, skipped and attempted certificates, scalar/prepared caches and isolated certificate failure.
`./start.sh --smoke --workers 1 --timeout 1800` builds and runs only those checks and is also supported locally on ARM.
It validates the package, not x86 performance.

## Retrieve the results

After completion, copy the automatically generated archive and checksum from the local repository:

```sh
mkdir -p .test-logs/residual-certificate-integration/returned
scp 'riccardo@100.64.144.53:/home/riccardo/crypto3-residual-combined-x86-20261104/residual-combined-results/*.tar.gz*' \
  .test-logs/residual-certificate-integration/returned/
```

This copies every run archive in that package's results directory. To retrieve just one, use the exact absolute
`RESULTS_ARCHIVE` printed at launch and completion, together with its `.sha256` sidecar.

The archive retains the complete frozen production sources, provenance, successful trials and partial work if another
case fails. Check `metadata.json` for the final outcome, `comparison.json` for paired results, per-case `commands.json`
for trial order and status, and `resources.jsonl` for memory pressure and worker accounting. A timeout or missing pair
is not a rejection or a speedup.
