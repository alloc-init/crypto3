# Production residual-certificate benchmark

This package measures the actual library continuation with certificates disabled and enabled. Both variants use the
same pinned clean production build. The default comparison covers **all remaining distinct-degree groups** from each
saved checkpoint. A composite residual pays for its failed certificate **and** its complete fallback scan.

The algorithm and caller controls are documented in [polynomial factorization](../../docs/polynomial_factorization.md).
The production default remains disabled. The package enables one attempt at the restored boundary; these states were
selected for the earlier experiment and do not calibrate a trigger for fresh grinding. No factor oracle is passed to
the production trigger. Captured factors are used only to validate inputs and outputs outside the measured interval.

## What this measures

| Variant | Measured work |
|---|---|
| `disabled` | Production continuation, with no certificate attempts |
| `prepared` | Production continuation, one certificate attempt with an additional prepared-cache allowance |
| `scalar` | The same continuation with zero prepared payload for the extra cache; ordinary powers still allocated |
| `certificate` | Isolated production certificate only; a failed result does not include a scan and is never reported as a continuation speedup |

The default variants are `disabled prepared`. Width is 246 on the five full-size states and 4 on small fixtures.
The existing composition budget is 8 GiB for full-size inputs. `prepared` gives the temporary cache the same cap;
this is a per-cache payload limit, not an RSS limit. Scalar mode preserves preparation on the existing giant-step cache.
Use `--scope next` to stop after the next group, or the default `--scope all` to return every remaining group.
`--minimum-excluded-degree` and `--maximum-attempts` configure the production policy; defaults here are 1 and 1.
The first eligible point is the restored checkpoint. Selecting a much larger threshold can suppress attempts.

The setup timer includes reading and validating the state and rebuilding the existing giant-step cache. The measured
interval includes all certificate planning, jump construction, temporary cache construction, arithmetic, fallback
scanning, and output-group copies. Serialized output validation happens afterward. Peak RSS covers the entire process,
including setup. There is no observer instrumentation in either variant. `plan.json` contains planned composition
counts for an isolated complete certificate, not observed counts when it exits early or discovery later changes R.

The results retain source and binary hashes, compiler commands/features, thread settings, CPU affinity, input hashes,
plans, complete factor-group bytes, wall/CPU time, process peak RSS, attempt counts, per-trial logs, host information,
and periodic resource/process-tree samples. The runner verifies state hashes before each trial and afterward, and
compares factor bytes and degree order across continuation variants. Trial order alternates by case and repetition.
`comparison.json` keeps individual results and signed savings; a regression remains a negative saving.

These are saved-residual component timings. They exclude earlier discovery, initial Frobenius construction,
individual-factor norm recovery and fresh-counter search. They do not establish a total grinding speedup. Broad
algorithm/kernel attribution remains in the earlier experiment; the permanent instrumentation branch must be adapted
separately before profiling this integrated implementation. Do not apply a universal profiling-overhead correction.

## Prior experiment and current status

The earlier clean prototype saved about 26.8–46.1 minutes on three irreducible residuals of degrees 11,614, 14,028 and
14,429. Attempts on two composite residuals added about 4.4 and 7.4 minutes before ordinary discovery could resume.
Those composite cases also belong to inputs with successful norm recovery. Peak RSS increased from roughly 6.4 GiB
to 11.6–11.8 GiB with the extra prepared cache. These figures describe the completed prototype run; integrated x86
latency and any online policy benefit remain unmeasured until this package is run.

Reference archive: `saved-residual-vo3g34sm.tar.gz`, SHA-256
`52ac775575974ea73e07cd93369e86905f3bbe84200811d0c5580734529fc717`.
The pinned reference package manifest has SHA-256
`de0ae157b1b4ca4b7b774b38775b4ba45f0b3abe2711454b73a6b5aa50eb8181`.
The package reuses that experiment's tested detachment, resource sampling and archive support. It includes pinned
small states and hashes of the full states, so it does not repeat their expensive construction.

## Build a package locally

From the repository root, after committing the production change:

```sh
python3 libs/math/test/bench_test/run_residual_certificate_integration_benchmark.py package \
  --reference-package .test-logs/residual-certificate/packages/crypto3-saved-residual-x86-20261103 \
  --saved-run .test-logs/save-residual/analysis-vo3g34sm/raw-complete/saved-residual-vo3g34sm \
  --output .test-logs/residual-certificate-integration/packages/crypto3-residual-certificate-integration-x86-20261104 \
  --checks .test-logs/residual-certificate-integration/checks
```

`package` exports the exact commit's headers and verifies the saved-state identities. It performs no remote execution.
For local development only, `--allow-dirty` records the three changed headers as an explicit source overlay; do not use
such a package as the final pinned distribution. A normal package requires committed benchmark sources and tests.

## Run on x86

The operator runs these commands on the benchmark machine after the tarball and checksum have been copied there:

```sh
cd ~
sha256sum -c crypto3-residual-certificate-integration-x86-20261104.tar.gz.sha256
tar -xzf crypto3-residual-certificate-integration-x86-20261104.tar.gz
cd crypto3-residual-certificate-integration-x86-20261104
./start.sh --workers 5 --timeout 14400
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

The default full panel is expected to take roughly two hours plus build, validation and archiving, based on the prior
component timings. This is a planning estimate, not a new measurement; the timeout is four hours. `--repetitions 2`
doubles the arithmetic work. Including `scalar` or `certificate` adds trials; scalar large-cache performance is not yet
established. Select variants or cases explicitly to narrow the experiment, for example:

```sh
./start.sh --variants disabled prepared scalar --scope next --workers 5 --timeout 21600
./start.sh --variants certificate --workers 5 --timeout 14400
```

`--cases` accepts case names from `panel.json`. Small correctness and lifecycle checks always run before full trials.
`./start.sh --smoke --workers 1 --timeout 1800` builds and runs only those checks and is also supported locally on ARM.
It validates the package, not x86 performance.

## Retrieve the results

After completion, use the absolute path printed as `RESULTS_ARCHIVE` in the persistent log. From the local repository,
replace the sample run suffix with the printed one:

```sh
mkdir -p .test-logs/residual-certificate-integration/returned
scp user@benchmark-host:/absolute/path/to/residual-integration-RUN.tar.gz \
  .test-logs/residual-certificate-integration/returned/
scp user@benchmark-host:/absolute/path/to/residual-integration-RUN.tar.gz.sha256 \
  .test-logs/residual-certificate-integration/returned/
```

The archive retains successful trials and partial work if another case fails. Check `metadata.json` for the final
outcome, `comparison.json` for paired results, per-case `commands.json` for trial order and status, and `resources.jsonl`
for memory pressure and worker accounting. A timeout or missing pair is not a rejection or a speedup.
