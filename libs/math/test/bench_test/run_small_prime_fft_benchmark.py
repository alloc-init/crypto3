#!/usr/bin/env python3
"""Package, verify and measure BN254 Fp12 polynomial products or coefficient combinations."""

import argparse
import csv
import datetime
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import statistics
import subprocess
import sys
import tarfile
import tempfile
import time
import traceback


BENCH = Path("libs/math/test/bench_test")
SOURCE = BENCH / "small_prime_fft_benchmark.cpp"
RUNNER = BENCH / "run_small_prime_fft_benchmark.py"
GUIDE = BENCH / "coefficient_combination_benchmark.md"
MANIFEST = "benchmark-snapshot.json"


def timestamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def digest(path):
    checksum = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def include_roots(repo):
    return sorted([*repo.glob("libs/*/include"), *repo.glob("libs/marshalling/*/include")])


def source_metadata(repo):
    snapshot = repo / MANIFEST
    if snapshot.is_file():
        metadata = json.loads(snapshot.read_text())
        for name, expected in metadata["files"].items():
            if digest(repo / name) != expected:
                raise RuntimeError(f"snapshot source changed: {name}")
        return metadata
    files = {repo / SOURCE, repo / RUNNER, repo / GUIDE, repo / BENCH / "CMakeLists.txt"}
    for root in include_roots(repo):
        files.update(path for path in root.rglob("*") if path.is_file())

    def git(*arguments):
        return subprocess.check_output(["git", "-C", str(repo), *arguments], text=True)

    if Path(git("rev-parse", "--show-toplevel").strip()).resolve() != repo:
        raise RuntimeError("--repo must be a Crypto3 checkout or an exported benchmark snapshot")
    return {
        "created_utc": timestamp(),
        "revision": git("rev-parse", "HEAD").strip(),
        "git_status": git("status", "--short", "--untracked-files=all"),
        "git_diff": git("diff", "--no-ext-diff", "--binary", "HEAD"),
        "files": {str(path.relative_to(repo)): digest(path) for path in sorted(files)},
    }


def package(repo, destination, coefficient_combination=False):
    metadata = source_metadata(repo)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation preserves previous handoffs. Every file is rooted under one directory.
    with destination.open("xb") as output, tarfile.open(fileobj=output, mode="w:gz") as archive:
        prefix = "crypto3-coefficient-combination-benchmark" if coefficient_combination else "crypto3-small-prime-benchmark"
        for name in metadata["files"]:
            archive.add(repo / name, arcname=f"{prefix}/{name}")
        if (repo / "LICENSE").is_file():
            archive.add(repo / "LICENSE", arcname=f"{prefix}/LICENSE")
        payload = (json.dumps(metadata, indent=2, sort_keys=True) + "\n").encode()
        entry = tarfile.TarInfo(f"{prefix}/{MANIFEST}")
        entry.size = len(payload)
        archive.addfile(entry, io.BytesIO(payload))
    # Do not produce a seemingly valid snapshot if a source changed during packaging.
    for name, expected in metadata["files"].items():
        if digest(repo / name) != expected:
            raise RuntimeError(f"source changed during packaging; discard {destination}: {name}")
    print(f"PACKAGE={destination}\nSHA256={digest(destination)}")


def summarize(samples, trials, warmups):
    with samples.open(newline="") as source:
        rows = list(csv.DictReader(source))
    results = {}
    common = {"multiply", "square", "multiply_low", "prepared_low", "divrem", "mulmod", "squaremod"}
    for backend in ("mixed_radix", "small_prime"):
        backend_rows = [row for row in rows if row["backend"] == backend]
        setup = [row for row in backend_rows if row["phase"] == "setup"]
        if {row["operation"] for row in setup} != {"backend", "divisor_context", "low_operand", "cyclic_divisor"}:
            raise RuntimeError(f"incomplete setup output for {backend}")
        cyclic = next(row for row in setup if row["operation"] == "cyclic_divisor")
        operations = common | ({"cyclic_remainder"} if int(cyclic["prepared_transform"]) else set())
        results[backend] = {"setup": setup, "operations": {}}
        for operation in sorted(operations):
            group = [row for row in backend_rows if row["phase"] == "trial" and row["operation"] == operation]
            if len(group) != trials or {int(row["trial"]) for row in group} != set(range(warmups, warmups + trials)):
                raise RuntimeError(f"incomplete samples: {backend}/{operation}")
            measurements = {}
            for metric in ("wall_ms", "cpu_ms"):
                values = [float(row[metric]) for row in group]
                if any(not math.isfinite(value) or value < 0 for value in values):
                    raise RuntimeError("invalid timing sample")
                measurements[metric] = {"median": statistics.median(values), "min": min(values), "max": max(values)}
            results[backend]["operations"][operation] = measurements
    speedups = {}
    for operation in sorted(results["mixed_radix"]["operations"].keys() & results["small_prime"]["operations"].keys()):
        baseline = results["mixed_radix"]["operations"][operation]["wall_ms"]["median"]
        candidate = results["small_prime"]["operations"][operation]["wall_ms"]["median"]
        speedups[operation] = baseline / candidate if candidate else None
        ratio = f"{speedups[operation]:.3f}x" if candidate else "below timer resolution"
        print(f"{operation:18} mixed_radix={baseline:.6f} ms  small_prime={candidate:.6f} ms  ratio={ratio}")
    return {"backends": results, "wall_time_speedups": speedups}


def summarize_combination(samples, trials, warmups):
    with samples.open(newline="") as source:
        rows = list(csv.DictReader(source))
    setup = [row for row in rows if row["phase"] == "setup"]
    if len(setup) != 4 or {row["operation"] for row in setup} != {
        "fixture_and_export", "backend", "divisor_context", "composition_cache"
    }:
        raise RuntimeError("incomplete coefficient-combination setup output")
    for phase, indices in (
        ("warmup", set(range(warmups))),
        ("trial", set(range(warmups, warmups + trials))),
        ("verification", set(range(warmups + trials))),
    ):
        group = [row for row in rows if row["phase"] == phase]
        if len(group) != len(indices) or {int(row["trial"]) for row in group} != indices:
            raise RuntimeError(f"incomplete coefficient-combination samples: {phase}")
        if any(row["operation"] != "coefficient_combination" for row in group):
            raise RuntimeError(f"unexpected coefficient-combination operation: {phase}")
    if len(rows) != 4 + 2 * (warmups + trials):
        raise RuntimeError("unexpected coefficient-combination rows")
    for row in rows:
        if any(not math.isfinite(float(row[key])) or float(row[key]) < 0 for key in ("wall_ms", "cpu_ms")):
            raise RuntimeError("invalid coefficient-combination timing sample")
    cache = next(row for row in setup if row["operation"] == "composition_cache")
    dimensions = {key: int(cache[key]) for key in (
        "degree", "outer_coefficients", "block_size", "block_count", "coefficient_products"
    )}
    for row in rows:
        if row["phase"] != "setup" and any(int(row[key]) != value for key, value in dimensions.items()):
            raise RuntimeError("coefficient-combination dimensions changed between samples")

    def statistics_for(phase):
        group = [row for row in rows if row["phase"] == phase and int(row["trial"]) >= warmups]
        result = {}
        for key in ("wall_ms", "cpu_ms"):
            values = [float(row[key]) for row in group]
            result[key] = {"median": statistics.median(values), "min": min(values), "max": max(values)}
        return result

    result = {
        "dimensions": dimensions, "setup": setup,
        "coefficient_combination": statistics_for("trial"),
        "verification": statistics_for("verification"),
        "memory": {key: max(int(row[key]) for row in rows) for key in (
            "process_peak_rss_bytes", "cache_coefficient_bytes", "output_and_reference_capacity_bytes"
        )},
        "verification_method": "all coefficients against reverse-order dot products; small fixtures also use schoolbook/Horner",
        "timing_scope": "sum of block-kernel calls; accumulator resets and verification excluded",
        "input_kind": "deterministic synthetic full-Fp12 inputs; real composition precomputation",
    }
    wall = result["coefficient_combination"]["wall_ms"]["median"] / 1000
    cpu = result["coefficient_combination"]["cpu_ms"]["median"] / 1000
    print(f"coefficient_combination: median wall={wall:.6f} s, CPU={cpu:.6f} s; "
          f"{dimensions['block_count']} blocks, {dimensions['coefficient_products']} Fp12 products")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[4])
    parser.add_argument("--package", type=Path, help="export source headers and this runner; do not build or measure")
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--boost-include", type=Path)
    parser.add_argument("--cpu", type=int, help="Linux CPU affinity; default is the highest permitted CPU")
    parser.add_argument("--coefficient-combination", action="store_true",
                        help="time the current scalar composition block kernel using real cached powers")
    parser.add_argument("--degree", type=int, default=14938)
    parser.add_argument("--trials", type=int, help="default: 3 for coefficient combinations, 5 for polynomial products")
    parser.add_argument("--warmups", type=int, help="default: 1 for coefficient combinations, 2 for polynomial products")
    parser.add_argument("--timeout", type=int, default=3600, help="limit for each command, in seconds")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.trials is None:
        args.trials = 3 if args.coefficient_combination else 5
    if args.warmups is None:
        args.warmups = 1 if args.coefficient_combination else 2
    repo = args.repo.resolve()
    if not (repo / SOURCE).is_file():
        parser.error("--repo does not contain the benchmark source")
    if not 1 <= args.degree <= 17487 or not 1 <= args.trials <= 1000 or not 0 <= args.warmups <= 1000 or args.timeout < 1:
        parser.error("degree must be 1..17487, trials 1..1000, warmups 0..1000, timeout positive")
    if args.package:
        package(repo, args.package.resolve(), args.coefficient_combination)
        return 0
    permitted = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    cpu = args.cpu if args.cpu is not None else (permitted[-1] if permitted else None)
    if cpu is not None and (cpu not in permitted or not shutil.which("taskset")):
        parser.error("CPU is outside the permitted affinity or taskset is unavailable")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    parent = repo / "build" / ("coefficient-combination-benchmarks" if args.coefficient_combination else "small-prime-fft-benchmarks")
    parent.mkdir(parents=True, exist_ok=True)
    logs = Path(tempfile.mkdtemp(prefix="run-", dir=parent))
    print(f"LOG_DIRECTORY={logs}\nCPU={cpu if cpu is not None else 'not pinned'}", flush=True)
    metadata = {
        "started_utc": timestamp(), "exit_status": 1, "platform": platform.platform(),
        "invocation": sys.argv, "configuration": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "cpu": cpu, "permitted_cpus": permitted,
        "environment": {key: env.get(key) for key in ("PATH", "CPATH", "CPLUS_INCLUDE_PATH", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "instrumentation": "none", "commands": [],
        "component": "coefficient_combination" if args.coefficient_combination else "polynomial_products",
    }
    save(logs / "metadata.json", metadata)

    def command(name, argv):
        print(f"{name}: {shlex.join(argv)}", flush=True)
        record = {"name": name, "argv": argv, "started_utc": timestamp(), "exit_status": None}
        metadata["commands"].append(record)
        save(logs / "metadata.json", metadata)
        started = time.monotonic()
        try:
            with (logs / f"{name}.stdout.log").open("w") as out, (logs / f"{name}.stderr.log").open("w") as err:
                completed = subprocess.run(argv, cwd=repo, env=env, stdout=out, stderr=err, timeout=args.timeout)
            record["exit_status"] = completed.returncode
            if completed.returncode:
                raise RuntimeError(f"{name} failed; see {logs / (name + '.stderr.log')}")
        finally:
            record["wall_seconds"] = time.monotonic() - started
            save(logs / "metadata.json", metadata)
        return (logs / f"{name}.stdout.log").read_text()

    try:
        sources = source_metadata(repo)
        save(logs / "sources.json", sources)
        for source in (SOURCE, RUNNER, GUIDE, BENCH / "CMakeLists.txt"):
            shutil.copyfile(repo / source, logs / source.name)
        command("machine", ["uname", "-a"])
        if sys.platform.startswith("linux"):
            command("cpu", ["lscpu"])
        elif sys.platform == "darwin":
            command("cpu", ["sysctl", "-n", "machdep.cpu.brand_string"])
        compiler = shlex.split(args.cxx)
        version = command("compiler", compiler + ["--version"])
        binary = logs / "small_prime_fft_benchmark"
        compile_command = compiler + ["-O3", "-DNDEBUG", "-std=c++23", "-march=native"]
        compile_command += ["-fconstexpr-steps=2147483647" if "clang" in version.lower() else "-fconstexpr-ops-limit=2147483647"]
        if args.boost_include:
            compile_command += ["-isystem", str(args.boost_include.resolve())]
        for root in include_roots(repo):
            compile_command += ["-I", str(root)]
        compile_command += [str(logs / SOURCE.name), "-o", str(binary)]
        command("compile", compile_command)
        metadata["binary_sha256"] = digest(binary)
        prefix = ["taskset", "-c", str(cpu)] if cpu is not None else []
        component = ["--coefficient-combination"] if args.coefficient_combination else []
        command("checks", prefix + [str(binary), *component, "--check", "--input-dir", str(logs / "check-inputs")])
        if not args.check_only:
            command("benchmark", prefix + [str(binary), *component, "--degree", str(args.degree), "--trials", str(args.trials),
                    "--warmups", str(args.warmups), "--input-dir", str(logs / "inputs")])
            shutil.copyfile(logs / "benchmark.stdout.log", logs / "samples.csv")
            summarize_run = summarize_combination if args.coefficient_combination else summarize
            save(logs / "summary.json", summarize_run(logs / "samples.csv", args.trials, args.warmups))
        for name, expected in sources["files"].items():
            if digest(repo / name) != expected:
                raise RuntimeError(f"source changed during the run: {name}")
        metadata["exit_status"] = 0
        print("Correctness checks passed." if args.check_only else "All measured outputs verified.", flush=True)
    except BaseException:
        (logs / "failure.log").write_text(traceback.format_exc())
        raise
    finally:
        metadata["input_sha256"] = {
            str(path.relative_to(logs)): digest(path)
            for directory in (logs / "inputs", logs / "check-inputs") if directory.exists()
            for path in sorted(directory.glob("*.txt"))
        }
        metadata["finished_utc"] = timestamp()
        save(logs / "metadata.json", metadata)
        archive_path = logs.with_suffix(".tar.gz")
        with tarfile.open(archive_path, "w:gz") as archive:
            for path in sorted(logs.iterdir()):
                if path.name != "small_prime_fft_benchmark":
                    archive.add(path, arcname=f"{logs.name}/{path.name}")
        print(f"EXIT_STATUS={metadata['exit_status']}\nRESULT_ARCHIVE={archive_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
