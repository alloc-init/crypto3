#!/usr/bin/env python3
"""Compare the production residual-discovery continuation with bounded certificates enabled and disabled."""

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import signal
import statistics
import subprocess
import sys
import tarfile
import tempfile
import time
import traceback

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
BENCH = Path('libs/math/test/bench_test')
SOURCE = 'residual_certificate_integration_benchmark.cpp'
RUNNER = 'run_residual_certificate_integration_benchmark.py'
TEST = 'test_residual_certificate_integration_runner.py'
GUIDE = 'residual_certificate_integration.md'
SUPPORT_SHA256 = 'ffcbd34365170a1d58f8423de5bf4554007b5059fe7623077b795f9d64ef3a04'
REFERENCE_MANIFEST_SHA256 = 'de0ae157b1b4ca4b7b774b38775b4ba45f0b3abe2711454b73a6b5aa50eb8181'
BASE = '5f2f2316a19d2e76c7347f9569cf33ffd5d27b5a'
PREFIX = Path('libs/math/include/nil/crypto3/math/polynomial')
HEADERS = ('backends/polynomial_backend.hpp', 'quotient_ring/polynomial_composition.hpp',
           'factorization/kaltofen_shoup_distinct_degree_factorization.hpp')
TESTS = ('distinct_degree_factorization', 'complete_factorization',
         'polynomial_x_norm_reconstruction', 'polynomial_composition')
SMALL = ('prime5', 'prime7', 'prime_power9', 'degree15', 'mixed5', 'equal6')
support = None


def load_support(reference=None):
    # Reuse the completed experiment's supervisor, hashing, resource sampling, source export and archive code.
    path = reference / 'run.py' if reference else HERE / 'support/residual_support.py'
    if hashlib.sha256(path.read_bytes()).hexdigest() != SUPPORT_SHA256:
        raise ValueError('benchmark support differs from the pinned completed experiment')
    spec = importlib.util.spec_from_file_location('residual_support', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def package(args):
    repo, destination = args.repo.resolve(), args.output.resolve()
    reference, saved = args.reference_package.resolve(), args.saved_run.resolve()
    if support.digest(reference / 'package.json') != REFERENCE_MANIFEST_SHA256:
        raise ValueError('reference package identity mismatch')
    support.verify_package(reference)
    revision = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', args.revision], text=True).strip()
    status = subprocess.check_output(['git', '-C', str(repo), 'status', '--short'], text=True)
    required = [str(BENCH / name) for name in (SOURCE, RUNNER, TEST, GUIDE)]
    required += [str(Path('libs/math/test') / (name + '.cpp')) for name in TESTS]
    if not args.allow_dirty:
        subprocess.run(['git', '-C', str(repo), 'ls-files', '--error-unmatch', *required],
                       stdout=subprocess.DEVNULL, check=True)
        if subprocess.check_output(['git', '-C', str(repo), 'diff', revision, '--',
                                    *required, 'libs/math/include'], text=True):
            raise ValueError('final packages must use committed source; use --allow-dirty only for local development')
    destination.mkdir(parents=True, exist_ok=False)
    support.export_sources(repo, revision, destination / 'production')
    if args.allow_dirty:
        for header in HEADERS:
            shutil.copyfile(repo / PREFIX / header, destination / 'production' / PREFIX / header)
        identity = support.hashes(destination / 'production')
        identity.pop('benchmark-snapshot.json')
        support.save(destination / 'production/benchmark-snapshot.json', {
            'revision': revision, 'files': identity, 'uncommitted_overlay': list(HEADERS),
            'git_diff': subprocess.check_output(['git', '-C', str(repo), 'diff', revision], text=True)})
    for name, target in ((SOURCE, SOURCE), (RUNNER, 'run.py'), (TEST, TEST), (GUIDE, 'README.md')):
        shutil.copyfile(repo / BENCH / name, destination / target)
    (destination / 'support').mkdir()
    shutil.copyfile(reference / 'run.py', destination / 'support/residual_support.py')
    shutil.copyfile(reference / 'support' / support.COMMON, destination / 'support' / support.COMMON)
    (destination / 'tests').mkdir()
    for name in TESTS:
        shutil.copyfile(repo / 'libs/math/test' / (name + '.cpp'), destination / 'tests' / (name + '.cpp'))
    metadata = json.loads((saved / 'metadata.json').read_text())
    if metadata.get('outcome') != 'complete' or metadata.get('smoke'):
        raise ValueError('source states must come from a completed full reference run')
    panel = json.loads((saved / 'panel.json').read_text())
    if panel != json.loads((reference / 'panel.json').read_text()):
        raise ValueError('source run and reference package panels differ')
    (destination / 'states').mkdir()
    for case in panel:
        identity = json.loads((saved / 'cases' / case['name'] / 'state-sha256.json').read_text())
        support.verify_hashes(saved / 'cases' / case['name'] / 'state', identity)
        support.save(destination / 'states' / (case['name'] + '.json'), identity)
    for name in SMALL:
        source = saved / 'checks' / (name + '-prepared0')
        identity = json.loads((source / 'state-sha256.json').read_text())
        support.verify_hashes(source / 'state', identity)
        target = destination / 'smoke' / name
        shutil.copytree(source / 'state', target / 'state')
        support.save(target / 'state-sha256.json', identity)
    support.save(destination / 'panel.json', panel)
    support.save(destination / 'source-provenance.json', {
        'created_utc': support.timestamp(), 'production_revision': revision, 'base_revision': BASE,
        'development_overlay': args.allow_dirty, 'git_status': status,
        'source_run': str(saved), 'source_run_metadata': metadata,
        'reference_package_manifest_sha256': REFERENCE_MANIFEST_SHA256,
        'scope': 'clean production continuation and optional isolated certificate; no observer build',
        'remote_execution': False})
    if args.checks:
        (destination / 'local-checks').mkdir()
        for path in sorted(args.checks.iterdir()):
            if path.suffix in ('.json', '.log'):
                shutil.copyfile(path, destination / 'local-checks' / path.name)
    (destination / 'start.sh').write_text('#!/bin/sh\nset -eu\npackage_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)\n'
                                         'exec python3 "$package_dir/run.py" start --package "$package_dir" "$@"\n')
    (destination / 'start.sh').chmod(0o755)
    support.save(destination / 'package.json', {'created_utc': support.timestamp(), 'production_revision': revision,
                                               'files': support.hashes(destination)})
    (destination / 'CONTENTS.sha256').write_text(''.join(f'{sha}  {name}\n' for name, sha in support.hashes(destination).items()))
    archive = destination.with_suffix('.tar.gz')
    if archive.exists():
        raise FileExistsError(archive)
    with tarfile.open(archive, 'w:gz', compresslevel=1) as output:
        output.add(destination, arcname=destination.name)
    Path(str(archive) + '.sha256').write_text(f'{support.digest(archive)}  {archive.name}\n')
    print(f'PACKAGE={archive}\nSHA256={support.digest(archive)}\nREMOTE_EXECUTION=false', flush=True)
    return 0


def build(package, directory, budget, cxx, smoke):
    common = support.common_module(package)
    compiler = shlex.split(cxx)
    version = budget.command('compiler', compiler + ['--version'])
    features = common.fp12_build_flags(package / 'production', directory, compiler,
                                       'portable' if smoke else 'auto', budget.command)
    flags = ['-std=c++23', '-O3', '-DNDEBUG', '-march=native', *features['flags']]
    if platform.system() == 'Darwin' and Path('/opt/homebrew/include/boost').exists():
        flags += ['-I/opt/homebrew/include']
    includes = [value for path in common.include_roots(package / 'production') for value in ('-I', str(path))]
    output = directory / 'build'
    output.mkdir(parents=True, exist_ok=True)
    commands = {}
    for name in (*TESTS, 'residual-benchmark'):
        source = package / SOURCE
        if name != 'residual-benchmark':
            source = output / (name + '-driver.cpp')
            source.write_text('#include "' + str(package / 'tests' / (name + '.cpp')) + '"\n'
                              '#include <boost/test/included/unit_test.hpp>\n')
        binary = output / name
        commands[name] = compiler + flags + includes + [str(source), '-o', str(binary)]
        support.save(directory / 'build.json', {'commands': commands, 'compiler': version,
                                               'flags': flags, 'features': features})
        budget.command('compile-' + name, commands[name])
        if name != 'residual-benchmark':
            budget.command('test-' + name, [str(binary), '--report_level=detailed', '--color_output=no'])
    support.save(directory / 'binary-sha256.json', {name: support.digest(output / name) for name in commands})
    return str(output / 'residual-benchmark')


def compare_case(directory, variants, repetitions):
    pairs, evidence = [], {}
    for trial in range(repetitions):
        reports = {}
        for variant in variants:
            folder = directory / f'trial-{trial}-{variant}'
            report = json.loads((folder / 'measurement.json').read_text())
            if report['rng_draws'] != 0 or report['attempts'] > report['maximum_attempts']:
                raise ValueError('RNG use or certificate attempt budget changed')
            reports[variant] = report
            if variant == 'certificate':
                continue
            hashes = {p.name: support.digest(p) for p in sorted(folder.glob('group-*.bin'))}
            if not hashes or (evidence and hashes != evidence):
                raise ValueError('continued discovery changed factor bytes')
            evidence = hashes
            if variant == 'disabled' and report['attempts'] != 0:
                raise ValueError('disabled certificate performed work')
        for variant in variants:
            if variant in ('disabled', 'certificate') or 'disabled' not in reports:
                continue
            old, new = reports['disabled'], reports[variant]
            if old['groups'] != new['groups'] or old['mode'] != new['mode']:
                raise ValueError('group order or continuation scope changed')
            pairs.append({'trial': trial, 'variant': variant, 'baseline': old, 'enabled': new,
                          'seconds_saved': old['seconds'] - new['seconds'],
                          'speedup': old['seconds'] / new['seconds']})
    result = {'pairs': pairs, 'group_sha256': evidence, 'repetitions': repetitions,
              'median_seconds_saved': {variant: statistics.median(p['seconds_saved'] for p in pairs
                                         if p['variant'] == variant)
                                      for variant in sorted({p['variant'] for p in pairs})},
              'interpretation': 'Failures include the failed certificate plus resumed discovery. No norm-recovery or grinding timing.'}
    support.save(directory / 'comparison.json', result)
    return result


def run_case(args):
    case = json.loads(args.spec.read_text())
    directory = Path(case['directory'])
    directory.mkdir(parents=True, exist_ok=True)
    if case['cpu'] is not None:
        os.sched_setaffinity(0, {case['cpu']})
    identity = json.loads(Path(case['state_manifest']).read_text())
    state = directory / 'state'
    support.verify_hashes(Path(case['state']), identity)
    shutil.copytree(case['state'], state)
    support.save(directory / 'state-sha256.json', identity)
    commands = []
    for trial in range(case['repetitions']):
        order = list(case['variants'])
        if (trial + case['order_index']) % 2:
            order.reverse()
        for variant in order:
            support.verify_hashes(state, identity)
            label = f'trial-{trial}-{variant}'
            prepared = case['prepared_bytes'] if variant in ('prepared', 'certificate') else 0
            argv = [case['binary'], '--state', str(state), '--output', str(directory / label),
                    '--mode', 'certificate' if variant == 'certificate' else case['scope'],
                    '--width', str(case['width']), '--prepared-bytes', str(case['prepared_bytes']),
                    '--minimum-excluded-degree', '0' if variant == 'disabled' else str(case['minimum_excluded']),
                    '--maximum-attempts', str(case['maximum_attempts']), '--certificate-prepared-bytes', str(prepared)]
            record = {'phase': label, 'argv': argv, 'started_utc': support.timestamp()}
            commands.append(record)
            support.save(directory / 'commands.json', commands)
            support.save(directory / 'progress.json', {'phase': label, 'started_utc': record['started_utc']})
            print(f'PHASE={label} CASE={case["name"]}', flush=True)
            with (directory / (label + '.log')).open('w') as output:
                # Stay in the worker process group so the supervisor owns the complete arithmetic process tree.
                with subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=output, stderr=subprocess.STDOUT,
                                      env=dict(os.environ, **support.THREADS)) as process:
                    record['pid'] = process.pid
                    support.save(directory / 'commands.json', commands)
                    status = process.wait()
            record.update(status=status, finished_utc=support.timestamp())
            support.save(directory / 'commands.json', commands)
            if status:
                raise RuntimeError(f'{label} failed: {status}')
    support.verify_hashes(state, identity)
    result = compare_case(directory, case['variants'], case['repetitions'])
    support.save(directory / 'progress.json', {'phase': 'complete', 'finished_utc': support.timestamp()})
    print(f'CASE_COMPLETE={case["name"]} SECONDS_SAVED={result["median_seconds_saved"]}', flush=True)
    return 0


def preflight(package, directory, budget, binary):
    for index, name in enumerate(SMALL):
        case = {'name': name, 'directory': str(directory / 'checks' / name), 'cpu': None, 'binary': binary,
                'state': str(package / 'smoke' / name / 'state'),
                'state_manifest': str(package / 'smoke' / name / 'state-sha256.json'),
                'variants': ['disabled', 'prepared', 'scalar', 'certificate'], 'scope': 'all',
                'width': 4, 'prepared_bytes': 1024 * 1024, 'minimum_excluded': 1, 'maximum_attempts': 1,
                'repetitions': 1, 'order_index': index}
        path = directory / (name + '-spec.json')
        support.save(path, case)
        budget.command('check-' + name, [sys.executable, str(package / 'run.py'), '_case', '--spec', str(path)])
        for variant in ('prepared', 'scalar'):
            report = json.loads((Path(case['directory']) / f'trial-0-{variant}/measurement.json').read_text())
            if report['attempts'] != 1:
                raise ValueError('small-fixture preflight did not execute the integrated certificate')
    support.save(directory / 'preflight.json', {'status': 'passed', 'cases': len(SMALL), 'trials': 4 * len(SMALL),
        'coverage': 'prime, prime-power, composite degree; failed full identity and proper-degree GCD; prepared and scalar; full continuation',
        'integrated_test_suites': TESTS})


def resolve_saved_run(args, package, panel):
    root = Path.home() / 'crypto3-saved-residual-x86-20261103/saved-residual-results'
    candidates = [args.saved_run.resolve()] if args.saved_run else sorted(root.glob('saved-residual-*'))
    for path in reversed(candidates):
        if path.is_dir() and all((path / 'cases' / case['name'] / 'state').is_dir() for case in panel):
            for case in panel:
                identity = json.loads((package / 'states' / (case['name'] + '.json')).read_text())
                support.verify_hashes(path / 'cases' / case['name'] / 'state', identity)
            return path
    raise ValueError('pass --saved-run /absolute/path/to/the/completed/saved-residual-run')


def run(args):
    package, directory = args.package.resolve(), args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    archive = directory.with_suffix('.tar.gz')
    print(f'RUN_DIR={directory}\nFOLLOW_COMMAND=tail -f {shlex.quote(str(directory / "run.log"))}\nRESULTS_ARCHIVE={archive}', flush=True)
    metadata = {'started_utc': support.timestamp(), 'outcome': 'failed', 'exit_status': 1, 'argv': sys.argv,
                'platform': platform.platform(), 'machine': platform.machine(), 'threads': support.THREADS,
                'timeout_seconds': args.timeout, 'smoke': args.smoke, 'repetitions': args.repetitions,
                'variants': args.variants, 'scope': args.scope, 'supervisor_pid': os.getpid(), 'results_archive': str(archive)}
    budget = support.Budget(args.timeout, directory)
    try:
        manifest = support.verify_package(package)
        support.save(directory / 'package.json', manifest)
        for name in (SOURCE, 'run.py', TEST, 'README.md', 'source-provenance.json', 'panel.json'):
            shutil.copyfile(package / name, directory / name)
        for name in ('support', 'tests', 'states'):
            shutil.copytree(package / name, directory / name)
        shutil.copyfile(package / 'production/benchmark-snapshot.json', directory / 'production-sources.json')
        for header in HEADERS:
            target = directory / 'production-headers' / header
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(package / 'production' / PREFIX / header, target)
        choices = support.selection(args.workers, args.smoke)
        metadata['selection'] = choices
        host = {'uname': list(os.uname()), 'python': sys.version}
        for name in ('cpuinfo', 'meminfo', 'version', 'self/cgroup'):
            path = Path('/proc') / name
            if path.exists():
                host[name] = path.read_text()
        support.save(directory / 'host.json', host)
        panel = json.loads((package / 'panel.json').read_text())
        if args.cases:
            selected = set(args.cases)
            panel = [case for case in panel if case['name'] in selected]
            if len(panel) != len(selected):
                raise ValueError('unknown or duplicate case selection')
        saved_run = None
        if not args.smoke:
            if shutil.disk_usage(directory).free < 15 * support.GIB:
                raise ValueError('need 15 GiB free for retained states and results archives')
            saved_run = resolve_saved_run(args, package, panel)
            metadata['saved_run'] = str(saved_run)
        support.save(directory / 'metadata.json', metadata)
        print(f'WORKERS={choices["workers"]} CPUS={choices["cpus"]} REPETITIONS={args.repetitions} VARIANTS={args.variants}', flush=True)
        budget.command('test-runner', [sys.executable, str(package / TEST)])
        binary = build(package, directory, budget, args.cxx, args.smoke)
        preflight(package, directory, budget, binary)
        if not args.smoke:
            jobs = []
            for index, original in enumerate(panel):
                case = dict(original, directory=str(directory / 'cases' / original['name']), binary=binary,
                            cpu=choices['cpus'][index % choices['workers']], repetitions=args.repetitions,
                            state=str(saved_run / 'cases' / original['name'] / 'state'),
                            state_manifest=str(package / 'states' / (original['name'] + '.json')),
                            variants=args.variants, scope=args.scope, order_index=index,
                            minimum_excluded=args.minimum_excluded_degree, maximum_attempts=args.maximum_attempts)
                path = directory / (case['name'] + '-spec.json')
                support.save(path, case)
                jobs.append((case['name'], [sys.executable, str(package / 'run.py'), '_case', '--spec', str(path)]))
            statuses = {}
            for offset in range(0, len(jobs), choices['workers']):
                statuses.update(budget.run_many(jobs[offset:offset + choices['workers']]))
            comparisons = {}
            for case in panel:
                path = directory / 'cases' / case['name'] / 'comparison.json'
                if path.exists():
                    comparisons[case['name']] = json.loads(path.read_text())
            support.save(directory / 'comparison.json', {'cases': comparisons, 'statuses': statuses,
                'scope': 'clean saved-residual continuations with failed attempts retained; no full recovery or fresh search timing'})
            if any(statuses.values()):
                raise RuntimeError('one or more cases failed; completed cases and partial evidence retained')
        support.verify_package(package)
        metadata.update(outcome='checks_only' if args.smoke else 'complete', exit_status=0)
    except BaseException as error:
        if isinstance(error, (TimeoutError, subprocess.TimeoutExpired)):
            metadata.update(outcome='timeout', exit_status=124)
        elif isinstance(error, (KeyboardInterrupt, InterruptedError)):
            metadata.update(outcome='canceled', exit_status=130)
        (directory / 'failure.log').write_text(traceback.format_exc())
        print(f'ERROR={error}', flush=True)
    finally:
        budget.finish_work()
        budget.stop()
        metadata.update(finished_utc=support.timestamp(), elapsed_seconds=time.monotonic() - budget.started)
        support.save(directory / 'metadata.json', metadata)
        print(f'ARCHIVING={archive} OUTCOME={metadata["outcome"]}', flush=True)
        try:
            support.archive_results(directory)
            print(f'RESULTS_ARCHIVE={archive}\nEXIT_STATUS={metadata["exit_status"]}', flush=True)
        finally:
            budget.close()
    return metadata['exit_status']


def start(args):
    package = args.package.resolve()
    support.verify_package(package)
    root = args.results.resolve()
    root.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='residual-integration-', dir=root))
    argv = [sys.executable, str(package / 'run.py'), 'run', '--package', str(package), '--directory', str(directory),
            '--timeout', str(args.timeout), '--workers', str(args.workers), '--cxx', args.cxx,
            '--repetitions', str(args.repetitions), '--scope', args.scope,
            '--minimum-excluded-degree', str(args.minimum_excluded_degree), '--maximum-attempts', str(args.maximum_attempts),
            '--variants', *args.variants]
    if args.smoke:
        argv.append('--smoke')
    if args.saved_run:
        argv += ['--saved-run', str(args.saved_run.resolve())]
    if args.cases:
        argv += ['--cases', *args.cases]
    with (directory / 'run.log').open('w') as output:
        previous = signal.signal(signal.SIGHUP, signal.SIG_IGN)
        try:
            process = subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=output, stderr=subprocess.STDOUT,
                                       start_new_session=True, close_fds=True)
        finally:
            signal.signal(signal.SIGHUP, previous)
    (directory / 'supervisor.pid').write_text(str(process.pid) + '\n')
    print(f'PID={process.pid}\nRUN_DIR={directory}\nFOLLOW_COMMAND=tail -f {shlex.quote(str(directory / "run.log"))}\n'
          f'CANCEL_COMMAND=kill -TERM {process.pid}\nRESULTS_ARCHIVE={directory.with_suffix(".tar.gz")}', flush=True)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('package')
    p.add_argument('--repo', type=Path, default=Path.cwd())
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--reference-package', type=Path, required=True)
    p.add_argument('--saved-run', type=Path, required=True)
    p.add_argument('--revision', default='HEAD')
    p.add_argument('--allow-dirty', action='store_true')
    p.add_argument('--checks', type=Path)
    for name in ('start', 'run'):
        p = sub.add_parser(name)
        p.add_argument('--package', type=Path, default=HERE)
        p.add_argument('--timeout', type=float, default=14400)
        p.add_argument('--workers', type=int, choices=range(1, 6), default=5)
        p.add_argument('--repetitions', type=int, choices=range(1, 6), default=1)
        p.add_argument('--cxx', default=os.environ.get('CXX', 'c++'))
        p.add_argument('--smoke', action='store_true')
        p.add_argument('--saved-run', type=Path)
        p.add_argument('--scope', choices=('next', 'all'), default='all')
        p.add_argument('--variants', nargs='+', choices=('disabled', 'prepared', 'scalar', 'certificate'),
                       default=['disabled', 'prepared'])
        p.add_argument('--cases', nargs='+')
        p.add_argument('--minimum-excluded-degree', type=int, default=1)
        p.add_argument('--maximum-attempts', type=int, default=1)
        if name == 'start':
            p.add_argument('--results', type=Path, default=Path.cwd() / 'residual-integration-results')
        else:
            p.add_argument('--directory', type=Path, required=True)
    p = sub.add_parser('_case', help=argparse.SUPPRESS)
    p.add_argument('--spec', type=Path, required=True)
    args = parser.parse_args(argv)
    if hasattr(args, 'timeout') and not 1 <= args.timeout <= 7 * 24 * 3600:
        parser.error('timeout must be between one second and seven days')
    if hasattr(args, 'variants') and (len(set(args.variants)) != len(args.variants) or
            args.minimum_excluded_degree < 1 or args.maximum_attempts < 0):
        parser.error('variants must be unique, minimum excluded degree positive and attempt count nonnegative')
    global support
    support = load_support(args.reference_package if args.command == 'package' else None)
    return {'package': package, 'start': start, 'run': run, '_case': run_case}[args.command](args)


if __name__ == '__main__':
    sys.exit(main())
