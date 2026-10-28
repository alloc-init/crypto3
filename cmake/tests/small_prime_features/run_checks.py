#!/usr/bin/env python3
"""Check native CMake propagation and compile-only feature/fallback configurations."""
import argparse
import json
from pathlib import Path
import platform
import shlex
import subprocess
import tempfile
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log-root', type=Path)
    parser.add_argument('--cmake-arg', action='append', default=[])
    parser.add_argument('--x86-compile-only', action='store_true',
                        help='also use Apple Clang x86 objects on an ARM Mac; never run those objects')
    parser.add_argument('--run-small-prime-tests', action='store_true',
                        help='build and run the four focused suites in the native configuration')
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    repo = source.parents[2]
    parent = args.log_root or repo / '.test-logs/small-prime-cmake'
    parent.mkdir(parents=True, exist_ok=True)
    logs = Path(tempfile.mkdtemp(prefix='checks-', dir=parent)).resolve()
    commands = []
    results = {}
    print('LOG_DIRECTORY=' + str(logs), flush=True)

    def command(label, argv):
        print('START ' + label, flush=True)
        started = time.monotonic()
        with (logs / (label + '.stdout.log')).open('w') as stdout, (logs / (label + '.stderr.log')).open('w') as stderr:
            result = subprocess.run(argv, stdout=stdout, stderr=stderr, timeout=600)
        commands.append(dict(label=label, argv=argv, exit_status=result.returncode,
                             seconds=time.monotonic() - started))
        (logs / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
        if result.returncode:
            raise RuntimeError(f'{label} failed: see {logs / (label + ".stderr.log")}')
        print('PASS ' + label, flush=True)

    def check(label, build, options, expected=None):
        command(label + '-configure', ['cmake', '-S', str(source), '-B', str(build),
                                      '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',
                                      *args.cmake_arg, *options])
        selection = dict(line.split('=', 1) for line in (build / 'selection.txt').read_text().splitlines())
        actual = tuple(int(selection[name]) for name in ('avx2', 'avx512', 'crt_ifma'))
        if expected is not None and actual != expected:
            raise RuntimeError(f'{label}: expected {expected}, configured {actual}')
        command(label + '-build', ['cmake', '--build', str(build), '--target', 'feature_consumer', '--parallel', '2'])
        entries = json.loads((build / 'compile_commands.json').read_text())
        entry = next(item for item in entries if Path(item['file']).resolve() == source / 'feature_consumer.cpp')
        argv = entry.get('arguments') or shlex.split(entry['command'])
        # Explicit interface options must reach an external consumer even when the
        # configuring project's global flags already enabled those instructions.
        for enabled, required in ((actual[0], ['-mavx2', '-mfma']),
                                  (actual[1], ['-mavx512f', '-mavx512dq'])):
            if enabled and not all(flag in argv for flag in required):
                raise RuntimeError(f'{label}: missing transitive flags: {argv}')
        results[label] = {'selection': selection, 'compile_command': argv}
        (logs / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        return selection

    native = logs / 'native'
    native_options = ['-DFEATURE_CHECK_HOST=native', '-DFEATURE_CHECK_COMPILE_ONLY=OFF']
    native_selection = check('native-fresh', native, native_options)
    check('native-reconfigure', native, native_options)
    command('native-run', ['ctest', '--test-dir', str(native), '-R', '^small_prime_feature_consumer$',
                           '--output-on-failure'])

    x86 = platform.machine().lower() in ('x86_64', 'amd64')
    architecture = []
    if args.x86_compile_only:
        if platform.system() != 'Darwin' or x86:
            parser.error('--x86-compile-only is for an ARM Mac using Apple Clang')
        architecture = ['-DCMAKE_OSX_ARCHITECTURES=x86_64']
        x86 = True
    supported_compiler = (native_selection['compiler'] in ('GNU', 'Clang', 'AppleClang') and
                          native_selection['frontend'] != 'MSVC')
    if x86 and supported_compiler:
        compile_only = ['-DFEATURE_CHECK_COMPILE_ONLY=ON', *architecture]
        cases = [('full', (1, 1, 1)), ('avx2', (1, 0, 0)), ('no-fma', (0, 0, 0)),
                 ('no-dq', (1, 0, 0)), ('no-ifma', (1, 1, 0)), ('unsupported', (0, 0, 0))]
        for host, expected in cases:
            selection = check(host + '-fresh', logs / host, [*compile_only, '-DFEATURE_CHECK_HOST=' + host], expected)
            if host == 'no-dq' and selection['ifma_found'] != 'TRUE':
                raise RuntimeError('DQ detection incorrectly disabled independent IFMA support')

        # Reuse a directory with successful compiler probes, then remove a host
        # feature, reject instruction probes, and restore the original settings.
        reused = logs / 'full'
        check('no-dq-reconfigure', reused, [*compile_only, '-DFEATURE_CHECK_HOST=no-dq'], (1, 0, 0))
        check('compiler-rejects', reused, [*compile_only, '-DFEATURE_CHECK_HOST=full',
              '-DCMAKE_REQUIRED_FLAGS=-crypto3-invalid-feature-probe-option'], (0, 0, 0))
        check('compiler-restored', reused, [*compile_only, '-DFEATURE_CHECK_HOST=full',
              '-DCMAKE_REQUIRED_FLAGS='], (1, 1, 1))

        # Check the header guards themselves with one compiler feature removed.
        # Positive propagation above uses only the target's options; these
        # negative compile-only probes deliberately override the inherited flags.
        for name, override, expected in [('missing-fma-guard', '-mno-fma', (0, 0, 0)),
                                         ('missing-dq-guard', '-mno-avx512dq', (1, 0, 0))]:
            original = iter(results['compiler-restored']['compile_command'])
            argv = []
            definitions = dict(zip(('EXPECT_FFT_AVX2', 'EXPECT_FFT_AVX512', 'EXPECT_CRT_IFMA'), expected))
            for arg in original:
                if arg in ('-o', '-MF', '-MT', '-MQ'):
                    next(original)
                elif arg not in ('-c', '-MD', '-MMD'):
                    key = arg[2:].split('=', 1)[0] if arg.startswith('-D') else ''
                    argv.append('-D' + key + '=' + str(definitions[key]) if key in definitions else arg)
            command(name, [*argv, override, '-fsyntax-only'])

        toolchain = logs / 'cross.cmake'
        # Explicit system selection makes this a CMake cross configuration even
        # when the available compiler targets the same OS. No target code runs.
        toolchain.write_text(f'set(CMAKE_SYSTEM_NAME {platform.system()})\n'
                             'set(CMAKE_SYSTEM_PROCESSOR x86_64)\n')
        check('cross', logs / 'cross', [*compile_only, '-DFEATURE_CHECK_HOST=native',
              '-DCMAKE_TOOLCHAIN_FILE=' + str(toolchain),
              '-DHOST_CPU_AVX2_FOUND=TRUE', '-DHOST_CPU_FMA_FOUND=TRUE',
              '-DHOST_CPU_AVX512F_FOUND=TRUE', '-DHOST_CPU_AVX512DQ_FOUND=TRUE',
              '-DHOST_CPU_AVX512VL_FOUND=TRUE', '-DHOST_CPU_AVX512IFMA_FOUND=TRUE'], (0, 0, 0))
        if results['cross']['selection']['cross'] != 'TRUE':
            raise RuntimeError('cross-compilation check was not a cross configuration')

    if args.run_small_prime_tests:
        names = ['arithmetic', 'crt', 'transform', 'backend']
        targets = ['crypto3_math_small_prime_fft_' + name + '_test' for name in names]
        command('native-tests-build', ['cmake', '--build', str(native), '--target', *targets, '--parallel', '2'])
        pattern = '^crypto3_math_small_prime_fft_(' + '|'.join(names) + ')_test$'
        command('native-tests-list', ['ctest', '--test-dir', str(native), '-R', pattern, '--show-only=json-v1'])
        tests = json.loads((logs / 'native-tests-list.stdout.log').read_text())['tests']
        if {test['name'] for test in tests} != set(targets):
            raise RuntimeError('the focused tests were not registered')
        command('native-tests-run', ['ctest', '--test-dir', str(native), '-R', pattern, '--output-on-failure'])
    print('All requested checks passed. Simulated and cross targets were compile-only.', flush=True)


if __name__ == '__main__':
    main()
