#!/usr/bin/env python3
"""Focused comparison, detachment, cancellation, timeout and partial-archive checks."""

import argparse
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import tarfile
import tempfile
import time
import unittest
from unittest import mock

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
source = HERE / 'run_residual_certificate_integration_benchmark.py'
if not source.exists():
    source = HERE / 'run.py'
spec = importlib.util.spec_from_file_location('integration_runner', source)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
reference = None
if not (HERE / 'support/residual_support.py').exists():
    reference = HERE.parents[3] / '.test-logs/residual-certificate/packages/crypto3-saved-residual-x86-20261103'
runner.support = runner.load_support(reference)
support = runner.support


class RunnerChecks(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='certificate-runner-check-')
        self.root = Path(self.temporary.name)
        self.processes = []

    def tearDown(self):
        for process in self.processes:
            if process.poll() is None:
                process.send_signal(signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
        self.temporary.cleanup()

    def fake_package(self):
        package = self.root / 'package'
        package.mkdir()
        shutil.copyfile(source, package / 'run.py')
        (package / 'support').mkdir()
        shutil.copyfile(reference / 'run.py' if reference else HERE / 'support/residual_support.py',
                        package / 'support/residual_support.py')
        common = reference / 'support' / support.COMMON if reference else HERE / 'support' / support.COMMON
        shutil.copyfile(common, package / 'support' / support.COMMON)
        for name in (runner.SOURCE, 'README.md'):
            (package / name).write_text('fixture\n')
        # A minimal preflight program avoids recursively launching this test suite in lifecycle checks.
        (package / runner.TEST).write_text('raise SystemExit(0)\n')
        support.save(package / 'source-provenance.json', {})
        support.save(package / 'panel.json', [])
        for name in ('tests', 'states'):
            (package / name).mkdir()
        for name in runner.HEADERS:
            path = package / 'production' / runner.PREFIX / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('fixture\n')
        support.save(package / 'production/benchmark-snapshot.json', {})
        support.save(package / 'package.json', {'files': support.hashes(package)})
        return package

    def launch(self, behavior, timeout=15):
        package = self.fake_package()
        compiler = self.root / 'compiler.py'
        compiler.write_text(behavior)
        args = argparse.Namespace(package=package, results=self.root / 'results', timeout=timeout, workers=1,
            cxx=shlex.quote(sys.executable) + ' ' + shlex.quote(str(compiler)), repetitions=1,
            scope='next', minimum_excluded_degree=1, maximum_attempts=1,
            variants=['disabled', 'scalar'], smoke=True, saved_run=None, cases=None)
        original_popen = subprocess.Popen
        def capture(*args, **kwargs):
            process = original_popen(*args, **kwargs)
            self.processes.append(process)
            return process
        output = io.StringIO()
        with mock.patch.object(runner.subprocess, 'Popen', side_effect=capture), contextlib.redirect_stdout(output):
            self.assertEqual(runner.start(args), 0)
        printed = dict(line.split('=', 1) for line in output.getvalue().splitlines())
        self.assertIn('FOLLOW_COMMAND', printed)
        self.assertIn('RESULTS_ARCHIVE', printed)
        process = self.processes[-1]
        self.assertEqual(int(printed['PID']), process.pid)
        self.assertEqual(os.getsid(process.pid), process.pid)
        self.assertIn('--scope', process.args)
        self.assertIn('scalar', process.args)
        return process, Path(printed['RUN_DIR']), Path(printed['RESULTS_ARCHIVE'])

    def check_archive(self, process, directory, archive, expected_status, outcome):
        self.assertEqual(process.wait(timeout=20), expected_status)
        metadata = json.loads((directory / 'metadata.json').read_text())
        self.assertEqual(metadata['outcome'], outcome)
        self.assertTrue(archive.is_file())
        self.assertTrue(Path(str(archive) + '.sha256').is_file())
        self.assertIn(f'RESULTS_ARCHIVE={archive}', (directory / 'run.log').read_text())
        with tarfile.open(archive) as result:
            archived = json.load(result.extractfile(directory.name + '/metadata.json'))
            self.assertEqual(archived['outcome'], outcome)
            self.assertIn(directory.name + '/failure.log', result.getnames())

    def waiting_compiler(self):
        heartbeat = self.root / 'heartbeat'
        # Both the compiler and its child stay in the supervisor-owned process group.
        child = 'from pathlib import Path; import time\np=Path(' + repr(str(heartbeat)) + ')\n' \
                'while True:\n p.write_text(str(time.monotonic()))\n time.sleep(0.02)\n'
        compiler = 'import subprocess, sys, time\nsubprocess.Popen([sys.executable, "-c", ' + repr(child) + '])\ntime.sleep(60)\n'
        return compiler, heartbeat

    def wait_for(self, path):
        deadline = time.monotonic() + 8
        while not path.exists() and time.monotonic() < deadline:
            time.sleep(0.02)
        self.assertTrue(path.exists(), 'child did not start')

    def check_child_stopped(self, heartbeat):
        self.wait_for(heartbeat)
        before = heartbeat.read_text()
        time.sleep(0.15)
        self.assertEqual(before, heartbeat.read_text(), 'arithmetic descendant survived cancellation')

    def test_detached_job_ignores_hangup_and_archives_compiler_failure(self):
        process, directory, archive = self.launch('import time\ntime.sleep(0.5)\nraise SystemExit(3)\n')
        os.kill(process.pid, signal.SIGHUP)
        self.check_archive(process, directory, archive, 1, 'failed')
        self.assertTrue((directory / 'compiler.log').is_file())

    def test_cancellation_kills_descendants_and_archives_partial_results(self):
        compiler, heartbeat = self.waiting_compiler()
        process, directory, archive = self.launch(compiler)
        self.wait_for(heartbeat)
        os.kill(process.pid, signal.SIGHUP)
        os.kill(process.pid, signal.SIGTERM)
        self.check_archive(process, directory, archive, 130, 'canceled')
        self.check_child_stopped(heartbeat)

    def test_timeout_kills_descendants_and_archives_partial_results(self):
        compiler, heartbeat = self.waiting_compiler()
        process, directory, archive = self.launch(compiler, timeout=3)
        self.check_archive(process, directory, archive, 124, 'timeout')
        self.check_child_stopped(heartbeat)

    def test_mutated_state_is_rejected_before_execution(self):
        state = self.root / 'state'
        state.mkdir()
        (state / 'state.txt').write_text('original')
        manifest = self.root / 'state.json'
        support.save(manifest, support.hashes(state))
        (state / 'state.txt').write_text('changed')
        spec = self.root / 'spec.json'
        support.save(spec, {'directory': str(self.root / 'case'), 'cpu': None,
                           'state': str(state), 'state_manifest': str(manifest)})
        with mock.patch.object(runner.subprocess, 'Popen') as spawn:
            with self.assertRaisesRegex(ValueError, 'identity mismatch'):
                runner.run_case(argparse.Namespace(spec=spec))
            spawn.assert_not_called()

    def test_failed_certificate_penalty_is_retained_in_comparison(self):
        for trial in range(2):
            for variant, seconds in (('disabled', 10), ('prepared', 12)):
                directory = self.root / f'trial-{trial}-{variant}'
                directory.mkdir()
                (directory / 'group-0.bin').write_bytes(b'identical factor')
                support.save(directory / 'measurement.json', {
                    'mode': 'all', 'seconds': seconds, 'rng_draws': 0, 'maximum_attempts': 1,
                    'attempts': int(variant != 'disabled'), 'groups': [{'factor_degree': 3, 'group_degree': 3}]})
        comparison = runner.compare_case(self.root, ['disabled', 'prepared'], 2)
        self.assertEqual(comparison['median_seconds_saved']['prepared'], -2)
        self.assertEqual(len(comparison['pairs']), 2)
        (self.root / 'trial-1-prepared/group-0.bin').write_bytes(b'different factor')
        with self.assertRaisesRegex(ValueError, 'factor bytes'):
            runner.compare_case(self.root, ['disabled', 'prepared'], 2)

    def test_disabled_mode_must_not_attempt_a_certificate(self):
        directory = self.root / 'trial-0-disabled'
        directory.mkdir()
        (directory / 'group-0.bin').write_bytes(b'factor')
        support.save(directory / 'measurement.json', {'rng_draws': 0, 'maximum_attempts': 1, 'attempts': 1})
        with self.assertRaisesRegex(ValueError, 'disabled certificate'):
            runner.compare_case(self.root, ['disabled'], 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
