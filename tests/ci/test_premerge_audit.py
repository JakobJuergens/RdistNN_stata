"""Behavioral checks for Stata verification wrappers; never launches Stata.

The PowerShell process launcher is replaced in a child PowerShell session.
The real wrappers still create their do-files, manage logs, and decide success.
"""

import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
PWSH = shutil.which("pwsh")

PROCESS_PROBE = r'''
Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# No executable is launched: model only process result and output log.
function global:Start-Process {
    param($FilePath, $ArgumentList, $WorkingDirectory, [switch]$Wait,
          [switch]$PassThru, $WindowStyle)
    if ($env:HARNESS_WRITE_LOG -eq "1") {
        Set-Content -LiteralPath (Join-Path $WorkingDirectory $env:HARNESS_LOG_NAME) `
            -Value $env:HARNESS_LOG_CONTENT
    }
    return [pscustomobject]@{ ExitCode = [int]$env:HARNESS_EXIT_CODE }
}

$arguments = @{ StataExe = $env:HARNESS_FAKE_EXE; RepoRoot = $env:HARNESS_ROOT }
if ($env:HARNESS_KEEP_LOG -eq "1") { $arguments.KeepLog = $true }
& $env:HARNESS_ENTRY_POINT @arguments
'''

PARSER_PROBE = r'''
Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
. $env:HARNESS_COMMON
Assert-StataSuccessLog -LogPath $env:HARNESS_LOG
'''


@unittest.skipUnless(PWSH, "PowerShell 7 is required for mocked wrapper checks")
class StataVerificationHarnessTests(unittest.TestCase):
    wrappers = {
        "runtime": ("stata_runtime_check.ps1", "tests.log", "All 123 tests passed."),
        "install": ("stata_install_check.ps1", "install_check.log", "INSTALL_CHECK_PASSED"),
        "documentation": (
            "stata_documentation_check.ps1", "documentation_examples.log",
            "DOCUMENTATION_EXAMPLES_PASSED 6",
        ),
    }

    def run_wrapper(self, name, *, exit_code=0, write_log=True, keep_log=False):
        entry, log_name, marker = self.wrappers[name]
        with tempfile.TemporaryDirectory(prefix="rdistnn-harness-") as directory:
            sandbox = Path(directory)
            (sandbox / "tests").mkdir()
            (sandbox / "tests/tests.do").write_text("// never executed\n")
            probe = sandbox / "probe.ps1"
            probe.write_text(PROCESS_PROBE)
            if keep_log:
                (sandbox / log_name).write_text("All 456 tests passed.\n")
            env = os.environ.copy()
            env.update({
                "HARNESS_ENTRY_POINT": str(ROOT / "tools/ci" / entry),
                "HARNESS_FAKE_EXE": sys.executable,
                "HARNESS_ROOT": str(sandbox),
                "HARNESS_LOG_NAME": log_name,
                "HARNESS_LOG_CONTENT": marker,
                "HARNESS_EXIT_CODE": str(exit_code),
                "HARNESS_WRITE_LOG": "1" if write_log else "0",
                "HARNESS_KEEP_LOG": "1" if keep_log else "0",
                "TMPDIR": directory,
                "TEMP": directory,
                "TMP": directory,
            })
            result = subprocess.run(
                [PWSH, "-NoProfile", "-File", str(probe)], env=env,
                capture_output=True, text=True, timeout=30,
            )
            previous_log = sandbox / "tests.previous.log"
            current_log = sandbox / log_name
            result.previous_log_content = previous_log.read_text() if previous_log.exists() else None
            result.current_log_content = current_log.read_text() if current_log.exists() else None
            return result

    def run_parser(self, content):
        with tempfile.TemporaryDirectory(prefix="rdistnn-log-parser-") as directory:
            sandbox = Path(directory)
            log = sandbox / "tests.log"
            log.write_text(content)
            probe = sandbox / "parse.ps1"
            probe.write_text(PARSER_PROBE)
            env = os.environ.copy()
            env.update({
                "HARNESS_COMMON": str(ROOT / "tools/ci/stata_common.ps1"),
                "HARNESS_LOG": str(log),
            })
            return subprocess.run(
                [PWSH, "-NoProfile", "-File", str(probe)], env=env,
                capture_output=True, text=True, timeout=30,
            )

    def test_wrappers_accept_fresh_successful_process(self):
        for name in self.wrappers:
            with self.subTest(wrapper=name):
                result = self.run_wrapper(name)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_wrappers_reject_nonzero_process_even_with_fresh_success_marker(self):
        for name in self.wrappers:
            with self.subTest(wrapper=name):
                result = self.run_wrapper(name, exit_code=2)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_wrappers_reject_process_that_writes_no_completion_log(self):
        for name in self.wrappers:
            with self.subTest(wrapper=name):
                result = self.run_wrapper(name, write_log=False)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_runtime_keep_log_does_not_accept_previous_run(self):
        result = self.run_wrapper("runtime", keep_log=True, write_log=False)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_runtime_keep_log_preserves_previous_log_and_checks_new_output(self):
        result = self.run_wrapper("runtime", keep_log=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.previous_log_content, "All 456 tests passed.\n")
        self.assertEqual(result.current_log_content, "All 123 tests passed.\n")
        self.assertIn("All 123 tests passed.", result.stdout)

    def test_runtime_parser_accepts_rendered_success_line(self):
        result = self.run_parser("All 123 tests passed.\n")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stdout.strip(), "123")

    def test_runtime_parser_rejects_echoed_command_without_success_output(self):
        result = self.run_parser('. di as result "All 123 tests passed."\nr(198);\n')
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_runtime_parser_rejects_zero_assertions(self):
        result = self.run_parser("All 0 tests passed.\n")
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_runtime_parser_rejects_failure_and_unexpected_abort(self):
        for failure in (
            "1 of 124 tests FAILED.",
            "  [runner] test-example.do aborted unexpectedly (rc=198)",
        ):
            with self.subTest(failure=failure):
                result = self.run_parser("All 123 tests passed.\n" + failure + "\n")
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)


@unittest.skipUnless(shutil.which("bash") and shutil.which("git"), "bash and git are required")
class StataVersionFloorTests(unittest.TestCase):
    def test_package_gate_rejects_old_integer_and_decimal_versions(self):
        # Copy only the static check's required package inputs, not a checkout
        # or build tree. An empty temporary Git index supports hygiene checks.
        manifest = (ROOT / "rdistnn.pkg").read_text()
        inputs = set(re.findall(r"^[fF]\s+(\S+)", manifest, re.M))
        inputs.update({
            ".gitignore", "README.md", "src/plugin/Makefile", "src/plugin/README.md",
            "tools/ci/pkg_check.sh", "tools/ci/plugin_build_check.ps1",
        })
        inputs.update(str(p.relative_to(ROOT)) for p in (ROOT / ".github/workflows").glob("*.yml"))
        with tempfile.TemporaryDirectory(prefix="rdistnn-version-gate-") as directory:
            sandbox = Path(directory)
            for name in inputs:
                destination = sandbox / name
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(ROOT / name, destination)
            (sandbox / "tests").mkdir()
            (sandbox / "tests/version-probe.do").write_text("version 17\n")
            subprocess.run(
                ["git", "init", "--quiet", str(sandbox)], check=True,
                capture_output=True, text=True,
            )
            target = sandbox / "src/np_dnn.ado"
            original = target.read_text()
            for version, should_pass in (("17", True), ("15", False), ("16.1", False), ("17.1", True)):
                with self.subTest(version=version):
                    mutated, count = re.subn(
                        r"(?m)^(\s*version\s+)[0-9.]+", lambda m: m.group(1) + version,
                        original,
                    )
                    self.assertGreater(count, 0, "fixture has no version declaration")
                    target.write_text(mutated)
                    result = subprocess.run(
                        ["bash", "tools/ci/pkg_check.sh"], cwd=sandbox,
                        capture_output=True, text=True, timeout=30,
                    )
                    if should_pass:
                        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    else:
                        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                        self.assertIn("  OLD      src/np_dnn.ado", result.stdout)


if __name__ == "__main__":
    unittest.main()
