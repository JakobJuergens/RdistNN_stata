"""Exercise real Ruff diagnostics in disposable independent checkouts."""

from contextlib import redirect_stderr
import importlib.util
import io
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]


class ScriptToolingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="rdistnn-script-tools-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "tools/ci").mkdir(parents=True)
        (self.root / "tests").mkdir()
        for source in (".ruff.toml", "tools/ci/check_scripts.py"):
            shutil.copy2(ROOT / source, self.root / source)
        self.runner = self.root / "tools/ci/check_scripts.py"

    def run_check(self, *args):
        return subprocess.run([sys.executable, str(self.runner), *args],
                              cwd=self.root.parent, capture_output=True, text=True)

    def fixture(self, name, content):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_real_lint_and_selected_formatting_outside_checkout(self):
        selected = self.fixture("tools/selected.py", "value= 1\nprint(value)\n")
        other = self.fixture("tools/other.py", "other= 2\nprint(other)\n")
        self.assertEqual(self.run_check().returncode, 0)
        result = self.run_check("--format-check", "tools/selected.py")
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.run_check("--format", "tools/selected.py").returncode, 0)
        self.assertEqual(selected.read_text(), "value = 1\nprint(value)\n")
        self.assertEqual(other.read_text(), "other= 2\nprint(other)\n")
        self.assertNotEqual(self.run_check("--format").returncode, 0)
        selected.write_text("import math\n", encoding="utf-8")
        result = self.run_check("tools/selected.py")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("F401", result.stdout + result.stderr)
        self.assertFalse((self.root / ".ruff_cache").exists())

    def test_excluded_sources_and_empty_selection_fail(self):
        for name in ("tools/vendor/broken.py", "tools/generated/broken.py",
                     "tools/build/broken.py", "tools/archives/broken.py",
                     "inst/paper-replication/affine/source/tools/broken.py"):
            path = self.fixture(name, "invalid python (\n")
            # Direct Ruff invocation must also respect force-exclude.
            result = subprocess.run(["ruff", "check", "--no-cache", "--config",
                                     str(self.root / ".ruff.toml"), str(path)],
                                    cwd=self.root.parent, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            if name.startswith("tools/"):
                result = self.run_check(name)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("zero eligible", result.stderr)
        self.assertEqual(self.run_check().returncode, 0)
        empty = self.root / "tests/empty"
        empty.mkdir()
        self.assertNotEqual(self.run_check("tests/empty").returncode, 0)
        self.assertNotEqual(self.run_check("tools/missing.py").returncode, 0)
        self.fixture("tools/not_python.ps1", "Write-Output 'valid'\n")
        self.assertNotEqual(self.run_check("tools/not_python.ps1").returncode, 0)
        self.assertNotEqual(self.run_check("../outside.py").returncode, 0)

    def test_format_directory_selection_cannot_rewrite_files(self):
        selected = self.fixture("tools/selected.py", "value= 1\nprint(value)\n")
        other = self.fixture("tests/other.py", "other= 2\nprint(other)\n")
        original = {path: path.read_bytes() for path in (selected, other)}
        for option in ("--format", "--format-check"):
            for directory in ("tools", "tests"):
                result = self.run_check(option, "tools/selected.py", directory)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("not directories", result.stderr)
                for path, content in original.items():
                    self.assertEqual(path.read_bytes(), content)

    def test_missing_and_wrong_ruff_fail_before_analysis(self):
        spec = importlib.util.spec_from_file_location("script_tooling_runner", self.runner)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with patch.object(module.shutil, "which", return_value=None), redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as failure:
                module.main([])
            self.assertEqual(failure.exception.code, 2)
        wrong = subprocess.CompletedProcess([], 0, "ruff 0.0.1\n", "")
        with patch.object(module.shutil, "which", return_value="ruff"):
            with patch.object(module.subprocess, "run", return_value=wrong), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as failure:
                    module.main([])
                self.assertEqual(failure.exception.code, 2)


class PowerShellToolingTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(shutil.which("pwsh"), "Install PowerShell 7 to run analyzer tests")
        self.temp = tempfile.TemporaryDirectory(prefix="rdistnn-powershell-tools-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "tools/ci").mkdir(parents=True)
        for name in ("check_powershell.ps1", "PSScriptAnalyzerSettings.psd1"):
            shutil.copy2(ROOT / "tools/ci" / name, self.root / "tools/ci" / name)
        self.runner = self.root / "tools/ci/check_powershell.ps1"
        self.fixture = self.root / "tools/fixtures/sample.ps1"
        self.fixture.parent.mkdir()

    def analyze(self, content, *options):
        self.fixture.write_text(content, encoding="utf-8")
        return subprocess.run(["pwsh", "-NoProfile", "-File", str(self.runner),
                               "-Paths", "tools/fixtures", *options],
                              cwd=self.root.parent, capture_output=True, text=True)

    def test_selected_directory_links_cannot_escape_or_bypass_exclusions(self):
        outside_temp = tempfile.TemporaryDirectory(prefix="rdistnn-outside-ps-")
        self.addCleanup(outside_temp.cleanup)
        outside = Path(outside_temp.name)
        (outside / "nested").mkdir()
        (outside / "nested/broken.ps1").write_text("function Broken {\n", encoding="utf-8")
        excluded = self.root / "build/excluded"
        excluded.mkdir(parents=True)
        (excluded / "broken.ps1").write_text("function Broken {\n", encoding="utf-8")
        (self.root / "tools/outside-link").symlink_to(outside, target_is_directory=True)
        (self.root / "tools/excluded-link").symlink_to(excluded, target_is_directory=True)
        for selected in ("tools/outside-link", "tools/outside-link/nested", "tools/excluded-link"):
            result = subprocess.run(["pwsh", "-NoProfile", "-File", str(self.runner),
                                     "-Paths", selected], cwd=self.root.parent,
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("must not traverse symbolic links or reparse points", result.stdout + result.stderr)
            self.assertNotIn("MissingEndCurlyBrace", result.stdout + result.stderr)

    def test_real_warning_parse_and_compatibility_diagnostics(self):
        result = self.analyze("Write-Output 'valid'\n")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        result = self.analyze("$unused = 1\nWrite-Output 'warning'\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PSUseDeclaredVarsMoreThanAssignments", result.stdout + result.stderr)
        result = self.analyze("function Broken {\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("MissingEndCurlyBrace", result.stdout + result.stderr)
        compatibility = ("Test-Json -Json '{}'\n"
                         "[System.IO.Path]::GetRelativePath('a', 'b')\n"
                         "Write-Output ($true ? 'yes' : 'no')\n")
        result = self.analyze(compatibility)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PSUseCompatibleSyntax", result.stdout + result.stderr)
        self.assertIn("PSUseCompatibleCommands", result.stdout + result.stderr)
        self.assertIn("PSUseCompatibleTypes", result.stdout + result.stderr)
        result = self.analyze(compatibility,
                              "-Compatibility", "PowerShell7")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_exclusions_empty_selection_and_module_dependency(self):
        self.fixture.write_text("Write-Output 'valid'\n", encoding="utf-8")
        for folder in ("vendor", "generated", "build", "archives"):
            excluded = self.fixture.parent / folder
            excluded.mkdir()
            (excluded / "broken.ps1").write_text("function Broken {\n", encoding="utf-8")
        result = subprocess.run(["pwsh", "-NoProfile", "-File", str(self.runner),
                                 "-Paths", "tools/fixtures"], cwd=self.root.parent,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        result = subprocess.run(["pwsh", "-NoProfile", "-File", str(self.runner),
                                 "-Paths", "tools/fixtures/vendor"], cwd=self.root.parent,
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Zero eligible", result.stdout + result.stderr)
        # Change module lookup only in this disposable process; no profile edits.
        runner = str(self.runner).replace("'", "''")
        result = subprocess.run(["pwsh", "-NoProfile", "-Command",
                                 "$env:PSModulePath = ''; & '" + runner + "'"],
                                cwd=self.root.parent, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Install-Module PSScriptAnalyzer", result.stdout + result.stderr)
        # A discoverable module with the wrong version must also fail at loading.
        old_module = self.root / "modules/PSScriptAnalyzer/0.0.1"
        old_module.mkdir(parents=True)
        (old_module / "PSScriptAnalyzer.psd1").write_text(
            "@{ ModuleVersion = '0.0.1'; RootModule = 'PSScriptAnalyzer.psm1' }\n",
            encoding="utf-8")
        (old_module / "PSScriptAnalyzer.psm1").write_text("", encoding="utf-8")
        modules = str(self.root / "modules").replace("'", "''")
        result = subprocess.run(["pwsh", "-NoProfile", "-Command",
                                 "$env:PSModulePath = '" + modules + "'; & '" + runner + "'"],
                                cwd=self.root.parent, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PSScriptAnalyzer 1.25.0 is required", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
