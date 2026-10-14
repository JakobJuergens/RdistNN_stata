import re
import subprocess
import tempfile
import unittest
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class PackageMaintenanceTests(unittest.TestCase):
    def test_master_runner_registers_every_regression_file_once(self):
        text = (ROOT / "tests/tests.do").read_text(encoding="utf-8")
        match = re.search(r"local test_files\s+///(.*?)// .*run each file", text, re.S)
        self.assertIsNotNone(match, "master runner's test inventory is missing")
        registered = re.findall(r"\btest-[a-z0-9-]+\b", match.group(1))
        expected = {p.stem for p in (ROOT / "tests").glob("test-*.do")}
        self.assertEqual(set(registered), expected)
        self.assertEqual(len(registered), len(set(registered)))

    def test_source_archive_has_public_docs_and_installation_payload(self):
        with tempfile.TemporaryDirectory(prefix="rdistnn-archive-test-") as directory:
            archive = Path(directory) / "source.zip"
            subprocess.run(
                ["bash", "tools/ci/build_source_archive.sh", str(archive)],
                cwd=ROOT, check=True, capture_output=True, text=True,
            )
            with zipfile.ZipFile(archive) as source:
                names = source.namelist()
                self.assertEqual(len(names), len(set(names)))
                required = {
                    "README.md", "RELEASE_NOTES.md", "rdistnn.pkg", "stata.toc",
                    "docs/command-guide.md", "examples/getting-started.do",
                    "src/rdistnn.plugin", "src/plugin/rdistnn_plugin.cpp",
                }
                required.update(re.findall(
                    r"^f\s+(\S+)", (ROOT / "rdistnn.pkg").read_text(), re.M,
                ))
                self.assertTrue(required.issubset(names), required.difference(names))
                for name in names:
                    self.assertFalse(name.startswith((
                        ".git/", ".github/", "_dist/", "archives/", "GH_logs/",
                        "markdown_notes/", "tests/", "tools/",
                    )), name)
                    self.assertFalse(name.endswith((".log", ".o", ".obj", ".mlib", ".mo")), name)
                    self.assertNotEqual(name, "src/plugin/rdistnn.plugin")
                    self.assertEqual(source.read(name), (ROOT / name).read_bytes())


if __name__ == "__main__":
    unittest.main()
