import hashlib
import posixpath
import re
import shutil
import subprocess
import tempfile
import unittest
import zipfile
from pathlib import Path
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[2]


class PackageMaintenanceTests(unittest.TestCase):
    def test_source_archive_respects_unstaged_deletions(self):
        with tempfile.TemporaryDirectory(prefix="rdistnn-deleted-source-test-") as directory:
            checkout = Path(directory)
            for name in ("src", "docs", "examples", "tools/ci"):
                shutil.copytree(ROOT / name, checkout / name)
            for name in (
                "README.md", "RELEASE_NOTES.md", "LICENSE.md", "LINKING_EXCEPTION.md",
                "rdistnn.pkg", "stata.toc", ".gitignore",
            ):
                shutil.copyfile(ROOT / name, checkout / name)
            subprocess.run(["git", "init", "--quiet"], cwd=checkout, check=True)
            retired = checkout / "src/retired_test_helper.ado"
            retired.write_text("version 17\n")
            subprocess.run(["git", "add", "--", str(retired)], cwd=checkout, check=True)
            retired.unlink()
            archive = checkout / "source.zip"
            subprocess.run(
                ["bash", "tools/ci/build_source_archive.sh", str(archive)],
                cwd=checkout, check=True, capture_output=True, text=True,
            )
            with zipfile.ZipFile(archive) as source:
                self.assertNotIn("src/retired_test_helper.ado", source.namelist())
                self.assertEqual(source.read("src/cate_dnn.ado"),
                                 (ROOT / "src/cate_dnn.ado").read_bytes())

    def test_release_source_download_has_verifiable_checksum_and_fresh_payload(self):
        with tempfile.TemporaryDirectory(prefix="rdistnn-release-test-") as directory:
            output = Path(directory)
            archive = output / "rdistnn-source.zip"
            checksums = output / "SHA256SUMS-source.txt"
            archive.write_bytes(b"stale release archive")
            checksums.write_text("stale checksum\n")
            subprocess.run(
                ["bash", "tools/ci/build_release_assets.sh", str(output)],
                cwd=ROOT, check=True, capture_output=True, text=True,
            )
            self.assertEqual(
                checksums.read_text(),
                hashlib.sha256(archive.read_bytes()).hexdigest()
                + "  rdistnn-source.zip\n",
            )
            subprocess.run(
                ["sha256sum", "--check", checksums.name],
                cwd=output, check=True, capture_output=True, text=True,
            )
            with zipfile.ZipFile(archive) as source:
                self.assertEqual(
                    source.read("rdistnn.pkg"), (ROOT / "rdistnn.pkg").read_bytes(),
                )
                self.assertEqual(
                    source.read("src/cate_dnn.ado"),
                    (ROOT / "src/cate_dnn.ado").read_bytes(),
                )

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
                    "README.md", "RELEASE_NOTES.md", "LICENSE.md", "rdistnn.pkg", "stata.toc",
                    "LINKING_EXCEPTION.md",
                    "docs/command-guide.md", "examples/getting-started.do",
                    "src/rdistnn.plugin", "src/plugin/rdistnn_plugin.cpp",
                }
                required.update(re.findall(
                    r"^[fF]\s+(\S+)", (ROOT / "rdistnn.pkg").read_text(), re.M,
                ))
                self.assertTrue(required.issubset(names), required.difference(names))
                # A relative contributor-guide link can work in a checkout
                # while breaking in the public archive, which omits that file.
                for name in names:
                    if not name.endswith(".md"):
                        continue
                    content = source.read(name).decode("utf-8")
                    for target in re.findall(r"\[[^\]]*\]\(([^\s)]+)\)", content):
                        parsed = urlsplit(target)
                        if parsed.scheme or parsed.netloc or not parsed.path:
                            continue
                        resolved = posixpath.normpath(posixpath.join(
                            posixpath.dirname(name), unquote(parsed.path),
                        ))
                        self.assertIn(resolved, names, f"{name}: broken archive link {target}")
                for name in names:
                    self.assertFalse(name.startswith((
                        ".git/", ".github/", "_dist/", "archives/", "GH_logs/",
                        "markdown_notes/", "tests/", "tools/",
                    )), name)
                    self.assertFalse(name.endswith((".log", ".o", ".obj", ".mlib", ".mo")), name)
                    self.assertNotEqual(name, "src/plugin/rdistnn.plugin")
                    self.assertEqual(source.read(name), (ROOT / name).read_bytes())

    def test_license_notices_are_installed_with_package_specific_names(self):
        manifest = (ROOT / "rdistnn.pkg").read_text()
        for original, installed in (
            ("LICENSE.md", "rdistnn_license.md"),
            ("LINKING_EXCEPTION.md", "rdistnn_linking_exception.md"),
        ):
            self.assertIn("F src/" + installed, manifest.splitlines())
            expected = (ROOT / original).read_bytes().replace(
                b"(LICENSE.md)", b"(rdistnn_license.md)",
            ).replace(
                b"(LINKING_EXCEPTION.md)", b"(rdistnn_linking_exception.md)",
            )
            self.assertEqual(
                (ROOT / "src" / installed).read_bytes(),
                expected,
            )


if __name__ == "__main__":
    unittest.main()
