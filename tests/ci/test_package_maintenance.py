import hashlib
import posixpath
import re
import shutil
import struct
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

    def test_archive_checker_uses_archived_manifest_and_rejects_nonpublic_paths(self):
        with tempfile.TemporaryDirectory(prefix="rdistnn-archive-validation-") as directory:
            root = Path(directory)
            baseline = root / "baseline.zip"
            subprocess.run(
                ["bash", "tools/ci/build_source_archive.sh", str(baseline)],
                cwd=ROOT, check=True, capture_output=True, text=True,
            )
            with zipfile.ZipFile(baseline) as source:
                contents = {info.filename: source.read(info) for info in source.infolist()}
            with self.subTest(case="relative-archive-outside-checkout"):
                result = subprocess.run(
                    ["bash", str(ROOT / "tools/ci/archive_check.sh"), baseline.name],
                    cwd=root, capture_output=True, text=True,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            manifest = contents["rdistnn.pkg"]
            cases = (
                ("unchanged", {}, True),
                # An archive can legitimately have a different installation
                # inventory from the checkout used to run the verifier.
                ("archive-only-helper", {
                    "rdistnn.pkg": manifest + b"f src/_rdistnn_archive_only.ado\n",
                    "src/_rdistnn_archive_only.ado": b"version 17\n",
                }, True),
                ("missing-declared-input", {
                    "rdistnn.pkg": manifest + b"f src/missing_required.ado\n",
                }, False),
                ("nonlocal-declaration", {
                    "rdistnn.pkg": manifest + b"f ../unrelated.ado\n",
                }, False),
                ("unregistered-helper", {
                    "src/_rdistnn_unregistered.ado": b"version 17\n",
                }, False),
                ("parent-path", {"../unrelated.ado": b"version 17\n"}, False),
                ("absolute-path", {"/tmp/unrelated.ado": b"version 17\n"}, False),
                ("backslash-path", {"src\\unrelated.ado": b"version 17\n"}, False),
                ("aliased-path", {"src/./unrelated.ado": b"version 17\n"}, False),
                ("private-notes", {".codex/private-note.md": b"private notes\n"}, False),
                ("checkout-tests", {"tests/test-extra.do": b"version 17\n"}, False),
                ("generated-object", {"src/plugin/temporary.o": b"native object"}, False),
            )
            for label, changes, should_pass in cases:
                with self.subTest(case=label):
                    archive = root / (label + ".zip")
                    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as source:
                        for name, content in (contents | changes).items():
                            source.writestr(name, content)
                    result = subprocess.run(
                        ["bash", "tools/ci/archive_check.sh", str(archive)], cwd=ROOT,
                        capture_output=True, text=True,
                    )
                    self.assertEqual(result.returncode == 0, should_pass,
                                     result.stdout + result.stderr)
            with self.subTest(case="corrupt-native-payload"):
                archive = root / "corrupt-payload.zip"
                with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_STORED) as source:
                    for name, content in contents.items():
                        source.writestr(name, content)
                with zipfile.ZipFile(archive) as source:
                    info = source.getinfo("src/rdistnn.plugin")
                payload = bytearray(archive.read_bytes())
                name_length, extra_length = struct.unpack_from("<HH", payload, info.header_offset + 26)
                data_offset = info.header_offset + 30 + name_length + extra_length
                # Change file data, preserving the inventory and stored CRC.
                payload[data_offset + 10] ^= 1
                archive.write_bytes(payload)
                result = subprocess.run(
                    ["bash", "tools/ci/archive_check.sh", str(archive)], cwd=ROOT,
                    capture_output=True, text=True,
                )
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("corrupt archive input: src/rdistnn.plugin", result.stderr)

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
