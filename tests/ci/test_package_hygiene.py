import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('package_hygiene', ROOT / 'tools/ci/package_hygiene.py')


class PackageHygieneTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.checker = importlib.util.module_from_spec(SPEC)
        SPEC.loader.exec_module(cls.checker)

    def test_internal_ado_files_are_required_installation_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'src').mkdir()
            (root / 'src/_private.ado').write_text('version 17\n')
            (root / 'rdistnn.pkg').write_text('v 3\n')
            self.assertTrue(any('_private.ado' in x for x in self.checker.package_issues(root)))
            (root / 'rdistnn.pkg').write_text('v 3\nf src/_private.ado\n')
            self.assertEqual(self.checker.package_issues(root), [])

    def test_manifest_rejects_duplicates_and_nonlocal_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'src').mkdir()
            (root / 'src/a.ado').write_text('version 17\n')
            (root / 'rdistnn.pkg').write_text('f src/a.ado\nf src/a.ado\nf ../outside\n')
            errors = self.checker.package_issues(root)
            self.assertTrue(any('duplicate' in x for x in errors))
            self.assertTrue(any('local' in x for x in errors))
            (root / 'rdistnn.pkg').write_text('f src/a.ado\nf src/./a.ado\n')
            errors = self.checker.package_issues(root)
            self.assertTrue(any('duplicate' in x for x in errors))
            self.assertTrue(any('canonical' in x for x in errors))

    def test_mata_functions_use_the_package_prefix(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'src').mkdir()
            helper = root / 'src/_weights.ado'
            helper.write_text('version 17\nmata:\n  real scalar lncomb(real scalar n)\n{ return(n) }\nend\n')
            (root / 'rdistnn.pkg').write_text('f src/_weights.ado\n')
            self.assertTrue(any('lncomb' in x for x in self.checker.package_issues(root)))
            helper.write_text(helper.read_text().replace('lncomb', 'rdistnn_fixture_log_combination'))
            self.assertEqual(self.checker.package_issues(root), [])

    def test_mata_functions_have_one_source_owner(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'src').mkdir()
            declaration = 'version 17\nmata:\nreal scalar rdistnn_probe() { return(1) }\nend\n'
            (root / 'src/a.ado').write_text(declaration)
            other = root / 'src/b.ado'
            other.write_text(declaration)
            (root / 'rdistnn.pkg').write_text('f src/a.ado\nf src/b.ado\n')
            errors = self.checker.package_issues(root)
            self.assertTrue(any('duplicate Mata function' in x and 'rdistnn_probe' in x
                                and 'a.ado' in x and 'b.ado' in x for x in errors))
            other.write_text(declaration.replace('rdistnn_probe', 'rdistnn_other'))
            self.assertEqual(self.checker.package_issues(root), [])

    def test_shared_headers_require_both_copies_and_identical_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            stata, r = base / 'stata', base / 'R'
            (stata / 'src/plugin').mkdir(parents=True)
            (r / 'inst/include').mkdir(parents=True)
            for name in self.checker.SHARED_HEADERS:
                (stata / 'src/plugin' / name).write_bytes(b'canonical\n')
                (r / 'inst/include' / name).write_bytes(b'canonical\n')
            self.assertEqual(self.checker.shared_header_issues(stata, r), [])
            for name in self.checker.SHARED_HEADERS:
                with self.subTest(header=name):
                    path = r / 'inst/include' / name
                    path.write_bytes(b'drift\n')
                    self.assertEqual(self.checker.shared_header_issues(stata, r),
                                     ['mirrored header differs: ' + name])
                    path.unlink()
                    self.assertEqual(self.checker.shared_header_issues(stata, r),
                                     ['missing mirrored header: ' + name])
                    path.write_bytes(b'canonical\n')

    def test_mata_modules_declare_the_supported_version_floor(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'src').mkdir()
            helper = root / 'src/_weights.ado'
            (root / 'rdistnn.pkg').write_text('f src/_weights.ado\n')
            for version in ('', 'version 16.1\n'):
                helper.write_text(version + 'mata:\nend\n')
                self.assertTrue(any('version' in x for x in self.checker.package_issues(root)))
            helper.write_text('version 17.1\nmata:\nend\n')
            self.assertEqual(self.checker.package_issues(root), [])

    def test_current_installation_inventory_is_complete(self):
        self.assertEqual(self.checker.package_issues(ROOT), [])
        self.assertEqual(self.checker.package_version_issues(ROOT), [])

    def write_version_fixture(self, root):
        (root / 'src').mkdir()
        (root / 'docs').mkdir()
        (root / 'rdistnn.pkg').write_text(
            'v 3\nd Version: 0.1.3.9001\nf src/rdistnn.sthlp\n')
        (root / 'README.md').write_text('Version **0.1.3.9001** (development).\n')
        (root / 'RELEASE_NOTES.md').write_text(
            '# Release notes\n\n## 0.1.3.9001 — development\n\n## 0.1.3 — historical\n')
        (root / 'src/rdistnn.sthlp').write_text(
            '{smcl}\n{* *! version 0.1.3.9001  2026-10-21}{...}\n'
            'Stata 17 or later. Package {bf:0.1.3.9001} provides commands.\n')

    def run_hygiene(self, root):
        return subprocess.run(
            [sys.executable, str(ROOT / 'tools/ci/package_hygiene.py'),
             '--repo-root', str(root)], capture_output=True, text=True,
        )

    def test_current_package_version_owners_must_agree(self):
        for owner in ('rdistnn.pkg', 'README.md', 'RELEASE_NOTES.md', 'src/rdistnn.sthlp'):
            with self.subTest(owner=owner), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                self.write_version_fixture(root)
                path = root / owner
                path.write_text(path.read_text().replace('0.1.3.9001', '0.1.3.9000', 1))
                result = self.run_hygiene(root)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn('version', result.stderr.lower())

    def test_stale_current_help_and_article_claims_are_reported(self):
        for target, text in (
            ('src/weights.sthlp', 'This help describes package {bf:0.1.3.9000}.\n'),
            ('docs/guide.md', 'This article documents package **0.1.3**.\n'),
            ('docs/crosswalk.md', 'This map covers the current **0.1.3.9000 development** interfaces.\n'),
        ):
            with self.subTest(target=target), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                self.write_version_fixture(root)
                (root / target).write_text(text)
                if target.startswith('src/'):
                    with (root / 'rdistnn.pkg').open('a') as manifest:
                        manifest.write('f ' + target + '\n')
                result = self.run_hygiene(root)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn(target, result.stderr)
                self.assertIn('0.1.3.9001', result.stderr)

    def test_manifest_version_and_first_release_heading_are_unambiguous(self):
        cases = (
            ('rdistnn.pkg', 'd Version: 0.1.3.9001\nd Version: invalid\n'),
            ('RELEASE_NOTES.md', '# Release notes\n\n## upcoming\n\n## 0.1.3.9001 — historical\n'),
        )
        for owner, text in cases:
            with self.subTest(owner=owner), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                self.write_version_fixture(root)
                if owner == 'rdistnn.pkg':
                    text = 'v 3\n' + text + 'f src/rdistnn.sthlp\n'
                (root / owner).write_text(text)
                result = self.run_hygiene(root)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn(owner, result.stderr)

    def test_all_installed_source_stamps_must_match_the_package_version(self):
        for name, valid in (
            ('src/command.ado', '*! version 0.1.3.9001  2026-10-23\nversion 17\n'),
            ('src/_helper.ado', '*! version 0.1.3.9001  2026-10-23\nversion 17\n'),
            ('src/command.sthlp', '{smcl}\n{* *! version 0.1.3.9001  2026-10-23}{...}\n'),
        ):
            for defect in ('stale', 'missing', 'duplicate', 'malformed', 'misplaced'):
                with self.subTest(name=name, defect=defect), tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    self.write_version_fixture(root)
                    lines = valid.splitlines(keepends=True)
                    index = 1 if name.endswith('.sthlp') else 0
                    if defect == 'stale':
                        lines[index] = lines[index].replace('0.1.3.9001', '0.2.1')
                    elif defect == 'missing':
                        del lines[index]
                    elif defect == 'duplicate':
                        lines.insert(index + 1, lines[index])
                    elif defect == 'malformed':
                        lines[index] = lines[index].replace('2026-10-23', '2026-02-30')
                    else:
                        lines.insert(index, '\n')
                    (root / name).write_text(''.join(lines))
                    with (root / 'rdistnn.pkg').open('a') as manifest:
                        manifest.write('f ' + name + '\n')
                    result = self.run_hygiene(root)
                    self.assertNotEqual(result.returncode, 0, result.stdout)
                    self.assertIn(name, result.stderr)
                    self.assertIn('stamp', result.stderr.lower())

    def test_version_hygiene_preserves_history_language_schema_and_pins(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_version_fixture(root)
            (root / 'src/current_header.ado').write_text(
                '*! version 0.1.3.9001  2026-10-23\nversion 17.0\n'
                'program define example\n    ereturn scalar model_version=4\nend\n')
            with (root / 'rdistnn.pkg').open('a') as manifest:
                manifest.write('f src/current_header.ado\n')
            (root / 'docs/guide.md').write_text(
                'Production version **0.1.3** remains pinned at v0.1.3.\n'
                'Refit schema 3 from version 0.1.3 before reuse.\n')
            (root / 'tests').mkdir()
            (root / 'tests/provenance.do').write_text('// R source version: 0.1.2.9000\n')
            result = self.run_hygiene(root)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_stamp_sync_preserves_bodies_newlines_and_other_notes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_version_fixture(root)
            command = root / 'src/command.ado'
            helper = root / 'src/_helper.ado'
            command_body = (b'*! Diagnostic command\r\nversion 17.0\r\n'
                            b'program define example\r\n    ereturn scalar model_version=4\r\nend\r\n')
            helper_body = b'version 17\nmata:\nvoid rdistnn_helper() {}\nend\n'
            command.write_bytes(b'*! version 0.0.9  12jul2026  Original Author\r\n' + command_body)
            helper.write_bytes(helper_body)
            with (root / 'rdistnn.pkg').open('a') as manifest:
                manifest.write('f src/command.ado\nf src/_helper.ado\n')
            help_path = root / 'src/rdistnn.sthlp'
            help_before = help_path.read_bytes()
            help_path.write_bytes(help_before.replace(b'0.1.3.9001  2026-10-21', b'0.2.1  17oct2026'))
            paths = (command, helper, help_path)
            before = {p: p.read_bytes() for p in paths}
            self.assertNotEqual(self.run_hygiene(root).returncode, 0)
            self.assertEqual({p: p.read_bytes() for p in paths}, before, 'Checking must be read-only')
            args = [sys.executable, str(ROOT / 'tools/ci/package_hygiene.py'),
                    '--repo-root', str(root), '--sync-version-stamps', '--stamp-date', '2026-10-23']
            result = subprocess.run(args, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(command.read_bytes(),
                             b'*! version 0.1.3.9001  2026-10-23  Original Author\r\n' + command_body)
            self.assertEqual(helper.read_bytes(), b'*! version 0.1.3.9001  2026-10-23\n' + helper_body)
            self.assertEqual(help_path.read_bytes(), help_before.replace(b'2026-10-21', b'2026-10-23'))
            after = {p: p.read_bytes() for p in paths}
            result = subprocess.run(args, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('Synchronized 0', result.stdout)
            self.assertEqual({p: p.read_bytes() for p in paths}, after)

    def test_ambiguous_stamp_sync_changes_no_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_version_fixture(root)
            helper = root / 'src/_helper.ado'
            helper.write_text('version 17\n')
            with (root / 'rdistnn.pkg').open('a') as manifest:
                manifest.write('f src/_helper.ado\n')
            help_path = root / 'src/rdistnn.sthlp'
            help_path.write_text(help_path.read_text() + '{* *! version 0.0.9  2026-10-23}{...}\n')
            before = {p: p.read_bytes() for p in (helper, help_path)}
            result = subprocess.run(
                [sys.executable, str(ROOT / 'tools/ci/package_hygiene.py'), '--repo-root', str(root),
                 '--sync-version-stamps', '--stamp-date', '2026-10-23'], capture_output=True, text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('ambiguous', result.stderr)
            self.assertEqual({p: p.read_bytes() for p in before}, before)

    def test_source_stamps_accept_production_and_development_versions(self):
        for version in ('0.1.4', '0.1.4.9000'):
            with self.subTest(version=version), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                self.write_version_fixture(root)
                for path in root.rglob('*'):
                    if path.is_file():
                        path.write_text(path.read_text().replace('0.1.3.9001', version))
                result = self.run_hygiene(root)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_custom_error_rendering_uses_the_existing_host_owners(self):
        cases = (
            ('command.ado', 'program define command\n    di as error "bad option"\n    exit 198\nend\n'),
            ('_helper.ado', 'mata:\nvoid rdistnn_probe() { errprintf("bad input\\n"); exit(198) }\nend\n'),
        )
        for name, body in cases:
            with self.subTest(source=name), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                self.write_version_fixture(root)
                (root / 'src' / name).write_text('*! version 0.1.3.9001  2026-10-23\nversion 17\n' + body)
                with (root / 'rdistnn.pkg').open('a') as manifest:
                    manifest.write('f src/' + name + '\n')
                result = self.run_hygiene(root)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn(name, result.stderr)
                self.assertIn('render', result.stderr.lower())

    def test_native_validation_and_standard_errors_remain_supported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_version_fixture(root)
            (root / 'src/command.ado').write_text(
                '*! version 0.1.3.9001  2026-10-23\nversion 17\n'
                'program define command\n    syntax varlist(numeric)\n'
                '    confirm matrix targets\n    if _N==0 error 2000\nend\n')
            with (root / 'rdistnn.pkg').open('a') as manifest:
                manifest.write('f src/command.ado\n')
            result = self.run_hygiene(root)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_error_policy_distinguishes_owners_bootstrap_and_literal_examples(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'src').mkdir()
            owner = root / 'src/_rdistnn_messages.ado'
            owner.write_text('display as error _asis `"`message\'"\'\n')
            bootstrap = root / 'src/_rdistnn_mata.ado'
            bootstrap.write_text('display as error _asis "Loader failed (`plugin_rc\')."\n'
                                 'exit `plugin_rc\'\n')
            caller = root / 'src/command.ado'
            caller.write_text('_rdistnn_messages error, message(`"Use display as error or errprintf()."\')\n'
                              '// display as error "comment"\n')
            self.assertEqual(self.checker.error_rendering_issues(root), [])
            bootstrap.write_text(bootstrap.read_text().replace('exit `plugin_rc\'', 'exit 198'))
            errors = self.checker.error_rendering_issues(root)
            self.assertEqual(len(errors), 1)
            self.assertIn('_rdistnn_mata.ado', errors[0])
