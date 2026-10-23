import importlib.util
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
            helper.write_text(helper.read_text().replace('lncomb', 'rdistnn_log_binomial'))
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
