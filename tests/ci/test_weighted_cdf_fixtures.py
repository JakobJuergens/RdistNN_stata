"""Independently rounded boundary references and checked-in Mata output."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "weighted_cdf_fixtures", ROOT / "tools/generate_weighted_cdf_fixtures.py")
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)


class WeightedCdfFixtureTests(unittest.TestCase):
    def test_rational_references_and_mata_fixture(self):
        with tempfile.TemporaryDirectory() as temporary:
            products = GENERATOR.products(ROOT, Path(temporary))
        mata = ROOT / "tests/fixtures/weighted_cdf_fixtures.do"
        self.assertEqual(mata.read_text(), products[mata])
        self.assertEqual(products[mata].count("test_assert"), 3)

    def test_literal_boundary_oracle(self):
        self.assertEqual(GENERATOR.compensated([2**54, 1, -2**54]), 1)
