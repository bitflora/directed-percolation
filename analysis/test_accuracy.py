"""Accuracy regressions: python -m unittest discover -s analysis -p 'test_*.py'."""
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import flint

import crt
import da
import gate
from series_input import load_record
import wall_analysis


ROOT = Path(__file__).resolve().parents[1]


class SeriesInputTests(unittest.TestCase):
    def test_excess_order_and_short_inputs_are_rejected(self):
        ser = json.loads((ROOT / "results/series/bond171.json").read_text())["bond_171"]
        self.assertEqual(len(gate.derived(ser, 171)["R1"]), 171)
        with self.assertRaises(ValueError):
            gate.derived(ser, 200)
        ser["mu20"] = ser["mu20"][:100]
        with self.assertRaises(ValueError):
            gate.derived(ser, 171)
        for a, b in (([1], [1, 2]), ([1, 2], [1])):
            with self.assertRaises(ValueError):
                gate.series_div(a, b, 2)

    def test_wall_ratios_respect_available_order(self):
        ser = json.loads((ROOT / "results/series/wall67.json").read_text())["bond_wall_67"]
        self.assertEqual(len(wall_analysis.all_series(ser, 67)["Rx"]), 67)
        with self.assertRaises(ValueError):
            wall_analysis.all_series(ser, 68)

    def test_ratio_coefficients_are_exact(self):
        self.assertEqual(gate.series_div([1, 0, 0, 0], [1, -2, 0, 0], 4), [1, 2, 4, 8])

    def test_unverified_reconstruction_requires_explicit_overrides(self):
        with tempfile.TemporaryDirectory() as td:
            runs = Path(td) / "runs.jsonl"
            output = Path(td) / "series.json"
            runs.write_text("\n".join(json.dumps({"model": "bond", "n": 1,
                            "prime": p, "S": [1, 100 % p]}) for p in (5, 7)))
            output.write_text("existing result")
            cmd = [sys.executable, str(ROOT / "analysis/crt.py"), str(runs), "-o", str(output)]
            failed = subprocess.run(cmd, capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("CRT verification failed", failed.stderr)
            self.assertEqual(output.read_text(), "existing result")
            subprocess.run(cmd + ["--allow-unverified"], check=True, capture_output=True)
            with self.assertRaises(ValueError):
                load_record(output, "bond_1")
            self.assertFalse(load_record(output, "bond_1", allow_unverified=True)["verified"])

    def test_verified_reconstruction_still_writes_and_loads(self):
        with tempfile.TemporaryDirectory() as td:
            runs, output = Path(td) / "runs.jsonl", Path(td) / "series.json"
            runs.write_text("\n".join(json.dumps({"model": "bond", "n": 1,
                            "prime": p, "S": [1, 100]}) for p in (1009, 1013)))
            subprocess.run([sys.executable, str(ROOT / "analysis/crt.py"), str(runs),
                            "-o", str(output)], check=True, capture_output=True)
            self.assertEqual(load_record(output, "bond_1")["S"], [1, 100])


class ApproximantTests(unittest.TestCase):
    def setUp(self):
        flint.ctx.prec = da.PREC_BITS

    def test_known_singularity_and_exponent(self):
        # (1 - 4p)^(-1/2): an independent exact reference.
        series = [math.comb(2*n, n) for n in range(12)]
        estimate = da.analyse_one(series, 1, 0, (1, 1), 0.25)
        self.assertIsNotNone(estimate)
        self.assertEqual(estimate.xc, 0.25)
        self.assertEqual(estimate.exponent, 0.5)

    def test_root_encloses_exact_irrational_value(self):
        root = da._refine_root(flint.arb_poly([-2, 0, 1]), 1.4)
        self.assertIsNotNone(root)
        self.assertTrue(root.overlaps(flint.arb(2).sqrt()))
        self.assertTrue(da._rel_ok(root))

    def test_nonreal_and_multiple_roots_are_rejected(self):
        self.assertIsNone(da._refine_root(flint.arb_poly([1, 0, 1]), 0.1))
        self.assertIsNone(da._refine_root(flint.arb_poly([1, -2, 1]), 0.99))

    def test_uncertain_root_is_rejected(self):
        poly = flint.arb_poly([flint.arb(1, "1e-12"), -2])
        self.assertIsNone(da._refine_root(poly, 0.5))

    def test_uncertain_exponent_is_rejected(self):
        # A precise Q_K does not guarantee a precise Q_{K-1} or exponent.
        Q = [[flint.arb(-1, "1e-5")], [flint.arb(1), flint.arb(-2)]]
        with patch.object(da, "approximant", return_value=(Q, [flint.arb(0)])):
            self.assertIsNone(da.analyse_one([1]*4, 1, 0, (0, 1), 0.5))


if __name__ == "__main__":
    unittest.main()
