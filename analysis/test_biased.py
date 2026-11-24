"""Numerical regression references captured before solve reuse was introduced."""
from fractions import Fraction
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import flint
import da
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


class BiasedTests(unittest.TestCase):
    def test_bond_site_and_wall_ensembles_keep_values_and_acceptance(self):
        cases = [
            ("bond171", "bond_171", "S", 2, ["0.644700180", "0.644700190"], [
                [2.2780483405856033, 2.2782627993547426, 2.2773286973976314, 2.277401302300671,
                 2.2776446541964823, 2.2777390866047114, 2.277730144293846, 2.277835343084708],
                [2.27805404571924, 2.278271523755124, 2.2773294097145578, 2.277403419628594,
                 2.2776502355311323, 2.277744300038565, 2.277735394034006, 2.2778407928963142]]),
            ("site158", "site_158", "mu20", 3, ["0.7054851", "0.7054853"], [
                [4.471254784110488, 4.471054288737057, 4.47107941035045, 4.470791186903231,
                 4.471047174623538, 4.4718394456136, 4.470670789984452, 4.471164871153595,
                 4.470863530685389, 4.471092529373154, 4.471190650650733, 4.471270523494312,
                 4.471127587741765, 4.4711874890751195],
                [4.471357733301024, 4.4711538068530805, 4.471178077141527, 4.470877346668767,
                 4.4711475666132365, 4.4719028552151645, 4.470760387345885, 4.4712758526277625,
                 4.470993051629257, 4.471208327058925, 4.471299058278803, 4.471370444745012,
                 4.471240485713192, 4.471296151307831]]),
            ("wall110", "bond_wall_110", "T", 2, ["0.644700180", "0.644700190"], [
                [1.0003634808466824, 1.0003429668809731, 1.0003414403782842, 1.000301966037848,
                 1.0003744906176528, 1.0003705384012336, 1.0004235068332332, 1.0004128954391085],
                [1.0003645753278017, 1.0003439604117534, 1.0003424604548525, 1.0003029162613217,
                 1.0003756711354148, 1.0003717328803146, 1.0004247667736688, 1.0004141019243828]])]
        for file, key, name, K, grid, expected in cases:
            series = json.loads((ROOT / f"results/series/{file}.json").read_text())[key][name]
            actual = da.biased_ensembles(series, K, [0, 2], [60, 62], map(Fraction, grid))
            for xc, values, reference in zip(grid, actual, expected):
                with self.subTest(series=file, xc=xc):
                    self.assertEqual(len(values), len(reference))
                    for v, r in zip(values, reference):
                        self.assertAlmostEqual(v, r, delta=1e-12)

    def test_bias_and_matched_series_equations_hold(self):
        series = json.loads((ROOT / "results/series/bond171.json").read_text())["bond_171"]["S"]
        K, L, degs = 2, 2, (18, 18, 19)
        xc = Fraction("0.644700185")
        theta = da._theta_powers(series, K)
        direct = da.approximant(series, K, L, degs, theta, bias=xc)
        basis = da._prepare_biased(series, K, L, degs, theta)
        self.assertIsNotNone(basis)
        updated = da._apply_bias(basis, K, L, degs, xc)
        self.assertIsNotNone(updated)
        for Q, P in (direct, updated):
            self.assertTrue(flint.arb_poly(Q[K])(da._to_arb(xc)).contains(0))
            for n in range(sum(d+1 for d in degs)+L-1):
                residual = sum(Q[i][j] * theta[i][n-j] for i in range(K+1)
                               for j in range(min(degs[i], n)+1)) - (P[n] if n <= L else 0)
                self.assertTrue(residual.contains(0), n)
        for qs, qr in zip(updated[0], direct[0]):
            for c, r in zip(qs, qr):
                self.assertTrue(c.overlaps(r))

    def test_known_biased_exponent_and_singular_system(self):
        series = [math.comb(2*n, n) for n in range(20)]
        self.assertEqual(da.biased_exponent(series, 1, 0, (1, 1), Fraction(1, 4)), 0.5)
        # Multiples of the exact ODE make the unbiased system singular.
        self.assertIsNone(da.approximant(series, 1, 0, (2, 2)))
        self.assertIsNotNone(da.approximant(series, 1, 0, (2, 2), bias=Fraction(3, 10)))
        self.assertIsNone(da._prepare_biased(series, 1, 0, (2, 2), da._theta_powers(series, 1)))

    def test_fallback_preserves_direct_results(self):
        series = json.loads((ROOT / "results/series/bond171.json").read_text())["bond_171"]["S"]
        grid = [Fraction("0.644700180"), Fraction("0.644700190")]
        expected = [da.biased_ensemble(series, 2, [0], [60], xc) for xc in grid]
        for name, result in (("_prepare_biased", None), ("_apply_bias", None), ("_biased_value", None)):
            if name == "_biased_value":
                # Only reject the stricter update, leaving the direct gate intact.
                original = da._biased_value
                replacement = lambda res, K, xc, digits=10: None if digits == 30 else original(res, K, xc, digits)
                context = patch.object(da, name, side_effect=replacement)
            else:
                context = patch.object(da, name, return_value=result)
            with self.subTest(fallback=name), context:
                self.assertEqual(da.biased_ensembles(series, 2, [0], [60], grid), expected)

    def test_empty_single_and_repeated_grids(self):
        series = [math.comb(2*n, n) for n in range(20)]
        xc = Fraction(1, 4)
        expected = da.biased_ensemble(series, 1, [0], [4], xc)
        self.assertEqual(da.biased_ensembles(series, 1, [0], [4], []), [])
        self.assertEqual(da.biased_ensembles(series, 1, [0], [4], [xc]), [expected])
        self.assertEqual(da.biased_ensembles(series, 1, iter([0]), iter([4]), [xc, xc]), [expected, expected])
        # Constant Q_K cannot satisfy a bias; its update denominator is zero.
        self.assertEqual(da.biased_ensembles(series, 1, [0], [2], [xc, xc]), [[], []])

    def test_full_length_bond_updates_match_direct_solves(self):
        series = json.loads((ROOT / "results/series/bond171.json").read_text())["bond_171"]["S"]
        grid = [Fraction("0.644700180"), Fraction("0.644700190")]
        actual = da.biased_ensembles(series, 2, [0], [166, 172], grid)
        expected = [da.biased_ensemble(series, 2, [0], [166, 172], xc) for xc in grid]
        for values, reference in zip(actual, expected):
            self.assertEqual(len(values), len(reference))
            for v, r in zip(values, reference):
                self.assertAlmostEqual(v, r, delta=1e-12)

    def test_scan_commands_preserve_ensemble_summaries(self):
        # Shortened verified inputs keep this an ordinary correctness test.
        for wall in (False, True):
            file, key, name = ("wall110", "bond_wall_110", "T") if wall else ("bond171", "bond_171", "S")
            original = json.loads((ROOT / f"results/series/{file}.json").read_text())[key]
            ser = {k: v[:61] if isinstance(v, list) else v for k, v in original.items()}
            ser["order"] = 60
            grid = [Fraction("0.644700180"), Fraction("0.644700190")]
            with tempfile.TemporaryDirectory() as td:
                inp, output = Path(td) / "series.json", Path(td) / "scan.json"
                inp.write_text(json.dumps({key: ser}))
                script = "wall_analysis.py" if wall else "biased_scan.py"
                args = [sys.executable, str(ROOT / "analysis" / script), str(inp), key,
                        "--names", name, "--K", "2", "--span", "0", "-o", str(output)]
                if wall:
                    args += ["--terms", "60", "--Lmax", "0", "--xc"] + [str(float(x)) for x in grid]
                else:
                    args += ["--xc", str(float(grid[0])), str(float(grid[1])), "2"]
                subprocess.run(args, check=True, capture_output=True, text=True)
                result = json.loads(output.read_text())
            for i, xc in enumerate(grid):
                values = da.biased_ensemble(ser[name], 2, [0] if wall else range(0, 11, 2), [61], xc)
                med = float(np.median(values))
                spread = float((np.quantile(values, 0.84)-np.quantile(values, 0.16))/2)
                row = result[i] if wall else result[0]
                self.assertAlmostEqual(row["exponent"] if wall else row["exponent"][i], med, delta=1e-12)
                self.assertAlmostEqual(row["spread"] if wall else row["spread"][i], spread, delta=1e-12)
                if wall:
                    self.assertEqual(row["count"], len(values))


if __name__ == "__main__":
    unittest.main()
