"""Biased differential-approximant scan: exponent estimates as a function of an
assumed critical point.

Assuming universality, the site exponents must equal the bond exponents. Jensen
1999 (section 4.2) estimated the site critical point by finding where the
biased site exponents cross the bond values; this script reproduces that
procedure.

Usage:
  python biased_scan.py results/series/site136.json site_136 --xc 0.7054850 0.7054856 7 \
      --target gamma=2.277730 gamma+2nu_perp=4.471425 -o results/analysis/site_biased.json
"""
import argparse
import json
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
import da  # noqa: E402
from gate import TARGET, derived  # noqa: E402
from series_input import load_record  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("series_json")
    ap.add_argument("key")
    ap.add_argument("--xc", nargs=3, required=True, help="start stop count")
    ap.add_argument("--names", nargs="*", default=["S", "mu20", "R1", "R2"])
    ap.add_argument("--K", type=int, default=2)
    ap.add_argument("--span", type=int, default=10)
    ap.add_argument("--target", nargs="*", default=[], help="quantity=value pairs to locate crossings")
    ap.add_argument("-o")
    ap.add_argument("--allow-unverified", action="store_true", help="analyse exploratory, unverified series")
    a = ap.parse_args()
    ser = load_record(a.series_json, a.key, a.allow_unverified)
    T = ser["order"]
    der = derived(ser, T)
    allser = {"S": ser["S"], "mu01": ser["mu01"], "mu02": ser["mu02"], "R1": der["R1"]}
    if "mu20" in ser:  # absent for --no-x runs
        allser["mu20"] = ser["mu20"]
        allser["R2"] = der["R2"]
    lo, hi, cnt = Fraction(a.xc[0]), Fraction(a.xc[1]), int(a.xc[2])
    grid = [lo + (hi - lo) * Fraction(i, cnt - 1) for i in range(cnt)]
    targets = dict(kv.split("=") for kv in a.target)
    out = []
    for name in a.names:
        s = allser[name]
        Tn = len(s)
        row = {"series": name, "quantity": TARGET[name], "xc": [], "exponent": [], "spread": []}
        ensembles = da.biased_ensembles(s, a.K, range(0, 11, 2), range(Tn - a.span, Tn + 1), grid)
        for xc, v in zip(grid, ensembles):
            med = float(np.median(v)) if v else float("nan")
            spr = float((np.quantile(v, 0.84) - np.quantile(v, 0.16)) / 2) if v else float("nan")
            row["xc"].append(float(xc))
            row["exponent"].append(med)
            row["spread"].append(spr)
            print(f"{name} xc={float(xc):.10f} exponent={med:.7f} ±{spr:.1e} ({len(v)})", flush=True)
        q = TARGET[name]
        if q in targets:
            tv = float(targets[q])
            xs, es = row["xc"], row["exponent"]
            for i in range(len(xs) - 1):
                if (es[i] - tv) * (es[i + 1] - tv) <= 0 and es[i] != es[i + 1]:
                    x = xs[i] + (tv - es[i]) * (xs[i + 1] - xs[i]) / (es[i + 1] - es[i])
                    row["crossing"] = x
                    print(f"  {name}: biased exponent equals {q} = {tv} at x_c = {x:.10f}")
        out.append(row)
    if a.o:
        json.dump(out, open(a.o, "w"), indent=1)


if __name__ == "__main__":
    main()
