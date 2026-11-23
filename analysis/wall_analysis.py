"""Differential-approximant analysis of the wall series (dpwall output).

Series and expected singularities at p -> pc-:
  T   = sum_t P(t)         ~ (pc - p)^-tau_1          (the mean cluster length)
  S                        ~ (pc - p)^-gamma_1
  mt1 = sum t C            ~ (pc - p)^-(gamma_1 + nu_par)
  mt2 = sum t^2 C          ~ (pc - p)^-(gamma_1 + 2 nu_par)
  mx1 = sum x C            ~ (pc - p)^-(gamma_1 + nu_perp)
  mx2 = sum x^2 C          ~ (pc - p)^-(gamma_1 + 2 nu_perp)
  Rt  = mt2 / mt1          ~ (pc - p)^-nu_par
  Rx  = mx2 / mx1          ~ (pc - p)^-nu_perp        (Jensen 1999 II's best pc)

Two modes:
  unbiased: estimate pc and the exponent from each approximant (as gate.py);
  biased:   impose pc = each value in --xc and estimate the exponent only.

Usage:
  python wall_analysis.py series.json bond_wall_100 --terms 90 100 [--names T S Rx]
         [--xc 0.644700180 0.644700185 0.644700190] -o out.json
"""
import argparse
import json
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
import da  # noqa: E402
from gate import series_div  # noqa: E402

TARGET = {"T": "tau_1", "S": "gamma_1", "mt1": "gamma_1+nu_par", "mt2": "gamma_1+2nu_par",
          "mx1": "gamma_1+nu_perp", "mx2": "gamma_1+2nu_perp", "Rt": "nu_par", "Rx": "nu_perp"}


def all_series(ser, T):
    out = {k: ser[k][:T + 1] for k in ("T", "S", "mt1", "mt2", "mx1", "mx2")}
    # mt*, mx* start at p^1: divide by p before forming ratios
    out["Rt"] = series_div(ser["mt2"][1:T + 1], ser["mt1"][1:T + 1], T)
    out["Rx"] = series_div(ser["mx2"][1:T + 1], ser["mx1"][1:T + 1], T)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("series_json")
    ap.add_argument("key")
    ap.add_argument("--terms", type=int, nargs="+", required=True)
    ap.add_argument("--names", nargs="*", default=list(TARGET))
    ap.add_argument("--K", type=int, nargs="*", default=[2, 3])
    ap.add_argument("--Lmax", type=int, default=10)
    ap.add_argument("--Lstep", type=int, default=2)
    ap.add_argument("--span", type=int, default=10)
    ap.add_argument("--mstep", type=int, default=1)
    ap.add_argument("--max-per", type=int, default=None)
    ap.add_argument("--guess", type=float, default=0.6447)
    ap.add_argument("--xc", nargs="*", default=None, help="biased mode: critical points to impose")
    ap.add_argument("-o")
    a = ap.parse_args()
    ser = json.load(open(a.series_json))[a.key]
    Ls = range(0, a.Lmax + 1, a.Lstep)
    results = []
    for T in a.terms:
        allser = all_series(ser, T)
        for name in a.names:
            s = allser[name]
            Tn = len(s)
            Ms = range(Tn - a.span, Tn + 1, a.mstep)
            for K in a.K:
                t0 = time.time()
                if a.xc is None:
                    ests = da.ensemble(s, K, Ls, Ms, a.guess, max_per=a.max_per, defect_radius=0.002)
                    summ = da.summarize(ests, f"{name} K={K}")
                    if not summ:
                        continue
                    summ.update({"series": name, "K": K, "terms": T, "quantity": TARGET[name],
                                 "seconds": round(time.time() - t0, 1)})
                    results.append(summ)
                    print(json.dumps(summ), flush=True)
                else:
                    for xc in a.xc:
                        v = da.biased_ensemble(s, K, Ls, Ms, Fraction(xc))
                        if not v:
                            continue
                        med = float(np.median(v))
                        spr = float((np.quantile(v, 0.84) - np.quantile(v, 0.16)) / 2)
                        r = {"series": name, "K": K, "terms": T, "quantity": TARGET[name], "xc": float(xc),
                             "exponent": med, "spread": spr, "count": len(v),
                             "seconds": round(time.time() - t0, 1)}
                        results.append(r)
                        print(json.dumps(r), flush=True)
            if a.o:
                json.dump(results, open(a.o, "w"), indent=1)


if __name__ == "__main__":
    main()
