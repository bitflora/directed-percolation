"""Differential-approximant analysis of the DP moment series, and the M2
precision-gain study (estimates and their spread as a function of the number
of series terms used).

Usage:
  python gate.py series.json bond_146 --terms 90 110 130 146 -o results/analysis/bond146.json

Series analysed (Jensen 1999, section 4.1):
  S                   ~ (pc - p)^-gamma
  mu01                ~ (pc - p)^-(gamma + nu_par)
  mu02                ~ (pc - p)^-(gamma + 2 nu_par)
  mu20                ~ (pc - p)^-(gamma + 2 nu_perp)
  R1 = mu02 / mu01    ~ (pc - p)^-nu_par
  R2 = mu20 mu02 / mu01^2 ~ (pc - p)^-(2 nu_perp)
"""
import argparse
import json
import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
import da  # noqa: E402


def series_div(a, b, n):
    """a / b as power series to n terms (Fractions); b[0] != 0."""
    out = []
    a = [Fraction(x) for x in a[:n]] + [Fraction(0)] * max(0, n - len(a))
    for k in range(n):
        s = a[k] - sum(out[i] * b[k - i] for i in range(max(0, k - len(b) + 1), k))
        out.append(s / b[0])
    return out


def series_mul(a, b, n):
    out = [0] * n
    for i, x in enumerate(a[:n]):
        if x == 0:
            continue
        for j, y in enumerate(b[: n - i]):
            out[i + j] += x * y
    return out


def derived(ser, n):
    """Ratio series. mu01 and mu02 start at p^1, so divide both by p first."""
    m01 = ser["mu01"][1:n + 1]
    m02 = ser["mu02"][1:n + 1]
    m20 = ser["mu20"][1:n + 1]
    r1 = series_div(m02, m01, n)
    r2 = series_div(series_mul(m20, m02, n), series_mul(m01, m01, n), n)
    return {"R1": r1, "R2": r2}


TARGET = {"S": "gamma", "mu01": "gamma+nu_par", "mu02": "gamma+2nu_par",
          "mu20": "gamma+2nu_perp", "R1": "nu_par", "R2": "2nu_perp"}


def analyse(ser, T, guess, Ks=(2, 3), Ls=range(0, 11, 2), span=10, names=None, mstep=1, max_per=None):
    out = []
    der = derived(ser, T)
    allser = {"S": ser["S"][:T + 1], "mu01": ser["mu01"][:T + 1], "mu02": ser["mu02"][:T + 1],
              "mu20": ser["mu20"][:T + 1], "R1": der["R1"][:T], "R2": der["R2"][:T]}
    for name in names or allser:
        s = allser[name]
        Tn = len(s)
        for K in Ks:
            t0 = time.time()
            ests = da.ensemble(s, K, Ls, range(Tn - span, Tn + 1, mstep), guess, max_per=max_per, defect_radius=0.002)
            summ = da.summarize(ests, f"{name} K={K}")
            if summ:
                summ.update({"series": name, "K": K, "terms": T, "quantity": TARGET[name],
                             "seconds": round(time.time() - t0, 1)})
                out.append(summ)
                print(json.dumps(summ), flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("series_json")
    ap.add_argument("key")
    ap.add_argument("--terms", type=int, nargs="+", required=True)
    ap.add_argument("--guess", type=float, default=None)
    ap.add_argument("--names", nargs="*")
    ap.add_argument("--K", type=int, nargs="*", default=[2, 3])
    ap.add_argument("--Lmax", type=int, default=10)
    ap.add_argument("--Lstep", type=int, default=2)
    ap.add_argument("--span", type=int, default=10)
    ap.add_argument("--mstep", type=int, default=1)
    ap.add_argument("--max-per", type=int, default=None, help="degree assignments per (M, L)")
    ap.add_argument("-o")
    a = ap.parse_args()
    ser = json.load(open(a.series_json))[a.key]
    guess = a.guess or (0.6447 if ser["model"] == "bond" else 0.70549)
    results = []
    for T in a.terms:
        results += analyse(ser, T, guess, Ks=a.K, names=a.names, Ls=range(0, a.Lmax + 1, a.Lstep),
                           span=a.span, mstep=a.mstep, max_per=a.max_per)
        if a.o:
            json.dump(results, open(a.o, "w"), indent=1)


if __name__ == "__main__":
    main()
