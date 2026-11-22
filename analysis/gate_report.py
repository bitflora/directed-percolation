"""Summarize gate.py output for the M2 precision-gain gate.

For every series and number of terms T, estimates from K = 2 and K = 3 are
pooled (median of the two medians). Two error proxies are reported:

* spread: half the 16-84% range of the approximant ensemble at fixed T, and
* drift: the change of the central estimate since the previous T.

The approximant spread at fixed T understates the true error when the
estimates still drift with T (as they do for this problem), so the gate uses
the drift. log10(drift per term) is fitted linearly against T and projected.

Usage: python gate_report.py results/analysis/bond146_gate.json [...] -o report.md
"""
import argparse
import json
import math
from collections import defaultdict

# Jensen 1999 (S01) bond estimates
J99 = {
    "xc": (0.644700185, 5e-9),
    "gamma": (2.277730, 5e-6),
    "nu_par": (1.733847, 6e-6),
    "2nu_perp": (2.193708, 4e-6),
    "gamma+nu_par": (4.01156, 1e-5),
    "gamma+2nu_par": (5.74539, 1e-5),
    "gamma+2nu_perp": (4.471425, 1.5e-5),
}


def fit_line(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return my - b * mx, b


def pooled(rows):
    """{(series, T): dict(xc, exp, spread_x, spread_e, quantity, count)} pooling K."""
    by = defaultdict(list)
    for r in rows:
        by[(r["series"], r["terms"])].append(r)
    out = {}
    for k, rs in by.items():
        xs = sorted(r["xc"] for r in rs)
        es = sorted(r["exponent"] for r in rs)
        mid = lambda v: v[len(v) // 2] if len(v) % 2 else (v[len(v) // 2 - 1] + v[len(v) // 2]) / 2
        out[k] = {
            "xc": mid(xs), "exp": mid(es),
            "sx": max(r["xc_spread"] for r in rs), "se": max(r["exponent_spread"] for r in rs),
            "kdiff_x": max(xs) - min(xs), "kdiff_e": max(es) - min(es),
            "quantity": rs[0]["quantity"], "count": sum(r["count"] for r in rs),
        }
    return out


def main():
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("-o")
    a = ap.parse_args()
    rows = []
    for p in a.inputs:
        rows += json.load(open(p))
    P = pooled(rows)
    names = ["S", "mu01", "mu02", "mu20", "R1", "R2"]
    Ts = sorted({t for (_, t) in P})
    out = ["## Central estimates against the number of terms (K = 2 and 3 pooled)", ""]
    out.append("| Series | Quantity | T | approximants | x_c | exponent | spread | K2–K3 difference | drift since previous T | exponent − J99 |")
    out.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    fits = {}
    for name in names:
        prev = None
        drifts = []
        for T in Ts:
            r = P.get((name, T))
            if not r:
                continue
            q = r["quantity"]
            drift = abs(r["exp"] - prev["exp"]) / (T - prev["T"]) if prev else None
            if drift:
                drifts.append((T, drift))
            dj = f"{r['exp'] - J99[q][0]:+.1e}" if q in J99 else ""
            out.append(
                f"| {name} | {q} | {T} | {r['count']} | {r['xc']:.10f} | {r['exp']:.7f} | {r['se']:.1e} | "
                f"{r['kdiff_e']:.1e} | {'' if drift is None else f'{drift:.1e}/term'} | {dj} |"
            )
            prev = dict(r, T=T)
        if len(drifts) >= 2:
            f = fit_line([t for t, _ in drifts], [math.log10(d) for _, d in drifts])
            fits[name] = (f, drifts)
    out.append("")
    out.append("## Drift fit and projection")
    out.append("")
    out.append("log10(drift per term) = a + b·T. If the drift keeps decaying geometrically, the remaining error after T terms is about drift(T)/(1 − 10^b) per term summed to infinity, i.e. drift(T)·10^b/(1 − 10^b).")
    out.append("")
    out.append("| Series | b (per term) | error factor per 10 terms | projected remaining error at T = 146 / 171 / 200 / 230 | gain 171 → 230 |")
    out.append("| --- | --- | --- | --- | --- |")
    for name, (f, drifts) in fits.items():
        a0, b = f
        if b >= 0:
            out.append(f"| {name} | {b:+.4f} | not decaying | — | — |")
            continue
        r = 10 ** b
        rem = lambda T: 10 ** (a0 + b * T) * r / (1 - r)
        proj = " / ".join(f"{rem(T):.1e}" for T in (146, 171, 200, 230))
        out.append(f"| {name} | {b:+.4f} | {10 ** (-10 * b):.2f}× | {proj} | {rem(171) / rem(230):.1f}× |")
    text = "\n".join(out)
    print(text)
    if a.o:
        open(a.o, "w", encoding="utf-8").write(text + "\n")


if __name__ == "__main__":
    main()
