"""Summarize gate.py output: estimates against the number of terms, comparison
with Jensen 1999, and a fit of the approximant spread against the number of
terms, projected to longer series.

Usage: python gate_report.py results/analysis/bond146.json [more.json ...] -o report.md
"""
import argparse
import json
import math
from collections import defaultdict

# Jensen 1999 (S01) bond estimates, eq. (15) and section 5
J99 = {
    "xc": (0.644700185, 5e-9),
    "gamma": (2.277730, 5e-6),
    "nu_par": (1.733847, 6e-6),
    "2nu_perp": (2.193708, 4e-6),
    "gamma+nu_par": (4.01156, 1e-5),
    "gamma+2nu_par": (5.74539, 1e-5),
    "gamma+2nu_perp": (4.471425, 1.5e-5),
}
J99_SITE = {"xc": (0.70548522, 4e-8)}


def fit_line(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return my - b * mx, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("-o")
    ap.add_argument("--model", default="bond")
    a = ap.parse_args()
    rows = []
    for p in a.inputs:
        rows += json.load(open(p))
    ref = J99 if a.model == "bond" else J99_SITE
    by = defaultdict(list)
    for r in rows:
        by[(r["series"], r["K"])].append(r)
    out = []
    out.append("| Series | Quantity | K | Terms | Approximants | x_c | spread | Exponent | spread | Exponent − J99 |")
    out.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for (name, K), rs in sorted(by.items()):
        for r in sorted(rs, key=lambda r: r["terms"]):
            q = r["quantity"]
            dref = ""
            if q in ref:
                dref = f"{r['exponent'] - ref[q][0]:+.2e}"
            out.append(
                f"| {name} | {q} | {K} | {r['terms']} | {r['count']} | {r['xc']:.10f} | {r['xc_spread']:.1e} | "
                f"{r['exponent']:.7f} | {r['exponent_spread']:.1e} | {dref} |"
            )
    out.append("")
    # spread against terms: fit log10(spread) = a + b T, per series, pooled over K
    out.append("Fits of log10(spread) against number of terms T, pooled over K (projection only):")
    out.append("")
    out.append("| Series | Quantity | slope per 10 terms (x_c) | slope per 10 terms (exponent) | projected exponent spread at 171 / 200 / 230 |")
    out.append("| --- | --- | --- | --- | --- |")
    names = sorted({r["series"] for r in rows})
    for name in names:
        rs = [r for r in rows if r["series"] == name and r["exponent_spread"] > 0 and r["xc_spread"] > 0]
        if len({r["terms"] for r in rs}) < 2:
            continue
        T = [r["terms"] for r in rs]
        fx = fit_line(T, [math.log10(r["xc_spread"]) for r in rs])
        fe = fit_line(T, [math.log10(r["exponent_spread"]) for r in rs])
        proj = " / ".join(f"{10 ** (fe[0] + fe[1] * t):.1e}" for t in (171, 200, 230)) if fe else ""
        out.append(
            f"| {name} | {rs[0]['quantity']} | {10 * fx[1]:+.2f} | {10 * fe[1]:+.2f} | {proj} |"
        )
    text = "\n".join(out)
    print(text)
    if a.o:
        open(a.o, "w", encoding="utf-8").write(text + "\n")


if __name__ == "__main__":
    main()
