"""Complete the order-171 bond series with mu20.

Input: results/runs/bond171_x.jsonl, all-moment runs (15-bit primes), and the
verified S, mu01, mu02 series in results/series/bond171.json.

* Every per-prime run must reproduce S, mu01 and mu02 of the verified series
  modulo its prime: an independent check of each run.
* mu20 is reconstructed by CRT from all primes but the last and accepted only
  if adding the last prime changes nothing (crt.py's leave-last-out rule).

Usage: python combine_bond171_mu20.py [-o results/series/bond171.json]
"""
import argparse
import json
import sys

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
from crt import crt_signed  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="results/runs/bond171_x.jsonl")
    ap.add_argument("--series", default="results/series/bond171.json")
    ap.add_argument("-o", help="write the completed series here")
    a = ap.parse_args()
    ser_all = json.load(open(a.series))
    key = next(iter(ser_all))
    ser = ser_all[key]
    N = ser["order"]
    recs = {}
    for line in open(a.runs):
        if line.strip():
            r = json.loads(line)
            assert r["model"] == "bond" and r["n"] == N and r["moments"] == 3
            recs[r["prime"]] = r  # a rerun of a prime replaces the earlier line
    recs = list(recs.values())
    primes = [r["prime"] for r in recs]
    for r in recs:
        p = r["prime"]
        for name in ("S", "mu01", "mu02"):
            exp = [x % p for x in map(int, ser[name])]
            if r[name] != exp:
                bad = next(k for k in range(len(exp)) if r[name][k] != exp[k])
                raise SystemExit(f"prime {p}: {name} differs from the verified series at order {bad}")
    print(f"{len(recs)} runs; every run reproduces the verified S, mu01, mu02 modulo its prime")
    bits = sum(p.bit_length() - 1 for p in primes)
    full, fewer = [], []
    for k in range(N + 1):
        res = [r["mu20"][k] for r in recs]
        full.append(crt_signed(res, primes)[0])
        fewer.append(crt_signed(res[:-1], primes[:-1])[0])
    ok = len(primes) > 1 and full == fewer
    maxbits = max(abs(x).bit_length() for x in full)
    print(f"mu20: {len(primes)} primes (product above 2^{bits}), largest coefficient {maxbits} bits, verified={ok}")
    old = None
    for f in ("results/series/bond156.json",):
        try:
            o = json.load(open(f))
            old = next(iter(o.values()))
        except FileNotFoundError:
            pass
    if old is not None and "mu20" in old:
        m = old["order"] + 1
        pre = [int(x) for x in old["mu20"][:m]] == full[:m]
        print(f"mu20 prefix matches the verified order-{old['order']} series: {pre}")
        ok = ok and pre
    if not ok:
        raise SystemExit("mu20 not verified; not writing")
    if a.o:
        ser["mu20"] = full
        ser["mu20_primes"] = len(primes)
        ser["verified"] = ser.get("verified", True) and ok
        json.dump(ser_all, open(a.o, "w"))
        print("wrote", a.o)


if __name__ == "__main__":
    main()
