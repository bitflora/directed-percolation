"""Combine per-prime engine output (JSON lines) into exact integer series.

Usage: python crt.py runs.jsonl [--check fixture.json] [-o series.json]

Every series is reconstructed with the Chinese remainder theorem as a signed
integer in (-M/2, M/2]. The reconstruction is accepted only if leaving out the
last prime gives the same integers (the last prime acts as a verification
prime), so at least two primes are required.
"""
import argparse
import json
from collections import defaultdict

NAMES = ("S", "mu01", "mu02", "mu20")  # bulk engine; other engines name their own series


def series_names(rec):
    return [k for k, v in rec.items() if isinstance(v, list)]


def crt_signed(residues, primes):
    M = 1
    x = 0
    for r, p in zip(residues, primes):
        # solve x' = x mod M, x' = r mod p
        t = ((r - x) * pow(M, -1, p)) % p
        x += M * t
        M *= p
    if x > M // 2:
        x -= M
    return x, M


def combine(records):
    primes = [r["prime"] for r in records]
    if len(set(primes)) != len(primes):
        raise SystemExit("duplicate primes in input")
    out = {}
    for name in series_names(records[0]):
        n1 = len(records[0][name])
        full, fewer = [], []
        for k in range(n1):
            res = [r[name][k] for r in records]
            full.append(crt_signed(res, primes)[0])
            fewer.append(crt_signed(res[:-1], primes[:-1])[0] if len(primes) > 1 else None)
        out[name] = full
        out[name + "_verified"] = len(primes) > 1 and full == fewer
    bits = sum(p.bit_length() for p in primes)
    return out, bits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs")
    ap.add_argument("--check")
    ap.add_argument("-o")
    a = ap.parse_args()
    groups = defaultdict(list)
    with open(a.runs) as f:
        for line in f:
            if line.strip():
                r = json.loads(line)
                groups[(r["model"], r["n"])].append(r)
    result = {}
    for (model, n), recs in sorted(groups.items()):
        ser, bits = combine(recs)
        ok = all(v for k, v in ser.items() if k.endswith("_verified"))
        maxbits = max(abs(x).bit_length() for k, v in ser.items() if not k.endswith("_verified") for x in v)
        print(f"{model} n={n}: {len(recs)} primes ({bits} bits), largest coefficient {maxbits} bits, verified={ok}")
        result[f"{model}_{n}"] = {"model": model, "order": n, "primes": len(recs), "verified": ok,
                                  **{k: v for k, v in ser.items() if not k.endswith('_verified')}}
    if a.check:
        fx = json.load(open(a.check))
        N = fx["order"]
        for key, r in result.items():
            ref = fx[r["model"]] if r["model"] in fx else fx  # wall fixtures are flat
            for name in series_names(ref):
                if name in r and name in ref:
                    m = min(N, r["order"]) + 1
                    assert r[name][:m] == ref[name][:m], (key, name, r[name][:m], ref[name][:m])
            print(f"{key}: matches fixture to order {min(N, r['order'])}")
    if a.o:
        with open(a.o, "w") as f:
            json.dump(result, f)
        print("wrote", a.o)


if __name__ == "__main__":
    main()
