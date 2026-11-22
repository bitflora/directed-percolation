"""Combine the order-166 (S, mu01, mu02) runs: the 31-bit primes first and the
15-bit u16 prime last, so that crt.py's leave-last-out check verifies the
reconstruction made from the six 31-bit primes."""
import json

recs = []
for f in ("results/runs/bond166.jsonl", "results/runs/bond166_r.jsonl"):
    try:
        recs += [json.loads(l) for l in open(f) if l.strip()]
    except FileNotFoundError:
        pass
seen = {}
for r in recs:
    seen[r["prime"]] = r
big = sorted((r for p, r in seen.items() if p >= 1 << 20), key=lambda r: -r["prime"])
small = [r for p, r in seen.items() if p < 1 << 20]
out = big + small
with open("results/runs/bond166_all.jsonl", "w") as f:
    for r in out:
        f.write(json.dumps(r) + "\n")
print(len(big), "31-bit primes,", len(small), "15-bit primes")
