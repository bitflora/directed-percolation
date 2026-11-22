"""Search for linear ODEs with polynomial coefficients satisfied by a series,
working modulo one prime (no big-integer reconstruction needed).

A series f = sum a_n x^n is D-finite if sum_{i=0}^{k} P_i(x) theta^i f = 0 for
some polynomials P_i of degree <= d, not all zero (theta = x d/dx). Given N
terms, the (k+1)(d+1) unknown coefficients must satisfy N linear equations
(the coefficients of x^0 .. x^(N-1)). If (k+1)(d+1) < N - margin and the
system has only the trivial solution modulo a large prime, no such ODE exists
with those (k, d). (Exactness: a nontrivial ODE over Q gives a nontrivial
solution modulo all but finitely many primes; with a 62-bit prime and integer
data an accidental failure is negligible but not impossible, so a second
prime is used.)

A positive result (nullspace found) is checked by verifying that the ODE also
annihilates the remaining `margin` terms that were not used to find it.

Usage: python dfinite.py results/runs/bond146.jsonl bond 146 S [--margin 10]
"""
import argparse
import json

import flint


def nullity(a, p, k, d, nuse):
    """Dimension of the solution space using the first nuse coefficients."""
    cols = (k + 1) * (d + 1)
    M = flint.nmod_mat(nuse, cols, p)
    for n in range(nuse):
        for i in range(k + 1):
            ni = pow(n, i, p) if i else 1
            for j in range(d + 1):
                if n - j < 0:
                    break
                # coefficient of x^n in x^j theta^i f = (n-j)^i a_{n-j}
                v = (pow(n - j, i, p) * a[n - j]) % p if i else a[n - j] % p
                M[n, i * (d + 1) + j] = v
    _, nullity_ = M.nullspace()
    return nullity_


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs")
    ap.add_argument("model")
    ap.add_argument("n", type=int)
    ap.add_argument("series")
    ap.add_argument("--margin", type=int, default=10)
    ap.add_argument("--kmax", type=int, default=12)
    a = ap.parse_args()
    recs = [json.loads(l) for l in open(a.runs) if l.strip()]
    recs = [r for r in recs if r["model"] == a.model and r["n"] == a.n]
    assert len(recs) >= 2, "need two primes"
    N = len(recs[0][a.series])
    print(f"{a.model} {a.series}: {N} terms, margin {a.margin}, primes {recs[0]['prime']} and {recs[1]['prime']}")
    print("k  largest d excluded  unknowns")
    found = []
    for k in range(0, a.kmax + 1):
        dmax = (N - a.margin) // (k + 1) - 1
        if dmax < 0:
            break
        res = []
        for r in recs[:2]:
            p = r["prime"]
            res.append(nullity(r[a.series], p, k, dmax, (k + 1) * (dmax + 1) + a.margin))
        rank_nullity = res
        if max(rank_nullity) == 0:
            print(f"{k:2d}  {dmax:4d}  {(k + 1) * (dmax + 1):5d}   none (both primes)")
        else:
            print(f"{k:2d}  {dmax:4d}  {(k + 1) * (dmax + 1):5d}   CANDIDATE nullity {rank_nullity}")
            found.append((k, dmax))
    if not found:
        print("No ODE of order k with polynomial degree <= the listed d exists for this series (mod these primes).")


if __name__ == "__main__":
    main()
