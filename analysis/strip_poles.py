"""Exploration for RESEARCH_PLAN M6 (haruspicy): poles of the mean cluster size
of bond directed percolation on a ring of L sites.

On a ring (cylinder) of circumference L, S_L(p) = sum_t E[number of occupied
sites at time t] is a rational function of p, because the Domany-Kinzel
transfer matrix T_L(p) is finite: S_L = u^T (I - T_L)^-1 n. Its coefficients
agree with the bulk series up to order about L. The haruspicy method of
Rechnitzer (S16) proves non-D-finiteness by showing that the poles of such
width-restricted rational functions accumulate on a curve. This script
computes the poles exactly for small L:

* states are nonempty subsets of Z_L up to rotation (the observable and the
  dynamics are rotation invariant);
* det(I - T_L(p)) is computed exactly: evaluated modulo several 62-bit primes
  at enough points, interpolated, and combined by CRT;
* its zeros are found with python-flint's rigorous complex root finder.

Output: results/analysis/strip_poles_L{L}.json and a plot.
Usage: python strip_poles.py 3 4 5 6 7
"""
import json
import sys

import flint


def necklaces(L):
    """Canonical (minimal-rotation) representatives of nonempty subsets of Z_L."""
    full = (1 << L) - 1
    rep = {}
    for m in range(1, 1 << L):
        r = min(((m << k) | (m >> (L - k))) & full for k in range(L))
        rep[m] = r
    reps = sorted(set(rep.values()))
    return reps, rep


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def transfer(L):
    """T[A_index][B_index] as integer polynomials in p (lists, low order first)."""
    reps, rep = necklaces(L)
    index = {r: i for i, r in enumerate(reps)}
    occ = {1: [0, 1], 2: [0, 2, -1]}      # 1 - (1-p)^k
    emp = {1: [1, -1], 2: [1, -2, 1]}     # (1-p)^k
    T = [dict() for _ in reps]
    for ia, A in enumerate(reps):
        # k_c = A[c-1] + A[c] (parents of site c are c-1 and c on the ring)
        ks = []
        for c in range(L):
            k = ((A >> c) & 1) + ((A >> ((c - 1) % L)) & 1)
            ks.append(k)
        cand = [c for c in range(L) if ks[c] > 0]
        for bits in range(1, 1 << len(cand)):
            poly = [1]
            B = 0
            for t, c in enumerate(cand):
                if (bits >> t) & 1:
                    poly = poly_mul(poly, occ[ks[c]])
                    B |= 1 << c
                else:
                    poly = poly_mul(poly, emp[ks[c]])
            ib = index[rep[B]]
            cur = T[ia].get(ib, [])
            if len(cur) < len(poly):
                cur = cur + [0] * (len(poly) - len(cur))
            for i, v in enumerate(poly):
                cur[i] += v
            T[ia][ib] = cur
    return reps, T


def det_poly(L, rank_one=False):
    """det(I - T_L(p)) as exact integer coefficients, by fraction-free (Bareiss)
    elimination over Z[p]. I - T(0) = I, so every pivot has constant term 1 and
    no pivoting is needed. With rank_one=True, returns det(I - T + n e_s^T),
    where n is the occupation count and e_s the singleton state; by the matrix
    determinant lemma, S_L = det(I - T + n e_s^T) / det(I - T) - 1."""
    reps, T = transfer(L)
    N = len(reps)
    M = [[flint.fmpz_poly([0]) for _ in range(N)] for _ in range(N)]
    for i in range(N):
        M[i][i] = flint.fmpz_poly([1])
        for jx, poly in T[i].items():
            M[i][jx] = M[i][jx] - flint.fmpz_poly(poly)
    if rank_one:
        js = reps.index(1)  # the singleton orbit
        for i in range(N):
            M[i][js] = M[i][js] + bin(reps[i]).count("1")
    prev = flint.fmpz_poly([1])
    for k in range(N - 1):
        piv = M[k][k]
        for i in range(k + 1, N):
            for j in range(k + 1, N):
                num = piv * M[i][j] - M[i][k] * M[k][j]
                q, r = divmod(num, prev)
                assert r == 0
                M[i][j] = q
        prev = piv
    det = M[N - 1][N - 1]
    return N, [int(c) for c in det.coeffs()]


def main():
    out = {}
    for L in map(int, sys.argv[1:]):
        N, coeffs = det_poly(L)
        P = flint.fmpz_poly(coeffs)
        # S_L = (det B - det A) / det A; cancel common factors to get true poles
        _, cb = det_poly(L, rank_one=True)
        num = flint.fmpz_poly(cb) - P
        g = P.gcd(num)
        P = P // g
        roots = []
        for r, mult in P.complex_roots():
            z = complex(r)
            roots.append([z.real, z.imag, mult])
        out[L] = {"states": N, "degree": P.degree(), "cancelled_degree": g.degree(), "roots": roots}
        near = sorted(roots, key=lambda r: abs(complex(r[0], r[1]) - 0.6447))[:3]
        print(f"L={L}: {N} orbit states, reduced denominator degree {P.degree()} (cancelled {g.degree()}), {len(roots)} distinct poles; "
              f"closest to p_c: {[f'{a:.4f}{b:+.4f}i' for a, b, _ in near]}", flush=True)
        json.dump(out, open("results/analysis/strip_poles_reduced.json", "w"))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.5, 6.5))
        for L, d in sorted(out.items()):
            xs = [r[0] for r in d["roots"]]
            ys = [r[1] for r in d["roots"]]
            ax.scatter(xs, ys, s=8, label=f"L = {L}")
        ax.scatter([0.644700185], [0], marker="x", s=80, color="k", label="p_c")
        for z in (-0.51670, complex(-0.22605, 0.4400), complex(-0.22605, -0.4400)):
            z = complex(z)
            ax.scatter([z.real], [z.imag], marker="+", s=80, color="r")
        ax.set_xlim(-1.2, 1.6)
        ax.set_ylim(-1.4, 1.4)
        ax.set_aspect("equal")
        ax.axhline(0, lw=0.5, color="0.6")
        ax.axvline(0, lw=0.5, color="0.6")
        ax.set_xlabel("Re p")
        ax.set_ylabel("Im p")
        ax.set_title("Poles of S_L(p) (common factors cancelled), bond DP on a ring of L sites\n"
                     "red +: bulk non-physical singularities (Jensen 1999)")
        ax.legend(loc="upper left", fontsize=8)
        fig.savefig("results/analysis/strip_poles_reduced.png", dpi=130, bbox_inches="tight")
        print("wrote results/analysis/strip_poles_reduced.png")
    except Exception as e:  # plotting is optional
        print("plot skipped:", e)


if __name__ == "__main__":
    main()
