"""M0 ground truth for directed percolation on the directed square lattice.

Three independent ways of computing the low-density series, in exact integer
arithmetic, cross-checked against each other:

1. Inclusion-exclusion over sets of directed paths (tiny t only). This checks
   the path-union sign rule d(g) = (-1)^(number of vertices with two incoming
   bonds) that the fast engine relies on.
2. The Domany-Kinzel transfer matrix: exact probabilities of each row
   configuration, swept one site at a time. Independent of the graph expansion.
3. The path-union expansion summed row by row, both directly and through the
   non-nodal reduction with S = 1/(1 - S_N) and the moment formulas.

Coordinates: row t has columns c = 0..t; site (t, c) sits at x = 2c - t and has
bonds to (t+1, c) (x-1) and (t+1, c+1) (x+1).

Output: JSON with the series S, mu01, mu02, mu20 for bond and site to order N,
used as a test fixture by the Rust engine.
"""
import itertools
import json
import sys
from functools import lru_cache


# ---------------------------------------------------------------- polynomials
def padd(a, b, N):
    out = list(a) + [0] * max(0, len(b) - len(a))
    for i, v in enumerate(b):
        out[i] += v
    return out[: N + 1]


def pmul(a, b, N):
    out = [0] * min(N + 1, len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0 or i > N:
            continue
        for j, y in enumerate(b):
            if i + j > N:
                break
            out[i + j] += x * y
    return out


def pscale(a, k):
    return [k * v for v in a]


def pinv(a, N):
    """Inverse of a power series with a[0] = 1 (or -1), truncated at N."""
    assert a[0] in (1, -1)
    inv = [0] * (N + 1)
    inv[0] = a[0]
    for n in range(1, N + 1):
        s = 0
        for k in range(1, min(n, len(a) - 1) + 1):
            s += a[k] * inv[n - k]
        inv[n] = -s * a[0]
    return inv


def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


# ------------------------------------------- 1. inclusion-exclusion over paths
def paths_to(x, t):
    """All directed paths from the origin to (x, t) as frozensets of bonds."""
    c_end = (x + t) // 2
    out = []
    for steps in itertools.product((0, 1), repeat=t):  # 0: col stays (x-1), 1: col+1 (x+1)
        if sum(steps) != c_end:
            continue
        c = 0
        bonds = []
        for u, s in enumerate(steps):
            bonds.append((u, c, s))
            c += s
        out.append(frozenset(bonds))
    return out


def pair_connectedness_ie(x, t, site=False):
    """C(x,t) by inclusion-exclusion over nonempty sets of paths (exact polynomial)."""
    paths = paths_to(x, t)
    poly = {}
    for r in range(1, len(paths) + 1):
        for subset in itertools.combinations(paths, r):
            union = frozenset().union(*subset)
            if site:
                sites = {(0, 0)}
                for (u, c, s) in union:
                    sites.add((u + 1, c + s))
                deg = len(sites) - 1  # origin is given occupied
            else:
                deg = len(union)
            poly[deg] = poly.get(deg, 0) + (-1) ** (r + 1)
    m = max(poly) if poly else 0
    return trim([poly.get(i, 0) for i in range(m + 1)])


# ------------------------------------------------- 2. Domany-Kinzel transfer matrix
def dk_rows(T, N, site=False):
    """Yield (t, {frozenset_of_cols: poly}) for t = 0..T, polys truncated at N.

    Bond DP: a site is occupied iff at least one occupied parent has its bond to
    it open, P = 1 - (1-p)^k.  Site DP: occupied with probability p iff at least
    one parent is occupied.  The origin is occupied with probability 1.
    """
    # (1-p)^k for k = 1, 2
    one_minus_p = {1: [1, -1], 2: [1, -2, 1]}
    occ = {k: padd([1], pscale(one_minus_p[k], -1), N) for k in (1, 2)}  # 1-(1-p)^k
    if site:
        occ = {1: [0, 1], 2: [0, 1]}
        emp = {1: [1, -1], 2: [1, -1]}
    else:
        emp = one_minus_p
    rows = {frozenset([0]): [1]}
    yield 0, rows
    for t in range(1, T + 1):
        new = {}
        for A, w in rows.items():
            # sites of the new row with at least one occupied parent
            cand = sorted({c for a in A for c in (a, a + 1)})
            ks = [(a_ in A) + ((a_ - 1) in A) for a_ in cand]
            for choice in itertools.product((0, 1), repeat=len(cand)):
                if not any(choice):
                    continue  # empty row: the cluster died; contributes nothing later
                poly = w
                for ch, k in zip(choice, ks):
                    poly = pmul(poly, occ[k] if ch else emp[k], N)
                    if not any(poly):
                        break
                if not any(poly):
                    continue
                B = frozenset(c for c, ch in zip(cand, choice) if ch)
                new[B] = padd(new.get(B, []), poly, N)
        rows = new
        yield t, rows


def dk_series(N, site=False):
    S = [0] * (N + 1)
    m01 = [0] * (N + 1)
    m02 = [0] * (N + 1)
    m20 = [0] * (N + 1)
    for t, rows in dk_rows(N, N, site):
        for B, w in rows.items():
            for c in B:
                x = 2 * c - t
                S = padd(S, w, N)
                m01 = padd(m01, pscale(w, t), N)
                m02 = padd(m02, pscale(w, t * t), N)
                m20 = padd(m20, pscale(w, x * x), N)
    return {"S": S, "mu01": m01, "mu02": m02, "mu20": m20}


# -------------------------------------------------- 3. path-union row expansion
def path_union_rows(T, N, site=False, nonnodal=False):
    """Yield (t, {frozenset_of_cols: poly}) for the path-union graph sum.

    Each occupied vertex of row t-1 chooses a nonempty subset of its two
    outgoing bonds.  Weight p per bond (bond) or p per newly occupied vertex
    (site), times -1 for every vertex receiving two bonds.  With nonnodal=True,
    singleton rows at 0 < t are removed after being reported.
    """
    rows = {frozenset([0]): [1]}
    for t in range(1, T + 1):
        new = {}
        for A, w in rows.items():
            A = sorted(A)
            for choice in itertools.product((1, 2, 3), repeat=len(A)):  # 1:L 2:R 3:LR
                indeg = {}
                nb = 0
                for a, ch in zip(A, choice):
                    if ch & 1:
                        indeg[a] = indeg.get(a, 0) + 1
                        nb += 1
                    if ch & 2:
                        indeg[a + 1] = indeg.get(a + 1, 0) + 1
                        nb += 1
                deg = len(indeg) if site else nb
                if deg > N:
                    continue
                sign = (-1) ** sum(1 for v in indeg.values() if v == 2)
                poly = [0] * deg + [sign]
                B = frozenset(indeg)
                new[B] = padd(new.get(B, []), pmul(w, poly, N), N)
        new = {B: w for B, w in new.items() if any(w)}
        yield t, new
        if nonnodal:
            new = {B: w for B, w in new.items() if len(B) > 1}
        rows = new


def pu_direct_series(N, site=False):
    S = [1] + [0] * N
    m01 = [0] * (N + 1)
    m02 = [0] * (N + 1)
    m20 = [0] * (N + 1)
    for t, rows in path_union_rows(N, N, site):
        for B, w in rows.items():
            if len(B) != 1:
                continue
            (c,) = B
            x = 2 * c - t
            S = padd(S, w, N)
            m01 = padd(m01, pscale(w, t), N)
            m02 = padd(m02, pscale(w, t * t), N)
            m20 = padd(m20, pscale(w, x * x), N)
    return {"S": S, "mu01": m01, "mu02": m02, "mu20": m20}


def assemble(SN, N01, N02, N20, N):
    """Full moments from non-nodal sums (Jensen 1999, eqs. 6-9)."""
    one_minus = padd([1], pscale(SN, -1), N)
    S = pinv(one_minus, N)
    S2 = pmul(S, S, N)
    mu01 = pmul(N01, S2, N)
    mu02 = pmul(padd(N02, pscale(pmul(pmul(N01, N01, N), S, N), 2), N), S2, N)
    mu20 = pmul(N20, S2, N)
    return {"S": S, "mu01": mu01, "mu02": mu02, "mu20": mu20}


def nonnodal_series(N, site=False):
    # non-nodal pieces at row t have >= 2t bonds (bond) or >= 2t-1 new sites (site)
    tmax = max(1, (N + 1) // 2 if site else N // 2)
    SN = [0] * (N + 1)
    N01 = [0] * (N + 1)
    N02 = [0] * (N + 1)
    N20 = [0] * (N + 1)
    for t, rows in path_union_rows(tmax, N, site, nonnodal=True):
        for B, w in rows.items():
            if len(B) != 1:
                continue
            (c,) = B
            x = 2 * c - t
            SN = padd(SN, w, N)
            N01 = padd(N01, pscale(w, t), N)
            N02 = padd(N02, pscale(w, t * t), N)
            N20 = padd(N20, pscale(w, x * x), N)
    return assemble(SN, N01, N02, N20, N)


def pad(d, N):
    return {k: (v + [0] * (N + 1 - len(v)))[: N + 1] for k, v in d.items()}


def main():
    # 1. sign rule: inclusion-exclusion vs path-union expansion, full polynomials
    for site in (False, True):
        big = 10 ** 6
        pu = dict(path_union_rows(4, big, site))
        for t in range(1, 5):
            for c in range(t + 1):
                x = 2 * c - t
                ie = pair_connectedness_ie(x, t, site)
                got = trim(pu[t].get(frozenset([c]), []))
                # path-union rows also hold multi-vertex rows; C(x,t) is the sum
                # over all graphs ENDING at (x,t), i.e. singleton row states.
                assert ie == got, (site, x, t, ie, got)
    print("inclusion-exclusion == path-union sign rule for t <= 4 (bond, site)")

    N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    out = {"order": N}
    for name, site in (("bond", False), ("site", True)):
        dk = pad(dk_series(N, site), N)
        pud = pad(pu_direct_series(N, site), N)
        nn = pad(nonnodal_series(N, site), N)
        for k in dk:
            assert dk[k] == pud[k], (name, k, dk[k], pud[k])
            assert dk[k] == nn[k], (name, k, dk[k], nn[k])
        print(f"{name}: DK == path-union direct == non-nodal to order {N}")
        print("  S    =", dk["S"])
        out[name] = dk
    path = sys.argv[2] if len(sys.argv) > 2 else None
    if path:
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print("wrote", path)


if __name__ == "__main__":
    main()
