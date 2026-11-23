"""Ground truth for bond directed percolation next to an impenetrable wall.

Geometry (Essam et al. 1996, Jensen 1999 II): directed square lattice, sites
(t, x) with t + x even, bonds (t, x) -> (t+1, x-1) and (t+1, x+1). The source
is at (0, 0) and bonds leading to sites with x < 0 are always absent.

Series (coefficients of p^n, exact integers):
  S     = sum C(t,x)                 mean cluster size,        exponent gamma_1
  mx1   = sum x C(t,x)               first x moment            gamma_1 + nu_perp
  mx2   = sum x^2 C(t,x)                                       gamma_1 + 2 nu_perp
  mt1   = sum t C(t,x)                                         gamma_1 + nu_par
  mt2   = sum t^2 C(t,x)                                       gamma_1 + 2 nu_par
  T     = sum_{t>=1} P(t)            mean cluster length,      exponent tau_1
          (P(t): probability that the cluster reaches row t)

Two independent methods:
1. Domany-Kinzel transfer matrix: exact probabilities of each wet row.
2. Path-union expansion. C(t,x) is the sum over unions of directed paths ending
   at (t,x) of (-1)^(joins) p^(bonds); P(t) is the same sum over unions of
   paths ending anywhere on row t, with an extra sign (-1)^(k-1) for k
   endpoints (inclusion-exclusion over the set of wet sites of row t, i.e. the
   d-weight with a super-sink joined to row t).

Usage: python wall_ground_truth.py N [--check] [-o fixture.json]
  --check compares with the table of Essam, Guttmann, Jensen, TanlaKishani,
  J. Phys. A 29 (1996) 1619, rows 0..N.
"""
import argparse
import itertools
import json

# Essam et al. 1996, Table 1, rows 0..30: T, S, mt1, mt2, mx1, mx2
# (their mu_{m,n} = sum t^m x^n C).
EGJT96 = [
    (0, 1, 0, 0, 0, 0), (1, 1, 1, 1, 1, 1), (2, 2, 4, 8, 2, 4), (2, 3, 9, 27, 5, 11),
    (5, 6, 24, 96, 10, 28), (5, 9, 47, 241, 21, 65), (11, 17, 108, 672, 40, 144),
    (13, 26, 201, 1499, 77, 303), (28, 47, 424, 3676, 142, 624),
    (25, 72, 762, 7644, 262, 1240), (75, 129, 1538, 17398, 470, 2438),
    (56, 194, 2675, 34369, 843, 4661), (188, 348, 5258, 74512, 1486, 8872),
    (112, 516, 8915, 141615, 2609, 16487), (458, 929, 17233, 296939, 4529, 30635),
    (319, 1351, 28518, 546394, 7846, 55734), (1157, 2456, 54636, 1119562, 13448, 101618),
    (312, 3506, 88459, 2004015, 23027, 181751), (3389, 6471, 169004, 4043156, 39096, 326608),
    (562, 8929, 266670, 7047626, 66320, 575790), (9193, 17029, 512651, 14102481, 111795, 1022909),
    (-2419, 22579, 786932, 23956166, 187946, 1781314),
    (24689, 44707, 1530464, 47809422, 314844, 3135130),
    (-6090, 55969, 2270857, 79011279, 526367, 5402999),
    (83997, 117836, 4516598, 158359672, 876362, 9435440),
    (-80845, 137313, 6439085, 254037643, 1455579, 16106911),
    (219791, 311654, 13207919, 514524887, 2415059, 27970523),
    (-95543, 324989, 17852082, 796972392, 3989542, 47305236),
    (653560, 833496, 38438680, 1646320650, 6597538, 81807186),
    (-1015961, 756309, 48640815, 2447308375, 10834513, 137158135),
    (2302634, 2242031, 111440275, 5201705453, 17869253, 236510661),
]
NAMES = ("T", "S", "mt1", "mt2", "mx1", "mx2")


def padd_into(acc, a, scale=1):
    for i, v in enumerate(a):
        if i < len(acc):
            acc[i] += scale * v


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


def empty(N):
    return {k: [0] * (N + 1) for k in NAMES}


def accumulate(ser, t, B, w):
    """Add the contribution of row t with wet set B (x values) and weight w."""
    for x in B:
        padd_into(ser["S"], w)
        padd_into(ser["mt1"], w, t)
        padd_into(ser["mt2"], w, t * t)
        padd_into(ser["mx1"], w, x)
        padd_into(ser["mx2"], w, x * x)


# ------------------------------------------------------ 1. Domany-Kinzel
def dk_series(N):
    """Rows up to t = N (P(t) = O(p^t)), probabilities truncated at p^N."""
    ser = empty(N)
    occ = {1: [0, 1], 2: [0, 2, -1]}          # 1 - (1-p)^k
    emp = {1: [1, -1], 2: [1, -2, 1]}         # (1-p)^k
    rows = {frozenset([0]): [1]}
    accumulate(ser, 0, [0], [1])
    for t in range(1, N + 1):
        new = {}
        for A, w in rows.items():
            cand = sorted({y for a in A for y in (a - 1, a + 1) if y >= 0})
            ks = [((y - 1) in A) + ((y + 1) in A) for y in cand]
            for choice in itertools.product((0, 1), repeat=len(cand)):
                if not any(choice):
                    continue
                poly = w
                for ch, k in zip(choice, ks):
                    poly = pmul(poly, occ[k] if ch else emp[k], N)
                    if not any(poly):
                        break
                if any(poly):
                    B = frozenset(y for y, ch in zip(cand, choice) if ch)
                    acc = new.setdefault(B, [0] * (N + 1))
                    padd_into(acc, poly)
        rows = {B: w for B, w in new.items() if any(w)}
        alive = [0] * (N + 1)
        for B, w in rows.items():
            padd_into(alive, w)
            accumulate(ser, t, B, w)
        padd_into(ser["T"], alive)
    return ser


# --------------------------------------------------- 2. path-union expansion
def pu_series(N):
    ser = empty(N)
    rows = {frozenset([0]): [1]}
    accumulate(ser, 0, [0], [1])
    for t in range(1, N + 1):
        new = {}
        for A, w in rows.items():
            A = sorted(A)
            opts = []
            for a in A:
                o = []
                for ch in (1, 2, 3):           # 1: x-1, 2: x+1, 3: both
                    ends = [y for bit, y in ((1, a - 1), (2, a + 1)) if ch & bit]
                    if all(y >= 0 for y in ends):
                        o.append(ends)
                opts.append(o)
            for choice in itertools.product(*opts):
                indeg = {}
                nb = 0
                for ends in choice:
                    for y in ends:
                        indeg[y] = indeg.get(y, 0) + 1
                        nb += 1
                sign = (-1) ** sum(1 for v in indeg.values() if v == 2)
                B = frozenset(indeg)
                acc = new.setdefault(B, [0] * (N + 1))
                for i, v in enumerate(w):
                    if v and i + nb <= N:
                        acc[i + nb] += sign * v
        rows = {B: w for B, w in new.items() if any(w)}
        alive = [0] * (N + 1)
        for B, w in rows.items():
            padd_into(alive, w, (-1) ** (len(B) - 1))
        padd_into(ser["T"], alive)
        # pair connectedness: C(t,x) = sum over graphs whose only endpoint is x
        for B, w in rows.items():
            if len(B) == 1:
                accumulate(ser, t, B, w)
    return ser


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("-o")
    a = ap.parse_args()
    dk = dk_series(a.N)
    pu = pu_series(a.N)
    for k in NAMES:
        same = dk[k] == pu[k]
        print(f"{k:4s} DK == path-union: {same}")
        if not same:
            bad = next(i for i in range(a.N + 1) if dk[k][i] != pu[k][i])
            print("   first difference at order", bad, dk[k][bad], pu[k][bad])
    if a.check:
        for j, k in enumerate(NAMES):
            ref = [row[j] for row in EGJT96[: a.N + 1]]
            ok = dk[k][: len(ref)] == ref
            print(f"{k:4s} matches Essam et al. 1996 rows 0..{len(ref) - 1}: {ok}")
            if not ok:
                bad = next(i for i in range(len(ref)) if dk[k][i] != ref[i])
                print("   first difference at order", bad, dk[k][bad], ref[bad])
    if a.o:
        json.dump({"model": "bond_wall", "order": a.N, **dk}, open(a.o, "w"))
        print("wrote", a.o)


if __name__ == "__main__":
    main()
