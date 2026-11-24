"""Inhomogeneous differential approximants for exact integer series.

A K-th order approximant to f(x) = sum a_n x^n is

    sum_{i=0}^{K} Q_i(x) (x d/dx)^i f(x) = P(x),

with deg Q_i = N_i, deg P = L, normalization Q_K(0) = 1, and the polynomial
coefficients fixed by matching the first M = sum(N_i + 1) + L coefficients
(the form used by Jensen 1999, eq. 14).

Singularities are zeros x_c of Q_K. If f ~ (x_c - x)^mu there, then
mu = K - 1 - Q_{K-1}(x_c) / (x_c Q_K'(x_c)), and the critical exponent of a
divergent quantity, f ~ (x_c - x)^(-lambda), is lambda = -mu.

The linear system is solved in arbitrary-precision ball arithmetic (python-flint
arb_mat), so ill-conditioning shows up as an error ball, not a silent wrong
answer. Approximants whose result is not determined to at least 30 digits are
dropped.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass

import flint
import numpy as np
from series_input import load_record

PREC_BITS = 400


@dataclass
class Estimate:
    M: int          # number of series terms used
    K: int
    L: int
    degs: tuple     # (N_0, ..., N_K)
    xc: float       # critical point (as float; xc_str has full precision)
    xc_str: str
    exponent: float  # lambda, with f ~ (xc - x)^(-lambda)
    others: list    # other zeros of Q_K (complex), as python complex


def _to_arb(v):
    if isinstance(v, int):
        return flint.arb(flint.fmpz(v))
    if hasattr(v, "numerator"):
        return flint.arb(flint.fmpq(v.numerator, v.denominator))
    return flint.arb(v)


def _theta_powers(a, K):
    """[n^i a_n] for i = 0..K, as arb (exact for ints and fractions up to PREC_BITS)."""
    flint.ctx.prec = PREC_BITS
    base = [_to_arb(v) for v in a]
    return [[base[n] * (n ** i) for n in range(len(a))] for i in range(K + 1)]


def approximant(a, K, L, degs, theta=None, bias=None):
    """Solve one DA. `a` are exact ints or fractions, len(a) >= M.
    With `bias` (a Fraction or arb), Q_K is forced to vanish there and one
    fewer series coefficient is matched. Returns (Q list, P) as arb lists or None."""
    flint.ctx.prec = PREC_BITS
    M = sum(d + 1 for d in degs) + L  # unknowns after fixing Q_K[0] = 1
    if M > len(a):
        return None
    if theta is None:
        theta = _theta_powers(a, K)
    # unknown ordering: Q_0[0..N_0], ..., Q_K[1..N_K], P[0..L]
    cols = []
    for i in range(K + 1):
        start = 1 if i == K else 0
        for k in range(start, degs[i] + 1):
            cols.append(("q", i, k))
    for k in range(L + 1):
        cols.append(("p", k))
    assert len(cols) == M, (len(cols), M)
    zero = flint.arb(0)
    rows = []
    nmatch = M - 1 if bias is not None else M
    for n in range(nmatch):
        row = []
        for col in cols:
            if col[0] == "q":
                _, i, k = col
                row.append(theta[i][n - k] if n - k >= 0 else zero)
            else:
                row.append(flint.arb(-1) if col[1] == n else zero)
        rows.append(row)
    rhs = [-theta[K][n] for n in range(nmatch)]
    if bias is not None:
        xb = _to_arb(bias)
        # Q_K(xb) = 1 + sum_{k>=1} q_{K,k} xb^k = 0
        row = []
        for col in cols:
            if col[0] == "q" and col[1] == K:
                row.append(xb ** col[2])
            else:
                row.append(zero)
        rows.append(row)
        rhs.append(flint.arb(-1))
    A = flint.arb_mat(rows)
    # the fixed Q_K[0] = 1 term goes to the right-hand side
    b = flint.arb_mat([[v] for v in rhs])
    try:
        sol = A.solve(b)
    except (ZeroDivisionError, ValueError):
        return None
    vals = [sol[c, 0] for c in range(M)]
    Q = [[None] * (degs[i] + 1) for i in range(K + 1)]
    Q[K][0] = flint.arb(1)
    P = [None] * (L + 1)
    for v, col in zip(vals, cols):
        if col[0] == "q":
            Q[col[1]][col[2]] = v
        else:
            P[col[1]] = v
    return Q, P


def _rel_ok(x, digits=30):
    if not x.is_finite():
        return False
    return x.rad() <= abs(x.mid()) * flint.arb(f"1e-{digits}")


def _refine_root(qpoly, r0):
    """Find a simple real root, then certify its enclosure by interval Newton."""
    dq = qpoly.derivative()
    x = flint.arb(r0)
    for _ in range(60):
        derivative = dq(x)
        if not derivative.is_finite() or derivative.contains(0):
            return None
        step = qpoly(x) / derivative
        if not step.is_finite():
            return None
        x = (x - step.mid()).mid()
        if abs(step.mid()) < flint.arb("1e-45"):
            break
    else:
        return None
    # Point iteration is only a seed: its midpoint alone carries no root
    # error bound. Strict interval-Newton inclusion certifies a unique root.
    interval = flint.arb(x, abs(x) * flint.arb("1e-32") + flint.arb("1e-60"))
    derivative = dq(interval)
    if not derivative.is_finite() or derivative.contains(0):
        return None
    root = x - qpoly(x) / derivative
    if not interval.contains_interior(root) or not _rel_ok(root):
        return None
    if not qpoly(root).contains(0):
        return None
    return root


def analyse_one(a, K, L, degs, guess, theta=None, window=0.05):
    res = approximant(a, K, L, degs, theta)
    if res is None:
        return None
    Q, P = res
    qk = Q[K]
    if not all(_rel_ok(c, 25) for c in qk):
        return None
    # rough roots with numpy (coefficients in increasing order)
    coeffs = [float(c.mid()) for c in qk]
    if not all(math.isfinite(c) for c in coeffs):
        return None
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
    if len(coeffs) < 2:
        return None
    try:
        roots = np.roots(coeffs[::-1])
    except np.linalg.LinAlgError:
        return None
    real_near = [r for r in roots if abs(r.imag) < 1e-6 and abs(r.real - guess) < window]
    if not real_near:
        return None
    r0 = min(real_near, key=lambda r: abs(r.real - guess)).real
    # Certify the root and propagate its uncertainty into the exponent.
    qpoly = flint.arb_poly(qk)
    dq = qpoly.derivative()
    x = _refine_root(qpoly, r0)
    if x is None or abs(float(x.mid()) - guess) >= window or abs(float(x.mid()) - r0) > 1e-6:
        return None
    qkm1 = flint.arb_poly(Q[K - 1])
    mu = K - 1 - qkm1(x) / (x * dq(x))
    lam = -mu
    if not _rel_ok(lam) or not math.isfinite(float(lam.mid())):
        return None
    others = [complex(r) for r in roots if abs(r - r0) > 1e-8]
    return Estimate(M=sum(d + 1 for d in degs) + L, K=K, L=L, degs=tuple(degs),
                    xc=float(x.mid()), xc_str=x.mid().str(40, radius=False),
                    exponent=float(lam.mid()), others=others)


def _biased_value(res, K, xc, digits=10):
    if res is None:
        return None
    Q, P = res
    x = _to_arb(xc)
    qk = flint.arb_poly(Q[K])
    dq = qk.derivative()
    mu = K - 1 - flint.arb_poly(Q[K - 1])(x) / (x * dq(x))
    lam = -mu
    if not _rel_ok(lam, digits):
        return None
    return float(lam.mid())


def biased_exponent(a, K, L, degs, xc, theta=None):
    """Exponent at a fixed critical point xc (Fraction), from a biased DA."""
    return _biased_value(approximant(a, K, L, degs, theta, bias=xc), K, xc)


def biased_ensemble(a, K, Ls, Ms, xc):
    theta = _theta_powers(a, K)
    out = []
    for M in Ms:
        for L in Ls:
            for degs in degree_sets(M, K, L):
                v = biased_exponent(a, K, L, degs, xc, theta)
                if v is not None and abs(v) < 50:
                    out.append(v)
    return out


def _prepare_biased(a, K, L, degs, theta):
    """Solve A*x=b and A*v=e_last together, before replacing A's last row.

    All but the last coefficient equation are shared by a biased family.
    The unbiased final equation supplies a convenient nonsingular base; it
    need not exist for every family, so callers retain the direct solver.
    """
    M = sum(d + 1 for d in degs) + L
    if M > len(a):
        return None
    cols = [("q", i, j) for i in range(K+1)
            for j in range(1 if i == K else 0, degs[i]+1)]
    cols += [("p", j) for j in range(L+1)]
    zero = flint.arb(0)
    rows = [[(theta[c[1]][n-c[2]] if n >= c[2] else zero) if c[0] == "q"
             else flint.arb(-1 if c[1] == n else 0) for c in cols] for n in range(M)]
    rhs = flint.arb_mat([[-theta[K][n], flint.arb(int(n == M-1))] for n in range(M)])
    try:
        sol = flint.arb_mat(rows).solve(rhs)
    except (ZeroDivisionError, ValueError):
        return None
    return cols, [sol[j, 0] for j in range(M)], [sol[j, 1] for j in range(M)]


def _apply_bias(basis, K, L, degs, xc):
    """Replace the last equation with Q_K(xc)=0 using ball arithmetic."""
    cols, base, direction = basis
    x = _to_arb(xc)
    weights = [(j, x ** c[2]) for j, c in enumerate(cols) if c[0] == "q" and c[1] == K]
    denom = sum((direction[j] * w for j, w in weights), flint.arb(0))
    if not denom.is_finite() or denom.contains(0):
        return None
    alpha = (-1 - sum(base[j] * w for j, w in weights)) / denom
    values = [b + d * alpha for b, d in zip(base, direction)]
    Q = [[None] * (d+1) for d in degs]
    Q[K][0] = flint.arb(1)
    P = [None] * (L+1)
    for value, c in zip(values, cols):
        if c[0] == "q":
            Q[c[1]][c[2]] = value
        else:
            P[c[1]] = value
    return Q, P


def biased_ensembles(a, K, Ls, Ms, xcs):
    """One ensemble per bias, sharing a solve for each degree assignment.

    Keep degree order and filtering identical to biased_ensemble. Updates
    must determine the exponent to 30 digits; singular or poorly resolved
    families fall back to the original biased solve and its accuracy gate.
    Only one family's basis is retained at a time.
    """
    xcs = list(xcs)
    if not xcs:
        return []
    if len(xcs) == 1:
        return [biased_ensemble(a, K, Ls, Ms, xcs[0])]
    theta = _theta_powers(a, K)
    Ls = tuple(Ls)
    out = [[] for _ in xcs]
    for M in Ms:
        for L in Ls:
            for degs in degree_sets(M, K, L):
                basis = _prepare_biased(a, K, L, degs, theta)
                for xc, values in zip(xcs, out):
                    v = None
                    if basis is not None:
                        v = _biased_value(_apply_bias(basis, K, L, degs, xc), K, xc, digits=30)
                    if v is None:
                        v = biased_exponent(a, K, L, degs, xc, theta)
                    if v is not None and abs(v) < 50:
                        values.append(v)
    return out


def degree_sets(M, K, L):
    """All (N_0..N_K) with sum(N_i+1) = M - L and max-min <= 1."""
    total = M - L
    base, extra = divmod(total, K + 1)
    if base < 1:
        return []
    # choose which `extra` polynomials get one more coefficient
    out = []
    import itertools
    for idx in itertools.combinations(range(K + 1), extra):
        degs = [base - 1 + (1 if i in idx else 0) for i in range(K + 1)]
        out.append(tuple(degs))
    return out


def ensemble(a, K, Ls, Ms, guess, max_per=None, defect_radius=0.01):
    """Run approximants over (L, M, degree assignment). Returns list of Estimate,
    dropping defective ones (another zero of Q_K within defect_radius of x_c,
    or a real zero in (x_c - 0.05, x_c))."""
    theta = _theta_powers(a, K)
    out = []
    for M in Ms:
        for L in Ls:
            ds = degree_sets(M, K, L)
            if max_per:
                ds = ds[:max_per]
            for degs in ds:
                e = analyse_one(a, K, L, degs, guess, theta)
                if e is None:
                    continue
                bad = False
                for r in e.others:
                    if abs(r - e.xc) < defect_radius:
                        bad = True
                    # a spurious real singularity just below x_c would
                    # dominate the approach to x_c; ones far away are harmless
                    if abs(r.imag) < 1e-6 and e.xc - 0.05 < r.real < e.xc:
                        bad = True
                if not bad:
                    out.append(e)
    return out


def summarize(ests, label=""):
    if not ests:
        return None
    xs = np.array([e.xc for e in ests])
    ls = np.array([e.exponent for e in ests])
    med_x, med_l = float(np.median(xs)), float(np.median(ls))
    # robust spread: half the 68% interquantile range
    sx = float((np.quantile(xs, 0.84) - np.quantile(xs, 0.16)) / 2)
    sl = float((np.quantile(ls, 0.84) - np.quantile(ls, 0.16)) / 2)
    return {"label": label, "count": len(ests), "xc": med_x, "xc_spread": sx,
            "exponent": med_l, "exponent_spread": sl}


def load_series(path, key, name, allow_unverified=False):
    return load_record(path, key, allow_unverified)[name]


if __name__ == "__main__":
    import argparse
    import time
    ap = argparse.ArgumentParser()
    ap.add_argument("series_json")
    ap.add_argument("key")
    ap.add_argument("name")
    ap.add_argument("--K", type=int, default=3)
    ap.add_argument("--terms", type=int, default=None, help="use only the first T terms")
    ap.add_argument("--guess", type=float, default=0.6447)
    ap.add_argument("--Lmax", type=int, default=10)
    ap.add_argument("--Lstep", type=int, default=2)
    ap.add_argument("--span", type=int, default=12, help="M ranges over [T-span, T]")
    ap.add_argument("--allow-unverified", action="store_true", help="analyse exploratory, unverified series")
    a = ap.parse_args()
    s = load_series(a.series_json, a.key, a.name, a.allow_unverified)
    T = a.terms if a.terms is not None else len(s)
    if not 1 <= T <= len(s):
        ap.error(f"requested {T} terms; available terms: {len(s)}")
    s = s[:T]
    t0 = time.time()
    ests = ensemble(s, a.K, range(0, a.Lmax + 1, a.Lstep), range(T - a.span, T + 1), a.guess)
    print(json.dumps(summarize(ests, f"{a.key} {a.name} K={a.K} T={T}")), f"{time.time()-t0:.1f}s")
