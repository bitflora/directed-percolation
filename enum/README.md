# dpenum: low-density series for 1+1 directed percolation

`dpenum` computes the moment series S, μ₀₁, μ₀₂ and μ₂₀ of the pair connectedness for bond or site directed percolation on the directed square lattice, modulo 62-bit primes. [`../analysis/crt.py`](../analysis/crt.py) combines the per-prime outputs into exact integers.

## Build and run

```
cargo build --release
target/release/dpenum -n 100 --bond --primes 3 -o out.jsonl      # 3 primes, one after another
target/release/dpenum -n 100 --site --first-prime 3 --primes 1   # just the 4th prime
python ../analysis/crt.py out.jsonl -o series.json               # exact series, checked with a spare prime
```

Options:

| Option | Meaning |
| --- | --- |
| `-n N` | Series order |
| `--bond` / `--site` | Model (bond is the default) |
| `--primes k`, `--first-prime i` | Run k primes starting at index i, from the largest primes below 2^62 |
| `--threads T` | T = 1 uses streaming steps that free each source shard once it is processed. T > 1 (default 4) runs sharded parallel blocks. Under heavy CPU contention T = 1 is faster. |
| `--block K` | Parallel mode only: steps per synchronization (default 8). A block of K column steps touches only K + 1 key bits, so states are sharded by the other bits and each thread runs the block independently. |
| `--no-x` | Store only M₀; outputs S, μ₀₁ and μ₀₂, using a third of the memory |
| `--no-mirror` | Disable the mirror reduction (testing) |
| `--loose` | Disable order pruning (testing only) |
| `-v` | Per-row progress: state count, peak memory, time |

Each prime's line is appended to the output file as soon as it finishes, so independent primes can run in separate processes and be combined afterwards.

Build variants:
- `cargo build --release` stores residues modulo 62-bit primes in u64.
- `cargo build --release --features u32` stores 31-bit residues in u32. That is 42% less memory per run, at the cost of twice as many primes. Prefer it when memory binds.
- `cargo build --release --features u16` stores 15-bit residues in u16, as Jensen 1999 did. It uses the least memory and needs four times the primes of u64, so use it when memory is the only constraint.
- Memory at n = 110, single prime:

| Build | All moments | `--no-x` |
| --- | --- | --- |
| u64 | 41.6 MiB | 18.1 MiB |
| u32 | 23.9 MiB | 12.2 MiB |
| u16 | 15.1 MiB | 9.3 MiB |

Diagnostic environment variables:
- `DP_SLACK=1`: window slack against the bound.
- `DP_NOSIGN=1`: count graphs without signs.
- `DP_CHEAP_LO=1`: the rejected cheap lower bound.

With `-v`, a coarse profile is also printed: step-loop time, finish time, identity copies, and time in the bound.

## Algorithm

The engine follows Jensen 1999's non-nodal path-union transfer matrix:
- C(x,t) is a sum over unions of directed paths with weight (−1)^(number of vertices with two incoming bonds) p^(bonds), or p^(sites) for site percolation.
- The boundary moves one site at a time, so each step rewrites two bits.
- Non-nodal graphs are kept by deleting single-vertex rows.
- The full series come from S = 1/(1 − S_N) and Jensen's eqs. 7–9.

The update rules were re-derived from scratch. They agree with the parts of Jensen's eq. 12 that survived text extraction.

Changes relative to 1999:

1. **Translation reduction.** States are cluster shapes, normalized to start at column 0 at each row end. Absolute position is carried as three moment polynomials M_k = Σ_o o^k P, which is exactly what x² needs. Translates share one entry.
2. **Mirror reduction.** A shape and its mirror are stored once, with weight 1/2 for symmetric shapes. Moments transform as M₁ → L M₀ − M₁ and M₂ → L² M₀ − 2L M₁ + M₂, with L = t − width.
3. **Interval-cover bounds** (`bound` in `src/main.rs`).
   - Every vertex has an ancestor d rows back and a descendant d rows ahead, within d columns.
   - So the minimum number of windows of length d covering the occupied columns bounds the vertex count of those rows.
   - Summing gives a lower bound `lo` on the bonds (sites) already placed, and a lower bound on those still needed. Partially processed columns are treated as fuzzy intervals.
   - Each state stores only coefficients `lo..hi`, and a state with lo > hi is dropped. Tightening `lo` alone cuts the state count by 4.6× at n = 100.
4. **No hashing for unique targets.** Targets that only one source can produce are appended without hashing, and identity steps reuse the source's window.
5. **Exact zeros and trimming.** All-zero states (exact cancellations) are dropped, and every window is trimmed to its nonzero range after each step. A target that only one source produces inherits the source's exact lower end.
6. **Mask-based bounds.** The interval covers run on 128-bit masks, a few bit operations per window.
7. **62-bit primes** (or 31-bit with `--features u32`). One run carries 62 bits; Jensen 1999 stored 16-bit residues and needed up to 12 primes.

The state count still grows as 2^(n/8), the same rate as Jensen's algorithm; the reductions only lower the prefactor.

## Tests

[`../analysis/m0_ground_truth.py`](../analysis/m0_ground_truth.py) computes the series three independent ways, in exact integers:
- inclusion–exclusion over sets of paths (t ≤ 4), which checks the sign rule;
- the Domany–Kinzel probability transfer matrix;
- a direct row-by-row path-union sum, with and without the non-nodal reduction.

It writes fixtures (`m0_order12.json`, `m0_order15.json`). `crt.py --check` compares engine output against them.

Further checks used during development:
- tight pruning agrees exactly with `--loose` at n = 36–40;
- mirror and no-mirror runs agree, and 1 and 4 threads agree;
- u32 (31-bit) and u64 (62-bit) runs agree modulo 2^31 − 1;
- the order-50 series is an exact prefix of the order-70 series.
