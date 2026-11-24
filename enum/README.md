# dpenum: low-density series for 1+1 directed percolation

`dpenum` computes the moment series S, μ₀₁, μ₀₂ and μ₂₀ of the pair connectedness for bond or site directed percolation on the directed square lattice, modulo 62-bit primes. [`../analysis/crt.py`](../analysis/crt.py) combines the per-prime outputs into exact integers.

## Build and run

```
cargo build --release
target/release/dpenum -n 100 --bond --primes 3 -o out.jsonl      # 3 primes, one after another
target/release/dpenum -n 100 --site --first-prime 3 --primes 1   # just the 4th prime
python ../analysis/crt.py out.jsonl -o series.json               # exact series, checked with a spare prime
```

CRT output is refused if the spare-prime check fails. `--allow-unverified`
permits exploratory output with `verified: false`; the analysis commands also
require `--allow-unverified` to consume such data. Requested analysis orders
must fit the available coefficients. Unbiased differential approximants use
interval Newton root certification and require 30-digit precision in both
the critical point and exponent.

Options:

| Option | Meaning |
| --- | --- |
| `-n N` | Series order |
| `--bond` / `--site` | Model (bond is the default) |
| `--primes k`, `--first-prime i` | Run k primes starting at index i, from the largest primes below 2^62 |
| `--threads T` | T = 1 uses streaming steps that free each source shard once it is processed. T > 1 (default 4) runs sharded parallel blocks. Under heavy CPU contention T = 1 is faster. |
| `--park` | Single-thread mode only. Each state waits in the table of the next column that changes it, so untouched states are never copied. Half the column steps and about 11% faster, but about 46% more peak memory (n = 120). Use it when time, not memory, binds. |
| `--block K` | Parallel mode only: steps per synchronization (default 8). A block of K column steps touches only K + 1 key bits, so states are sharded by the other bits and each thread runs the block independently. |
| (parallel memory) | Since v14, each parallel block first moves the old layer into per-target pieces, freeing every source shard once split and every piece once stepped. The old and new layers are therefore never both held in full: at n = 150, all moments, u16, 6 threads, the peak working set fell from 694 to 491 MiB, with identical output. Freed heap memory is returned to Windows after every block (`DP_NO_TRIM=1` disables this). |
| `--no-x` | Store only M₀; outputs S, μ₀₁ and μ₀₂, using a third of the memory |
| `--no-mirror` | Disable the mirror reduction (testing) |
| `--loose` | Disable order pruning (testing only) |
| `-v` | Per-row progress: state count, estimated table peak, measured heap peak (counting allocator), time |

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
- `DP_NO_TRIM=1`: do not return freed heap memory to the OS after each block.

`analysis/peakmem.ps1 <exe> <args>` runs the engine and reports its peak working set, sampled every 100 ms. The heap peak printed by `-v` excludes heap fragmentation; the working set includes it.

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

Run the Rust regression tests with `cargo test`, and the analysis accuracy
regressions from the repository root with
`python -m unittest discover -s analysis -p 'test_*.py'` (requires numpy and python-flint).

## dpwall: the wall series

`src/wall.rs` builds a second binary, `dpwall`, for bond DP next to an impenetrable wall (sites with x < 0 removed, source at x = 0; Essam et al. 1996, Jensen 1999 II). It computes, modulo primes:

- T = Σ_{t≥1} P(t), the mean cluster length, where P(t) is the probability that the cluster reaches row t (exponent τ₁);
- S, Σ x C, Σ x² C, Σ t C and Σ t² C (exponents γ₁ + k ν⊥ and γ₁ + k ν∥).

```
cargo build --release --bin dpwall
target/release/dpwall -n 110 --threads 6 --block 6 --primes 3 -v -o wall110.jsonl
python ../analysis/crt.py wall110.jsonl -o wall110.json
```

**Method.** The path-union expansion, as in the bulk engine, but every state is kept: P(t) is the sum over row-t states B of (−1)^(|B|−1) times the weight of B (inclusion–exclusion over the wet sites of row t), so paths never have to merge again and the bulk engine's collapse pruning does not apply. The wall also removes translation invariance, so there is no nodal or translation reduction. Pruning is by order only: a state on the boundary between rows t and t + 1 with k occupied sites has weight O(p^(t+k−1)), and each pending occupied parent still needs a bond.

**Sweep.** One parent at a time. Step j changes only two adjacent key bits (the pending child and parent j), so a step maps groups of at most four states to themselves. Blocks of steps run in parallel on shards keyed by the untouched bits, as in the bulk engine. Targets whose pending-child bit is clear have a single source and are appended without hashing (about 1.4× faster).

**Cost** (v1: u64 residues, u128 keys, 6 threads, one prime): order 67 in 4.5 s, 80 in 45 s, 90 in 135 s, 100 in about 10 min with 8 M states and 0.6 GB, 110 in 37 min with 28 M states and 2.2 GB. States grow by about 1.14 per order and time by about 1.16.

**Memory (v2).** Tighter order bounds cannot shrink the state count: a wall state contributes to P(t) on every row, so no future bonds are owed, and a state whose weight vanishes below order n + 1 is already dropped when its window is trimmed. The savings come from bytes per state instead:
- Keys are u64 (row t uses bits up to t/2 + 3, so n ≤ 121).
- Layer tables are read only in order, so they keep no per-state arena offsets (`Table::drop_offsets`, `Table::iter`).
- Per-state overhead falls from 22 to 10 bytes. Coefficients (about 8 per state at the peak) dominate the rest, so the 31-bit build (`--features u32`) roughly halves the total.

At order 100 (8.0 M states): v1 u64 about 616 MiB, v2 u64 525 MiB, v2 u32 298 MiB (and 239 s against 393 s, from the halved memory traffic). Order 120 with v2 u32 is projected at about 4 GB and needs five 31-bit primes ([analysis/wall120_queue.sh](../analysis/wall120_queue.sh)).

**Why not swap.** Each block reads the whole layer once and writes a new one, a cyclic scan. Under paging that is the worst case: the least recently used page is the next one needed, faults are small and synchronous, and Windows wastes time compressing incompressible residues.

**Spilling to disk (`--spill DIR`).** Each block reads its input one target shard at a time, runs its steps, and splits the output by the next block's shard grouping straight into the next generation ([src/spill.rs](src/spill.rs)). A generation larger than `--spill-above` MiB (default 256) is written to disk as one file per shard, appended in chunks; reading a shard deletes its file. Smaller generations stay in memory through the same code path, so the early rows create no files. Each spilled block is one sequential read and one sequential write of the layer.
- Each disk generation reserves a private subdirectory atomically, so independent runs can safely share the spill root. Its subdirectory is removed when the generation is released.
- Memory holds only the shards in flight, about threads × 2 × layer ÷ shards. `--shards` defaults to 64 in spill mode; raise it for large layers. Every shard costs one file create and delete per block, so very many shards make small layers slow.
- Disk needs about two layers' worth of space.
- Beyond n = 121 the keys no longer fit in u64: build with `--features wide` (u128 keys, 8 more bytes per state). Features combine, for example `--features u32,wide`.
- Order 100, u32, 4 threads, one prime, on a busy machine: in memory, peak working set 600 MiB and 695 s; spilled above 32 MiB, peak working set 155 MiB, 761 s (+9.5%) and 32 GiB of disk traffic. The outputs are identical. All-disk and `wide` runs reproduce wall110 exactly at order 80.

```
cargo build --release --bin dpwall --features u32,wide --target-dir target_u32_wide
target_u32_wide/release/dpwall -n 130 --threads 6 --spill D:/dpspill --shards 256 --primes 6 -v -o wall130.jsonl
```

**Tests.** `analysis/wall_ground_truth.py N --check` computes the same series by the exact Domany–Kinzel transfer matrix and by the path-union expansion in Python, and compares both with Table 1 of Essam et al. (1996). The engine reproduces the fixture `analysis/wall_order30.json` and every entry of that table to order 67 (`analysis/egjt96_wall_table.txt`; `analysis/check_wall_table.py` does the comparison).
