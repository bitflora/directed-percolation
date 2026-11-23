# Findings

Running record of results for the directed percolation study. Each claim lists its evidence. Numbers marked *provisional* come from runs or analyses still in progress.

## Summary of the first night (November 21–22, 2026)

- **Engine.** An independent reimplementation of Jensen's 1999 algorithm, checked against exact ground truth (F1, F2). It adds translation and mirror reductions, exact window trimming and faster bounds (F3). Per prime it needs roughly 6× less memory than the 1999 implementation.
- **Verified exact series:**
  - bond to order 156 for all four moment series (S, μ₀₁, μ₀₂, μ₂₀) (F4);
  - bond to order 171 for all four moment series, Jensen's own length (F10, F11, F12; μ₂₀ added on the second night);
  - site to order 158 for all four moment series, Jensen's own length (F8, F13).

  Jensen reached 171 (bond) and 158 (site).
- **Reproduction of the 1999 results.**
  - At 171 terms (F11), the unbiased p_c, γ, ν∥ and the μ₀₁ and μ₀₂ exponents all agree with Jensen's values within his error bars.
  - At Jensen's p_c, biased approximants give γ, ν∥, 2ν⊥ and γ + 2ν⊥ within his error bars.
  - The ν∥ and 2ν⊥ estimates are stable from 146 to 156 terms.
  - The site critical point from universality is 0.70548522, as he found (F5, F8).
- **Gate (M2) passed.** Estimates of ν∥ and ν⊥ improve about 1.6–2× per 10 extra terms (fits over 90–171 terms), so 200–230 terms should give about one more significant digit (F5, F10, F11).
- **Solvability.** There is no small linear ODE for any of the series (F6). On rings, the poles of the mean cluster size pinch toward p_c and cluster around the bulk non-physical singularities: a first look at the haruspicy route (F9).
- **Second night (November 22–23).** μ₂₀ at 171 terms completes Jensen's set; every unbiased estimate, including 2ν⊥ = 2.1937077 and γ + 2ν⊥ = 4.4714139, agrees with his. The μ₂₀μ₀₂/μ₀₁² critical point (0.6447001905) sits 5 × 10⁻⁹ above the S and μ₀₂/μ₀₁ ones, inside his error bar; this gap is the clearest target for longer series. A split-layer parallel step cut peak memory by about 30% (F12). The site series then reached 158 terms in 1 h 48 min; every site estimate agrees with Jensen's, and biased approximants give p_c(site) = 0.70548521, against his 0.70548522(4) (F13).
- **Compute decision (M4).** Measured: order 171 for S, μ₀₁ and μ₀₂ takes 12 min per prime with 6 threads and 2.2 GB on this laptop (F11), and parallel efficiency is about 60–75% (F7). Orders 200 and above need a large-memory machine; projections are in F7.

## F1. Independent ground truth for small orders (M0)

[analysis/m0_ground_truth.py](analysis/m0_ground_truth.py) computes the four moment series (S, μ₀₁, μ₀₂, μ₂₀) for bond and site DP in exact integers, three independent ways:

1. inclusion–exclusion over sets of directed paths, for t ≤ 4;
2. the Domany–Kinzel probability transfer matrix;
3. the path-union graph expansion summed row by row, both directly and through the non-nodal reduction.

All three agree to order 15 for bond and site. The fixtures are [analysis/m0_order12.json](analysis/m0_order12.json) and [analysis/m0_order15.json](analysis/m0_order15.json).

The bond series begins S = 1 + 2p + 4p² + 8p³ + 15p⁴ + 28p⁵ + 50p⁶ + 90p⁷ + 156p⁸ + 274p⁹ + ….

Two details of the method had to be pinned down:
- **Sign rule.** The path-union sign is (−1)^(number of vertices with two incoming bonds), for both bond and site, with site weights p per occupied vertex. Inclusion–exclusion confirms this.
- **Site row count.** For site percolation, a non-nodal piece ending at row t can have as few as 2t − 1 new sites. The site calculation therefore needs rows up to (n + 1)/2, one more than the bond rule n/2. Using n/2 rows gives a wrong coefficient at order 9, which the cross-check caught.

## F2. Reimplemented and improved transfer-matrix engine (M1)

[enum/](enum/README.md) is a Rust implementation of Jensen's non-nodal path-union algorithm, with the update rules re-derived independently.

**Additions over 1999:**
- translation reduction, with position moments carried as three polynomials;
- mirror reduction;
- per-state coefficient windows from an interval-cover lower bound;
- 62-bit modular arithmetic.

**Validation:**
- exact agreement with F1 to order 15 (bond and site);
- tight pruning equals no pruning at n = 36–40;
- mirror and no-mirror runs agree, and 1 and 4 threads agree;
- 31-bit and 62-bit runs agree modulo 2³¹ − 1;
- the order-50 series is an exact prefix of the order-70 series;
- the cover-function speedup matches the linear version on 20,000 random point sets.

**Growth:** the peak state count grows by 1.09 per term (single prime, bond):

| n | 60 | 70 | 80 | 90 | 100 | 110 |
| --- | --- | --- | --- | --- | --- | --- |
| peak states | 1,622 | 4,218 | 10,853 | 26,668 | 65,551 | 155,778 |

The rate is 2^(1/8) = 1.0905 per term, the same as Jensen 1999's algorithm. The reductions lower the prefactor, not the growth rate.

**Memory compared with 1999:**
- Order 126 needs 210 MiB, storing three moments as 8-byte residues. Extrapolating at 1.09 per term gives about 10 GB at order 171.
- In Jensen's storage format (one polynomial, 2-byte residues), the same states would need roughly 10 GB / 12 ≈ 0.85 GB, against the 5.5 GB he reports for order 171.
- So the shape, mirror and bound reductions are worth roughly a factor of 6, or about 20 extra terms at fixed memory. This is *provisional*: it extrapolates our measurements and assumes his memory figure is comparable.
- Measured at order 171 (F11): S, μ₀₁ and μ₀₂ with 4-byte residues peaked at 2.2 GB in parallel mode. That is roughly 1.6 GB in streaming mode, allowing for the 40% parallel overhead. Jensen's 5.5 GB covered all four series with 2-byte residues, so the comparison is indicative only.

**Bounds matter most.** Replacing the interval-cover lower bound on bonds already placed with the simple bound 2s multiplies the state count by 4.6 at n = 100.

**Compute environment.** All timings were taken while another project's job (`mixlab`) occupied all 12 logical cores. Only about 3–4 GB of the 16 GB RAM was free. Multithreaded steps ran slower than single-threaded ones under this contention, because every column step waits for the slowest thread. Production therefore runs one process per prime. Absolute times overstate what the machine can do.

## F3. Engine optimization results (M3, partly)

Each change below was verified to give identical series: the order-15 fixture, and exact agreement with the earlier binary at order 70 for bond and site.

| Change | Effect at n = 100 (single prime) |
| --- | --- |
| Interval-cover lower bound on degree already placed, instead of 2s | 4.6× fewer states than the simple bound |
| Mirror reduction | Bound calls 13.7M → 5.8M; time about 1.6–2.8× lower |
| Streaming steps (free each source shard after use) | Peak memory about −30% (no full double buffer) |
| Drop all-zero states | Peak states −15% (bond), −26% (site) |
| Trim each window to its nonzero range at the end of every step | Peak states −21% (bond), −31% (site) combined with the previous row; memory −8% / −19% |
| `--features u32` (31-bit primes) | Memory −42% per run; twice as many primes |
| Binary-search cover | No measurable change under contention |
| 128-bit mask cover (each window a few bit operations; checked against the point-list version on 400,000 random states) | Bound time −25%, total about −10% |
| Single-source targets inherit the source's trimmed lower end plus the degree shift, skipping the past-cover loop | Bound time 5.2 s → 3.3–4.0 s, total about −20%, memory −5% |
| Rejected: lower end 2s + fixed for hashed targets, relying on trimming | Peak memory +62%, no time gain |
| Per-state bookkeeping: 8-bit window ends, 32-bit per-shard offsets | Memory −4% (u32), −9% (u32 `--no-x`) at n = 110 |
| `--features u16` (15-bit primes, Jensen's choice) | At n = 110 with `--no-x`: 9.3 MiB, against 12.2 MiB for u32 and 18.1 MiB for u64; four times the primes of u64 |
| Drop the redundant zero check; copy identity blocks without zero-filling | About −3% (n = 120, idle machine: 30.6 s → 29.7 s) |
| Rejected: `-C target-cpu=native` (AVX2) | About −2–3%; not worth the loss of portable binaries |
| `--park` (opt-in): each state waits in the table of the next column that changes it, so untouched states are never copied | At n = 120: state-steps 54M → 26M, time −11% (30.4 s → 27.2 s), peak memory +46% (104 → 152 MiB). Pending tables keep hash slots and untrimmed windows until they are processed. Use it when time binds. |

**Profile on the idle machine** (bond, n = 120, one prime, single thread, 30.0 s in total):

| Part | Time | Share |
| --- | --- | --- |
| Bound computation (32M calls, 0.26 µs each) | 8.4 s | 28% |
| Rest of the column steps: hashing, accumulation, copies (54M state-steps, 28M of them identity copies) | 17.3 s | 58% |
| Window trimming and compaction (`finish`) | 2.9 s | 10% |
| Row ends | about 1.3 s | 4% |

The next speed-ups worth trying:
- **Make `--park` memory-neutral** by trimming pending tables as they fill.
- **Fewer cache misses** in the target hash.
- **Cheaper bounds** for hashed targets.

**Bound tightness.** A diagnostic mode counts graphs without signs, so nothing cancels.
- Even then, the interval-cover lower bound on degree placed is 1–3 below the true structural minimum for 74% of states at n = 90. The looseness is structural, not cancellation.
- Window trimming after every step already replaces that bound with the exact minimum, so the states *kept* between steps are already near-minimal. A tighter static bound would still save construction-time memory and work, because fewer doomed states would be created and windows would be allocated shorter. It would not change the state count between steps.
- The completion bound is exact on worked examples.
- Further large savings would need a different algorithm, for example the Zbarsky-type sparsification in REVIEW §10.1, not tighter bounds.

**Profile.** With bounds computed twice instead of once, run time rises from 12.2 s to 15.4 s. So hashed-target bounds take about a quarter of the time, and all bounds roughly a third. Hashing, memory traffic and accumulation take the rest.

## F4. Verified bond series to orders 146 and 156

[results/series/bond146.json](results/series/bond146.json) holds the exact S, μ₀₁, μ₀₂ and μ₂₀ series for bond DP on the directed square lattice to order 146 (Jensen 1999 reached 171).

**Verification:**
- The coefficients were reconstructed from four 62-bit primes and checked with the fourth.
- The largest coefficient has 140 bits.
- The independently computed 120-term series (three other primes) is an exact prefix.

**Cost:** each prime took 31 minutes on a single thread with all cores contended, and peaked at 1.19 GB with the v4 binary.

The coefficients cannot be compared term by term with Jensen's, because his series files are not available online. The checks against his series are therefore the critical-point and exponent estimates in F5.

**Extended to order 156.** [results/series/bond156.json](results/series/bond156.json) is the same four series to order 156.
- Reconstructed from seven 31-bit primes and verified (the largest coefficient has 149 bits).
- The 146-term series is an exact prefix.
- No ODE with (k + 1)(d + 1) ≤ 147 exists for any of the four series ([results/analysis/dfinite_bond156.txt](results/analysis/dfinite_bond156.txt)).

## F5. Reproduction of the 1999 estimates and the precision-gain gate (M2)

The 146-term bond series was analysed by [analysis/gate.py](analysis/gate.py) at T = 80, 90, …, 146 terms. Each estimate is the median over inhomogeneous differential approximants of order K = 2 and 3:
- L = 0, 3, 6, 9;
- M = T−10 … T in steps of 2;
- up to two degree assignments per (M, L).

Defective approximants are dropped: another zero of Q_K within 0.002 of x_c, or a real zero in (x_c − 0.05, x_c). The full table is in [results/analysis/bond146_gate_report.md](results/analysis/bond146_gate_report.md).

**Reproduction.** At 146 terms the estimates agree with Jensen's 171-term values:

| Quantity | Series | This work, T = 146 | Jensen 1999, T = 171 |
| --- | --- | --- | --- |
| ν∥ | μ₀₂/μ₀₁ | 1.7338449 (spread 1.6e-6) | 1.733847(6) |
| 2ν⊥ | μ₂₀μ₀₂/μ₀₁² | 2.1937084 (spread 1.9e-6) | 2.193708(4) |
| γ | S | 2.2777184 (spread 4.8e-6) | 2.277730(5) |
| γ + 2ν⊥ | μ₂₀ | 4.471404 (spread 1.1e-5) | 4.471425(15) |
| p_c | μ₀₂/μ₀₁, μ₂₀μ₀₂/μ₀₁² | 0.6447001840, 0.6447001916 | 0.644700185(5) |
| p_c | S | 0.6447001775 | |

For ν∥, ν⊥ and p_c the 146-term values already sit inside Jensen's error bars. The S, μ₀₁, μ₀₂ and μ₂₀ estimates are still drifting upward toward his values, by about 1–5 × 10⁻⁷ per term at T = 146. Jensen reported the same upward drift.

**Error measure.** For every series the approximant spread at fixed T is smaller than the change between successive T. So the spread understates the error, and the gate uses the drift (change of the central estimate per term) instead.

**Gate result** (provisional: the fits pass through noisy medians). log10(drift per term) was fitted linearly in T over 90–146 and extrapolated as a geometric tail:

| Quantity | Error reduction per 10 terms | Projected gain 171 → 230 terms |
| --- | --- | --- |
| ν∥ (R1) | 1.8× | about 30× |
| 2ν⊥ (R2) | 1.9× | about 40× |
| γ (S) | 1.4× | about 6× |
| γ + 2ν⊥ (μ₂₀) | 1.1× | about 2× |

The plan's criterion was a 3× gain at 230 terms or 2× at 200. It is met comfortably for ν∥, ν⊥, γ and p_c; μ₂₀ alone converges too slowly to matter.

**Verdict: the M2 gate is passed.** Extending the series to about 200–230 terms should sharpen ν∥ and ν⊥ by roughly one order of magnitude beyond the 1999 values.

**Caveats:**
- The extrapolations assume the drift keeps decaying geometrically.
- Our projected remaining error at 171 terms (about 2 × 10⁻⁷ for ν∥) is much smaller than Jensen's quoted ±6 × 10⁻⁶. Either his error bars were conservative or the drift proxy is optimistic; a longer series is what will tell.

**Update with 151 and 156 terms** ([results/analysis/bond156_gate_report.md](results/analysis/bond156_gate_report.md)). The ratio-series estimates have stabilized, and γ and γ + 2ν⊥ keep rising toward Jensen's values:

| Quantity | T = 146 | T = 151 | T = 156 | Jensen 1999 (171 terms) |
| --- | --- | --- | --- | --- |
| ν∥ (R1) | 1.7338449 | 1.7338449 | 1.7338448 (spread 7e-7) | 1.733847(6) |
| 2ν⊥ (R2) | 2.1937084 | 2.1937081 | 2.1937080 (spread 8e-7) | 2.193708(4) |
| γ (S) | 2.2777184 | 2.2777220 | 2.2777224 | 2.277730(5) |
| γ + 2ν⊥ (μ₂₀) | 4.4714038 | 4.4714069 | 4.4714110 | 4.471425(15) |
| x_c (R1 / R2) | 0.6447001840 / …1916 | 0.6447001840 / …1911 | 0.6447001839 / …1908 | 0.644700185(5) |

With the 151- and 156-term points included, the drift fits give a 60× (ν∥) and 40× (2ν⊥) projected gain from 171 to 230 terms, consistent with the 146-term analysis.

**Biased approximants.** Q_K is forced to vanish at x_c (K = 2, L = 0, 2, …, 10, M = 136 … 146; [analysis/biased_scan.py](analysis/biased_scan.py), output [results/analysis/bond146_biased.json](results/analysis/bond146_biased.json)). At Jensen's p_c = 0.644700185 the 146-term series gives:

| Quantity | Biased estimate, 146 terms | Jensen 1999 | Difference |
| --- | --- | --- | --- |
| γ | 2.2777279 ± 0.6e-6 | 2.277730(5) | −2.1e-6 |
| ν∥ | 1.7338458 ± 0.1e-6 | 1.733847(6) | −1.2e-6 |
| 2ν⊥ | 2.1937038 ± 1.1e-6 | 2.193708(4) | −4.2e-6 |
| γ + 2ν⊥ | 4.4714237 ± 3.7e-6 | 4.471425(15) | −1.3e-6 |

The ± values are the spread over approximants. All four agree with Jensen within his error bars. Repeating the scan with the 156-term series ([results/analysis/bond156_biased.json](results/analysis/bond156_biased.json)) gives, at p_c = 0.644700185: γ = 2.2777275, ν∥ = 1.7338459, 2ν⊥ = 2.1937016 and γ + 2ν⊥ = 4.4714237. These are unchanged from 146 terms within 1 × 10⁻⁶ (2ν⊥ within 2.2 × 10⁻⁶), so the biased estimates have converged at this length.

Sensitivity to the assumed critical point is dγ/dp_c ≈ 1.5 × 10³ and dν∥/dp_c ≈ 1.0 × 10³. So Jensen's ±5 × 10⁻⁹ on p_c alone allows ±5–8 × 10⁻⁶ in these exponents. The biased S and μ₀₂/μ₀₁ estimates equal Jensen's γ and ν∥ at x_c = 0.64470018615 and 0.64470018620, inside his p_c interval.


## F6. No small linear ODE (D-finiteness search)

[analysis/dfinite.py](analysis/dfinite.py) solves the linear system for an ODE Σ_{i≤k} P_i(x) θⁱ f = 0 with deg P_i ≤ d modulo two 62-bit primes. It uses 10 more equations than unknowns, and was checked on a known D-finite series (central binomial coefficients, found at order 1, degree 1).

For all four 146-term bond series there is no such ODE whenever (k + 1)(d + 1) ≤ 137, for k = 1 … 10. Examples of excluded sizes:

| Order k | Largest degree d excluded |
| --- | --- |
| 1 | 67 |
| 3 | 33 |
| 6 | 18 |
| 10 | 11 |

This is a weak exclusion. Known D-finite lattice series often need ODEs far larger than this (the Ising χ̃⁽⁵⁾ ODE has order 33). It is reported as a baseline for the haruspicy approach (RESEARCH_PLAN M6), not as evidence of non-D-finiteness. Output: [results/analysis/dfinite_bond146.txt](results/analysis/dfinite_bond146.txt).

## F7. Measured cost and projections for the compute decision (M4)

**Measured** (single thread, one prime per process, all 12 logical cores also busy with another job):

| Run | Binary | Peak memory | Time per prime |
| --- | --- | --- | --- |
| bond, n = 126 | u64, v4 | 210 MiB | 2.8 min |
| bond, n = 146 | u64, v4 | 1.19 GB | 31 min |
| bond, n = 156 | u32, v8 | 1.33 GB | 62–65 min with 4–5 of our processes running at once; 32–41 min with fewer |
| bond, n = 156 | u32, v10/v11 | 1.31 GB | 40 min with 4 of our processes running; 21–22 min with 3 |

**Growth:** states and memory grow by 1.09 per term, time by about 1.10 per term. Since these runs, engine v10 is about 20% faster and 5% leaner, and `--no-x` (dropping μ₂₀) roughly halves memory.

**Projections** (rough; they scale the measured numbers by 1.09–1.10 per term). The bond series needs about 0.96 bits per term, so 7 primes of 31 bits or 4 of 62 bits including a verification prime.

| Target order | Memory per prime (u32, all moments) | Single-thread time per prime, contended | Feasible where |
| --- | --- | --- | --- |
| 171 (reproduce 1999) | about 4.6 GB (about 2.5 GB with `--no-x`) | about 1.5 h with v11 (scaled from 22 min at n = 156) | This PC: about 10 h for 7 primes with all moments in one process, or about 4–5 h with `--no-x` in two processes, plus a separate μ₂₀ run |
| 200 | about 57 GB | about 70 h (much less multithreaded and uncontended) | A 128 GB node with the sharded multithreaded mode (`--threads`), roughly 3 h per prime at 16 threads (estimate) |
| 230 | about 0.8 TB | far beyond single-thread | A 1 TB node, roughly 15–20 h per prime at 32 threads (estimate), so several days for 7 primes |

**Parallel scaling** (measured at 03:29, after the other job had stopped). Bond, n = 120, one prime, u64, block-parallel mode (`--threads T --block 8`). The 6-core / 12-thread laptop was still running four other single-threaded processes of ours.

| Threads | Time | Speed-up | Peak memory |
| --- | --- | --- | --- |
| 1 (streaming) | 51.4 s | 1.0× | 104 MiB |
| 2 | 31.2 s | 1.65× | 146 MiB |
| 4 | 17.4 s | 2.95× | 146 MiB |
| 6 | 14.1 s | 3.65× | 146 MiB |
| 6, `--block 16` | 29.5 s | 1.7× | 145 MiB |

Parallel efficiency is about 75% at 4 threads and 60% at 6 here, limited by the busy physical cores. The parallel mode lacks streaming, so it needs about 40% more memory. The block size matters. At 6 threads, blocks of 2, 4, 6, 8 and 12 took 15.8, 13.9, 14.0, 14.2 and 17.2 s, so 4–8 is best (default 8). In the first measurement, a block of 16 took twice as long as 8. Scaling on a many-core server is still unmeasured; it is the first thing to check before renting one.

## F8. Site series to orders 136 and 146, and the site critical point from universality

[results/series/site136.json](results/series/site136.json) holds the exact site series to order 136 (Jensen 1999 reached 158). It was reconstructed from four 62-bit primes and verified with the fourth; the largest coefficient has 152 bits, and the series matches the M0 ground truth to order 15.

**Unbiased analysis** ([results/analysis/site136_gate_report.md](results/analysis/site136_gate_report.md)). The site series converges more slowly and less regularly than the bond series, as Jensen also found. At T = 136:

| Quantity | This work, T = 136 | Jensen site, T = 158 | Bond value |
| --- | --- | --- | --- |
| x_c | 0.7054851–0.7054853 | 0.70548515(20) | |
| γ | 2.27760 | 2.27765(6) | 2.277730 |
| ν∥ (R1) | 1.73378 | 1.73381(4) | 1.733847 |
| 2ν⊥ (R2) | 2.19378 | 2.19377(5) | 2.193708 |

**Biased analysis.** Following Jensen 1999 §4.2, biased approximants were scanned over x_c and the points found where each biased site exponent equals the bond value. Output: [results/analysis/site136_biased.json](results/analysis/site136_biased.json).

| Series | Bond exponent imposed | Site x_c where it is matched |
| --- | --- | --- |
| S | γ = 2.277730 | 0.70548526 |
| μ₀₂/μ₀₁ | ν∥ = 1.733847 | 0.70548523 |
| μ₂₀μ₀₂/μ₀₁² | 2ν⊥ = 2.193708 | 0.70548516 |
| μ₂₀ | γ + 2ν⊥ = 4.471425 | 0.70548525 |

The mean of the four is 0.70548522 with a spread of about ±4 × 10⁻⁸, matching Jensen's biased estimate p_c(site) = 0.70548522(4) from 158 terms. The two lattices are consistent with a common set of exponents (universality) at this precision.

**Extended to order 146** ([results/series/site146.json](results/series/site146.json)):
- seven 31-bit primes, verified; the largest coefficient has 163 bits; the 136-term series is an exact prefix;
- 17–18 min per prime once the other job had stopped.

Unbiased estimates at 146 terms ([results/analysis/site146_gate_a.json](results/analysis/site146_gate_a.json)) barely moved from 136:

| Quantity | T = 146 |
| --- | --- |
| x_c | 0.7054851–0.7054853 |
| γ | 2.27759–2.27764 |
| ν∥ | 1.73379–1.73380 |
| 2ν⊥ | 2.19376–2.19379 |

Biased crossings at 146 terms ([results/analysis/site146_biased.json](results/analysis/site146_biased.json)):
- μ₀₂/μ₀₁ at 0.70548523 and μ₂₀μ₀₂/μ₀₁² at 0.70548519, both consistent with 0.70548522(4).
- The biased S and μ₂₀ estimates became erratic above x_c ≈ 0.7054852 (spread up to 3 × 10⁻⁴) and did not cross the bond values in the scanned range, so they are not used.

The site series converges slowly. As in 1999, the site critical point is best fixed by imposing the bond exponents.

## F9. Exploratory: zeros of the ring transfer matrix (start of M6, haruspicy)

[analysis/strip_poles.py](analysis/strip_poles.py) computes det(I − T_L(p)) exactly for bond DP on a ring of L sites. It uses states up to rotation and fraction-free elimination over Z[p], then finds all complex zeros with python-flint's rigorous root finder. These zeros are the candidate poles of the ring mean cluster size S_L(p), which is a rational function of p. On rings, S_L agrees with the bulk series through order 2L − 1 (checked for L = 4, 6, 8).

| L | Orbit states | det degree | Zero nearest p_c |
| --- | --- | --- | --- |
| 5 | 7 | 32 | 0.7163 ± 0.5275i |
| 6 | 13 | 58 | 0.6749 ± 0.4064i |
| 7 | 19 | 102 | 0.6577 ± 0.3348i |
| 8 | 35 | 198 | 0.6489 ± 0.2866i |
| 9 | 59 | 362 | 0.6439 ± 0.2516i |
| 10 | 107 | 710 | 0.6408 ± 0.2249i |

**Observations** (exploratory; see the plot [results/analysis/strip_poles.png](results/analysis/strip_poles.png)):
- **Pinching at p_c.** A pair of zeros pinches toward p_c as L grows, in the way Yang–Lee edges do. The effective exponent from successive imaginary parts falls toward 1/ν⊥ ≈ 0.91: it is about 1.1 between L = 8 and 9.
- **Clustering at the bulk non-physical singularities.** Zeros cluster around −0.517 on the negative axis and −0.226 ± 0.440i, the singularities Jensen located in the bulk series. With growing L they fill out bands and arcs rather than isolated points.

**Cancellation.** By the matrix determinant lemma, S_L = det(I − T + n e_sᵀ)/det(I − T) − 1, and the script removes the common factor of numerator and denominator. Nothing cancels for L ≤ 6. The cancelled degree is 6, 22, 65 and 162 for L = 7, 8, 9 and 10, so the true pole counts are 83, 161, 280 and 529. The reduced rational functions reproduce the bulk series through order 2L − 1, and the poles nearest p_c are unchanged. The plot of the true poles is [results/analysis/strip_poles_reduced.png](results/analysis/strip_poles_reduced.png).

**What this does not show.**
- A non-D-finiteness argument would need the limit of these poles as L → ∞ and a proof that they accumulate densely. Neither is established here.

L = 10 took about ten minutes. An attempt at L = 11 (187 orbit states) was stopped after 3.2 CPU-hours and 2.8 GB without finishing even the first determinant: Python-level Bareiss elimination on polynomial matrices does not scale past L ≈ 10. Going to L ≈ 14–16 for M6 will need determinants computed modulo primes at many evaluation points in compiled code (for example inside the Rust engine), followed by interpolation and CRT.

## F10. Bond S, μ₀₁, μ₀₂ to order 166

[results/series/bond166.json](results/series/bond166.json) holds the exact S, μ₀₁ and μ₀₂ series to order 166. It was computed with `--no-x`, so μ₂₀ is not included (μ₂₀ remains at order 156).
- Reconstructed from six 31-bit primes and verified with a seventh, 15-bit prime from the u16 build. Mixed prime sizes work because the reconstruction is by the Chinese remainder theorem.
- The largest coefficient has 159 bits.
- The 156-term series is an exact prefix.
- Each prime took about 37–41 min (u32 and u16, three processes at once, `mixlab` busy).

Estimates (full table: [results/analysis/bond166_gate_report.md](results/analysis/bond166_gate_report.md); biased: [results/analysis/bond166_biased.json](results/analysis/bond166_biased.json)):

| Quantity | T = 156 | T = 161 | T = 166 | Jensen 1999 (171) |
| --- | --- | --- | --- | --- |
| γ (S, unbiased) | 2.2777224 | 2.2777244 | 2.2777279 | 2.277730(5) |
| x_c (S, unbiased) | 0.6447001808 | 0.6447001824 | 0.6447001839 | 0.644700185(5) |
| ν∥ (R1, unbiased) | 1.7338448 | 1.7338450 | 1.7338462 | 1.733847(6) |
| x_c (R1, unbiased) | 0.6447001839 | 0.6447001841 | 0.6447001852 | 0.644700185(5) |
| γ (biased at p_c = 0.644700185) | 2.2777275 | | 2.2777278 | 2.277730(5) |
| ν∥ (biased at p_c = 0.644700185) | 1.7338459 | | 1.7338459 | 1.733847(6) |

**Reading.**
- At 166 terms the unbiased critical point from μ₀₂/μ₀₁ is 0.6447001852, Jensen's central value.
- The unbiased γ and ν∥ estimates are within 2.1 × 10⁻⁶ and 0.8 × 10⁻⁶ of his.
- The biased estimates at his p_c have not changed by more than 3 × 10⁻⁷ between 146 and 166 terms.

These numbers come from an independently written program, and they are consistent with Jensen's 1999 series and analysis.

## F11. Bond S, μ₀₁, μ₀₂ at Jensen's own length, 171 terms

[results/series/bond171.json](results/series/bond171.json) holds the exact S, μ₀₁ and μ₀₂ series to order 171, the length of Jensen's 1999 bond series. μ₂₀ is not included (`--no-x`).
- Reconstructed from seven 31-bit primes and verified; the largest coefficient has 164 bits; the 166-term series is an exact prefix.
- The analysis below ran on the 6-prime reconstruction while the 7th prime was computing. That reconstruction proved identical to the verified series.

**Cost.** One process with 6 threads in block-parallel mode (`--threads 6 --block 6`, u32, `--no-x`) took 11.0–12.1 min per prime, about 80 min for all seven, with a 2.2 GB peak. This was on the laptop after the other job had stopped.

**Estimates at equal length** ([results/analysis/bond171_gate_report.md](results/analysis/bond171_gate_report.md), [results/analysis/bond171_biased.json](results/analysis/bond171_biased.json)):

| Quantity | This work, 171 terms | Jensen 1999, 171 terms | Difference |
| --- | --- | --- | --- |
| p_c from μ₀₂/μ₀₁ (unbiased) | 0.6447001855 | 0.644700185(5) | +5e-10 |
| p_c from S (unbiased) | 0.6447001847 | | −3e-10 |
| ν∥ (μ₀₂/μ₀₁, unbiased) | 1.7338467 | 1.733847(6) | −3.5e-7 |
| γ (S, unbiased) | 2.2777288 | 2.277730(5) | −1.2e-6 |
| γ + ν∥ (μ₀₁) | 4.0115503 | 4.01156(1) | −9.7e-6 |
| γ + 2ν∥ (μ₀₂) | 5.7453796 | 5.74539(1) | −1.0e-5 |
| γ (biased at 0.644700185) | 2.2777276 | 2.277730(5) | −2.4e-6 |
| ν∥ (biased at 0.644700185) | 1.7338459 | 1.733847(6) | −1.1e-6 |

**Conclusions.**
- Every estimate agrees with Jensen's within his quoted error, from series computed by an independently written program and analysed with an independently written approximant code. This is the strongest available evidence, short of comparing his coefficient files, that the two series are the same and that his 1999 numbers are sound.
- The biased estimates are identical to the 146-, 156- and 166-term values within about 3 × 10⁻⁷, so at fixed p_c the limiting uncertainty is p_c itself.
- dν∥/dp_c ≈ 1.0 × 10³, so pinning p_c to ±1 × 10⁻¹⁰ would fix ν∥ to about ±1 × 10⁻⁷. That is the case for the longer series in M5.

## F12. All four bond moment series at 171 terms: μ₂₀ completes Jensen's 1999 set

[results/series/bond171.json](results/series/bond171.json) now holds μ₂₀ to order 171, alongside the verified S, μ₀₁ and μ₀₂ (F11). This completes the set of series Jensen 1999 analysed.

**Verification.**
- μ₂₀ was reconstructed from twelve 15-bit primes ([analysis/combine_bond171_mu20.py](analysis/combine_bond171_mu20.py)).
- Its largest coefficient has 157 bits; eleven primes give 165 bits, and the twelfth prime changes nothing (leave-last-out check).
- The order-156 μ₂₀ series is an exact prefix.
- Every one of the twelve runs also recomputed S, μ₀₁ and μ₀₂ and reproduced the 171-term series of F11 modulo its prime. That makes twelve independent checks of the earlier result, by a different build (u16 against u32) with different moment storage (all three moments against `--no-x`).

**Engine change (v14).** Parallel blocks now split the old layer into per-target pieces and free each piece once it has been stepped. Previously the whole old layer and the whole new layer were held at the same time. Freed heap memory is returned to Windows after each block. Output is identical.

Peak working set, order 150, all moments, u16, 6 threads, measured with [analysis/peakmem.ps1](analysis/peakmem.ps1):

| Engine | Peak working set | Time |
| --- | --- | --- |
| v13 | 694 MiB | 272 s |
| v14 | 486–491 MiB | 235–256 s |

**Cost at order 171** (all moments, u16, `--threads 6 --block 6`, with another job using all 12 logical cores):
- 22–24 min per prime, 4.5 h for twelve;
- heap peak 1.80 GB and peak working set 2.1 GB, against about 4.4 GB projected in F7 for the 31-bit build with the old engine.

**Estimates** ([results/analysis/bond171_gate_report.md](results/analysis/bond171_gate_report.md), [results/analysis/bond171_biased_x.json](results/analysis/bond171_biased_x.json)):

| Quantity | 156 terms | 171 terms | Jensen 1999 |
| --- | --- | --- | --- |
| 2ν⊥ (μ₂₀μ₀₂/μ₀₁², unbiased) | 2.1937080 | 2.1937077 (spread 3e-7) | 2.193708(4) |
| p_c from μ₂₀μ₀₂/μ₀₁² (unbiased) | 0.6447001908 | 0.6447001905 | 0.644700185(5) |
| γ + 2ν⊥ (μ₂₀, unbiased) | 4.4714110 | 4.4714139 | 4.471425(15) |
| γ + 2ν⊥ (biased at 0.644700185) | 4.4714237 | 4.4714267 | |
| 2ν⊥ (biased at 0.644700185) | 2.1937016 | 2.1936998 | |

**Drift.** The unbiased 2ν⊥ estimate drifts by only about 2 × 10⁻⁸ per term from 156 to 171 terms. With the 161–171-term points included, the drift fit gives an error factor of 1.9× per 10 terms. The fit projects a 40× gain from 171 to 230 terms, unchanged from F5. γ + 2ν⊥ from μ₂₀ is still creeping up by about 1 × 10⁻⁷ per term, as Jensen also found.

**One tension, inside Jensen's error bars.** The μ₂₀μ₀₂/μ₀₁² ratio puts the critical point at 0.6447001905, about 5 × 10⁻⁹ above the p_c from S and μ₀₂/μ₀₁ (0.6447001847–0.6447001855). That is at the upper edge of Jensen's 0.644700185(5).
- The 2ν⊥ estimate is steep in the assumed p_c: about 1.9 × 10⁻⁵ per 10⁻⁸.
- Biased at Jensen's central p_c, 2ν⊥ = 2.1936998, 2.0 of his error bars low.
- Unbiased, 2ν⊥ = 2.1937077, in agreement with his value.
- So the 2ν⊥ value depends on which critical point it is paired with. Jensen's quoted 2ν⊥ evidently used the series' own critical point.
- Consistency check: biased at 0.644700185, (γ + 2ν⊥) − γ = 4.4714267 − 2.2777276 = 2.1936991, which matches the biased R2 value. So the two routes agree once the same p_c is used.
- The 5 × 10⁻⁹ gap between the two series' critical points is the most useful target for longer series: shrinking it pins p_c and both ν's together.

**Conclusion.** All four 1999 series are reproduced at full length, and every unbiased estimate agrees with Jensen's within his quoted errors. That finishes the reproduction goal for the bond lattice (plan goal 1). The site series was completed to 158 next (F13).

## F13. Site series at Jensen's own length, 158 terms

[results/series/site158.json](results/series/site158.json) holds all four site series (S, μ₀₁, μ₀₂, μ₂₀) to order 158, the length of Jensen's 1999 site series.

**Verification.**
- Fourteen 15-bit primes (u16 build, [analysis/site158_queue.sh](analysis/site158_queue.sh)) give a product above 2²¹⁰; the largest coefficient has 177 bits.
- The leave-last-out check passes (crt.py).
- All four order-146 site series (F8), computed with 31-bit primes by the older engine, are exact prefixes.

**Cost** (all moments, `--threads 6 --block 6`):
- Heap peak 1.34 GB per prime, the same for every prime.
- 390–400 s per prime with the machine otherwise idle, and 800–890 s while the other job was running.
- 1 h 48 min for all fourteen primes.

**Unbiased estimates** ([results/analysis/site158_gate_report.md](results/analysis/site158_gate_report.md), K = 2 and 3, same settings as the bond analysis):

| Quantity | This work, T = 146 | This work, T = 158 | Jensen site, T = 158 | Bond value (J99) |
| --- | --- | --- | --- | --- |
| x_c (S, μ₀₂/μ₀₁) | 0.7054851–0.7054853 | 0.7054851–0.7054852 | 0.70548515(20) | |
| γ (S) | 2.27762 | 2.27758–2.27765 | 2.27765(6) | 2.277730 |
| ν∥ (μ₀₂/μ₀₁) | 1.73380 | 1.73380–1.73382 | 1.73381(4) | 1.733847 |
| 2ν⊥ (μ₂₀μ₀₂/μ₀₁²) | 2.19377 | 2.19376–2.19378 | 2.19377(5) | 2.193708 |
| γ + 2ν⊥ (μ₂₀) | 4.47124 | 4.47124–4.47130 | | 4.471425 |

Every site estimate agrees with Jensen's 158-term site values. The site exponents are still visibly drifting toward the bond values, as he noted. Measured from the bond values at 158 terms:
- γ is 1.2 × 10⁻⁴ low and has nearly stopped moving;
- ν∥ is 3.7 × 10⁻⁵ low and closes by about 1.3 × 10⁻⁶ per term;
- 2ν⊥ is 6.1 × 10⁻⁵ high and closes by about 4 × 10⁻⁷ per term.

The drift fits give error factors of 1.4× (ν∥) and 1.6× (2ν⊥) per 10 terms, against 1.9× for bond 2ν⊥ (F12). The site lattice is the slower route to the exponents.

**The site critical point from universality** ([results/analysis/site158_biased.json](results/analysis/site158_biased.json)): biased approximants, scanning x_c and locating where each site exponent equals the bond value:

| Series | Bond exponent imposed | 136 terms | 146 terms | 158 terms |
| --- | --- | --- | --- | --- |
| μ₀₂/μ₀₁ | ν∥ = 1.733847 | 0.70548523 | 0.70548523 | 0.70548523 |
| μ₂₀μ₀₂/μ₀₁² | 2ν⊥ = 2.193708 | 0.70548516 | 0.70548519 | 0.70548520 |

- Mean at 158 terms: 0.70548521 with a half-spread of 1.5 × 10⁻⁸. Jensen's biased estimate is p_c(site) = 0.70548522(4).
- As at 146 terms, the biased S and μ₂₀ exponents peak just below the bond values (2.277715 and 4.471385, both near x_c = 0.7054852–0.7054853) and do not cross them, so they are not used.

**Same pattern as bond.** The unbiased critical point from μ₂₀μ₀₂/μ₀₁² (0.70548526) sits above the ones from S and μ₀₂/μ₀₁ (0.70548515–0.70548518), the same ordering as the 5 × 10⁻⁹ gap on the bond lattice (F12), but about 20× wider here.

**Conclusion.** Both of Jensen's 1999 series sets are now reproduced at full length (bond 171, site 158), with every estimate inside his quoted errors. Plan goal 1 (an independent, public reproduction) is complete. The series files are in [results/series/](results/series/).
