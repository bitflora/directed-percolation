# Findings

Running record of results for the directed percolation study. Each claim lists its evidence. Numbers marked *provisional* come from runs or analyses still in progress.

## Summary of the first night (November 21–22, 2026)

- **Engine.** An independent reimplementation of Jensen's 1999 algorithm, checked against exact ground truth (F1, F2). It adds translation and mirror reductions, exact window trimming and faster bounds (F3). Per prime it needs roughly 6× less memory than the 1999 implementation.
- **Verified exact series:**
  - bond to order 156 for all four moment series (S, μ₀₁, μ₀₂, μ₂₀) (F4);
  - bond to order 166 for S, μ₀₁ and μ₀₂ (F10);
  - site to order 136 (F8).

  Jensen reached 171 (bond) and 158 (site).
- **Reproduction of the 1999 results.**
  - At Jensen's p_c, biased approximants give γ, ν∥, 2ν⊥ and γ + 2ν⊥ within his error bars.
  - The ν∥ and 2ν⊥ estimates are stable from 146 to 156 terms.
  - The site critical point from universality is 0.70548522, as he found (F5, F8).
- **Gate (M2) passed.** Estimates of ν∥ and ν⊥ improve about 1.8–2× per 10 extra terms, so 200–230 terms should give about one more significant digit (F5).
- **Solvability.** There is no small linear ODE for any of the series (F6). On rings, the poles of the mean cluster size pinch toward p_c and cluster around the bulk non-physical singularities: a first look at the haruspicy route (F9).
- **Compute decision (M4).** Order 171 is feasible on this PC. Orders 200 and above need a large-memory machine; projections are in F7.

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

**Untested at scale:**
- The sharded multithreaded mode has been verified only for correctness (order 15 and 40). Its parallel efficiency without contention has not been measured, because this machine was fully loaded all night.
- Measuring it on an idle machine is the first step before renting anything.

## F8. Site series to order 136, and the site critical point from universality

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

The computation is cheap: L = 10 took under ten minutes. So M6 can continue to L ≈ 14–16.

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
