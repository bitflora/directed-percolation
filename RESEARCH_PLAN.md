# Research plan: extending the 1+1 directed percolation series

Status, morning of November 22, 2026:
- M0 and the M2 gate are done.
- M1 is done to order 171 (bond, all moments) and 146 (site). (Updated November 23.)
- M3 is partly done.
- M6 has started. Results are in [FINDINGS.md](FINDINGS.md). The plan follows the prior-art assessment in [REVIEW.md](REVIEW.md); source IDs refer to [SOURCES.md](SOURCES.md).

| Milestone | Status |
| --- | --- |
| M0 ground truth | Done: three independent exact methods agree to order 15 (F1) |
| M1 reimplementation | Done in [enum/](enum/README.md). Verified bond series to 171 terms, all four moments (Jensen's full 1999 set); site to 146 terms. At 171 terms every unbiased estimate agrees with Jensen's within his errors (FINDINGS F11, F12). Still to run: site at 158 |
| M2 analysis and gate | Done. **Gate passed**: ν∥ and ν⊥ improve about 1.6–2× per 10 terms (FINDINGS F5) |
| M3 engineering | Partly done (F2, F3): shape, mirror, zero-state and window-trimming reductions; mask bounds and exact lower ends for single-source targets; u32 and u16 builds, streaming and block-parallel steps. Not yet done: SIMD, out-of-core, Zbarsky sparsification, and measuring parallel efficiency on an idle machine |
| M4 compute decision | Open; measured projections in FINDINGS |
| M5 | Not started (wall series for τ₁) |
| M6 | Started: D-finiteness search (no ODE with (k + 1)(d + 1) ≤ 147); exact ring rational functions up to L = 10 (F6, F9) |
| M7 | FINDINGS.md kept current; no paper draft yet |

## Goals, in priority order

1. Reproduce Jensen's 1999 square-lattice series exactly (bond 171, site 158) with an independent, public implementation.
2. Extend the bond series, then the site series, far enough to gain at least one significant digit on p_c, ν∥ and ν⊥.
3. Settle whether the wall exponent τ₁ equals 1, using a direct series for the average cluster length near a wall (S02).
4. Test proposed exact forms, and test whether the series are D-finite, at the new precision. Report negative results plainly.
5. Attempt a proof that a width-restricted DP generating function is not D-finite, using haruspicy (S16). This is independent of the long series and can run in parallel with M1–M3.

A milestone is complete only when its artifacts are in this folder and its numbers are written down. The two decision points are M2, a go/no-go on whether longer series are worth it, and M4, the compute spend.

## Working estimates (to be replaced by measurements)

Memory scales as M(n) ≈ 5.5 GB × 2^((n − 171)/8) (S01). This ignores the slow growth in coefficients per state, which makes it optimistic.

| Memory available | Order reachable | Setting |
| --- | --- | --- |
| about 11 GB | about 179 | This PC (16 GB RAM), in memory |
| about 160 GB | about 210 | This PC with out-of-core streaming to the free SSD space, if run time allows |
| 1 TB | about 231 | Cloud high-memory node |
| 4 TB | about 247 | Largest common cloud nodes |

Run time grows at least as fast as the state count, and is multiplied by the number of primes. Up to 12 primes of 15 bits were needed at order 171, and a few more will be needed at order 230. Time, not memory, may bind on the local machine.

## M0: Ground truth for small orders

- **Brute force.** Write a brute-force enumerator: for t ≤ 10 or so, sum over all bond (or site) configurations in the backward light cone of (x,t), and compute C(x,t) as an exact polynomial in p.
- **Simple transfer matrix.** Write the plain Domany–Kinzel transfer matrix of S01 section 2, which costs 2^(t/2) states. Use it to compute S(t) and X(t) exactly to about t = 30.
- **Assembly check.** Form the moment series and the non-nodal assembly formulas (REVIEW section 3). Check the first 20–40 terms against each other and against published low-order coefficients (S10, S11, S08, S06 tables; OEIS if present).
- **Update rules.** Re-derive the path-union update rules of S01 figure 3 and equation 12, which did not survive text extraction. Confirm them against brute force on every move type.

Deliverable: a test suite of exact small-order polynomials that every later implementation must match.

## M1: Reimplement the 1999 algorithm

The code is a Rust crate at `directed-percolation/enum/`; Rust is installed on this machine. Components:

- Sparse state hash keyed by a compact boundary encoding, with truncated polynomial coefficients per state.
- Non-nodal zeroing at horizontal boundary positions, the pruning rules (width and N_cont), and per-row coefficient trimming.
- Arithmetic modulo primes below 2¹⁵, or below 2³¹ if that proves no worse in bytes per bit. There is one run per prime, with CRT recombination and a verification prime that is not used in the reconstruction.
- Bond and site variants sharing one engine.
- Instrumentation: number of live states per row, peak resident memory, and wall time, logged to CSV.
- Also log, per row, the histogram of N_cont over stored states and the number of stored coefficients that can never reach order n. This measures the per-state window saving of REVIEW section 10.1 before any redesign.

Targets:

- Bond to 171 and site to 158. The published coefficients are not online, so the check is indirect:
  - agreement with the M0 suite at low order;
  - stability of the CRT reconstruction (an extra prime must agree);
  - recovery of the S01 critical estimates in M2.
- If the user decides to ask Jensen for the original files, an exact coefficient comparison becomes possible. That is an outward-facing step for the user to take, not something this plan does.
- Fit the configuration-growth factor and the memory and time laws against n, and compare with S01 figure 4 (about 2^(1/8)).

## M2: Analysis toolchain and the precision-gain gate

The code is a Python package at `directed-percolation/analysis/`, with a venv holding numpy, mpmath and sympy. Python 3.10 is present but currently has none of these.

- Implement inhomogeneous differential approximants of order K = 1–4, with varied L and the S01 rule that the Q_i degrees differ by at most 1. Solve the defining linear systems exactly, over the rationals or modulo primes with reconstruction, before any floating-point root finding; the series have 170+ terms with very large integer coefficients. Use mpmath only for the roots and the indicial equation.
- Add biased approximants that fix the non-physical singularities of S01 table 1, at their located positions and with the simple exponents where those look exact. Also try a conformal change of variable that moves the complex pair outside the physical disc. Compare both with the unbiased analysis.
- Reproduce the S01 analysis: unbiased p_c and exponent estimates, plots against number of terms used, the biased site analysis, the confluent-exponent checks (Baker–Hunter, Adler–Moshe–Privman), and the non-physical singularity table. The S01 values should come out within their error bars.
- **Gate.** Analyze the bond series truncated at 112, 130, 150 and 171 terms. Measure how the spread of the estimates shrinks with n, fit a law (for example, error ∝ e^(−an)), and project to 200, 230 and 250 terms. Cross-check the projection with Guttmann's coefficient-prediction approach (M02).
  - Continue to M3 only if 230 terms would give at least a 3× gain in ν∥ and ν⊥, or 200 terms at least 2×.
  - Otherwise, stop and write up a reproduction, the public series release, and the wall and D-finiteness tests on 171 terms.

## M3: Algorithm engineering on the local machine

These are ordered by expected payoff per effort. Measure each change against the M1 baseline.

1. **Memory per state.**
   - Compact boundary keys.
   - Open addressing with no per-entry pointers.
   - Residue width chosen for minimum bytes per useful bit.
   - Storing only the live coefficient window.
2. **Exact contribution orders.** Following Clisby and Jensen (M06, M09), compute for each configuration the exact minimum number of bonds needed to complete it into a non-nodal path-union graph ending at row t, instead of the S01 bound N_cont. Reformulating the state as the required future connectivity may make this O(width). Use it both to prune and to give each state its own coefficient window. Do this first: it has the best expected gain per unit of effort, and the M1 histograms say in advance how large the gain can be.
3. **Parallelism.** Six cores and 12 threads. Partition the state table by the boundary bits away from the update kink, following M05.
4. **Out-of-core.** Stream sorted state partitions from the SSD so that memory holds only the working partitions. Measure throughput; the 167 GB free sets the ceiling.
5. **Vectorization.** Rewrite the coefficient update as SIMD shifted-add kernels with AVX2, and choose between 16-bit lanes and 64-bit Barrett or Montgomery arithmetic by measured bytes per useful bit and throughput.
6. **Research option.** Determine whether a Zbarsky-style sparse-cut inclusion–exclusion (M04) can apply to path-union graphs.
   - First measure, from the M1 run, the distribution of occupied boundary sites per row over contributing configurations.
   - Then check on paper whether inclusion–exclusion over cut classes commutes with the ±1 path-union weights and the non-nodal zeroing.
   - Pursue an implementation only if both are favourable.

Deliverable: the measured order reachable on this machine, the per-order time, and updated projections for larger machines.

## M4: Compute decision (the user decides)

Present a table of order reached, estimated wall time and cost for:

- local out-of-core;
- a cloud high-memory node of about 1–4 TB;
- a university or institutional cluster, if the user has access.

No paid or external resources are used without explicit approval.

**State after the first night (November 22, 2026).** The table is in FINDINGS F7. In short:
- The gate passed: ν∥ and ν⊥ improve about 1.8× per 10 extra terms.
- Order 171 is feasible on this PC in about 16–30 h of compute.
- Order 200 needs a node with about 128 GB.
- Order 230 needs about 1 TB.

Recommended next steps, cheapest first:
1. Parallel efficiency is about 75% at 4 threads and 60% at 6 on the laptop (FINDINGS F7). Measure it on a many-core server before renting one for long.
2. Order 171 is done for all four series (F11, F12; μ₂₀ took twelve 15-bit primes at about 23 min each and 2.1 GB). An exact coefficient check still needs his series files.
3. Only then decide on a 128 GB to 1 TB node for 200–230 terms.

## M5: Production runs and analysis

- **Bond series** to the agreed order, then **site**.
- **Wall series** (S02 geometry): compute the moments needed for γ₁ and ν∥,1, and a direct series for the average cluster length that gives τ₁ without scaling relations.
- **Optional:** the triangular lattice (S03) as a universality cross-check, if the engine generalizes cheaply.
- **Optional independent check:** exact diagonalization of the 1D contact-process generator, using translation and reflection symmetry and a sparse eigensolver, for L up to about 28–30. Extrapolate z and ν⊥ by finite-size scaling with Bulirsch–Stoer (as in N01), and compare with the series values. This is feasible on this machine and has different systematic errors from series.
- **Analysis:** repeat the M2 analysis on every series.
  - Quote final values with error bars from the spread of approximants, and say plainly that this is the convention.
  - Check the scaling relations (β from different combinations) and bond–site universality.

## M6: Exact-form and solvability tests

- Re-test the S08 rationals, Sire's β (N08) and τ₁ = 1 at the new precision.
- **Simple-form search.** Search small-denominator rationals and quadratic irrationals over a declared class, sized in advance. Report the expected number of chance matches at the achieved precision (the look-elsewhere correction), and claim a match only if it is far beyond that.
- **Non-physical singularities.** Test whether their exponents are exactly the simple values S01 suggests. These are a more plausible place for exact structure than the physical exponents.
- **D-finiteness.**
  - Guess linear ODEs with polynomial coefficients modulo single primes, directly from the per-prime enumeration residues (M12, M13), and validate any candidate on held-out terms and on other primes.
  - Use the order–degree trade-off (M14) to reach the largest excludable class for the available length.
  - Do the same for algebraic and D-algebraic relations (M07).
  - Without a proof, record the largest ODE order and degree excluded.
- **Haruspicy on strips (S16).** This can start any time after M0, because it needs only small transfer matrices.
  - Compute exactly the rational generating functions of the DP pair connectedness (or survival probability) on strips of width w = 1, 2, …, as far as exact arithmetic allows.
  - Track how their poles move and accumulate as w grows.
  - If they accumulate densely on a curve, attempt Rechnitzer's argument to prove that the corresponding unrestricted generating function is not D-finite.
  - Report a proof, a conjecture with evidence, or an obstruction, whichever results.

## M7: Write-up

- Write FINDINGS.md, recording each claim with its evidence.
- Release the series files and the code in this folder.
- Draft a paper for J. Phys. A or Phys. Rev. E. Lead with the new series and estimates, then the τ₁ verdict, then the exact-form and D-finiteness exclusions, then the algorithmic changes.
- Before submission, read N04 and N05 in full and check them for competing estimates.

## Risks

- **Reproduction drift.** A subtle error in the update rules can produce plausible but wrong high-order terms. M0's brute-force suite and the extra verification prime are the defences.
- **No precision gain.** The DA estimates in S01 were still drifting (μ₀₁ and μ₀₂). More terms may expose drift rather than shrink errors. The M2 gate exists for this.
- **Time.** Out-of-core may be memory-feasible but take months. Measure before committing.
- **Overclaiming.** Exact-form "matches" at 7 digits are cheap. Use the look-elsewhere accounting.
