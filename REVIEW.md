# Literature review: directed percolation in 1+1 dimensions

Search date: November 21, 2026.

Source identifiers refer to [the annotated bibliography](SOURCES.md), where inspection depth is recorded. Only Jensen 1999 (S01) was read in full. The other entries were checked at abstract, HTML-summary or metadata level. Where this review proposes an opening, that is an inference from the retrieved literature, not a published open-problem statement.

## Findings that change the original proposal

**The 27-year-old record is real.** Jensen's 1999 low-density series for the directed square lattice still define the state of the art: bond to order 171 and site to order 158 (S01). A Semantic Scholar sweep of every paper citing S01, through 2026, found no extension. Targeted searches for newer series work also found none. Jensen's later papers extended other cases: a wall (S02, 173 terms), other lattices (S03) and temporal disorder (S04). Monte Carlo work in 1+1 does not compete; the main high-precision study covers only d ≥ 2 (N03).

**The "simple exact form" angle is mostly closed already.** The one published rational conjecture for the bulk exponents (S08, 1988) is excluded by the 1999 numbers themselves, by 10 to 90 standard errors (table below). The best-known closed-form guess for β (N08) is an approximation that misses by about 12 standard errors. No newer bulk conjecture was found. Longer series would therefore mostly tighten numbers rather than test live hypotheses. A positive "hidden solvability" result is a long shot.

**A sharper near-exact question survives at a wall.** Near an impenetrable wall, scaling relations give τ₁ = 1.00014(2), which is 7 standard errors from the conjectured τ₁ = 1 (S02). Jensen called this tantalizingly close and asked for a direct series for the average cluster length. That calculation uses the same machinery, needs no new lattice, and is the cleanest yes/no exact-value test available.

**The experimental motivation in the idea note needs correcting.** The liquid-crystal turbulence experiments (E05, E06) are (2+1)-dimensional DP, not 1+1. The 1+1 realizations are quasi-one-dimensional Couette flow (E01, E02) and pipe flow (E03, E07). Their precision is about two digits, for example β = 0.29 ± 0.05 (E02). No experiment needs a seventh digit. The real consumers of better numbers are approximate theories and numerical methods (section 6).

**The user's hardware is the binding constraint.** The 1999 run needed about 5.5 GB for order 171, and memory grows by a factor of 2 for every 8 extra terms (S01). The 16 GB machine available here reaches only about order 179 in RAM. A publishable extension needs algorithmic gains, out-of-core streaming, or rented memory; see [RESEARCH_PLAN.md](RESEARCH_PLAN.md).

## 1. The problem

Bond DP on the directed square lattice is equivalent to the Domany–Kinzel cellular automaton with W(0|l,r) = (1 − p)^(l+r) (S18). Above p_c an infinite directed cluster exists. Below p_c the moments of the pair connectedness C(x,t),

μ_{n,m}(p) = Σ_t Σ_x xⁿ tᵐ C(x,t) ~ (p_c − p)^−(γ + nν⊥ + mν∥),

diverge at p_c. The average cluster size is S = μ₀₀. The Janssen–Grassberger conjecture (T01, T02) says any transition into a single absorbing state, with short-range interactions and no extra symmetry or conservation law, is in this universality class. The class includes epidemics without immunity, the contact process, catalytic surface reactions, and the laminar–turbulent transition in shear flows (R01, R02).

In 1+1 dimensions no exponent is known exactly (R01, T03, T08). This is unlike isotropic 2D percolation, where conformal field theory gives exact rational exponents. DP is anisotropic (z ≠ 1), so ordinary conformal methods do not apply. Known exact results concern neighbouring problems: directed animals, which are equivalent to hard squares (S17); compact DP on a special line of the Domany–Kinzel model (S18); exact correction-term structure in the percolation-probability series (S06, S13, S14, S15); and exact survival probabilities at unphysical complex p (N07).

## 2. History of the low-density series record (square bond)

| Year | Terms | Method | Source |
| --- | --- | --- | --- |
| 1977 | 21 | Transfer matrix, growth √2 per term | S10 |
| 1983–1986 | 25, 35 | Transfer matrix plus weak-subgraph expansion | S11 |
| 1988 | 49 | Non-nodal graph resummation (doubles the order), growth 2^(1/4) | S08 |
| 1996 | 112 | t_max = 47 plus a correction-term extrapolation procedure | S05 |
| 1999 | 171 (site 158) | Direct non-nodal path-union transfer matrix, growth ≈ 2^(1/8) | S01 |

High-density series for the percolation probability P(q) form a separate line (S09, S06, S07). That line is shorter, but it gives β directly instead of through a scaling relation.

## 3. The 1999 algorithm, as needed for reimplementation

From S01, sections 2 and 3:

- **Expansion.** C(x,t) = Σ_g d(g) p^|g| over graphs that are unions of directed paths from the origin to (x,t) (S12). The weight d(g) is (−1)^(number of path splits), so d(g) = ±1, and graphs with dangling ends have weight zero.
- **Non-nodal reduction.** A graph is nodal if every path passes through some intermediate point. Computing only the non-nodal parts S_N, μᴺ₀₁, μᴺ₀₂ and μᴺ₂₀ to t_max ≈ n/2 is enough. The full moments follow from

  S = 1/(1 − S_N), μ₀₁ = μᴺ₀₁ S², μ₀₂ = [μᴺ₀₂ + 2(μᴺ₀₁)² S] S², μ₂₀ = μᴺ₂₀ S².

  Nodal configurations are removed by zeroing the states with a single incoming bond (encoded as integers 2^k) whenever the boundary is horizontal.
- **Transfer matrix.** A boundary line is moved through the lattice one site at a time. Each site has two states: whether it has an incoming bond from the row above (bond) or whether it is occupied (site). Updates involve only the kink site and its two neighbours below, so they are local and parallelize naturally. The bond update rules are given pictorially in S01 figure 3. The site rules (S01 equation 12) are the bond rules with the weight divided by p wherever a bond enters a site that is already occupied. The figure symbols did not survive text extraction, so M0 re-derives the rules and checks them against brute force.
- **Pruning.**
  - A configuration whose outermost occupied sites are r apart needs at least r − 1 more rows to collapse, so it is dropped if t′ + r > t_max.
  - It is also dropped if its minimum contributing order N_cont = 2N_min + 2t′ − 4t (bond) or 2N_min + 2t′ − N_occ − 4t (site) exceeds n.
  - Polynomials store only the n − 2t coefficients that can be nonzero, and the two lowest are discarded at each new row.
- **Storage and arithmetic.**
  - Sparse states live in a hash table.
  - Coefficients are computed modulo primes p_i = 2¹⁵ − r_i, stored as 16-bit residues, and recombined with the Chinese remainder theorem. Up to 12 primes were needed at order 171.
  - Intermediate values need not be representable, only the final coefficients.
- **Cost.** The number of configurations grows as about 2^(n/8) ≈ 1.091ⁿ, against 2^(n/4) for the 1996 method. Doubling memory buys about 8 more terms. Order 171 took about 5.5 GB.

The paper reports neither run time nor the number of separate runs (one per prime is implied by the 16-bit storage).

## 4. Current best values in 1+1 dimensions

From S01 unless noted. The bond series is the precise one; site estimates are an order of magnitude less precise and only marginally consistent with bond.

| Quantity | Value | Note |
| --- | --- | --- |
| p_c, square bond | 0.644700185(5) | |
| p_c, square site | 0.70548522(4) | Biased by assuming the bond exponents; unbiased value 0.70548515(20) |
| γ | 2.277730(5) | |
| ν∥ | 1.733847(6) | |
| ν⊥ | 1.096854(4) | From 2ν⊥ = 2.193708(4) |
| β | 0.276486(8) | From the scaling relation β = (ν∥ + 2ν⊥ − γ)/2, not measured directly |
| z = ν∥/ν⊥ | 1.580745(10) | Derived |
| δ = β/ν∥ | 0.159464(6) | Derived |
| Δ₁ | consistent with 1 | No sign of non-analytic corrections |
| λ_c, contact process | 1.64896(2) | Exact diagonalization (N01) |
| Wall exponents | ν∥,1 = 1.733820(12), γ₁ = 1.82051(1), τ₁ = 1.00014(2) | S02 |
| p_c, other lattices | triangular bond 0.47802525(5), kagomé bond 0.65896910(8), (4.8²) bond 0.76765116(7) | S03, from a summary; verify |

The 1999 analysis also located non-physical singularities: a negative-axis singularity at −0.51670(3) for bond and −0.451952165(8) for site, plus complex-conjugate pairs. The site singularities carry exponents consistent with simple exact values (1/2, −1/2, −3/2 at the negative-axis point; 3, 2, 1 at the complex pair). The bond negative-axis point suggests a logarithmic divergence in S. Such exact non-physical exponents are a reasonable place to look for hidden structure.

## 5. Other numerical routes and their precision

- **Exact diagonalization and DMRG of the generator** (N01, N02). These give z and ν⊥ to about 5 digits, within 1 to 2 of their own error bars of the series values. They are limited by system size.
- **Strip transfer matrices and Yang–Lee zeros** (N04, N05, N06). Recent and precise, but the numbers were not retrieved this session. They are the most likely recent competitor for ν⊥ and should be checked first.
- **Tensor networks** (N10). These measure entropy-related exponents, not the standard set.
- **Monte Carlo.** No 1+1 study competes with series. N03 covers d ≥ 2; N09 covers the parity-conserving class, a different universality class.
- **Machine learning.** Several 2023–2026 papers learn DP transitions. None is a precision benchmark.

## 6. Field theory and bootstrap: who uses better numbers

- **ε-expansion around d_c = 4.** Now complete to three loops (T03), with a 2026 extension to scaling functions (T04). At d = 1 the expansion parameter is ε = 3, so raw truncations are useless. Testing resummations needs precise 1+1 targets.
- **Non-perturbative or functional RG** (T05, T06, T07). Reported to compare accurately at d = 1 at LPA and LPA′; the numbers still need to be tabulated.
- **1D-specific analytic schemes**, for example a fermionic field theory (T09).
- **Positivity bootstrap** (T08). It gives rigorous but currently loose bounds: contact process λ_c ≥ 1.46 against 1.649. It is a long-term consumer, and the gap shows how far rigorous methods are from series precision.
- **Local scale invariance tests** (T10). Their verdict for the 1D contact process was limited by exponent precision.
- **Quantum-hardware simulations of DP dynamics** (E08). These use the classical values as the reference.

## 7. Status of exact-value proposals

Distances below use the S01 error bars. Several of these entries are approximations, not proposals of exact values; they are included because they are the forms most likely to be cited.

| Proposal | Value | Best estimate | Distance |
| --- | --- | --- | --- |
| γ = 41/18 (S08) | 2.277778 | 2.277730(5) | ≈ 10σ |
| ν⊥ = 79/72 (S08) | 1.097222 | 1.096854(4) | ≈ 92σ |
| ν∥ = 26/15 (S08) | 1.733333 | 1.733847(6) | ≈ 86σ |
| β = (1 − 1/√5)/2 (N08, approximation) | 0.276393 | 0.276486(8) | ≈ 12σ |
| p_c = (√5 − 1)/2 (N08, approximation) | 0.618034 | 0.644700185(5) | far |
| Wall τ₁ = 1 (S02) | 1 | 1.00014(2) | ≈ 7σ, from scaling relations only |

The golden-ratio results of N07 are exact, but they hold at unphysical p and do not constrain exponents. The 2005 "exact solution" (N11) concerns a different row-alternating model.

Two structural results bear on solvability. The directed bond-animal generating function is proven non-D-finite (S16), and the DP generating function is not believed to be D-finite. No paper found tests the DP moment series themselves for D-finiteness or D-algebraicity with modern guessing tools (M07, M08). With 171 or more exact terms, such a test is cheap. A negative result is a modest, publishable observation, and a positive one would be a major one.

The main methodological risk for any exact-form claim is look-elsewhere: at 6–7 digits, many "simple" closed forms will fit some exponent by chance. N09 is a recent example of how near-rational claims are argued in a sister class.

## 8. Experiments

- **1+1 realizations.**
  - Circular Couette flow (E01, E02): β = 0.29 ± 0.05.
  - Pipe flow, where puffs behave as DP with a jammed phase above criticality (E03, E07).
  - Trapped-ion quantum simulations of a Floquet contact process (E08).
- **2+1 realizations**, not relevant to the 1+1 numbers: turbulent liquid crystals (E05, E06) and channel flow (E04).

Experiments confirm the universality class at the 10% level. They are not a reason to want more digits, and the write-up should not claim they are.

## 9. Novelty assessment

| Contribution | Status | Expected value |
| --- | --- | --- |
| Reimplementation and public release of the 171/158-term series | Coefficients were only ever distributed by e-mail from a defunct web page (S04) | Modest. Enables independent reanalysis. |
| Extension of the square bond and site series | Unchallenged since 1999 | The core result. About one more digit if 220–250 terms are reached; less if only about 200. |
| Direct wall series for the average cluster length | Requested in S02, apparently never done | A clean verdict on τ₁ = 1. Small but sharp. |
| Systematic exact-form and D-finiteness tests | Bulk rational forms already excluded; no guessing study found | Low cost. Negative results are publishable as part of the main paper. |
| Algorithmic improvement (per-configuration contribution orders, Zbarsky-type sparsification M04, out-of-core) | Open | Potentially the most original part, if any of it works. See section 10. |
| Proof that a width-restricted DP generating function is not D-finite (haruspicy, S16) | Done for directed bond animals, apparently not for DP | A theorem rather than a numerical hint about solvability. Possibly the most interesting mathematical outcome. |

Realistic outcome: a solid J. Phys. A or Phys. Rev. E paper combining the extended series, the wall test and a careful exact-form exclusion. A positive exact form for a bulk exponent is unlikely, given section 7.

## 10. Techniques developed since 1999

This section lists methods that did not exist in 1999, or were not used on DP then, ranked by expected payoff. Statements marked *inference* are this review's reasoning, not published results, and each is tied to a measurement in [RESEARCH_PLAN.md](RESEARCH_PLAN.md).

### 10.1 Enumeration

**Exact per-configuration contribution orders (M06, M09).** This is the most likely real gain.
- Clisby and Jensen rebuilt their polygon transfer matrix so that each state records how its boundary edges must be connected in the future, not how they were connected in the past. The cheapest completion of each configuration is then unique and costs O(W) to compute. This sharpened pruning, and polygons went from 110 to 130 steps.
- S01 uses a coarser bound for DP, N_cont = 2N_min + 2t′ − 4t, and stores the same n − 2t coefficients for every state in a row.
- *Inference:* an exact minimum completion cost for a path-union configuration, given row, width and the positions of occupied sites, would prune more states. It would also let each state keep only the coefficients that can still reach order n, a window of about n − N_cont + 1 instead of n − 2t.
- Memory is dominated by coefficients, not keys. The number of states grows exponentially with order, so most stored states lie near the cutoff and would need very few coefficients. Each factor of 2 saved buys about 8 terms.
- The size of the gain is unknown until the M1 code reports the distribution of N_cont over stored states.

**Sparse-cut inclusion–exclusion (M04).** This is the only idea found that could change the asymptotic cost.
- For self-avoiding walks and polygons, Zbarsky reduces the cost from exp(Cn) to exp(C√(n log n)). The argument is a pigeonhole one: some residue class of evenly spaced cut lines must be crossed only a few times. Inclusion–exclusion over those classes then restricts the transfer matrix to narrow states at the cuts.
- *Inference:* a contributing DP graph at order n has about n bonds spread over about n/2 rows, so most rows are crossed by few paths, which is the sparsity the argument needs.
- What is unknown is whether the inclusion–exclusion bookkeeping is compatible with the ±1 path-union weights and the non-nodal reduction. If it is, this would be the most original part of a paper.

**Hardware and engineering.**
- The core update is a shifted add of residue arrays, which vectorizes well with AVX2.
- 64-bit modular arithmetic with Barrett or Montgomery reduction is cheap.
- NVMe streaming makes out-of-core sweeps practical, following the M05 partitioning.
- Cloud nodes offer terabytes of RAM; S01 had 5.5 GB.
- *Inference:* per-core throughput is 50 to 100 times that of 1999 hardware, so memory rather than time should bind on the local machine. M1 measures this.

**Considered and set aside.**
- GPUs: the computation is memory-bound, and the local card has 4 GB.
- Frontier-based search with ZDDs: equivalent to the existing hash-table transfer matrix for weighted counting.
- The critical-polynomial method (M10, M11): it gives 10 to 15 digit thresholds for isotropic 2D percolation, but it relies on structure DP lacks. No DP application was found.

### 10.2 Series analysis

- **Exact or modular differential approximants.** The coefficients are integers, so the linear systems defining an approximant can be solved exactly (modulo primes, or over the rationals) before any floating-point step. This avoids loss of precision at 200+ terms and makes large approximant ensembles cheap. Fast Hermite–Padé (order-basis) algorithms suit this.
- **Coefficient prediction (M02).** Differential approximants predict further coefficients approximately, and the predicted terms improve ratio-type estimates. The extra gain is modest; it also helps project the benefit of longer series in M2.
- **Built-in non-physical singularities.** The radius of convergence is set by the complex pair at |p| ≈ 0.49, inside p_c ≈ 0.645 (S01 table 1). For site percolation their exponents look exactly simple. Approximants that fix these singularities, or a change of variable that moves them away, should sharpen the physical estimates. The idea predates 1999, but only became usable once the singularities were located to 4 to 9 digits.

### 10.3 Solvability tests

- **ODE guessing modulo a prime (M12, M13).** The Ising-susceptibility programme guessed a Fuchsian ODE of order 33 from a series computed modulo one prime, then lifted factors to exact arithmetic. DP enumeration produces per-prime residues anyway, so a mod-p guessing test needs no full reconstruction.
  - The limit is length: those computations used thousands of terms. With about 250 DP terms, only ODEs whose (order + 1) × (degree + 1) is below about 250 can be excluded, even after using the order–degree trade-off (M14).
  - D-algebraic guessing (M07) extends the excluded class.
- **Haruspicy (S16).** Rechnitzer proved the directed bond-animal generating function is not D-finite. He solved width-restricted versions exactly and showed that their singularities accumulate densely on a curve, which a D-finite function cannot have.
  - *Inference:* the same strategy may apply to a width-restricted DP generating function, for example the pair connectedness on strips, where the transfer matrix is finite and each strip generating function is rational.
  - If the poles of the strip functions can be controlled as the width grows, the result is a theorem that the corresponding DP generating function is not D-finite. That is stronger than any numerical exclusion, and the strip transfer matrices are small enough to experiment with on this machine.

### 10.4 Cross-checks and consumers

- **Exact diagonalization with modern sparse eigensolvers.** N01 used 18 sites in 1999 and obtained z and ν⊥ to about 5 digits. *Inference:* with translation and reflection symmetry, roughly 28 to 30 sites fit in 16 GB. That gives an independent check of z and ν⊥ on this machine, with a different source of systematic error from series.
- **Tensor networks** (TRG and its successors since 2007) are highly precise for equilibrium 2D models. They are hampered by DP's absorbing, non-Hermitian dynamics, and current DP results (N10) are not competitive.
- **The functional RG** at high derivative-expansion order now matches bootstrap precision for the 3D Ising model (T11). Applying the same programme to DP needs precise 1+1 targets. The positivity bootstrap (T08) is in the same position: a consumer of better values, not a competitor.
- **Monte Carlo** remains uncompetitive with series in 1+1.
