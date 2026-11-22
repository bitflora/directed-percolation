# Findings

Running record of results for the directed percolation study. Each claim lists its evidence. Numbers marked *provisional* come from runs or analyses still in progress.

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
