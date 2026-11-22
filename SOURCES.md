# Annotated bibliography for 1+1 directed percolation

Search date: November 21, 2026. The entries map the relevant literature; this is not an exhaustive citation census. "Inspection" records what was actually consulted. Only S01 was read in full. Entries marked "background knowledge" or "search snippet" were not opened this session, and their claims should be checked before they are cited in a paper. A link labelled as a Crossref search marks a DOI that could not be verified. [sources.json](sources.json) holds the same records in machine-readable form.

## Series expansions and exact combinatorics

### S01 Low-density series expansions for directed percolation I: A new efficient algorithm with applications to the square lattice

Iwan Jensen. 1999. J. Phys. A 32, 5233 (1999).

[Link](https://arxiv.org/abs/cond-mat/9906036). Inspection: full text (arXiv v1, extracted with pdftotext; figure symbols in the update rules did not survive extraction).

The current record and the reproduction target. Square bond series to order 171 and site series to order 158 for the moments of the pair connectedness. Non-nodal path-union transfer matrix, configuration growth about 2^(n/8), about 5.5 GB at order 171, residues modulo primes near 2^15 with up to 12 primes. Final estimates: bond p_c = 0.644700185(5), site p_c = 0.70548522(4), gamma = 2.277730(5), nu_par = 1.733847(6), nu_perp = 1.096854(4), beta = 0.276486(8) from scaling, Delta_1 consistent with 1.

### S02 Low-density series expansions for directed percolation II: The square lattice with a wall

Iwan Jensen. 1999. J. Phys. A 32 (1999); DOI listed as 10.1088/0305-4470/32/33/304.

[Link](https://arxiv.org/abs/cond-mat/9906037). Inspection: primary abstract and HTML text via fetch summary.

Series to order 173 near an impenetrable wall. Surface exponents nu_par,1 = 1.733820(12), nu_perp,1 = 1.0968(1), gamma_1 = 1.82051(1). Scaling gives tau_1 = 1.00014(2), which appears to exclude the conjecture tau_1 = 1; the paper calls the value tantalizingly close and asks for a direct series for the average cluster length. This is the sharpest surviving exact-value question found.

### S03 Low-density series expansions for directed percolation III. Some two-dimensional lattices

Iwan Jensen. 2004. J. Phys. A 37, 6899 (2004).

[Link](https://arxiv.org/abs/cond-mat/0405504). Inspection: primary abstract; table values from a fetch summary of the truncated HTML text.

Triangular bond 90 and site 82 terms, honeycomb bond 301, kagome 173/169, (4.8^2) 255/239. Exponents agree across lattices to the fifth digit. The square lattice was not extended. Reported p_c values (summary, not verified verbatim): triangular bond 0.47802525(5), kagome bond 0.65896910(8), (4.8^2) bond 0.76765116(7).

### S04 Low-density series expansions for directed percolation IV. Temporal disorder

Iwan Jensen. 2005. J. Phys. A (2005).

[Link](https://arxiv.org/abs/cond-mat/0509713). Inspection: search snippet.

Extends the same algorithm to temporally disordered DP. States that series were available by e-mail or from the author's Melbourne web page; that page is no longer reachable.

### S05 Low-density series expansions for directed percolation on square and triangular lattices

Iwan Jensen. 1996. J. Phys. A 29, 7013 (1996).

[Link](https://search.crossref.org/search/works?q=Jensen+Low-density+series+expansions+for+directed+percolation+on+square+and+triangular+lattices+1996&from_ui=yes). Inspection: metadata; described in S01 section 2; DOI not verified.

The previous record: transfer matrix to t_max = 47 plus a correction-term extension procedure, reaching 112 terms (111 for mu_20). Earlier estimate p_c = 0.64470015(5). Volume and page from S01 reference 27; the link is a Crossref search because the DOI could not be verified.

### S06 Series expansions of the percolation probability for directed square and honeycomb lattices

Iwan Jensen and Anthony J. Guttmann. 1995. J. Phys. A 28, 4813 (1995).

[Link](https://arxiv.org/abs/cond-mat/9509121). Inspection: primary abstract.

High-density series for the percolation probability P(q). Square bond q_c = 0.177143(2), beta = 0.2763(2). The first correction term is a Catalan number and higher corrections are rational functions of Catalan numbers, which drives the Baxter-Guttmann style extension. An independent observable for cross-checking beta.

### S07 Series expansions of the percolation probability on the directed triangular lattice

Iwan Jensen and Anthony J. Guttmann. 1996. J. Phys. A (1996).

[Link](https://arxiv.org/abs/cond-mat/9511084). Inspection: search snippet.

Triangular bond and site P(q) series with Catalan-type correction formulas. Lattice cross-check for universality.

### S08 On two-dimensional directed percolation

J. W. Essam, A. J. Guttmann and K. De'Bell. 1988. J. Phys. A 21, 3815 (1988).

[Link](https://iopscience.iop.org/article/10.1088/0305-4470/21/19/018). Inspection: primary abstract.

Introduced the non-nodal resummation that doubles the series length (49 terms) and conjectured gamma = 41/18, nu_perp = 79/72, nu_par = 26/15. The S01 estimates exclude all three by many error bars.

### S09 Series expansion of the percolation probability for the directed square lattice

R. J. Baxter and A. J. Guttmann. 1988. J. Phys. A 21, 3193 (1988).

[Link](https://search.crossref.org/search/works?q=Baxter+Guttmann+Series+expansion+of+the+percolation+probability+for+the+directed+square+lattice+1988&from_ui=yes). Inspection: metadata; cited as S01 reference 28.

Origin of the finite-lattice extrapolation idea used to extend series beyond the directly computed order. The DOI guessed earlier resolved to an unrelated paper, so the link is a Crossref search.

### S10 Series expansions for the directed-bond percolation problem

J. Blease. 1977. J. Phys. A 10, 917 and 3461 (1977).

[Link](https://search.crossref.org/search/works?q=Blease+directed+bond+percolation+series+1977+J.+Phys.+A+10+917&from_ui=yes). Inspection: metadata; cited as S01 reference 22.

First transfer-matrix approach, 21 low-density terms, growth sqrt(2) per term. Useful low-order check values. Title not verified; link is a Crossref search.

### S11 Analysis of extended series for bond percolation on the directed square lattice

J. W. Essam, K. De'Bell, J. Adler and F. M. Bhatti. 1986. Phys. Rev. B 33, 1982 (1986).

[Link](https://doi.org/10.1103/PhysRevB.33.1982). Inspection: metadata verified via Crossref; cited as S01 reference 24.

35-term series using a weak subgraph expansion. Low-order check values.

### S12 Percolation theory on directed graphs

D. K. Arrowsmith and J. W. Essam. 1977. J. Math. Phys. 18, 235 (1977).

[Link](https://doi.org/10.1063/1.523262). Inspection: metadata verified via Crossref; cited as S01 reference 29.

Source of the pair-connectedness expansion over unions of directed paths with weights d(g) = +-1, which is the basis of the S01 algorithm.

### S13 Catalan numbers in a series expansion of the directed percolation probability on a square lattice

N. Inui and M. Katori. 1996. J. Phys. A 29, 4347 (1996).

[Link](https://search.crossref.org/search/works?q=Inui+Katori+Catalan+numbers+series+expansion+directed+percolation+probability+square+lattice&from_ui=yes). Inspection: bibliographic lead from reference lists; DOI not verified (link is a Crossref search).

Exact structure in the correction terms of the P(q) series. Relevant to any claim of hidden exact structure, and possibly to checking computed coefficients.

### S14 Hypergeometric series in a series expansion of the directed-bond percolation probability on the square lattice

M. Katori, N. Inui, G. Komatsu and K. Kameoka. 1997. J. Stat. Phys. (1997).

[Link](https://link.springer.com/article/10.1007/BF02180198). Inspection: search result title.

Follow-up to S13: higher correction terms written as hypergeometric series.

### S15 New enumerative results on two-dimensional directed animals

Mireille Bousquet-Melou. 1998. Discrete Math. (1998).

[Link](https://www.sciencedirect.com/science/article/pii/S0012365X9700109X). Inspection: search snippet.

Directed-animal results used to derive leading correction terms in DP series. Exact combinatorics adjacent to DP.

### S16 Haruspicy 3: The directed bond-animal generating function is not D-finite

Andrew Rechnitzer. 2004. arXiv math/0408054.

[Link](https://arxiv.org/abs/math/0408054). Inspection: search snippet.

Proves the anisotropic directed bond-animal generating function is not D-finite, and so cannot be guessed by GFUN or differential approximants. Supports the expectation that DP is not D-finite; the method is a template for a non-D-finiteness argument.

### S17 Equivalence of the two-dimensional directed-site animal problem to Baxter's hard-square lattice-gas model

Deepak Dhar. 1982. Phys. Rev. Lett. 49, 959 (1982).

[Link](https://doi.org/10.1103/PhysRevLett.49.959). Inspection: metadata from background knowledge; not retrieved this session.

Exact solution of directed animals, the classic example of hidden solvability next to DP. Background only.

### S18 Equivalence of cellular automata to Ising models and directed percolation

E. Domany and W. Kinzel. 1984. Phys. Rev. Lett. 53, 311 (1984).

[Link](https://doi.org/10.1103/PhysRevLett.53.311). Inspection: metadata; cited as S01 reference 11.

Domany-Kinzel cellular automaton, which contains bond DP via W(0|l,r) = (1-p)^(l+r) and the exactly solvable compact-DP line. Basis of the simple transfer matrix used for M0 checks.

## Enumeration, series-analysis and guessing methods

### M01 Asymptotic analysis of power-series expansions

Anthony J. Guttmann. 1989. Phase Transitions and Critical Phenomena, Vol. 13, eds. C. Domb and J. L. Lebowitz (Academic Press, 1989).

[Link](https://www.worldcat.org/search?q=Phase+Transitions+and+Critical+Phenomena+volume+13). Inspection: metadata; cited as S01 reference 20.

Standard reference for differential approximants and ratio methods as used in S01.

### M02 Series extension: predicting approximate series coefficients from a finite number of exact coefficients

Anthony J. Guttmann. 2016. J. Phys. A 49 (2016).

[Link](https://arxiv.org/abs/1604.08254). Inspection: search snippet.

Differential approximants predict later coefficients well, and the predicted terms improve ratio-method estimates. Relevant to the M2 precision-gain projection and as a cheap complement to exact extension.

### M03 Analysis of series expansions for non-algebraic singularities

Anthony J. Guttmann. 2015. J. Phys. A 48 (2015).

[Link](https://arxiv.org/abs/1405.5327). Inspection: search snippet.

Methods for singularities outside the algebraic class. Background for interpreting the logarithmic singularity S01 finds on the negative axis of the bond series.

### M04 Asymptotically faster algorithm for counting self-avoiding walks and self-avoiding polygons

Samuel Zbarsky. 2019. J. Phys. A 52, 505001 (2019).

[Link](https://arxiv.org/abs/1903.04054). Inspection: primary abstract and HTML text via fetch summary.

Improves transfer-matrix enumeration from exp(Cn) to exp(C sqrt(n log n)) in 2D, using pigeonhole over residue classes of cut columns plus inclusion-exclusion, for objects that are sparse. The author suggests it may generalize; whether DP path-union graphs are sparse in the required sense is an open question for M3.

### M05 A parallel algorithm for the enumeration of self-avoiding polygons on the square lattice

Iwan Jensen. 2003. J. Phys. A 36 (2003); DOI 10.1088/0305-4470/36/21/304 as listed on arXiv.

[Link](https://arxiv.org/abs/cond-mat/0301468). Inspection: primary abstract.

Parallel transfer-matrix design: the state vector is partitioned so that each task updates its own piece and data is redistributed when the update position moves across the partition. Template for multithreading and out-of-core streaming in M3.

### M06 A new transfer-matrix algorithm for exact enumerations: self-avoiding polygons on the square lattice

Nathan Clisby and Iwan Jensen. 2012. J. Phys. A 45, 115202 (2012).

[Link](https://arxiv.org/abs/1111.5877). Inspection: primary abstract and HTML text sections 1 to 2.2 via fetch agent.

Square-lattice self-avoiding polygons extended from 110 to 130 steps. Key change: instead of tracking how boundary edges are connected to the left (the past), the state records how they must be connected to the right (the future). Pruning then has a unique cheapest completion, and its cost drops from exponential in the width to O(W). The update rules become more complicated. This is the model for computing an exact per-configuration contribution order in the DP code (REVIEW section 10).

### M07 D-algebraic guessing

Bertrand Teguia Tabuguia. 2025. arXiv 2510.26869.

[Link](https://arxiv.org/abs/2510.26869). Inspection: primary abstract.

Guessing polynomial differential equations from finitely many terms. Possible tool for the M6 solvability tests beyond D-finite guessing.

### M08 The art of algorithmic guessing in gfun

see paper. 2022. arXiv 2209.03059.

[Link](https://arxiv.org/abs/2209.03059). Inspection: search result title.

Practical guide to guessing linear ODEs and recurrences (D-finite) from series. Tooling reference for M6.

### M09 A new transfer-matrix algorithm for exact enumerations: self-avoiding walks on the square lattice

Iwan Jensen. 2013. arXiv 1309.6709.

[Link](https://arxiv.org/abs/1309.6709). Inspection: arXiv metadata and search snippet.

Applies the M06 approach to walks. States that the gain comes from making it simple to calculate the order at which a configuration will contribute, at the cost of more complicated update rules.

### M10 Transfer matrix computation of generalised critical polynomials in percolation

Christian R. Scullard and Jesper Lykke Jacobsen. 2012. arXiv 1209.1451.

[Link](https://arxiv.org/abs/1209.1451). Inspection: arXiv metadata and search snippet.

Critical-polynomial method for isotropic percolation thresholds, computed by transfer matrix. No application to directed percolation was found. The method relies on isotropic structure that DP lacks, so it is listed as a technique considered and set aside.

### M11 High-precision percolation thresholds and Potts-model critical manifolds from graph polynomials

Jesper Lykke Jacobsen. 2014. arXiv 1401.7847.

[Link](https://arxiv.org/abs/1401.7847). Inspection: arXiv metadata and search snippet.

Graph-polynomial thresholds to about 10 to 15 digits for 2D isotropic lattices. A benchmark for what large transfer-matrix computations achieve now, not a DP method.

### M12 High order Fuchsian equations for the square lattice Ising model: chi-tilde^(5)

A. Bostan, S. Boukraa, A. J. Guttmann, S. Hassani, I. Jensen, J.-M. Maillard and N. Zenine. 2009. J. Phys. A 42, 275209 (2009) (DOI listed on arXiv).

[Link](https://arxiv.org/abs/0904.1601). Inspection: arXiv metadata and search snippet.

Found an order-33 Fuchsian ODE from a series computed modulo a single prime, without full CRT reconstruction. Template for testing DP series for D-finiteness directly from per-prime enumeration output. The Ising work used thousands of terms, so with about 250 DP terms only small ODEs can be excluded.

### M13 High order Fuchsian equations for the square lattice Ising model: chi^(6)

S. Boukraa, S. Hassani, I. Jensen, J.-M. Maillard and N. Zenine. 2010. J. Phys. A 43, 115201 (2010) (DOI listed on arXiv).

[Link](https://arxiv.org/abs/0912.4968). Inspection: arXiv metadata and search snippet.

ODE factors guessed modulo several primes, then reconstructed in exact arithmetic. Shows the mod-prime-guess-then-lift workflow.

### M14 Differential equations for algebraic functions

A. Bostan, F. Chyzak, G. Lecerf, B. Salvy and E. Schost. 2007. ISSAC 2007.

[Link](https://search.crossref.org/search/works?q=Bostan+Chyzak+Lecerf+Salvy+Schost+Differential+equations+for+algebraic+functions+ISSAC+2007&from_ui=yes). Inspection: metadata from background knowledge; not retrieved this session (link is a Crossref search).

Source of the order-degree trade-off: a non-minimal-order ODE can have much lower degree, so fewer series terms are needed to guess it. Relevant to how many terms a D-finiteness test needs.

## Other numerical routes and exact-form claims

### N01 Precise critical exponents for the basic contact process

J. Ricardo G. de Mendonca. 1999. J. Phys. A 32, L467 (1999).

[Link](https://arxiv.org/abs/cond-mat/9906455). Inspection: primary abstract and HTML text via fetch summary.

Exact diagonalization of the master operator for L up to 18 with Bulirsch-Stoer extrapolation: lambda_c = 1.64896(2), z = 1.58077(2), nu_perp = 1.09681(2), delta = 0.162(2). Within 1 to 2 of its own error bars of the S01 values, but 2 to 5 times less precise than series.

### N02 Density matrix renormalization group and reaction-diffusion processes

E. Carlon, M. Henkel and U. Schollwoeck. 1999. Eur. Phys. J. B 12, 99 (1999).

[Link](https://search.crossref.org/search/works?q=Carlon+Henkel+Schollwock+density+matrix+renormalization+group+reaction-diffusion+processes+1999&from_ui=yes). Inspection: secondary description in the DMRG review cond-mat/0409292; venue from memory and not verified (link is a Crossref search).

DMRG for 1D reaction-diffusion including the contact process, with exponents consistent with DP. Benchmark consumer of improved values.

### N03 High-precision Monte Carlo study of directed percolation in (d+1) dimensions

Junfeng Wang, Zongzheng Zhou, Qingquan Liu, Timothy M. Garoni and Youjin Deng. 2013. Phys. Rev. E 88, 042102 (2013).

[Link](https://arxiv.org/abs/1201.3006). Inspection: primary abstract; journal reference verified on arXiv.

Covers 2 <= d <= 7 only; 1+1 is left to series. Confirms that series is the reference method in 1+1.

### N04 A transfer-matrix study of directed lattice animals and directed percolation on a square lattice

Dragica Knezevic and Milan Knezevic. 2016. J. Phys. A 49, 115001 (2016).

[Link](https://doi.org/10.1088/1751-8113/49/11/115001). Inspection: catalog abstract only.

Finite-strip transfer matrices give sequences of effective p_c and correlation-length exponent estimates. Numbers not retrieved; check precision against S01.

### N05 Transverse-size critical exponent of directed percolation from Yang-Lee zeros of survival probability

Milan Knezevic and Milos Knezevic. 2020. Phys. Rev. E 101, 012107 (2020).

[Link](https://doi.org/10.1103/PhysRevE.101.012107). Inspection: catalog abstract only.

Strip transfer matrix for survival probability plus Yang-Lee zeros gives nu_perp and p_c precisely. Value not retrieved; it is the most recent dedicated 1+1 nu_perp estimate found.

### N06 On directed interacting animals and directed percolation

Milan Knezevic and Jean Vannimenus. 2002. J. Phys. A (2002).

[Link](https://arxiv.org/abs/cond-mat/0203367). Inspection: search snippet.

Strip transfer matrix linking interacting directed animals and DP. Precursor of N04.

### N07 Directed percolation and the golden ratio

Stephan M. Dammer, Silvio R. Dahmen and Haye Hinrichsen. 2001. arXiv cond-mat/0106396.

[Link](https://arxiv.org/abs/cond-mat/0106396). Inspection: primary abstract and HTML text via fetch summary.

Exact survival probability at unphysical complex p = 2 and p = (1 +- sqrt 5)/2. No exponent conjecture. Quotes beta = 0.27649(4) and p_c = 0.6447001(1) from Jensen 1996 PRL. Shows exact structure exists away from the physical critical point.

### N08 A precise approximation for directed percolation in d=1+1

Clement Sire. 2002. Phys. Rev. E 66, 066128 (2002).

[Link](https://arxiv.org/abs/cond-mat/0205490). Inspection: primary abstract and HTML text via fetch summary.

Closure approximation giving p_c = (sqrt 5 - 1)/2 = 0.618 and beta = (1 - 1/sqrt 5)/2 = 0.276393. The beta differs from 0.276486(8) by about 12 error bars, so it is an approximation, not an exact value.

### N09 High-precision simulations of the parity conserving directed percolation universality class in 1+1 dimensions

Peter Grassberger. 2025. Phys. Rev. E (2026); arXiv 2512.02241.

[Link](https://arxiv.org/abs/2512.02241). Inspection: primary abstract and search snippets.

Different universality class, but a recent model of how near-rational exponent claims are argued (beta = 1.000 and nu = 2 there, z not rational). Useful precedent for framing the M6 exact-form tests.

### N10 Entropy governed by the absorbing state of directed percolation

Kenji Harada and Naoki Kawashima. 2019. Phys. Rev. Lett. 123, 090601 (2019).

[Link](https://arxiv.org/abs/1902.10479). Inspection: search snippet.

Tensor-network study of 1+1 DP. Its exponents are for entropy relaxation (delta_bond = 0.1593(1) etc.), not the standard set. A 2020 companion (arXiv 2008.10807) studies spectrum structure at criticality.

### N11 Directed percolation in two dimensions: An exact solution

Shiwei Chen and Fa-Yueh Wu (per fetch summary). 2006. Nankai Tracts in Mathematics vol. 10 (World Scientific, 2006), pp. 160-168.

[Link](https://arxiv.org/abs/cond-mat/0511296). Inspection: primary abstract.

Despite the title, solves a different row-alternating model whose criticality is in aspect ratio. Not a solution of square-lattice 1+1 DP; listed to prevent misreading.

## Field theory, bootstrap and scaling theory

### T01 On the nonequilibrium phase transition in reaction-diffusion systems with an absorbing stationary state

H. K. Janssen. 1981. Z. Phys. B 42, 151 (1981).

[Link](https://search.crossref.org/search/works?q=Janssen+nonequilibrium+phase+transition+reaction-diffusion+absorbing+stationary+state+1981&from_ui=yes). Inspection: metadata from background knowledge; DOI not verified (link is a Crossref search).

Half of the DP conjecture: generic absorbing-state transitions with a single absorbing state are DP. Background.

### T02 On phase transitions in Schloegl's second model

P. Grassberger. 1982. Z. Phys. B 47, 365 (1982).

[Link](https://search.crossref.org/search/works?q=Grassberger+phase+transitions+Schlogl+second+model+1982&from_ui=yes). Inspection: metadata from background knowledge; DOI not verified (link is a Crossref search).

Other half of the Janssen-Grassberger DP conjecture. Background.

### T03 Field-theoretic analysis of directed percolation: Three-loop approximation

L. Ts. Adzhemyan, M. Hnatic, E. V. Ivanova, M. V. Kompaniets, T. Lucivjansky and L. Mizisin. 2023. Phys. Rev. E 107, 064138 (2023).

[Link](https://arxiv.org/abs/2306.17057). Inspection: primary abstract and HTML text via fetch summary; author list verified on arXiv.

All DP exponents to O(epsilon^3), epsilon = 4 - d. Not evaluated at d = 1 (epsilon = 3), where the truncated series is useless without resummation. A consumer of better 1+1 numbers as a resummation test.

### T04 Renormalization group analysis of directed percolation process: Towards multiloop calculation of scaling functions

M. Hnatic, M. Kecer, T. Lucivjansky and L. Mizisin. 2026. Theor. Math. Phys. 227, 984 (2026).

[Link](https://arxiv.org/abs/2602.11369). Inspection: arXiv abstract page metadata.

Extends the three-loop DP analysis to the equation of state. Indicates active field-theory interest in 2026.

### T05 Nonperturbative renormalization-group study of reaction-diffusion processes

L. Canet, B. Delamotte, O. Deloubriere and N. Wschebor. 2004. Phys. Rev. Lett. 92, 195703 (2004).

[Link](https://arxiv.org/abs/cond-mat/0309504). Inspection: search snippet.

NPRG recovers DP exponents in low dimensions. Consumer of reference values.

### T06 The nonperturbative functional renormalization group and its applications

N. Dupuis, L. Canet, A. Eichhorn, W. Metzner, J. M. Pawlowski, M. Tissier and N. Wschebor. 2021. Phys. Rep. 910, 1 (2021).

[Link](https://arxiv.org/abs/2006.04853). Inspection: search snippet; author list from memory, not verified.

Review stating that DP exponents in d = 1, 2, 3 have been computed at LPA and LPA' and compare accurately with other methods. Points to the FRG d = 1 numbers to be tabulated in REVIEW.

### T07 Background field functional renormalization group for absorbing state phase transitions

Michael Buchhold and Sebastian Diehl. 2016. Phys. Rev. E 94, 012138 (2016).

[Link](https://arxiv.org/abs/1604.07612). Inspection: search snippet.

FRG variant reporting good DP exponents even in low dimensions. Consumer of reference values.

### T08 Bootstrapping nonequilibrium stochastic processes

Minjae Cho. 2025. SciPost Phys. 19, 124 (2025).

[Link](https://arxiv.org/abs/2505.13609). Inspection: primary HTML text, first part verbatim and remainder via fetch summary.

Positivity bootstrap with LP and SDP gives rigorous bounds for the contact process and Domany-Kinzel model. Current bounds are loose (contact-process lambda_c >= 1.46 at L = 10 against 1.6489). Frames 1+1 DP as unsolved. A potential future consumer of reference values.

### T09 Fermionic field theory for directed percolation in (1+1) dimensions

V. Brunel, K. Oerding and F. van Wijland. 1999. arXiv cond-mat/9911095.

[Link](https://arxiv.org/abs/cond-mat/9911095). Inspection: search snippet; authors verified on arXiv.

An approximation scheme specific to d = 1. Example of a 1D-specific analytic approach that needs sharp benchmarks.

### T10 On logarithmic extensions of local scale-invariance

Malte Henkel. 2013. Nucl. Phys. B 869, 282 (2013); arXiv 1009.4139.

[Link](https://arxiv.org/abs/1009.4139). Inspection: search snippet; metadata verified on arXiv.

Local scale invariance tested against transfer-matrix renormalization group data for the 1D critical contact process; the non-logarithmic form shows deviations (cond-mat/0603369). The tests are limited by how precisely DP exponents are known.

### T11 Convergence of non-perturbative approximations to the renormalization group

Ivan Balog, Hugues Chate, Bertrand Delamotte, Maroje Marohnic and Nicolas Wschebor. 2019. Phys. Rev. Lett. 123, 240604 (2019).

[Link](https://arxiv.org/abs/1907.01829). Inspection: arXiv metadata.

The derivative expansion of the functional RG, carried to high order, reaches conformal-bootstrap-level precision for the 3D Ising model. Makes precise 1+1 DP values more valuable as a target for the same programme applied to DP.

## Experiments and physical realizations

### E01 Directed percolation phase transition to sustained turbulence in Couette flow

Gregoire Lemoult, Liang Shi, Kerstin Avila, Shreyas V. Jalikop, Marc Avila and Bjoern Hof. 2016. Nature Physics 12, 254 (2016).

[Link](https://www.nature.com/articles/nphys3675). Inspection: search snippets; author list from memory, verify.

Quasi-1D circular Couette experiment, the main 1+1 DP realization in turbulence. Fitted exponent values not retrieved.

### E02 Scale invariance at the onset of turbulence in Couette flow

Liang Shi, Marc Avila and Bjoern Hof. 2013. Phys. Rev. Lett. 110, 204502 (2013).

[Link](https://doi.org/10.1103/PhysRevLett.110.204502). Inspection: secondary quote in search snippet.

beta = 0.29 +- 0.05 against the 1+1 DP value 0.276. Shows experimental precision is about two digits.

### E03 Directed percolation and puff jamming near the transition to pipe turbulence

Gregoire Lemoult, Vasudevan Mukund, Hong-Yan Shih, Gaute Linga, Joachim Mathiesen, Nigel Goldenfeld and Bjoern Hof. 2024. Nature Physics 20, 1339 (2024).

[Link](https://doi.org/10.1038/s41567-024-02513-0). Inspection: search snippets and publisher listing.

Pipe flow transition is in the DP class (effectively 1+1), with a jammed puff phase above criticality. Exponent values not retrieved.

### E04 A universal transition to turbulence in channel flow

Masaki Sano and Keiichi Tamai. 2016. Nature Physics 12, 249 (2016).

[Link](https://www.nature.com/articles/nphys3659). Inspection: search snippet.

Channel flow, 2+1 DP. Not a 1+1 test; listed to separate dimensionalities.

### E05 Directed percolation criticality in turbulent liquid crystals

Kazumasa A. Takeuchi, Masafumi Kuroda, Hugues Chate and Masaki Sano. 2007. Phys. Rev. Lett. 99, 234503 (2007).

[Link](https://arxiv.org/abs/0706.4151). Inspection: primary abstract.

The liquid-crystal experiment cited in the idea note. It is (2+1)-dimensional DP, not 1+1.

### E06 Experimental realization of directed percolation criticality in turbulent liquid crystals

Kazumasa A. Takeuchi, Masafumi Kuroda, Hugues Chate and Masaki Sano. 2009. Phys. Rev. E 80, 051116 (2009).

[Link](https://arxiv.org/abs/0907.4297). Inspection: primary abstract.

Twelve exponents, five scaling functions and eight scaling relations, all 2+1 DP. Its comparison table lists earlier 1+1 experiments.

### E07 Directed percolation describes lifetime and growth of turbulent puffs and slugs

Maksim Sipos and Nigel Goldenfeld. 2011. Phys. Rev. E 84, 035304(R) (2011).

[Link](https://arxiv.org/abs/1101.5835). Inspection: arXiv metadata verified.

DP model of pipe-flow puff lifetimes and splitting, consistent with super-exponential lifetime data.

### E08 Characterizing a non-equilibrium phase transition on a quantum computer

Eli Chertkov and 15 coauthors. 2022. Nature Physics 19, 1799 (2023).

[Link](https://arxiv.org/abs/2209.12889). Inspection: search snippet.

Floquet quantum contact process on 73 trapped-ion qubits; data and simulations consistent with classical 1+1 DP. A recent consumer of reference exponents.

## Reviews and books

### R01 Non-equilibrium critical phenomena and phase transitions into absorbing states

Haye Hinrichsen. 2000. Adv. Phys. 49, 815 (2000).

[Link](https://arxiv.org/abs/cond-mat/0001070). Inspection: primary abstract.

Standard review. States that DP is one of very few critical phenomena whose 1+1 exponents are not known exactly, and that they do not seem rational.

### R02 On possible experimental realizations of directed percolation

Haye Hinrichsen. 2000. Braz. J. Phys. 30, 69 (2000).

[Link](http://www.scielo.br/j/bjp/a/XdPLNw7wCkz4PwhM7MDdnzr/?lang=en). Inspection: search snippet.

Why DP was hard to observe experimentally before 2007. Context for the experiments section.

### R03 Universality classes in nonequilibrium lattice systems

Geza Odor. 2004. Rev. Mod. Phys. 76, 663 (2004).

[Link](https://arxiv.org/abs/cond-mat/0205644). Inspection: metadata.

Broad review of absorbing-state universality classes with exponent tables.

### R04 Universal scaling behavior of non-equilibrium phase transitions

Sven Luebeck. 2004. Int. J. Mod. Phys. B 18, 3977 (2004).

[Link](https://arxiv.org/abs/cond-mat/0501259). Inspection: search snippet.

Review with universal scaling functions and amplitude ratios for DP. Notes that DP exponents do not seem rational, unlike isotropic percolation.

### R05 Non-Equilibrium Phase Transitions, Volume 1: Absorbing Phase Transitions

Malte Henkel, Haye Hinrichsen and Sven Luebeck. 2008. Springer (2008).

[Link](https://doi.org/10.1007/978-1-4020-8765-3). Inspection: metadata from background knowledge; not retrieved this session.

Textbook reference for DP scaling theory and exponent tables.

### R06 The field theory approach to percolation processes

Hans-Karl Janssen and Uwe C. Taeuber. 2005. Ann. Phys. 315, 147 (2005).

[Link](https://arxiv.org/abs/cond-mat/0409670). Inspection: arXiv metadata verified.

Field-theory review for DP and dynamic isotropic percolation.
