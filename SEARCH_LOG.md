# Literature search log

Search date: November 21, 2026.

## Scope and method

The questions were:

- whether Jensen's 1999 series record for 1+1 directed percolation has been broken;
- what the current best values are and how they were obtained;
- which exact or closed-form proposals exist and whether they survive;
- which other methods compete or would use better numbers;
- what the algorithm needs in order to be reproduced.

Methods used:

- **Web searches**, about 30 queries in standard and extended modes, listed below.
- **A citation sweep.** The Semantic Scholar API listed every paper citing cond-mat/9906036, which returned a single page with no continuation. Titles from 2015 to 2026 were screened for series extensions or new 1+1 precision estimates. None were found.
- **Full-text reading.** Jensen 1999 (S01) was read in full from the arXiv PDF, extracted with pdftotext. The figure symbols in the update rules did not survive extraction.
- **Abstract or HTML reading** through a fetch agent for S02, S03, M04, N01, N07, N08, N11, T03 and T08.
- **Identifier verification.** Crossref and arXiv were checked for 20 identifiers. Two guessed identifiers resolved to unrelated papers and were replaced by Crossref search links. Five Crossref lookups hit rate limits and remain unverified; they are labelled as such in SOURCES.md.

This is a broad map, not a systematic review with a complete denominator. That no extension was found is evidence, not proof, that none exists.

## Query batches

### Batch 1: record and core values

- Jensen low-density series expansions directed percolation square lattice order 171 critical exponents
- directed percolation 1+1 dimensions series expansion new record longer series 2020..2026 (extended)
- high-precision Monte Carlo directed percolation (1+1) dimensions critical exponents Wang Zhou Garoni Deng
- Grassberger directed percolation 1+1 dimensions high precision simulation 2024 2025 critical point exponents (extended)
- contact process series expansion Dickman Jensen time-dependent perturbation theory critical point

### Batch 2: exact-form proposals and solvability

- directed percolation exponents conjectured rational values … exact conjecture ruled out (extended)
- "directed percolation" exponents rational conjecture Essam Guttmann De'Bell 1988 … 26/15 (extended)
- directed percolation 1+1 critical exponents "not rational" OR "exact value" conjecture 2010..2026 (extended)
- directed percolation series not D-finite unsolvability Guttmann Bousquet-Mélou directed animals exact solution (extended)
- directed percolation generating function D-finite test differential equation guessing … haruspicy
- Inui Katori Catalan numbers series expansions directed percolation exact coefficients
- directed percolation local scale invariance conformal invariance 1+1 exact exponent conjecture Henkel

### Batch 3: competing methods and theory

- directed percolation four-loop renormalization group exponents 2023 Adzhemyan Kompaniets
- non-perturbative renormalization group directed percolation 1+1 exponents Canet Delamotte derivative expansion 2024
- functional renormalization group directed percolation d=1 exponents … 2020..2026
- directed percolation 1+1 tensor network DMRG critical exponents contact process precise estimate
- Knežević 2016 "transfer-matrix study of directed lattice animals and directed percolation on a square lattice"
- Knežević 2020 Yang-Lee zeros survival probability directed percolation transverse exponent
- "Bootstrapping nonequilibrium stochastic processes" directed percolation bounds
- directed percolation 1+1 Monte Carlo GPU large scale simulation critical exponents 2021–2023 (extended)

### Batch 4: experiments

- directed percolation turbulence experiment Couette flow 1+1 Lemoult Nature Physics exponents
- liquid crystal electroconvection turbulence directed percolation quasi-one-dimensional 1+1 Takeuchi exponents
- pipe flow transition turbulence directed percolation 1+1 universality Hof Barkley experiment exponents 2023 2024
- Lemoult 2024 Nature Physics directed percolation puff jamming pipe exponents
- "directed percolation" quantum computer OR Rydberg experiment 1+1 critical exponents 2022

### Batch 5: algorithms, data availability, analysis methods

- Jensen "directed percolation" series data website Melbourne series coefficients download
- Jensen transfer matrix parallel algorithm enumeration memory pruning self-avoiding polygons 2003 2004
- differential approximants modern series analysis methods Guttmann review 2015 2020
- Hinrichsen 2000 review directed percolation table exponents d=1

## Access limitations

- **Jensen's series files.** S01 and S04 say they were available by e-mail or from a Melbourne web page. No current copy was found online. Asking the author is possible, but it is an outward-facing step for the user to decide.
- **Paywalled papers.** Several journal versions are paywalled. Where an arXiv version exists, it is linked.
- **Formatting.** The pdftotext extraction of S01 lost the pictorial symbols in figure 3 and equation 12, so the update rules must be re-derived (RESEARCH_PLAN M0).

## Remaining checks before writing a paper

1. Read N04 and N05 (Knežević, 2016 and 2020) in full for their p_c and ν⊥ precision. These are the most recent dedicated 1+1 estimates.
2. Record the FRG d = 1 exponent values (T05–T07) for a comparison table.
3. Extract the fitted 1+1 exponents from E01 and E03.
4. Verify the five DOIs left unverified (S05, S10, S13, T01, T02) and the S03 table values, which came from a summary.
5. Confirm via Google Scholar, ADS and arXiv full-text search that no square-lattice DP series beyond 171 terms appeared in a venue Semantic Scholar does not index, such as theses or conference proceedings.
6. Check the 2016 "Semi-directed percolation in two dimensions" and the 2022 "Direct measurement of correlation length in one-dimensional contact process", both from the citation sweep, for any new 1+1 estimates.
