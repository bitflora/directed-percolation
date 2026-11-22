# Directed percolation in 1+1 dimensions

The goal is to extend the low-density series for bond and site directed percolation on the directed square lattice past Jensen's 1999 record of 171 and 158 terms, then sharpen the critical point and exponents and test proposed exact forms.

Literature search completed November 21, 2026. This folder contains the findings and a proposed study. No code has been written and nothing has been computed.

Main findings:

- The 1999 record still stands.
- The published rational conjecture for the bulk exponents is already excluded by 10 to 90 error bars, so the "hidden exact solution" angle is a long shot.
- A sharper near-exact question remains at a wall: τ₁ = 1.00014(2) against the conjecture τ₁ = 1.
- The 1+1 experiments are Couette and pipe flow. The liquid-crystal experiments are 2+1. All measure exponents to about two digits.
- On this 16 GB machine, the 1999 algorithm reaches only about order 179. A worthwhile extension needs algorithmic gains, out-of-core streaming, or rented memory.

Read these files in order:

1. [Literature review](REVIEW.md): the algorithm, current best values, status of exact-form proposals, competing methods, experiments, novelty, and (section 10) techniques developed since 1999 that could help.
2. [Research plan](RESEARCH_PLAN.md): milestones M0–M7, with a go/no-go gate at M2 and a compute decision at M4.
3. [Annotated bibliography](SOURCES.md): 68 records, each with its inspection level.
4. [Search log](SEARCH_LOG.md): queries, citation sweep, access limits, and remaining checks.
5. [Machine-readable sources](sources.json): the same bibliography.

Recommended first step: M0 and M1, an independent reimplementation of the 1999 non-nodal transfer matrix in Rust, checked against brute force at low order. Then M2: analyze the series truncated at 112, 130, 150 and 171 terms to measure how much precision more terms would actually buy. Only then spend effort on extension.

Key prior art: [Jensen 1999, square lattice to 171 terms](https://arxiv.org/abs/cond-mat/9906036), [Jensen 1999, square lattice with a wall](https://arxiv.org/abs/cond-mat/9906037), [Jensen 2004, other lattices](https://arxiv.org/abs/cond-mat/0405504), [Essam, Guttmann and De'Bell 1988, rational conjecture](https://iopscience.iop.org/article/10.1088/0305-4470/21/19/018), [Hinrichsen review](https://arxiv.org/abs/cond-mat/0001070).
