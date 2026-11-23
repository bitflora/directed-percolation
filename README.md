# Directed percolation in 1+1 dimensions

The goal is to extend the low-density series for bond and site directed percolation on the directed square lattice past Jensen's 1999 record of 171 and 158 terms, then sharpen the critical point and exponents and test proposed exact forms.

Literature search completed November 21, 2026. Implementation started the same night:
- an exact ground-truth checker (M0);
- a reimplemented and improved transfer-matrix engine in [enum/](enum/README.md) (M1);
- differential-approximant analysis and the precision-gain gate (M2).

Results so far are in [FINDINGS.md](FINDINGS.md), and the milestone status is at the top of [RESEARCH_PLAN.md](RESEARCH_PLAN.md).

Headline: an independent implementation reproduced Jensen's 1999 results at his own lengths. All four moment series are verified to 171 terms (bond) and 158 terms (site), and p_c, γ, ν∥, 2ν⊥ and the moment exponents all agree within his error bars, as does the site critical point from universality. The precision-gain gate passed, so extending to 200–230 terms is worth doing. Next to a wall, a new direct series for the mean cluster length (110 terms) gives τ₁ = 1.00010(2), against the conjecture τ₁ = 1.

Main findings from the literature search (computed results are in [FINDINGS.md](FINDINGS.md), summarized at its top):

- The 1999 record still stands.
- The published rational conjecture for the bulk exponents is already excluded by 10 to 90 error bars, so the "hidden exact solution" angle is a long shot.
- A sharper near-exact question remains at a wall: τ₁ = 1.00014(2) against the conjecture τ₁ = 1. A direct 110-term series computed here gives τ₁ = 1.00010(2), which disfavours τ₁ = 1 (FINDINGS F14).
- The 1+1 experiments are Couette and pipe flow. The liquid-crystal experiments are 2+1. All measure exponents to about two digits.
- With the 1999 algorithm, this 16 GB machine would reach only about order 179. The new engine needs roughly 6× less memory for the same series, which put order 171 within reach locally: about 12 min per prime for S, μ₀₁, μ₀₂ (u32, `--no-x`, 2.2 GB) and about 23 min per prime for all moments (u16, 2.1 GB working set, while another job used every core). A worthwhile extension beyond 171 still needs rented memory or a better algorithm.

Read these files in order:

0. [Findings](FINDINGS.md): what has been computed and verified, with evidence.
1. [Literature review](REVIEW.md): the algorithm, current best values, status of exact-form proposals, competing methods, experiments, novelty, and (section 10) techniques developed since 1999 that could help.
2. [Research plan](RESEARCH_PLAN.md): milestones M0–M7, with a go/no-go gate at M2 and a compute decision at M4.
3. [Annotated bibliography](SOURCES.md): 68 records, each with its inspection level.
4. [Search log](SEARCH_LOG.md): queries, citation sweep, access limits, and remaining checks.
5. [Machine-readable sources](sources.json): the same bibliography.

Code and data:

| Path | Contents |
| --- | --- |
| [enum/](enum/README.md) | `dpenum`, the Rust series engine |
| [analysis/m0_ground_truth.py](analysis/m0_ground_truth.py) | Exact ground truth by three independent methods |
| [analysis/crt.py](analysis/crt.py) | Combine per-prime engine output into verified exact series |
| [analysis/da.py](analysis/da.py), [gate.py](analysis/gate.py), [gate_report.py](analysis/gate_report.py) | Differential approximants, the analysis ensemble, the gate report |
| [analysis/dfinite.py](analysis/dfinite.py) | Search for linear ODEs modulo primes (D-finiteness test) |
| [results/series/](results/series/) | Verified exact series (JSON) |
| [results/analysis/](results/analysis/) | Analysis outputs |

The analysis tools need the venv at `analysis/.venv` (numpy, mpmath, sympy, python-flint). Create it with `python -m venv analysis/.venv` and `analysis/.venv/Scripts/python -m pip install numpy mpmath sympy python-flint`.

Key prior art: [Jensen 1999, square lattice to 171 terms](https://arxiv.org/abs/cond-mat/9906036), [Jensen 1999, square lattice with a wall](https://arxiv.org/abs/cond-mat/9906037), [Jensen 2004, other lattices](https://arxiv.org/abs/cond-mat/0405504), [Essam, Guttmann and De'Bell 1988, rational conjecture](https://iopscience.iop.org/article/10.1088/0305-4470/21/19/018), [Hinrichsen review](https://arxiv.org/abs/cond-mat/0001070).
