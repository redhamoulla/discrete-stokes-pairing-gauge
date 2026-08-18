# Discrete Stokes Pairing Gauge

This repository accompanies the paper:

> Redha Moulla, *Gauge Freedom in Discrete Stokes Structures: Affine Pairings, Intrinsic Nodes, and the Rigidity of Gauss-Legendre Locality*, 2026.

The paper studies finite-dimensional Stokes structures before a quadrature rule, metric, or nodal grid is chosen. Its main results include:

- the affine classification of mixed pairings compatible with a discrete Stokes identity;
- an exact full-rank existence criterion and a classification of admissible radical lines;
- the image, fibres, and topology of the polynomial kernel map;
- the hyperbolic balance condition for admissible interior real nodes;
- rigidity of Gauss-Legendre nodes and weights under scalar locality;
- a canonical minimum-gauge representative and its rank stratification;
- the global synthesis obstruction, the exceptional orders `N = 2, 4, 8`, and an explicit two-chart construction;
- unavoidable radial blow-up near the deleted constant class;
- connections with rectangular summation-by-parts identities and continuous port-Hamiltonian gauges.

## Repository layout

```text
paper/          Final manuscript source and compiled PDF
experiments/    Reproducible finite-dimensional Stokes checks
results/        Reference JSON output from the experiments
```

## Manuscript

- [PDF](paper/moulla_2026_stokes_pairing_gauge.pdf)
- [LaTeX source](paper/moulla_2026_stokes_pairing_gauge.tex)

Compile with:

```bash
make paper
```

This requires a standard LaTeX installation with `latexmk` and `pdflatex`.

## Reproducible checks

The numerical script constructs discrete Stokes data for transformed polynomial, periodic trigonometric, exponential, Gaussian radial-basis, and neural `tanh` trial spaces. It reports:

- the rank of the mixed pairing matrix;
- the dimension lost under differentiation;
- the Stokes identity residual;
- the minimal endpoint-trace completion;
- the residual of a balanced rank factorization.

These checks illustrate the abstract construction; the mathematical results of the paper are proved analytically and do not depend on the experiments.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
make experiments
```

The deterministic reference output is stored in [`results/automatic_stokes_examples.json`](results/automatic_stokes_examples.json).

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## Author

Redha Moulla  
AXIA, France  
<redha.moulla@axia-conseil.com>

