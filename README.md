# Orbifold Zeta Functions

Code associated to [arXiv:2602.23173](https://arxiv.org/abs/2602.23173).

## Setup

Requires [SageMath](https://www.sagemath.org/). Run scripts with:

```
sage path/to/script.sage
```

## Modules

Code is loaded using `load('module.sage')` statements for simplicity.

### `potentials.sage`

Matrix potential decomposition and Koszul cohomology basis computation.
Decomposes invertible integer matrices into Fermat, chain, and loop
basic potentials, and computes cohomology bases for each type before
combining them again.

### `orbifold.sage`

Orbifold cohomology and trace formula. Computes symmetry groups,
lambda-sector bases, and Frobenius trace contributions for weighted
Fermat hypersurfaces.

### `orbifold_display.sage`

Display and visualization utilities. LaTeX rendering of trace
contributions, Hodge tables, zeta function assembly, etc.

### `diagonal_k3.sage`

Computes point counts and zeta functions for weighted diagonal K3
surfaces via Jacobi sums and cyclic quotient singularity resolution,
following Goto, "Arithmetic of weighted diagonal surfaces over finite
fields," *J. Number Theory* **59**(1), 1996, Theorem 5.2.

### `deformed_diagonal_k3.sage`

Computes point counts and zeta functions for weighted deformed diagonal
K3 surfaces, following Goto, "Arithmetic of weighted diagonal surfaces
and weighted deformed diagonal surfaces over finite fields," Ph.D.
thesis, Queen's University, 1994, Theorem 4.3.4 / Corollary 4.3.7.

### `quintic.sage`

Point counting on the Dwork quintic threefold over finite fields via
Gauss sum formulas from Candelas, de la Ossa, and Rodriguez Villegas,
"Calabi-Yau manifolds over finite fields, I," arXiv:hep-th/0012233,
2000, Section 9.

### `quintic_factorization.sage`

Zeta function factorization for the Dwork quintic at psi = 0, rho = 1.
Replicates Table 12.1 from Candelas, de la Ossa, and Rodriguez Villegas,
"Calabi-Yau manifolds over finite fields. II," *Fields Inst. Commun.*
**38**, 2003, pp. 121-157.

### `orbifold_explorer.ipynb`

Interactive Jupyter notebook for exploring orbifold cohomology and zeta
functions. Configure a matrix A and symmetry group, and the notebook
finds the first eligible prime, displays all cohomology elements in a
table, computes point counts over $F_p$ and $F_{p^2}$, renders the
orbifold zeta function, and expands the log zeta to extract point counts.
Requires a SageMath Jupyter kernel.

## Tests

Tests live in `tests/` and are run via `bash tests/run_tests.sh`.

- `potentials_test.sage` — unit tests for matrix decomposition and basis computation.
- `orbifold_test.sage` — unit tests for lambda-sector basis computation.
- `diagonal_k3_test.sage` — integration tests for select weighted diagonal K3 surfaces.
- `deformed_diagonal_k3_test.sage` — integration tests for select deformed diagonal K3 surfaces.

## Examples

The `examples/` directory contains worked examples that use the core
modules for various purposes. Each `.sage` script
produces a companion `.md` file with formatted output. Examples cover
Fermat, chain, and loop potentials, diagonal and deformed K3 surfaces,
the Fermat quintic, and Greene-Plesser mirror pairs.

## Calculations

The `calculations/` directory contains scripts that generate tables
and verification data:

- `diagonal_k3_zetas.sage` — zeta functions for diagonal K3 surfaces.
- `deformed_diagonal_k3_zetas.sage` — zeta functions for deformed diagonal K3 surfaces.
- `quintic_point_counts.sage` — point counts for the Dwork quintic.
- `quintic_zeta.sage` — zeta function verification for the quintic.
- `per_element_contributions/` — per-element trace contribution TeX tables for diagonal and deformed diagonal K3 surfaces.
