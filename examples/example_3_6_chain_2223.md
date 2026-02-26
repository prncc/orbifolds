# Example 3.6: chain(2,2,2,3) Non-Calabi-Yau Chain

$W_A = x_1^2 x_2 + x_2^2 x_3 + x_3^2 x_4 + x_4^3$, $A = \text{chain}(2,2,2,3)$
  $n = 4$, $d = 3$, $\det(A) = 24$, weights $Q = (1, 1, 1, 1)$
  CY condition: $J A^{-T} J^T = 4/3$ (NOT integer, CY FAILS)
  $|G| = 3$

## Orbifold Cohomology Basis: 8 elements
  Integer-age elements: 6
  Fractional-age elements: 2

| Element | $\lambda$ | $\gamma$ | age($\lambda$) | age($\gamma$) | Hodge | Integer ages? |
|---------|--------|-------|----------|----------|-------|---------------|
| $x_{1}^{1} x_{2}^{1} x_{3}^{1} x_{4}^{3}  e_{1} e_{2} e_{3} e_{4}$ | (0, 0, 0, 0) | (1, 1, 1, 3) | 0 | 2 | (1,1) | yes |
| $x_{1}^{1} x_{2}^{1} x_{3}^{2} x_{4}^{2}  e_{1} e_{2} e_{3} e_{4}$ | (0, 0, 0, 0) | (1, 1, 2, 2) | 0 | 2 | (1,1) | yes |
| $x_{1}^{1} x_{2}^{2} x_{3}^{1} x_{4}^{2}  e_{1} e_{2} e_{3} e_{4}$ | (0, 0, 0, 0) | (1, 2, 1, 2) | 0 | 2 | (1,1) | yes |
| $x_{1}^{1} x_{2}^{2} x_{3}^{2} x_{4}^{1}  e_{1} e_{2} e_{3} e_{4}$ | (0, 0, 0, 0) | (1, 2, 2, 1) | 0 | 2 | (1,1) | yes |
| $x_{1}^{2} x_{2}^{1} x_{3}^{1} x_{4}^{2}  e_{1} e_{2} e_{3} e_{4}$ | (0, 0, 0, 0) | (2, 1, 1, 2) | 0 | 2 | (1,1) | yes |
| $x_{1}^{2} x_{2}^{1} x_{3}^{2} x_{4}^{1}  e_{1} e_{2} e_{3} e_{4}$ | (0, 0, 0, 0) | (2, 1, 2, 1) | 0 | 2 | (1,1) | yes |
| $y_{1}^{1} y_{2}^{1} y_{3}^{1} y_{4}^{1}$ | (1, 1, 1, 1) | (0, 0, 0, 0) | 4/3 | 0 | (1/3,1/3) | NO |
| $y_{1}^{2} y_{2}^{2} y_{3}^{2} y_{4}^{2}$ | (2, 2, 2, 2) | (0, 0, 0, 0) | 8/3 | 0 | (5/3,5/3) | NO |

## Computations at $p = 73$ ($\det(A) = 24 \mid 72 = p-1$)

### Integer-age element eigenvalues

| Element | Hodge | Eigenvalue |
|---------|-------|------------|
| $x_{1}^{1} x_{2}^{1} x_{3}^{1} x_{4}^{3}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | 73 |
| $x_{1}^{1} x_{2}^{1} x_{3}^{2} x_{4}^{2}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | 73 |
| $x_{1}^{1} x_{2}^{2} x_{3}^{1} x_{4}^{2}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | 73 |
| $x_{1}^{1} x_{2}^{2} x_{3}^{2} x_{4}^{1}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | 73 |
| $x_{1}^{2} x_{2}^{1} x_{3}^{1} x_{4}^{2}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | 73 |
| $x_{1}^{2} x_{2}^{1} x_{3}^{2} x_{4}^{1}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | 73 |

### Partial count (integer-age elements only): $N_{\text{partial}} = 438$
  $H^{1,1}$ has $6$ integer-age elements, each with eigenvalue $p$
  Contribution from $H^{1,1}$: $6 \cdot 73 = 438$

### True point count (brute force): $N_{\text{true}} = 5841$

### Comparison
  $N_{\text{true}} - N_{\text{partial}} = 5841 - 438 = 5403$
  $1 + p + p^2 = 1 + 73 + 73^2 = 5403$

PASS: $N_{\text{true}} - N_{\text{partial}} = 1 + p + p^2 = 5403$
  (The fractional-age sectors would contribute exactly $|\mathbb{P}^2(\mathbb{F}_p)|$ if CY held)
