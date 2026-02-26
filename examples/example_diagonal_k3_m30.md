# Example: Diagonal K3 Surface (m=30, Q=(15,10,3,2))

Surface: $x_0^2 + x_1^3 + x_2^{10} + x_3^{15} = 0$ in $\mathbb{P}^3(15,10,3,2)$
  $A = \text{diag}(2, 3, 10, 15)$, $\det(A) = 900$, $d = 30$
  $|G| = 30$

## Orbifold Cohomology: 24 elements

  $h^{0,0} = 1$
  $h^{0,2} = 1$
  $h^{1,1} = 20$
  $h^{2,0} = 1$
  $h^{2,2} = 1$

## Computations at $p = 1801$

### Eigenvalue table

| Element | Hodge | $\gamma A^{-1}$ | Eigenvalue formula | tc value |
|---------|-------|--------------|-------------------|----------|
| $x_{1}^{1} x_{2}^{1} x_{3}^{1} x_{4}^{1}  e_{1} e_{2} e_{3} e_{4}$ | (0,2) | (1/2, 1/3, 1/10, 1/15) | $- \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{1}{3} \right) \Gamma_p \left( \frac{1}{10} \right) \Gamma_p \left( \frac{1}{15} \right)$ | $72 + 845 \cdot 1801 + 1550 \cdot 1801^2 + 721 \cdot 1801^3 + \cdots$ |
| $x_{1}^{1} x_{2}^{1} x_{3}^{3} x_{4}^{13}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | (1/2, 1/3, 3/10, 13/15) | $p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{1}{3} \right) \Gamma_p \left( \frac{3}{10} \right) \Gamma_p \left( \frac{13}{15} \right)$ | $1416 \cdot 1801 + 391 \cdot 1801^2 + 37 \cdot 1801^3 + 915 \cdot 1801^4 + \cdots$ |
| $x_{1}^{1} x_{2}^{1} x_{3}^{5} x_{4}^{10}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | (1/2, 1/3, 1/2, 2/3) | $p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{1}{3} \right) \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{2}{3} \right)$ | $p$ |
| $x_{1}^{1} x_{2}^{1} x_{3}^{7} x_{4}^{7}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | (1/2, 1/3, 7/10, 7/15) | $p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{1}{3} \right) \Gamma_p \left( \frac{7}{10} \right) \Gamma_p \left( \frac{7}{15} \right)$ | $1081 \cdot 1801 + 438 \cdot 1801^2 + 783 \cdot 1801^3 + 164 \cdot 1801^4 + \cdots$ |
| $x_{1}^{1} x_{2}^{1} x_{3}^{9} x_{4}^{4}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | (1/2, 1/3, 9/10, 4/15) | $p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{1}{3} \right) \Gamma_p \left( \frac{9}{10} \right) \Gamma_p \left( \frac{4}{15} \right)$ | $409 \cdot 1801 + 1590 \cdot 1801^2 + 1437 \cdot 1801^3 + 612 \cdot 1801^4 + \cdots$ |
| $x_{1}^{1} x_{2}^{2} x_{3}^{1} x_{4}^{11}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | (1/2, 2/3, 1/10, 11/15) | $p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{2}{3} \right) \Gamma_p \left( \frac{1}{10} \right) \Gamma_p \left( \frac{11}{15} \right)$ | $1026 \cdot 1801 + 1383 \cdot 1801^2 + 580 \cdot 1801^3 + 1654 \cdot 1801^4 + \cdots$ |
| $x_{1}^{1} x_{2}^{2} x_{3}^{3} x_{4}^{8}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | (1/2, 2/3, 3/10, 8/15) | $p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{2}{3} \right) \Gamma_p \left( \frac{3}{10} \right) \Gamma_p \left( \frac{8}{15} \right)$ | $903 \cdot 1801 + 410 \cdot 1801^2 + 194 \cdot 1801^3 + 913 \cdot 1801^4 + \cdots$ |
| $x_{1}^{1} x_{2}^{2} x_{3}^{5} x_{4}^{5}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | (1/2, 2/3, 1/2, 1/3) | $p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{2}{3} \right) \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{1}{3} \right)$ | $p$ |
| $x_{1}^{1} x_{2}^{2} x_{3}^{7} x_{4}^{2}  e_{1} e_{2} e_{3} e_{4}$ | (1,1) | (1/2, 2/3, 7/10, 2/15) | $p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{2}{3} \right) \Gamma_p \left( \frac{7}{10} \right) \Gamma_p \left( \frac{2}{15} \right)$ | $1525 \cdot 1801 + 1463 \cdot 1801^2 + 303 \cdot 1801^3 + 112 \cdot 1801^4 + \cdots$ |
| $x_{1}^{1} x_{2}^{2} x_{3}^{9} x_{4}^{14}  e_{1} e_{2} e_{3} e_{4}$ | (2,0) | (1/2, 2/3, 9/10, 14/15) | $- p^{2} \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{2}{3} \right) \Gamma_p \left( \frac{9}{10} \right) \Gamma_p \left( \frac{14}{15} \right)$ | $1776 \cdot 1801^2 + 1343 \cdot 1801^3 + 473 \cdot 1801^4 + 1577 \cdot 1801^5 + \cdots$ |
| $y_{1}^{1} y_{2}^{1} y_{3}^{1} y_{4}^{1}$ | (0,0) | (0, 0, 0, 0) | $1$ | $1$ |
| $x_{1}^{1} x_{3}^{5} y_{2}^{1} y_{4}^{10} e_{1} e_{3}$ | (1,1) | (1/2, 0, 1/2, 0) | $- p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{1}{2} \right)$ | $p$ |
| $x_{1}^{1} x_{3}^{5} y_{2}^{2} y_{4}^{5} e_{1} e_{3}$ | (1,1) | (1/2, 0, 1/2, 0) | $- p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{1}{2} \right)$ | $p$ |
| $x_{2}^{1} x_{4}^{10} y_{1}^{1} y_{3}^{5} e_{2} e_{4}$ | (1,1) | (0, 1/3, 0, 2/3) | $- p \, \Gamma_p \left( \frac{1}{3} \right) \Gamma_p \left( \frac{2}{3} \right)$ | $p$ |
| $x_{2}^{2} x_{4}^{5} y_{1}^{1} y_{3}^{5} e_{2} e_{4}$ | (1,1) | (0, 2/3, 0, 1/3) | $- p \, \Gamma_p \left( \frac{2}{3} \right) \Gamma_p \left( \frac{1}{3} \right)$ | $p$ |
| $y_{1}^{1} y_{2}^{2} y_{3}^{5} y_{4}^{5}$ | (1,1) | (0, 0, 0, 0) | $p \,$ | $p$ |
| $y_{1}^{1} y_{2}^{2} y_{3}^{7} y_{4}^{2}$ | (1,1) | (0, 0, 0, 0) | $p \,$ | $p$ |
| $y_{1}^{1} y_{2}^{1} y_{3}^{9} y_{4}^{4}$ | (1,1) | (0, 0, 0, 0) | $p \,$ | $p$ |
| $y_{1}^{1} y_{2}^{1} y_{3}^{3} y_{4}^{13}$ | (1,1) | (0, 0, 0, 0) | $p \,$ | $p$ |
| $y_{1}^{1} y_{2}^{1} y_{3}^{5} y_{4}^{10}$ | (1,1) | (0, 0, 0, 0) | $p \,$ | $p$ |
| $y_{1}^{1} y_{2}^{2} y_{3}^{1} y_{4}^{11}$ | (1,1) | (0, 0, 0, 0) | $p \,$ | $p$ |
| $y_{1}^{1} y_{2}^{1} y_{3}^{7} y_{4}^{7}$ | (1,1) | (0, 0, 0, 0) | $p \,$ | $p$ |
| $y_{1}^{1} y_{2}^{2} y_{3}^{3} y_{4}^{8}$ | (1,1) | (0, 0, 0, 0) | $p \,$ | $p$ |
| $y_{1}^{1} y_{2}^{2} y_{3}^{9} y_{4}^{14}$ | (2,2) | (0, 0, 0, 0) | $p^{2} \,$ | $p^2$ |

### Orbifold vs Goto comparison

Orbifold point count: $N = 3270689$
Orbifold zeta: $(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{14} \cdot (1 - 1873 t + 4068p t^{2} - 5981p^{2} t^{3} + 4595p^{3} t^{4} - 5981p^{4} t^{5} + 4068p^{5} t^{6} - 1873p^{6} t^{7} + p^{8} t^{8})$

Goto point count: $N = 3270689$
Goto zeta: $(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{14} \cdot (1 - 1873 t + 4068p t^{2} - 5981p^{2} t^{3} + 4595p^{3} t^{4} - 5981p^{4} t^{5} + 4068p^{5} t^{6} - 1873p^{6} t^{7} + p^{8} t^{8})$

PASS: Point counts match ($N = 3270689$)
PASS: Zeta denominators match

### Singularity structure

| Pair $(i,j)$ | $(q_i, q_j)$ | $d_{ij}$ | $e_{ij}$ | $f_{ij}$ | $\omega_{ij}$ | $r_{ij}$ |
|------------|------------|------|------|------|----------|------|
| (0,1) | (15,10) | 5 | 30 | 1 | 1 | 4 |
| (0,2) | (15,3) | 3 | 15 | 1 | 2 | 2 |
| (1,3) | (10,2) | 2 | 10 | 1 | 3 | 1 |

### Middle cohomology $P_2(t)$
  Orbifold: $(1 - p t)^{14} \cdot (1 - 1873 t + 4068p t^{2} - 5981p^{2} t^{3} + 4595p^{3} t^{4} - 5981p^{4} t^{5} + 4068p^{5} t^{6} - 1873p^{6} t^{7} + p^{8} t^{8})$
  Goto:     $(1 - p t)^{14} \cdot (1 - 1873 t + 4068p t^{2} - 5981p^{2} t^{3} + 4595p^{3} t^{4} - 5981p^{4} t^{5} + 4068p^{5} t^{6} - 1873p^{6} t^{7} + p^{8} t^{8})$
PASS: Middle cohomology polynomials match

Riemann hypothesis: PASS
