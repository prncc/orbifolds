# Example 3.3: chain(2,2,3) Calabi-Yau Curve

$W_A = x_1^2 x_2 + x_2^2 x_3 + x_3^3$, $A = \text{chain}(2,2,3)$
  $n = 3$, $d = 3$, $\det(A) = 12$, weights $Q = (1, 1, 1)$
  CY condition: $J A^{-T} J^T = 1$ (integer, CY holds)
  $|G| = 3$

## Orbifold Cohomology Basis: 4 elements

  $h^{0,0} = 1$
  $h^{0,1} = 1$
  $h^{1,0} = 1$
  $h^{1,1} = 1$

| Element | $\lambda$ | $\gamma$ | age($\lambda$) | age($\gamma$) | Hodge | Eigenvalue |
|---------|--------|-------|----------|----------|-------|------------|
| $x_{1}^{1} x_{2}^{1} x_{3}^{1}  e_{1} e_{2} e_{3}$ | (0, 0, 0) | (1, 1, 1) | 0 | 1 | (0,1) | $- \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{1}{4} \right) \Gamma_p \left( \frac{1}{4} \right)$ |
| $x_{1}^{1} x_{2}^{2} x_{3}^{3}  e_{1} e_{2} e_{3}$ | (0, 0, 0) | (1, 2, 3) | 0 | 2 | (1,0) | $p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{3}{4} \right) \Gamma_p \left( \frac{3}{4} \right)$ |
| $y_{1}^{1} y_{2}^{1} y_{3}^{1}$ | (1, 1, 1) | (0, 0, 0) | 1 | 0 | (0,0) | $1$ |
| $y_{1}^{2} y_{2}^{2} y_{3}^{2}$ | (2, 2, 2) | (0, 0, 0) | 2 | 0 | (1,1) | $p \,$ |

## Good prime: $p = 13$ ($\det(A) = 12 \mid 12 = p-1$)

### Trace contributions

| Element | Hodge | Contribution | $\gamma A^{-1}$ | Value |
|---------|-------|--------------|--------------|-------|
| $x_{1}^{1} x_{2}^{1} x_{3}^{1}  e_{1} e_{2} e_{3}$ | (0,1) | $\Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{1}{4} \right) \Gamma_p \left( \frac{1}{4} \right)$ | (1/2, 1/4, 1/4) | $6 + 2 \cdot 13 + 10 \cdot 13^2 + 12 \cdot 13^3 + 2 \cdot 13^4 + \cdots$ |
| $x_{1}^{1} x_{2}^{2} x_{3}^{3}  e_{1} e_{2} e_{3}$ | (1,0) | $- p \, \Gamma_p \left( \frac{1}{2} \right) \Gamma_p \left( \frac{3}{4} \right) \Gamma_p \left( \frac{3}{4} \right)$ | (1/2, 3/4, 3/4) | $11 \cdot 13 + 2 \cdot 13^2 + 10 \cdot 13^4 + 10 \cdot 13^5 + \cdots$ |
| $y_{1}^{1} y_{2}^{1} y_{3}^{1}$ | (0,0) | $1$ | (0, 0, 0) | 1 |
| $y_{1}^{2} y_{2}^{2} y_{3}^{2}$ | (1,1) | $p \,$ | (0, 0, 0) | 13 |

### Point count (orbifold): $N = 20$
### Point count (brute force): $N = 20$

PASS: orbifold = brute force = $20$

### Zeta function denominator: $169*t^4 - 104*t^3 - 58*t^2 - 8*t + 1$
### Middle cohomology $P_1(t) = 13*t^2 + 6*t + 1$
### Factored: $(1 - t) \cdot (1 - p t) \cdot (1 + 6 t + p t^{2})$

PASS: $P_1(t) = 1 + 6t + pt^2$

Riemann hypothesis: PASS

## Bad prime: $p = 7$ ($d=3 \mid 6$ but $\det(A)=12$ does not divide $6$)

Orbifold count: rational_reconstruction FAILS (as expected)
  Trace formula value: $6*7^2 + 7^3 + 5*7^4 + 5*7^5 + 5*7^6 + 4*7^7 + 2*7^8 + 7^10 + 5*7^11 + 7^12 + 3*7^13 + 2*7^14 + 6*7^15 + 6*7^16 + 6*7^17 + 6*7^18 + 2*7^19 + O(7^20)$
Brute force count: $N = 8$
