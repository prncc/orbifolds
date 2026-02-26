# Example: Greene-Plesser Mirror Pair of the Fermat Quintic

Fermat quintic: $x_1^5 + x_2^5 + x_3^5 + x_4^5 + x_5^5 = 0$ in $\mathbb{P}^4$
  $A = 5 I_5$, $d = 5$, $p = 11$

## Orbifold 1: $|G| = 125$

  $|G| = 125$
  Cohomology elements: $80$

  $h^{1,1} = 17$
  $h^{2,1} = 21$
  $h^{1,2} = 21$
  $\chi = -8$

PASS: $h^{1,1}=17$, $h^{2,1}=21$, $\chi=-8$

  Numerator $P_3(t)$: degree $44$
  Factored: $(1 - 89 t + 351p t^{2} - 89p^{3} t^{3} + p^{6} t^{4}) \cdot (1 + p t - 9p^{2} t^{2} + p^{4} t^{3} + p^{6} t^{4})^{10}$

  Denominator: $(1 - t) \cdot (1 - p^{3} t) \cdot (1 - p t)^{17} \cdot (1 - p^{2} t)^{17}$

  Riemann hypothesis ($|\alpha|=p^{3/2}$): PASS

## Orbifold 2: $|G^T| = 25$ (mirror)

  $|G^T| = 25$
  Cohomology elements: $80$

  $h^{1,1} = 21$
  $h^{2,1} = 17$
  $h^{1,2} = 17$
  $\chi = 8$

PASS: $h^{1,1}=21$, $h^{2,1}=17$, $\chi=8$

  Numerator $P_3(t)$: degree $36$
  Factored: $(1 - 89 t + 351p t^{2} - 89p^{3} t^{3} + p^{6} t^{4}) \cdot (1 + p t - 9p^{2} t^{2} + p^{4} t^{3} + p^{6} t^{4})^{8}$

  Denominator: $(1 - t) \cdot (1 - p^{3} t) \cdot (1 - p t)^{21} \cdot (1 - p^{2} t)^{21}$

  Riemann hypothesis ($|\alpha|=p^{3/2}$): PASS

## Mirror Symmetry

  $h^{1,1}(G) = 17 = h^{2,1}(G^T) = 17$: PASS
  $h^{2,1}(G) = 21 = h^{1,1}(G^T) = 21$: PASS

PASS: Hodge numbers swap under mirror symmetry

## Zeta Factorization Structure

  $\deg(P_3(G)) = 44$ (expected $44 = 4 + 10 \cdot 4$)
  $\deg(P_3(G^T)) = 36$ (expected $36 = 4 + 8 \cdot 4$)

  $R_1(t) = (1 - 89 t + 351p t^{2} - 89p^{3} t^{3} + p^{6} t^{4})$
  $R_A(pt)^2 = (1 + p t - 9p^{2} t^{2} + p^{4} t^{3} + p^{6} t^{4})$

  PASS: $P_3(G) = R_1 \cdot [R_A(pt)^2]^{10}$
  PASS: $P_3(G^T) = R_1 \cdot [R_A(pt)^2]^8$
