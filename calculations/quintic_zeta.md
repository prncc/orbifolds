# Zeta Function Verification: Fermat Quintic Threefold

Surface: $x_1^5 + x_2^5 + x_3^5 + x_4^5 + x_5^5 = 0$ in $\mathbb{P}^4$ ($\psi = 0$)

Matrix: $A = 5 I_5$, $\det(A) = 5^5 = 3125$, $d = 5$

Primes: $37501$, $62501$, $112501$, $118751$ (first 4 with $\det(A) \mid (p-1)$)

## Zeta function structure

For a Calabi--Yau threefold:
$$\zeta_p(t) = \frac{P_3(t)}{(1-t)(1-pt)(1-p^2 t)(1-p^3 t)}$$
where $P_3(t) = R_1(t) \cdot \bigl[R_A(pt, 0)^2\bigr]^{50}$ has degree $4 + 50 \cdot 4 = 204$.

---

## $p = 37501$

### Method 1: Gauss sum formula

Computation time: 5.1s

$R_1(t)$: degree 4 polynomial

$R_A(pt, 0)^2$: degree 4 polynomial

| Parameter | Value |
|-----------|-------|
| $a_1$ | $-8414879$ |
| $b_1$ | $1287051631$ |
| $c$ | $271$ |
| $d$ | $93331$ |
| $\alpha_+$ | $\frac{-8414879 + 7741525\sqrt{5}}{2}$ |
| $\alpha_-$ | $\frac{-8414879 -7741525\sqrt{5}}{2}$ |
| $\gamma_+$ | $\frac{271 + 5\sqrt{5}}{2}$ |
| $\gamma_-$ | $\frac{271 -5\sqrt{5}}{2}$ |

#### Factored $P_3(t)$

$$(1 - 8414879 t + 1287051631p t^{2} - 8414879p^{3} t^{3} + p^{6} t^{4}) \cdot (1 + 271p t + 93331p^{2} t^{2} + 271p^{4} t^{3} + p^{6} t^{4})^{50}$$

#### Complete zeta function

$$\zeta_p(t) = \frac{(1 - 8414879 t + 1287051631p t^{2} - 8414879p^{3} t^{3} + p^{6} t^{4}) \cdot (1 + 271p t + 93331p^{2} t^{2} + 271p^{4} t^{3} + p^{6} t^{4})^{50}}{(1 - t)(1 - p\,t)(1 - p^{2} t)(1 - p^{3} t)}$$

### Method 2: Orbifold trace formula

Computation time: 74.2s

Number of cohomology elements: 208

Orbifold numerator $P_3(t)$: degree 204

Orbifold denominator: degree 4

#### Factored orbifold $P_3(t)$

$$(1 - 8414879 t + 1287051631p t^{2} - 8414879p^{3} t^{3} + p^{6} t^{4}) \cdot (1 + 271p t + 93331p^{2} t^{2} + 271p^{4} t^{3} + p^{6} t^{4})^{50}$$

Orbifold denominator (factored): $(1 - t) \cdot (1 - p t) \cdot (1 - p^{2} t) \cdot (1 - p^{3} t)$

### Comparison

| | Gauss sum | Orbifold |
|---|-----------|----------|
| $\deg P_3$ | 204 | 204 |
| Match | yes | yes |

### Point counts

Extracted from the zeta function expansion $\log Z(t) = \sum_{n \geq 1} N_n \, t^n / n$:

| $n$ | $N_n$ (points over $\mathbb{F}_{37501^n}$) |
|-----|-------|
| 1 | $52740499948675$ |
| 2 | $2781359284579565153342704675$ |
| 3 | $146684977590415830796619713088162291370925$ |

---

## $p = 62501$

### Method 1: Gauss sum formula

Computation time: 8.6s

$R_1(t)$: degree 4 polynomial

$R_A(pt, 0)^2$: degree 4 polynomial

| Parameter | Value |
|-----------|-------|
| $a_1$ | $30690371$ |
| $b_1$ | $9257997381$ |
| $c$ | $71$ |
| $d$ | $99981$ |
| $\alpha_+$ | $\frac{30690371 + 10775725\sqrt{5}}{2}$ |
| $\alpha_-$ | $\frac{30690371 -10775725\sqrt{5}}{2}$ |
| $\gamma_+$ | $\frac{71 + 145\sqrt{5}}{2}$ |
| $\gamma_-$ | $\frac{71 -145\sqrt{5}}{2}$ |

#### Factored $P_3(t)$

$$(1 + 30690371 t + 9257997381p t^{2} + 30690371p^{3} t^{3} + p^{6} t^{4}) \cdot (1 + 71p t + 99981p^{2} t^{2} + 71p^{4} t^{3} + p^{6} t^{4})^{50}$$

#### Complete zeta function

$$\zeta_p(t) = \frac{(1 + 30690371 t + 9257997381p t^{2} + 30690371p^{3} t^{3} + p^{6} t^{4}) \cdot (1 + 71p t + 99981p^{2} t^{2} + 71p^{4} t^{3} + p^{6} t^{4})^{50}}{(1 - t)(1 - p\,t)(1 - p^{2} t)(1 - p^{3} t)}$$

### Method 2: Orbifold trace formula

Computation time: 130.3s

Number of cohomology elements: 208

Orbifold numerator $P_3(t)$: degree 204

Orbifold denominator: degree 4

#### Factored orbifold $P_3(t)$

$$(1 + 30690371 t + 9257997381p t^{2} + 30690371p^{3} t^{3} + p^{6} t^{4}) \cdot (1 + 71p t + 99981p^{2} t^{2} + 71p^{4} t^{3} + p^{6} t^{4})^{50}$$

Orbifold denominator (factored): $(1 - t) \cdot (1 - p t) \cdot (1 - p^{2} t) \cdot (1 - p^{3} t)$

### Comparison

| | Gauss sum | Orbifold |
|---|-----------|----------|
| $\deg P_3$ | 204 | 204 |
| Match | yes | yes |

### Point counts

Extracted from the zeta function expansion $\log Z(t) = \sum_{n \geq 1} N_n \, t^n / n$:

| $n$ | $N_n$ (points over $\mathbb{F}_{62501^n}$) |
|-----|-------|
| 1 | $244156502943925$ |
| 2 | $59610367065473834056333248175$ |
| 3 | $14554010838275253898527781005005357802359425$ |

---

## $p = 112501$

### Method 1: Gauss sum formula

Computation time: 16.3s

$R_1(t)$: degree 4 polynomial

$R_A(pt, 0)^2$: degree 4 polynomial

| Parameter | Value |
|-----------|-------|
| $a_1$ | $17212621$ |
| $b_1$ | $-915109619$ |
| $c$ | $-479$ |
| $d$ | $278581$ |
| $\alpha_+$ | $\frac{17212621 + 49191475\sqrt{5}}{2}$ |
| $\alpha_-$ | $\frac{17212621 -49191475\sqrt{5}}{2}$ |
| $\gamma_+$ | $\frac{-479 + 55\sqrt{5}}{2}$ |
| $\gamma_-$ | $\frac{-479 -55\sqrt{5}}{2}$ |

#### Factored $P_3(t)$

$$(1 + 17212621 t - 915109619p t^{2} + 17212621p^{3} t^{3} + p^{6} t^{4}) \cdot (1 - 479p t + 278581p^{2} t^{2} - 479p^{4} t^{3} + p^{6} t^{4})^{50}$$

#### Complete zeta function

$$\zeta_p(t) = \frac{(1 + 17212621 t - 915109619p t^{2} + 17212621p^{3} t^{3} + p^{6} t^{4}) \cdot (1 - 479p t + 278581p^{2} t^{2} - 479p^{4} t^{3} + p^{6} t^{4})^{50}}{(1 - t)(1 - p\,t)(1 - p^{2} t)(1 - p^{3} t)}$$

### Method 2: Orbifold trace formula

Computation time: 259.1s

Number of cohomology elements: 208

Orbifold numerator $P_3(t)$: degree 204

Orbifold denominator: degree 4

#### Factored orbifold $P_3(t)$

$$(1 + 17212621 t - 915109619p t^{2} + 17212621p^{3} t^{3} + p^{6} t^{4}) \cdot (1 - 479p t + 278581p^{2} t^{2} - 479p^{4} t^{3} + p^{6} t^{4})^{50}$$

Orbifold denominator (factored): $(1 - t) \cdot (1 - p t) \cdot (1 - p^{2} t) \cdot (1 - p^{3} t)$

### Comparison

| | Gauss sum | Orbifold |
|---|-----------|----------|
| $\deg P_3$ | 204 | 204 |
| Match | yes | yes |

### Point counts

Extracted from the zeta function expansion $\log Z(t) = \sum_{n \geq 1} N_n \, t^n / n$:

| $n$ | $N_n$ (points over $\mathbb{F}_{112501^n}$) |
|-----|-------|
| 1 | $1423876073488675$ |
| 2 | $2027394654052389497109812802175$ |
| 3 | $2886738507011079685634568411080155451821512175$ |

---

## $p = 118751$

### Method 1: Gauss sum formula

Computation time: 17.5s

$R_1(t)$: degree 4 polynomial

$R_A(pt, 0)^2$: degree 4 polynomial

| Parameter | Value |
|-----------|-------|
| $a_1$ | $11436371$ |
| $b_1$ | $-14628025869$ |
| $c$ | $521$ |
| $d$ | $275331$ |
| $\alpha_+$ | $\frac{11436371 + 63993725\sqrt{5}}{2}$ |
| $\alpha_-$ | $\frac{11436371 -63993725\sqrt{5}}{2}$ |
| $\gamma_+$ | $\frac{521 + 155\sqrt{5}}{2}$ |
| $\gamma_-$ | $\frac{521 -155\sqrt{5}}{2}$ |

#### Factored $P_3(t)$

$$(1 + 11436371 t - 14628025869p t^{2} + 11436371p^{3} t^{3} + p^{6} t^{4}) \cdot (1 + 521p t + 275331p^{2} t^{2} + 521p^{4} t^{3} + p^{6} t^{4})^{50}$$

#### Complete zeta function

$$\zeta_p(t) = \frac{(1 + 11436371 t - 14628025869p t^{2} + 11436371p^{3} t^{3} + p^{6} t^{4}) \cdot (1 + 521p t + 275331p^{2} t^{2} + 521p^{4} t^{3} + p^{6} t^{4})^{50}}{(1 - t)(1 - p\,t)(1 - p^{2} t)(1 - p^{3} t)}$$

### Method 2: Orbifold trace formula

Computation time: 269.9s

Number of cohomology elements: 208

Orbifold numerator $P_3(t)$: degree 204

Orbifold denominator: degree 4

#### Factored orbifold $P_3(t)$

$$(1 + 11436371 t - 14628025869p t^{2} + 11436371p^{3} t^{3} + p^{6} t^{4}) \cdot (1 + 521p t + 275331p^{2} t^{2} + 521p^{4} t^{3} + p^{6} t^{4})^{50}$$

Orbifold denominator (factored): $(1 - t) \cdot (1 - p t) \cdot (1 - p^{2} t) \cdot (1 - p^{3} t)$

### Comparison

| | Gauss sum | Orbifold |
|---|-----------|----------|
| $\deg P_3$ | 204 | 204 |
| Match | yes | yes |

### Point counts

Extracted from the zeta function expansion $\log Z(t) = \sum_{n \geq 1} N_n \, t^n / n$:

| $n$ | $N_n$ (points over $\mathbb{F}_{118751^n}$) |
|-----|-------|
| 1 | $1674620058737425$ |
| 2 | $2804294711853468324003533172175$ |
| 3 | $4696079921757156471284481243293773424762344675$ |

