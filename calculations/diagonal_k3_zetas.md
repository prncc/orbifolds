# Diagonal K3 Zeta Function Report

For each of the 14 diagonal K3 surfaces (Goto Table 7),
we compute the zeta function using both Goto's formula and the orbifold trace formula,
then verify they match.

**Goto zeta denominator:** Computed from Goto's Theorem 5.2 using Jacobi sums and singularity corrections.

**Orbifold zeta denominator:** Computed from the orbifold trace formula (Conjecture 4.1).

**Goto point counts:** Extracted from the power series expansion of the Goto zeta function.

**Orbifold point counts:** Computed directly from the orbifold trace formula.

## Surface $\#1$: $Q=(21, 14, 6, 1)$, $m=42$

Equation: $x_0^{2} + x_1^{3} + x_2^{7} + x_3^{42} = 0$

### $p = 3529$ ($\det A = 1764$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{3529}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{10} \cdot (1 - 4675 t + 13699p t^{2} - 16339p^{2} t^{3} + 26774p^{3} t^{4} - 25606p^{4} t^{5} + 33125p^{5} t^{6} - 25606p^{6} t^{7} + 26774p^{7} t^{8} - 16339p^{8} t^{9} + 13699p^{9} t^{10} - 4675p^{10} t^{11} + p^{12} t^{12})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $12493807$ | $12493807$ | yes |
| $2$ | $155098205359775$ | $155098205359775$ | yes |
| $3$ | $1931567770373312334109$ | $1931567770373312334109$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#2$: $Q=(15, 10, 3, 2)$, $m=30$

Equation: $x_0^{2} + x_1^{3} + x_2^{10} + x_3^{15} = 0$

### $p = 1801$ ($\det A = 900$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{1801}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{14} \cdot (1 - 1873 t + 4068p t^{2} - 5981p^{2} t^{3} + 4595p^{3} t^{4} - 5981p^{4} t^{5} + 4068p^{5} t^{6} - 1873p^{6} t^{7} + p^{8} t^{8})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $3270689$ | $3270689$ | yes |
| $2$ | $10520981712809$ | $10520981712809$ | yes |
| $3$ | $34125755766076001084$ | $34125755766076001084$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#3$: $Q=(12, 8, 3, 1)$, $m=24$

Equation: $x_0^{2} + x_1^{3} + x_2^{8} + x_3^{24} = 0$

### $p = 1153$ ($\det A = 1152$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{1153}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{10} \cdot (1 + p^{2} t^{2})^{2} \cdot (1 - 2912 t + 2080p t^{2} + 2464p^{2} t^{3} - 5246p^{3} t^{4} + 2464p^{4} t^{5} + 2080p^{5} t^{6} - 2912p^{6} t^{7} + p^{8} t^{8})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $1343852$ | $1343852$ | yes |
| $2$ | $1767339949000$ | $1767339949000$ | yes |
| $3$ | $2349502142967849260$ | $2349502142967849260$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#4$: $Q=(10, 5, 4, 1)$, $m=20$

Equation: $x_0^{2} + x_1^{4} + x_2^{5} + x_3^{20} = 0$

### $p = 1601$ ($\det A = 800$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{1601}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{14} \cdot (1 + 6232 t + 13148p t^{2} + 19304p^{2} t^{3} + 21830p^{3} t^{4} + 19304p^{4} t^{5} + 13148p^{5} t^{6} + 6232p^{6} t^{7} + p^{8} t^{8})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $2579384$ | $2579384$ | yes |
| $2$ | $6570031989144$ | $6570031989144$ | yes |
| $3$ | $16840229006482409144$ | $16840229006482409144$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#5$: $Q=(9, 6, 2, 1)$, $m=18$

Equation: $x_0^{2} + x_1^{3} + x_2^{9} + x_3^{18} = 0$

### $p = 2917$ ($\det A = 972$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{2917}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{16} \cdot (1 + 4530 t + 987p t^{2} - 1316p^{2} t^{3} + 987p^{3} t^{4} + 4530p^{4} t^{5} + p^{6} t^{6})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $8551032$ | $8551032$ | yes |
| $2$ | $72401342919288$ | $72401342919288$ | yes |
| $3$ | $616053706694430767160$ | $616053706694430767160$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#6$: $Q=(6, 4, 1, 1)$, $m=12$

Equation: $x_0^{2} + x_1^{3} + x_2^{12} + x_3^{12} = 0$

### $p = 2593$ ($\det A = 864$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{2593}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{6} \cdot (1 + p t)^{8} \cdot (1 + p t + p^{2} t^{2})^{2} \cdot (1 + 850 t - 3405p t^{2} + 850p^{2} t^{3} + p^{4} t^{4})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $6712428$ | $6712428$ | yes |
| $2$ | $45207554939820$ | $45207554939820$ | yes |
| $3$ | $303959065482434371464$ | $303959065482434371464$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#7$: $Q=(6, 3, 2, 1)$, $m=12$

Equation: $x_0^{2} + x_1^{4} + x_2^{6} + x_3^{12} = 0$

### $p = 577$ ($\det A = 576$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{577}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{18} \cdot (1 + 92 t + 966p t^{2} + 92p^{2} t^{3} + p^{4} t^{4})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $343224$ | $343224$ | yes |
| $2$ | $110846605464$ | $110846605464$ | yes |
| $3$ | $36902426197572024$ | $36902426197572024$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#8$: $Q=(4, 4, 3, 1)$, $m=12$

Equation: $x_0^{3} + x_1^{3} + x_2^{4} + x_3^{12} = 0$

### $p = 433$ ($\det A = 432$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{433}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{18} \cdot (1 + 68 t + 294p t^{2} + 68p^{2} t^{3} + p^{4} t^{4})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $195216$ | $195216$ | yes |
| $2$ | $35155249944$ | $35155249944$ | yes |
| $3$ | $6590638235507856$ | $6590638235507856$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#9$: $Q=(4, 3, 3, 2)$, $m=12$

Equation: $x_0^{3} + x_1^{4} + x_2^{4} + x_3^{6} = 0$

### $p = 577$ ($\det A = 288$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{577}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{18} \cdot (1 + 92 t + 966p t^{2} + 92p^{2} t^{3} + p^{4} t^{4})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $343224$ | $343224$ | yes |
| $2$ | $110846605464$ | $110846605464$ | yes |
| $3$ | $36902426197572024$ | $36902426197572024$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#10$: $Q=(5, 2, 2, 1)$, $m=10$

Equation: $x_0^{2} + x_1^{5} + x_2^{5} + x_3^{10} = 0$

### $p = 3001$ ($\det A = 500$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{3001}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{10} \cdot (1 - 9799 t + 14001p t^{2} - 9799p^{2} t^{3} + p^{4} t^{4}) \cdot (1 + p t + p^{2} t^{2} + p^{3} t^{3} + p^{4} t^{4})^{2}$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $9039809$ | $9039809$ | yes |
| $2$ | $81108138046409$ | $81108138046409$ | yes |
| $3$ | $730459215726830632409$ | $730459215726830632409$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#11$: $Q=(4, 2, 1, 1)$, $m=8$

Equation: $x_0^{2} + x_1^{4} + x_2^{8} + x_3^{8} = 0$

### $p = 7681$ ($\det A = 512$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{7681}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 + p t)^{4} \cdot (1 - p t)^{14} \cdot (1 - 4300 t - 5466p t^{2} - 4300p^{2} t^{3} + p^{4} t^{4})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $59078872$ | $59078872$ | yes |
| $2$ | $3480736967431512$ | $3480736967431512$ | yes |
| $3$ | $205355619016224986824792$ | $205355619016224986824792$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#12$: $Q=(2, 2, 1, 1)$, $m=6$

Equation: $x_0^{3} + x_1^{3} + x_2^{6} + x_3^{6} = 0$

### $p = 1297$ ($\det A = 324$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{1297}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{16} \cdot (1 + 478 t + p^{2} t^{2}) \cdot (1 + p t + p^{2} t^{2})^{2}$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $1699890$ | $1699890$ | yes |
| $2$ | $2829847534674$ | $2829847534674$ | yes |
| $3$ | $4760360695111029144$ | $4760360695111029144$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#13$: $Q=(3, 1, 1, 1)$, $m=6$

Equation: $x_0^{2} + x_1^{6} + x_2^{6} + x_3^{6} = 0$

### $p = 433$ ($\det A = 432$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{433}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{20} \cdot (1 + 862 t + p^{2} t^{2})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $195288$ | $195288$ | yes |
| $2$ | $35156242968$ | $35156242968$ | yes |
| $3$ | $6590638254808536$ | $6590638254808536$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

## Surface $\#14$: $Q=(1, 1, 1, 1)$, $m=4$

Equation: $x_0^{4} + x_1^{4} + x_2^{4} + x_3^{4} = 0$

### $p = 257$ ($\det A = 256$ divides $p-1$)

Minimal field of definition: $\mathbb{F}_{257}$

#### Zeta Denominator

$$(1 - t) \cdot (1 - p^{2} t) \cdot (1 - p t)^{20} \cdot (1 + 510 t + p^{2} t^{2})$$

#### Point Counts

| $\nu$ | $N(\mathbb{F}_{q^\nu})$ (Goto) | Orbifold | Match |
|-------|-------------------------------|----------|-------|
| $1$ | $70680$ | $70680$ | yes |
| $2$ | $4363919384$ | $4363919384$ | yes |
| $3$ | $288137115411480$ | $288137115411480$ | yes |

Zeta denominators match: yes

Zeta numerators are 1: yes

Point counts match: yes

---

All 14 surfaces passed.
