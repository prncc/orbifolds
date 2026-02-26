"""
Example 3.6 from the paper: chain(2,2,2,3) non-Calabi-Yau chain.

The chain(2,2,2,3) matrix defines a smooth cubic in P^3. The Calabi-Yau
condition fails (J*A^{-T}*J^T = 4/3), causing fractional ages in twisted
sectors. Shows:
  - 8 cohomology elements: 6 with integer ages, 2 with fractional ages
  - At p=73: partial count from integer-age elements = 6*73 = 438
  - True count = 5841, difference = 5403 = 1 + 73 + 73^2 = |P^2(F_73)|

Usage: sage examples/example_3_6_chain_2223.sage
"""

import os
import itertools

load("orbifold.sage")
load("orbifold_display.sage")


# --- Brute-force point counting ---

def apply_W_to_point(A, point):
    total = 0
    for row in A:
        term = 1
        for i, x in enumerate(point):
            term *= x ** row[i]
        total += term
    return total


def count_points_projective(A, p):
    n = A.nrows()
    count = 0
    coords = list(range(p))
    for first_nonzero in range(n):
        prefix = [0] * first_nonzero + [1]
        n_free = n - first_nonzero - 1
        if n_free == 0:
            if apply_W_to_point(A, prefix) % p == 0:
                count += 1
        else:
            for suffix in itertools.product(*([coords] * n_free)):
                point = prefix + list(suffix)
                if apply_W_to_point(A, point) % p == 0:
                    count += 1
    return count


# --- Setup ---

A = build_chain_matrix([2, 2, 2, 3])
n = A.nrows()
d = get_d(A)
det_A = det(A)
weights = tuple(get_weights(A))

assert n == 4
assert d == 3
assert det_A == 24
assert weights == (1, 1, 1, 1)

J = vector([1] * n)
cy_value = J * A.inverse().T * J
assert cy_value not in ZZ, f"CY condition should fail, but J*A^{{-T}}*J^T = {cy_value}"

# Build cohomology basis
G = generate_symmetry_group([[1] * n], A)
all_elements = []
for lam in G:
    for element in get_lambda_sector_basis(A, G, lam, d):
        all_elements.append(element)

sorted_elements = sorted(all_elements, key=lambda e: (e.age_lam, e.age_gamma))
hodge_to_elements = group_elements_by_hodge(all_elements)

assert len(all_elements) == 8

# Classify by integer vs fractional ages
integer_age_elements = [e for e in sorted_elements if e.age_lam in ZZ and e.age_gamma in ZZ]
fractional_age_elements = [e for e in sorted_elements if e.age_lam not in ZZ or e.age_gamma not in ZZ]

assert len(integer_age_elements) == 6
assert len(fractional_age_elements) == 2

# --- Report ---

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_3_6_chain_2223.md")
report_lines = []


def out(s=""):
    report_lines.append(s)
    print(s)


out("# Example 3.6: chain(2,2,2,3) Non-Calabi-Yau Chain")
out()
out(f"$W_A = x_1^2 x_2 + x_2^2 x_3 + x_3^2 x_4 + x_4^3$, $A = \\text{{chain}}(2,2,2,3)$")
out(f"  $n = {n}$, $d = {d}$, $\\det(A) = {det_A}$, weights $Q = {weights}$")
out(f"  CY condition: $J A^{{-T}} J^T = {cy_value}$ (NOT integer, CY FAILS)")
out(f"  $|G| = {len(G)}$")
out()

# --- Cohomology basis ---

out(f"## Orbifold Cohomology Basis: {len(all_elements)} elements")
out(f"  Integer-age elements: {len(integer_age_elements)}")
out(f"  Fractional-age elements: {len(fractional_age_elements)}")
out()

out("| Element | $\\lambda$ | $\\gamma$ | age($\\lambda$) | age($\\gamma$) | Hodge | Integer ages? |")
out("|---------|--------|-------|----------|----------|-------|---------------|")
for element in sorted_elements:
    s, r = element.hodge_number
    int_ages = element.age_lam in ZZ and element.age_gamma in ZZ
    flag = "yes" if int_ages else "NO"

    def render_rational(r):
        if r in ZZ:
            return str(Integer(r))
        return f"{r.numerator()}/{r.denominator()}"

    out(f"| ${latex(element)}$ | {element.y_lambda} | {tuple(element.x_gamma)} | "
        f"{render_rational(element.age_lam)} | {render_rational(element.age_gamma)} | "
        f"({render_rational(s)},{render_rational(r)}) | {flag} |")
out()

# --- Computations at p=73 ---

p = 73
assert (p - 1) % det_A == 0, f"det(A)={det_A} should divide {p-1}"
R = Zp(p, 200)

out(f"## Computations at $p = {p}$ ($\\det(A) = {det_A} \\mid {p-1} = p-1$)")
out()

# Integer-age eigenvalue contributions
out("### Integer-age element eigenvalues")
out()
out("| Element | Hodge | Eigenvalue |")
out("|---------|-------|------------|")
for element in integer_age_elements:
    s, r = element.hodge_number
    alpha = element.twisted_eigenvalue(p, R)
    alpha_int = Integer(alpha.rational_reconstruction())
    out(f"| ${latex(element)}$ | ({Integer(s)},{Integer(r)}) | {alpha_int} |")
out()

# Partial count from integer-age elements
partial_total = R(0)
for element in integer_age_elements:
    tc = element.trace_contribution(p, R)
    partial_total += tc
N_partial = Integer(partial_total.rational_reconstruction())

out(f"### Partial count (integer-age elements only): $N_{{\\text{{partial}}}} = {N_partial}$")

# Check that H^{1,1} contributes 6*p
h11_elements = hodge_to_elements.get((1, 1), [])
h11_integer = [e for e in h11_elements if e.age_lam in ZZ and e.age_gamma in ZZ]
out(f"  $H^{{1,1}}$ has ${len(h11_integer)}$ integer-age elements, each with eigenvalue $p$")
out(f"  Contribution from $H^{{1,1}}$: ${len(h11_integer)} \\cdot {p} = {len(h11_integer) * p}$")
out()

# Brute force
print(f"  Brute force at p={p} ...")
N_true = count_points_projective(A, p)
out(f"### True point count (brute force): $N_{{\\text{{true}}}} = {N_true}$")
out()

# Comparison
diff = N_true - N_partial
trivial = 1 + p + p^2
out(f"### Comparison")
out(f"  $N_{{\\text{{true}}}} - N_{{\\text{{partial}}}} = {N_true} - {N_partial} = {diff}$")
out(f"  $1 + p + p^2 = 1 + {p} + {p}^2 = {trivial}$")
out()

assert diff == trivial, f"Expected difference {trivial}, got {diff}"
out(f"PASS: $N_{{\\text{{true}}}} - N_{{\\text{{partial}}}} = 1 + p + p^2 = {trivial}$")
out(f"  (The fractional-age sectors would contribute exactly $|\\mathbb{{P}}^2(\\mathbb{{F}}_p)|$ if CY held)")
out()

# --- Write report ---

with open(OUTPUT_FILE, "w") as f:
    f.write("\n".join(report_lines))
print(f"\nReport written to {OUTPUT_FILE}")
