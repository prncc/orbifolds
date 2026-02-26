"""
Example 3.3 from the paper: chain(2,2,3) Calabi-Yau curve.

Orbifold trace formula for the CY curve W_A = x_1^2*x_2 + x_2^2*x_3 + x_3^3
in P^2 with G = <J>. Shows:
  - 4 cohomology elements with eigenvalue table
  - Good prime p=13 (det(A)=12 | 12=p-1): orbifold count matches brute force (N=20)
  - Zeta function: (1 + 6t + pt^2) / ((1-t)(1-pt))
  - Bad prime p=7 (det(A)=12 does not divide 6=p-1): formula fails

Usage: sage examples/example_3_3_chain_223.sage
"""

import os
import itertools

load("orbifold.sage")
load("orbifold_display.sage")


# --- Brute-force point counting ---

def apply_W_to_point(A, point):
    """Evaluate W_A(x) = sum of monomials defined by rows of A."""
    total = 0
    for row in A:
        term = 1
        for i, x in enumerate(point):
            term *= x ** row[i]
        total += term
    return total


def count_points_projective(A, p):
    """Count projective points on {W_A = 0} in P^{n-1}(F_p)."""
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

A = build_chain_matrix([2, 2, 3])
d = get_d(A)
n = A.nrows()
det_A = det(A)
weights = tuple(get_weights(A))
A_inv = A.inverse()

assert n == 3
assert d == 3
assert det_A == 12
assert weights == (1, 1, 1)

J = vector([1] * n)
cy_value = J * A.inverse().T * J
assert cy_value in ZZ, "CY condition should hold"

# Build cohomology basis
G = generate_symmetry_group([[1, 1, 1]], A)
all_elements = []
for lam in G:
    for element in get_lambda_sector_basis(A, G, lam, d):
        all_elements.append(element)

sorted_elements = sorted(all_elements, key=lambda e: (e.age_lam, e.age_gamma))
hodge_to_elements = group_elements_by_hodge(all_elements)

assert len(all_elements) == 4

# --- Report generation ---

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_3_3_chain_223.md")
report_lines = []


def out(s=""):
    report_lines.append(s)
    print(s)


out("# Example 3.3: chain(2,2,3) Calabi-Yau Curve")
out()
out(f"$W_A = x_1^2 x_2 + x_2^2 x_3 + x_3^3$, $A = \\text{{chain}}(2,2,3)$")
out(f"  $n = {n}$, $d = {d}$, $\\det(A) = {det_A}$, weights $Q = {weights}$")
out(f"  CY condition: $J A^{{-T}} J^T = {cy_value}$ (integer, CY holds)")
out(f"  $|G| = {len(G)}$")
out()

# --- Cohomology basis ---

out(f"## Orbifold Cohomology Basis: {len(all_elements)} elements")
out()
for hodge in sorted(hodge_to_elements):
    s, r = hodge
    out(f"  $h^{{{s},{r}}} = {len(hodge_to_elements[hodge])}$")
out()

out("| Element | $\\lambda$ | $\\gamma$ | age($\\lambda$) | age($\\gamma$) | Hodge | Eigenvalue |")
out("|---------|--------|-------|----------|----------|-------|------------|")
for element in sorted_elements:
    s, r = element.hodge_number
    eigenval = render_eigenvalue(element, simplify=True)
    gamma_entries = [Rational(g) for g in vector(element.x_gamma) * A_inv]
    gamma_str = "(" + ", ".join(str(g) for g in gamma_entries) + ")"
    out(f"| ${latex(element)}$ | {element.y_lambda} | {tuple(element.x_gamma)} | "
        f"{element.age_lam} | {element.age_gamma} | ({s},{r}) | ${eigenval}$ |")
out()

# --- Good prime p=13 ---

p = 13
assert (p - 1) % det_A == 0, f"det(A)={det_A} should divide {p-1}"
R = Zp(p, 200)

out(f"## Good prime: $p = {p}$ ($\\det(A) = {det_A} \\mid {p-1} = p-1$)")
out()

# Trace contributions
out("### Trace contributions")
out()
out("| Element | Hodge | Contribution | $\\gamma A^{-1}$ | Value |")
out("|---------|-------|--------------|--------------|-------|")
for element in sorted_elements:
    s, r = element.hodge_number
    tc = element.trace_contribution(p, R)
    gamma_entries = [Rational(g) for g in vector(element.x_gamma) * A_inv]
    gamma_str = "(" + ", ".join(str(g) for g in gamma_entries) + ")"
    contrib = render_contribution(element, simplify=True)
    try:
        val = Integer(tc.rational_reconstruction())
        val_str = str(val)
    except (ArithmeticError, TypeError):
        v = tc.ordp()
        coeffs = [Integer(c) for c in tc.expansion()]
        nterms = min(5, len(coeffs))
        parts = []
        for idx in range(nterms):
            c = coeffs[idx]
            k = v + idx
            if c == 0:
                continue
            if k == 0:
                parts.append(f"{c}")
            elif k == 1:
                parts.append(f"{c} \\cdot {p}" if c != 1 else f"{p}")
            else:
                parts.append(f"{c} \\cdot {p}^{k}" if c != 1 else f"{p}^{k}")
        val_str = "$" + " + ".join(parts) + " + \\cdots$"
    out(f"| ${latex(element)}$ | ({s},{r}) | ${contrib}$ | {gamma_str} | {val_str} |")
out()

# Orbifold point count
N_orb = count_points_orbifold(A, d, p, 1, R)
out(f"### Point count (orbifold): $N = {N_orb}$")

# Brute force
N_brute = count_points_projective(A, p)
out(f"### Point count (brute force): $N = {N_brute}$")
out()

assert N_orb == 20, f"Expected N=20 at p=13, got {N_orb}"
assert N_brute == 20, f"Expected N=20 at p=13, got {N_brute}"
assert N_orb == N_brute
out(f"PASS: orbifold = brute force = ${N_orb}$")
out()

# Zeta function
zeta_denom = compute_orbifold_zeta_denominator(all_elements, p, R)
out(f"### Zeta function denominator: ${zeta_denom}$")

Qt = PolynomialRing(QQ, 't')
t = Qt.gen()
trivial = (1 - t) * (1 - p * t)
middle = zeta_denom / trivial
out(f"### Middle cohomology $P_1(t) = {middle}$")
out(f"### Factored: ${render_factored_poly_latex(zeta_denom, p)}$")
out()

# Verify P_1(t) = 1 + 6t + pt^2
expected_middle = 1 + 6*t + p*t^2
assert middle == expected_middle, f"Expected {expected_middle}, got {middle}"
out(f"PASS: $P_1(t) = 1 + 6t + pt^2$")
out()

# Riemann hypothesis
rh_results = check_riemann_hypothesis(zeta_denom, p)
rh_ok = all(r[3] for r in rh_results)
out(f"Riemann hypothesis: {'PASS' if rh_ok else 'FAIL'}")
out()

# --- Bad prime p=7 ---

p_bad = 7
assert (p_bad - 1) % d == 0, f"d={d} should divide {p_bad-1}"
assert (p_bad - 1) % det_A != 0, f"det(A)={det_A} should NOT divide {p_bad-1}"

out(f"## Bad prime: $p = {p_bad}$ ($d={d} \\mid {p_bad-1}$ but $\\det(A)={det_A}$ does not divide ${p_bad-1}$)")
out()

R_bad = Zp(p_bad, 20)

# Compute supertrace
orbifold_total = R_bad(0)
for element in all_elements:
    tc = element.trace_contribution(p_bad, R_bad)
    orbifold_total += tc

try:
    orbifold_count = Integer(orbifold_total.rational_reconstruction())
    out(f"Orbifold count: ${orbifold_count}$ (unexpectedly succeeded)")
except (ArithmeticError, TypeError, ValueError):
    out(f"Orbifold count: rational_reconstruction FAILS (as expected)")
    out(f"  Trace formula value: ${orbifold_total}$")

N_brute_bad = count_points_projective(A, p_bad)
out(f"Brute force count: $N = {N_brute_bad}$")
out()

# --- Write report ---

with open(OUTPUT_FILE, "w") as f:
    f.write("\n".join(report_lines))
print(f"\nReport written to {OUTPUT_FILE}")
