"""
Example from the paper: chain(3,3,3,4) Calabi-Yau K3 surface.

The chain(3,3,3,4) matrix defines a smooth K3 surface in P^3 with
det(A)=108, d=4. Shows:
  - 24 cohomology elements with eigenvalue table
  - Good prime p=109 (det=108 | 108=p-1): orbifold count = brute force = 12147
  - Bad prime p=5 (d=4 | 4 but det=108 does not divide 4): formula fails

Usage: sage examples/example_chain_3334.sage
"""

import os

load("orbifold.sage")
load("orbifold_display.sage")


# --- Setup ---

A = build_chain_matrix([3, 3, 3, 4])
n = A.nrows()
d = get_d(A)
det_A = det(A)
weights = tuple(get_weights(A))
A_inv = A.inverse()

assert n == 4
assert d == 4
assert det_A == 108
assert weights == (1, 1, 1, 1)

J = vector([1] * n)
cy_value = J * A.inverse().T * J
assert cy_value in ZZ, f"CY condition should hold, got {cy_value}"

# Build cohomology basis
G = generate_symmetry_group([[1] * n], A)
all_elements = []
for lam in G:
    for element in get_lambda_sector_basis(A, G, lam, d):
        all_elements.append(element)

sorted_elements = sorted(all_elements, key=lambda e: (e.age_lam, e.age_gamma))
hodge_to_elements = group_elements_by_hodge(all_elements)

assert len(all_elements) == 24

# --- Report ---

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_chain_3334.md")
report_lines = []


def out(s=""):
    report_lines.append(s)
    print(s)


out("# Example: chain(3,3,3,4) Calabi-Yau K3 Surface")
out()
out(f"$A = \\text{{chain}}(3,3,3,4)$")
out(f"  $W_A = x_1^3 x_2 + x_2^3 x_3 + x_3^3 x_4 + x_4^4$")
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
    out(f"| ${latex(element)}$ | {element.y_lambda} | {tuple(element.x_gamma)} | "
        f"{element.age_lam} | {element.age_gamma} | ({s},{r}) | ${eigenval}$ |")
out()

# --- Good prime p=109 ---

p = 109
assert (p - 1) % det_A == 0, f"det(A)={det_A} should divide {p-1}"
R = Zp(p, 200)

out(f"## Good prime: $p = {p}$ ($\\det(A) = {det_A} \\mid {p-1} = p-1$)")
out()

# Trace contributions
out("### Trace contributions")
out()
out("| Element | Hodge | $\\gamma A^{-1}$ | Value |")
out("|---------|-------|--------------|-------|")
for element in sorted_elements:
    s, r = element.hodge_number
    tc = element.trace_contribution(p, R)
    gamma_entries = [Rational(g) for g in vector(element.x_gamma) * A_inv]
    gamma_str = "(" + ", ".join(str(g) for g in gamma_entries) + ")"
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
    out(f"| ${latex(element)}$ | ({s},{r}) | {gamma_str} | {val_str} |")
out()

# Orbifold point count
N_orb = count_points_orbifold(A, d, p, 1, R)
out(f"### Point count (orbifold): $N = {N_orb}$")

# Zeta function
zeta_denom = compute_orbifold_zeta_denominator(all_elements, p, R)
Qt = PolynomialRing(QQ, 't')
t = Qt.gen()
trivial = (1 - t) * (1 - p^2 * t)
middle = zeta_denom / trivial
out(f"### Zeta denominator factored: ${render_factored_poly_latex(zeta_denom, p)}$")
out(f"### Middle $P_2(t)$ factored: ${render_factored_poly_latex(middle, p)}$")
out()

# Riemann hypothesis
rh_results = check_riemann_hypothesis(zeta_denom, p)
rh_ok = all(r[3] for r in rh_results)
out(f"Riemann hypothesis: {'PASS' if rh_ok else 'FAIL'}")
out()

# Brute force
print(f"  Brute force at p={p} (this may take ~1 minute) ...")
N_brute = count_points_projective(A, p)
out(f"### Point count (brute force): $N = {N_brute}$")
out()

assert N_orb == 12147, f"Expected N=12147 at p=109, got {N_orb}"
assert N_brute == 12147, f"Expected N=12147 at p=109, got {N_brute}"
assert N_orb == N_brute
out(f"PASS: orbifold = brute force = ${N_orb}$")
out()

# --- Bad prime p=5 ---

p_bad = 5
assert (p_bad - 1) % d == 0, f"d={d} should divide {p_bad-1}"
assert (p_bad - 1) % det_A != 0, f"det(A)={det_A} should NOT divide {p_bad-1}"

out(f"## Bad prime: $p = {p_bad}$ ($d={d} \\mid {p_bad-1}$ but $\\det(A)={det_A}$ does not divide ${p_bad-1}$)")
out()

R_bad = Zp(p_bad, 20)

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
