"""
Example from the paper: Diagonal K3 surface with m=30, Q=(15,10,3,2).

Verifies the orbifold trace formula against Goto's explicit Jacobi sum formula
for the weighted diagonal K3 surface x_0^2 + x_1^3 + x_2^10 + x_3^15 = 0
in P^3(15,10,3,2). Shows:
  - 24 cohomology elements, h^{1,1}=20
  - Orbifold zeta matches Goto zeta at p=1801
  - Singularity structure: P_{01}=A_{5,4}, P_{02}=A_{3,2}, P_{13}=A_{2,1}
  - middle cohomology P_2(t) and Riemann hypothesis

Usage: sage examples/example_diagonal_k3_m30.sage
"""

import os

load("orbifold.sage")
load("orbifold_display.sage")
load("diagonal_k3.sage")


# --- Setup ---

m = 30
WEIGHTS = (15, 10, 3, 2)
A = build_fermat_matrix([m / w for w in WEIGHTS])

assert tuple(get_weights(A)) == WEIGHTS
assert get_d(A) == m
assert det(A) == 900

A_inv = A.inverse()

# Build cohomology basis
G = generate_symmetry_group([[1] * A.nrows()], A)
all_elements = []
for lam in G:
    for element in get_lambda_sector_basis(A, G, lam, m):
        all_elements.append(element)

sorted_elements = sorted(all_elements, key=lambda e: (e.age_lam, e.age_gamma))
hodge_to_elements = group_elements_by_hodge(all_elements)

assert len(all_elements) == 24

# --- Report ---

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_diagonal_k3_m30.md")
report_lines = []


def out(s=""):
    report_lines.append(s)
    print(s)


out("# Example: Diagonal K3 Surface (m=30, Q=(15,10,3,2))")
out()
out(f"Surface: $x_0^2 + x_1^3 + x_2^{{10}} + x_3^{{15}} = 0$ in $\\mathbb{{P}}^3(15,10,3,2)$")
out(f"  $A = \\text{{diag}}(2, 3, 10, 15)$, $\\det(A) = {det(A)}$, $d = {m}$")
out(f"  $|G| = {len(G)}$")
out()

# --- Cohomology ---

out(f"## Orbifold Cohomology: {len(all_elements)} elements")
out()
for hodge in sorted(hodge_to_elements):
    s, r = hodge
    out(f"  $h^{{{s},{r}}} = {len(hodge_to_elements[hodge])}$")
out()

# --- Eigenvalue table ---

p = 1801
assert (p - 1) % Integer(det(A)) == 0
R = Zp(p, 256)
CC = ComplexField(1024)

out(f"## Computations at $p = {p}$")
out()

out("### Eigenvalue table")
out()
out("| Element | Hodge | $\\gamma A^{-1}$ | Eigenvalue formula | tc value |")
out("|---------|-------|--------------|-------------------|----------|")
for element in sorted_elements:
    s, r = element.hodge_number
    eigenval = render_eigenvalue(element, simplify=True)
    gamma_entries = [Rational(g) for g in vector(element.x_gamma) * A_inv]
    gamma_str = "(" + ", ".join(str(g) for g in gamma_entries) + ")"
    tc = element.trace_contribution(p, R)
    try:
        val = Integer(tc.rational_reconstruction())
        if val == p:
            val_str = "$p$"
        elif val == -p:
            val_str = "$-p$"
        elif val == p^2:
            val_str = "$p^2$"
        elif val == 1:
            val_str = "$1$"
        else:
            val_str = f"${val}$"
    except (ArithmeticError, TypeError):
        v = tc.ordp()
        coeffs = [Integer(c) for c in tc.expansion()]
        nterms = min(4, len(coeffs))
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
    out(f"| ${latex(element)}$ | ({s},{r}) | {gamma_str} | ${eigenval}$ | {val_str} |")
out()

# --- Orbifold vs Goto comparison ---

out("### Orbifold vs Goto comparison")
out()

# Orbifold
N_orb = count_points_orbifold(A, m, p, 1, R)
zeta_orb = compute_orbifold_zeta_denominator(all_elements, p, R)
out(f"Orbifold point count: $N = {N_orb}$")
out(f"Orbifold zeta: ${render_factored_poly_latex(zeta_orb, p)}$")
out()

# Goto
N_goto = Integer(round(CC(count_points_goto(p, 1, m, WEIGHTS, CC)).real()))
zeta_goto = compute_goto_zeta_denominator(p, m, WEIGHTS, CC)
out(f"Goto point count: $N = {N_goto}$")
out(f"Goto zeta: ${render_factored_poly_latex(zeta_goto, p)}$")
out()

assert N_orb == N_goto, f"Point count mismatch: orbifold={N_orb}, Goto={N_goto}"
assert zeta_orb == zeta_goto, f"Zeta mismatch!"
out(f"PASS: Point counts match ($N = {N_orb}$)")
out(f"PASS: Zeta denominators match")
out()

# --- Singularity structure ---

out("### Singularity structure")
out()

sing_terms = collect_k3_singularity_terms(p, m, WEIGHTS, CC)
out("| Pair $(i,j)$ | $(q_i, q_j)$ | $d_{ij}$ | $e_{ij}$ | $f_{ij}$ | $\\omega_{ij}$ | $r_{ij}$ |")
out("|------------|------------|------|------|------|----------|------|")
for sing in sing_terms:
    i, j = sing.ij
    out(f"| ({i},{j}) | ({WEIGHTS[i]},{WEIGHTS[j]}) | {sing.d_ij} | "
        f"{sing.e_ij} | {sing.f_ij} | {sing.omega_ij} | {sing.r_ij} |")
out()

# --- Middle cohomology ---

Qt = PolynomialRing(QQ, 't')
t = Qt.gen()
trivial = (1 - t) * (1 - p^2 * t)
middle_orb = zeta_orb / trivial
middle_goto = zeta_goto / trivial
out(f"### Middle cohomology $P_2(t)$")
out(f"  Orbifold: ${render_factored_poly_latex(middle_orb, p)}$")
out(f"  Goto:     ${render_factored_poly_latex(middle_goto, p)}$")
assert middle_orb == middle_goto
out(f"PASS: Middle cohomology polynomials match")
out()

# Riemann hypothesis
rh_results = check_riemann_hypothesis(zeta_orb, p)
rh_ok = all(r[3] for r in rh_results)
out(f"Riemann hypothesis: {'PASS' if rh_ok else 'FAIL'}")
out()

# --- Write report ---

with open(OUTPUT_FILE, "w") as f:
    f.write("\n".join(report_lines))
print(f"\nReport written to {OUTPUT_FILE}")
