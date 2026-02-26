"""
Example from the paper: Deformed diagonal K3 surface #16.

Surface: x_0^10 + x_1^10 + x_2^2 + x_0*x_3^3 = 0 in P^3(1,1,5,3)
  Q = (1,1,5,3), m = 10, det(A) = 600

Shows:
  - 24 cohomology elements, h^{1,1}=20
  - Isolated singularity of type A_{3,2} at [0:0:0:1]
  - Orbifold zeta matches deformed Goto zeta at p=601
  - P_2(t) factorization

Usage: sage examples/example_deformed_k3_surface16.sage
"""

import os

load("orbifold.sage")
load("orbifold_display.sage")
load("deformed_diagonal_k3.sage")
load("deformed_diagonal_k3_enumeration.sage")


def build_deformed_matrix(weights, m):
    """Build the 4x4 potential matrix for a deformed diagonal surface."""
    exps = get_exponents(weights, m)
    m_0, m_1, m_2, m_3 = exps
    return matrix(ZZ, [
        [m_0, 0, 0, 0],
        [0, m_1, 0, 0],
        [0, 0, m_2, 0],
        [1, 0, 0, m_3],
    ])


# --- Setup ---

weights, m = get_surface(16)
assert weights == (1, 1, 5, 3)
assert m == 10

exps = get_exponents(weights, m)
A = build_deformed_matrix(weights, m)
det_A = det(A)

CC = ComplexField(1024)

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

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_deformed_k3_surface16.md")
report_lines = []


def out(s=""):
    report_lines.append(s)
    print(s)


eq_str = surface_equation_string(weights, m)
out(f"# Example: Deformed Diagonal K3 Surface #16")
out()
out(f"Surface: ${eq_str} = 0$ in $\\mathbb{{P}}^3{weights}$")
out(f"  Exponents: ${exps}$")
out(f"  $Q = {weights}$, $m = {m}$, $\\det(A) = {det_A}$")
out(f"  $|G| = {len(G)}$")
out()

# --- Cohomology ---

out(f"## Orbifold Cohomology: {len(all_elements)} elements")
out()
for hodge in sorted(hodge_to_elements):
    s, r = hodge
    out(f"  $h^{{{s},{r}}} = {len(hodge_to_elements[hodge])}$")
out()

# --- Find test prime ---

p = None
for candidate in Primes():
    if candidate > 10000:
        break
    if (candidate - 1) % Integer(det_A) == 0 and (candidate - 1) % m == 0:
        p = candidate
        break

assert p == 601, f"Expected first valid prime 601, got {p}"
R = Zp(p, 256)

f_min = compute_minimal_field_degree(p, m, weights)

out(f"## Computations at $p = {p}$ ($\\det(A) = {det_A} \\mid {p-1}$)")
out(f"  Minimal field degree: $f = {f_min}$")
out()

# --- Singularity structure ---

out("### Singularity structure")
out()
isolated_term, curve_terms = collect_all_singularity_terms(p, m, weights, CC)
if isolated_term is not None:
    out(f"  Isolated singularity: $\\alpha={isolated_term.alpha_3}$, "
        f"$q_3={isolated_term.q_3}$, $r={isolated_term.r_3}$")
for sing in curve_terms:
    i, j = sing.ij
    out(f"  Pair $({i},{j})$: $d={sing.d_ij}$, $\\alpha={sing.alpha_ij}$, "
        f"$e={sing.e_ij}$, $f={sing.f_ij}$, $\\omega={sing.omega_ij}$, $r={sing.r_ij}$")
out()

# --- Orbifold vs Goto comparison ---

out("### Orbifold vs Deformed Goto comparison")
out()

# Orbifold
N_orb = count_points_orbifold(A, m, p, 1, R)
zeta_orb = compute_orbifold_zeta_denominator(all_elements, p, R)
out(f"Orbifold point count: $N = {N_orb}$")
out(f"Orbifold zeta: ${render_factored_poly_latex(zeta_orb, p)}$")
out()

# Goto
N_goto = Integer(round(CC(count_points_deformed_goto(p, 1, m, weights, CC)).real()))
zeta_goto = compute_deformed_goto_zeta_denominator(p, m, weights, CC)
out(f"Deformed Goto point count: $N = {N_goto}$")
out(f"Deformed Goto zeta: ${render_factored_poly_latex(zeta_goto, p)}$")
out()

assert N_orb == N_goto, f"Point count mismatch: orbifold={N_orb}, Goto={N_goto}"
assert zeta_orb == zeta_goto, f"Zeta mismatch!"
out(f"PASS: Point counts match ($N = {N_orb}$)")
out(f"PASS: Zeta denominators match")
out()

# --- Middle cohomology ---

Qt = PolynomialRing(QQ, 't')
t = Qt.gen()
trivial = (1 - t) * (1 - p^2 * t)
middle = zeta_orb / trivial
out(f"### Middle cohomology $P_2(t)$")
out(f"  ${render_factored_poly_latex(middle, p)}$")
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
