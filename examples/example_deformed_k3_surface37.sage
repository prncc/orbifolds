"""
Example from the paper: Deformed diagonal K3 surface #37.

Surface: x_0^5 + x_1^15 + x_2^3 + x_0*x_3^2 = 0 in P^3(3,1,5,6)
  Q = (3,1,5,6), m = 15, det(A) = 450

Shows:
  - 24 cohomology elements, h^{1,1}=20
  - Orbifold vs deformed Goto at p=1801 (f=1) and p=2251 (f=2, with base change)
  - Singularity structure at different primes
  - Zeta function over F_p and F_{p^2}

Usage: sage examples/example_deformed_k3_surface37.sage
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

weights, m = get_surface(37)
assert weights == (3, 1, 5, 6)
assert m == 15

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

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_deformed_k3_surface37.md")
report_lines = []


def out(s=""):
    report_lines.append(s)
    print(s)


eq_str = surface_equation_string(weights, m)
out(f"# Example: Deformed Diagonal K3 Surface #37")
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

# --- Test at two primes ---

# Find primes where det(A) | (p-1) and m | (p-1)
test_primes = []
for candidate in Primes():
    if len(test_primes) >= 2:
        break
    if candidate > 10000:
        break
    if (candidate - 1) % Integer(det_A) == 0 and (candidate - 1) % m == 0:
        test_primes.append(candidate)

for p in test_primes:
    R = Zp(p, 256)
    f_min = compute_minimal_field_degree(p, m, weights)

    out(f"## $p = {p}$, minimal field degree $f = {f_min}$")
    out()

    # Singularity structure
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

    # --- f=1: over F_p ---

    out(f"### Over $\\mathbb{{F}}_{{{p}}}$ ($\\nu=1$)")
    out()

    N_orb = count_points_orbifold(A, m, p, 1, R)
    N_goto = Integer(round(CC(count_points_deformed_goto(p, 1, m, weights, CC)).real()))
    out(f"  Orbifold: $N = {N_orb}$")
    out(f"  Goto:     $N = {N_goto}$")
    assert N_orb == N_goto, f"Point count mismatch at p={p}, nu=1"
    out(f"  PASS: match")
    out()

    zeta_orb = compute_orbifold_zeta_denominator(all_elements, p, R)
    zeta_goto = compute_deformed_goto_zeta_denominator(p, m, weights, CC)
    out(f"  Orbifold zeta: ${render_factored_poly_latex(zeta_orb, p)}$")
    out(f"  Goto zeta:     ${render_factored_poly_latex(zeta_goto, p)}$")
    assert zeta_orb == zeta_goto, f"Zeta mismatch at p={p}"
    out(f"  PASS: zeta match")
    out()

    # Middle cohomology
    Qt = PolynomialRing(QQ, 't')
    t = Qt.gen()
    trivial = (1 - t) * (1 - p^2 * t)
    middle = zeta_orb / trivial
    out(f"  $P_2(t)$: ${render_factored_poly_latex(middle, p)}$")
    out()

    # Riemann hypothesis
    rh_results = check_riemann_hypothesis(zeta_orb, p)
    rh_ok = all(r[3] for r in rh_results)
    out(f"  Riemann hypothesis: {'PASS' if rh_ok else 'FAIL'}")
    out()

    # --- Base change to F_{p^f} if f > 1 ---

    if f_min > 1:
        f = f_min
        q = p^f
        out(f"### Over $\\mathbb{{F}}_{{p^{f}}} = \\mathbb{{F}}_{{{q}}}$ ($\\nu={f}$)")
        out()

        N_orb_f = count_points_orbifold(A, m, p, f, R)
        N_goto_f = Integer(round(CC(count_points_deformed_goto(p, f, m, weights, CC)).real()))
        out(f"  Orbifold: $N = {N_orb_f}$")
        out(f"  Goto:     $N = {N_goto_f}$")
        assert N_orb_f == N_goto_f, f"Point count mismatch at p={p}, nu={f}"
        out(f"  PASS: match")
        out()

        zeta_orb_f = compute_orbifold_zeta_denominator_extension(all_elements, p, R, f)
        zeta_goto_f = compute_deformed_goto_zeta_denominator_extension(p, f, m, weights, CC)
        out(f"  Orbifold zeta: ${render_factored_poly_latex(zeta_orb_f, q)}$")
        out(f"  Goto zeta:     ${render_factored_poly_latex(zeta_goto_f, q)}$")
        assert zeta_orb_f == zeta_goto_f, f"Zeta mismatch at p={p}, f={f}"
        out(f"  PASS: zeta match")
        out()

# --- Write report ---

with open(OUTPUT_FILE, "w") as f:
    f.write("\n".join(report_lines))
print(f"\nReport written to {OUTPUT_FILE}")
