"""
Example from the paper: L2L2 K3 surface.

The L2L2 matrix is a block diagonal of two loop L2(3,3) matrices:
  W_A = x_0^3*x_1 + x_0*x_1^3 + x_2^3*x_3 + x_2*x_3^3 = 0 in P^3(1,1,1,1)
with det(A)=64, d=4, G = <J>. Shows:
  - 24 cohomology elements: h^{0,0}=h^{2,2}=h^{0,2}=h^{2,0}=1, h^{1,1}=20
  - All 20 H^{1,1} elements have eigenvalue p
  - Point counts at p=193 (N=40920) and p=257 (N=70680)
  - Zeta function at p=193: (1-t)(1-pt)^20 * (1+190t+p^2*t^2) * (1-p^2*t)
  - p-adic expansions of alpha_{0,2} and alpha_{2,0}

Usage: sage examples/example_L2L2.sage
"""

import os

load("orbifold.sage")
load("orbifold_display.sage")


# --- Setup ---

A = block_diagonal_matrix(
    build_loop_matrix([3, 3]),
    build_loop_matrix([3, 3]),
)
n = A.nrows()
d = get_d(A)
det_A = det(A)
weights = tuple(get_weights(A))

assert n == 4
assert d == 4
assert det_A == 64

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

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_L2L2.md")
report_lines = []


def out(s=""):
    report_lines.append(s)
    print(s)


out("# Example: L2L2 K3 Surface")
out()
out(f"$W_A = x_0^3 x_1 + x_0 x_1^3 + x_2^3 x_3 + x_2 x_3^3$ in $\\mathbb{{P}}^3(1,1,1,1)$")
out(f"  $\\det(A) = {det_A}$, $d = {d}$, $|G| = {len(G)}$")
out()

# --- Cohomology ---

out(f"## Orbifold Cohomology: {len(all_elements)} elements")
out()
for hodge in sorted(hodge_to_elements):
    s, r = hodge
    out(f"  $h^{{{s},{r}}} = {len(hodge_to_elements[hodge])}$")
out()

# --- Point counts and zeta at multiple primes ---

EXPECTED = {193: 40920, 257: 70680, 449: 209880, 577: 343320, 641: 424920}

for p in [193, 257]:
    R = Zp(p, 256)

    out(f"## $p = {p}$")
    out()

    # Eigenvalue table (only at p=193 for brevity)
    if p == 193:
        out("### Eigenvalue table")
        out()
        out("| Element | Hodge | Eigenvalue |")
        out("|---------|-------|------------|")
        for element in sorted_elements:
            s, r = element.hodge_number
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
                # Show p-adic expansion
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
            out(f"| ${latex(element)}$ | ({s},{r}) | {val_str} |")
        out()

        # Verify all H^{1,1} eigenvalues are p
        h11_elements = hodge_to_elements[(1, 1)]
        for el in h11_elements:
            alpha = el.twisted_eigenvalue(p, R)
            assert Integer(alpha.rational_reconstruction()) == p, \
                f"Expected eigenvalue p for H^{{1,1}} element, got {alpha.rational_reconstruction()}"
        out(f"PASS: All ${len(h11_elements)}$ $H^{{1,1}}$ eigenvalues equal $p$")
        out()

        # p-adic expansions of alpha_{0,2} and alpha_{2,0}
        out("### $p$-adic eigenvalue expansions")
        out()
        for hodge in [(0, 2), (2, 0)]:
            el = hodge_to_elements[hodge][0]
            alpha = el.twisted_eigenvalue(p, R)
            v = alpha.ordp()
            coeffs = [Integer(c) for c in alpha.expansion()]
            nterms = min(6, len(coeffs))
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
            out(f"  $\\alpha_{{{hodge[0]},{hodge[1]}}} = {' + '.join(parts)} + \\cdots$")
        out()

    # Point count
    N_orb = count_points_orbifold(A, d, p, 1, R)
    expected = EXPECTED[p]
    out(f"### Point count: $N = {N_orb}$ (expected ${expected}$)")
    assert N_orb == expected, f"Expected N={expected} at p={p}, got {N_orb}"
    out(f"PASS")
    out()

    # Zeta function
    zeta_denom = compute_orbifold_zeta_denominator(all_elements, p, R)
    out(f"### Zeta denominator factored: ${render_factored_poly_latex(zeta_denom, p)}$")
    out()

    # Middle cohomology
    Qt = PolynomialRing(QQ, 't')
    t = Qt.gen()
    trivial = (1 - t) * (1 - p^2 * t)
    middle = zeta_denom / trivial
    out(f"### Middle $P_2(t)$: ${render_factored_poly_latex(middle, p)}$")
    out()

    # Riemann hypothesis
    rh_results = check_riemann_hypothesis(zeta_denom, p)
    rh_ok = all(r[3] for r in rh_results)
    out(f"Riemann hypothesis: {'PASS' if rh_ok else 'FAIL'}")
    out()

# --- Write report ---

with open(OUTPUT_FILE, "w") as f:
    f.write("\n".join(report_lines))
print(f"\nReport written to {OUTPUT_FILE}")
