"""
Example from the paper: Greene-Plesser mirror pair of the Fermat quintic.

Two quotient orbifolds of the Fermat quintic x_1^5+...+x_5^5 = 0 at p=11:
  (1) G = <J, (0,1,2,3,4), (0,1,1,4,4)>, |G|=125
      80 elements, h^{1,1}=17, h^{2,1}=21, chi=-8
      P_3(t) = R_1(t) * [R_A(pt)^2]^10, numerator degree 44
  (2) G^T = <J, (0,4,1,1,4)>, |G^T|=25
      80 elements, h^{1,1}=21, h^{2,1}=17, chi=8
      P_3(t) = R_1(t) * [R_A(pt)^2]^8, numerator degree 36
Mirror symmetry: Hodge numbers swap.

Usage: sage examples/example_greene_plesser.sage
"""

import os

load("orbifold.sage")
load("orbifold_display.sage")


# --- Setup ---

A = build_fermat_matrix([5, 5, 5, 5, 5])
m = get_d(A)
p = 11
R = Zp(p, 256)

assert m == 5
assert (p - 1) % m == 0

# --- Report ---

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_greene_plesser.md")
report_lines = []


def out(s=""):
    report_lines.append(s)
    print(s)


def compute_orbifold_data(generators):
    """Compute orbifold cohomology and zeta data for a symmetry group."""
    G = generate_symmetry_group(generators, A)
    all_elements = []
    for lam in G:
        for element in get_lambda_sector_basis(A, G, lam, m):
            all_elements.append(element)
    hodge_to_elements = group_elements_by_hodge(all_elements)
    numerator, denominator = compute_orbifold_zeta_numerator_denominator(all_elements, p, R)
    return {
        "G": G,
        "all_elements": all_elements,
        "hodge_to_elements": hodge_to_elements,
        "numerator": numerator,
        "denominator": denominator,
    }


def hodge_count(hodge_to_elements, s, r):
    return len(hodge_to_elements.get((s, r), []))


def euler_characteristic(hodge_to_elements):
    chi = 0
    for (s, r), elements in hodge_to_elements.items():
        chi += (-1)**(s + r) * len(elements)
    return chi


out("# Example: Greene-Plesser Mirror Pair of the Fermat Quintic")
out()
out(f"Fermat quintic: $x_1^5 + x_2^5 + x_3^5 + x_4^5 + x_5^5 = 0$ in $\\mathbb{{P}}^4$")
out(f"  $A = 5 I_5$, $d = {m}$, $p = {p}$")
out()

# ================================================================
# Orbifold 1: |G| = 125
# ================================================================

out("## Orbifold 1: $|G| = 125$")
out()

gen1 = [[1, 1, 1, 1, 1], [0, 1, 2, 3, 4], [0, 1, 1, 4, 4]]
print("Computing |G|=125 orbifold...")
data1 = compute_orbifold_data(gen1)

G1 = data1["G"]
h2e_1 = data1["hodge_to_elements"]
P3_1 = data1["numerator"]
denom_1 = data1["denominator"]

assert len(G1) == 125, f"Expected |G|=125, got {len(G1)}"
out(f"  $|G| = {len(G1)}$")
out(f"  Cohomology elements: ${len(data1['all_elements'])}$")
out()

h11_1 = hodge_count(h2e_1, 1, 1)
h21_1 = hodge_count(h2e_1, 2, 1)
h12_1 = hodge_count(h2e_1, 1, 2)
chi_1 = euler_characteristic(h2e_1)

out(f"  $h^{{1,1}} = {h11_1}$")
out(f"  $h^{{2,1}} = {h21_1}$")
out(f"  $h^{{1,2}} = {h12_1}$")
out(f"  $\\chi = {chi_1}$")
out()

assert h11_1 == 17, f"Expected h^{{1,1}}=17, got {h11_1}"
assert h21_1 == 21, f"Expected h^{{2,1}}=21, got {h21_1}"
assert chi_1 == -8, f"Expected chi=-8, got {chi_1}"
out(f"PASS: $h^{{1,1}}=17$, $h^{{2,1}}=21$, $\\chi=-8$")
out()

out(f"  Numerator $P_3(t)$: degree ${P3_1.degree()}$")
out(f"  Factored: ${render_factored_poly_latex(P3_1, p)}$")
out()
out(f"  Denominator: ${render_factored_poly_latex(denom_1, p)}$")
out()

# Riemann hypothesis on numerator (check each irreducible factor)
rh_results_1 = check_riemann_hypothesis(P3_1, p, max_weight=3)
rh_ok_1 = all(r[2] is not None for r in rh_results_1)
out(f"  Riemann hypothesis ($|\\alpha|=p^{{3/2}}$): {'PASS' if rh_ok_1 else 'FAIL'}")
out()

# ================================================================
# Orbifold 2: |G^T| = 25
# ================================================================

out("## Orbifold 2: $|G^T| = 25$ (mirror)")
out()

gen2 = [[1, 1, 1, 1, 1], [0, 4, 1, 1, 4]]
print("Computing |G^T|=25 orbifold...")
data2 = compute_orbifold_data(gen2)

G2 = data2["G"]
h2e_2 = data2["hodge_to_elements"]
P3_2 = data2["numerator"]
denom_2 = data2["denominator"]

assert len(G2) == 25, f"Expected |G^T|=25, got {len(G2)}"
out(f"  $|G^T| = {len(G2)}$")
out(f"  Cohomology elements: ${len(data2['all_elements'])}$")
out()

h11_2 = hodge_count(h2e_2, 1, 1)
h21_2 = hodge_count(h2e_2, 2, 1)
h12_2 = hodge_count(h2e_2, 1, 2)
chi_2 = euler_characteristic(h2e_2)

out(f"  $h^{{1,1}} = {h11_2}$")
out(f"  $h^{{2,1}} = {h21_2}$")
out(f"  $h^{{1,2}} = {h12_2}$")
out(f"  $\\chi = {chi_2}$")
out()

assert h11_2 == 21, f"Expected h^{{1,1}}=21, got {h11_2}"
assert h21_2 == 17, f"Expected h^{{2,1}}=17, got {h21_2}"
assert chi_2 == 8, f"Expected chi=8, got {chi_2}"
out(f"PASS: $h^{{1,1}}=21$, $h^{{2,1}}=17$, $\\chi=8$")
out()

out(f"  Numerator $P_3(t)$: degree ${P3_2.degree()}$")
out(f"  Factored: ${render_factored_poly_latex(P3_2, p)}$")
out()
out(f"  Denominator: ${render_factored_poly_latex(denom_2, p)}$")
out()

# Riemann hypothesis
rh_results_2 = check_riemann_hypothesis(P3_2, p, max_weight=3)
rh_ok_2 = all(r[2] is not None for r in rh_results_2)
out(f"  Riemann hypothesis ($|\\alpha|=p^{{3/2}}$): {'PASS' if rh_ok_2 else 'FAIL'}")
out()

# ================================================================
# Mirror symmetry check
# ================================================================

out("## Mirror Symmetry")
out()
out(f"  $h^{{1,1}}(G) = {h11_1} = h^{{2,1}}(G^T) = {h21_2}$: {'PASS' if h11_1 == h21_2 else 'FAIL'}")
out(f"  $h^{{2,1}}(G) = {h21_1} = h^{{1,1}}(G^T) = {h11_2}$: {'PASS' if h21_1 == h11_2 else 'FAIL'}")
out()

assert h11_1 == h21_2
assert h21_1 == h11_2
out("PASS: Hodge numbers swap under mirror symmetry")
out()

# ================================================================
# Zeta factorization structure
# ================================================================

out("## Zeta Factorization Structure")
out()

# Extract R_1 and R_A(pt)^2 from factorizations
# P3_1 = R_1 * (R_A^2)^10, degree 4 + 40 = 44
# P3_2 = R_1 * (R_A^2)^8, degree 4 + 32 = 36

assert P3_1.degree() == 44, f"Expected deg(P3_1)=44, got {P3_1.degree()}"
assert P3_2.degree() == 36, f"Expected deg(P3_2)=36, got {P3_2.degree()}"
out(f"  $\\deg(P_3(G)) = {P3_1.degree()}$ (expected $44 = 4 + 10 \\cdot 4$)")
out(f"  $\\deg(P_3(G^T)) = {P3_2.degree()}$ (expected $36 = 4 + 8 \\cdot 4$)")
out()

# Factor P3_1: should be (unit) * R_1^1 * (R_A^2)^10
f1 = P3_1.factor()
f2 = P3_2.factor()

# Collect quartic factors and multiplicities from P3_1
quartics_1 = {fac: mult for fac, mult in f1 if fac.degree() == 4}
quartics_2 = {fac: mult for fac, mult in f2 if fac.degree() == 4}

# Identify R_1 (appears with mult 1 in P3_1) and R_A^2 (appears with mult 10)
R1 = None
RA2_single = None
for fac, mult in quartics_1.items():
    if mult == 1:
        R1 = fac
    elif mult == 10:
        RA2_single = fac

assert R1 is not None, f"Could not find R_1: {f1}"
assert RA2_single is not None, f"Could not find R_A^2: {f1}"

# Verify R_A^2 appears with mult 8 in P3_2
assert RA2_single in quartics_2, "R_A^2 not found in P3_2"
assert quartics_2[RA2_single] == 8, f"Expected mult 8 in P3_2, got {quartics_2[RA2_single]}"
assert quartics_2.get(R1, 0) == 1, f"Expected R_1 with mult 1 in P3_2"

out(f"  $R_1(t) = {render_factored_poly_latex(R1, p)}$")
out(f"  $R_A(pt)^2 = {render_factored_poly_latex(RA2_single, p)}$")
out()
out(f"  PASS: $P_3(G) = R_1 \\cdot [R_A(pt)^2]^{{10}}$")
out(f"  PASS: $P_3(G^T) = R_1 \\cdot [R_A(pt)^2]^8$")
out()

# --- Write report ---

with open(OUTPUT_FILE, "w") as f:
    f.write("\n".join(report_lines))
print(f"\nReport written to {OUTPUT_FILE}")
