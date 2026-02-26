"""
Example from the paper: Fermat quintic threefold.

Surface: x_1^5 + x_2^5 + x_3^5 + x_4^5 + x_5^5 = 0 in P^4
  A = 5*I_5, d=5, det(A)=3125, G=<J>, |G|=5
  208 cohomology elements, b_3 = 204

Verifies the orbifold trace formula against the Gauss sum method
(Candelas-de la Ossa-Rodriguez Villegas) at primes with det(A) | (p-1).

Shows:
  - P_3(t) = R_1(t) * [R_A(pt)^2]^50 has degree 204
  - Coefficients (a_1, b_1, c, d) at each prime
  - Point counts N_1, N_2, N_3 from zeta expansion
  - Orbifold P_3 matches Gauss sum P_3

Set FULL_TEST = True to test all 4 primes (~1 hour).
Default tests only the first prime p=37501 (~5-15 minutes).

Usage: sage examples/example_fermat_quintic.sage
"""

import os
import time

load("orbifold.sage")
load("orbifold_display.sage")
load("quintic_factorization.sage")

# Large primes need more PARI stack
pari.allocatemem(10**10)

FULL_TEST = False  # Set True to test all 4 primes (takes ~1 hour)


def point_counts_from_zeta_cy3(P3_poly, p, num_terms=3):
    """Extract point counts N_1, ..., N_{num_terms} from the CY3 zeta function."""
    PS = PowerSeriesRing(QQ, 't', default_prec=num_terms + 5)
    t = PS.gen()
    trivial = (1 - t) * (1 - p*t) * (1 - p**2*t) * (1 - p**3*t)
    log_zeta = PS(P3_poly).log() - PS(trivial).log()
    counts = []
    for n in range(1, num_terms + 1):
        N_n = n * log_zeta[n]
        counts.append(Integer(N_n))
    return counts


# --- Setup ---

A = build_fermat_matrix([5, 5, 5, 5, 5])
m = get_d(A)
det_A = 5**5  # = 3125

assert m == 5
assert det(A) == det_A

# Find primes
NUM_PRIMES = 4 if FULL_TEST else 1
primes_to_test = []
for q in Primes():
    if len(primes_to_test) >= NUM_PRIMES:
        break
    if q > 10**6:
        break
    if (q - 1) % det_A == 0:
        primes_to_test.append(q)

# Build orbifold elements (depends only on A and m)
print("Building orbifold cohomology basis...")
G = generate_symmetry_group([[1] * A.nrows()], A)
all_elements = []
for lam in G:
    for element in get_lambda_sector_basis(A, G, lam, m):
        all_elements.append(element)

hodge_to_elements = group_elements_by_hodge(all_elements)

assert len(all_elements) == 208

# Expected values from the paper (Table 3.13)
EXPECTED = {
    37501: {"a1": -8414879, "b1": 1287051631, "c": 271, "d": 93331,
            "N1": 52740499948675},
    62501: {"a1": 30690371, "b1": 9257997381, "c": 71, "d": 99981,
            "N1": 244156502943925},
    112501: {"a1": 17212621, "b1": -915109619, "c": -479, "d": 278581,
             "N1": 1423876073488675},
    118751: {"a1": 11436371, "b1": -14628025869, "c": 521, "d": 275331,
             "N1": 1674620058737425},
}

# --- Report ---

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_fermat_quintic.md")
report_lines = []


def out(s=""):
    report_lines.append(s)
    print(s)


out("# Example: Fermat Quintic Threefold")
out()
out(f"Surface: $x_1^5 + x_2^5 + x_3^5 + x_4^5 + x_5^5 = 0$ in $\\mathbb{{P}}^4$")
out(f"  $A = 5 I_5$, $d = {m}$, $\\det(A) = {det_A}$")
out(f"  Cohomology elements: ${len(all_elements)}$")
out(f"  Primes to test: ${primes_to_test}$")
out()

for hodge in sorted(hodge_to_elements):
    s, r = hodge
    out(f"  $h^{{{s},{r}}} = {len(hodge_to_elements[hodge])}$")
out()

n_odd = sum(len(e) for (s, r), e in hodge_to_elements.items() if (s + r) % 2 == 1)
n_even = sum(len(e) for (s, r), e in hodge_to_elements.items() if (s + r) % 2 == 0)
out(f"  Numerator elements (odd $s+r$): ${n_odd} = b_3$")
out(f"  Denominator elements (even $s+r$): ${n_even}$")
out()

for p in primes_to_test:
    out(f"## $p = {p}$")
    out()

    # --- Gauss sum method ---
    print(f"  p={p}: Gauss sum method...")
    t0 = time.time()
    res = extract_candelas_factors_psi_zero(p)
    t_gauss = time.time() - t0

    R1_poly = res['R1']
    RA2_poly = res['RA2']
    P3_gauss = R1_poly * RA2_poly**50
    out(f"### Gauss sum method ({t_gauss:.1f}s)")
    out(f"  $a_1 = {res['a1']}$")
    out(f"  $b_1 = {res['b1']}$")
    out(f"  $c = {res['c']}$")
    out(f"  $d = {res['d']}$")
    out(f"  $\\deg(P_3) = {P3_gauss.degree()}$")
    out()

    # --- Orbifold method ---
    print(f"  p={p}: Orbifold method...")
    t0 = time.time()
    R_padic = Zp(p, 256)
    numerator, denominator = compute_orbifold_zeta_numerator_denominator(all_elements, p, R_padic)
    t_orb = time.time() - t0

    P3_orb = numerator
    out(f"### Orbifold method ({t_orb:.1f}s)")
    out(f"  $\\deg(P_3) = {P3_orb.degree()}$")
    out()

    # --- Comparison ---
    match = P3_gauss == P3_orb
    out(f"### Comparison")
    out(f"  $P_3$(Gauss) $= P_3$(Orbifold): {'PASS' if match else 'FAIL'}")
    assert match, f"P_3 mismatch at p={p}!"
    out()

    # --- Verify against paper values ---
    if p in EXPECTED:
        exp = EXPECTED[p]
        assert res['a1'] == exp['a1'], f"a1 mismatch at p={p}"
        assert res['b1'] == exp['b1'], f"b1 mismatch at p={p}"
        assert res['c'] == exp['c'], f"c mismatch at p={p}"
        assert res['d'] == exp['d'], f"d mismatch at p={p}"
        out(f"  PASS: $(a_1, b_1, c, d)$ match expected values")
        out()

    # --- Point counts ---
    counts = point_counts_from_zeta_cy3(P3_gauss, p, num_terms=3)
    out(f"### Point counts from zeta expansion")
    for i, N in enumerate(counts, 1):
        out(f"  $N_{i} = {N}$")
    out()

    if p in EXPECTED:
        assert counts[0] == EXPECTED[p]['N1'], f"N_1 mismatch at p={p}"
        out(f"  PASS: $N_1$ matches expected value")
        out()

# --- Write report ---

with open(OUTPUT_FILE, "w") as f:
    f.write("\n".join(report_lines))
print(f"\nReport written to {OUTPUT_FILE}")
