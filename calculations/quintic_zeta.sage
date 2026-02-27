"""
Zeta function verification for the Fermat quintic threefold at psi = 0.

Computes the zeta function P_3(t) via two independent methods:
  (1) Gauss sum formula (Candelas--de la Ossa--Rodriguez Villegas):
      P_3(t) = R_1(t) * [R_A(pt, 0)^2]^{50}
  (2) Orbifold trace formula (supertrace of twisted Frobenius on Borisov complex)

Then extracts the first 3 point counts N_1, N_2, N_3 from the zeta function.

Usage: sage calculations/quintic_zeta.sage

Output: calculations/quintic_zeta.md
"""

import os
import time

load("orbifold.sage")
load("orbifold_display.sage")
load("quintic_factorization.sage")

# Large primes need more PARI stack for p-adic gamma computations
pari.allocatemem(10**10)

OUTPUT_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "quintic_zeta.md",
)


def write(text, mode='a'):
    with open(OUTPUT_FILE, mode) as f:
        f.write(text)
        f.flush()


def format_sqrt5_latex(elt):
    """Format a Q(sqrt5) element as LaTeX with proper sqrt{5} notation."""
    coeffs = elt.list()
    c0, c1 = QQ(coeffs[0]), QQ(coeffs[1])

    if c1 == 0:
        return latex(c0)

    d0 = c0.denominator()
    d1 = c1.denominator()
    common_denom = lcm(d0, d1)
    a = ZZ(c0 * common_denom)
    b = ZZ(c1 * common_denom)

    if b == 1:
        sqrt_part = r"\sqrt{5}"
    elif b == -1:
        sqrt_part = r"-\sqrt{5}"
    else:
        sqrt_part = r"%d\sqrt{5}" % b

    if a == 0:
        numerator_str = sqrt_part
    elif b > 0:
        numerator_str = f"{a} + {sqrt_part}"
    else:
        numerator_str = f"{a} {sqrt_part}"

    if common_denom == 1:
        return numerator_str
    else:
        return r"\frac{%s}{%s}" % (numerator_str, common_denom)


def get_orbifold_elements(A, m):
    """Get all orbifold cohomology basis elements."""
    G = generate_symmetry_group([[1] * A.nrows()], A)
    all_elements = []
    for lam in G:
        for element in get_lambda_sector_basis(A, G, lam, m):
            all_elements.append(element)
    return all_elements


def point_counts_from_zeta_cy3(P3_poly, p, num_terms=3):
    """Extract point counts N_1, ..., N_{num_terms} from the CY3 zeta function.

    For a CY3: Z(t) = P_3(t) / [(1-t)(1-pt)(1-p^2 t)(1-p^3 t)]
    so log Z(t) = log P_3(t) - log[(1-t)(1-pt)(1-p^2 t)(1-p^3 t)]
    and N_n = n * [t^n] log Z(t).
    """
    R = PowerSeriesRing(QQ, 't', default_prec=num_terms + 5)
    t = R.gen()

    trivial = (1 - t) * (1 - p*t) * (1 - p**2*t) * (1 - p**3*t)
    log_zeta = R(P3_poly).log() - R(trivial).log()

    counts = []
    for n in range(1, num_terms + 1):
        N_n = n * log_zeta[n]
        counts.append(Integer(N_n))
    return counts


def main():
    # Setup
    A = build_fermat_matrix([5, 5, 5, 5, 5])
    M = get_d(A)
    detA = 5**5  # = 3125

    # Find first 4 primes with det(A) | (p - 1)
    NUM_PRIMES = 4
    primes_to_test = []
    for q in Primes():
        if len(primes_to_test) >= NUM_PRIMES:
            break
        if q > 10**7:
            break
        if (q - 1) % detA == 0:
            primes_to_test.append(q)
    assert len(primes_to_test) == NUM_PRIMES, \
        f"Only found {len(primes_to_test)} primes"
    print(f"Testing primes: {primes_to_test}")

    header = rf"""# Zeta Function Verification: Fermat Quintic Threefold

Surface: $x_1^5 + x_2^5 + x_3^5 + x_4^5 + x_5^5 = 0$ in $\mathbb{{P}}^4$ ($\psi = 0$)

Matrix: $A = 5 I_5$, $\det(A) = 5^5 = 3125$, $d = {M}$

Primes: {', '.join(f'${p}$' for p in primes_to_test)} (first {NUM_PRIMES} with $\det(A) \mid (p-1)$)

## Zeta function structure

For a Calabi--Yau threefold:
$$\zeta_p(t) = \frac{{P_3(t)}}{{(1-t)(1-pt)(1-p^2 t)(1-p^3 t)}}$$
where $P_3(t) = R_1(t) \cdot \bigl[R_A(pt, 0)^2\bigr]^{{50}}$ has degree $4 + 50 \cdot 4 = 204$.

"""
    write(header, mode='w')

    # Orbifold elements depend only on A and M, compute once
    orbifold_elements = get_orbifold_elements(A, M)
    print(f"Found {len(orbifold_elements)} cohomology elements")

    for p in primes_to_test:
        print(f"\n{'='*60}")
        print(f"Processing p = {p}")
        print(f"{'='*60}")

        write(f"---\n\n## $p = {p}$\n\n")

        # =====================================================================
        # Method 1: Gauss sum formula
        # =====================================================================
        print("--- Method 1: Gauss sum formula ---")
        t0 = time.time()
        res = extract_candelas_factors_psi_zero(p)
        t_gauss = time.time() - t0

        R1_poly = res['R1']
        QA_poly = res['RA2']  # This is R_A(pt, 0)^2
        P3_gauss = R1_poly * QA_poly**50
        print(f"  Gauss sum: deg(P3) = {P3_gauss.degree()} ({t_gauss:.1f}s)")

        write(f"### Method 1: Gauss sum formula\n\n")
        write(f"Computation time: {t_gauss:.1f}s\n\n")
        write(f"$R_1(t)$: degree {R1_poly.degree()} polynomial\n\n")
        write(f"$R_A(pt, 0)^2$: degree {QA_poly.degree()} polynomial\n\n")
        write(f"| Parameter | Value |\n")
        write(f"|-----------|-------|\n")
        write(f"| $a_1$ | ${res['a1']}$ |\n")
        write(f"| $b_1$ | ${res['b1']}$ |\n")
        write(f"| $c$ | ${res['c']}$ |\n")
        write(f"| $d$ | ${res['d']}$ |\n")
        write(f"| $\\alpha_+$ | ${format_sqrt5_latex(res['alpha_pm'][0])}$ |\n")
        write(f"| $\\alpha_-$ | ${format_sqrt5_latex(res['alpha_pm'][1])}$ |\n")
        write(f"| $\\gamma_+$ | ${format_sqrt5_latex(res['gamma_pm'][0])}$ |\n")
        write(f"| $\\gamma_-$ | ${format_sqrt5_latex(res['gamma_pm'][1])}$ |\n")
        write(f"\n")

        P3_factored = render_factored_poly_latex(P3_gauss, p)
        write(f"#### Factored $P_3(t)$\n\n")
        write(f"$${P3_factored}$$\n\n")

        denom_latex = r"(1 - t)(1 - p\,t)(1 - p^{2} t)(1 - p^{3} t)"
        write(f"#### Complete zeta function\n\n")
        write(f"$$\\zeta_p(t) = \\frac{{{P3_factored}}}{{{denom_latex}}}$$\n\n")

        # =====================================================================
        # Method 2: Orbifold trace formula
        # =====================================================================
        print("--- Method 2: Orbifold trace formula ---")
        t0 = time.time()
        R_padic = Zp(p, 256)

        numerator, denominator = compute_orbifold_zeta_numerator_denominator(
            orbifold_elements, p, R_padic)
        t_orb = time.time() - t0

        P3_orbifold = numerator
        print(f"  Orbifold: deg(P3) = {P3_orbifold.degree()} ({t_orb:.1f}s)")

        write(f"### Method 2: Orbifold trace formula\n\n")
        write(f"Computation time: {t_orb:.1f}s\n\n")
        write(f"Number of cohomology elements: {len(orbifold_elements)}\n\n")
        write(f"Orbifold numerator $P_3(t)$: degree {P3_orbifold.degree()}\n\n")
        write(f"Orbifold denominator: degree {denominator.degree()}\n\n")

        orb_factored = render_factored_poly_latex(P3_orbifold, p)
        write(f"#### Factored orbifold $P_3(t)$\n\n")
        write(f"$${orb_factored}$$\n\n")

        orb_denom_factored = render_factored_poly_latex(denominator, p)
        write(f"Orbifold denominator (factored): ${orb_denom_factored}$\n\n")

        # =====================================================================
        # Comparison
        # =====================================================================
        match = P3_gauss == P3_orbifold
        print(f"\n=== P_3(t) match: {match} ===")

        write(f"### Comparison\n\n")
        write(f"| | Gauss sum | Orbifold |\n")
        write(f"|---|-----------|----------|\n")
        write(f"| $\\deg P_3$ | {P3_gauss.degree()} | {P3_orbifold.degree()} |\n")
        write(f"| Match | {'yes' if match else '**NO**'} | {'yes' if match else '**NO**'} |\n")
        write(f"\n")

        # =====================================================================
        # Point counts from zeta expansion
        # =====================================================================
        print("\n--- Extracting point counts ---")
        counts = point_counts_from_zeta_cy3(P3_gauss, p, num_terms=3)
        for i, N in enumerate(counts, 1):
            print(f"  N_{i} = {N}")

        write(f"### Point counts\n\n")
        write(f"Extracted from the zeta function expansion ")
        write(f"$\\log Z(t) = \\sum_{{n \\geq 1}} N_n \\, t^n / n$:\n\n")
        write(f"| $n$ | $N_n$ (points over $\\mathbb{{F}}_{{{p}^n}}$) |\n")
        write(f"|-----|-------|\n")
        for i, N in enumerate(counts, 1):
            write(f"| {i} | ${N}$ |\n")
        write(f"\n")

    print(f"\nReport written to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
