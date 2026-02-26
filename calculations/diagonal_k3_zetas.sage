"""
Generate zeta function report for all 14 diagonal K3 surfaces.

Compares Goto's formula (Theorem 5.2) against the orbifold trace formula
for point counts over the minimal field of definition.

Usage: sage analysis/diagonal_zeta_report.sage

Output: analysis/diagonal_zeta_report.md
"""

import os
import sys

OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "diagonal_k3_zetas.md")

load("diagonal_k3.sage")
load("orbifold.sage")
load("orbifold_display.sage")

# All 14 diagonal K3 surfaces from Goto Proposition 8.1 / Table 7
ALL_K3_SURFACES = [
    (42, (21, 14, 6, 1)),
    (30, (15, 10, 3, 2)),
    (24, (12, 8, 3, 1)),
    (20, (10, 5, 4, 1)),
    (18, (9, 6, 2, 1)),
    (12, (6, 4, 1, 1)),
    (12, (6, 3, 2, 1)),
    (12, (4, 4, 3, 1)),
    (12, (4, 3, 3, 2)),
    (10, (5, 2, 2, 1)),
    (8,  (4, 2, 1, 1)),
    (6,  (2, 2, 1, 1)),
    (6,  (3, 1, 1, 1)),
    (4,  (1, 1, 1, 1)),
]


def select_primes_for_surface(weights, m, A, M):
    """Select primes to test for a surface.

    Returns a list containing the first prime p with det(A) | (p-1) and M | (p-1).
    """
    det_prime = find_first_prime_dividing_det(A, M)
    if det_prime is not None:
        return [det_prime]
    return []


def find_first_prime_dividing_det(A, M):
    """Find the first prime p with det(A) | (p - 1) and M | (p - 1)."""
    det_A = abs(A.det())
    target = lcm(det_A, M)
    for p in Primes():
        if p > 50000:
            return None
        if (p - 1) % target == 0:
            return p
    return None


def get_orbifold_elements(A, m):
    """Get all orbifold cohomology basis elements for a surface."""
    G = generate_symmetry_group([[1] * A.nrows()], A)
    all_elements = []
    for lam in G:
        for element in get_lambda_sector_basis(A, G, lam, m):
            all_elements.append(element)
    return all_elements


def point_counts_from_zeta(zeta_poly, q, num_terms=3):
    """Extract point counts N_1, N_2, N_3 from the zeta denominator.

    The zeta function Z(T) = exp(sum_{n>=1} N_n T^n / n).
    So log(Z(T)) = sum N_n T^n / n, which means
    N_n = n * [T^n] log(Z(T)).

    Since Z(T) = 1/denom(T) for our surfaces, log(Z(T)) = -log(denom(T)).
    """
    R = PowerSeriesRing(QQ, 't', default_prec=num_terms + 5)
    t = R.gen()

    # Convert polynomial to power series
    denom_series = R(zeta_poly)

    # log(1/denom) = -log(denom)
    log_zeta = -denom_series.log()

    counts = []
    for n in range(1, num_terms + 1):
        # N_n = n * coefficient of T^n in log(Z)
        N_n = n * log_zeta[n]
        counts.append(Integer(N_n))

    return counts


def format_equation_latex(weights, m):
    """Format the surface equation in LaTeX."""
    exponents = [m // w for w in weights]
    terms = " + ".join(f"x_{i}^{{{e}}}" for i, e in enumerate(exponents))
    return terms + " = 0"


def format_zeta_latex(poly, p):
    """Format a zeta polynomial in factored LaTeX form."""
    return render_factored_poly_latex(poly, p)


def generate_prime_report(idx, weights, m, p, A, M, CC):
    """Generate report section for one prime of a surface."""
    lines = []

    # Compute minimal field of definition
    min_field_deg = compute_minimal_field_degree_diagonal(p, m, weights)
    q = p^min_field_deg

    if min_field_deg > 1:
        field_str = f"$\\mathbb{{F}}_{{{p}^{{{min_field_deg}}}}}$"
    else:
        field_str = f"$\\mathbb{{F}}_{{{p}}}$"

    det_A = abs(A.det())
    det_divides = (p - 1) % det_A == 0
    det_note = f" ($\\det A = {det_A}$ divides $p-1$)" if det_divides else ""

    lines.append(f"### $p = {p}${det_note}")
    lines.append("")
    lines.append(f"Minimal field of definition: {field_str}")
    lines.append("")

    R = Zp(p, 256)
    orbifold_elements = get_orbifold_elements(A, M)

    # Compute zeta denominators
    if min_field_deg > 1:
        goto_zeta = compute_goto_zeta_denominator_extension(
            p, min_field_deg, m, weights, CC)
        orbifold_zeta = compute_orbifold_zeta_denominator_extension(
            orbifold_elements, p, R, min_field_deg)
    else:
        goto_zeta = compute_goto_zeta_denominator(p, m, weights, CC)
        orbifold_zeta = compute_orbifold_zeta_denominator(orbifold_elements, p, R)

    # Compute zeta numerators (should be 1 for K3 surfaces)
    odd_elements = [e for e in orbifold_elements if sum(e.hodge_number) % 2 == 1]
    if min_field_deg > 1:
        orbifold_numerator = _accumulate_zeta_product(odd_elements, p, R, f=min_field_deg)
    else:
        orbifold_numerator = _accumulate_zeta_product(odd_elements, p, R)

    Qt = PolynomialRing(QQ, 't')
    goto_numerator = Qt(1)  # Goto's formula gives numerator = 1 for K3

    lines.append("#### Zeta Denominator")
    lines.append("")
    lines.append(f"$${format_zeta_latex(goto_zeta, p)}$$")
    lines.append("")

    # Check matches
    zeta_match = (goto_zeta == orbifold_zeta)
    numerator_match = (goto_numerator == orbifold_numerator)

    # Compute point counts from zeta power series
    goto_counts = point_counts_from_zeta(goto_zeta, q, num_terms=3)

    # Compute point counts directly
    nus = [min_field_deg * k for k in [1, 2, 3]]
    orbifold_counts = []
    for nu in nus:
        count = count_points_orbifold(A, M, p, nu, R)
        orbifold_counts.append(Integer(count))

    # Check if counts match
    counts_match = (goto_counts == orbifold_counts)

    lines.append("#### Point Counts")
    lines.append("")
    lines.append("| $\\nu$ | $N(\\mathbb{F}_{q^\\nu})$ (Goto) | Orbifold | Match |")
    lines.append("|-------|-------------------------------|----------|-------|")
    for i, (nu, gc, oc) in enumerate(zip(nus, goto_counts, orbifold_counts)):
        match_str = "yes" if gc == oc else "**NO**"
        lines.append(f"| ${nu}$ | ${gc}$ | ${oc}$ | {match_str} |")
    lines.append("")

    all_match = zeta_match and counts_match and numerator_match
    lines.append(f"Zeta denominators match: {'yes' if zeta_match else '**NO**'}")
    lines.append("")
    lines.append(f"Zeta numerators are 1: {'yes' if numerator_match else '**NO**'}")
    lines.append("")
    lines.append(f"Point counts match: {'yes' if counts_match else '**NO**'}")
    lines.append("")

    return lines, all_match


def write_to_file(text, mode='a'):
    """Write text to the output file."""
    with open(OUTPUT_FILE, mode) as f:
        f.write(text)
        f.flush()


def main():
    """Generate full report for all 14 surfaces."""
    CC = ComplexField(1024)

    # Initialize the file with header
    header = """# Diagonal K3 Zeta Function Report

For each of the 14 diagonal K3 surfaces (Goto Table 7),
we compute the zeta function using both Goto's formula and the orbifold trace formula,
then verify they match.

**Goto zeta denominator:** Computed from Goto's Theorem 5.2 using Jacobi sums and singularity corrections.

**Orbifold zeta denominator:** Computed from the orbifold trace formula (Conjecture 4.1).

**Goto point counts:** Extracted from the power series expansion of the Goto zeta function.

**Orbifold point counts:** Computed directly from the orbifold trace formula.

"""
    write_to_file(header, mode='w')
    print(f"Writing report to {OUTPUT_FILE}")

    all_passed = True
    failed_surfaces = []
    num_surfaces = len(ALL_K3_SURFACES)

    for idx, (m, weights) in enumerate(ALL_K3_SURFACES, 1):
        exponents = [m // w for w in weights]
        M = lcm(exponents)
        A = build_fermat_matrix(exponents)

        # Select primes: first valid, one with different min field (if any), det(A) | (p-1)
        primes = select_primes_for_surface(weights, m, A, M)

        if not primes:
            no_prime_text = f"## Surface #{idx}: Q={weights}, m={m}\n\nNo valid primes found for M={M}\n\n"
            write_to_file(no_prime_text)
            print(f"[{idx}/{num_surfaces}] Q={weights}, m={m} - No valid primes for M={M}")
            continue

        # Write surface header
        surface_header = f"## Surface $\\#{idx}$: $Q={weights}$, $m={m}$\n\n"
        surface_header += f"Equation: ${format_equation_latex(weights, m)}$\n\n"
        write_to_file(surface_header)

        for p in primes:
            print(f"[{idx}/{num_surfaces}] Q={weights}, m={m}, p={p} ... ", end="", flush=True)

            lines, passed = generate_prime_report(idx, weights, m, p, A, M, CC)

            if passed:
                print("OK")
            else:
                print("FAIL")
                all_passed = False
                if idx not in failed_surfaces:
                    failed_surfaces.append(idx)

            # Write this prime's results immediately
            write_to_file("\n".join(lines) + "\n")

    # Write footer
    footer = "---\n\n"
    if all_passed:
        footer += "All 14 surfaces passed.\n"
    else:
        footer += f"FAILURES: Surfaces {failed_surfaces}\n"
    write_to_file(footer)

    print(f"\nReport complete: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
