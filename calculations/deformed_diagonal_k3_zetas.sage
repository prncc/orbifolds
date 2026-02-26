"""
Goto vs Orbifold zeta function comparison for all 85 deformed diagonal K3 surfaces.

For each surface, compares Goto's zeta denominator (Theorem 4.3.4) against
the orbifold trace formula (Conjecture 4.1) at selected (p, f) pairs:

  - First prime p with det(A) | (p-1): test at (p, 1)
  - First prime p with M | (p-1) and f > 1: test at (p, 1) and (p, f)

Usage: sage calculations/deformed_diagonal_zetas.sage

Output: calculations/deformed_diagonal_zetas.md
"""

import os

OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "deformed_diagonal_k3_zetas.md")

load("deformed_diagonal_k3.sage")
load("deformed_diagonal_k3_enumeration.sage")
load("orbifold.sage")
load("orbifold_display.sage")


def build_deformed_matrix(weights, m):
    """Build the 4x4 potential matrix for a deformed diagonal surface."""
    exps = get_exponents(weights, m)
    m_0, m_1, m_2, m_3 = exps
    return matrix([
        [m_0, 0, 0, 0],
        [0, m_1, 0, 0],
        [0, 0, m_2, 0],
        [1, 0, 0, m_3]
    ])


def get_orbifold_elements(A, M):
    """Get all orbifold cohomology basis elements for a surface."""
    G = generate_symmetry_group([[1] * A.nrows()], A)
    all_elements = []
    for lam in G:
        for element in get_lambda_sector_basis(A, G, lam, M):
            all_elements.append(element)
    return all_elements


def has_fgt1_prime(m, weights, M, max_p=100000):
    """Check if there exists a prime p with M | (p-1) and f > 1."""
    for p in Primes():
        if p > max_p:
            return False
        if (p - 1) % M != 0:
            continue
        if compute_minimal_field_degree(p, m, weights) > 1:
            return True
    return False


def find_test_primes(m, weights, M, det_A, max_p=10000):
    """Yield distinct (p, f) pairs to test.

    - First prime with det(A) | (p-1): yield (p, 1)
    - First prime with M | (p-1) and f > 1: yield (p, 1) and (p, f)
    - Deduplicates: if a pair was already yielded, skip it.
    - Stops once both categories are found.
    """
    found_det = False
    found_fgt1 = False
    seen = set()
    for p in Primes():
        if p > max_p or (found_det and found_fgt1):
            return
        if (p - 1) % det_A != 0:
            continue
        if (p - 1) % M != 0:
            continue
        f = compute_minimal_field_degree(p, m, weights)
        if not found_det:
            if (p, Integer(1)) not in seen:
                seen.add((p, Integer(1)))
                yield (p, Integer(1))
            found_det = True
        if not found_fgt1 and f > 1:
            if (p, Integer(1)) not in seen:
                seen.add((p, Integer(1)))
                yield (p, Integer(1))
            if (p, Integer(f)) not in seen:
                seen.add((p, Integer(f)))
                yield (p, Integer(f))
            found_fgt1 = True


def compute_zeta_pair(p, f, m, weights, A, M, CC):
    """Compute Goto and Orbifold zeta denominators over F_{p^f}.

    Returns (goto_denom, orbifold_denom) as polynomials in QQ['t'].
    """
    R = Zp(p, 256)
    orbifold_elements = get_orbifold_elements(A, M)

    if f > 1:
        goto_denom = compute_deformed_goto_zeta_denominator_extension(
            p, f, m, weights, CC)
        orbifold_denom = compute_orbifold_zeta_denominator_extension(
            orbifold_elements, p, R, f)
    else:
        goto_denom = compute_deformed_goto_zeta_denominator(
            p, m, weights, CC)
        orbifold_denom = compute_orbifold_zeta_denominator(
            orbifold_elements, p, R)

    return goto_denom, orbifold_denom


def extract_P2(zeta_denominator, q):
    """Extract P_2(t) from the full zeta denominator.

    The denominator = (1-t) * (1-q^2*t) * P_2(t).
    """
    Qt = PolynomialRing(QQ, 't')
    t = Qt.gen()
    trivial = (1 - t) * (1 - q**2 * t)
    P2, remainder = zeta_denominator.quo_rem(trivial)
    if remainder != 0:
        raise ValueError("Denominator doesn't factor as expected, remainder=%s" % remainder)
    return P2


def format_test_latex(goto_denom, orbifold_denom, q, p):
    """Format LaTeX output for a test result.

    Returns list of markdown lines with factored and expanded forms
    for both the Goto and Orbifold zeta denominators.
    """
    lines = []

    for label, denom in [("Goto", goto_denom), ("Orbifold", orbifold_denom)]:
        P2 = extract_P2(denom, q)
        p2_factored = render_factored_poly_latex(P2, p)
        full_factored = render_factored_poly_latex(denom, p)

        lines.append("**%s:**" % label)
        lines.append("")
        lines.append("$P_2(t) = %s$" % p2_factored)
        lines.append("")
        lines.append("$\\text{denom}(t) = %s$" % full_factored)
        lines.append("")
        lines.append("$\\text{denom}(t) = %s$" % latex(denom))
        lines.append("")

    return lines


def write_to_file(text, mode='a'):
    with open(OUTPUT_FILE, mode) as fh:
        fh.write(text)
        fh.flush()


def write_summary(results):
    """Write summary table and statistics."""
    lines = []
    lines.append("---")
    lines.append("")
    lines.append("## Summary")
    lines.append("")

    all_tests = []
    for r in results:
        all_tests.extend(r['tests'])

    total = len(all_tests)
    passed = sum(1 for t in all_tests if t['passed'])
    failed = sum(1 for t in all_tests if not t['passed'])

    lines.append(f"**{passed}/{total} tests passed, {failed} failed.**")
    lines.append("")

    failures = [(r, t) for r in results for t in r['tests'] if not t['passed']]
    if failures:
        lines.append("### Failures")
        lines.append("")
        for r, t in failures:
            lines.append(f"- Surface #{r['idx']}: Q={r['weights']}, "
                         f"m={r['m']}, p={t['p']}, f={t['f']}")
        lines.append("")

    lines.append("### Full Results")
    lines.append("")
    lines.append("| # | Q | m | det(A) | tests |")
    lines.append("|---|---|---|--------|-------|")
    for r in results:
        if r['tests']:
            test_strs = []
            for t in r['tests']:
                status = "pass" if t['passed'] else "FAIL"
                test_strs.append("(p=%s, f=%s): %s" % (t['p'], t['f'], status))
            tests_col = "; ".join(test_strs)
        else:
            tests_col = "-"
        lines.append(f"| {r['idx']} | {r['weights']} | {r['m']} | "
                     f"{r['det_A']} | {tests_col} |")
    lines.append("")

    write_to_file("\n".join(lines))


def main():
    CC = ComplexField(1024)

    header = """# Deformed Diagonal K3: Goto vs Orbifold Zeta over F_p and F_{p^f}

For each of the 85 deformed diagonal K3 surfaces, we compare the
Goto zeta denominator (Theorem 4.3.4) with the Orbifold zeta
denominator (Conjecture 4.1).

For each surface we select (p, f) pairs:
- First prime p with det(A) | (p-1): test over F_p
- First prime p with M | (p-1) and f > 1: test over F_p and F_{p^f}
- Duplicate pairs are skipped.

"""
    write_to_file(header, mode='w')
    print(f"Writing report to {OUTPUT_FILE}")

    results = []

    for idx in range(1, 86):
        weights, m = get_surface(idx)
        exps = get_exponents(weights, m)
        M = lcm(exps)
        A = build_deformed_matrix(weights, m)
        det_A = abs(A.det())

        # # Skip surfaces with no f>1 prime
        # if not has_fgt1_prime(m, weights, M):
        #     print("[%s/85] Q=%s, m=%s  skipped (no f>1 prime)" % (idx, weights, m))
        #     continue

        result = {
            'idx': idx,
            'weights': weights,
            'm': m,
            'det_A': det_A,
            'M': M,
            'tests': [],
        }

        lines = []
        eq_str = surface_equation_string(weights, m)
        lines.append(f"## Surface #{idx}: Q={weights}, m={m}")
        lines.append("")
        lines.append(f"Equation: {eq_str}")
        lines.append("")
        lines.append(f"det(A) = {det_A}, M = {M}")
        lines.append("")

        test_pairs = list(find_test_primes(m, weights, M, det_A))
        if not test_pairs:
            lines.append("No suitable test primes found (up to 100000)")
            lines.append("")
            print("[%s/85] Q=%s, m=%s  no primes found" % (idx, weights, m))

        for i, (p, f) in enumerate(test_pairs):
            q = p**f
            print("[%s/85] Q=%s, m=%s  p=%s, f=%s ... " % (idx, weights, m, p, f),
                  end="", flush=True)
            goto_denom, orbifold_denom = compute_zeta_pair(
                p, f, m, weights, A, M, CC)
            match = (goto_denom == orbifold_denom)
            result['tests'].append({'p': p, 'f': f, 'passed': match})

            lines.append("### Test %s: p=%s, f=%s" % (i + 1, p, f))
            lines.append("")
            lines.append("Match: %s" % ("yes" if match else "**NO**"))
            lines.append("")
            lines.extend(format_test_latex(goto_denom, orbifold_denom, q, p))
            print("OK" if match else "FAIL")

        write_to_file("\n".join(lines) + "\n\n")
        results.append(result)

    write_summary(results)
    print(f"\nReport complete: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
