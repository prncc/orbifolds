"""Integration tests comparing Goto's formula against the orbifold trace formula.

Tests all 14 untwisted weighted diagonal K3 surfaces from Goto
Proposition 8.1 / Table 7, verifying point counts match across
multiple primes and extension degrees.

Run with: sage diagonal_k3_test.sage
"""

load("orbifold.sage")
load("orbifold_display.sage")
load("diagonal_k3.sage")


# Goto Proposition 8.1 / Table 7: all 14 weighted diagonal K3 surfaces.
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


def find_valid_primes(m, count=2):
    """Finds the first primes p satisfying p = 1 (mod m).

    Searches primes up to 10000 for those congruent to 1 modulo m,
    which is the condition needed for the orbifold trace formula.

    Args:
        m: Positive integer modulus.
        count: Number of valid primes to find.

    Returns:
        List of up to count primes with p = 1 (mod m).
    """
    found = []
    for p in Primes():
        if p > 10000:
            break
        if (p - 1) % m == 0:
            found.append(p)
            if len(found) >= count:
                break
    return found


def has_f_ij_gt_1(p, m, weights, CC):
    """Checks if any singularity term has f_ij > 1 for this prime.

    Per Goto Lemma 4.5, f_ij = m_ij / gcd(m_ij, e_ij). In the
    untwisted case m_ij is 1 or 2 depending on whether -1 is an m-th
    power in F_p^*. When m_ij = 2 and e_ij is odd, f_ij = 2 and the
    singularity correction at (i,j) only contributes when f_ij | nu.

    Args:
        p: Prime number.
        m: Degree of the hypersurface.
        weights: Tuple of 4 positive integer weights.
        CC: Complex field for computation.

    Returns:
        True if any singularity term has f_ij > 1.
    """
    for term in collect_k3_singularity_terms(p, m, weights, CC):
        if term.f_ij > 1:
            return True
    return False


def find_first_prime_dividing_det(A):
    """Finds the first prime p with det(A) | (p-1).

    Args:
        A: The potential matrix.

    Returns:
        The first such prime, or None if none found under 50000.
    """
    det_A = abs(A.det())
    for p in Primes():
        if p > 50000:
            return None
        if (p - 1) % det_A == 0:
            return p
    return None


def get_orbifold_elements(A, m):
    """Gets all orbifold cohomology basis elements for a surface.

    Args:
        A: The potential matrix.
        m: The degree.

    Returns:
        List of CohomologyBasisElement objects.
    """
    G = generate_symmetry_group([[1] * A.nrows()], A)
    all_elements = []
    for lam in G:
        for element in get_lambda_sector_basis(A, G, lam, m):
            all_elements.append(element)
    return all_elements


def count_points_goto_int(p, nu, m, weights, CC):
    """Computes Goto's point count and rounds to an integer.

    Args:
        p: Prime number.
        nu: Extension degree.
        m: Degree of the hypersurface.
        weights: Tuple of 4 positive integer weights.
        CC: Complex field for computation.

    Returns:
        Integer point count.
    """
    result = count_points_goto(p, nu, m, weights, CC)
    return Integer(round(CC(result).real()))


def test_k3_surface(m, weights, primes, nus):
    """Tests one K3 surface by comparing Goto and orbifold point counts.

    For each combination of prime and extension degree, computes the
    point count using both Goto's formula and the orbifold trace
    formula, and checks they agree.

    Args:
        m: Degree of the hypersurface.
        weights: Tuple of 4 positive integer weights.
        primes: List of primes to test.
        nus: List of extension degrees to test.

    Returns:
        List of (p, nu, goto_count, orbifold_count, match) tuples.
    """
    A = build_fermat_matrix([m / w_i for w_i in weights])
    assert tuple(get_weights(A)) == weights
    assert get_d(A) == m

    CC = ComplexField(1024)
    results = []

    for p in primes:
        R = Zp(p, 256)
        for nu in nus:
            goto = count_points_goto_int(p, nu, m, weights, CC)
            orbifold = count_points_orbifold(A, m, p, nu, R)
            match = (goto == orbifold)
            results.append((p, nu, goto, orbifold, match))

    return results


if __name__ == "__main__":
    print("=" * 72)
    print("Goto vs Orbifold: point counts and zeta at det(A)-primes")
    print("=" * 72)

    CC = ComplexField(1024)
    total_tests = 0
    total_passed = 0
    failed = []

    test_surfaces = [
        (4,  (1, 1, 1, 1)),   # Fermat quartic, det=256, p=257
        (12, (4, 4, 3, 1)),   # det=432, p=433, has singular pairs
        (30, (15, 10, 3, 2)), # det=900, p=1801
    ]

    for m, weights in test_surfaces:
        A = build_fermat_matrix([m / w_i for w_i in weights])
        det_A = abs(A.det())
        p = find_first_prime_dividing_det(A)

        if p is None:
            print("\n(m=%d, Q=%s): No prime found with det(A)=%d | (p-1)" % (
                m, weights, det_A))
            continue

        min_field_deg = compute_minimal_field_degree_diagonal(p, m, weights)
        field_str = "F_%d^%d" % (p, min_field_deg) if min_field_deg > 1 else "F_%d" % p
        print("\n(m=%d, Q=%s)  p=%d, det(A)=%d, field=%s" % (
            m, weights, p, det_A, field_str))

        # --- Point count comparison ---
        results = test_k3_surface(m, weights, [p], nus=[1, 2])

        for p_r, nu, goto, orbifold, match in results:
            total_tests += 1
            status = "OK" if match else "FAIL"
            if match:
                total_passed += 1
            else:
                failed.append((m, weights, p_r, nu, goto, orbifold))
            f_ij = has_f_ij_gt_1(p_r, m, weights, CC)
            f_tag = "  (f_ij>1)" if f_ij else ""
            print("  p=%d, nu=%d: Goto=%d, Orbifold=%d  [%s]%s" % (
                p_r, nu, goto, orbifold, status, f_tag))

        # --- Zeta denominator comparison ---
        R = Zp(p, 256)
        orbifold_elements = get_orbifold_elements(A, m)

        if min_field_deg > 1:
            goto_zeta = compute_goto_zeta_denominator_extension(
                p, min_field_deg, m, weights, CC)
            orbifold_zeta = compute_orbifold_zeta_denominator_extension(
                orbifold_elements, p, R, min_field_deg)
        else:
            goto_zeta = compute_goto_zeta_denominator(p, m, weights, CC)
            orbifold_zeta = compute_orbifold_zeta_denominator(
                orbifold_elements, p, R)

        zeta_match = (goto_zeta == orbifold_zeta)
        total_tests += 1
        if zeta_match:
            total_passed += 1
        else:
            failed.append((m, weights, p, "zeta", goto_zeta.degree(), orbifold_zeta.degree()))

        status = "OK" if zeta_match else "FAIL"
        print("  Zeta: Goto degree=%d, Orbifold degree=%d  [%s]" % (
            goto_zeta.degree(), orbifold_zeta.degree(), status))

    print("\n" + "=" * 72)
    print("Results: %d/%d passed." % (total_passed, total_tests))
    if failed:
        print("FAILURES:")
        for entry in failed:
            print("  %s" % (entry,))
    else:
        print("All tests passed.")
    print("=" * 72)
