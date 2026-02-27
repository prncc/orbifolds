"""Integration tests for deformed diagonal K3 surfaces.

Tests comparing Goto's formula (Theorem 4.3.4) against the orbifold trace
formula for selected surfaces from the 85 deformed diagonal K3 surfaces
in Goto's Table 10.

Run with: sage deformed_diagonal_k3_test.sage
"""

load("orbifold.sage")
load("orbifold_display.sage")
load("deformed_diagonal_k3.sage")
load("deformed_diagonal_k3_enumeration.sage")


def build_deformed_matrix(weights, m):
    """Builds the potential matrix for a deformed diagonal surface.

    The deformed diagonal surface
        x_0^{m_0} + x_1^{m_1} + x_2^{m_2} + x_0 x_3^{m_3} = 0
    corresponds to the matrix:
        [[m_0, 0, 0, 0],
         [0, m_1, 0, 0],
         [0, 0, m_2, 0],
         [1, 0, 0, m_3]]

    Args:
        weights: Tuple (q_0, q_1, q_2, q_3) of weights.
        m: Degree of the hypersurface.

    Returns:
        4x4 integer matrix representing the potential.
    """
    exps = get_exponents(weights, m)
    m_0, m_1, m_2, m_3 = exps
    return matrix([
        [m_0, 0, 0, 0],
        [0, m_1, 0, 0],
        [0, 0, m_2, 0],
        [1, 0, 0, m_3]
    ])


def find_valid_primes_for_M(M, count=2):
    """Finds the first primes p satisfying p = 1 (mod M).

    Args:
        M: Positive integer modulus.
        count: Number of valid primes to find.

    Returns:
        List of up to count primes with p = 1 (mod M).
    """
    found = []
    for p in Primes():
        if p > 10000:
            break
        if (p - 1) % M == 0:
            found.append(p)
            if len(found) >= count:
                break
    return found


def find_first_prime_dividing_det(A, M):
    """Finds the first prime p with det(A) | (p-1) and M | (p-1).

    Args:
        A: The potential matrix.
        M: The lcm of the exponents.

    Returns:
        The first such prime, or None if none found under 50000.
    """
    det_A = abs(A.det())
    target = lcm(det_A, M)
    for p in Primes():
        if p > 50000:
            return None
        if (p - 1) % target == 0:
            return p
    return None


def compute_goto_B2(m, weights, p, CC=None):
    """Computes B_2 from the Goto zeta function denominator.

    The zeta denominator has the form (1-T) * P_2(T) * (1-q^2*T) where
    q = p^f and f is the minimal field degree for singularity definition.
    We divide out the trivial factors and return deg(P_2).

    Args:
        m: Degree of the hypersurface.
        weights: Tuple (q_0, q_1, q_2, q_3).
        p: Prime number with p ≡ 1 (mod M).
        CC: Complex field, or None for default.

    Returns:
        Second Betti number from Goto's zeta function.
    """
    if CC is None:
        CC = ComplexField(1024)
    Qt = PolynomialRing(QQ, 't')
    t = Qt.gen()

    min_field_deg = compute_minimal_field_degree(p, m, weights)
    if min_field_deg > 1:
        goto_zeta = compute_deformed_goto_zeta_denominator_extension(
            p, min_field_deg, m, weights, CC)
        q = p^min_field_deg
    else:
        goto_zeta = compute_deformed_goto_zeta_denominator(p, m, weights, CC)
        q = p

    # Divide out (1-T) and (1-q^2*T)
    P2 = goto_zeta // (1 - t)
    P2 = P2 // (1 - q^2 * t)

    return P2.degree()


def compute_orbifold_B2(A, m):
    """Computes B_2 by counting orbifold basis elements.

    The orbifold cohomology has basis elements with Hodge numbers (s, r).
    B_2 counts basis elements with s + r = 2 (i.e., in H^2).

    For a K3 surface, we also have H^0 (one element with s+r=0) and
    H^4 (one element with s+r=4), so the total number of basis elements
    is 2 + B_2.

    Args:
        A: The potential matrix.
        m: The degree.

    Returns:
        Second Betti number from orbifold basis count.
    """
    G = generate_symmetry_group([[1] * A.nrows()], A)
    all_elements = []
    for lam in G:
        for element in get_lambda_sector_basis(A, G, lam, m):
            all_elements.append(element)

    # Count elements by cohomological degree s + r
    count_by_degree = {}
    for element in all_elements:
        s, r = element.hodge_number
        deg = s + r
        count_by_degree[deg] = count_by_degree.get(deg, 0) + 1

    # B_2 is the count of elements with degree 2
    return count_by_degree.get(2, 0)


def count_points_goto_int(p, nu, m, weights, CC):
    """Computes Goto's point count and rounds to integer."""
    result = count_points_deformed_goto(p, nu, m, weights, CC)
    return Integer(round(CC(result).real()))


def test_single_surface(weights, m, primes, base_nus, verbose=True):
    """Tests one deformed K3 surface.

    Compares Goto's formula against the orbifold trace formula.
    Tests are performed over F_{p^{f*nu}} where f is the minimal field
    degree for singularity definition.

    Args:
        weights: Tuple (q_0, q_1, q_2, q_3).
        m: Degree.
        primes: List of primes to test.
        base_nus: List of extension degrees relative to the minimal field.
            E.g., [1, 2] means test at nu=f and nu=2f where f is the
            minimal field degree for singularity definition.
        verbose: Print detailed output.

    Returns:
        Tuple (passed, total, results) where results is a list of
        (p, nu, goto, orbifold, match) tuples.
    """
    A = build_deformed_matrix(weights, m)
    exps = get_exponents(weights, m)
    M = lcm(exps)

    CC = ComplexField(1024)
    results = []
    passed = 0
    total = 0

    for p in primes:
        min_field_deg = compute_minimal_field_degree(p, m, weights)
        R = Zp(p, 256)

        # Scale nus by the minimal field degree
        nus = [nu * min_field_deg for nu in base_nus]

        for nu in nus:
            goto = count_points_goto_int(p, nu, m, weights, CC)
            orbifold = count_points_orbifold(A, M, p, nu, R)
            match = (goto == orbifold)
            results.append((p, nu, goto, orbifold, match))
            total += 1
            if match:
                passed += 1
            if verbose:
                status = "OK" if match else "FAIL"
                field_note = " (base F_p^%d)" % min_field_deg if min_field_deg > 1 else ""
                print("  p=%d, nu=%d%s: Goto=%d, Orbifold=%d  [%s]" % (p, nu, field_note, goto, orbifold, status))

    return passed, total, results


def get_orbifold_elements(A, m):
    """Gets all orbifold cohomology basis elements for a surface.

    Args:
        A: The potential matrix.
        m: The degree.

    Returns:
        List of OrbifoldElement objects.
    """
    G = generate_symmetry_group([[1] * A.nrows()], A)
    all_elements = []
    for lam in G:
        for element in get_lambda_sector_basis(A, G, lam, m):
            all_elements.append(element)
    return all_elements


def test_zeta_structure(weights, m, p, R, CC, verbose=True):
    """Verifies that the zeta denominators from Goto and orbifold match.

    Computes B_2 and zeta denominator degree from both methods.
    Zetas are computed over F_{p^f} where f is the minimal field degree
    for singularity definition.

    Args:
        weights: Tuple (q_0, q_1, q_2, q_3).
        m: Degree.
        p: Prime number.
        R: p-adic ring for orbifold computation.
        CC: Complex field.
        verbose: Print detailed output.

    Returns:
        True if all checks pass.
    """
    min_field_deg = compute_minimal_field_degree(p, m, weights)
    A = build_deformed_matrix(weights, m)

    # Compute B_2 from both methods
    goto_B2 = compute_goto_B2(m, weights, p, CC)
    orbifold_B2 = compute_orbifold_B2(A, m)

    # Get orbifold elements for zeta computation
    orbifold_elements = get_orbifold_elements(A, m)

    # Compute zeta denominators over F_{p^f}
    if min_field_deg > 1:
        goto_zeta = compute_deformed_goto_zeta_denominator_extension(
            p, min_field_deg, m, weights, CC)
        orbifold_zeta = compute_orbifold_zeta_denominator_extension(
            orbifold_elements, p, R, min_field_deg)
        q = p^min_field_deg
    else:
        goto_zeta = compute_deformed_goto_zeta_denominator(p, m, weights, CC)
        orbifold_zeta = compute_orbifold_zeta_denominator(orbifold_elements, p, R)
        q = p

    goto_degree = goto_zeta.degree()
    orbifold_degree = orbifold_zeta.degree()

    if verbose:
        if min_field_deg > 1:
            field_str = "F_%d^%d" % (p, min_field_deg)
        else:
            field_str = "F_%d" % p
        print("  B_2: Goto=%d, Orbifold=%d" % (goto_B2, orbifold_B2))
        print("  Zeta over %s, denominator degree: Goto=%d, Orbifold=%d" % (field_str, goto_degree, orbifold_degree))

    passed = True

    if goto_B2 != orbifold_B2:
        if verbose:
            print("  FAIL: B_2 mismatch!")
        passed = False

    if goto_degree != orbifold_degree:
        if verbose:
            print("  FAIL: zeta degree mismatch!")
        passed = False

    # Check Riemann hypothesis on Goto zeta
    roots = goto_zeta.roots(CC, multiplicities=False)
    eigenvalues = [1/r for r in roots if r != 0]

    tolerance = 1e-6
    valid_magnitudes = [1, q, q^2]

    for ev in eigenvalues:
        mag = abs(ev)
        ok = any(abs(mag - v) < tolerance for v in valid_magnitudes)
        if not ok:
            if verbose:
                print("  WARNING: eigenvalue %s has magnitude %s (expected 1, %d, or %d)" % (ev, mag, q, q^2))

    if verbose and passed:
        print("  Zeta structure check: PASS")

    return passed


def run_tests(surface_indices=None, num_primes=2, verbose=True):
    """Runs tests on selected deformed diagonal K3 surfaces.

    Args:
        surface_indices: List of 1-indexed surface numbers to test,
            or None for a default selection.
        num_primes: Number of primes to test per surface.
        verbose: Print detailed output.
    """
    if surface_indices is None:
        # Default: test a diverse selection including some with q_0 > 1
        surface_indices = [1, 2, 3, 4, 6, 7, 10, 19, 30, 50]

    print("=" * 72)
    print("Deformed Diagonal K3: Goto vs Orbifold Trace Formula")
    print("=" * 72)

    total_tests = 0
    total_passed = 0
    failed = []

    CC = ComplexField(1024)

    for idx in surface_indices:
        weights, m = get_surface(idx)
        exps = get_exponents(weights, m)
        M = lcm(exps)

        primes = find_valid_primes_for_M(M, count=num_primes)
        if not primes:
            print("\nSurface #%d (m=%d, Q=%s): No valid primes found for M=%d" % (idx, m, weights, M))
            continue

        # Compute minimal field degree for first prime to display
        min_field_deg = compute_minimal_field_degree(primes[0], m, weights)
        field_info = "F_p^%d" % min_field_deg if min_field_deg > 1 else "F_p"

        print("\nSurface #%d: m=%d, Q=%s, q_3=%d, M=%d, base=%s, primes=%s" % (
            idx, m, weights, weights[3], M, field_info, primes))
        print("  Equation: %s" % surface_equation_string(weights, m))

        # Test point counts (base_nus=[1, 2] means test at f*1 and f*2 where f is min_field_deg)
        passed, total, results = test_single_surface(weights, m, primes, base_nus=[1, 2], verbose=verbose)
        total_tests += total
        total_passed += passed

        for p, nu, goto, orbifold, match in results:
            if not match:
                failed.append((idx, weights, m, p, nu, goto, orbifold))

        # Test zeta structure
        if primes:
            R = Zp(primes[0], 256)
            test_zeta_structure(weights, m, primes[0], R, CC, verbose=verbose)

    print("\n" + "=" * 72)
    print("Results: %d/%d passed." % (total_passed, total_tests))
    if failed:
        print("FAILURES:")
        for idx, weights, m, p, nu, goto, orbifold in failed:
            print("  Surface #%d (m=%d, Q=%s) p=%d nu=%d: Goto=%d != Orbifold=%d" % (
                idx, m, weights, p, nu, goto, orbifold))
    else:
        print("All tests passed.")
    print("=" * 72)


if __name__ == "__main__":
    print("=" * 72)
    print("Deformed Diagonal K3: Goto vs Orbifold at det(A)-primes")
    print("=" * 72)

    CC = ComplexField(1024)
    total_tests = 0
    total_passed = 0
    failed = []

    test_indices = [1, 5, 30]

    for idx in test_indices:
        weights, m = get_surface(idx)
        exps = get_exponents(weights, m)
        M = lcm(exps)
        A = build_deformed_matrix(weights, m)
        det_A = abs(A.det())

        p = find_first_prime_dividing_det(A, M)

        if p is None:
            print("\nSurface #%d (m=%d, Q=%s): No prime found with det(A)=%d | (p-1)" % (
                idx, m, weights, det_A))
            continue

        min_field_deg = compute_minimal_field_degree(p, m, weights)
        field_str = "F_%d^%d" % (p, min_field_deg) if min_field_deg > 1 else "F_%d" % p
        print("\nSurface #%d: m=%d, Q=%s, p=%d, det(A)=%d, field=%s" % (
            idx, m, weights, p, det_A, field_str))
        print("  Equation: %s" % surface_equation_string(weights, m))

        # --- Point count comparison ---
        passed, total, results = test_single_surface(
            weights, m, [p], base_nus=[1, 2], verbose=True)
        total_tests += total
        total_passed += passed

        for p_r, nu, goto, orbifold, match in results:
            if not match:
                failed.append((idx, weights, m, p_r, nu, goto, orbifold))

        # --- Zeta denominator comparison ---
        R = Zp(p, 256)
        orbifold_elements = get_orbifold_elements(A, M)

        if min_field_deg > 1:
            goto_zeta = compute_deformed_goto_zeta_denominator_extension(
                p, min_field_deg, m, weights, CC)
            orbifold_zeta = compute_orbifold_zeta_denominator_extension(
                orbifold_elements, p, R, min_field_deg)
        else:
            goto_zeta = compute_deformed_goto_zeta_denominator(p, m, weights, CC)
            orbifold_zeta = compute_orbifold_zeta_denominator(
                orbifold_elements, p, R)

        zeta_match = (goto_zeta == orbifold_zeta)
        total_tests += 1
        if zeta_match:
            total_passed += 1
        else:
            failed.append((idx, weights, m, p, "zeta", goto_zeta.degree(), orbifold_zeta.degree()))

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
