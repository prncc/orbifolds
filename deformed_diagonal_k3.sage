"""Goto's formula for weighted deformed diagonal K3 surfaces.

Computes point counts and zeta functions via Jacobi sums and cyclic
quotient singularity resolution, following Theorem 4.3.4 and
Corollary 4.3.7 from Goto's thesis.

A deformed diagonal surface has equation:
    c_0 x_0^{m_0} + c_1 x_1^{m_1} + c_2 x_2^{m_2} + c_3 x_0 x_3^{m_3} = 0

where m_i = m/q_i for i=0,1,2 and m_3 = (m - q_0)/q_3.

Reference: Goto, "Arithmetic of weighted diagonal surfaces and weighted
deformed diagonal surfaces over finite fields" (1994), Chapter 4. Ph.D. Thesis.
"""

load("diagonal_k3.sage")


def calculate_deformed_exponents(m, weights):
    """Computes the exponents m_i for a deformed diagonal surface.

    For a deformed diagonal surface:
        m_i = m/q_i for i = 0, 1, 2
        m_3 = (m - q_0)/q_3

    Args:
        m: Degree of the hypersurface.
        weights: Tuple (q_0, q_1, q_2, q_3) of positive integer weights.

    Returns:
        Tuple (m_0, m_1, m_2, m_3) of exponents.
    """
    q_0, q_1, q_2, q_3 = weights
    m_0 = m // q_0
    m_1 = m // q_1
    m_2 = m // q_2
    m_3 = (m - q_0) // q_3
    return (m_0, m_1, m_2, m_3)


def iterate_deformed_coordinate_pairs():
    """Yields coordinate pairs in I_1 = I \ {(1,3), (2,3)}.

    For deformed diagonal surfaces, the sets P_{13} and P_{23} are empty
    (cf. Goto eq 4.3), so we only consider pairs (0,1), (0,2), (0,3), (1,2).

    Returns:
        Generator of (i, j) tuples with 0 <= i < j <= 3.
    """
    excluded = {(1, 3), (2, 3)}
    for i in range(4):
        for j in range(i + 1, 4):
            if (i, j) not in excluded:
                yield (i, j)


def calculate_alpha_3(weights):
    """Computes alpha_3 for the isolated singularity at [0,0,0,1].

    Finds the unique alpha_3 in [1, q_3) satisfying:
        q_1 * alpha_3 ≡ q_2 (mod q_3)

    Per Goto Proposition 4.2.7.

    Args:
        weights: Tuple (q_0, q_1, q_2, q_3) of positive integer weights.

    Returns:
        Positive integer alpha_3 in [1, q_3), or None if q_3 = 1.

    Raises:
        ValueError: If no valid alpha_3 is found (should not happen).
    """
    q_0, q_1, q_2, q_3 = weights
    if q_3 == 1:
        return None
    for alpha_3 in range(1, q_3):
        if (q_1 * alpha_3 - q_2) % q_3 == 0:
            return alpha_3
    raise ValueError(f"Unable to calculate alpha_3 for weights {weights}.")


def calculate_r_3(q_3, alpha_3):
    """Computes r_3, the length of the continued fraction of q_3/alpha_3.

    This equals the number of exceptional divisors in the minimal
    resolution of the cyclic quotient singularity A_{q_3, alpha_3}.

    Args:
        q_3: The weight q_3.
        alpha_3: The alpha_3 parameter from calculate_alpha_3.

    Returns:
        Non-negative integer r_3. Returns 0 if q_3 = 1.
    """
    if q_3 == 1 or alpha_3 is None:
        return 0
    return len(get_resolution_continued_fraction(q_3 / alpha_3))


def _has_solution_gamma_pow_eq_minus_one(p, f, m_eff):
    """Check if gamma^{m_eff} = -1 has a solution in F_{p^f}.

    The equation gamma^m = -1 has solutions in F_{p^f}^* iff -1 is an m-th power.
    Since F_{p^f}^* is cyclic of order p^f - 1, -1 is an m-th power iff
    (-1)^{(p^f - 1) / gcd(m, p^f - 1)} = 1, i.e., (p^f - 1) / gcd(m, p^f - 1) is even.
    """
    order = p^f - 1
    g = gcd(m_eff, order)
    exponent = order // g
    return exponent % 2 == 0


def compute_minimal_field_degree(p, m, weights):
    """Computes the minimal extension degree f such that all singularities are defined over F_{p^f}.

    For singularities in P_ij, coordinates satisfy gamma^{m_eff} = -1.
    We find the smallest f such that this equation has solutions in F_{p^f} for all
    singular pairs (i, j).

    The isolated singularity [0,0,0,1] is always defined over F_p.

    Args:
        p: Prime number.
        m: Degree of the hypersurface.
        weights: Tuple (q_0, q_1, q_2, q_3) of positive integer weights.

    Returns:
        Positive integer f such that all P_ij singularities are defined over F_{p^f}.
    """
    # Collect all m_eff values we need to satisfy
    m_eff_values = []

    for i, j in iterate_deformed_coordinate_pairs():
        d_ij = calculate_goto_d_ij(weights, i, j)
        if d_ij < 2:
            continue

        # For (0,3), use m - q_0 instead of m (Proposition 4.2.4)
        if (i, j) == (0, 3):
            m_eff = m - weights[0]
        else:
            m_eff = m

        m_eff_values.append(m_eff)

    if not m_eff_values:
        return Integer(1)  # No singular pairs, defined over F_p

    # Find minimal f such that gamma^{m_eff} = -1 has solutions for all m_eff
    for f in range(1, 101):
        all_have_solutions = True
        for m_eff in m_eff_values:
            if not _has_solution_gamma_pow_eq_minus_one(p, f, m_eff):
                all_have_solutions = False
                break
        if all_have_solutions:
            return Integer(f)

    raise ValueError("Could not find suitable field extension degree")


def generate_deformed_V_terms(m, weights):
    """Generates index pairs (a_1, a_2) in the set V.

    The set V consists of pairs (a_1, a_2) such that:
        - a_i ∈ M_i * Z/MZ, a_i ≠ 0 for i = 1, 2
        - a_1 + a_2 = 0

    where M = lcm(m_0, m_1, m_2, m_3) and M_i = M/m_i.

    Per Goto Proposition 4.3.2.

    Args:
        m: Degree of the hypersurface.
        weights: Tuple (q_0, q_1, q_2, q_3) of positive integer weights.

    Returns:
        List of (a_1, a_2) tuples.
    """
    exponents = calculate_deformed_exponents(m, weights)
    m_0, m_1, m_2, m_3 = exponents
    M = lcm([m_0, m_1, m_2, m_3])
    M_1 = M // m_1
    M_2 = M // m_2

    results = []
    # a_1 is a nonzero multiple of M_1 in [1, M-1]
    for k_1 in range(1, m_1):  # k_1 = 1, ..., m_1 - 1
        a_1 = k_1 * M_1
        # a_2 = -a_1 mod M = M - a_1
        a_2 = M - a_1
        # Check a_2 is a nonzero multiple of M_2
        if a_2 > 0 and a_2 % M_2 == 0:
            results.append((a_1, a_2))
    return results


def generate_deformed_W_terms(m, weights):
    """Generates index tuples (a_0, a_1, a_2, a_3) in the set W.

    The set W consists of tuples (a_0, a_1, a_2, a_3) such that:
        - a_0 ∈ Z/MZ, a_0 ≠ 0
        - a_i ∈ M_i * Z/MZ, a_i ≠ 0 for i = 1, 2, 3
        - sum(a_i) = 0 mod M
        - m_0 * a_0 + a_3 = 0 mod M

    where M = lcm(m_0, m_1, m_2, m_3) and M_i = M/m_i.

    Per Goto Proposition 4.3.2.

    Args:
        m: Degree of the hypersurface.
        weights: Tuple (q_0, q_1, q_2, q_3) of positive integer weights.

    Returns:
        List of (a_0, a_1, a_2, a_3) tuples.
    """
    exponents = calculate_deformed_exponents(m, weights)
    m_0, m_1, m_2, m_3 = exponents
    M = lcm([m_0, m_1, m_2, m_3])
    M_1 = M // m_1
    M_2 = M // m_2
    M_3 = M // m_3

    results = []

    # a_0 ranges over all nonzero elements of Z/MZ
    for a_0 in range(1, M):
        # a_3 is determined by: m_0 * a_0 + a_3 ≡ 0 (mod M)
        a_3 = (-m_0 * a_0) % M
        # Check a_3 is nonzero and a multiple of M_3
        if a_3 == 0 or a_3 % M_3 != 0:
            continue

        # a_1 is a nonzero multiple of M_1
        for k_1 in range(1, m_1):
            a_1 = k_1 * M_1

            # a_2 is determined by: a_0 + a_1 + a_2 + a_3 ≡ 0 (mod M)
            a_2 = (-(a_0 + a_1 + a_3)) % M
            # Check a_2 is nonzero and a multiple of M_2
            if a_2 == 0 or a_2 % M_2 != 0:
                continue

            results.append((a_0, a_1, a_2, a_3))

    return results


def calculate_direct_fij(p, f, m, e_ij):
    """Computes f_ij directly as the multiplicative order of (-1)^{(p^f-1)/N_ij}.

    For non-(0,3) pairs, f_ij is determined by evaluating (-1) raised to
    the power (p^f - 1) / N_ij in GF(p^f) and taking its multiplicative order,
    where N_ij = m / e_ij.

    Args:
        p: Prime number.
        f: Extension degree (singularities defined over F_{p^f}).
        m: Degree of the hypersurface.
        e_ij: lcm(q_i, q_j) for the coordinate pair.

    Returns:
        Positive integer giving the multiplicative order.
    """
    field_size = p^f
    F = GF(field_size)
    N_ij = m // e_ij
    exponent = (field_size - 1) // N_ij
    return (F(-1)^exponent).multiplicative_order()


class DeformedK3SingularityTerm(SageObject):
    """A singularity correction term for a coordinate pair (i, j) in I_1.

    Similar to K3SingularityTerm but only for pairs in
    I_1 = I \ {(1,3), (2,3)}.

    Attributes:
        m: Degree of the hypersurface.
        weights: Tuple of 4 weights.
        ij: Tuple (i, j) of coordinate indices.
        m_ij: 1 or 2, depending on -1 being an m-th power in F_p^*.
        d_ij: gcd(q_i, q_j).
        e_ij: lcm(q_i, q_j).
        f_ij: m_ij / gcd(m_ij, e_ij).
        omega_ij: m / (e_ij * f_ij).
        alpha_ij: Singularity type parameter.
        eta_ij: Primitive f_ij-th root of unity.
        r_ij: Length of continued fraction of d_ij / alpha_ij.
    """

    def __init__(self, p, m, weights, i, j, CC=None):
        """Initializes a singularity term for coordinate pair (i, j).

        Args:
            p: Prime number.
            m: Degree of the hypersurface.
            weights: Tuple of 4 positive integer weights.
            i: First coordinate index.
            j: Second coordinate index.
            CC: Complex field, or None for default 1024-bit precision.
        """
        self.CC = CC or ComplexField(1024)
        self.m = m
        self.weights = weights
        self.ij = (i, j)

        # For (0,3), use (m - q_0) instead of m (Proposition 4.2.8)
        # This affects m_ij, f_ij, and omega_ij
        if (i, j) == (0, 3):
            m_effective = m - weights[0]
        else:
            m_effective = m

        self.m_ij = get_goto_m_ij(p, m_effective)
        self.d_ij = calculate_goto_d_ij(weights, i, j)
        self.e_ij = calculate_goto_e_ij(weights, i, j)
        if (i, j) == (0, 3):
            self.f_ij = calculate_direct_fij(p, 1, m_effective, self.e_ij)
        else:
            self.f_ij = calculate_goto_f_ij(self.m_ij, self.e_ij)
        self.omega_ij = calculate_goto_omega_ij(m_effective, self.e_ij, self.f_ij)

        self.alpha_ij = calculate_goto_alpha_ij(self.d_ij, weights, i, j)
        self.eta_ij = calculate_goto_eta_ij(self.f_ij, self.CC)
        self.r_ij = calculate_goto_r_ij(self.d_ij, self.alpha_ij)

    def _latex_(self):
        """Returns LaTeX representation for Sage display."""
        return f"A_{{{self.d_ij}, {self.alpha_ij}}} ({self.alpha_ij}, {self.f_ij})"

    def _repr_(self):
        """Returns plain text representation for Sage display."""
        return f"A_({self.d_ij}, {self.alpha_ij}) ({self.omega_ij})"

    def get_zeta_contribution(self, p, nu, var='u'):
        """Computes the singularity contribution to the zeta function denominator.

        Returns prod_{k=0}^{f-1} (1 - eta^k * p^nu * var)^{r * omega}.

        Args:
            p: Prime number.
            nu: Extension degree.
            var: Variable name for the polynomial ring.

        Returns:
            Polynomial in CC[var].
        """
        term = 1
        CRng = PolynomialRing(self.CC, var)
        u = CRng.gen()
        for f_ij_power in range(self.f_ij):
            term *= (1 - (self.eta_ij^f_ij_power * p)^nu * u)
        return term^(self.r_ij * self.omega_ij)

    def get_point_count(self, p, nu):
        """Computes the singularity correction to the point count.

        Returns omega * r * sum_{k=0}^{f-1} (eta^k * p)^nu if f | nu,
        else 0.

        Args:
            p: Prime number.
            nu: Extension degree.

        Returns:
            Complex number (real-valued) giving the correction.
        """
        if not Integer(self.f_ij).divides(nu):
            return 0
        total = 0
        for a in range(0, self.f_ij):
            total += (self.eta_ij^a * p)^nu
        return total * self.omega_ij * self.r_ij


class IsolatedSingularityTerm(SageObject):
    """The isolated singularity at [0,0,0,1] for deformed diagonal surfaces.

    When q_3 >= 2, the point [0,0,0,1] is a cyclic quotient singularity
    of type A_{q_3, alpha_3}. This contributes (1-qT)^{r_3} to P_2.

    Attributes:
        q_3: The weight q_3.
        alpha_3: Singularity type parameter.
        r_3: Length of continued fraction of q_3/alpha_3.
    """

    def __init__(self, weights, CC=None):
        """Initializes the isolated singularity term.

        Args:
            weights: Tuple of 4 positive integer weights.
            CC: Complex field, or None for default 1024-bit precision.
        """
        self.CC = CC or ComplexField(1024)
        self.weights = weights
        self.q_3 = weights[3]
        self.alpha_3 = calculate_alpha_3(weights)
        self.r_3 = calculate_r_3(self.q_3, self.alpha_3)

    def _latex_(self):
        """Returns LaTeX representation for Sage display."""
        if self.q_3 == 1:
            return "\\text{(no isolated singularity)}"
        return f"A_{{{self.q_3}, {self.alpha_3}}} \\text{{ at }} [0,0,0,1]"

    def _repr_(self):
        """Returns plain text representation for Sage display."""
        if self.q_3 == 1:
            return "(no isolated singularity)"
        return f"A_({self.q_3}, {self.alpha_3}) at [0,0,0,1]"

    def get_zeta_contribution(self, p, nu, var='u'):
        """Computes the contribution to the zeta function denominator.

        Returns (1 - p^nu * var)^{r_3}.

        Args:
            p: Prime number.
            nu: Extension degree.
            var: Variable name for the polynomial ring.

        Returns:
            Polynomial in CC[var].
        """
        CRng = PolynomialRing(self.CC, var)
        u = CRng.gen()
        return (1 - p^nu * u)^self.r_3

    def get_point_count(self, p, nu):
        """Computes the point count contribution from the isolated singularity.

        Returns r_3 * p^nu.

        Args:
            p: Prime number.
            nu: Extension degree.

        Returns:
            Integer giving the correction.
        """
        return self.r_3 * p^nu


def collect_all_singularity_terms(p, m, weights, CC):
    """Collects all singularity correction terms for a deformed K3 surface.

    Per Proposition 4.2.3, the singular locus is:
        X_sing = union_{(i,j) in I_1} P_ij                  if q_3 = 1
        X_sing = {[0,0,0,1]} union_{(i,j) in I_1} P_ij      if q_3 >= 2

    where I_1 = {(i,j) in I \ {(1,3), (2,3)} | d_ij >= 2}.

    Args:
        p: Prime number.
        m: Degree of the hypersurface.
        weights: Tuple of 4 positive integer weights.
        CC: Complex field for computation.

    Returns:
        Tuple (isolated_term, curve_terms) where:
            isolated_term: IsolatedSingularityTerm or None (for [0,0,0,1])
            curve_terms: List of DeformedK3SingularityTerm objects (for P_ij)
    """
    q_3 = weights[3]

    # Isolated singularity at [0,0,0,1] when q_3 >= 2
    if q_3 >= 2:
        isolated_term = IsolatedSingularityTerm(weights, CC)
    else:
        isolated_term = None

    # Singularities from I_1 = {(i,j) in I \ {(1,3), (2,3)} | d_ij >= 2}
    curve_terms = []
    for i, j in iterate_deformed_coordinate_pairs():
        d_ij = calculate_goto_d_ij(weights, i, j)
        if d_ij <= 1:
            continue
        curve_terms.append(DeformedK3SingularityTerm(p, m, weights, i, j, CC))

    return isolated_term, curve_terms


def count_points_deformed_goto(p, nu, m, weights, CC=None):
    """Counts F_{p^nu}-points on the minimal resolution of a deformed diagonal K3.

    Uses Goto's Theorem 4.3.4 and Corollary 4.3.7 for the surface
        x_0^{m_0} + x_1^{m_1} + x_2^{m_2} + x_0 x_3^{m_3} = 0
    in P^3(Q) with twist c = (1,1,1,1).

    Computes:
        N_nu(X~) = 1 + q^{2nu}
                   + (1 + r_3) * q^nu
                   + sum_{(a_1,a_2) in V} q^nu * chi^{-nu}((-1)^{a_1})
                   + sum_{a in W} j(a)^nu
                   + sum_{(i,j) in I_1} singularity corrections

    where q = p.

    Args:
        p: Prime number with p ≡ 1 (mod M).
        nu: Extension degree (count points over F_{p^nu}).
        m: Degree of the hypersurface.
        weights: Tuple of 4 positive integer weights.
        CC: Complex field, or None for default 1024-bit precision.

    Returns:
        Complex number (should be a real integer) giving the point count.
    """
    if CC is None:
        CC = ComplexField(1024)

    K = GF(p)
    q = p

    # Base: 1 + q^{2nu}
    total = 1 + q^(2*nu)

    # The factor (1-qT) in P_2(X_k, T) contributes q^nu to the point count
    total += q^nu

    # V terms: sum over (a_1, a_2) in V of q^nu * chi^{-nu}((-1)^{a_1})
    # In the untwisted case, chi^{-1}((-1)^{a_1}) is:
    #   1 if a_1 is even or -1 is an M-th power in F_p^*
    #   -1 if a_1 is odd and -1 is not an M-th power
    exponents = calculate_deformed_exponents(m, weights)
    M = lcm(exponents)

    # Check if -1 is an M-th power in F_p^*
    minus_one_is_Mth_power = any(K(x^M) == K(-1) for x in K.unit_group())

    v_terms = generate_deformed_V_terms(m, weights)
    for a_1, a_2 in v_terms:
        if minus_one_is_Mth_power:
            chi_factor = 1
        else:
            # chi^{-1}((-1)^{a_1}) = (-1)^{a_1} when chi(-1) = -1
            chi_factor = 1 if a_1 % 2 == 0 else -1
        total += q^nu * chi_factor^nu

    # W terms: sum over a in W of j(a)^nu
    w_terms = generate_deformed_W_terms(m, weights)
    for a in w_terms:
        character_fractions = [a_i / M for a_i in a]
        js = calculate_jacobi_sum(K, character_fractions)
        total += js^nu

    # Singularity corrections from resolution (Proposition 4.2.3)
    isolated_term, curve_terms = collect_all_singularity_terms(p, m, weights, CC)

    # Isolated singularity at [0,0,0,1] when q_3 >= 2
    if isolated_term is not None:
        total += isolated_term.get_point_count(p, nu)

    # Curve singularities from I_1
    for sing in curve_terms:
        total += sing.get_point_count(p, nu)

    return total


def compute_deformed_goto_zeta_denominator(p, m, weights, CC=None):
    """Computes the zeta function denominator using Goto's Theorem 4.3.4.

    Assembles the denominator as (1-t)(1-q^2 t) * P_2(X_k, t) where P_2
    includes:
        - (1-qt)^{1+r_3}
        - prod over V: (1 - q * chi^{-1}((-1)^a_1) * t)
        - prod over W: (1 - j(a) * t)
        - singularity terms from I_1

    Args:
        p: Prime number with p ≡ 1 (mod M).
        m: Degree of the hypersurface.
        weights: Tuple of 4 positive integer weights.
        CC: Complex field, or None for default 1024-bit precision.

    Returns:
        Polynomial in QQ['t'] giving the zeta function denominator.
    """
    if CC is None:
        CC = ComplexField(1024)

    Ct = PolynomialRing(CC, 'u')
    u = Ct.gen()

    K = GF(p)

    # H^0 and H^4: (1 - t)(1 - p^2 t)
    denom = (1 - u) * (1 - p^2 * u)

    # The (1-qt) factor from P_2
    denom *= (1 - p * u)

    # V terms
    exponents = calculate_deformed_exponents(m, weights)
    M = lcm(exponents)

    minus_one_is_Mth_power = any(K(x^M) == K(-1) for x in K.unit_group())

    v_terms = generate_deformed_V_terms(m, weights)
    for a_1, a_2 in v_terms:
        if minus_one_is_Mth_power:
            chi_factor = 1
        else:
            # chi^{-1}((-1)^{a_1}) = (-1)^{a_1} when chi(-1) = -1
            chi_factor = 1 if a_1 % 2 == 0 else -1
        eigenvalue = p * chi_factor
        denom *= (1 - eigenvalue * u)

    # W terms: prod over a in W of (1 - j(a) * t)
    w_terms = generate_deformed_W_terms(m, weights)
    for a in w_terms:
        character_fractions = [a_i / M for a_i in a]
        js = CC(calculate_jacobi_sum(K, character_fractions))
        denom *= (1 - js * u)

    # Singularity corrections (Proposition 4.2.3)
    isolated_term, curve_terms = collect_all_singularity_terms(p, m, weights, CC)

    # Isolated singularity at [0,0,0,1]: (1 - pt)^{r_3}
    if isolated_term is not None:
        denom *= isolated_term.get_zeta_contribution(p, 1)

    # Curve singularities from I_1
    for sing in curve_terms:
        for k in range(sing.f_ij):
            denom *= (1 - sing.eta_ij^k * p * u) ^ (sing.r_ij * sing.omega_ij)

    # Round complex coefficients to integers
    Qt = PolynomialRing(QQ, 't')
    t = Qt.gen()
    result = Qt(0)
    for deg, coef in enumerate(list(denom)):
        result += Integer(coef.real().round()) * t^deg

    return result


def compute_deformed_goto_zeta_denominator_extension(p, f, m, weights, CC=None):
    """Computes the zeta function denominator over F_{p^f} using Goto's formula.

    Base change raises all Frobenius eigenvalues to the f-th power.

    Args:
        p: Prime number with p ≡ 1 (mod M).
        f: Extension degree (positive integer).
        m: Degree of the hypersurface.
        weights: Tuple of 4 positive integer weights.
        CC: Complex field, or None for default 1024-bit precision.

    Returns:
        Polynomial in QQ['t'] giving the zeta denominator over F_{p^f}.
    """
    if CC is None:
        CC = ComplexField(1024)

    Ct = PolynomialRing(CC, 'u')
    u = Ct.gen()

    K = GF(p)
    q = p^f

    # H^0 and H^4: (1 - T)(1 - q^2 T)
    denom = (1 - u) * (1 - q^2 * u)

    # The (1-qT) factor from P_2
    denom *= (1 - q * u)

    # V terms with eigenvalues raised to f-th power
    exponents = calculate_deformed_exponents(m, weights)
    M = lcm(exponents)

    minus_one_is_Mth_power = any(K(x^M) == K(-1) for x in K.unit_group())

    v_terms = generate_deformed_V_terms(m, weights)
    for a_1, a_2 in v_terms:
        if minus_one_is_Mth_power:
            chi_factor = 1
        else:
            # chi^{-1}((-1)^{a_1}) = (-1)^{a_1} when chi(-1) = -1
            chi_factor = 1 if a_1 % 2 == 0 else -1
        eigenvalue = (p * chi_factor)^f
        denom *= (1 - eigenvalue * u)

    # W terms with Jacobi sums raised to f-th power
    w_terms = generate_deformed_W_terms(m, weights)
    for a in w_terms:
        character_fractions = [a_i / M for a_i in a]
        js = CC(calculate_jacobi_sum(K, character_fractions))
        denom *= (1 - js^f * u)

    # Singularity corrections (Proposition 4.2.3)
    isolated_term, curve_terms = collect_all_singularity_terms(p, m, weights, CC)

    # Isolated singularity at [0,0,0,1]: (1 - q T)^{r_3}
    if isolated_term is not None:
        denom *= isolated_term.get_zeta_contribution(p, f)

    # Curve singularities from I_1: eigenvalues (eta^k * p)^f
    for sing in curve_terms:
        for k in range(sing.f_ij):
            eigenvalue = (sing.eta_ij^k * p)^f
            denom *= (1 - eigenvalue * u) ^ (sing.r_ij * sing.omega_ij)

    # Round complex coefficients to integers
    Qt = PolynomialRing(QQ, 't')
    t = Qt.gen()
    result = Qt(0)
    for deg, coef in enumerate(list(denom)):
        result += Integer(coef.real().round()) * t^deg

    return result
