"""Goto's formula for weighted diagonal K3 surfaces.

Computes point counts and zeta functions via Jacobi sums and cyclic
quotient singularity resolution, following Goto's Theorem 5.2.
"""

import itertools


def generate_diagonal_k3_jacobi_terms(m, weights):
    """Generates Jacobi sum argument tuples for a diagonal K3 surface.

    Enumerates all tuples a = (a_0, ..., a_3) with each a_i a nonzero
    multiple of weights[i], such that sum(a_i) is divisible by m.
    These tuples parameterize the Jacobi sums appearing in the smooth
    part of the point-count formula.

    Args:
        m: The degree of the weighted projective hypersurface.
        weights: Tuple of 4 positive integer weights (q_0, q_1, q_2, q_3).

    Returns:
        List of 4-tuples of integers, each a valid Jacobi sum argument.
    """
    a_slots = [[] for _ in range(len(weights))]
    for i in range(len(weights)):
        for s_i in range(m / weights[i]):
            if s_i == 0:
                continue
            a_slots[i].append(s_i * weights[i])
    results = []
    for a in itertools.product(*a_slots):
        a_sum = sum([a_i for a_i in a])
        if a_sum % m == 0:
            results.append(a)
    return results


@cached_function
def calculate_gauss_sum(K, a_fraction):
    """Computes a Gauss sum over the finite field K.

    Evaluates g(chi_a) = sum_{t in K^*} psi(w^t) chi_a(w^t), where
    w is a multiplicative generator of K^*, psi is the standard
    additive character, and chi_a is determined by a_fraction.

    Args:
        K: A finite field GF(p).
        a_fraction: Rational number determining the multiplicative
            character chi_a via chi_a(w) = e^(2 pi i * a_fraction).

    Returns:
        Complex number giving the Gauss sum.
    """
    p = K.characteristic()
    w = K.multiplicative_generator()

    def _psi(t):
        """Additive character (K, +) -> CC.

        For x \in K, we'd use e^(2 * pi * I * Integer(x) / p),
        but since we're summing over K^* jointly with _chi,
        we pass in t = Integer(x) directly and encounter the same
        terms over the sum.
        """
        return e^(2 * pi * I * Integer(w^t) / p)

    def _chi(t):
        """Multiplicative character determined by a_fraction \in Q \ Z.

        If K^* = <w> and x = x^t, this is the character that sends
        w -> e^(2 * pi * I * a_fraction) and x -> e^(2 * pi * I * t * a_fraction).
        However, since we're summing over K^* jointly with _psi,
        we pass t in directly and encounter the same terms over the sum.
        """
        return e^(2 * pi * I * t * a_fraction)

    total = 0
    for t in range(p - 1):
        term = _psi(t) * _chi(t)
        total += term
    return total


def calculate_jacobi_sum(K, character_fractions):
    """Computes a Jacobi sum as a product of Gauss sums.

    Evaluates J(chi_{a_1}, ..., chi_{a_n}) = (1/p) * prod g(chi_{a_i}).

    Args:
        K: A finite field GF(p).
        character_fractions: List of rational numbers determining the
            multiplicative characters.

    Returns:
        Complex number giving the Jacobi sum.

    Raises:
        NotImplementedError: If K is not a prime field (degree > 1).
    """
    p = K.characteristic()
    if K.degree() > 1:
        raise NotImplementedError("Function assumes K = F_p.")
    total = 1 / p
    for a_fraction in character_fractions:
        total *= calculate_gauss_sum(K, a_fraction)
    return total


def iterate_coordinate_pairs():
    """Yields all 6 coordinate pairs (i, j) with 0 <= i < j <= 3.

    Returns:
        Generator of (i, j) tuples.
    """
    for i in range(4):
        for j in range(i + 1, 4):
            yield (i, j)


def _has_solution_gamma_pow_eq_minus_one(p, f, m):
    """Check if gamma^m = -1 has a solution in F_{p^f}.

    The equation gamma^m = -1 has solutions in F_{p^f}^* iff -1 is an m-th power.
    Since F_{p^f}^* is cyclic of order p^f - 1, -1 is an m-th power iff
    (-1)^{(p^f - 1) / gcd(m, p^f - 1)} = 1, i.e., (p^f - 1) / gcd(m, p^f - 1) is even.

    Args:
        p: Prime number.
        f: Extension degree.
        m: Exponent in the equation gamma^m = -1.

    Returns:
        True if a solution exists, False otherwise.
    """
    order = p^f - 1
    g = gcd(m, order)
    exponent = order // g
    return exponent % 2 == 0


def compute_minimal_field_degree_diagonal(p, m, weights):
    """Computes the minimal extension degree f for diagonal surface singularities.

    For the untwisted case (c = (1,1,1,1)), computes the smallest f such that
    all singularities P_ij are defined over F_{p^f}.

    By Proposition 4.6 and Lemma 4.5, for untwisted diagonal surfaces:
      - If -1 is an m-th power in F_p^*: f_ij = 1 for all (i,j)
      - Otherwise: f_ij = 2 if e_ij is odd, else f_ij = 1

    We return max_{(i,j)} f_ij, which equals 1 or 2.

    Args:
        p: Prime number.
        m: Degree of the hypersurface.
        weights: Tuple (q_0, q_1, q_2, q_3) of positive integer weights.

    Returns:
        Positive integer f such that all P_ij singularities are defined over F_{p^f}.
    """
    # Check if there are any singular pairs
    has_singularities = False
    for i, j in iterate_coordinate_pairs():
        d_ij = calculate_goto_d_ij(weights, i, j)
        if d_ij >= 2:
            has_singularities = True
            break

    if not has_singularities:
        return Integer(1)  # No singular pairs, defined over F_p

    # For diagonal surfaces, all singularities use m_eff = m
    # Find minimal f such that gamma^m = -1 has solutions in F_{p^f}
    for f in range(1, 101):
        if _has_solution_gamma_pow_eq_minus_one(p, f, m):
            return Integer(f)

    raise ValueError("Could not find suitable field extension degree")


def get_complement(i, j):
    """Returns the two coordinate indices not equal to i or j.

    Args:
        i: First coordinate index (0-3).
        j: Second coordinate index (0-3), with j != i.

    Returns:
        Tuple of the two remaining indices from {0, 1, 2, 3}.
    """
    return tuple(
        index for index in range(4)
        if (index != i and index != j)
    )


def get_resolution_continued_fraction(fraction):
    """Computes the Hirzebruch-Jung continued fraction expansion.

    Given a rational number > 1, computes its continued fraction
    expansion [b_1, b_2, ...] where each b_i = ceil of the current
    value. This parameterizes the minimal resolution of a cyclic
    quotient singularity.

    Args:
        fraction: Rational number greater than 1.

    Returns:
        List of integers giving the continued fraction expansion.
    """
    continued_fraction = []
    while fraction > 1:
        b = ceil(fraction)
        continued_fraction.append(b)
        remaining = b - fraction
        if remaining == 0:
            break
        fraction = 1/remaining
    return continued_fraction


def calculate_goto_d_ij(weights, i, j):
    """Computes d_ij = gcd(q_i, q_j) for a coordinate pair.

    Args:
        weights: Tuple of 4 positive integer weights.
        i: First coordinate index.
        j: Second coordinate index.

    Returns:
        The gcd of weights[i] and weights[j].
    """
    return gcd(weights[i], weights[j])


def calculate_goto_e_ij(weights, i, j):
    """Computes e_ij = lcm(q_i, q_j) for a coordinate pair.

    Args:
        weights: Tuple of 4 positive integer weights.
        i: First coordinate index.
        j: Second coordinate index.

    Returns:
        The lcm of weights[i] and weights[j].
    """
    return lcm(weights[i], weights[j])


def calculate_goto_f_ij(m_ij, e_ij):
    """Computes f_ij = m_ij / gcd(m_ij, e_ij).

    This determines the order of the eta root of unity appearing in
    the singularity correction.

    Args:
        m_ij: The value m_ij (1 or 2, depending on whether -1 is an
            m-th power in F_p^*).
        e_ij: The lcm of the two weights.

    Returns:
        Positive integer f_ij.
    """
    return m_ij / gcd(m_ij, e_ij)


def get_goto_m_ij(p, m):
    """Computes m_ij: 1 if -1 is an m-th power in F_p^*, else 2.

    In the untwisted case (c = (1,1,1,1)), m_ij depends only on
    whether -1 lies in the subgroup of m-th powers.

    Args:
        p: Prime number.
        m: Degree of the hypersurface.

    Returns:
        1 or 2.
    """
    K = GF(p)
    for x in K.unit_group():
        if K(-1) == K(x^m):
            return 1
    return 2


def calculate_goto_r_ij(d_ij, alpha_ij):
    """Computes r_ij, the length of the continued fraction of d_ij / alpha_ij.

    This equals the number of exceptional divisors in the minimal
    resolution of the cyclic quotient singularity A_{d_ij, alpha_ij}.

    Args:
        d_ij: The gcd of the two weights.
        alpha_ij: The Goto alpha parameter for this pair.

    Returns:
        Non-negative integer r_ij.
    """
    return len(get_resolution_continued_fraction(d_ij / alpha_ij))


def calculate_goto_eta_ij(f_ij, CC=None):
    """Computes eta_ij, a primitive f_ij-th root of unity.

    Args:
        f_ij: Positive integer order of the root.
        CC: Complex field for computation, or None for default
            1024-bit precision.

    Returns:
        A complex f_ij-th root of unity.
    """
    if CC is None:
        CC = ComplexField(1024)
    return cyclotomic_polynomial(f_ij).roots(CC, multiplicities=False)[0]


def calculate_goto_omega_ij(m, e_ij, f_ij):
    """Computes omega_ij = m / (e_ij * f_ij).

    This multiplicity factor counts how many singularity points on
    the surface correspond to a given coordinate pair.

    Args:
        m: Degree of the hypersurface.
        e_ij: The lcm of weights[i] and weights[j].
        f_ij: The f_ij parameter.

    Returns:
        Positive integer omega_ij.
    """
    return m / (e_ij * f_ij)


def calculate_goto_alpha_ij(d_ij, weights, i, j):
    """Computes alpha_ij, the Goto parameter for the singularity type.

    Finds the smallest alpha in [1, d_ij) such that
    q_{i*} * alpha - q_{j*} is divisible by d_ij, where i*, j* are
    the complement coordinates.

    Args:
        d_ij: The gcd of weights[i] and weights[j].
        weights: Tuple of 4 positive integer weights.
        i: First coordinate index.
        j: Second coordinate index.

    Returns:
        Positive integer alpha_ij in [1, d_ij).

    Raises:
        ValueError: If no valid alpha is found.
    """
    i_star, j_star = get_complement(i, j)
    q_is = weights[i_star]
    q_js = weights[j_star]
    for alpha_ij in range(1, d_ij):
        if (q_is * alpha_ij - q_js) % d_ij == 0:
            return alpha_ij
    raise ValueError(f"Unable to calculate alpha_ij for {weights}, i = {i}, j = {j}.")


class K3SingularityTerm(SageObject):
    """A singularity correction term for a coordinate pair (i, j).

    Represents the contribution of the cyclic quotient singularity
    A_{d_ij, alpha_ij} to the point count and zeta function of the
    minimal resolution of a weighted diagonal K3 surface.

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

        self.m_ij = get_goto_m_ij(p, m)
        self.d_ij = calculate_goto_d_ij(weights, i, j)
        self.e_ij = calculate_goto_e_ij(weights, i, j)
        self.f_ij = calculate_goto_f_ij(self.m_ij, self.e_ij)
        self.omega_ij = calculate_goto_omega_ij(m, self.e_ij, self.f_ij)

        self.alpha_ij = calculate_goto_alpha_ij(self.d_ij, weights, i, j)
        self.eta_ij = calculate_goto_eta_ij(self.f_ij, self.CC)
        self.r_ij = calculate_goto_r_ij(self.d_ij, self.alpha_ij)

    def _latex_(self):
        """Returns LaTeX representation for Sage display."""
        return f"A_{{{self.d_ij}, {self.alpha_ij}}} ({self.alpha_ij}, {self.f_ij})"

    def _repr_(self):
        """Returns plain text representation for Sage display."""
        return f"A_({self.d_ij}, {self.alpha_ij}) ({self.omega_ij})"

    def get_zeta_contribution(self, p, nu):
        """Computes the singularity contribution to the zeta function denominator.

        Returns prod_{k=0}^{f-1} (1 - eta^k * p^nu * x)^{r * omega}.

        Args:
            p: Prime number.
            nu: Extension degree.

        Returns:
            Polynomial in CC['x'].
        """
        term = 1
        CRng.<x> = PolynomialRing(self.CC)
        for f_ij_power in range(self.f_ij):
            term *= (1 - (self.eta_ij^f_ij_power * p)^nu * x)
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


def collect_k3_singularity_terms(p, m, weights, CC):
    """Collects all singularity correction terms for a K3 surface.

    Iterates over coordinate pairs (i, j) and creates K3SingularityTerm
    objects for those with d_ij = gcd(q_i, q_j) >= 2.

    Args:
        p: Prime number.
        m: Degree of the hypersurface.
        weights: Tuple of 4 positive integer weights.
        CC: Complex field for computation.

    Returns:
        List of K3SingularityTerm objects.
    """
    singularity_terms = []

    for i, j in iterate_coordinate_pairs():
        d_ij = calculate_goto_d_ij(weights, i, j)
        if d_ij <= 1:
            continue
        singularity_terms.append(K3SingularityTerm(p, m, weights, i, j, CC))

    return singularity_terms


def count_points_goto(p, nu, m, weights, CC=None):
    """Counts F_{p^nu}-points on the minimal resolution of a weighted diagonal K3.

    Uses Goto's Theorem 5.2 for the surface
    x_0^{m/q_0} + ... + x_3^{m/q_3} = 0 in P^3(Q) with twist c = (1,1,1,1).

    Computes:
        N_nu(X~) = 1 + q^nu + q^{2*nu}
                   + sum_{a in A} j(a)^nu
                   + sum_{(i,j) in J_0} omega_ij * r_ij * {q^nu + (eta*q)^nu + ...}

    where q = p and J_0 = {(i,j) : gcd(q_i, q_j) >= 2}.

    Args:
        p: Prime number with p = 1 (mod m).
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

    # Smooth part: 1 + q^nu + q^{2*nu} + sum of Jacobi sum powers
    total = 1 + q^nu + q^(2*nu)

    jacobi_terms = generate_diagonal_k3_jacobi_terms(m, weights)
    for a in jacobi_terms:
        character_fractions = [a_i / m for a_i in a]
        js = calculate_jacobi_sum(K, character_fractions)
        total += js^nu

    # Singularity correction from resolution of cyclic quotient singularities
    sing_terms = collect_k3_singularity_terms(p, m, weights, CC)
    for sing in sing_terms:
        total += sing.get_point_count(p, nu)

    return total


def compute_goto_zeta_denominator(p, m, weights, CC=None):
    """Computes the zeta function denominator using Goto's Theorem 5.2.

    Assembles the denominator as prod(1 - alpha_i * t) over all
    Frobenius eigenvalues: trivial cohomology contributions (H^0, H^4),
    Jacobi sum eigenvalues (H^2 smooth part), and singularity correction
    eigenvalues (exceptional divisors).

    Args:
        p: Prime number with p = 1 (mod m).
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

    # Trivial cohomology: (1 - t)(1 - pt)(1 - p^2 t)
    denom = (1 - u) * (1 - p * u) * (1 - p^2 * u)

    # Jacobi sum eigenvalues: prod_{a in A} (1 - j(a) t)
    jacobi_terms = generate_diagonal_k3_jacobi_terms(m, weights)
    for a in jacobi_terms:
        character_fractions = [a_i / m for a_i in a]
        js = CC(calculate_jacobi_sum(K, character_fractions))
        denom *= (1 - js * u)

    # Singularity correction: prod_{(i,j)} prod_{k=0}^{f-1} (1 - eta^k p t)^{r*omega}
    sing_terms = collect_k3_singularity_terms(p, m, weights, CC)
    for sing in sing_terms:
        for k in range(sing.f_ij):
            denom *= (1 - sing.eta_ij^k * p * u) ^ (sing.r_ij * sing.omega_ij)

    # Round complex coefficients to integers
    Qt = PolynomialRing(QQ, 't')
    t = Qt.gen()
    result = Qt(0)
    for deg, coef in enumerate(list(denom)):
        result += Integer(coef.real().round()) * t^deg

    return result


def compute_goto_zeta_denominator_extension(p, f, m, weights, CC=None):
    """Computes the zeta function denominator over F_{p^f} using Goto's formula.

    Base change raises all Frobenius eigenvalues to the f-th power:
    trivial eigenvalues become 1, p^f, p^{2f}; Jacobi sums become
    j(a)^f; singularity eigenvalues become (eta^k * p)^f.

    Args:
        p: Prime number with p = 1 (mod m).
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

    # Trivial cohomology: (1 - T)(1 - q T)(1 - q^2 T)
    denom = (1 - u) * (1 - q * u) * (1 - q^2 * u)

    # Jacobi sum eigenvalues raised to f-th power
    jacobi_terms = generate_diagonal_k3_jacobi_terms(m, weights)
    for a in jacobi_terms:
        character_fractions = [a_i / m for a_i in a]
        js = CC(calculate_jacobi_sum(K, character_fractions))
        denom *= (1 - js^f * u)

    # Singularity correction: eigenvalues (eta^k * p)^f
    sing_terms = collect_k3_singularity_terms(p, m, weights, CC)
    for sing in sing_terms:
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
