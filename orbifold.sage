"""Orbifold cohomology and trace formula computation.

Computes symmetry groups G_A, lambda-sector bases, and Frobenius trace
contributions for weighted Fermat hypersurfaces.
"""

load("potentials.sage")

import functools

from sage.structure.sage_object import SageObject


def _render_monomial(variable_name, exponents, latex):
    """Renders a monomial as a string with subscripted variables.

    Args:
        variable_name: Base name for variables (e.g. "x" or "y").
        exponents: Tuple of integer exponents for each variable.
        latex: If True, use LaTeX formatting; otherwise plain text.

    Returns:
        String representation of the monomial, omitting zero-exponent terms.
    """
    parts = []
    for index, one_exponent in enumerate(exponents):
        if one_exponent == 0:
            continue
        subscript = index + 1
        if latex:
            parts.append(f"{variable_name}_{{{subscript}}}^{{{one_exponent}}}")
        else:
            parts.append(f"{variable_name}_{subscript}^{one_exponent}")
    return " ".join(parts)


def _render_wedge(indices, latex):
    """Renders an exterior product of basis vectors.

    Args:
        indices: Sequence of coordinate indices for the wedge product.
        latex: If True, use LaTeX formatting; otherwise plain text.

    Returns:
        String representation of the exterior product.
    """
    parts = []
    for one_index in indices:
        subscript = one_index + 1
        if latex:
            parts.append(f"e_{{{subscript}}}")
        else:
            parts.append(f"e_{subscript}")
    return " ".join(parts)


@functools.lru_cache(maxsize=None)
def _p_adic_gamma_value(gi, p, prec):
    """Cached p-adic Gamma evaluation for a single rational value.

    Args:
        gi: Rational number at which to evaluate Gamma_p.
        p: Prime number.
        prec: p-adic precision.

    Returns:
        p-adic Gamma value Gamma_p(gi) in Zp(p, prec).
    """
    R = Zp(p, prec)
    return R(gi).gamma("pari")


class CohomologyBasisElement(SageObject):
    """A single element of the orbifold cohomology basis.

    Represents a basis element x^gamma * y^lambda * e_I in the
    lambda-sector of the orbifold cohomology, where x^gamma is the
    Koszul part, y^lambda is the twisted sector label, and e_I is the
    exterior algebra part.

    Attributes:
        x_gamma: Tuple of integer exponents for the Koszul monomial.
        y_lambda: Tuple of integer exponents for the twisted sector.
        ext_coordinates: Tuple of coordinate indices fixed by lambda.
        full_matrix: The potential matrix A.
    """

    def __init__(self, x_gamma, y_lambda, ext_coordinates, full_matrix):
        """Initializes a cohomology basis element.

        Args:
            x_gamma: Exponent tuple for the Koszul monomial x^gamma.
            y_lambda: Exponent tuple for the twisted sector label y^lambda.
            ext_coordinates: Tuple of coordinate indices fixed by lambda,
                determining the exterior algebra factor.
            full_matrix: The potential matrix A defining the hypersurface.
        """
        self.x_gamma = x_gamma
        self.y_lambda = y_lambda
        self.ext_coordinates = ext_coordinates
        self.full_matrix = full_matrix

    def to_tuple(self):
        """Converts this element to a hashable tuple representation.

        Returns:
            Triple of tuples (x_gamma, y_lambda, ext_coordinates).
        """
        return (tuple(self.x_gamma), tuple(self.y_lambda), tuple(self.ext_coordinates))

    def _render(self, latex):
        """Renders this element as a string.

        Args:
            latex: If True, use LaTeX formatting; otherwise plain text.

        Returns:
            String representation combining x, y, and exterior parts.
        """
        x_part = _render_monomial("x", self.x_gamma, latex)
        y_part = _render_monomial("y", self.y_lambda, latex)
        ext_part = _render_wedge(self.ext_coordinates, latex)
        return f"{x_part} {y_part} {ext_part}".strip()

    def _latex_(self):
        """Returns LaTeX representation for Sage display."""
        return self._render(latex=True)

    def _repr_(self):
        """Returns plain text representation for Sage display."""
        return self._render(latex=False)

    @functools.cached_property
    def y_coordinates(self):
        """Coordinate indices not fixed by lambda.

        Returns:
            Generator of indices i not in ext_coordinates.
        """
        return (
            i for i in range(self.full_matrix.nrows()) if i not in self.ext_coordinates
        )

    @functools.cached_property
    def lambda_matrix(self):
        """Sub-matrix of A restricted to non-fixed coordinates.

        Returns:
            Square sub-matrix of full_matrix at the y_coordinates.
        """
        return self.full_matrix.matrix_from_rows_and_columns(
            self.y_coordinates, self.y_coordinates
        )

    @functools.cached_property
    def age_lam(self):
        """Age of the twisted sector element y^lambda.

        Returns:
            Rational number giving the age of lambda w.r.t. A.
        """
        return get_lambda_age(self.y_lambda, self.full_matrix)

    @functools.cached_property
    def age_gamma(self):
        """Age of x^gamma.

        Returns:
            Rational number giving the age of gamma w.r.t. A.
        """
        return get_gamma_age(self.x_gamma, self.full_matrix)

    @functools.cached_property
    def dim_lam(self):
        """Dimension of the lambda-sector fixed locus.

        Counts how many coordinates of d * lambda * A^{-T} are divisible
        by d.

        Returns:
            Non-negative integer.
        """
        d = get_d(self.full_matrix)
        scaling = d * vector(self.y_lambda) * ~self.full_matrix.T
        return len([z for z in scaling if (z % d) == 0])

    def p_adic_gamma(self, R):
        """Computes the p-adic Gamma product for this element.

        Evaluates prod_i Gamma_p(gamma_i) where gamma_i are the entries
        of gamma * A^{-1}.

        Args:
            R: A p-adic ring (e.g. Zp(p, prec)).

        Returns:
            Product of p-adic Gamma values in R.
        """
        p = R.prime()
        prec = R.precision_cap()
        total = R(1)
        for gi in vector(self.x_gamma) * ~self.full_matrix:
            total *= _p_adic_gamma_value(Rational(gi), p, prec)
        return total

    def trace_contribution(self, p, R):
        """Computes this element's contribution to the orbifold trace formula.

        The contribution is:
            (-1)^{dim(lambda)} * p^{age(lambda)-1} * (-p)^{age(gamma)}
            * prod Gamma_p(gamma * A^{-1}).

        Args:
            p: Prime number.
            R: A p-adic ring (e.g. Zp(p, prec)).

        Returns:
            The trace contribution as a p-adic number.
        """
        trace = (-1) ** (self.dim_lam)
        trace *= p ** (self.age_lam - 1)
        trace *= (-p) ** self.age_gamma
        trace *= self.p_adic_gamma(R)
        return trace

    def twisted_eigenvalue(self, p, R):
        """Eigenvalue of p^{age+n-1} H(Fr_A) on this cohomology class.

        From the paper's Proposition 3.4 and eq 3.19, H(Fr_A) acts
        diagonally with eigenvalue epsilon (eq 3.20). The twisted
        operator p^{age+n-1} H(Fr_A) then has eigenvalue:

            alpha_i = p^{age(lambda)-1} * (-p)^{age(gamma)} * Gamma_p(gamma A^{-1})

        For nu > 1, the nu-th iterate [p^{age+n-1} H(Fr_A)]^nu has
        eigenvalue alpha_i^nu on each basis element.

        The trace contribution (eq 3.18) bakes in the supertrace sign:

            tc_i = (-1)^{dim(lambda)} * alpha_i

        so alpha_i = (-1)^{dim(lambda)} * tc_i = (-1)^{s+r} * tc_i.

        Args:
            p: Prime number.
            R: A p-adic ring (e.g. Zp(p, prec)).

        Returns:
            The twisted Frobenius eigenvalue as a p-adic number.
        """
        return (-1) ** self.dim_lam * self.trace_contribution(p, R)

    @property
    def hodge_number(self):
        """Hodge number (s, r) for this cohomology element.

        Returns:
            Tuple (s, r) where s = age_lam + age_gamma - 1 and
            r = dim_lam - 1 + age_lam - age_gamma.
        """
        return (
            self.age_lam + self.age_gamma - 1,
            self.dim_lam - 1 + self.age_lam - self.age_gamma,
        )


def reduce_modulo_lattice(v, A, A_inverse=None):
    """Reduces a vector to the fundamental domain of the lattice G_A.

    Here, G_A = Z^n / Z^n A^T. Computes v mod the lattice by subtracting
    the floor of A^{-1} v times A.

    Args:
        v: Integer vector or tuple to reduce.
        A: Square integer matrix defining the lattice.
        A_inverse: Precomputed inverse of A, or None to compute it.

    Returns:
        Tuple of integers representing v in the fundamental domain.
    """
    _A_inverse = A.inverse() if A_inverse == None else A_inverse
    vv = vector(v)
    floors = [floor(m) for m in _A_inverse * vv]
    return tuple(vv - A * vector(floors))


def get_lambda_basis(A):
    """Returns symmetries of A represented as elements of G_A.

    Enumerates all elements of the quotient group G_A = Z^n / Z^n A^T
    by brute-force reduction of all tuples modulo the lattice.

    Args:
        A: Square integer matrix defining the potential.

    Returns:
        List of tuples, each an element of G_A.
    """
    lambdas = set()
    d = get_d(A)
    for one_lambda in Tuples(range(d), A.nrows()).list():  # Slow...
        lambdas.add(tuple(reduce_modulo_lattice(one_lambda, A, ~A)))
    return list(lambdas)


def generate_symmetry_group(group_generators, A):
    """Computes the subgroup of G_A generated by the given generators.

    Elements are represented as lambda uniquely written as
    zeta_d^(d * lambda * A^{-T}), where d = get_d(A).

    Args:
        group_generators: List of tuples, each a generator of the
            desired subgroup of G_A.
        A: Square integer matrix defining the potential.

    Returns:
        Set of tuples representing the generated subgroup.
    """
    n = A.nrows()
    d = get_d(A)
    group = set()
    all_tuples = [[]]
    sources = [range(d)] * len(group_generators)
    for source in sources:
        all_tuples = [t + [s] for t in all_tuples for s in source]
    # Scale all basis vectors in all ways and add to set (ensures uniqueness).
    for one_tuple in all_tuples:
        g = vector([0] * n)
        for i in range(len(one_tuple)):
            g = (g + one_tuple[i] * vector(group_generators[i])) % d
        group.add(tuple(reduce_modulo_lattice(g, A)))
    return group


def generate_dual_group(symmetry_group, A):
    """Computes the dual of a symmetry group with respect to A^T.

    An element mu of G_{A^T} is in the dual if for all lambda in the
    symmetry group, lambda * A^{-T} * mu is integral.

    Args:
        symmetry_group: Set of tuples representing the symmetry group.
        A: Square integer matrix defining the potential.

    Returns:
        Set of tuples representing the dual symmetry group.
    """
    A_T = A.T
    A_T_inverse = ~A_T
    G_AT = get_lambda_basis(A_T)

    def is_dual(mu):
        for lam in symmetry_group:
            if not is_zero_mod_one(vector(lam) * A_T_inverse * vector(mu)):
                return False
        return True

    dual_symmetry_group = set()
    for mu in G_AT:
        if is_dual(mu):
            dual_symmetry_group.add(mu)
    return dual_symmetry_group


def is_zero_mod_one(value):
    """Checks whether a rational number is an integer.

    Args:
        value: A rational number.

    Returns:
        True if value equals floor(value), i.e. value is integral.
    """
    return value == floor(value)


def get_fixed_coordinates(A, y_lambda, d):
    """Finds coordinates fixed by the lambda-sector element.

    A coordinate i is fixed if the i-th entry of d * lambda * A^{-T}
    is divisible by d.

    Args:
        A: Square integer matrix defining the potential.
        y_lambda: Tuple of integers representing the lambda-sector element.
        d: The minimal degree d = get_d(A).

    Returns:
        List of integer indices of the fixed coordinates.
    """
    mu = d * vector(y_lambda) * A.inverse().transpose()
    a = [i for i in range(A.nrows()) if is_zero_mod_one(mu[i]) and (mu[i] % d) == 0]
    return a


def is_invariant(x_gamma, symmetry_group, full_matrix, d=None):
    """Checks if x^gamma is invariant under the symmetry group.

    Tests whether d * gamma * A^{-1} * lambda is divisible by d for
    every lambda in the symmetry group.

    Args:
        x_gamma: Integer vector of exponents for x^gamma.
        symmetry_group: Set of tuples representing the symmetry group.
        full_matrix: The potential matrix A.
        d: The minimal degree, or None to compute it.

    Returns:
        True if x^gamma is invariant under all group elements.
    """
    d = d or get_d(full_matrix)
    for y_lambda in symmetry_group:
        scaling = d * vector(x_gamma) * full_matrix.inverse() * vector(y_lambda)
        if not (scaling % d) == 0:
            return False
    return True


def get_gamma_age(x_gamma, A):
    """Computes the age of x^gamma with respect to A.

    The age is gamma * A^{-1} * epsilon, where epsilon = (1,...,1).

    Args:
        x_gamma: Integer vector of exponents for x^gamma.
        A: Square integer matrix defining the potential.

    Returns:
        Rational number giving the age.
    """
    epsilon = vector([1] * int(A.nrows()))
    return vector(x_gamma) * A.inverse() * epsilon


def get_lambda_age(y_lambda, A):
    """Computes the age of y^lambda with respect to A.

    The age is lambda * A^{-T} * epsilon, where epsilon = (1,...,1).

    Args:
        y_lambda: Integer vector of exponents for y^lambda.
        A: Square integer matrix defining the potential.

    Returns:
        Rational number giving the age.
    """
    epsilon = vector([1] * int(A.nrows()))
    return vector(y_lambda) * A.T.inverse() * epsilon


def expand_gamma(gamma, fixed_coordinates, n):
    """Embeds a reduced gamma vector back into full dimension.

    Places the entries of gamma at the positions specified by
    fixed_coordinates, with zeros elsewhere.

    Args:
        gamma: Integer vector of exponents in the reduced (fixed) coordinates.
        fixed_coordinates: List of coordinate indices where gamma entries
            should be placed.
        n: Total dimension of the full vector.

    Returns:
        List of integers of length n with gamma entries at
        fixed_coordinates and zeros elsewhere.
    """
    gamma_expanded = [0] * n
    for from_index, to_index in enumerate(fixed_coordinates):
        gamma_expanded[to_index] = gamma[from_index]
    return gamma_expanded


def get_lambda_sector_basis(full_matrix, symmetry_group, y_lambda, d):
    """Computes the cohomology basis elements for a single lambda-sector.

    For a given twisted sector y^lambda, restricts A to the fixed
    coordinates, computes the Koszul basis of the restricted matrix,
    and filters for invariance and integrality of the age.

    Args:
        full_matrix: The potential matrix A.
        symmetry_group: Set of tuples representing the symmetry group G.
        y_lambda: Tuple of integers for the lambda-sector element.
        d: The minimal degree d = get_d(A).

    Returns:
        List of CohomologyBasisElement objects for this sector.
    """
    n = full_matrix.nrows()
    fixed_coordinates = get_fixed_coordinates(full_matrix, y_lambda, d)
    A_lambda = full_matrix.matrix_from_rows_and_columns(
        fixed_coordinates, fixed_coordinates
    )

    if not fixed_coordinates:
        return [
            CohomologyBasisElement(
                x_gamma=[0] * n,
                y_lambda=y_lambda,
                ext_coordinates=fixed_coordinates,
                full_matrix=full_matrix,
            )
        ]

    elements = []
    for x_gamma in get_gamma_basis(A_lambda):
        expanded_gamma = expand_gamma(x_gamma, fixed_coordinates, n)
        if not is_invariant(
            expanded_gamma,
            symmetry_group,
            full_matrix,
            d,
        ):
            continue
        if not is_zero_mod_one(get_gamma_age(x_gamma, A_lambda)):
            continue
        elements.append(
            CohomologyBasisElement(
                x_gamma=expanded_gamma,
                y_lambda=y_lambda,
                ext_coordinates=fixed_coordinates,
                full_matrix=full_matrix,
            )
        )

    return elements


def count_points_orbifold(A, m, p, nu, R):
    """Computes the orbifold point count over F_{p^nu}.

    Computes the supertrace of [p^{age+n-1} H(Fr_A)]^nu acting on
    H(B_A, d_A + d_A^v):

        N_{p^nu} = sum_i (-1)^{s_i+r_i} * alpha_i^nu

    where alpha_i = twisted_eigenvalue() is the eigenvalue of the
    twisted Frobenius p^{age+n-1} H(Fr_A) (see eq 3.19 and Proposition
    3.4). The sign (-1)^{s+r} is the supertrace grading; the only
    nu-dependent term is the exponent on alpha_i.

    At nu=1 this is exactly the supertrace of p^{age+n-1} H(Fr_A)
    from Theorem 3.6 (eq 3.18). Since tc_i = (-1)^{dim(lambda)} *
    alpha_i and dim(lambda) = s+r mod 2, the sum reduces to
    sum_i tc_i = ST_p = N_p.

    Args:
        A: Square integer matrix defining the potential.
        m: The minimal degree d = get_d(A).
        p: Prime number with p = 1 (mod m).
        nu: Extension degree (count points over F_{p^nu}).
        R: A p-adic ring (e.g. Zp(p, prec)).

    Returns:
        Integer point count N_{p^nu}.
    """
    G = generate_symmetry_group([[1] * A.nrows()], A)
    all_elements = []
    for lam in G:
        for element in get_lambda_sector_basis(A, G, lam, m):
            all_elements.append(element)
    total = R(0)
    for element in all_elements:
        alpha = element.twisted_eigenvalue(p, R)
        s, r = element.hodge_number
        total += (-1)^(s + r) * alpha^nu
    return Integer(total)
