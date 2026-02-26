"""Display and zeta function utilities for orbifold cohomology.

Provides LaTeX rendering of trace contributions, Hodge tables, zeta
function assembly from cohomology data, and Riemann hypothesis
verification (Weil conjecture checks).
"""

load("orbifold.sage")

import collections


def _render_gamma_factors(element):
    """Renders the Gamma_p factors for gamma * A^{-1} as LaTeX strings.

    Args:
        element: A CohomologyBasisElement.

    Returns:
        List of LaTeX strings, one per non-zero entry of gamma * A^{-1}.
    """
    x_gamma = element.x_gamma
    if not x_gamma:
        return []
    gamma_entries = vector(x_gamma) * element.full_matrix.inverse()
    parts = []
    for gamma_entry in gamma_entries:
        if gamma_entry != 0:
            gnum = gamma_entry.numerator()
            gdenom = gamma_entry.denominator()
            parts.append(
                "\\Gamma_p \\left( \\frac{%s}{%s} \\right)" % (gnum, gdenom)
            )
    return parts


def render_contribution(element, simplify=False):
    """Renders the trace formula contribution of a cohomology element as LaTeX.

    The trace contribution is:
        (-1)^{dim(lambda)} p^{age(lambda) - 1} (-p)^{age(gamma)}
        prod Gamma_p(gamma A^{-1})

    Args:
        element: A CohomologyBasisElement whose contribution to render.
        simplify: If True, collapse signs and p-powers into a single
            combined term.

    Returns:
        LaTeX string representing the trace contribution.
    """
    parts = []
    if simplify:
        sign_total = element.dim_lam + element.age_gamma
        if (sign_total % 2) == 1:
            parts.append("-")

        p_power = element.age_lam + element.age_gamma - 1
        if p_power == 1:
            parts.append("p")
        elif p_power > 1:
            parts.append(f"p^{{{p_power}}}")

        if parts:
            parts.append("\\,")
    else:
        parts.extend([
            "(-1)^{%i}" % element.dim_lam,
            "p^{%i - 1}" % element.age_lam,
            "(-p)^{%i}" % element.age_gamma,
        ])

    parts.extend(_render_gamma_factors(element))
    if not parts:
        return "1"
    return " ".join(parts)


def render_eigenvalue(element, simplify=False):
    """Renders the twisted Frobenius eigenvalue of a cohomology element as LaTeX.

    The twisted eigenvalue is (eq 3.19):
        p^{age(lambda) - 1} (-p)^{age(gamma)} prod Gamma_p(gamma A^{-1})

    This omits the (-1)^{dim(lambda)} supertrace sign present in the
    trace contribution (render_contribution).

    Args:
        element: A CohomologyBasisElement whose eigenvalue to render.
        simplify: If True, collapse signs and p-powers into a single
            combined term.

    Returns:
        LaTeX string representing the twisted eigenvalue.
    """
    parts = []
    if simplify:
        sign_total = element.age_gamma
        if (sign_total % 2) == 1:
            parts.append("-")

        p_power = element.age_lam + element.age_gamma - 1
        if p_power == 1:
            parts.append("p")
        elif p_power > 1:
            parts.append(f"p^{{{p_power}}}")

        if parts:
            parts.append("\\,")
    else:
        parts.extend([
            "p^{%i - 1}" % element.age_lam,
            "(-p)^{%i}" % element.age_gamma,
        ])

    parts.extend(_render_gamma_factors(element))
    if not parts:
        return "1"
    return " ".join(parts)


def group_elements_by_hodge(all_elements):
    """Groups cohomology elements by their Hodge number.

    Args:
        all_elements: List of CohomologyBasisElement objects.

    Returns:
        Dict mapping (s, r) Hodge number tuples to sorted lists of
        elements with that Hodge number.
    """
    hodge_to_elements = collections.defaultdict(list)
    for element in all_elements:
        hodge_to_elements[element.hodge_number].append(element)
    for hodge in hodge_to_elements:
        hodge_to_elements[hodge] = sorted(
            hodge_to_elements[hodge],
            key=lambda element: (element.x_gamma, element.y_lambda)
        )
    return hodge_to_elements


def build_contribution_table(hodge_to_elements, p, R, simplify=True):
    """Builds a table of trace contributions grouped by Hodge number.

    Args:
        hodge_to_elements: Dict mapping Hodge numbers to lists of
            CohomologyBasisElement objects (from group_elements_by_hodge).
        p: Prime number.
        R: A p-adic ring (e.g. Zp(p, prec)).
        simplify: If True, render simplified contributions.

    Returns:
        List of rows suitable for Sage's table() display, with header
        row followed by data rows.
    """
    header = ["Hodge Number", "Element", "Contribution", "TC p-adic"]
    rows = [header]
    seen_hodges = set()
    for hodge_number in sorted(hodge_to_elements):
        for element in hodge_to_elements[hodge_number]:
            tc = element.trace_contribution(p, R)
            rows.append([
                hodge_number if hodge_number not in seen_hodges else "",
                element,
                f"${render_contribution(element, simplify=simplify)}$",
                tc,
            ])
            seen_hodges.add(hodge_number)
    return rows


def compute_orbifold_zeta_denominator(all_elements, p, R):
    """Computes the orbifold zeta function denominator as a rational polynomial.

    Builds prod_i (1 - alpha_i * t) over all cohomology elements, using
    twisted Frobenius eigenvalues alpha_i (eq 3.19). Rational integer
    eigenvalues are accumulated exactly; the rest are accumulated p-adically
    and then reconstructed to rationals.

    Args:
        all_elements: List of CohomologyBasisElement objects.
        p: Prime number.
        R: A p-adic ring (e.g. Zp(p, prec)).

    Returns:
        Polynomial in QQ['t'] giving the zeta function denominator.
    """
    return _accumulate_zeta_product(all_elements, p, R)


def _accumulate_zeta_product(elements, p, R, f=1):
    """Accumulates prod_i (1 - alpha_i^f * t) for a list of elements.

    Uses twisted_eigenvalue() (= alpha_i, the twisted Frobenius
    eigenvalue from eq 3.19) rather than trace_contribution (= tc_i),
    since the standard zeta function is built from eigenvalues:

        zeta(t) = prod_{even s+r} (1 - alpha_i t)^{-1}
                * prod_{odd s+r}  (1 - alpha_i t)

    Handles rational and p-adic eigenvalues separately, then combines
    them into a single rational polynomial. When f > 1, each eigenvalue
    is raised to the f-th power (base change to F_{p^f}).

    Args:
        elements: List of CohomologyBasisElement objects.
        p: Prime number.
        R: A p-adic ring (e.g. Zp(p, prec)).
        f: Extension degree (default 1). Eigenvalues are raised to
            the f-th power for base change to F_{p^f}.

    Returns:
        Polynomial in QQ['t'].
    """
    Qt = PolynomialRing(QQ, 't')
    t = Qt.gen()
    Rx = PolynomialRing(R, 'x')
    x = Rx.gen()

    rational_terms = Qt(1)
    p_adic_terms = Rx(1)

    for element in elements:
        alpha = element.twisted_eigenvalue(p, R)
        alpha_f = alpha^f
        try:
            alpha_rational = alpha_f.rational_reconstruction()
            alpha_integer = Integer(alpha_rational)
            rational_terms *= (1 - alpha_integer * t)
        except (ArithmeticError, TypeError):
            p_adic_terms *= (1 - alpha_f * x)

    # Convert p-adic polynomial to rational.
    p_adic_as_rational = Qt(0)
    for deg, coef in enumerate(list(p_adic_terms)):
        p_adic_as_rational += coef.rational_reconstruction() * t**deg

    return rational_terms * p_adic_as_rational


def compute_orbifold_zeta_numerator_denominator(all_elements, p, R):
    """Computes orbifold zeta numerator and denominator separately.

    Builds the zeta function from twisted Frobenius eigenvalues alpha_i
    (see twisted_eigenvalue()), split by Hodge parity:

        numerator   = prod_{odd s+r}  (1 - alpha_i * t)
        denominator = prod_{even s+r} (1 - alpha_i * t)

    so that zeta(t) = numerator / denominator^{-1} in the standard
    convention where odd-cohomology factors appear in the numerator.

    Args:
        all_elements: List of CohomologyBasisElement objects.
        p: Prime number.
        R: A p-adic ring (e.g. Zp(p, prec)).

    Returns:
        Tuple (numerator, denominator) as polynomials in QQ['t'].
    """
    odd_elements = []
    even_elements = []
    for element in all_elements:
        s, r = element.hodge_number
        if (s + r) % 2 == 1:
            odd_elements.append(element)
        else:
            even_elements.append(element)

    numerator = _accumulate_zeta_product(odd_elements, p, R)
    denominator = _accumulate_zeta_product(even_elements, p, R)
    return numerator, denominator


def compute_orbifold_zeta_denominator_extension(all_elements, p, R, f):
    """Computes the orbifold zeta denominator over F_{p^f}.

    Each twisted eigenvalue alpha_i (computed over F_p) is raised to
    the f-th power before forming the factor (1 - alpha_i^f * T). This
    implements base change: the Frobenius over F_{p^f} is Frob_p^f,
    whose eigenvalues are the f-th powers of the F_p eigenvalues.

    Args:
        all_elements: List of CohomologyBasisElement objects.
        p: Prime number.
        R: A p-adic ring (e.g. Zp(p, prec)).
        f: Extension degree (positive integer).

    Returns:
        Polynomial in QQ['t'] giving the zeta denominator over F_{p^f}.
    """
    return _accumulate_zeta_product(all_elements, p, R, f=f)


def build_hodge_table(hodge_to_elements):
    """Builds a table of Hodge numbers h^{s,r} and their multiplicities.

    Args:
        hodge_to_elements: Dict mapping Hodge numbers to lists of
            CohomologyBasisElement objects.

    Returns:
        List of rows suitable for Sage's table() display.
    """
    rows = [["$h^{s,r}$", "Count"]]
    for hodge_number in sorted(hodge_to_elements):
        s, r = hodge_number
        rows.append([
            f"$h^{{{s},{r}}}$",
            len(hodge_to_elements[hodge_number]),
        ])
    return rows


def build_cohomology_table(hodge_to_elements, p, R, simplify=True):
    """Builds a table of contributions grouped by total cohomology degree H^n.

    Groups elements by n = s + r (total degree), then shows each
    H^{s,r} component within that group with its contributions.

    Args:
        hodge_to_elements: Dict mapping Hodge numbers to lists of
            CohomologyBasisElement objects.
        p: Prime number.
        R: A p-adic ring (e.g. Zp(p, prec)).
        simplify: If True, render simplified contributions.

    Returns:
        List of rows suitable for Sage's table() display.
    """
    # Group by total degree n = s + r
    degree_groups = collections.defaultdict(list)
    for hodge_number, elements in hodge_to_elements.items():
        s, r = hodge_number
        n = s + r
        degree_groups[n].append((hodge_number, elements))

    header = ["$H^n$", "$H^{s,r}$", "Element", "Contribution", "Value at $p$"]
    rows = [header]
    seen_degrees = set()
    for n in sorted(degree_groups):
        hodge_components = sorted(degree_groups[n])
        seen_hodges = set()
        for hodge_number, elements in hodge_components:
            s, r = hodge_number
            for element in elements:
                tc = element.trace_contribution(p, R)
                try:
                    tc_val = str(Integer(tc.rational_reconstruction()))
                except (ArithmeticError, TypeError):
                    tc_val = str(tc)
                rows.append([
                    f"$H^{n}$" if n not in seen_degrees else "",
                    f"$H^{{{s},{r}}}$" if hodge_number not in seen_hodges else "",
                    element,
                    f"${render_contribution(element, simplify=simplify)}$",
                    tc_val,
                ])
                seen_degrees.add(n)
                seen_hodges.add(hodge_number)
    return rows


def render_zeta_latex(zeta_denominator, p):
    """Renders the zeta function as a LaTeX product formula.

    Factors the denominator and writes Z(X, t) as 1 / prod(factors),
    grouping factors by cohomological weight determined from root norms.

    Args:
        zeta_denominator: Polynomial in QQ['t'] giving the zeta
            function denominator.
        p: Prime number.

    Returns:
        LaTeX string suitable for display with Math().
    """
    Qt = zeta_denominator.parent()
    t = Qt.gen()
    factored = zeta_denominator.factor()

    # Collect factors and classify by the root's absolute value
    CC = ComplexField(53)
    weight_factors = collections.defaultdict(list)

    for poly, mult in factored:
        # Find roots to classify weight
        roots = poly.roots(CC, multiplicities=False)
        if roots:
            norm = abs(1 / roots[0])
            if abs(norm - 1) < 0.01:
                w = 0
            elif abs(norm - float(CC(p))) < 0.01:
                w = 2
            elif abs(norm - float(CC(p)**2)) < 0.01:
                w = 4
            elif abs(norm - float(sqrt(CC(p)))) < 0.01:
                w = 1
            else:
                w = -1  # unknown
        else:
            w = -1
        weight_factors[w].append((poly, mult))

    parts = []
    parts.append("Z(X_{\\mathbb{F}_{%d}}, t) = " % p)

    # Build numerator / denominator pieces
    # For a surface: Z = P_2(t)^{(-1)} / ((1-t)(1-pt)(1-p^2 t))
    # But we have the full denominator, so Z = 1 / denom
    # Display as product of factors grouped by weight
    factor_strs = []
    for w in sorted(weight_factors):
        for poly, mult in weight_factors[w]:
            # Render polynomial factor
            poly_latex = _render_factor_latex(poly, mult, p)
            factor_strs.append(poly_latex)

    parts.append("\\frac{1}{")
    parts.append(" \\cdot ".join(factor_strs))
    parts.append("}")

    return "".join(parts)


def render_factored_poly_latex(poly, p):
    """Renders a factored polynomial in normalized (1 - alpha t) form.

    Drops the leading scalar unit from factorization and renders each
    factor using p-notation where possible.

    Args:
        poly: Polynomial in QQ['t'] to factor and render.
        p: Prime number for p-notation.

    Returns:
        LaTeX string suitable for display with Math().
    """
    factored = poly.factor()
    factor_strs = []
    for f, mult in factored:
        factor_strs.append(_render_factor_latex(f, mult, p))
    return " \\cdot ".join(factor_strs)


def _render_factor_latex(poly, mult, p):
    """Renders a single polynomial factor as LaTeX.

    Converts factors like (t - 1/5) to (1 - 5t) form, rendering
    eigenvalues in terms of p where possible.

    Args:
        poly: An irreducible polynomial factor.
        mult: Multiplicity of this factor.
        p: Prime number for p-notation.

    Returns:
        LaTeX string for this factor, e.g. "(1 - p^2 t)^{3}".
    """
    coeffs = list(poly)

    if poly.degree() == 1:
        # Linear: a + b*t, root at -a/b, eigenvalue tc = -b/a
        a, b = coeffs[0], coeffs[1]
        tc = Rational(-b / a)
        inner = _render_linear_inner(tc, p)
        if mult == 1:
            return f"({inner})"
        return f"({inner})^{{{mult}}}"

    else:
        # General degree: normalize so constant term is 1, render with p-notation
        a0 = coeffs[0]
        if a0 != 0:
            parts = ["1"]
            for deg in range(1, poly.degree() + 1):
                bk = Rational(coeffs[deg] / a0)
                coeff_latex = _signed_coeff_latex(bk, p)
                if not coeff_latex:
                    continue
                if deg == 1:
                    parts.append(f"{coeff_latex} t")
                else:
                    parts.append(f"{coeff_latex} t^{{{deg}}}")
            inner = " ".join(parts)
        else:
            inner = latex(poly / poly.leading_coefficient())
        if mult == 1:
            return f"({inner})"
        return f"({inner})^{{{mult}}}"


def _render_linear_inner(tc, p):
    """Renders the inside of a linear factor (1 - tc * t).

    Args:
        tc: Rational eigenvalue coefficient.
        p: Prime number for p-notation.

    Returns:
        LaTeX string like "1 - p^2 t" or "1 + 3p t".
    """
    tc_latex = _coeff_latex(tc, p)
    sep = " " if tc_latex else ""
    if tc > 0:
        return f"1 - {tc_latex}{sep}t"
    elif tc < 0:
        return f"1 + {_coeff_latex(-tc, p)}{sep}t"
    else:
        return "1"


def _coeff_latex(r, p):
    """Renders a rational coefficient using p-notation where possible.

    Returns the coefficient string without the variable t. For
    coefficient 1, returns empty string. Uses p-notation: p, p^2,
    -p, 2p, etc.

    Args:
        r: Rational coefficient value.
        p: Prime number for p-notation.

    Returns:
        LaTeX string for the coefficient, or empty string for 1.
    """
    r = Rational(r)
    if r == 0:
        return "0"
    if r == 1:
        return ""
    if r == -1:
        return "-"

    # Check if r is an integer
    if r in ZZ:
        n = Integer(r)
        return _integer_as_p_power(n, p)

    # Check if r is a ratio involving powers of p
    num = r.numerator()
    den = r.denominator()

    # Check denominator as power of p
    d = den
    p_exp = 0
    while d > 1 and d % p == 0:
        d //= p
        p_exp += 1

    if d == 1 and p_exp > 0:
        # Denominator is p^k
        num_latex = _integer_as_p_power(abs(num), p) if abs(num) != 1 else ""
        sign = "-" if num < 0 else ""
        if p_exp == 1:
            return f"{sign}{num_latex}p^{{-1}}"
        return f"{sign}{num_latex}p^{{-{p_exp}}}"

    return f"\\frac{{{num}}}{{{den}}}"


def _integer_as_p_power(n, p):
    """Renders an integer in terms of p, factoring out the maximal p-power.

    Args:
        n: Non-negative integer to render.
        p: Prime number for p-notation.

    Returns:
        String like "p^3", "151p^2", or the plain integer string.
    """
    n = Integer(abs(n))
    if n == 1:
        return "1"

    # Check if n is a power of p
    if n == p:
        return "p"
    k = 2
    pk = p * p
    while pk <= n:
        if pk == n:
            return f"p^{{{k}}}"
        k += 1
        pk *= p

    # Factor out maximal power of p: n = c * p^k
    remaining = n
    k = 0
    while remaining > 1 and remaining % p == 0:
        remaining //= p
        k += 1
    if k > 0:
        if k == 1:
            return f"{remaining}p"
        return f"{remaining}p^{{{k}}}"

    return str(n)


def _signed_coeff_latex(r, p):
    """Renders a rational coefficient with explicit sign for use in sums.

    Args:
        r: Rational coefficient value.
        p: Prime number for p-notation.

    Returns:
        String like "+ p^2" or "- 3p", or empty string for zero.
    """
    r = Rational(r)
    if r == 0:
        return ""
    if r > 0:
        coeff = _coeff_latex(r, p)
        return f"+ {coeff}" if coeff else "+ "
    coeff = _coeff_latex(-r, p)
    return f"- {coeff}" if coeff else "- "


def build_rh_table(rh_results, p):
    """Builds a table of Riemann hypothesis check results.

    Args:
        rh_results: List of (root, norm, weight, multiplicity) tuples
            from check_riemann_hypothesis.
        p: Prime number.

    Returns:
        List of rows suitable for Sage's table() display, showing
        each root's absolute value, expected weight, cohomology group,
        multiplicity, and pass/fail status.
    """
    header = ["$|\\alpha|$", "Weight", "Cohomology", "Mult", "Status"]
    rows = [header]

    # Group by (weight, norm) for cleaner display
    grouped = collections.defaultdict(lambda: [0, []])
    for root, norm, weight, mult in rh_results:
        key = (weight, round(norm, 4))
        grouped[key][0] += mult
        grouped[key][1].append(root)

    weight_labels = {
        0: "$H^0$",
        1: "$H^1$",
        2: "$H^2$",
        3: "$H^3$",
        4: "$H^4$",
    }

    for (weight, norm), (total_mult, _) in sorted(grouped.items(), key=lambda x: (x[0][0] is None, x[0][0] or 0, x[0][1])):
        if weight is not None:
            expected = float(p ** (weight / Integer(2)))
            status = "ok" if abs(norm - expected) < 0.01 else "FAIL"
            label = weight_labels.get(weight, f"$H^{{{weight}}}$")
            norm_str = f"$p^{{{weight}/2}}={norm:.4f}$" if weight % 2 == 1 else f"$p^{{{weight // 2}}}={Integer(round(norm))}$"
        else:
            status = "FAIL"
            label = "?"
            norm_str = f"${norm:.6f}$"
        rows.append([norm_str, weight, label, total_mult, status])

    return rows


def check_riemann_hypothesis(zeta_denominator, p, CC=None, max_weight=4):
    """Checks the roots of the zeta denominator.

    For a smooth projective variety of dimension d, eigenvalues of
    Frobenius on H^i have absolute value p^{i/2}. The zeta denominator
    roots satisfy |1/root| = p^{w/2} for weight w.

    Args:
        zeta_denominator: Polynomial in QQ['t'] giving the zeta
            function denominator.
        p: Prime number.
        CC: Complex field for root computation, or None for default
            256-bit precision.
        max_weight: Maximum cohomological weight to check against.

    Returns:
        List of (root, |1/root|, weight, multiplicity) tuples. Weight
        is None if the root doesn't match any expected value.
    """
    if CC is None:
        CC = ComplexField(256)
    Cv = PolynomialRing(CC, 'v')
    v = Cv.gen()

    complex_poly = Cv(0)
    for deg, coef in enumerate(list(zeta_denominator)):
        complex_poly += CC(coef) * v**deg

    results = []
    for root, multiplicity in complex_poly.roots():
        norm = abs(1/root)
        weight = None
        for w in range(max_weight + 1):
            expected = CC(p) ** (w / Integer(2))
            if abs(norm - expected) < 0.001:
                weight = w
                break
        results.append((root, float(norm), weight, multiplicity))
    return results
