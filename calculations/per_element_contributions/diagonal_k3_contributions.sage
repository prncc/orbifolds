"""
Generate LaTeX document with per-element contribution tables for all 14 diagonal K3 surfaces.

Usage: sage calculations/per_element_contributions/diagonal_k3_contributions.sage

Output: calculations/per_element_contributions/diagonal_k3_contributions.tex
"""

import os

OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "diagonal_k3_contributions.tex")

load("diagonal_k3.sage")
load("orbifold.sage")
load("orbifold_display.sage")

# All 14 diagonal K3 surfaces from Goto Table 7
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


def find_prime(A):
    """Find the first prime p with det(A) | (p - 1)."""
    det_A = abs(A.det())
    for p in Primes():
        if p > 50000:
            return None
        if (p - 1) % det_A == 0:
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


def format_padic_expansion(padic_val, p, num_terms=3):
    """Format p-adic number as a + b*p + c*p^2 + ... in LaTeX."""
    try:
        # Try rational reconstruction first
        rat = padic_val.rational_reconstruction()
        int_val = Integer(rat)
        return f"${int_val}$"
    except (ArithmeticError, TypeError, ValueError):
        pass

    # Extract p-adic digits
    try:
        coeffs = []
        val = padic_val
        for i in range(num_terms):
            coeff = Integer(val.residue())
            coeffs.append(coeff)
            val = (val - coeff) / p

        # Format as a + b*p + c*p^2 + ...
        terms = []
        for i, c in enumerate(coeffs):
            if c == 0:
                continue
            if i == 0:
                terms.append(f"{c}")
            elif i == 1:
                if c == 1:
                    terms.append(f"{p}")
                else:
                    terms.append(f"{c}{{\\cdot}}{p}")
            else:
                if c == 1:
                    terms.append(f"{p}^{{{i}}}")
                else:
                    terms.append(f"{c}{{\\cdot}}{p}^{{{i}}}")

        if not terms:
            return "$0 + \\cdots$"
        return "$" + "{+}".join(terms) + "{+}\\cdots$"
    except Exception:
        return "$\\cdots$"


def format_equation_latex(weights, m):
    """Format the surface equation in LaTeX."""
    exponents = [m // w for w in weights]
    terms = " + ".join(f"x_{i}^{{{e}}}" for i, e in enumerate(exponents))
    return terms + " = 0"


def generate_surface_table(idx, weights, m, p, A, R):
    """Generate LaTeX table for one surface."""
    lines = []

    # Get all elements and group by Hodge number
    all_elements = get_orbifold_elements(A, m)
    hodge_to_elements = group_elements_by_hodge(all_elements)

    # Table header
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(f"\\caption{{Trace contributions for $" + format_equation_latex(weights, m) + f"$ at $p={p}$.}}\\label{{tab:diagonal-{idx}}}")
    lines.append(r"\renewcommand{\arraystretch}{1.1}")
    lines.append(r"{\scriptsize")
    lines.append(r"\begin{tabular}{@{}lllll@{}}")
    lines.append(r"\toprule")
    lines.append(r"$H^{s,r}$ & Element & $\mathrm{tc}_i$ & $\alpha_i$ & $\alpha_i$ expanded \\")
    lines.append(r"\midrule")

    # Group by Hodge number
    prev_hodge = None
    for hodge in sorted(hodge_to_elements.keys()):
        s, r = hodge
        elements = hodge_to_elements[hodge]

        if prev_hodge is not None:
            lines.append(r"\midrule")

        for i, element in enumerate(elements):
            # Hodge label (only on first element of group)
            hodge_str = f"$H^{{{s},{r}}}$" if i == 0 else ""

            # Element in LaTeX
            elem_latex = element._latex_()

            # Symbolic trace contribution
            tc_symbolic = render_contribution(element, simplify=True)

            # Symbolic eigenvalue
            alpha_symbolic = render_eigenvalue(element, simplify=True)

            # Expanded eigenvalue
            alpha_val = element.twisted_eigenvalue(p, R)
            alpha_expanded = format_padic_expansion(alpha_val, p)

            lines.append(f"{hodge_str} & ${elem_latex}$ & ${tc_symbolic}$ & ${alpha_symbolic}$ & {alpha_expanded} \\\\")

        prev_hodge = hodge

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"}")
    lines.append(r"\renewcommand{\arraystretch}{1.0}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


def main():
    """Generate full LaTeX document."""
    R_prec = 256

    # LaTeX preamble
    preamble = r"""\documentclass[11pt]{article}
\usepackage{booktabs}
\usepackage{amsmath,amssymb}
\usepackage[margin=0.75in]{geometry}
\usepackage{longtable}

\begin{document}

\title{Per-Element Trace Contributions for Diagonal K3 Surfaces}
\author{Generated by diagonal\_k3\_contributions.sage}
\maketitle

For each of the 14 diagonal K3 surfaces (Goto Table 7), we list all orbifold
cohomology basis elements with their twisted Frobenius eigenvalues $\alpha_i$.

"""

    with open(OUTPUT_FILE, 'w') as f:
        f.write(preamble)

    print(f"Writing to {OUTPUT_FILE}")

    num_surfaces = len(ALL_K3_SURFACES)
    for idx, (m, weights) in enumerate(ALL_K3_SURFACES, 1):
        exponents = [m // w for w in weights]
        A = build_fermat_matrix(exponents)

        p = find_prime(A)
        if p is None:
            print(f"[{idx}/{num_surfaces}] Q={weights}, m={m} - No valid prime found")
            continue

        R = Zp(p, R_prec)

        print(f"[{idx}/{num_surfaces}] Q={weights}, m={m}, p={p} ... ", end="", flush=True)

        # Section header for surface
        section = f"\n\\section*{{Surface \\#{idx}: $Q={weights}$, $m={m}$}}\n\n"
        section += f"Equation: ${format_equation_latex(weights, m)}$\n\n"
        with open(OUTPUT_FILE, 'a') as f:
            f.write(section)

        # Generate table
        table = generate_surface_table(idx, weights, m, p, A, R)

        with open(OUTPUT_FILE, 'a') as f:
            f.write(f"\\subsection*{{$p = {p}$}}\n\n")
            f.write(table)
            f.write("\n\n")

        print("OK")

        with open(OUTPUT_FILE, 'a') as f:
            f.write("\\newpage\n")

    # LaTeX footer
    footer = r"""
\end{document}
"""
    with open(OUTPUT_FILE, 'a') as f:
        f.write(footer)

    print(f"\nDone. Output: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
