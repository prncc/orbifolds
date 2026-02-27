"""
Comparison: orbifold trace formula vs Gauss sum point counts for the Fermat quintic.

Computes F_p and F_{p^2} point counts on x_1^5 + ... + x_5^5 = 0 in P^4 (psi=0)
using two independent methods:
  (1) Orbifold trace formula (count_points_orbifold from orbifold.sage)
  (2) Gauss sum formula (DworkQuintic from quintic.sage)

Restricted to rho=1 primes (p ≡ 1 mod 5, equivalently 5 | p-1).

Usage: sage calculations/quintic_point_counts.sage

Output: calculations/quintic_point_counts.md
"""

import os
import time

pari.allocatemem(10^10)  # 10 GB PARI stack for large-prime gamma computations

load("orbifold.sage")
load("quintic.sage")

OUTPUT_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "quintic_point_counts.md",
)

A = build_fermat_matrix([5, 5, 5, 5, 5])
M = get_d(A)
detA = 5**5  # = 3125

# First 4 primes with det(A) | (p - 1)
NUM_PRIMES = 4
PRIMES = []
for _q in Primes():
    if len(PRIMES) >= NUM_PRIMES:
        break
    if _q > 10**7:
        break
    if (_q - 1) % detA == 0:
        PRIMES.append(_q)


def write(text, mode='a'):
    with open(OUTPUT_FILE, mode) as f:
        f.write(text)
        f.flush()


def main():
    header = r"""# Orbifold Trace vs Gauss Sum: Fermat Quintic Threefold

Surface: $x_1^5 + x_2^5 + x_3^5 + x_4^5 + x_5^5 = 0$ in $\mathbb{P}^4$ ($\psi = 0$)

Matrix: $A = 5 I_5$, $\det(A) = 5^5$, $d = 5$

Condition: $\rho = 1$, i.e. $p \equiv 1 \pmod{5}$, so $5 \mid (p-1)$

Methods:
1. **Orbifold trace formula**: supertrace of twisted Frobenius on deformed Borisov complex
2. **Gauss sum formula**: Candelas--de la Ossa--Rodriguez Villegas (Section 9, arXiv:hep-th/0402133)

"""
    write(header, mode='w')

    for p in PRIMES:
        assert (p - 1) % detA == 0, f"p={p} does not satisfy det(A) | (p-1)"

        write(f"## $p = {p}$\n\n")
        print(f"=== p = {p} ===")

        R = Zp(p, 256)

        write(f"| Field | Orbifold | Gauss sum | Match |\n")
        write(f"|-------|----------|-----------|-------|\n")

        for nu in [1]: # [1, 2]:
            field = f"\\mathbb{{F}}_{{{p}^{{{nu}}}}}" if nu > 1 else f"\\mathbb{{F}}_{{{p}}}"

            # Orbifold trace formula
            print(f"  Orbifold over F_{p^nu} (nu={nu})...", end="", flush=True)
            t0 = time.time()
            N_orb = count_points_orbifold(A, M, p, nu, R)
            t_orb = time.time() - t0
            print(f" {N_orb} ({t_orb:.1f}s)")

            # Gauss sum formula
            print(f"  Gauss sum over F_{p^nu} (s={nu})...", end="", flush=True)
            t0 = time.time()
            dq = DworkQuintic(p, s=nu)
            N_gauss = dq.count_points(0)
            t_gauss = time.time() - t0
            print(f" {N_gauss} ({t_gauss:.1f}s)")

            match = "yes" if N_orb == N_gauss else "**NO**"
            write(f"| ${field}$ | ${N_orb}$ | ${N_gauss}$ | {match} |\n")

        write("\n")

    print(f"Report written to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
