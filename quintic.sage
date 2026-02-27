"""
Point Counting on the Dwork Quintic Threefold over Finite Fields.

Implements the Gauss sum formulas from Section 9 of:
"Calabi-Yau Manifolds Over Finite Fields, I"
by P. Candelas, X. de la Ossa, and F. Rodriguez Villegas (2000).

The family of Calabi-Yau manifolds is defined in P^4 by:
x_1^5 + x_2^5 + x_3^5 + x_4^5 + x_5^5 - 5 * psi * x_1 * x_2 * x_3 * x_4 * x_5 = 0
"""

import itertools
from sage.all import (
    GF, ComplexField, Integer, exp, pi, I, discrete_log, PowerSeriesRing, QQ, prod
)

class DworkQuintic:
    def __init__(self, p, s=1, prec=1024):
        """
        Initializes the point counter for the Dwork Quintic over GF(p^s).

        Args:
            p (int): Prime characteristic (p != 5).
            s (int): Extension degree (q = p^s).
            prec (int): Precision for the ComplexField to evaluate Gauss sums
                        without rounding errors (1024 bits is safe for huge q).
        """
        if p == 5:
            raise ValueError("The quintic is degenerate at characteristic p=5.")
        if p == 2:
            raise NotImplementedError("Characteristic 2 requires modified Gauss sum signs.")

        self.p = p
        self.s = s
        self.q = p**s
        self.K = GF(self.q, 'a')
        self.CC = ComplexField(prec)
        self.g = self.K.multiplicative_generator()

        self._precompute_characters()
        self._gauss_cache = {0: self.CC(-1)}

    def _precompute_characters(self):
        """
        Precomputes the additive and multiplicative character arrays used
        to evaluate individual Gauss sums on demand.
        """
        q, p, g, CC = self.q, self.p, self.g, self.CC

        # Additive character Theta(x) = exp(2*pi*i*Tr(g^t)/p)
        self._theta = [CC(exp(2 * pi * I * Integer((g**t).trace()) / p)) for t in range(q - 1)]

        # Multiplicative character roots zeta^t
        zeta = CC(exp(2 * pi * I / (q - 1)))
        self._zeta_pows = [zeta**t for t in range(q - 1)]

    def _gauss_sum(self, m):
        """
        Returns the Gauss sum G_m, computing and caching on first access.

        G_m = sum_{t=0}^{q-2} theta[t] * zeta^{m*t}
        """
        m = m % (self.q - 1)
        if m not in self._gauss_cache:
            total = self.CC(0)
            for t in range(self.q - 1):
                total += self._theta[t] * self._zeta_pows[(m * t) % (self.q - 1)]
            self._gauss_cache[m] = total
        return self._gauss_cache[m]

    def count_points(self, psi_val):
        """
        Returns the exact number of projective F_q-rational points.
        """
        psi = self.K(psi_val)
        if psi == 0:
            return self._count_points_psi_zero()
        else:
            return self._count_points_psi_nonzero(psi)

    def _count_points_psi_zero(self):
        """Handles the purely diagonal Fermat hypersurface at psi = 0."""
        q = self.q
        trivial_pts = 1 + q + q**2 + q**3

        if (q - 1) % 5 != 0:
            return Integer(trivial_pts)

        k = (q - 1) // 5
        sum_jacobi = self.CC(0)

        for c in itertools.product(range(1, 5), repeat=5):
            if sum(c) % 5 == 0:
                term = self.CC(1)
                for ci in c:
                    term *= self._gauss_sum(ci * k)
                sum_jacobi += term

        pts = trivial_pts + Integer((sum_jacobi / self.CC(q)).real().round())
        return pts

    def _count_points_psi_nonzero(self, psi):
        """Handles the general deformation psi != 0."""
        q = self.q
        lam = (self.K(5) * psi)**(-5)
        L = discrete_log(lam, self.g)

        def teich(m):
            return self.CC(exp(2 * pi * I * ((m * L) % (q - 1)) / (q - 1)))

        # CASE 1: 5 does not divide (q - 1). Mapping x -> x^5 is a bijection.
        # (Follows Equation 9.8)
        if (q - 1) % 5 != 0:
            nu = self.CC(1 + q**4)
            for m in range(1, q - 1):
                term = (self._gauss_sum(m)**5 / self._gauss_sum(5 * m)) * teich(-m)
                nu += term
            nu_int = Integer(nu.real().round())
            return (nu_int - 1) // (q - 1)

        # CASE 2: 5 divides (q - 1). 5th roots of unity introduce conifold/toric periods.
        # (Follows Equations 9.13 and 9.19)
        else:
            k = (q - 1) // 5

            # (Monomial representation, Permutation Multiplicity gamma_v)
            V_list = [
                ((0,0,0,0,0), 1),
                ((4,1,0,0,0), 20),
                ((3,2,0,0,0), 20),
                ((3,1,1,0,0), 30),
                ((2,2,1,0,0), 30),
                ((4,3,2,1,0), 24)
            ]

            nu = self.CC(q**4)
            for v, gamma in V_list:
                for m in range(q - 1):
                    # Eq 9.19a: Smooth periods
                    if m % k != 0:
                        prod_G = self.CC(1)
                        for vj in v:
                            prod_G *= self._gauss_sum(-(m + k * vj))
                        sign = -1 if m % 2 != 0 else 1
                        beta = (sign / self.CC(q)) * self._gauss_sum(5 * m) * prod_G

                    # Eq 9.19b: Toric/Boundary periods
                    else:
                        a = m // k
                        z = sum(1 for vj in v if (vj + a) % 5 == 0)

                        delta_z0 = 1 if z == 0 else 0
                        num = -self.CC(q**(5 - z - delta_z0))

                        prod_G = self.CC(1)
                        for vj in v:
                            prod_G *= self._gauss_sum((vj + a) * k)
                        beta = num / prod_G

                    nu += gamma * beta * teich(m)

            nu_int = Integer(nu.real().round())
            return (nu_int - 1) // (q - 1)
