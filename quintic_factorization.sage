"""
Zeta Function Factorization for the Dwork Quintic at psi = 0, rho = 1.
Replicates Table 12.1 from Candelas, de la Ossa, Rodriguez-Villegas (2000).
"""


def extract_candelas_factors_psi_zero(p):
    """
    Computes R_1(t) and R_A(pt, 0)^2 directly from Jacobi sums over F_p,
    and extracts the parameters a1, b1, c, d, alpha_pm, gamma_pm.
    """
    if p % 5 != 1:
        raise ValueError(f"This script handles rho=1 primes (p = 1 mod 5).")

    print(f"=== Extracting Zeta Factors for p={p}, psi=0 ===")

    # 1. Setup exact complex field and roots of unity
    # High precision is used to ensure perfect rounding back to integers
    CC = ComplexField(2048)
    zeta_5 = CC(exp(2 * pi * I / 5))
    zeta_p = CC(exp(2 * pi * I / p))

    # 2. Find a primitive root and map elements to their discrete logarithm
    K = GF(p)
    g = K.multiplicative_generator()
    ind = {Integer(g**t): t for t in range(p - 1)}

    # 3. Compute the 4 fundamental Gauss sums for the character of order 5
    # G[m] = sum_{x=1}^{p-1} zeta_5^{m * ind[x]} * zeta_p^x
    G = {}
    for m in range(1, 5):
        G[m] = sum(zeta_5**((m * ind[x]) % 5) * zeta_p**x for x in range(1, p))

    # 4. Form the roots (eigenvalues) for the two distinct orbits
    # By the Hasse-Davenport relation, the eigenvalues of the Frobenius operator
    # on the middle cohomology are lambda = - (1/p) * G(c1)G(c2)G(c3)G(c4)G(c5)
    # The polynomials are prod(1 - lambda * t) = prod(1 + u * t) where u = -lambda

    u1 = [] # Roots for R_1(t) -> The invariant tuple (1, 1, 1, 1, 1)
    for m in range(1, 5):
        u1.append((G[m]**5) / p)

    uA = [] # Roots for R_A(pt, 0)^2 -> Any asymmetric tuple summing to 0 mod 5
    for m in range(1, 5):
        c1, c3, c4 = (1*m) % 5, (3*m) % 5, (4*m) % 5
        uA.append((G[c1]**3 * G[c3] * G[c4]) / p)

    # 5. Form polynomials from the roots: P(t) = prod(1 + u*t)
    R = PolynomialRing(QQ, 't')
    t = R.gen()

    def poly_from_roots(u_roots):
        poly = R(1)
        for u in u_roots:
            poly *= (1 + u * t)
        # Round the complex coefficients to perfect integers
        return sum(Integer(coef.real().round()) * t**i for i, coef in enumerate(poly.list()))

    R1_poly = poly_from_roots(u1)
    QA_poly = poly_from_roots(uA) # This forms the single quartic R_A(pt, 0)^2

    # 6. Extract integer coefficients based on the paper's definitions
    # R_1(t) = 1 + a_1*t + b_1*p*t^2 + a_1*p^3*t^3 + p^6*t^4
    a1 = R1_poly[1]
    b1 = R1_poly[2] // p

    # R_A(pt, 0)^2 = 1 + c*p*t + d*p^2*t^2 + c*p^4*t^3 + p^6*t^4
    c_val = QA_poly[1] // p
    d_val = QA_poly[2] // (p**2)

    # 7. Solve for alpha_pm and gamma_pm over Q[sqrt(5)]
    x = polygen(QQ, 'x')
    K_sqrt5 = NumberField(x**2 - 5, 'sqrt5')

    # alpha^2 - a1*alpha + p*(b1 - 2*p^2) = 0
    disc_alpha = a1**2 - 4 * p * (b1 - 2*p**2)
    alpha_plus = (a1 + K_sqrt5(disc_alpha).sqrt()) / 2
    alpha_minus = (a1 - K_sqrt5(disc_alpha).sqrt()) / 2

    # gamma^2 - c*gamma + (d - 2*p) = 0
    disc_gamma = c_val**2 - 4 * (d_val - 2*p)
    gamma_plus = (c_val + K_sqrt5(disc_gamma).sqrt()) / 2
    gamma_minus = (c_val - K_sqrt5(disc_gamma).sqrt()) / 2

    return {
        'R1': R1_poly,
        'RA2': QA_poly,
        'a1': a1, 'b1': b1,
        'c': c_val, 'd': d_val,
        'alpha_pm': (alpha_plus, alpha_minus),
        'gamma_pm': (gamma_plus, gamma_minus)
    }

if __name__ == '__main__':
    # A = diag(5,5,5,5,5), det(A) = 5^5 = 3125
    # Conditions: det(A) | (p - 1) and rho = 1 (i.e. p ≡ 1 mod 5)
    detA = 5**5
    primes = []
    for p in Primes():
        if p <= 100000 and ((p - 1) % detA) == 0:
            primes.append(p)
            break
    results = []
    for p in primes:
        results.append(extract_candelas_factors_psi_zero(p))

    # Print table of parameters
    print("\nZeta Function Factorization for the Dwork Quintic at psi = 0")
    print("(Table 12.1 from Candelas, de la Ossa, Rodriguez-Villegas)")
    print()
    header = f"{'p':>5}  {'a1':>8}  {'b1':>8}  {'c':>8}  {'d':>8}  {'alpha+':>20}  {'alpha-':>20}  {'gamma+':>20}  {'gamma-':>20}"
    print(header)
    print("-" * len(header))
    for p, res in zip(primes, results):
        print(f"{p:>5}  {str(res['a1']):>8}  {str(res['b1']):>8}  {str(res['c']):>8}  {str(res['d']):>8}  {str(res['alpha_pm'][0]):>20}  {str(res['alpha_pm'][1]):>20}  {str(res['gamma_pm'][0]):>20}  {str(res['gamma_pm'][1]):>20}")

    # Print polynomials
    print()
    print(f"{'p':>5}  {'R_1(t)':>20}  {'R_A(pt,0)^2':>20}")
    print("-" * 50)
    for p, res in zip(primes, results):
        print(f"{p:>5}  {res['R1']}")
        print(f"{'':>5}  {res['RA2']}")
        print()