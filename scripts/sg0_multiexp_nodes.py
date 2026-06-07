#!/usr/bin/env python3
"""Generate Gauss quadrature nodes/weights on (0,1) for weight ln^2(x),
used by the MultiExp radial grid of SG-0 (P. Gill, S. Chien,
J. Comput. Chem. 24 (2003) 732; S. Chien, P. Gill, J. Comput. Chem. 27
(2006) 730).  Moments: m_k = int_0^1 x^k ln^2 x dx = 2/(k+1)^3.

Method (same as Psi4's cubature.cc, but at high precision): Hankel
moment matrix -> Cholesky -> three-term recurrence coefficients
(Golub-Welsch) -> symmetric tridiagonal Jacobi matrix -> eigenvalues
(nodes) and weights m_0 * q_1i^2.

Radial map: r_i = -R ln(x_i), w_i = R^3 omega_i / x_i (w includes the
r^2 Jacobian; OpenQP's getSliceNonZero multiplies rAtm^3 = 1).
Validation: sum_i w_i exp(-r_i^2) = int_0^inf r^2 exp(-r^2) dr
          = sqrt(pi)/4 for R = 1.

Output: 17-digit Fortran parameter arrays, sorted by ascending r
(descending x).
"""
from mpmath import mp, mpf, matrix, cholesky, eigsy, sqrt, pi, exp, log

mp.dps = 200


def gauss_ln2(n):
    # Hankel moment matrix, (n+1) x (n+1)
    M = matrix(n + 1, n + 1)
    for i in range(n + 1):
        for j in range(n + 1):
            M[i, j] = mpf(2) / mpf(i + j + 1) ** 3
    L = cholesky(M)        # M = L L^T, L lower triangular
    R = L.T                # upper triangular, as in Golub-Welsch

    def r(i, j):
        return R[i, j]

    alpha = [None] * n
    beta = [None] * (n - 1)
    alpha[0] = r(0, 1) / r(0, 0)
    for j in range(1, n):
        alpha[j] = r(j, j + 1) / r(j, j) - r(j - 1, j) / r(j - 1, j - 1)
    for j in range(n - 1):
        beta[j] = r(j + 1, j + 1) / r(j, j)

    # Jacobi matrix
    J = matrix(n, n)
    for i in range(n):
        J[i, i] = alpha[i]
    for i in range(n - 1):
        J[i, i + 1] = beta[i]
        J[i + 1, i] = beta[i]

    E, Q = eigsy(J)        # ascending eigenvalues
    m0 = mpf(2)            # int_0^1 ln^2 x dx
    x = [E[i] for i in range(n)]
    w = [m0 * Q[0, i] ** 2 for i in range(n)]
    return x, w


def validate(n, x, w):
    # Gauss exactness: integrals of x^k ln^2 x for k = 0..2n-1
    maxerr_mom = mpf(0)
    for k in range(2 * n):
        s = sum(wi * xi ** k for xi, wi in zip(x, w))
        maxerr_mom = max(maxerr_mom, abs(s - mpf(2) / mpf(k + 1) ** 3))
    # Radial test: sum w_i/x_i exp(-r_i^2) = sqrt(pi)/4   (R = 1)
    s = sum((wi / xi) * exp(-log(xi) ** 2) for xi, wi in zip(x, w))
    err_rad = abs(s - sqrt(pi) / 4)
    return maxerr_mom, err_rad


def emit_fortran(name, vals, per_line=3):
    out = [f"  real(kind=dp), parameter :: {name}({len(vals)}) = [ &"]
    line = []
    body = []
    for v in vals:
        body.append(mp.nstr(v, 17, strip_zeros=False).replace('e', 'd')
                    if 'e' in mp.nstr(v, 17) else mp.nstr(v, 17, strip_zeros=False) + 'd0')
    for i in range(0, len(body), per_line):
        chunk = ", ".join(body[i:i + per_line])
        sep = ", &" if i + per_line < len(body) else "]"
        out.append(f"        {chunk}{sep}")
    return "\n".join(out)


for n in (23, 26):
    x, w = gauss_ln2(n)
    em, er = validate(n, x, w)
    print(f"n = {n}: max moment error = {mp.nstr(em, 3)}, "
          f"radial sqrt(pi)/4 error = {mp.nstr(er, 3)}")
    # Moments must be exact (rule correctness); the exp(-r^2) test is
    # limited by the intrinsic accuracy of the n-point rule itself
    # (e^{-r^2} is not a polynomial in x), a few 1e-9 for n = 23.
    assert em < mpf(10) ** -30 and er < mpf(10) ** -7
    # double-precision re-validation of the radial test
    import math
    xd = [float(v) for v in x]
    wd = [float(v) for v in w]
    s = sum((wi / xi) * math.exp(-math.log(xi) ** 2) for xi, wi in zip(xd, wd))
    print(f"  double precision: sum w_i exp(-r_i^2) = {s:.16f}, "
          f"target {float(sqrt(pi)/4):.16f}, err = {abs(s - float(sqrt(pi)/4)):.2e}")
    # ascending r = descending x
    xr = list(reversed(x))
    wr = list(reversed(w))
    print(emit_fortran(f"me{n}_x", xr))
    print(emit_fortran(f"me{n}_w", wr))
    print()
