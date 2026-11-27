"""Shared helpers for the verification scripts.

Notation follows the manuscript:

    Q12(q)   = (q, q^11; q^12)_inf / (q^5, q^7; q^12)_inf
    P(q)     = 1 / Q12(q) = sum_{r >= 0} p(r) q^r
    T(tau)   = q^{-1} P(q) = q^{-1} + 1 + sum_{n >= 1} a(n) q^n,   a(n) = p(n+1).

Everything in this module that produces coefficients uses exact integer
arithmetic only.
"""

from math import gcd


def rho(j):
    """Exponent of (1 - q^j) in P(q) = prod_j (1 - q^j)^rho(j)."""
    r = j % 12
    if r in (1, 11):
        return -1
    if r in (5, 7):
        return 1
    return 0


def p_by_recurrence(N):
    """p(0..N) from r p(r) = sum_{k=1}^r B(k) p(r-k), B(k) = -sum_{d|k} d rho(d)."""
    B = [0] * (N + 1)
    for d in range(1, N + 1):
        w = d * rho(d)
        if w:
            for k in range(d, N + 1, d):
                B[k] -= w
    p = [0] * (N + 1)
    p[0] = 1
    for r in range(1, N + 1):
        s = 0
        for k in range(1, r + 1):
            s += B[k] * p[r - k]
        if s % r:
            raise AssertionError(f"non-exact division at r={r}")
        p[r] = s // r
    return p


def p_by_product(N):
    """p(0..N) by multiplying out the product directly (independent check)."""
    c = [0] * (N + 1)
    c[0] = 1
    for j in range(1, N + 1):
        e = rho(j)
        if e == 1:  # multiply by (1 - q^j)
            for i in range(N, j - 1, -1):
                c[i] -= c[i - j]
        elif e == -1:  # divide by (1 - q^j): multiply by 1/(1-q^j)
            for i in range(j, N + 1):
                c[i] += c[i - j]
    return c


def a_coeffs(N):
    """a(1..N) as a dict, via the recurrence."""
    p = p_by_recurrence(N + 1)
    return {n: p[n + 1] for n in range(1, N + 1)}


def predicted_sign(n):
    """Sign of cos(pi (n-1)/6), computed exactly from (n-1) mod 12."""
    k = (n - 1) % 12
    if k in (0, 1, 2, 10, 11):
        return 1
    if k in (3, 9):
        return 0
    return -1


def units(c):
    return [d for d in range(c) if gcd(d, c) == 1]


def inverse(d, c):
    return pow(d, -1, c)


def tau(n):
    """Number of divisors."""
    t, m, p = 1, n, 2
    while p * p <= m:
        e = 0
        while m % p == 0:
            m //= p
            e += 1
        t *= e + 1
        p += 1
    if m > 1:
        t *= 2
    return t
