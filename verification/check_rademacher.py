"""Numerical checks of the exact formula (Theorem 1.1) and of the tail bound.

Part A (floating point, numpy):
  Truncate the series (1.6) at c <= CMAX and compare with the exact integers
  a(n) for 1 <= n <= NMAX.  Also check that the asymptotic expansion (1.8)
  is accurate.

Part B (high precision, decimal):
  Compute E(n) = a(n) - M(n) exactly enough, with
      M(n) = pi/(3 sqrt n) cos(pi (n-1)/6) I_1(pi sqrt(n)/3),
  and check |E(n)| < 42 exp(pi sqrt(n)/6)   (Proposition 4.2)
  and |M(n)| > |E(n)| for n >= 576 with cos != 0 (Section 5),
  for every n up to NB.

Usage: python3 check_rademacher.py [NMAX] [CMAX] [NB]
"""

import math
import sys
from decimal import Decimal, getcontext
from math import gcd

import numpy as np
from scipy.special import iv

from common import a_coeffs, predicted_sign


# ---------------------------------------------------------------- Part A
def truncated_rademacher(NMAX, CMAX):
    ns = np.arange(1, NMAX + 1, dtype=float)
    total = np.zeros(NMAX)
    for c in range(12, CMAX + 1, 12):
        d = np.array([x for x in range(c) if x % 12 in (1, 11) and gcd(x, c) == 1])
        dbar = np.array([pow(int(x), -1, c) for x in d])
        # K12(n;c) for all n at once (real because the set is closed under d -> -d)
        phase = (np.outer(np.arange(1, NMAX + 1), d) - dbar[None, :]) % c
        K = np.cos(2 * np.pi * phase / c).sum(axis=1)
        total += K / c * iv(1, 4 * np.pi * np.sqrt(ns) / c)
    return 2 * np.pi / np.sqrt(ns) * total


def part_a(NMAX=150, CMAX=2400):
    a = a_coeffs(NMAX)
    approx = truncated_rademacher(NMAX, CMAX)
    worst = max(abs(approx[n - 1] - a[n]) for n in range(1, NMAX + 1))
    for n in range(1, NMAX + 1):
        assert round(approx[n - 1]) == a[n], (n, approx[n - 1], a[n])
    print(f"[ok] series (1.6) truncated at c <= {CMAX} rounds to a(n) "
          f"for 1 <= n <= {NMAX} (max abs error {worst:.3f})")
    # asymptotic (1.8): relative error of the bracket
    big = a_coeffs(3001)
    for n in (1001, 2003, 3001):
        A = big[n]
        x = math.pi * math.sqrt(n) / 3
        lead = (math.cos(math.pi * (n - 1) / 6) / math.sqrt(6) * n ** -0.75
                * math.exp(x)
                * (1 - 9 / (8 * math.pi * math.sqrt(n))
                   - 135 / (128 * math.pi ** 2 * n)))
        rel = abs(A / lead - 1)
        print(f"     n = {n}: a(n)/asymptotic - 1 = {rel:.2e}  "
              f"(n^-3/2 = {n ** -1.5:.2e})")
        assert rel < 10 * n ** -1.5


# ---------------------------------------------------------------- Part B
def dec_pi(prec):
    getcontext().prec = prec + 10
    # Chudnovsky would be overkill; Machin's formula is enough here.
    def arctan_inv(x):
        x = Decimal(x)
        s, term, k, x2 = Decimal(0), 1 / x, 1, x * x
        while term != 0:
            s += term / k if (k // 2) % 2 == 0 else -term / k
            term /= x2
            k += 2
        return s
    return 4 * (4 * arctan_inv(5) - arctan_inv(239))


def dec_I1(x):
    """I_1(x) by its power series (all terms positive, so no cancellation)."""
    half = x / 2
    term = half
    s = term
    j = 0
    while True:
        j += 1
        term = term * half * half / (j * (j + 1))
        s += term
        if term < s * Decimal(10) ** (-getcontext().prec):
            return s


COS = {0: Decimal(1), 1: None, 2: Decimal("0.5"), 3: Decimal(0),
       4: Decimal("-0.5"), 5: None, 6: Decimal(-1), 7: None,
       8: Decimal("-0.5"), 9: Decimal(0), 10: Decimal("0.5"), 11: None}


def part_b(NB=10000):
    PI = dec_pi(80)
    getcontext().prec = 80
    r3 = Decimal(3).sqrt() / 2
    COS[1], COS[11], COS[5], COS[7] = r3, r3, -r3, -r3
    a = a_coeffs(NB)
    worst = Decimal(0)
    for n in range(1, NB + 1):
        sn = Decimal(n).sqrt()
        x = PI * sn / 3
        M = PI / (3 * sn) * COS[(n - 1) % 12] * dec_I1(x)
        E = Decimal(a[n]) - M
        ratio = abs(E) / (x / 2).exp()       # |E| / exp(pi sqrt n / 6)
        worst = max(worst, ratio)
        assert ratio < 42, n
        if n >= 576 and predicted_sign(n) != 0:
            assert abs(M) > abs(E), n
    print(f"[ok] |E(n)| < 42 exp(pi sqrt n/6) for 1 <= n <= {NB} "
          f"(largest observed ratio {float(worst):.4f})")
    print(f"[ok] |M(n)| > |E(n)| for 576 <= n <= {NB}, n != 4 (mod 6)")


if __name__ == "__main__":
    args = [int(x) for x in sys.argv[1:]]
    NMAX = args[0] if len(args) > 0 else 150
    CMAX = args[1] if len(args) > 1 else 2400
    NB = args[2] if len(args) > 2 else 4000
    part_a(NMAX, CMAX)
    part_b(NB)
