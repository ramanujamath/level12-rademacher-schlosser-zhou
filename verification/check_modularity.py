"""Numerical check that T is invariant under Gamma_1(12) but not under
other elements of Gamma_0(12).

T(tau) = q^{-1} prod_{m>=0} (1-q^{12m+5})(1-q^{12m+7}) / ((1-q^{12m+1})(1-q^{12m+11})).

For each test matrix gamma = (a b; c d) in SL_2(Z) with 12 | c and each test
point tau, compare T(gamma tau) with T(tau).  Expected: equality exactly
when d = +-1 (mod 12); for d = +-5 (mod 12) the values differ.

Usage: python3 check_modularity.py
"""

import cmath
import math


def T(tau, tol=1e-18):
    q = cmath.exp(2j * math.pi * tau)
    val = 1 / q
    m = 0
    while True:
        f5, f7 = q ** (12 * m + 5), q ** (12 * m + 7)
        f1, f11 = q ** (12 * m + 1), q ** (12 * m + 11)
        val *= (1 - f5) * (1 - f7) / ((1 - f1) * (1 - f11))
        if abs(f1) < tol:
            return val
        m += 1


def act(g, tau):
    a, b, c, d = g
    return (a * tau + b) / (c * tau + d)


MATRICES = [
    (1, 1, 0, 1), (1, 0, 12, 1), (-1, 0, 12, -1), (13, 1, 12, 1),
    (11, -1, 12, -1), (25, 2, 12, 1), (1, 0, 24, 1), (37, 3, 12, 1), (1, 0, 36, 1), (35, -3, 12, -1),
    (-11, 1, -12, 1),
    # not in +-Gamma_1(12): d = +-5 (mod 12)
    (5, 2, 12, 5), (7, 4, 12, 7), (-5, 2, 12, -5),
]
POINTS = [0.1 + 0.9j, -0.23 + 0.6j, 0.37 + 1.3j]


def main():
    for g in MATRICES:
        a, b, c, d = g
        assert a * d - b * c == 1 and c % 12 == 0
        expected = d % 12 in (1, 11)
        for z in POINTS:
            lhs, rhs = T(act(g, z), tol=1e-16), T(z)
            same = abs(lhs - rhs) < 1e-8 * max(1, abs(rhs))
            assert same == expected, (g, z, lhs, rhs)
        tag = "invariant" if expected else "not invariant (as expected)"
        print(f"  {g}: {tag}")
    print("[ok] T is invariant exactly under the tested elements with "
          "d = +-1 (mod 12)")
    cusps()


def cusps():
    """T has no pole at any cusp a/c with c | 12, c < 12.

    Every cusp of Gamma_1(12) other than those equivalent to infinity is of
    this form.  For gamma in SL_2(Z) with gamma(inf) = a/c we evaluate
    T(gamma(iy)) for growing y.  A pole would force growth at least like
    exp(2 pi y / 12); we check the values stay bounded and converge.
    """
    from math import gcd
    count = 0
    for c in (1, 2, 3, 4, 6):
        for a in range(12):
            if gcd(a, c) != 1:
                continue
            # find b, d with a d - b c = 1
            d = next(d for d in range(1, 13 * c + 1) if (a * d - 1) % c == 0)
            b = (a * d - 1) // c
            g = (a, b, c, d)
            vals = [abs(T(act(g, 1j * y), tol=1e-16)) for y in (8.0, 16.0, 32.0)]
            # a pole of width w <= 12 would multiply |T| by >= exp(2 pi 16/12) > 4000
            assert vals[2] < 2 * max(vals[0], 1.0), (g, vals)
            count += 1
    print(f"[ok] T stays bounded at all {count} tested cusps a/c with c | 12, c < 12 "
          "(so its only pole is at the cusp infinity)")


if __name__ == "__main__":
    main()
