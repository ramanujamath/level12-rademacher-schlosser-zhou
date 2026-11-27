"""Checks on the restricted Kloosterman sums of Sections 3 and 4.

For 12 | c:
    K12(n; c)    = sum_{d mod c, (d,c)=1, d = +-1 (12)} e((n d - dbar)/c)
    S(A, B; c)   = sum_{d mod c, (d,c)=1} e((A d + B dbar)/c)
    S_chi(A,B;c) = same with weight chi12(d)

Checked for every c = 12, 24, ..., CMAX and every n mod c:
  (a) K12(n;12) = 2 cos(pi (n-1)/6)                      (eq. (1.7))
  (b) K12 = (S(n,-1) + S_chi(n,-1)) / 2                  (eq. (3.1))
  (c) S(n,-1;c) = 0 for n even                           (Lemma 3.1)
  (d) S_chi(n,-1;c) = 0 for n = 1 (mod 3)                (Lemma 3.3)
  (e) K12(n;c) = 0 for n = 4 (mod 6)                     (Theorem 1.2)
  (f) the Chinese-remainder factorisation (Lemma 3.2) of S_chi
  (g) K12 = (1/6) sum_h cos(pi h/6) S(n + h c/12, -1; c) (eq. (4.2))
  (h) |S(m,-1;c)| <= tau(c) sqrt(c)                      (Weil bound, used)
  (i) |K12(n;c)| <= ((2+sqrt3)/3) tau(c) sqrt(c)         (Lemma 4.1, eq. (4.1))

Floating point is used for the exponential sums; tolerances are 1e-8.

Usage: python3 check_kloosterman.py [CMAX]
"""

import cmath
import math
import sys

from common import inverse, tau, units

TOL = 1e-8


def chi12(d):
    r = d % 12
    return 1 if r in (1, 11) else (-1 if r in (5, 7) else 0)


def chi3(d):
    r = d % 3
    return 0 if r == 0 else (1 if r == 1 else -1)


def chi4(d):
    r = d % 4
    return 0 if r % 2 == 0 else (1 if r == 1 else -1)


def e(x):
    return cmath.exp(2j * math.pi * x)


def sums_for_modulus(c):
    """Return functions n -> K12, S, S_chi (as complex) for fixed c."""
    U = units(c)
    inv = {d: inverse(d, c) for d in U}

    def S(A, B, weight=lambda d: 1):
        return sum(weight(d) * e(((A * d + B * inv[d]) % c) / c) for d in U)

    def K(n):
        return sum(e(((n * d - inv[d]) % c) / c) for d in U if d % 12 in (1, 11))

    return S, K


def twisted(A, B, c, psi):
    U = units(c)
    return sum(psi(d) * e(((A * d + B * pow(d, -1, c)) % c) / c) for d in U)


def main(CMAX=360):
    C_const = (2 + math.sqrt(3)) / 3
    worst_ratio = 0.0
    worst_weil = 0.0
    for c in range(12, CMAX + 1, 12):
        S, K = sums_for_modulus(c)
        t = tau(c)
        # factorisation c = 3^beta * c'
        beta, cp = 0, c
        while cp % 3 == 0:
            cp //= 3
            beta += 1
        q3 = 3 ** beta
        lam = pow(cp, -1, q3)       # inverse of c' mod 3^beta
        mu = pow(q3, -1, cp)        # inverse of 3^beta mod c'
        for n in range(c):
            k = K(n)
            s0 = S(n, -1)
            s1 = S(n, -1, chi12)
            assert abs(k.imag) < TOL
            assert abs(k - (s0 + s1) / 2) < TOL, (c, n)
            if c == 12:
                assert abs(k - 2 * math.cos(math.pi * (n - 1) / 6)) < TOL
            if n % 2 == 0:
                assert abs(s0) < TOL, ("Lemma 3.1", c, n)
            if n % 3 == 1:
                assert abs(s1) < TOL, ("Lemma 3.3", c, n)
            if n % 6 == 4:
                assert abs(k) < TOL, ("Theorem 1.2", c, n)
            # (f) CRT factorisation, for a few n per c
            if n < 6:
                lhs = s1
                rhs = (twisted(lam * n % q3, -lam % q3, q3, chi3)
                       * twisted(mu * n % cp, -mu % cp, cp, chi4))
                assert abs(lhs - rhs) < 1e-7, ("CRT", c, n)
                if n % 3 == 1:
                    assert abs(twisted(lam * n % q3, -lam % q3, q3, chi3)) < TOL
            # (g) additive-character expansion
            g = sum(math.cos(math.pi * h / 6) * S(n + h * c // 12, -1)
                    for h in range(12)) / 6
            assert abs(k - g) < 1e-7, ("additive", c, n)
            # (h), (i)
            worst_weil = max(worst_weil, abs(s0) / (t * math.sqrt(c)))
            worst_ratio = max(worst_ratio, abs(k) / (t * math.sqrt(c)))
        print(f"  c = {c:4d}: all identities hold for every n mod c")
    assert worst_weil <= 1 + TOL
    assert worst_ratio <= C_const + TOL
    print(f"[ok] (1.7), (3.1), Lemmas 3.1-3.3, Theorem 1.2 and identity (4.2) "
          f"verified for 12 | c <= {CMAX}")
    print(f"[ok] max |S(n,-1;c)|/(tau(c) sqrt c)  = {worst_weil:.4f} <= 1")
    print(f"[ok] max |K12(n;c)|/(tau(c) sqrt c)   = {worst_ratio:.4f} "
          f"<= (2+sqrt3)/3 = {C_const:.4f}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 360)
