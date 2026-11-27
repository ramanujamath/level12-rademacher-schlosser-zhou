"""Exact checks on the coefficients a(n) and b(N).

1. The recurrence (Section 5) agrees with direct expansion of the product.
2. The first terms agree with OEIS A230536.
3. For 1 <= n <= N (default 10000):
     a(n) = 0  <=>  n = 5 or n = 4 (mod 6),
     otherwise sgn a(n) = sgn cos(pi (n-1)/6).
   The manuscript needs this only for n <= 575; the larger range is a
   safety margin.
4. Corollary 1.5: the sign pattern of b(N), the coefficients of 1/Q12.

Integer arithmetic only.

Usage: python3 check_coefficients.py [N]
"""

import sys

from common import a_coeffs, p_by_product, p_by_recurrence, predicted_sign

# OEIS A230536, offset -1: a(-1), a(0), a(1), ...
OEIS = [1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, 0, 1, 2, 2, 2, 1, 0,
        -2, -3, -4, -4, -2, 0, 3, 5, 7, 6, 4, 0]


def main(N=10000):
    # 1. recurrence vs product
    M = 1500
    assert p_by_recurrence(M) == p_by_product(M), "recurrence != product"
    print(f"[ok] recurrence agrees with direct product expansion up to q^{M}")

    # 2. OEIS
    p = p_by_recurrence(len(OEIS))
    assert p[:len(OEIS)] == OEIS
    print(f"[ok] first {len(OEIS)} terms match OEIS A230536")

    # 3. sign law and zero set
    a = a_coeffs(N)
    zeros = []
    for n in range(1, N + 1):
        if a[n] == 0:
            zeros.append(n)
        if n == 5 or n % 6 == 4:
            assert a[n] == 0, n
        else:
            s = (a[n] > 0) - (a[n] < 0)
            assert s == predicted_sign(n) and s != 0, n
    extra = [n for n in zeros if n % 6 != 4]
    assert extra == [5], extra
    print(f"[ok] sign law of Theorem 1.4 holds for 1 <= n <= {N}; "
          f"zeros outside n = 4 (mod 6): {extra}")

    # 4. Corollary 1.5 for b(N) = p(N)
    pb = p_by_recurrence(N)
    for Nn in range(0, N + 1):
        r, v = Nn % 12, pb[Nn]
        if r in (0, 1, 2, 3, 4):
            assert v > 0, Nn
        elif r in (5, 11):
            assert v == 0, Nn
        elif r == 6:
            assert (v == 0) if Nn == 6 else (v < 0), Nn
        else:
            assert v < 0, Nn
    print(f"[ok] Corollary 1.5 (sign pattern of 1/Q12) holds for 0 <= N <= {N}")
    print(f"     b(0..11) = {pb[:12]}")
    print(f"     a({N}) has {len(str(abs(a[N - 1])))} digits")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 10000)
