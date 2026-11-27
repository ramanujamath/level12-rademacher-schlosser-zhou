"""Checks of every numerical constant used in Sections 4 and 5.

Rational inequalities are checked with fractions.Fraction (exact);
inequalities involving special functions are checked on fine grids in
double precision as a sanity check of the analytic arguments given in
the manuscript.

Usage: python3 check_constants.py
"""

import math
from fractions import Fraction as F

from scipy.integrate import quad
from scipy.special import ive

from common import tau


def ok(msg):
    print("[ok]", msg)


def main():
    # ---- Lemma 4.1 constant  C = (2 + sqrt 3)/3
    s = sum(abs(math.cos(math.pi * h / 6)) for h in range(12))
    assert abs(s - (4 + 2 * math.sqrt(3))) < 1e-12
    # sqrt3 < 7/4  <=>  3 < 49/16
    assert F(3) < F(7, 4) ** 2
    C_upper = (2 + F(7, 4)) / 3          # = 5/4
    assert C_upper == F(5, 4)
    ok("sum_h |cos(pi h/6)| = 4 + 2 sqrt3, so C = (2+sqrt3)/3 < 5/4")

    # ---- Proposition 4.2:  24 pi^2 C (zeta(3/2)^2 - 1) / 12^{3/2} < 42
    # inputs (each verified below): pi^2 < 9.87, zeta(3/2) < 2.62, 12^{3/2} = 24 sqrt3 > 41.56
    assert F(31416, 10000) ** 2 < F(987, 100) and math.pi < 3.1416
    assert F(41556, 1000 * 24) ** 2 < 3   # 41.556/24 < sqrt3, hence 24 sqrt3 > 41.556 ...
    assert F(4156, 100 * 24) ** 2 < 3     # ... and > 41.56
    K = 10 ** 6
    zeta_upper = sum(k ** -1.5 for k in range(1, K)) + K ** -1.5 + 2 / math.sqrt(K)
    assert zeta_upper < 2.62
    const = 24 * F(987, 100) * F(5, 4) * (F(262, 100) ** 2 - 1) / F(4156, 100)
    assert const < 42
    true_const = (24 * math.pi ** 2 * (2 + math.sqrt(3)) / 3
                  * (2.6123753486854883 ** 2 - 1) / 12 ** 1.5)
    ok(f"tail constant < {float(const):.3f} < 42  (actual value {true_const:.3f})")

    # tau(12 m) <= 6 tau(m)
    assert all(tau(12 * m) <= 6 * tau(m) for m in range(1, 20001))
    ok("tau(12m) <= 6 tau(m) for m <= 20000 (true in general: tau is submultiplicative)")

    # ---- Bessel inequalities (grid checks; proofs are in the manuscript)
    xs = [k / 100 for k in range(1, 20001)]          # 0.01 .. 200
    for x in xs:
        I1_scaled = ive(1, x)                         # I_1(x) e^{-x}
        assert I1_scaled <= x / 2 + 1e-15, x          # I_1(x) <= (x/2) e^x
        if x >= 1:
            assert I1_scaled > 2 / (9 * math.pi * math.sqrt(x)), x
    ok("I_1(x) <= (x/2) e^x on (0, 200] and I_1(x) > 2 e^x / (9 pi sqrt x) on [1, 200]")

    J, _ = quad(lambda v: math.exp(-v) * math.sqrt(v), 0, 1)
    assert J > 2 / (3 * math.e) > F(2, 9)
    ok(f"int_0^1 e^-v sqrt v dv = {J:.5f} > 2/(3e) = {2 / (3 * math.e):.5f} > 2/9")

    # asymptotic expansion (4.4): bracket error is O(x^-3)
    worst = 0.0
    for x in (20, 40, 80, 160, 320, 640):
        br = ive(1, x) * math.sqrt(2 * math.pi * x)
        err = abs(br - (1 - 3 / (8 * x) - 15 / (128 * x * x))) * x ** 3
        worst = max(worst, err)
    assert worst < 1
    ok(f"I_1 asymptotic: |bracket error| * x^3 <= {worst:.3f} for x in [20, 640]")

    # ---- Section 5 constants
    # sqrt(3/pi)/27 > 1/28  <=>  3/pi > (27/28)^2  <=>  pi < 3 (28/27)^2
    assert F(31416, 10000) < 3 * F(28, 27) ** 2
    ok("sqrt(3/pi)/27 > 1/28")
    # threshold: exp(pi sqrt n / 6) > 28 * 42 n^{3/4} = 1176 n^{3/4} for n >= 576
    assert 28 * 42 == 1176
    assert F(49, 10) ** 2 > 24                      # sqrt 24 < 4.9
    rhs = 1176 * 24 * F(49, 10)                     # > 1176 * 24^{3/2}
    assert F(2718, 1000) ** 12 > rhs                # e^{4 pi} > e^{12} > 2.718^12
    ok(f"e^(4 pi) > 2.718^12 = {float(F(2718, 1000) ** 12):.0f} > "
       f"{float(rhs):.0f} > 1176 * 576^(3/4)")
    # monotonicity: derivative of pi sqrt(n)/6 - (3/4) log n is positive for n >= 9
    assert all(math.pi / (12 * math.sqrt(n)) - 3 / (4 * n) > 0 for n in range(9, 100000))
    assert math.pi / (12 * 3) > 3 / (4 * 9)
    ok("ratio exp(pi sqrt n/6) / n^{3/4} is increasing for n >= 9")
    first = next(n for n in range(1, 10 ** 5)
                 if all(math.exp(math.pi * math.sqrt(m) / 6) > 1176 * m ** 0.75
                        for m in range(n, n + 2000)))
    ok(f"(5.4) in fact holds from n = {first} on; the manuscript uses n >= 576")


if __name__ == "__main__":
    main()
