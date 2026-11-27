"""Run every verification script with its default parameters.

Usage: python3 run_all.py
"""

import time

import check_coefficients
import check_constants
import check_kloosterman
import check_modularity
import check_rademacher

STEPS = [
    ("Exact coefficients, sign law and zero set (Lemma 5.1, Thm 1.4, Cor 1.5)",
     lambda: check_coefficients.main(10000)),
    ("Modularity of T and its cusps (Section 2)", check_modularity.main),
    ("Kloosterman-sum identities and bounds (Sections 2-4)",
     lambda: check_kloosterman.main(360)),
    ("Exact formula and tail bound (Thm 1.1, Prop 4.2, Section 5)",
     lambda: (check_rademacher.part_a(150, 2400), check_rademacher.part_b(4000))),
    ("Numerical constants (Sections 4-5)", check_constants.main),
]

if __name__ == "__main__":
    t0 = time.time()
    for title, fn in STEPS:
        print(f"\n=== {title}")
        fn()
    print(f"\nAll checks passed in {time.time() - t0:.1f} s.")
