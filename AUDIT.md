# Audit of the original manuscript

This file records what was checked in the original draft
(`paper/original/level12_rademacher_schlosser_zhou_original.pdf`), what was
changed in the revised manuscript (`paper/level12_rademacher_schlosser_zhou.tex`),
and why. Equation and lemma numbers refer to the revised version unless stated.

## 1. Priority and novelty (most important)

The δ = −1 case of Schlosser–Zhou's Conjecture 24 has already been treated in
two preprints that the original draft did not cite:

- B. He and L. Li, arXiv:2509.10023 (Sept. 2025). They give the δ = −1 sign
  pattern for large n by the circle method, with the zero classes taken from the
  Andrews–Bressoud vanishing theorem.
- J. S. Padhi, arXiv:2609.22324 (Sept. 2026). This settles Conjecture 24 for
  δ = ±1 and 2 ≤ δ ≤ 3, and determines exactly which δ in (−1, 0) work. The
  δ = −1 case is proved by dominance for n ≥ 181, plus an exact check up to 400.

In addition, a(6m+4) = 0 is an instance of the Andrews–Bressoud vanishing
theorem (J. Austral. Math. Soc. 27 (1979)), and it was already recorded in OEIS A230536.

**Timeline (from the arXiv records):**

| Date | Event |
|---|---|
| 12 Sep 2025 | He–Li, arXiv:2509.10023 v1 (δ = −1 for large indices) |
| 16 Sep 2026 | Padhi, arXiv:2609.22324 v1 (complete treatment) |
| 22 Sep 2026 | This manuscript submitted to *The Ramanujan Journal* (per the author) |

Padhi's preprint appeared six days before submission, and the author was not
aware of it. The two works are therefore independent and essentially
simultaneous. The He–Li preprint predates the submission by a year, so it is
cited as earlier work on the large-index case.

**Change:** the abstract and introduction no longer present the result as the
first proof. A new subsection, "Relation to other work", cites He–Li, Padhi and
Andrews–Bressoud. It presents the paper as an independent proof by a different
method, with these features: an exact convergent formula, term-by-term
vanishing of the Rademacher series, fully explicit constants, and the complete
zero set.

The statement of Conjecture 24 was checked against the published paper
(Ramanujan J. 61 (2023)): for δ = −1 the pattern is `+++++0-----0`. This
matches the draft.

## 2. Mathematical content: checked, and correct

| Item | Status |
|---|---|
| T = q^{-1}/Q12 has the expansion q^{-1} + 1 + Σ a(n)qⁿ, with a(n) = b(n+1) | ✔ exact |
| T is Γ1(12)-invariant, and not invariant for d ≡ ±5 (mod 12) | ✔ numerical (`check_modularity.py`) |
| T is bounded at all cusps a/c with c \| 12, c < 12, so its only pole is at ∞ | ✔ numerical; T → 2+√3 at the cusp 0 |
| Double-coset parametrisation (Lemma 2.1) | ✔ proof re-derived |
| Exact formula (1.6), including its normalisation | ✔ truncated series reproduces a(n), n ≤ 150, to within 0.03 |
| K12(n;12) = 2cos(π(n−1)/6) | ✔ |
| Lemma 3.1 (2-adic involution) | ✔ proof and numerics |
| CRT factorisation | ✔ now stated precisely and proved (Lemma 3.2) |
| 3-adic vanishing (Lemma 3.3) | ✔ proof and numerics for all 12 \| c ≤ 360 |
| Asymptotic (1.8): constants 1/√6, 9/(8π), 135/(128π²) | ✔ re-derived; numerical relative error ≈ 10⁻⁶ at n = 3000 |
| Bessel inequalities (4.3), (5.3), integral representation (5.2) | ✔ |
| Recurrence (5.7) | ✔ agrees with direct product expansion |
| Sign law and zero set | ✔ exact, up to n = 10 000 |

## 3. Changes to the argument

1. **Twisted Kloosterman bound replaced.** The draft bounded
   S_χ12(n,−1;c) by citing Carmichael–Kimmel (Mathematika 2026, Lemma 2.5),
   with a conductor factor g(12) = √2. I could not verify that citation or its
   hypotheses. The revised Lemma 4.1 avoids twisted sums entirely: it detects
   d ≡ ±1 (mod 12) with additive characters, using identity (4.2), and then
   applies only the classical Weil bound (Iwaniec–Kowalski, Cor. 11.12). This
   gives the sharper constant (2+√3)/3 < 5/4 in place of 10/3. The
   Carmichael–Kimmel reference is dropped.
2. **New constants:** the tail bound is |E(n)| < 42 e^{π√n/6}
   (previously 160). The dominance threshold is n ≥ 576 (previously 784), so
   the finite check is now n ≤ 575. Every inequality is checked in
   `check_constants.py`, with rational ones checked in exact arithmetic.
3. **Lemma 3.2 (twisted vanishing) simplified.** The draft split this into
   the cases 24 | c and v₂(c) = 2. The 3-adic argument works for every
   12 | c, so the proof now has one case.
4. **Main-term lower bound:** now 1/28 (previously 1/33). Both are valid.

## 4. Minor issues fixed

- The CRT factorisation formula lost its overbars in the draft and did not say
  which inverses were meant. It is now stated precisely.
- "(n, −1, c) = 1" now appears as gcd notation.
- τ(c) (the divisor function) and τ ∈ ℍ are now explicitly distinguished.
- The definitions the draft assumed are now given: modular function, cusp,
  Hauptmodul, Kloosterman sum, χ12, I₁, and the meaning of "weak sign".
- The appendix code now checks n ≤ 575 and asserts a(n) ≠ 0 explicitly.
- The data-availability statement now points to this repository.
- The acknowledgment adds a sentence on Claude's role, alongside the existing
  ChatGPT statement. The author should confirm this wording.

## 5. Items relying on citation (not re-proved here)

- That T is a Hauptmodul for Γ1(12) (Yang 2004, Table 1, no. 13; also stated in
  OEIS). This is supported numerically.
- The Niebur/Zuckerman exact-formula framework. This is supported by the
  numerical agreement in (1).
- The Weil bound for composite moduli, |S(m,n;c)| ≤ τ(c)(m,n,c)^{1/2}c^{1/2}.
  The exact corollary number in Iwaniec–Kowalski should be confirmed against
  the printed book.
