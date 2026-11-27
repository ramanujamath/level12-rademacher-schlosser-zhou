# A Rademacher expansion for a level-12 Hauptmodul and a sign conjecture of Schlosser and Zhou

**Author:** Sumit Kumar Jha (sumitkumarjha75@gmail.com)

This repository contains a paper and its verification code. The paper gives an
exact formula for the coefficients of

$$Q_{12}(q)^{-1}=\frac{(q^5,q^7;q^{12})_\infty}{(q,q^{11};q^{12})_\infty}=\sum_{N\ge0} b(N)q^N .$$

It uses that formula to prove the sign pattern `+ + + + + 0 − − − − − 0`
(mod 12) that Schlosser and Zhou conjectured, which is the case δ = −1 of
Conjecture 24 in *Ramanujan J.* 61 (2023). It also shows that b(6) is the
only zero outside the residue classes 5 and 11 mod 12.

The main ingredients are:

- an exact convergent Rademacher-type series for the coefficients of the
  Γ1(12) Hauptmodul T = q⁻¹/Q12 (OEIS A230536), whose terms involve Kloosterman
  sums restricted to d ≡ ±1 (mod 12);
- a quadratic-character splitting of these sums, which shows that every term
  of the series vanishes when n ≡ 4 (mod 6);
- an explicit tail bound that uses only the classical Weil bound, together
  with an exact computer check of the first 575 coefficients.

The δ = −1 case has also been treated, by the circle method, by He–Li
(arXiv:2509.10023) and Padhi (arXiv:2609.22324). This paper gives an
independent proof by a different route.

## Contents

| Path | What it is |
|---|---|
| `paper/level12_rademacher_schlosser_zhou.tex`, `.pdf` | Revised manuscript |
| `paper/original/` | The original draft, kept for reference |
| `AUDIT.md` | What was checked and changed, and why |
| `verification/` | Verification scripts |

## Running the checks

```bash
pip install -r requirements.txt
cd verification
python3 run_all.py        # about 15 s
```

| Script | Checks |
|---|---|
| `check_coefficients.py [N]` | Exact a(n) by recurrence and by direct expansion, against OEIS; the sign law and zero set (Lemma 5.1, Thm 1.4, Cor 1.5) up to N = 10 000 |
| `check_modularity.py` | Γ1(12)-invariance of T; boundedness at the other cusps |
| `check_kloosterman.py [C]` | (1.7), (3.1), Lemmas 3.1–3.3, Thm 1.2, (4.2) and the bound (4.1) for all 12 \| c ≤ C and all n mod c |
| `check_rademacher.py` | The truncated exact formula against a(n); the asymptotic (1.8); with 80-digit precision, Prop. 4.2 and main-term dominance for 576 ≤ n ≤ 4000 |
| `check_constants.py` | Every numerical inequality in Sections 4–5 |

To build the PDF, run `cd paper && latexmk -pdf level12_rademacher_schlosser_zhou.tex`.
