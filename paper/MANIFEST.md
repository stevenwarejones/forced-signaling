# Manifest — "The signaling cost of finite-speed hidden influences" (draft v1.13)

Every circulated copy of main.pdf must be accompanied by this package.

Solver for all floating-point LPs: scipy linprog (HiGHS) with the primal and dual
feasibility tolerances set EXPLICITLY to 1e-9 (the HiGHS defaults are 1e-7). Every
solve is checked for success; a failed, missing or nonfinite result raises LPFailure
rather than being read as a zero, and each solution's equality residual, inequality
violation and bound violation are recomputed from the returned vector and reported in
each script's closing diagnostics line. Decimal agreement reported anywhere in this
package is OBSERVED agreement, not a certified error bound: the exact-arithmetic
certificates for Theorems 1 and 2 are the only certified numbers here.

## Verify integrity (run from the package root)
    sha256sum -c hashes.txt
Hashes cover: main.tex, main.pdf, all scripts, figure source, figures, certificates, this manifest.

## Dependencies
python3 (3.9+), numpy, scipy (HiGHS via linprog), sympy (verify_Sigma.py only).
Each reproduction script prints the exact python/numpy/scipy/sympy versions it ran under.
Typical runtimes (machine-dependent; measured on python 3.11 / numpy 2.4 / scipy 1.17):
verify_K8.py 0.1-0.3 s (expected: "CERTIFICATE VALID: S4^op <= 6 + 8*Delta_sig...");
verify_Sigma.py 0.9-1.5 s (expected: "CERTIFICATES VALID: Sigma_HIC(Q_LC4) = (sqrt2-1)/4 exactly").
Both verifiers EXIT NONZERO on any failed check or malformed certificate, so they may be
used directly as automated gates.

## Reproduction commands
    # all commands run from the package root
    python3 paper/verify_K8.py              # Theorem 1 exact certificate (independent verifier)
    python3 paper/verify_Sigma.py           # Theorem 2 exact certificates (independent verifier; sympy)
    python3 paper/verify_directional.py     # Proposition 1: dual split + two one-sided certificates
    python3 threadB/reproduce_core.py       # QUICK: core numerical claims, minutes
    FULL=1 python3 threadB/reproduce_core.py  # complete: 512-spectrum, parallel 2-copy, 500k MC
    python3 threadB/reproduce_extra.py      # 4 supersets, random+random; FULL=1 adds LC5, mixed-flavor
    python3 threadB/reproduce_theorem4.py   # 24-vertex check, 400 constrained distances, perturbations
    (cd paper && python3 figures.py)        # regenerates the three figures
    python3 tests/regression_checks.py      # tamper/failure gates (see below)

QUICK vs FULL: QUICK mode validates a documented subset. It does NOT exercise the
512-completion spectrum or its multiplicities (it evaluates a 1/8 stride PLUS the K=8
optimal completion explicitly, since the plain stride omits both K=8 and K=10). Claims
marked (FULL) below are not reproduced by a QUICK run and must not be reported as such.

## Failure-mode gates (tests/regression_checks.py)
Verification is only meaningful if it fails when it should. This script checks, against
temporary copies (never the shipped certificates), that: a 1e-30 perturbation survives
certificate parsing rather than being snapped onto the nearby exact value; perturbed,
equality-breaking, sign-violating and malformed certificates are all REJECTED with a
nonzero exit, including under `python -O`; an unsuccessful, nonfinite-objective or
nonfinite-solution LP raises LPFailure instead of yielding Sigma = 0; and no
fallback-to-zero solver pattern remains anywhere in threadB/.

## Per-claim coverage inventory
| claim | script | status |
|---|---|---|
| Thm 1 certificate (K=8) | verify_K8.py | independent exact verifier |
| Prop 1: dual split (2 per direction) | verify_directional.py | exact; read off K8_certificate.json, so it presupposes verify_K8.py passes |
| Prop 1: one-sided attainment (dA,dD)=(s,0),(0,s) | verify_directional.py | exact primal certificates over Q(sqrt2); these are NEW proof obligations, not consequences of the Thm 1 dual |
| Thm 4: 24-vertex inequality | reproduce_theorem4.py | exhaustive (finite proof step) |
| Thm 4: 400 random constrained distances | reproduce_theorem4.py (seed 3) | covered; compared against max{0,(S-2)/8} -- 376 of the 400 have S<2, where the max is what makes the claim true |
| perturbation sweep near cluster point | reproduce_theorem4.py (seeds 13/17) | 100+50 in-package; the 300/150 sweep cited in the paper was an external audit (not regenerable here) |
| Thm 2 certificates (Q(sqrt2)) | verify_Sigma.py | independent exact standalone verifier |
| completion spectrum {8..16} | reproduce_core.py (FULL) | machine-numerical, asserted in FULL only (K=8 optimality is EXACT via Corollary 1); QUICK checks the optimal completion but not the multiplicities |
| LC4 numeric | reproduce_core.py | covered |
| tilt identity Sigma = max{0,(S4-6)/8} | reproduce_core.py | covered at 4 tested theta; S4 computed directly from the state and compared (max deviation ~1e-15). NOT a claim for all theta, and NOT a uniqueness claim for the S4 facet |
| fixed cluster-point dual = (S4-6)/8 on tilted states | reproduce_core.py | covered at the same 4 theta (max deviation ~6e-16); the dual is extracted at the cluster point and applied to each tilted state's marginals, including the negative value at theta=0.85 |
| GHZ/W/random = <1e-9 | reproduce_core.py (seed 7) | covered |
| chained n=3,4 | chained.py / reproduce_core.py | covered |
| parallel flatness | reproduce_core.py (FULL) | covered; asserted, floating-point only (no exact certificate) |
| MC delay cover | reproduce_core.py (seed 1) | covered (QUICK: 50k trials) |
| event budgets | reproduce_core.py | covered (allocation follows Li et al.; counts are ours) |
| CGLMP ladder d=2..8 | — | matches Brito et al. Fig.3; script OPEN ITEM |
| LC5 saturation | reproduce_extra.py (FULL) | covered |
| locked mixed-flavor saturation | reproduce_extra.py (FULL) | covered |
| Barnea evaluation | barnea.py (library) | WITHHELD from manuscript (transcription unverified) |

## Certificate file schema (Sigma_LC4_certificates.json)
Entries are exact elements of Q(sqrt2) as strings. Index space: primal vector of
length 385 = 256 model weights t[(x,w,a,d,fb,fg)] in lexicographic order of
(x,w,a,d,fb,fg) with fb,fg in 0..3 encoding B/C response functions (bit y of fb
= response to setting y), then 128 slack variables (16 signaling contexts x 8
outcomes, contexts ordered x-change then w-change, each over (w,y,z) resp.
(x,y,z) lexicographic; outcomes (b,c,d) resp. (a,b,c) lexicographic), then Delta
at index 384. Dual: 'dual_lambda_nonzero' indexes the 272 inequality rows (256 absolute-value rows, two per slack variable-pairing, plus 16 TV rows, in construction order);
'dual_mu' lists the 132 equality duals (4 normalization, 64 ABD, 64 ACD).

## Certificate file schema (directional_certificates.json)
Two primal vectors of length 386: the 385-column layout above with the single Delta
column split in two -- delta_A at 384 (charged by the eight A-switch contexts) and
delta_D at 385 (the eight D-switch contexts). Entries are [r,t] pairs of rational
strings meaning r + t*sqrt(2); omitted indices are zero. Both certificates satisfy
delta_A + delta_D = (sqrt(2)-1)/2 with the other direction exactly zero.

## Known open items (mirror of the manuscript's Statement on AI use)
- Independent HUMAN expert verification (none yet; AI cross-verification by
  Codex (OpenAI) has independently reconstructed and confirmed Theorems 1-2,
  the completion spectrum, tilt identity, cover numbers, and event budgets).
- Script for the one remaining OPEN ITEM row above (CGLMP ladder; values match Brito et al. Fig. 3).
- v1.6 correction: the earlier flag-parity assemblage claim (0.9508/0.990) was WITHDRAWN;
  the assemblage has no joint no-signaling extension (audit finding). Its driver was removed.
- Barnea transcription discrepancy.
- Proposition 1 is stated for the two early parties A and D of this scenario only.
  Whether the directional test pays experimentally depends on how asymmetric the
  achievable bounds are, and it requires SIMULTANEOUS valid confidence bounds on
  delta_A and delta_D where the scalar test needs one bound on their maximum; that
  statistical cost is not priced here.
- The error-relaxed extension of Proposition 1 (a three-parameter budget region
  admitting marginal deviations) is deliberately NOT included: its error parameter
  is mathematically defined but has no established operational reading, and it is
  not needed for any claim made here.
- The tilted-family identity is verified at four points only; neither the identity for
  all theta nor uniqueness of the active constraint is established (the earlier
  uniqueness wording was withdrawn in v1.12 -- it is false at theta=0.85, where
  S4 = 5.9669 < 6 and the zero is enforced by nonnegativity, not by the S4 facet).
- The experimental proposal of Sec. 6 is an architecture with the design obligations of
  Sec. 6.3 open (confidence-bound construction, simultaneous comparisons, the continuous
  candidate region, full budgets including the 8*Delta_sig term, setting randomization,
  and non-i.i.d./memory-valid treatment). Nothing in this package establishes that the
  experiment is ready to run.
- Whether the encoded conditionally-local optimization fully captures the intended
  finite-speed physical interpretation is a modelling question for specialist assessment,
  not something any certificate here settles.
