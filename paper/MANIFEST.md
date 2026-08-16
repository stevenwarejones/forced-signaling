# Manifest — "The signaling cost of finite-speed hidden influences" (draft v1.10)

Every circulated copy of main.pdf must be accompanied by this package.
Solver for all floating-point LPs: scipy linprog (HiGHS), feasibility tol ~1e-9.

## Verify integrity (run from the package root)
    sha256sum -c hashes.txt
Hashes cover: main.tex, main.pdf, all scripts, figure source, figures, certificates, this manifest.

## Dependencies
python3 (3.9+), numpy, scipy (HiGHS via linprog), sympy (verify_Sigma.py only).
Typical runtimes: verify_K8.py < 1 s (expected: "CERTIFICATE VALID: S4^op <= 6 + 8*Delta_sig...");
verify_Sigma.py seconds-to-minutes depending on machine, dominated by symbolic marginal
computation (expected: "CERTIFICATES VALID: Sigma_HIC(Q_LC4) = (sqrt2-1)/4 exactly").

## Reproduction commands
    # all commands run from the package root
    python3 paper/verify_K8.py              # Theorem 1 exact certificate (independent verifier)
    python3 paper/verify_Sigma.py           # Theorem 2 exact certificates (independent verifier; sympy)
    python3 threadB/reproduce_core.py       # QUICK: core numerical claims, minutes
    FULL=1 python3 threadB/reproduce_core.py  # complete: 512-spectrum, parallel 2-copy, 500k MC
    python3 threadB/reproduce_extra.py      # 4 supersets, random+random; FULL=1 adds LC5, mixed-flavor
    python3 threadB/reproduce_theorem4.py   # 24-vertex check, 400 constrained distances, perturbations
    (cd paper && python3 figures.py)        # regenerates the three figures

## Per-claim coverage inventory
| claim | script | status |
|---|---|---|
| Thm 1 certificate (K=8) | verify_K8.py | independent exact verifier |
| Thm 4: 24-vertex inequality | reproduce_theorem4.py | exhaustive (finite proof step) |
| Thm 4: 400 random constrained distances | reproduce_theorem4.py (seed 3) | covered |
| perturbation sweep near cluster point | reproduce_theorem4.py (seeds 13/17) | 100+50 in-package; the 300/150 sweep cited in the paper was an external audit (not regenerable here) |
| Thm 2 certificates (Q(sqrt2)) | verify_Sigma.py | independent exact standalone verifier |
| completion spectrum {8..16} | reproduce_core.py (FULL) | machine-numerical (K=8 optimality is EXACT via Corollary 1) |
| LC4 numeric / tilt identity | reproduce_core.py | covered |
| GHZ/W/random = <1e-9 | reproduce_core.py (seed 7) | covered |
| chained n=3,4 | chained.py / reproduce_core.py | covered |
| parallel flatness | reproduce_core.py (FULL) | covered |
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

## Known open items (mirror of the manuscript's Statement on AI use)
- Independent HUMAN expert verification (none yet; AI cross-verification by
  Codex (OpenAI) has independently reconstructed and confirmed Theorems 1-2,
  the completion spectrum, tilt identity, cover numbers, and event budgets).
- Script for the one remaining OPEN ITEM row above (CGLMP ladder; values match Brito et al. Fig. 3).
- v1.6 correction: the earlier flag-parity assemblage claim (0.9508/0.990) was WITHDRAWN;
  the assemblage has no joint no-signaling extension (audit finding). Its driver was removed.
- Barnea transcription discrepancy.
