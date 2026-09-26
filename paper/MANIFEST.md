# Manifest — "The signaling cost of finite-speed hidden influences" (draft v1.15)

Every circulated copy of main.pdf must be accompanied by this package.

Solver for all floating-point LPs: scipy linprog (HiGHS) with the primal and dual
feasibility tolerances set EXPLICITLY to 1e-9 (the HiGHS defaults are 1e-7). Every
solve is checked for success; a failed, missing or nonfinite result raises LPFailure
rather than being read as a zero, and each solution's equality residual, inequality
violation and bound violation are recomputed from the returned vector and reported in
each script's closing diagnostics line. Decimal agreement reported anywhere in this
package is OBSERVED agreement, not a certified error bound: the exact-arithmetic
certificates for Theorems 1–2 and the attaining models for Propositions 1–2 are
checked exactly. The pinned Lean development additionally checks the results
listed below, including the entire white-noise curve. Numerical reproduction
rows do not inherit that proof status.

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
    python3 paper/verify_invisibility.py    # Proposition 2: pairwise-invisible attaining models
    # or all four in dependency order, failing loudly (run in a subshell/script):
    #   for v in K8 Sigma directional invisibility; do
    #     python3 paper/verify_$v.py || { echo "FAILED: verify_$v.py"; exit 1; }
    #   done
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
fallback-to-zero solver pattern remains anywhere in threadB/; and that
the support/coefficient validator is itself live and load-bearing (gate block D2): the
shipped certificate must report both structural facts; a broken primal must SUPPRESS the
reduced-program conclusion rather than merely exit nonzero; and mis-stating either the core
set or the 1/8 coefficient must REJECT the shipped, mathematically valid certificate, so the
gates fail if that validation is disabled.

## Numbered claims and formal verification scope

Lean revision: [`c10474e9b45cca2ea260eec2d9f348a12cf62677`](https://github.com/stevenwarejones/ontology-separation/tree/c10474e9b45cca2ea260eec2d9f348a12cf62677).
This is the merged main revision containing PRs #77–#81. All Lean references
below are pinned to it, and declaration names have the prefix
`OntologySeparation.`. The model class is finite stochastic conditional-local;
the extension to infinite hidden spaces and its finite-speed physical
interpretation are not kernel-checked. Values are expectation-level minima,
not finite-sample confidence bounds. No human expert has reviewed them and no
novelty or priority claim is made.

Each numbered assertion has **one** primary verification label. Where a proof
has a checked ingredient but an analytic conclusion, the overall assertion
keeps its analytic label; the ingredient appears separately below. Definitions
1–2 specify quantities, rather than assertions needing proof.

| Numbered claim | Verification label | Evidence and boundary |
|---|---|---|
| Theorem 1: measured-signaling tradeoff | Lean | `HiddenInfluence.sharp_tradeoff`, `HiddenInfluence.coefficient_optimal` in [SignalingTradeoff.lean](https://github.com/stevenwarejones/ontology-separation/blob/c10474e9b45cca2ea260eec2d9f348a12cf62677/OntologySeparation/Experiments/SignalingTradeoff.lean); independent `verify_K8.py` also retained. |
| Corollary 1: optimal slope over all completions | Lean | `HiddenInfluenceCompletion.coefficient_lower_bound_all_completions`, `optimal_completion_globally_sharp`, `stochastic_completion_globally_sharp`, `completion_count` in [ForcedSignalingCompletions.lean](https://github.com/stevenwarejones/ontology-separation/blob/c10474e9b45cca2ea260eec2d9f348a12cf62677/OntologySeparation/Experiments/ForcedSignalingCompletions.lean); 512 completions, including the stochastic lift. |
| Theorem 2: exact LC4 minimum | Lean | `ForcedSignalingTheorem2.exact_forced_signaling`, `exact_forced_signaling_stochastic_value` in [ForcedSignalingTheorem2.lean](https://github.com/stevenwarejones/ontology-separation/blob/c10474e9b45cca2ea260eec2d9f348a12cf62677/OntologySeparation/Experiments/ForcedSignalingTheorem2.lean); exact cluster-state target and finite stochastic bridge, plus `verify_Sigma.py`. |
| Theorem 3: quantum bound through the S4 facet | analytic (not machine-checked) | CHSH decomposition and Tsirelson's bound in the paper; not a bound on the full forced-signaling optimum. |
| Theorem 4: marginal-preserving local approximation | analytic (not machine-checked) | The paper supplies the construction and reduction; `reproduce_theorem4.py` exhaustively checks its 24-vertex inequality. That finite ingredient is not a machine proof of the full theorem. |
| Theorem 5: weak universal ceiling | analytic (not machine-checked) | Conditional application of Theorem 4 plus the triangle inequality; not a universal Lean theorem. |
| Theorem 6: sliced marginal-assemblage reduction | analytic (not machine-checked) | Model/selection equivalence and sliced-TV calculation in the paper. |
| Theorem 7: optimal delay cover | analytic (not machine-checked) | Arc-cover argument; `reproduce_core.py` supplies numerical/Monte Carlo checks, not a proof of the continuous statement. |
| Proposition 1: directional refinement and tightness | Lean | `ForcedSignalingPropositions.proposition1`, `proposition1_A_attains`, `proposition1_D_attains`, `stochastic_directional_bound`, `stochastic_lc4_directional_lower_bound` in [ForcedSignalingPropositions.lean](https://github.com/stevenwarejones/ontology-separation/blob/c10474e9b45cca2ea260eec2d9f348a12cf62677/OntologySeparation/Experiments/ForcedSignalingPropositions.lean); both exact directional certificates also checked by `verify_directional.py`. |
| Proposition 2: accessible signaling as a property of the layout | analytic (not machine-checked) | The geometric criterion and every-model branch use collectibility and Lemma 2. Lean checks its invisible attaining models only, as detailed below. |
| Lemma 1: one CHSH violation at a time | analytic (not machine-checked) | Elementary sign-pattern argument in the paper. |
| Lemma 2: only recipient sets containing the blind pair can carry signal | analytic (not machine-checked) | Marginalization of the fixed ABD/ACD families; FULL subset LPs are supporting numerical checks only. |
| Observation 1: support of every optimal dual | analytic (not machine-checked) | `verify_Sigma.py` checks the exhibited primal/dual support, coefficients and slack exactly; the extension to every optimal dual is the paper's complementary-slackness argument. |
| Observation 2: mechanism identity at four tilted points | numerical only | `reproduce_core.py`, QUICK or FULL; four specified theta values, not an all-theta theorem. |
| Conjecture 1: tight universal ceiling | conjecture | Adversarial numerical record below; the certified LC4 value does not prove the universal ceiling. |

### Checked ingredients and additional results

| Claim or ingredient | Verification label | Evidence and boundary |
|---|---|---|
| Proposition 2: all three pairwise-invisible attaining models and parity shift | Lean | `ForcedSignalingPropositions.proposition2_A`, `proposition2_D`, `proposition2_balanced`, `pairwise_invisible_optimum` in [ForcedSignalingPropositions.lean](https://github.com/stevenwarejones/ontology-separation/blob/c10474e9b45cca2ea260eec2d9f348a12cf62677/OntologySeparation/Experiments/ForcedSignalingPropositions.lean); `PairwiseInvisible` includes all six proper recipient projections in both early directions, B/C silence, and pure parity. No spacetime conclusion is included. |
| White-noise LC4 curve, all visibilities in [0,1] | Lean | `NoisyLC4.exact_curve`, `exact_curve_stochastic`, `sigma_ninety_percent`, `sigma_ninetyfive_percent` in [NoisyLC4ForcedSignaling.lean](https://github.com/stevenwarejones/ontology-separation/blob/c10474e9b45cca2ea260eec2d9f348a12cf62677/OntologySeparation/Experiments/NoisyLC4ForcedSignaling.lean); explicit zero/threshold/LC4 models and convex interpolation. The threshold is the known S4 noise tolerance. |
| Theorem 4: finite 24-vertex inequality | Python exact verifier | Part (1) of `reproduce_theorem4.py`: exhaustive over local and PR vertices and CHSH sign patterns. Although stored as NumPy floats, all inputs and intermediate sums in this finite part are exactly representable small dyadic numbers; no LP is used here. Parts (2)–(3) are floating-point numerical checks, not exact verifiers. The extension from vertices and the theorem's construction remain analytic. |
| Observation 1: exhibited certificate, row slack and reduced-program optimum | Python exact verifier | `verify_Sigma.py`, exact arithmetic in Q(sqrt2); does not mechanize the universal complementary-slackness argument. |
| Per-completion slope spectrum and multiplicities | numerical only | `reproduce_core.py (FULL)`; distinct from the Lean proof that the minimum possible slope is 8. |
| Collectibility and the explicit four-site geometries | analytic (not machine-checked) | Light-cone arguments in Section 6; not formalized by the invisible-model certificates. |

Section 6's experiment is a **proposal, not verified**. It is not an additional
theorem: confidence construction, the continuous candidate region, complete
event budgets, measurement-setting randomization, and memory-valid inference
remain open. A proof of Theorem 7 alone does not make the architecture ready.

The pinned [Tests/Audit.lean](https://github.com/stevenwarejones/ontology-separation/blob/c10474e9b45cca2ea260eec2d9f348a12cf62677/Tests/Audit.lean) and
[docs/AXIOM_AUDIT.txt](https://github.com/stevenwarejones/ontology-separation/blob/c10474e9b45cca2ea260eec2d9f348a12cf62677/docs/AXIOM_AUDIT.txt) register the public proof roots.
See the README's **Formal verification in Lean** section for exact reproduction
commands and the allowed-axiom policy. The certificate inputs for the
Propositions 1–2 Lean port are this repository's
`bae83b865ea60b1fbc4b808e66dc9bbaf7da2d4b` revision; values were not altered.

## Per-claim numerical and certificate coverage inventory
| claim | script | status |
|---|---|---|
| Thm 1 certificate (K=8) | verify_K8.py | independent exact verifier |
| Prop 1: dual split (2 per direction) | verify_directional.py | exact; read off K8_certificate.json, so it presupposes verify_K8.py passes |
| Prop 1: one-sided attainment (dA,dD)=(s,0),(0,s) | verify_directional.py | exact primal certificates over Q(sqrt2); these are NEW proof obligations, not consequences of the Thm 1 dual |
| Prop 2: optimum attainable with all singles/pairs blind | verify_invisibility.py | exact; three models (A_only, D_only, balanced), each reproducing every ABD/ACD marginal exactly, with all one- and two-party recipient marginals non-signaling and the triple difference a pure parity shift |
| Lemma 2: only {B,C} and {B,C,D} can carry a signal | (analytic) + reproduce_extra.py (FULL) | proof is one line -- any recipient set omitting B or C is fixed by the reproduced ABD/ACD families, which are no-signaling; the FULL script confirms it by LP, maximising a SINGLE-OUTCOME probability difference (not TV, which is not linear): five subsets give exactly 0; {B,C} reaches 0.323 and {B,C,D} 0.213, certifying compatible models with TV >= 0.323 and 0.213 (a singleton is an event, so TV >= max_o |dp(o)|) |
| Prop 2: co-located layout has zero accessible signaling FOR EVERY compatible model | (analytic) | needs BOTH conditions to fail: C_full from J_c^+(K_D) subset J_c^+(K_A); C_BC from the early parties lying on segment BC, which makes the two "closer to B/C than to the sender" half-spaces disjoint. The second is an explicit geometric HYPOTHESIS, not a consequence of co-location alone |
| Prop 2: separated 4-site example is collectible both ways | (analytic, exact rationals) | margins 1/20 and 3/5 at c=1, v=4; arithmetic stated in the text |
| Prop 2: four-site restoration (12 km line) | (analytic) | early parties outboard at +/-6 km, blind pair inboard at +/-5 km, 100 ns stagger; all pairs c-spacelike, both records collectible with ~3 us margin. AN EXAMPLE AT v=1e4 c (where the inclusions hold with 12-300x margin), valid down to kappa > 400.28 in the laboratory frame and kappa > 863.35 under the worst aligned boost -- it does NOT cover the whole c < v <= 1e4 c region; see the open item below. Collectibility is a LIGHT-cone condition, hence Lorentz invariant, so it adds no burden to the delay cover and its margin is ~70x the +/-42 ns delay span |
| Obs: support of the optimal dual | verify_Sigma.py | exact; all nonzero dual entries lie in comparisons {5,7,12,14} (A-changes at w=1 with z=1, D-changes at x=1 with z=0, both y), including all four of their total-variation rows, and all 36 nonzero entries equal 1/8 (the inequality dual is 1/8 times a 0/1 vector). The restriction holds for EVERY optimal dual, by optimality rather than feasibility: the twelve non-core TV rows are slack by exactly (sqrt2-1)/2 in the certified primal (checked exactly), so complementary slackness zeroes their multipliers, and the slack-column coupling lam[r1]+lam[r2] <= lam[TV] zeroes the associated absolute-value multipliers. A merely FEASIBLE dual can carry non-core weight (all equality multipliers zero, 1/8 on one non-core TV row, objective 0). The reduced statement is two-sided and both sides come from objects already certified here: dropping the twelve other TV constraints enlarges the primal feasible set (so the exhibited model still applies) and the dual carries no weight on them (so it stays dual feasible), giving Sigma over the four comparisons = (sqrt2-1)/4 exactly. This is complementary slackness at Q_LC4, NOT a bound: deleting constraints can only DECREASE Sigma, so an UPPER bound proved for the four-comparison restriction does not transfer to the full program (lower bounds do transfer). Whether the support has this shape at other behaviors is not established here |
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
- The four-site restoration is an example AT v=1e4 c, valid down to v ~ 8.6e2 c
  (kappa > 400.28 in the laboratory frame, kappa > 863.35 under the worst aligned
  boost). It does NOT inherit Theorem 7's coverage of the whole c < v <= 1e4 c region:
  collectibility forces tau < L/c while the v-inclusion forces tau > L/(kappa c), so a
  fixed early geometry covers only kappa > L/(c tau). Programming the early stagger to
  cover the full range is NOT worked out here.
- Proposition 2's universal (every-model) branch requires the early parties to lie between
  the two blind stations. Move them off that segment and {B,C} can become collectible, at
  which point only the weaker existential statement survives.
- Proposition 2 concerns a SINGLE-TRIAL, single-sender record-collection task with the
  other settings fixed in advance, no ancillas, no adaptivity and no repetition. It does
  not show that a finite-v theory evades the Bancal theorem, which asserts a channel in
  SOME suitable arrangement; it shows that the arrangement proposed here is not one of
  them, and that the statistical exclusion does not depend on that reading.
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
