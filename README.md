# The signaling cost of finite-speed hidden influences

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21962607.svg)](https://doi.org/10.5281/zenodo.21962607)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Manuscript, exact-arithmetic certificates, and independent verifiers for the result
that any finite-speed hidden-influence model reproducing the observable no-blind-pair
marginals of the four-qubit linear-cluster correlations of Li, Hu, Deng and Scarani
([arXiv:2608.05271](https://arxiv.org/abs/2608.05271)) must exhibit **at least
(√2−1)/4 ≈ 0.1036** of operational signaling — and that this minimum is attained
within the specified conditionally-local model class, so the bound is tight rather
than merely valid. Admissible models may signal more; what is fixed is the floor.
The data a model must reproduce are those marginals, not the entire quantum behavior.

> ### Status: research memorandum, not a submission
>
> **This work was produced by AI systems — Claude (Anthropic) and Codex (OpenAI) — in an
> iterative, human-directed session with adversarial cross-checking between them. The
> listed author is not a physicist. No claim in it has been verified by a human domain
> expert.** Its two central results carry exact-arithmetic certificates with standalone
> verifiers (below), and most numerical claims have a reproduction script — the
> per-claim inventory in `paper/MANIFEST.md` states exactly which do not, and every
> such gap is listed there rather than implied away. Even for the covered claims, that
> is a different and weaker thing than expert review. The manuscript states this on its title
> page and in a dedicated Statement on AI use, and is explicitly **not for submission in
> its current form**. It is circulated to invite exactly the scrutiny it has not yet had.
>
> If you find an error, that is the most useful possible outcome — please open an issue
> or email the address below.

## The result

Bancal, Pironio, Acín, Liang, Scarani and Gisin (Nat. Phys. **8**, 867 (2012)) proved that
any hidden-influence model propagating at finite speed *v > c* in a preferred frame must
become operationally signaling in order to reproduce quantum correlations in suitable
arrangements. That is a qualitative dilemma. This work makes it quantitative.

For a behavior *Q* in a "blind-pair" scenario, define the **forced signaling** Σ<sub>HIC</sub>(*Q*):
the minimum operational signaling, in total variation per setting change, that any
conditionally-local (finite-*v*) model must exhibit while reproducing every observable
no-blind-pair marginal of *Q*. Principal results:

- **Σ<sub>HIC</sub>(Q<sub>LC4</sub>) = (√2−1)/4 exactly** — primal and dual certificates over ℚ(√2).
- **S₄<sup>op</sup> ≤ 6 + 8·Δ<sub>sig</sub>**, an exact rational certificate — so an experiment can
  *measure* the signaling rather than assume no-signaling. The constant 8 is provably
  optimal, and is completion-dependent: the natural all-defaults estimation completion
  gives 16, a factor of two in sensitivity lost for free.
- A universal bound Σ<sub>HIC</sub> ≤ (√2−1)/2, a reduction of Σ<sub>HIC</sub> to steering assemblages, and a
  conjectured tight ceiling (√2−1)/4 with its full adversarial record.
- **Directional refinement:** the same certificate resolves by sender, giving
  S₄<sup>op</sup> ≤ 6 + 4δ<sub>A</sub> + 4δ<sub>D</sub> — tighter than 6 + 8Δ<sub>sig</sub> whenever the two
  directions differ — and either party alone can carry the whole minimum budget
  (√2−1)/2, so bounding the signaling out of one early party can never exclude the class.
- **What the experiment does and does not show:** the forced signaling can be carried
  entirely by the three-party joint record, with every one- and two-party recipient
  marginal exactly non-signaling. Whether it is *usable* then depends on the layout —
  and in the co-located three-site arrangement below *no* compatible model gives an
  accessible signal, because only recipient sets containing the whole blind pair can
  carry one and none of those is collectible there. The statistical exclusion is
  unaffected either way. A four-site variant
  (early parties outboard at ±6 km, blind pair inboard at ±5 km) restores it, with
  ~3 µs collectibility margins; collectibility is a light-cone condition and so is
  frame-independent, adding nothing to the delay-cover problem.
- An experimental architecture: three sites, one ~10 km baseline, ~16 programmed delays
  covering every preferred frame with |β| ≤ 1.34×10⁻³ and every hidden speed *v* ≤ 10⁴*c*.

## Verify the central claims in a couple of seconds

The verifiers rebuild the linear programs **from first principles** — integer and exact
symbolic arithmetic — and validate the stored certificates independently of the code that
generated them. Certificate entries are read by a restricted parser that accepts only
exact integer/rational literals and the token `sqrt(2)`: there is no floating-point atom,
no approximate number recognition, and no general expression evaluation anywhere in the
certificate path, so a supplied value cannot be silently replaced by a nearby simpler one.
All four **exit nonzero** if any check fails or the certificate is malformed, so they
can be used as automated gates. The audit surface is a few hundred lines of standalone
Python, not the whole package.

```bash
python3 paper/verify_K8.py           # Theorem 1
python3 paper/verify_Sigma.py        # Theorem 2
python3 paper/verify_directional.py  # Proposition 1 (directional refinement)
python3 paper/verify_invisibility.py # Proposition 2 (pairwise-invisible attaining models)
```

Each completes in a few seconds on an ordinary laptop (measured: 0.1–0.3 s, 0.9–1.5 s,
~3 s and ~1 s respectively); timings are machine-dependent. `verify_directional.py` reads
the Theorem 1 certificate, so run `verify_K8.py` first — or just run them all in order:

```bash
for v in K8 Sigma directional invisibility; do
  python3 paper/verify_$v.py || { echo "FAILED: verify_$v.py"; exit 1; }
done
```

(run that in a subshell or script so the `exit 1` is meaningful — a bare `|| echo`
would report success even when a verifier fails).

Expected output, verbatim:

```
dual feasible: True | min slack: 0
sum TV duals: 4 (need 4) | sum normalization duals: 6 (need 6)
VERDICT: CERTIFICATE VALID: S4^op <= 6 + 8*Delta_sig for all Delta_sig >= 0
```

```
PRIMAL: nonneg True, equalities True, inequalities True, value=(sqrt2-1)/4 True
DUAL: sign True, feasibility True, value=(sqrt2-1)/4 True
VERDICT: CERTIFICATES VALID: Sigma_HIC(Q_LC4) = (sqrt2-1)/4 exactly
```

**Dependencies:** python3 (3.9+), numpy, scipy (HiGHS via `linprog`), and sympy for
`verify_Sigma.py` only.

```bash
pip install numpy scipy sympy
```

## Reproducing the numerical record

All commands run from the repository root.

```bash
python3 threadB/reproduce_core.py          # QUICK: core numerical claims, minutes
FULL=1 python3 threadB/reproduce_core.py   # complete: 512-completion spectrum, parallel 2-copy, 500k MC
python3 threadB/reproduce_extra.py         # 4 setting supersets, random+random; FULL=1 adds LC5, mixed-flavor, subset pinning
python3 threadB/reproduce_theorem4.py      # 24-vertex check, 400 constrained distances, perturbation sweep
(cd paper && python3 figures.py)           # regenerates the three figures
```

Every LP is solved with the HiGHS primal and dual feasibility tolerances set explicitly to
`1e-9` (the library defaults are `1e-7`), a failed or nonfinite solve raises rather than
being read as a zero, and each run ends with a diagnostics line reporting the worst
equality residual, inequality violation and bound violation actually recomputed from the
returned solutions. Observed decimal agreement is reported as such; it is not a certified
error bound.

QUICK mode validates a documented subset — it does **not** exercise the full 512-completion
spectrum or its multiplicities (it does check the optimal completion explicitly). Claims
that need `FULL=1` are labelled in the manifest; do not report a FULL claim as reproduced
from a QUICK run.

```bash
python3 tests/regression_checks.py         # certificate-tamper and solver-failure gates
```

`paper/MANIFEST.md` carries a **per-claim coverage inventory** stating, for every numerical
claim in the manuscript, which script regenerates it and which claims are not yet covered.
Claims are labelled throughout as *proven*, *machine-verified*, *numerically established*,
or *conjectured*, and the distinction is meant literally.

## Contents

| Path | What it is |
|---|---|
| `paper/main.pdf`, `paper/main.tex` | the manuscript |
| `paper/verify_K8.py` | standalone exact verifier, Theorem 1 (integer/`Fraction` arithmetic) |
| `paper/verify_Sigma.py` | standalone exact verifier, Theorem 2 (sympy over ℚ(√2)) |
| `paper/K8_certificate.json` | rational dual vector, Theorem 1 |
| `paper/Sigma_LC4_certificates.json` | primal and dual certificates over ℚ(√2), Theorem 2 |
| `paper/verify_directional.py`, `paper/directional_certificates.json` | standalone exact verifier and the two one-sided certificates, Proposition 1 |
| `paper/verify_invisibility.py`, `paper/invisible_certificates.json` | standalone exact verifier and the three pairwise-invisible models, Proposition 2 |
| `paper/MANIFEST.md` | reproduction commands, per-claim coverage inventory, certificate schema, open items |
| `paper/figures.py`, `paper/fig_*.pdf` | figure source and output |
| `threadB/reproduce_*.py` | reproduction drivers for the numerical record |
| `threadB/adversary.py` | the Σ<sub>HIC</sub> linear program used by the adversarial constructions |
| `threadB/chained.py`, `qutrit.py`, `barnea.py` | scenario-specific libraries |
| `tests/regression_checks.py` | gates that the verifiers reject tampered certificates and that failed solves raise |
| `hashes.txt` | sha256 over the manuscript package |

## Integrity

```bash
sha256sum -c hashes.txt          # Linux
shasum -a 256 -c hashes.txt      # macOS (no sha256sum by default)
```

`hashes.txt` covers the manuscript package — manuscript, scripts, figures, certificates,
the manifest and the regression gates. It does not cover this README, the license, the
citation metadata, or other repository scaffolding.
The manuscript requires that any circulated copy of `main.pdf` be accompanied by this
package, without which its certificate claims cannot be assessed; a link to a tagged
release of this repository satisfies that.

## Known open items

Carried in full in the manuscript's Statement on AI use and in `paper/MANIFEST.md`. In brief:

- **No human domain expert has verified any claim.** AI cross-verification by a second
  system is not a substitute.
- Conjecture 1 (the tight ceiling) is supported by an adversarial numerical record only;
  the proof strategy — a steering-NPA moment hierarchy over jointly-extendible assemblages —
  is stated but not carried out.
- A transcription discrepancy in the Barnea *et al.* tripartite evaluation, whose numerical
  result is therefore withheld from the manuscript.
- The exact relation between input-maximized and input-averaged normalizations away from
  the extremal points.
- Literature priority of the parallel-flatness observation: a targeted sweep found nothing,
  but a systematic citation-graph pass is recommended before any posting.
- The statistical protocol of §6.3 is a sound outline, not a complete specification; the
  remaining design obligations — including confidence-bound construction, the continuous
  candidate region, full event budgets, setting randomization, and the non-i.i.d./memory
  structure induced by the timestamped protocol — are enumerated there explicitly rather
  than assumed away.
- The experiment is a proposed architecture with those obligations open. Nothing in this
  repository establishes that it is ready to run.

## License

Code (`*.py`) and certificates (`*.json`) are released under the MIT License — see
[`LICENSE`](LICENSE). The manuscript (`paper/main.tex`, `paper/main.pdf`) and the figures
are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Citation

Archived on Zenodo: [**10.5281/zenodo.21962607**](https://doi.org/10.5281/zenodo.21962607) (v1.10).
Machine-readable metadata is in [`CITATION.cff`](CITATION.cff); GitHub's "Cite this
repository" button reads it.

**The working tree is ahead of the last archived release.** The DOI above pins a specific
earlier snapshot; the draft version on the title page of `paper/main.pdf` is authoritative
for the manuscript in this tree. Cite the archived release you actually used, and check
the Zenodo record's "versions" list for the newest one.

Note that this is an unrefereed research memorandum whose claims have not been checked
by a human expert, and it should be cited as such.

## Contact

Steven W. Jones — stevenwarejones@gmail.com

Corrections, refutations and pointers to prior art are all welcome, and a refutation is
worth more to me than agreement.
