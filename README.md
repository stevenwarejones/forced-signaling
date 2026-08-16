# The signaling cost of finite-speed hidden influences

Manuscript, exact-arithmetic certificates, and independent verifiers for the result
that finite-speed hidden-influence models explaining the four-qubit linear-cluster
correlations of Li, Hu, Deng and Scarani ([arXiv:2608.05271](https://arxiv.org/abs/2608.05271))
must exhibit **exactly (√2−1)/4 ≈ 0.1036** of operational signaling.

> ### Status: research memorandum, not a submission
>
> **This work was produced by AI systems — Claude (Anthropic) and Codex (OpenAI) — in an
> iterative, human-directed session with adversarial cross-checking between them. The
> listed author is not a physicist. No claim in it has been verified by a human domain
> expert.** Its two central results carry exact-arithmetic certificates with standalone
> verifiers (below), and every numerical claim has a reproduction script; that is a
> different and weaker thing than expert review. The manuscript states this on its title
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
- An experimental architecture: three sites, one ~10 km baseline, ~16 programmed delays
  covering every preferred frame with |β| ≤ 1.34×10⁻³ and every hidden speed *v* ≤ 10⁴*c*.

## Verify the central claims in about two seconds

Both verifiers rebuild the linear programs **from first principles** — integer and exact
symbolic arithmetic — and validate the stored certificates independently of the code that
generated them. The audit surface is ~180 lines of standalone Python, not the whole
package.

```bash
python3 paper/verify_K8.py      # Theorem 1  (~0.1 s)
python3 paper/verify_Sigma.py   # Theorem 2  (~2 s)
```

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
python3 threadB/reproduce_extra.py         # 4 setting supersets, random+random; FULL=1 adds LC5, mixed-flavor
python3 threadB/reproduce_theorem4.py      # 24-vertex check, 400 constrained distances, perturbation sweep
(cd paper && python3 figures.py)           # regenerates the three figures
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
| `paper/MANIFEST.md` | reproduction commands, per-claim coverage inventory, certificate schema, open items |
| `paper/figures.py`, `paper/fig_*.pdf` | figure source and output |
| `threadB/reproduce_*.py` | reproduction drivers for the numerical record |
| `threadB/adversary.py` | the Σ<sub>HIC</sub> linear program used by the adversarial constructions |
| `threadB/chained.py`, `qutrit.py`, `barnea.py` | scenario-specific libraries |
| `hashes.txt` | sha256 over the manuscript package |

## Integrity

```bash
sha256sum -c hashes.txt
```

`hashes.txt` covers the manuscript package — manuscript, scripts, figures, certificates
and manifest. It does not cover this README, the license, or other repository scaffolding.
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
  five remaining design obligations are listed explicitly rather than assumed away.

## License

Code (`*.py`) and certificates (`*.json`) are released under the MIT License — see
[`LICENSE`](LICENSE). The manuscript (`paper/main.tex`, `paper/main.pdf`) and the figures
are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Citation

If you use the certificates or the Σ<sub>HIC</sub> linear program, please cite the manuscript; see
[`CITATION.cff`](CITATION.cff). Note that it is an unrefereed research memorandum whose
claims have not been checked by a human expert, and it should be cited as such.

## Contact

Steven W. Jones — stevenwarejones@gmail.com

Corrections, refutations and pointers to prior art are all welcome, and a refutation is
worth more to me than agreement.
