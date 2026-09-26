# Shot budget under i.i.d. sampling

This analysis applies the proved model-class inequalities to finite samples. It
is not a Lean proof of concentration and is not a loophole-free test. All 16
settings (x,y,z,w) receive n independent runs, N=16n in total, with a fixed law
within each setting. The six correlator terms use five distinct settings:
ABD and CD share (1,0,0,0). Simulation uses their shared counts and reuses the
same counts for the recipient marginals; these statistics are not assumed
independent.

`shot_budget.py --check` recomputes all 256 Born entries exactly in Q(sqrt(2)),
compares them with the pinned Lean table in `data/lc4_born_lean.json`, and checks
all entries against `threadB/adversary.py`. It verifies normalization,
nonnegativity, all four no-signaling directions and S4=4+2sqrt(2) exactly.
The fixture is a transcription of the kernel-checked `fullProb_eq_table`
certificate, whose revision is recorded in the JSON. White-noise probabilities
are evaluated with 80 decimal digits before conversion to NumPy doubles.
Floating-point package agreement is reported with its numerical tolerance.

## Confidence rule: rigorous under the stated i.i.d. assumptions

The two inequalities used, with natural logarithms, are

    Pr(|Ehat−E| ≥ t) ≤ 2 exp(−n t²/2),
    Pr(||phat−p||₁ ≥ ε) ≤ (2^k−2) exp(−n ε²/2).

The first is Hoeffding (1963), Theorem 2, equation (2.6), for variables in
[−1,1], applied to both tails. The second follows from Weissman et al. (2003),
Theorem 2.1, equation (8), and φ(πP)≥2 (Proposition 4.3; stated after equation
(10)). Their norm is L1, twice TV. For k=8 the prefactor is 254, not 14.
These forms were checked against the original papers:

- W. Hoeffding, *Probability Inequalities for Sums of Bounded Random Variables*,
  JASA 58(301), 13–30 (1963).
  [Original scan](https://www.cs.rpi.edu/academics/courses/spring06/random/hoefding.pdf),
  [DOI](https://doi.org/10.1080/01621459.1963.10500830).
- T. Weissman, E. Ordentlich, G. Seroussi, S. Verdú, M. J. Weinberger,
  *Inequalities for the L1 Deviation of the Empirical Distribution*,
  HPL-2003-97 (R.1), June 2003.
  [Original report mirror](https://shiftleft.com/mirrors/www.hpl.hp.com/techreports/2003/HPL-2003-97R1.pdf).

Allocate α/2 to six correlator events (α/12 each) and α/2 to 32 L1 events
(α/64 each): one BCD and one ABC marginal at every setting. The union bound
requires no independence between these 38 events. Define

    tα = sqrt(2 log(24/α)/n),
    εα = sqrt(2 log(16256/α)/n),
    L_S = Shat−8tα,
    U_A = min(1,deltahat_A+εα),
    U_D = min(1,deltahat_D+εα),
    U_Δ = max(U_A,U_D).

All confidence bounds hold jointly with probability at least 1−α. Each TV
comparison uses half of each endpoint's L1 deviation. The rejection rules are

    scalar:      L_S > 6+8 U_Δ,
    directional: L_S > 6+4 U_A+4 U_D.

Both have size at most α. The directional rejection region contains the scalar
region for every dataset with these same bounds. Choosing the directional test
therefore requires no additional multiplicity penalty for these two nested
rules. This does not justify optional stopping on experimental data or trying
unrelated tests without further error control.

## A: guaranteed power and sufficient budget

Apply the same 38-event bound with failure β=.05 to the quantum sampling law.
On that event, Shat≥S−8tβ and both empirical directional TVs are at most εβ,
although their true values are zero. Hence both tests reject with probability
at least 1−β whenever

    gap := p(4+2sqrt(2))−6 > 8(tα+tβ+εα+εβ).

This includes the positive empirical-TV bias and fluctuations. It does not
substitute deltahat=0. The α event is not union-bounded with β for power:
α defines the rule's deterministic radii; the β event alone ensures rejection
under the alternative.

Let C=8[sqrt(2log(24/α))+sqrt(2log(24/β))+
         sqrt(2log(16256/α))+sqrt(2log(16256/β))]. For gap>0,

    n_A = floor((C/gap)²)+1,    N_A=16 n_A.

The script checks the strict integer boundary. This is the smallest allocation
certified by this conservative sufficient condition, not the true minimum of
the test or an information-theoretic lower bound. Scalar and directional
worst-case budgets coincide. Since gap=(4+2sqrt(2))(p−p*), they diverge as
(p−p*)⁻². At or below p* this analysis certifies no finite sufficient budget.
Displayed integers numerically evaluate the rigorous formula; the arithmetic
is not Lean-certified.

## B: Monte Carlo estimated crossing

Run `python3 threadB/shot_budget.py --study` with `requirements-ci.txt`.
The study samples multinomial full-output counts at all 16 settings and applies
the rules above without modifying their confidence radii. Each search point
uses 5,000 trials with seed 260926, p, α, n and a stream number. Bisection locates
an empirical .95 crossing to 0.5% in n; because each point uses fresh draws,
empirical power need not be monotone. This is an estimated crossing, not a
certified minimum.

Each chosen n receives an independent, fixed 20,000-trial validation stream.
Its power estimate and exact Clopper–Pearson 95% binomial interval are reported.
These intervals describe pointwise simulation uncertainty, not simultaneous
coverage over the table. Validation power can fall slightly below .95; its
interval and the recorded search bracket make this uncertainty explicit.
Full points, successes, seeds and versions are in `data/shot_budget_results.json`;
[SHOT_BUDGET_RESULTS.md](SHOT_BUDGET_RESULTS.md) is the rendered table.

`--check` replays a stored point and checks the analytic budgets, monotonicity,
near-threshold divergence, simulated thresholds below the conservative budgets,
and the Lean-proved Sigma curve. Full Monte Carlo is excluded from CI.

## Scope

Not covered: memory effects across runs, detection / fair-sampling and
postselection, settings randomization, unknown-frame timing (the delay cover),
or systematic errors and calibration. The underlying geometry and model-class
premises still apply. A memory-robust test-martingale or prediction-based-ratio
test is future work. The baseline allocation is uniform; optional nonuniform
allocations must use setting-specific radii and preserve the same error budget.

## Optional nonuniform allocation

The same significance and power error splits apply with n replaced by each
setting's count n_s. Put c_S=sqrt(2log(24/α))+sqrt(2log(24/β)) and
c_T=sqrt(2log(16256/α))+sqrt(2log(16256/β)). Let w_s be the sum of absolute
score coefficients using setting s. A sufficient condition for both tests is

    gap > c_S sum_s w_s/sqrt(n_s)
          + 4 c_T max_(i,j in the 16 switch pairs) [1/sqrt(n_i)+1/sqrt(n_j)].

SLSQP numerically minimizes total counts using reciprocal-square-root
coordinates and a variable bounding every switch-pair sum. Every count is then
rounded upward, and the strict sufficient condition is checked again. Thus
feasibility is checked independently of optimizer success; global optimality
is not asserted. The optional table reports these sufficient budgets with
all 16 positive allocations in the JSON. It does not claim nonuniform Monte
Carlo power estimates. `--study` regenerates them; `--check` reproduces them.
