#!/usr/bin/env python3
"""Standalone verifier for the directional refinement (Proposition 1).

Two claims are checked, both in exact Q(sqrt2) arithmetic:

(a) DUAL SPLIT.  In K8_certificate.json the total-variation multipliers sum to 2
    over the eight A-switch contexts (rows 256-263 in verify_K8.py's row order)
    and to 2 over the eight D-switch contexts (rows 264-271).  Since each TV row's
    right-hand side is twice that direction's budget, the certificate already
    proves  S4^op <= 6 + 4*delta_A + 4*delta_D,  which is strictly tighter than
    S4^op <= 6 + 8*max(delta_A,delta_D) whenever the two directions differ.
    This check presupposes that verify_K8.py itself passes; run it first.

(b) ONE-SIDED TIGHTNESS.  Two primal certificates are validated against a
    from-first-principles reconstruction of the LP with the single Delta column
    split into delta_A (column 384) and delta_D (column 385): conditionally-local
    models reproducing every no-blind-pair marginal of Q_LC4 exactly, with
    (delta_A,delta_D) = ((sqrt2-1)/2, 0) and (0, (sqrt2-1)/2).  Together with (a)
    these give  min(delta_A + delta_D) = (sqrt2-1)/2 = 2*Sigma_HIC,  attained with
    either direction carrying the whole budget.

The exact-arithmetic model (Q2) and the restricted certificate parser are imported
from verify_Sigma.py rather than duplicated, so the package has ONE audited exact
number implementation; the LP reconstruction here is independent of it and of the
generating code.  No floating point is used anywhere in the checks.

Exit status: 0 if every check passes, 1 otherwise (including malformed input).
Run: python3 verify_directional.py"""
import json, os, sys
import numpy as np, sympy as sp
from fractions import Fraction
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from verify_Sigma import Q2, ZERO, ONE, parse_Q2, sym_to_Q2, ParseError, R2   # noqa: E402


def fail(msg):
    print(f"VERIFICATION ERROR: {msg}")
    print("VERDICT: FAILED")
    sys.exit(1)


def load(name):
    try:
        return json.load(open(os.path.join(HERE, name)))
    except (OSError, ValueError) as e:
        fail(f"cannot read {name}: {e}")


# ---------------------------------------------------------------- (a) dual split
k8 = load('K8_certificate.json')
if 'nonzero_dual_entries' not in k8 or not isinstance(k8['nonzero_dual_entries'], dict):
    fail("K8_certificate.json: missing or malformed 'nonzero_dual_entries'")
A_BLOCK, D_BLOCK = range(256, 264), range(264, 272)
sumA = sumD = Fraction(0)
for k, v in k8['nonzero_dual_entries'].items():
    if not (isinstance(k, str) and k.lstrip('-').isdigit()):
        fail(f"K8_certificate.json: index key {k!r} is not an integer literal")
    j = int(k)
    q = parse_Q2(v)
    if q.b != 0:
        fail(f"K8_certificate.json[{j}]: expected a rational, got {v!r}")
    if j in A_BLOCK: sumA += q.a
    elif j in D_BLOCK: sumD += q.a
print(f"DUAL SPLIT: A-switch TV multipliers sum to {sumA} (need 2), "
      f"D-switch to {sumD} (need 2)")
split_ok = (sumA == 2 and sumD == 2)
if split_ok:
    print("  => S4^op <= 6 + 4*delta_A + 4*delta_D  (given verify_K8.py passes)")

# ------------------------------------------- exact LC4 marginals (Li et al. settings)
def skron(*Ms):
    out = sp.Matrix([[1]])
    for M in Ms:
        out = sp.Matrix(np.kron(np.array(out.tolist(), dtype=object),
                                np.array(M.tolist(), dtype=object)))
    return out


X = sp.Matrix([[0, 1], [1, 0]]); Z = sp.Matrix([[1, 0], [0, -1]]); I = sp.eye(2)
plus = sp.Matrix([1, 1]) / R2
psi = plus
for _ in range(3):
    psi = sp.Matrix(np.kron(np.array(psi.tolist(), dtype=object),
                            np.array(plus.tolist(), dtype=object)))


def CZ(n, i, j):
    M = sp.zeros(2 ** n, 2 ** n)
    for s in range(2 ** n):
        M[s, s] = -1 if ((s >> (n - 1 - i)) & 1 and (s >> (n - 1 - j)) & 1) else 1
    return M


psi = CZ(4, 0, 1) * CZ(4, 1, 2) * CZ(4, 2, 3) * psi
pp = lambda O: {0: (I + O) / 2, 1: (I - O) / 2}
PA = [pp(X), pp(Z)]; PB = [pp((Z + X) / R2), pp((Z - X) / R2)]
PC = [pp(Z), pp(X)]; PD = [pp(X), pp(Z)]
qp = lambda x, y, z, w, a, b, c, d: sym_to_Q2(
    (psi.T * skron(PA[x][a], PB[y][b], PC[z][c], PD[w][d]) * psi)[0, 0])

# ------------------- LP reconstruction with the Delta column split in two (integers)
idx = {}; k = 0
for x, w, a, d, fb, fg in product(range(2), range(2), range(2), range(2), range(4), range(4)):
    idx[(x, w, a, d, fb, fg)] = k; k += 1
N = k; fB = lambda fn, y: (fn >> y) & 1; fC = lambda fn, z: (fn >> z) & 1

eq = []; rhs = []
for x, w in product(range(2), range(2)):
    r = np.zeros(N, dtype=int)
    for a, d, fb, fg in product(range(2), range(2), range(4), range(4)):
        r[idx[(x, w, a, d, fb, fg)]] = 1
    eq.append(r); rhs.append(ONE)
for x, y, w, a, b, d in product(range(2), repeat=6):
    r = np.zeros(N, dtype=int)
    for fb in range(4):
        if fB(fb, y) != b: continue
        for fg in range(4): r[idx[(x, w, a, d, fb, fg)]] += 1
    eq.append(r); rhs.append(qp(x, y, 0, w, a, b, 0, d) + qp(x, y, 0, w, a, b, 1, d))
for x, z, w, a, c, d in product(range(2), repeat=6):
    r = np.zeros(N, dtype=int)
    for fg in range(4):
        if fC(fg, z) != c: continue
        for fb in range(4): r[idx[(x, w, a, d, fb, fg)]] += 1
    eq.append(r); rhs.append(qp(x, 0, z, w, a, 0, c, d) + qp(x, 0, z, w, a, 1, c, d))

contexts = [('x',) + t for t in product(range(2), repeat=3)] + \
           [('w',) + t for t in product(range(2), repeat=3)]
total = N; cs = {}
for ci, _ in enumerate(contexts):
    cs[ci] = list(range(total, total + 8)); total += 8
dA = total; dD = total + 1; total += 2          # the split: 384 = delta_A, 385 = delta_D


def mcols(x, y, z, w, a, b, c, d):
    return [idx[(x, w, a, d, fb, fg)] for fb in range(4) if fB(fb, y) == b
            for fg in range(4) if fC(fg, z) == c]


Aub = []
for ci, ctx in enumerate(contexts):
    for oi, out in enumerate(product(range(2), repeat=3)):
        uv = cs[ci][oi]; r1 = np.zeros(total, dtype=int); r2 = np.zeros(total, dtype=int)
        if ctx[0] == 'x':
            _, w, y, z = ctx; b, c, d = out
            for a in range(2):
                for col in mcols(0, y, z, w, a, b, c, d): r1[col] += 1; r2[col] -= 1
                for col in mcols(1, y, z, w, a, b, c, d): r1[col] -= 1; r2[col] += 1
        else:
            _, x, y, z = ctx; a, b, c = out
            for d in range(2):
                for col in mcols(x, y, z, 0, a, b, c, d): r1[col] += 1; r2[col] -= 1
                for col in mcols(x, y, z, 1, a, b, c, d): r1[col] -= 1; r2[col] += 1
        r1[uv] -= 1; r2[uv] -= 1; Aub += [r1, r2]
    r = np.zeros(total, dtype=int)
    for uv in cs[ci]: r[uv] = 1
    r[dA if ctx[0] == 'x' else dD] = -2       # A-switch contexts charge delta_A, D-switch delta_D
    Aub.append(r)
Aub = np.array(Aub)
if Aub.shape[0] != 272 or len(eq) != 132 or total != 386:
    fail(f"reconstruction shape mismatch: {Aub.shape[0]} inequality rows (need 272), "
         f"{len(eq)} equality rows (need 132), {total} columns (need 386)")

# ------------------------------------------------- (b) the two one-sided certificates
cert = load('directional_certificates.json')
target_sum = Q2(Fraction(-1, 2), Fraction(1, 2))          # (sqrt(2)-1)/2
results = {}
for name, on, off in (('A_only', dA, dD), ('D_only', dD, dA)):
    if name not in cert or not isinstance(cert[name], dict):
        fail(f"directional_certificates.json: missing or malformed {name!r}")
    v = [ZERO] * total
    for i, pair in cert[name].items():
        if not (isinstance(i, str) and i.lstrip('-').isdigit()):
            fail(f"{name}: index key {i!r} is not an integer literal")
        j = int(i)
        if not 0 <= j < total:
            fail(f"{name}: index {j} out of range [0,{total})")
        if not (isinstance(pair, list) and len(pair) == 2):
            fail(f"{name}[{j}]: expected a [rational, rational] pair, got {pair!r}")
        try:
            r_, t_ = parse_Q2(pair[0]), parse_Q2(pair[1])
        except ParseError as e:
            fail(f"{name}[{j}]: {e}")
        if r_.b != 0 or t_.b != 0:
            fail(f"{name}[{j}]: components must be rationals, got {pair!r}")
        v[j] = Q2(r_.a, t_.a)
    nonneg = all(z.sign() >= 0 for z in v)
    eqs = all((sum((Q2(int(eq[i][j])) * v[j] for j in np.nonzero(eq[i])[0]), ZERO)
               - rhs[i]).is_zero() for i in range(132))
    ineq = all(sum((Q2(int(Aub[i, j])) * v[j] for j in np.nonzero(Aub[i])[0]), ZERO).sign() <= 0
               for i in range(272))
    val = (v[dA] + v[dD] - target_sum).is_zero()
    zero = v[off].is_zero()
    results[name] = all([nonneg, eqs, ineq, val, zero])
    print(f"{name}: nonneg {nonneg}, marginals {eqs}, signaling rows {ineq}, "
          f"delta_A+delta_D=(sqrt2-1)/2 {val}, other direction exactly 0 {zero}")

ok = split_ok and all(results.values())
print("VERDICT:", "DIRECTIONAL CERTIFICATES VALID: S4^op <= 6 + 4*delta_A + 4*delta_D, "
      "and min(delta_A+delta_D) = (sqrt2-1)/2 is attained one-sidedly in either direction"
      if ok else "FAILED")
sys.exit(0 if ok else 1)
