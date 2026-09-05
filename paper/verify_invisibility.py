#!/usr/bin/env python3
"""Standalone verifier for Proposition 2 (accessible signaling).

Checks, in exact Q(sqrt2) arithmetic, that the forced signaling of Q_LC4 can be
carried entirely by the THREE-PARTY joint record, with every one- and two-party
recipient marginal exactly non-signaling.  For each stored model it verifies:

  * nonnegativity and per-setting normalization of the 256 response weights;
  * that every ABD and ACD marginal equals the quantum target EXACTLY
    (reconstructed symbolically here from the state, not read from a file);
  * that for either early sender, every recipient subset of size 1 or 2 has
    identical marginals under the two sender settings, at every fixed setting
    of the remaining parties and with no postselection;
  * that B and C signal nowhere at all;
  * that the surviving three-party difference is a pure parity shift, i.e. of
    the form lambda*(-1)^(r1+r2+r3) -- the only form compatible with all proper
    marginals agreeing (expand in the binary Fourier basis: normalization kills
    the constant, and the one- and two-body agreements kill every remaining
    coefficient except the three-body one);
  * that (delta_A, delta_D) takes its claimed value.

Consequence, with Proposition 1's lower bound: the minimum forced signaling is
unchanged by demanding pairwise invisibility, so a model can saturate it while
no proper subset of recipients sees anything.  Whether that signal is usable
therefore depends on whether the full complementary record can be assembled
outside the sender's future light cone -- a property of the spacetime layout,
not of the correlations.

Exit status: 0 if every check passes, 1 otherwise (including malformed input).
Run: python3 verify_invisibility.py"""
import json, os, sys
import numpy as np, sympy as sp
from fractions import Fraction
from itertools import product, combinations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from verify_Sigma import Q2, ZERO, ONE, parse_Q2, sym_to_Q2, ParseError, R2   # noqa: E402


def fail(msg):
    print(f"VERIFICATION ERROR: {msg}")
    print("VERDICT: FAILED")
    sys.exit(1)


# ---------------------------------------------- exact LC4 target (Li et al. settings)
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
OBS = [[pp(X), pp(Z)], [pp((Z + X) / R2), pp((Z - X) / R2)],
       [pp(Z), pp(X)], [pp(X), pp(Z)]]
Q = {}
for x, y, z, w in product(range(2), repeat=4):
    for a, b, c, d in product(range(2), repeat=4):
        Q[(x, y, z, w, a, b, c, d)] = sym_to_Q2(
            (psi.T * skron(OBS[0][x][a], OBS[1][y][b],
                           OBS[2][z][c], OBS[3][w][d]) * psi)[0, 0])

# ------------------------------------------------------------- model reconstruction
KEYS = list(product(range(2), range(2), range(2), range(2), range(4), range(4)))
fB = lambda fn, y: (fn >> y) & 1
fC = lambda fn, z: (fn >> z) & 1
OUT4 = list(product(range(2), repeat=4))
SET4 = list(product(range(2), repeat=4))


def build(weights):
    """P(a,b,c,d | x,y,z,w) from the 256 response weights, exactly."""
    P = {st: {o: ZERO for o in OUT4} for st in SET4}
    for st in SET4:
        x, y, z, w = st
        for (kx, kw, a, d, fb, fg), t in zip(KEYS, weights):
            if kx == x and kw == w:
                o = (a, fB(fb, y), fC(fg, z), d)
                P[st][o] = P[st][o] + t
    return P


def marg(P, st, sub):
    return tuple(sum((P[st][o] for o in OUT4 if tuple(o[i] for i in sub) == v), ZERO)
                 for v in product(range(2), repeat=len(sub)))


def tv(p, q):
    return Q2(Fraction(1, 2)) * sum((abs_q2(a - b) for a, b in zip(p, q)), ZERO)


def abs_q2(v):
    return v if v.sign() >= 0 else Q2(-v.a, -v.b)


def qmax(xs):
    best = ZERO
    for v in xs:
        if (v - best).sign() > 0: best = v
    return best


# ------------------------------------------------------------------ load and check
path = os.path.join(HERE, 'invisible_certificates.json')
try:
    certs = json.load(open(path))
except (OSError, ValueError) as e:
    fail(f"cannot read invisible_certificates.json: {e}")

s = Q2(Fraction(-1, 2), Fraction(1, 2))          # (sqrt2-1)/2
half_s = Q2(Fraction(-1, 4), Fraction(1, 4))     # (sqrt2-1)/4
EXPECTED = {'A_only': (s, ZERO), 'D_only': (ZERO, s), 'balanced': (half_s, half_s)}

all_ok = True
for name in ('A_only', 'D_only', 'balanced'):
    if name not in certs or not isinstance(certs[name], dict):
        fail(f"invisible_certificates.json: missing or malformed {name!r}")
    w = [ZERO] * 256
    for i, pair in certs[name].items():
        if not (isinstance(i, str) and i.lstrip('-').isdigit()):
            fail(f"{name}: index key {i!r} is not an integer literal")
        j = int(i)
        if not 0 <= j < 256:
            fail(f"{name}: index {j} out of range [0,256)")
        if not (isinstance(pair, list) and len(pair) == 2):
            fail(f"{name}[{j}]: expected a [rational, rational] pair, got {pair!r}")
        try:
            r_, t_ = parse_Q2(pair[0]), parse_Q2(pair[1])
        except ParseError as e:
            fail(f"{name}[{j}]: {e}")
        if r_.b != 0 or t_.b != 0:
            fail(f"{name}[{j}]: components must be rationals, got {pair!r}")
        w[j] = Q2(r_.a, t_.a)

    nonneg = all(v.sign() >= 0 for v in w)
    P = build(w)
    norm = all((sum((P[st][o] for o in OUT4), ZERO) - ONE).is_zero() for st in SET4)

    # exact reproduction of every ABD and ACD marginal (epsilon = 0)
    exact = True
    for x, y, w_, a, b, d in product(range(2), repeat=6):
        lhs = sum((P[(x, y, 0, w_)][o] for o in OUT4 if (o[0], o[1], o[3]) == (a, b, d)), ZERO)
        rhs = sum((Q[(x, y, 0, w_, a, b, c, d)] for c in range(2)), ZERO)
        if not (lhs - rhs).is_zero(): exact = False; break
    if exact:
        for x, z, w_, a, c, d in product(range(2), repeat=6):
            lhs = sum((P[(x, 0, z, w_)][o] for o in OUT4 if (o[0], o[2], o[3]) == (a, c, d)), ZERO)
            rhs = sum((Q[(x, 0, z, w_, a, b, c, d)] for b in range(2)), ZERO)
            if not (lhs - rhs).is_zero(): exact = False; break

    # invisibility to every single and every pair; parity form of the triple difference
    invisible = True; parity = True; deltas = []
    for sender in range(4):
        rec = tuple(i for i in range(4) if i != sender)
        triple = []
        for others in product(range(2), repeat=3):
            st = [0] * 4
            for i, v in zip(rec, others): st[i] = v
            st0 = tuple(st); st[sender] = 1; st1 = tuple(st)
            for size in (1, 2):
                for sub in combinations(rec, size):
                    if marg(P, st0, sub) != marg(P, st1, sub):
                        invisible = False
            p0, p1 = marg(P, st0, rec), marg(P, st1, rec)
            triple.append(tv(p0, p1))
            diff = [a - b for a, b in zip(p0, p1)]
            for k, o in enumerate(product(range(2), repeat=3)):
                sgn = 1 if sum(o) % 2 == 0 else -1
                if not (diff[k] - Q2(sgn) * diff[0]).is_zero(): parity = False
        deltas.append(qmax(triple))
    bc_silent = deltas[1].is_zero() and deltas[2].is_zero()
    want = EXPECTED[name]
    values = (deltas[0] - want[0]).is_zero() and (deltas[3] - want[1]).is_zero()

    ok = all([nonneg, norm, exact, invisible, parity, bc_silent, values])
    all_ok = all_ok and ok
    print(f"{name}: nonneg {nonneg}, normalized {norm}, marginals exact {exact}, "
          f"singles/pairs blind {invisible}, triple difference is pure parity {parity}, "
          f"B,C silent {bc_silent}, (delta_A,delta_D) as claimed {values}")

print("VERDICT:", "INVISIBILITY CERTIFICATES VALID: the forced signaling of Q_LC4 is "
      "attainable with every one- and two-party recipient marginal exactly non-signaling"
      if all_ok else "FAILED")
sys.exit(0 if all_ok else 1)
