#!/usr/bin/env python3
"""Standalone verifier for the Theorem 2 certificates (Sigma_LC4_certificates.json).
Rebuilds the exact LC4 quantum marginals and the LP constraint system from first
principles (sympy exact arithmetic over Q(sqrt2)); does NOT import the generating code.
Checks: primal feasibility (nonneg, 132 equalities, 272 inequalities) with objective
(sqrt2-1)/4, and dual feasibility with the same objective.

All certificate arithmetic is carried out in an exact model of Q(sqrt2) built on
fractions.Fraction (class Q2 below).  Certificate entries are read by a restricted
recursive-descent parser that accepts only integers, +, -, *, /, parentheses and the
literal sqrt(2): no floating-point atoms, no approximate number recognition, and no
general expression evaluation.  In particular nothing here can silently replace a
supplied value by a nearby simpler one.

Exit status: 0 if every check passes, 1 otherwise (including malformed input).
Run: python3 verify_Sigma.py"""
import json, os, sys
import numpy as np, sympy as sp
from fractions import Fraction
from itertools import product

R2 = sp.sqrt(2)


# ---------------------------------------------------------------- exact Q(sqrt2)
class Q2:
    """Exact element a + b*sqrt(2) of Q(sqrt2), a,b in Q. No floats anywhere."""
    __slots__ = ('a', 'b')

    def __init__(self, a=0, b=0):
        if isinstance(a, float) or isinstance(b, float):
            raise TypeError("Q2 refuses float input")
        self.a = Fraction(a); self.b = Fraction(b)

    def __add__(s, o): o = _q2(o); return Q2(s.a + o.a, s.b + o.b)
    __radd__ = __add__

    def __neg__(s): return Q2(-s.a, -s.b)

    def __sub__(s, o): return s + (-_q2(o))

    def __rsub__(s, o): return _q2(o) + (-s)

    def __mul__(s, o):
        o = _q2(o)
        return Q2(s.a * o.a + 2 * s.b * o.b, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = _q2(o)
        den = o.a * o.a - 2 * o.b * o.b          # norm; zero only if o == 0
        if den == 0:
            raise ZeroDivisionError("division by zero in Q(sqrt2)")
        return s * Q2(o.a / den, -o.b / den)

    def __eq__(s, o):
        o = _q2(o); return s.a == o.a and s.b == o.b

    def __hash__(s): return hash((s.a, s.b))

    def sign(s):
        """Exact sign of a + b*sqrt(2), by comparing a^2 with 2b^2. No floats."""
        a, b = s.a, s.b
        if b == 0: return (a > 0) - (a < 0)
        if a == 0: return (b > 0) - (b < 0)
        if a > 0 and b > 0: return 1
        if a < 0 and b < 0: return -1
        t = a * a - 2 * b * b                     # sqrt(2) irrational => t != 0 here
        if a > 0: return 1 if t > 0 else -1       # a > 0 > b
        return -1 if t > 0 else 1                 # a < 0 < b

    def is_zero(s): return s.a == 0 and s.b == 0

    def __repr__(s): return f"Q2({s.a}, {s.b})"


def _q2(v):
    if isinstance(v, Q2): return v
    if isinstance(v, float): raise TypeError("float in exact arithmetic")
    return Q2(v)


ZERO = Q2(0); ONE = Q2(1); SQRT2 = Q2(0, 1)


# ------------------------------------------- restricted exact parser for Q(sqrt2)
class ParseError(ValueError):
    pass


def parse_Q2(text):
    """Parse an exact Q(sqrt2) literal. Grammar (whitespace ignored):
         expr   := term (('+'|'-') term)*
         term   := factor (('*'|'/') factor)*
         factor := ('+'|'-') factor | atom
         atom   := INTEGER | 'sqrt(2)' | '(' expr ')'
    Deliberately has no float literals, no exponent notation, no symbols and no
    function calls other than the exact token sqrt(2)."""
    if not isinstance(text, str):
        raise ParseError(f"certificate entry is {type(text).__name__}, expected a string")
    s = text.replace(' ', '')
    if not s:
        raise ParseError("empty certificate entry")
    pos = 0

    def peek():
        return s[pos] if pos < len(s) else ''

    def expr():
        nonlocal pos
        v = term()
        while peek() in ('+', '-'):
            op = s[pos]; pos += 1
            v = v + term() if op == '+' else v - term()
        return v

    def term():
        nonlocal pos
        v = factor()
        while peek() in ('*', '/'):
            op = s[pos]; pos += 1
            v = v * factor() if op == '*' else v / factor()
        return v

    def factor():
        nonlocal pos
        if peek() == '-':
            pos += 1; return -factor()
        if peek() == '+':
            pos += 1; return factor()
        return atom()

    def atom():
        nonlocal pos
        if peek() == '(':
            pos += 1; v = expr()
            if peek() != ')':
                raise ParseError(f"unbalanced parenthesis in {text!r}")
            pos += 1; return v
        if s.startswith('sqrt(2)', pos):
            pos += 7; return SQRT2
        start = pos
        while pos < len(s) and s[pos].isdigit():
            pos += 1
        if pos == start:
            raise ParseError(f"unexpected character {peek()!r} in {text!r}")
        return Q2(int(s[start:pos]))

    v = expr()
    if pos != len(s):
        raise ParseError(f"trailing input {s[pos:]!r} in {text!r}")
    return v


def sym_to_Q2(e):
    """Convert an exact sympy number known to lie in Q(sqrt2) into a Q2, exactly.
    Rejects floats and anything outside the declared field. No nsimplify: the
    coefficients are read off the expanded expression and then re-checked."""
    e = sp.expand(sp.radsimp(sp.expand(e)))
    if e.atoms(sp.Float):
        raise ParseError(f"floating-point atom in reconstructed value {e}")
    if not e.is_number or not e.is_real:
        raise ParseError(f"reconstructed value {e} is not a real number")
    a = Fraction(0); b = Fraction(0)
    for base, coef in e.as_coefficients_dict().items():
        if not coef.is_Rational:
            raise ParseError(f"non-rational coefficient {coef} in {e}")
        c = Fraction(int(coef.p), int(coef.q))
        if base == sp.S.One:
            a += c
        elif base == R2:
            b += c
        else:
            raise ParseError(f"value {e} is outside Q(sqrt2): basis element {base}")
    out = Q2(a, b)
    # independent re-check that the extracted pair really equals the input
    if sp.simplify(e - (sp.Rational(out.a.numerator, out.a.denominator)
                        + sp.Rational(out.b.numerator, out.b.denominator) * R2)) != 0:
        raise ParseError(f"field extraction failed for {e}")
    return out


def fail(msg):
    print(f"VERIFICATION ERROR: {msg}")
    print("VERDICT: FAILED")
    sys.exit(1)


def main():
    # ------------------------------- exact quantum marginals of |LC4> (Li et al. settings)
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


    def qp(x, y, z, w, a, b, c, d):
        """Exact Born-rule probability, returned as an element of Q(sqrt2)."""
        return sym_to_Q2((psi.T * skron(PA[x][a], PB[y][b], PC[z][c], PD[w][d]) * psi)[0, 0])


    # ----------------------------------------------- LP reconstruction (integers only)
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
    delta = total; total += 1


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
        r[delta] = -2; Aub.append(r)
    Aub = np.array(Aub)   # 272 rows
    if Aub.shape[0] != 272 or len(eq) != 132:
        fail(f"reconstruction shape mismatch: {Aub.shape[0]} inequality rows "
             f"(expected 272), {len(eq)} equality rows (expected 132)")

    # ------------------------------------------------------- load and validate certificate
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Sigma_LC4_certificates.json')
    try:
        cert = json.load(open(path))
    except (OSError, ValueError) as e:
        fail(f"cannot read {os.path.basename(path)}: {e}")

    for key in ('primal_t_nonzero', 'dual_lambda_nonzero', 'dual_mu'):
        if key not in cert:
            fail(f"certificate is missing required key {key!r}")


    def sparse_vec(entries, length, label):
        """Exact sparse -> dense, with strict key and range validation."""
        if not isinstance(entries, dict):
            fail(f"{label}: expected an object of index -> value, got {type(entries).__name__}")
        out = [ZERO] * length
        for i, v in entries.items():
            if not (isinstance(i, str) and (i.isdigit() or (i[:1] == '-' and i[1:].isdigit()))):
                fail(f"{label}: index key {i!r} is not an integer literal")
            j = int(i)
            if not 0 <= j < length:
                fail(f"{label}: index {j} out of range [0,{length})")
            try:
                out[j] = parse_Q2(v)
            except ParseError as e:
                fail(f"{label}[{j}]: {e}")
        return out


    tF = sparse_vec(cert['primal_t_nonzero'], total, 'primal_t_nonzero')
    lam = sparse_vec(cert['dual_lambda_nonzero'], 272, 'dual_lambda_nonzero')
    if not isinstance(cert['dual_mu'], list) or len(cert['dual_mu']) != 132:
        fail(f"dual_mu: expected a list of 132 entries, got "
             f"{len(cert['dual_mu']) if isinstance(cert['dual_mu'], list) else type(cert['dual_mu']).__name__}")
    try:
        mu = [parse_Q2(v) for v in cert['dual_mu']]
    except ParseError as e:
        fail(f"dual_mu: {e}")
    target = Q2(Fraction(-1, 4), Fraction(1, 4))          # (sqrt(2)-1)/4
    if 'value' in cert:
        try:
            if parse_Q2(cert['value']) != target:
                fail(f"certificate 'value' field {cert['value']!r} is not (sqrt(2)-1)/4")
        except ParseError as e:
            fail(f"value: {e}")

    # ------------------------------------------------------------------------- checks
    okPn = all(e.sign() >= 0 for e in tF)
    okPe = all((sum((Q2(int(eq[i][j])) * tF[j] for j in np.nonzero(eq[i])[0]), ZERO)
                - rhs[i]).is_zero() for i in range(132))
    okPu = all(sum((Q2(int(Aub[i, j])) * tF[j] for j in np.nonzero(Aub[i])[0]), ZERO).sign() <= 0
               for i in range(272))
    okPv = (tF[delta] - target).is_zero()
    print(f"PRIMAL: nonneg {okPn}, equalities {okPe}, inequalities {okPu}, value=(sqrt2-1)/4 {okPv}")

    okDl = all(l.sign() >= 0 for l in lam)
    okDf = True
    for j in range(total):
        s = sum((Q2(-int(Aub[i, j])) * lam[i] for i in np.nonzero(Aub[:, j])[0]), ZERO)
        if j < N:
            for i in range(132):
                if eq[i][j]: s = s + Q2(int(eq[i][j])) * mu[i]
        if (s - (ONE if j == delta else ZERO)).sign() > 0:
            okDf = False; break
    okDv = (sum((mu[i] * rhs[i] for i in range(132)), ZERO) - target).is_zero()
    print(f"DUAL: sign {okDl}, feasibility {okDf}, value=(sqrt2-1)/4 {okDv}")

    ok = okPn and okPe and okPu and okPv and okDl and okDf and okDv
    print("VERDICT:", "CERTIFICATES VALID: Sigma_HIC(Q_LC4) = (sqrt2-1)/4 exactly" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
