#!/usr/bin/env python3
"""Focused regression gates for the verification and reproduction machinery.

These do not re-derive any physics.  They check that the tooling FAILS when it
should -- the property that makes a passing run mean anything:

  A. verify_Sigma.py accepts the shipped certificate and exits 0.
  B. Certificate parsing is strictly exact: a perturbation of 1e-30 added to an
     entry survives parsing (it is not collapsed onto the nearby exact value),
     and the verifier rejects the perturbed certificate.
  C. An equality-breaking primal weight is rejected.
  D. A negative (sign-violating) dual entry is rejected.
  E. Malformed certificate data (float atom, bad index, wrong length) is rejected
     with a nonzero exit rather than coerced.
  F. verify_K8.py accepts the shipped certificate, and rejects a tampered one.
  G. Both verifiers keep their essential checks under `python -O` (asserts off).
  H. A failed LP solve raises LPFailure instead of being read as Sigma = 0, and a
     nonfinite objective is rejected.

All tamper tests run against temporary copies; the shipped certificates are never
modified.  Exit status 0 iff every gate holds.
Run: python3 tests/regression_checks.py
"""
import json, os, shutil, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(ROOT, 'paper')
sys.path.insert(0, os.path.join(ROOT, 'threadB'))

results = []


def check(name, ok, detail=''):
    results.append((name, ok, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    return ok


def run_verifier(script, workdir, opt=False):
    """Run a verifier with cwd=workdir so it reads that directory's certificate."""
    cmd = [sys.executable] + (['-O'] if opt else []) + [os.path.join(workdir, script)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def sandbox():
    """A temp copy of paper/ so tampering never touches the shipped files."""
    d = tempfile.mkdtemp(prefix='fs-regress-')
    for f in ('verify_Sigma.py', 'verify_K8.py', 'verify_directional.py',
              'verify_invisibility.py', 'Sigma_LC4_certificates.json',
              'K8_certificate.json', 'directional_certificates.json',
              'invisible_certificates.json'):
        shutil.copy(os.path.join(PAPER, f), d)
    return d


def tamper(workdir, certfile, mutate):
    path = os.path.join(workdir, certfile)
    cert = json.load(open(path))
    mutate(cert)
    json.dump(cert, open(path, 'w'))


print("== A. shipped certificates verify ==")
rc, out = run_verifier('verify_Sigma.py', PAPER)
check("verify_Sigma.py exits 0 on the shipped certificate", rc == 0, f"exit={rc}")
check("verify_Sigma.py reports CERTIFICATES VALID", 'CERTIFICATES VALID' in out)
rc, out = run_verifier('verify_K8.py', PAPER)
check("verify_K8.py exits 0 on the shipped certificate", rc == 0, f"exit={rc}")
check("verify_K8.py reports CERTIFICATE VALID", 'CERTIFICATE VALID' in out)

print("== B. parsing is strictly exact (no approximate number recognition) ==")
sys.path.insert(0, PAPER)
from verify_Sigma import parse_Q2, Q2, ParseError          # noqa: E402
from fractions import Fraction                              # noqa: E402
target = Q2(Fraction(-1, 4), Fraction(1, 4))
perturbed_literal = "(-1/4 + sqrt(2)/4) + 1/1000000000000000000000000000000"
parsed = parse_Q2(perturbed_literal)
check("a 1e-30 perturbation survives parsing (not snapped to the exact value)",
      parsed != target and (parsed - target).sign() > 0,
      f"difference = {(parsed - target).a}")

d = sandbox()
tamper(d, 'Sigma_LC4_certificates.json',
       lambda c: c['primal_t_nonzero'].__setitem__(
           '384', "(-1/4 + sqrt(2)/4) + 1/1000000000000000000000000000000"))
rc, out = run_verifier('verify_Sigma.py', d)
check("perturbed objective entry (index 384) is REJECTED", rc != 0, f"exit={rc}")
shutil.rmtree(d)

print("== C. equality-breaking primal weight is rejected ==")
d = sandbox()
tamper(d, 'Sigma_LC4_certificates.json',
       lambda c: c['primal_t_nonzero'].__setitem__('0', '1/4'))   # was 1/8
rc, out = run_verifier('verify_Sigma.py', d)
check("altered primal weight is REJECTED", rc != 0, f"exit={rc}")
shutil.rmtree(d)

print("== D. dual sign violation is rejected ==")
d = sandbox()
tamper(d, 'Sigma_LC4_certificates.json',
       lambda c: c['dual_lambda_nonzero'].__setitem__('85', '-1/8'))
rc, out = run_verifier('verify_Sigma.py', d)
check("negative dual entry is REJECTED", rc != 0, f"exit={rc}")
shutil.rmtree(d)

print("== D2. the certified dual's support is pinned to the four active comparisons ==")
# Observation obs:dualsupport rests on WHERE the dual weight sits, not only on its
# feasibility.  Moving weight onto a comparison outside {5,7,12,14} must be caught
# even when the moved entry is individually harmless (positive, same magnitude).
for label, row in [("an A-change comparison outside the core", 0 * 17 + 16),
                   ("a D-change comparison outside the core", 9 * 17 + 16)]:
    d = sandbox()
    tamper(d, 'Sigma_LC4_certificates.json',
           lambda c, r=row: c['dual_lambda_nonzero'].__setitem__(str(r), '1/8'))
    rc, out = run_verifier('verify_Sigma.py', d)
    check(f"dual weight added on {label} is REJECTED", rc != 0, f"exit={rc}")
    shutil.rmtree(d)
d = sandbox()
tamper(d, 'Sigma_LC4_certificates.json',
       lambda c: c['dual_lambda_nonzero'].pop('101'))
rc, out = run_verifier('verify_Sigma.py', d)
check("removing a core total-variation row is REJECTED", rc != 0, f"exit={rc}")
shutil.rmtree(d)
# The observation also states the COEFFICIENTS: every nonzero dual entry is exactly 1/8.
# Rescaling one of them keeps the support intact, so only the coefficient check catches it.
for label, row, val in [("a core total-variation row", '101', '1/4'),
                        ("a core non-TV row", '85', '1/16')]:
    d = sandbox()
    tamper(d, 'Sigma_LC4_certificates.json',
           lambda c, r=row, v=val: c['dual_lambda_nonzero'].__setitem__(r, v))
    rc, out = run_verifier('verify_Sigma.py', d)
    check(f"rescaling {label} away from 1/8 is REJECTED", rc != 0, f"exit={rc}")
    shutil.rmtree(d)

print("== E. malformed certificate data is rejected, not coerced ==")
for label, cf, mut in [
    ("float atom in a Sigma entry", 'Sigma_LC4_certificates.json',
     lambda c: c['primal_t_nonzero'].__setitem__('0', 0.125)),
    ("float-shaped string in a Sigma entry", 'Sigma_LC4_certificates.json',
     lambda c: c['primal_t_nonzero'].__setitem__('0', '0.125')),
    ("out-of-range Sigma index", 'Sigma_LC4_certificates.json',
     lambda c: c['primal_t_nonzero'].__setitem__('99999', '1/8')),
    ("wrong dual_mu length", 'Sigma_LC4_certificates.json',
     lambda c: c.__setitem__('dual_mu', c['dual_mu'][:-1])),
    ("float atom in the K8 certificate", 'K8_certificate.json',
     lambda c: c['nonzero_dual_entries'].__setitem__('80', 1.0)),
    ("out-of-range K8 index", 'K8_certificate.json',
     lambda c: c['nonzero_dual_entries'].__setitem__('99999', '1')),
]:
    d = sandbox(); tamper(d, cf, mut)
    script = 'verify_Sigma.py' if 'Sigma' in cf else 'verify_K8.py'
    rc, out = run_verifier(script, d)
    check(f"{label} is REJECTED", rc != 0, f"exit={rc}")
    shutil.rmtree(d)

print("== F. K8 certificate tampering is rejected ==")
d = sandbox()
tamper(d, 'K8_certificate.json',
       lambda c: c['normalization_duals'].__setitem__(1, '3'))   # breaks sum == 6
rc, out = run_verifier('verify_K8.py', d)
check("altered normalization dual is REJECTED", rc != 0, f"exit={rc}")
shutil.rmtree(d)

print("== F2. directional certificates (Proposition 1) ==")
rc, out = run_verifier('verify_directional.py', PAPER)
check("verify_directional.py exits 0 on the shipped certificates", rc == 0, f"exit={rc}")
check("verify_directional.py reports DIRECTIONAL CERTIFICATES VALID",
      'DIRECTIONAL CERTIFICATES VALID' in out)
for label, cf, mut in [
    # break the one-sided budget: delta_A no longer equals (sqrt2-1)/2
    ("altered delta_A budget", 'directional_certificates.json',
     lambda c: c['A_only'].__setitem__('384', ['-1/4', '1/4'])),
    # make the 'off' direction nonzero, so the model is no longer one-sided
    ("nonzero off-direction in D_only", 'directional_certificates.json',
     lambda c: c['D_only'].__setitem__('384', ['1/100', '0'])),
    # perturb a model weight, breaking the marginal equalities
    ("perturbed model weight", 'directional_certificates.json',
     lambda c: c['A_only'].__setitem__('0', ['1/4', '0'])),
    # malformed entry shape
    ("malformed entry shape", 'directional_certificates.json',
     lambda c: c['A_only'].__setitem__('0', '1/8')),
    # break the dual split: move an A-block TV multiplier so the block no longer sums to 2
    ("broken dual split (A block != 2)", 'K8_certificate.json',
     lambda c: c['nonzero_dual_entries'].__setitem__('261', '2')),
]:
    d = sandbox(); tamper(d, cf, mut)
    rc, out = run_verifier('verify_directional.py', d)
    check(f"{label} is REJECTED", rc != 0, f"exit={rc}")
    shutil.rmtree(d)

print("== F3. pairwise-invisibility certificates (Proposition 2) ==")
rc, out = run_verifier('verify_invisibility.py', PAPER)
check("verify_invisibility.py exits 0 on the shipped certificates", rc == 0, f"exit={rc}")
check("verify_invisibility.py reports INVISIBILITY CERTIFICATES VALID",
      'INVISIBILITY CERTIFICATES VALID' in out)
for label, mut in [
    # break a weight: the ABD/ACD marginals no longer match the quantum target exactly
    ("perturbed weight (epsilon no longer 0)",
     lambda c: c['balanced'].__setitem__('0', ['1/4', '0'])),
    # move weight between response functions so a PAIR marginal starts to signal
    ("pair marginal made to signal",
     lambda c: c['A_only'].__setitem__('1', ['1/16', '0'])),
    # zero out a weight so normalization fails
    ("broken normalization",
     lambda c: c['D_only'].__setitem__('0', ['0', '0'])),
    # negative weight
    ("negative weight", lambda c: c['balanced'].__setitem__('2', ['-1/8', '0'])),
    ("malformed entry shape", lambda c: c['balanced'].__setitem__('0', '1/8')),
]:
    d = sandbox(); tamper(d, 'invisible_certificates.json', mut)
    rc, out = run_verifier('verify_invisibility.py', d)
    check(f"{label} is REJECTED", rc != 0, f"exit={rc}")
    shutil.rmtree(d)

print("== G. essential checks survive `python -O` ==")
d = sandbox()
tamper(d, 'Sigma_LC4_certificates.json',
       lambda c: c['dual_lambda_nonzero'].__setitem__('85', '-1/8'))
rc, out = run_verifier('verify_Sigma.py', d, opt=True)
check("verify_Sigma.py -O still rejects a bad certificate", rc != 0, f"exit={rc}")
shutil.rmtree(d)
d = sandbox()
tamper(d, 'K8_certificate.json', lambda c: c['nonzero_dual_entries'].__setitem__('80', '-1'))
rc, out = run_verifier('verify_K8.py', d, opt=True)
check("verify_K8.py -O still rejects a negative dual", rc != 0, f"exit={rc}")
shutil.rmtree(d)

print("== H. solver failures raise instead of reading as Sigma = 0 ==")
import numpy as np                                            # noqa: E402
import adversary                                              # noqa: E402
from adversary import LPFailure                               # noqa: E402


class _FakeRes:
    def __init__(self, success=True, fun=0.0, x=None, status=4, message='injected'):
        self.success, self.fun, self.x = success, fun, x
        self.status, self.message = status, message


_real_linprog = adversary.linprog
try:
    adversary.linprog = lambda *a, **k: _FakeRes(success=False, fun=None, x=None)
    lp = adversary.SigmaLP(2, 2, 2, adversary.PA, adversary.PB, adversary.PC, adversary.PD)
    try:
        lp.solve(adversary.cluster4()); ok = False; detail = "no exception raised"
    except LPFailure as e:
        ok = True; detail = str(e)[:60]
    check("unsuccessful solve raises LPFailure", ok, detail)

    adversary.linprog = lambda *a, **k: _FakeRes(success=True, fun=float('nan'),
                                                 x=np.zeros(1))
    try:
        lp.solve(adversary.cluster4()); ok = False; detail = "no exception raised"
    except LPFailure as e:
        ok = True; detail = str(e)[:60]
    check("nonfinite objective raises LPFailure", ok, detail)

    adversary.linprog = lambda *a, **k: _FakeRes(success=True, fun=0.0,
                                                 x=np.full(1, np.inf))
    try:
        lp.solve(adversary.cluster4()); ok = False; detail = "no exception raised"
    except LPFailure as e:
        ok = True; detail = str(e)[:60]
    check("nonfinite solution vector raises LPFailure", ok, detail)
finally:
    adversary.linprog = _real_linprog

print("== I. no fallback-to-zero survives in the reproduction drivers ==")
import re                                                     # noqa: E402
offenders = []
for fn in sorted(os.listdir(os.path.join(ROOT, 'threadB'))):
    if not fn.endswith('.py'):
        continue
    src = open(os.path.join(ROOT, 'threadB', fn)).read()
    if re.search(r'\.solve\([^)]*\)\s*or\s*0', src):
        offenders.append(f"{fn}: `.solve(...) or 0`")
    if re.search(r'res\.fun\s+if\s+res\.success\s+else\s+None', src):
        offenders.append(f"{fn}: `res.fun if res.success else None`")
check("no `or 0` / silent-None solver fallbacks remain", not offenders, '; '.join(offenders))

failed = [n for n, ok, _ in results if not ok]
print()
print(f"{len(results) - len(failed)}/{len(results)} gates passed")
if failed:
    print("FAILED gates: " + "; ".join(failed))
sys.exit(1 if failed else 0)
