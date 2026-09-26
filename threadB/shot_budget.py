#!/usr/bin/env python3
"""i.i.d. sufficient shot budgets and seeded power estimates. Not loophole-free.
Run --check in CI; --study regenerates the full Monte Carlo results.
"""
from __future__ import annotations
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import sys
import numpy as np
from scipy.stats import beta as binomial_beta

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'paper'))
from verify_Sigma import Q2
from adversary import cluster4, PA, PB, PC, PD, kron

SETTINGS = list(product(range(2), repeat=4))
OUTCOMES = SETTINGS.copy()
TERMS = [(1,(0,0,0,1),(0,1)), (1,(0,1,0,1),(0,1)),
         (1,(1,0,0,0),(0,1,3)), (-1,(1,1,0,0),(0,1,3)),
         (2,(1,0,0,0),(2,3)), (2,(0,0,1,1),(0,2,3))]
VISIBILITIES = [.90,.92,.95,.98,1.00]
ALPHAS = [.05,.01,.001]
PSTAR = 3-3*math.sqrt(2)/2
SQUANTUM = 4+2*math.sqrt(2)
SEED = 260926
BETA = .05
COEFFICIENTS = np.zeros((16,16))
for c,s,parties in TERMS:
    COEFFICIENTS[SETTINGS.index(s)] += [c*(-1)**sum(o[k] for k in parties) for o in OUTCOMES]

def require(condition,message):
    if not condition: raise AssertionError(message)

def effect(party,setting,outcome,i,j):
    sign=(-1)**outcome
    if party==1:
        return Q2(Fraction(1,2) if i==j else 0,
                  Fraction(sign*((-1)**i if i==j else (-1)**setting),4))
    is_z=(party in (0,3) and setting==1) or (party==2 and setting==0)
    if is_z: return Q2(int(i==j==outcome))
    return Q2(Fraction(1,2) if i==j else Fraction(sign,2))

@lru_cache(None)
def exact_born():
    amplitudes=[Fraction((-1)**(a*b+b*c+c*d),4) for a,b,c,d in OUTCOMES]
    ans=[]
    for s in SETTINGS:
        row=[]
        for o in OUTCOMES:
            value=Q2(0)
            for ii,i in enumerate(OUTCOMES):
                for jj,j in enumerate(OUTCOMES):
                    term=Q2(amplitudes[ii]*amplitudes[jj])
                    for party in range(4):
                        term*=effect(party,s[party],o[party],i[party],j[party])
                        if term.is_zero(): break
                    value+=term
            require(value.sign()>=0,'STOP: negative exact Born entry')
            row.append(value)
        require(sum(row,Q2(0))==1,'STOP: Born normalization mismatch')
        ans.append(row)
    return ans

@lru_cache(None)
def probabilities(p):
    if not math.isfinite(p) or not 0<=p<=1: raise ValueError('visibility must lie in [0,1]')
    with localcontext() as ctx:
        ctx.prec=80
        rt=Decimal(2).sqrt(); pd=Decimal(str(p))
        def dec(f): return Decimal(f.numerator)/Decimal(f.denominator)
        q=np.array([[float(pd*(dec(v.a)+dec(v.b)*rt)+(1-pd)/16) for v in row]
                    for row in exact_born()])
    require(np.all(q>=0) and np.max(abs(q.sum(axis=1)-1))<2e-15,'invalid sampling law')
    return q

def cross_check():
    exact=exact_born()
    fixture=json.loads((ROOT/'threadB/data/lc4_born_lean.json').read_text())
    pinned=[Q2(Fraction(a),Fraction(b)) for a,b in fixture['coefficients_q_sqrt2']]
    require([v for row in exact for v in row]==pinned,
            'STOP: Born table contradicts pinned Lean certificate')
    psi=cluster4(); q=probabilities(1); maxerr=0.
    for si,s in enumerate(SETTINGS):
        for oi,o in enumerate(OUTCOMES):
            value=float(np.real(psi.conj()@kron(PA[s[0]][o[0]],PB[s[1]][o[1]],
                       PC[s[2]][o[2]],PD[s[3]][o[3]])@psi))
            maxerr=max(maxerr,abs(value-q[si,oi]))
    require(maxerr<3e-14,'STOP: Born law contradicts threadB/adversary.py')
    for sender in range(4):
        for fixed in product(range(2),repeat=3):
            s0=list(fixed);s0.insert(sender,0);s1=s0.copy();s1[sender]=1
            for recipients in product(range(2),repeat=3):
                def marginal(s):
                    return sum((exact[SETTINGS.index(tuple(s))][oi] for oi,o in enumerate(OUTCOMES)
                        if tuple(o[k] for k in range(4) if k!=sender)==recipients),Q2(0))
                require(marginal(s0)==marginal(s1),'STOP: exact no-signaling mismatch')
    score=sum((Q2(int(COEFFICIENTS[s,o]))*exact[s][o] for s in range(16) for o in range(16)),Q2(0))
    require(score==Q2(4,2),'STOP: score contradicts Lean')
    return {'all_256_lean_entries_equal_exactly':True,'all_sender_marginals_equal_exactly':True,
            'package_max_absolute_error':maxerr,'score_exact':'4+2*sqrt(2)'}

def radii(n,failure):
    if n<=0 or not 0<failure<1: raise ValueError('n>0 and failure in (0,1) required')
    return math.sqrt(2*math.log(24/failure)/n),math.sqrt(2*math.log(16256/failure)/n)

def sufficient_n(p,alpha,beta=BETA):
    if not 0<=p<=1 or not 0<alpha<1 or not 0<beta<1: raise ValueError('invalid parameter')
    gap=p*SQUANTUM-6
    if gap<=0:return None
    threshold=(8*sum((*radii(1,alpha),*radii(1,beta)))/gap)**2
    n=math.floor(threshold)+1
    def ok(k):return gap>8*sum((*radii(k,alpha),*radii(k,beta)))
    while not ok(n):n+=1
    while n>1 and ok(n-1):n-=1
    return n

def statistics(counts,n):
    empirical=counts/n
    shat=np.einsum('tso,so->t',empirical,COEFFICIENTS)
    bcd=empirical.reshape(-1,16,2,8).sum(axis=2)
    abc=empirical.reshape(-1,16,8,2).sum(axis=3)
    da=(abs(bcd[:,:8]-bcd[:,8:]).sum(axis=2)/2).max(axis=1)
    dd=(abs(abc[:,0::2]-abc[:,1::2]).sum(axis=2)/2).max(axis=1)
    return shat,da,dd

def decisions(counts,n,alpha):
    s,a,d=statistics(counts,n);t,e=radii(n,alpha)
    ua=np.minimum(1,a+e);ud=np.minimum(1,d+e);ls=s-8*t
    scalar=ls>6+8*np.maximum(ua,ud);directional=ls>6+4*ua+4*ud
    require(np.all(~scalar|directional),'directional test must dominate pointwise')
    return scalar,directional

def cp_interval(k,n):
    return [0. if k==0 else float(binomial_beta.ppf(.025,k,n-k+1)),
            1. if k==n else float(binomial_beta.ppf(.975,k+1,n-k))]

def simulate(p,alpha,n,trials,stream=0):
    key=[SEED,round(p*1000000),round(alpha*1000000),int(n),stream]
    rng=np.random.default_rng(np.random.SeedSequence(key))
    counts=np.stack([rng.multinomial(n,row,size=trials) for row in probabilities(p)],axis=1)
    s,d=decisions(counts,n,alpha)
    return {'n':int(n),'N':16*int(n),'trials':trials,'stream':stream,
        'scalar_successes':int(s.sum()),'directional_successes':int(d.sum()),
        'scalar_power':float(s.mean()),'directional_power':float(d.mean()),
        'scalar_cp95':cp_interval(int(s.sum()),trials),'directional_cp95':cp_interval(int(d.sum()),trials)}

def study(p,alpha,trials=5000,validation_trials=20000):
    bound=sufficient_n(p,alpha);require(bound is not None,'no sufficient finite budget')
    observations={};estimates={}
    def at(n):
        if n not in observations:observations[n]=simulate(p,alpha,n,trials)
        return observations[n]
    for test in ('scalar','directional'):
        lo=1;hi=bound
        require(at(hi)[test+'_power']>=.95,'inspect: MC power below .95 at rigorous bound')
        while hi-lo>max(1,math.ceil(hi*.005)):
            mid=(lo+hi)//2
            if at(mid)[test+'_power']>=.95:hi=mid
            else:lo=mid
        validation=simulate(p,alpha,hi,validation_trials,stream=1)
        estimates[test]={'n':hi,'N':16*hi,'search_bracket_N':[16*lo,16*hi],
            'power':validation[test+'_power'],'cp95':validation[test+'_cp95'],'validation':validation}
        require(hi<=bound,'estimated crossing exceeds sufficient bound')
    return {'visibility':p,'alpha':alpha,'beta':BETA,'rigorous_scalar_N':16*bound,
            'rigorous_directional_N':16*bound,'monte_carlo':estimates,
            'search_points':[observations[n] for n in sorted(observations)]}

def check():
    verified=cross_check()
    require(sum(abs(c) for c,_,_ in TERMS)==8,'coefficient sum')
    require(len({s for _,s,_ in TERMS})==5,'six terms share five settings')
    for p,want in { .90:.0181980515339464,.95:.0608757210636100,1.:.1035533905932738 }.items():
        require(abs(max(0,(p*SQUANTUM-6)/8)-want)<2e-14,'STOP: Sigma curve mismatch')
    for a in ALPHAS:
        ns=[sufficient_n(p,a) for p in VISIBILITIES]
        require(all(x>y for x,y in zip(ns,ns[1:])),'visibility monotonicity')
    require(sufficient_n(PSTAR,.01) is None,'threshold has no finite sufficient n')
    near=[sufficient_n(PSTAR+e,.01) for e in (.01,.001,.0001)]
    require(near[1]>99*near[0] and near[2]>99*near[1],'inverse-square divergence')
    for invalid in (-1,1.1,float('nan')):
        try:probabilities(invalid)
        except ValueError:pass
        else:raise AssertionError('invalid visibility accepted')
    path=ROOT/'threadB/data/shot_budget_results.json'
    if path.exists():
        obj=json.loads(path.read_text())
        for row in obj['rows']:
            n=sufficient_n(row['visibility'],row['alpha'])
            require(row['rigorous_scalar_N']==16*n,'stale analytic result')
            for test in ('scalar','directional'):
                e=row['monte_carlo'][test]
                require(e['N']<=16*n,'MC estimate exceeds sufficient bound')
                require(e['validation']['trials']>=2000,'too few validation trials')
            require(all(x['trials']>=2000 for x in row['search_points']),'too few search trials')
        first=obj['rows'][0];point=first['search_points'][len(first['search_points'])//2]
        require(simulate(first['visibility'],first['alpha'],point['n'],point['trials'])==point,
                'seeded replay mismatch; use pinned dependencies')
    print(json.dumps({'checks':'passed','born_cross_check':verified,'near_threshold_n':near},indent=2))

def render_table(obj):
    lines=['# Shot budgets under i.i.d. sampling','',
      'A: rigorous sufficient total N for at least 95% power; minimum integer allocation certified by this concentration bound, not the true minimum.',
      'B: seeded Monte Carlo estimated 95%-power crossing (0.5% search resolution); independent validation power and pointwise exact 95% binomial interval.','',
      '| p | α | A scalar | A directional | B scalar: N; power [CI] | B directional: N; power [CI] |',
      '|---|---|---:|---:|---|---|']
    for r in obj['rows']:
        vals=[]
        for t in ('scalar','directional'):
            e=r['monte_carlo'][t];lo,hi=e['cp95']
            vals.append(f"{e['N']:,}; {e['power']:.4f} [{lo:.4f}, {hi:.4f}]")
        lines.append(f"| {r['visibility']:.2f} | {r['alpha']:g} | {r['rigorous_scalar_N']:,} | {r['rigorous_directional_N']:,} | {vals[0]} | {vals[1]} |")
    lines+=['','The directional rejection region contains the scalar region for every dataset. Worst-case budgets coincide; Monte Carlo crossings have sampling uncertainty.','',
            '## Near the threshold','', '| p − p* | A total N (α=.01) |','|---|---:|']
    for row in obj['near_threshold']:lines.append(f"| {row['offset']:g} | {row['N']:,} |")
    lines+=['','p*=3−3√2/2. N grows as (p−p*)⁻²; no finite sufficient budget is certified at or below p*.','',
        f"Seed {SEED}; {obj['search_trials']:,} trials per search point; {obj['validation_trials']:,} independent validation trials per crossing.",
        'Intervals are pointwise Monte Carlo uncertainty, not simultaneous guarantees or experimental significance bounds.',
        'Full brackets, counts, environment: `data/shot_budget_results.json`.','']
    return '\n'.join(lines)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check',action='store_true');ap.add_argument('--study',action='store_true')
    ap.add_argument('--trials',type=int,default=5000);ap.add_argument('--validation-trials',type=int,default=20000)
    args=ap.parse_args()
    if args.check:check()
    if args.study:
        require(args.trials>=2000 and args.validation_trials>=2000,'at least 2,000 trials required')
        import scipy,platform
        obj={'method':'i.i.d. Hoeffding + Weissman; half error budget per family',
            'seed':SEED,'search_trials':args.trials,'validation_trials':args.validation_trials,
            'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'born_cross_check':cross_check(),'rows':[],
            'near_threshold':[{'offset':e,'N':16*sufficient_n(PSTAR+e,.01)} for e in (.01,.001,.0001)]}
        for p in VISIBILITIES:
            for alpha in ALPHAS:
                row=study(p,alpha,args.trials,args.validation_trials);obj['rows'].append(row)
                print(p,alpha,{t:row['monte_carlo'][t]['N'] for t in ('scalar','directional')},flush=True)
                (ROOT/'threadB/data/shot_budget_results.json').write_text(json.dumps(obj,indent=2)+'\n')
        (ROOT/'threadB/SHOT_BUDGET_RESULTS.md').write_text(render_table(obj))
    if not(args.check or args.study):ap.error('choose --check or --study')
if __name__=='__main__':main()
