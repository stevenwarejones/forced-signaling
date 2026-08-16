#!/usr/bin/env python3
"""Regenerates fig_geometry.pdf, fig_cover.pdf, fig_tradeoff.pdf (matplotlib)."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':8.5,'font.family':'serif','axes.linewidth':0.6})
INK='#1a1a1a'; BLUE='#2456a4'; ORANGE='#c46210'; GRAY='#8a8a8a'; LBLUE='#dce7f5'; LORANGE='#f7e8d4'
# geometry
fig,ax=plt.subplots(figsize=(3.4,2.6))
tA,tD,tBC=0.,0.35,3.; xB,xC=-5,5; cs=.55; vs=cs/8
for t0 in (tA,tD):
    xs=np.linspace(-6,6,10); ax.plot(xs,t0+np.abs(xs)*cs,color=GRAY,lw=.7,ls='--',zorder=1)
xs=np.linspace(-6,6,10); ax.fill_between(xs,tA+np.abs(xs)*vs,4.2,color=LBLUE,alpha=.5,zorder=0)
for x0 in (xB,xC):
    xs=np.linspace(x0-3,x0+3,10); ax.fill_between(xs,tBC+np.abs(xs-x0)*vs,4.2,color=LORANGE,alpha=.65,zorder=0)
ax.plot([0],[tA],'o',color=BLUE,ms=5,zorder=5); ax.text(.25,tA-.28,'$A$')
ax.plot([0],[tD],'o',color=BLUE,ms=5,zorder=5); ax.text(.25,tD+.1,'$D$')
ax.plot([xB],[tBC],'s',color=ORANGE,ms=5,zorder=5); ax.text(xB-.1,tBC+.25,'$B$',ha='center')
ax.plot([xC],[tBC+.12],'s',color=ORANGE,ms=5,zorder=5); ax.text(xC+.1,tBC+.37,'$C$',ha='center')
ax.plot([0,xB],[tD,tBC],color=GRAY,lw=.8); ax.plot([0,xC],[tD,tBC+.12],color=GRAY,lw=.8)
ax.annotate("$B\\sim C$: outside each other's $v$-cones",xy=(0,3.86),ha='center',fontsize=8)
ax.annotate('early $v$-cone reaches both',xy=(0,1.55),ha='center',color=BLUE,fontsize=7.5)
ax.annotate('$\\tau=t_B-t_C$ programmed',xy=(0,2.62),ha='center',color=ORANGE,fontsize=7.5)
ax.set_xlabel('space'); ax.set_ylabel('time (preferred frame)')
ax.set_xlim(-6.4,6.4); ax.set_ylim(-.4,4.3); ax.set_xticks([]); ax.set_yticks([])
for s in ('top','right'): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig('fig_geometry.pdf'); plt.close(fig)
# cover
B,K=1.34e-3,1e4
gB=1/np.sqrt(1-B*B); a=1/(gB*np.sqrt(K*K-1)); phi=np.arctan(a); U=np.arcsin(B*np.cos(phi))
N=int(np.floor(U/phi))+1; th=np.array([-U+(j+.5)*2*U/N for j in range(N)])
fig,ax=plt.subplots(figsize=(3.4,1.9))
ax.plot([-U,U],[0,0],color=INK,lw=2.4,solid_capstyle='butt',zorder=3)
for t in th:
    ax.axvspan(t-phi,t+phi,ymin=.30,ymax=.72,color=BLUE,alpha=.28,lw=0)
    ax.plot([t],[.21],marker='v',color=ORANGE,ms=4,zorder=5)
ax.annotate('covered interval $[-U,\\,U]$ of candidate frames',xy=(0,-.16),ha='center',fontsize=8)
ax.annotate(f'$N={N}$ overlapping arcs of half-width $\\phi$ (overlaps darker)',xy=(0,.52),ha='center',fontsize=8,color=BLUE)
ax.annotate('delay settings $\\tau_j=(d/c)\\sin\\theta_j$',xy=(0,.115),ha='center',fontsize=7.5,color=ORANGE)
ax.set_xlim(-U*1.15,U*1.15); ax.set_ylim(-.28,.62); ax.axis('off')
fig.tight_layout(); fig.savefig('fig_cover.pdf'); plt.close(fig)
# tradeoff
fig,ax=plt.subplots(figsize=(3.4,2.5))
x=np.linspace(0,.16,100)
ax.fill_between(x,5.4,6+8*x,color=LBLUE,alpha=.7,zorder=0)
ax.plot(x,6+8*x,color=BLUE,lw=1.6); ax.plot(x,6+16*x,color=BLUE,lw=1.,ls='--')
SQ=4+2*np.sqrt(2); ce=(np.sqrt(2)-1)/4
ax.plot([0],[SQ],marker='*',ms=11,color=ORANGE,zorder=6)
ax.annotate(r'quantum point $(0,\ 4+2\sqrt{2})$',xy=(.004,SQ+.015),fontsize=8)
ax.plot([ce],[SQ],'o',ms=4,color=INK,zorder=6)
ax.annotate(r'$\Delta_{\rm sig}\geq(\sqrt{2}-1)/4$',xy=(ce+.004,SQ-.09),fontsize=8.5)
ax.annotate('',xy=(ce-.003,SQ),xytext=(.004,SQ),arrowprops=dict(arrowstyle='->',color=ORANGE,lw=1.1))
ax.annotate('finite-$v$ models:\n$S_4^{\\rm op}\\leq 6+8\\Delta_{\\rm sig}$',xy=(.105,6.32),color=BLUE,fontsize=8,ha='center')
ax.annotate('default completion (slope 16)',xy=(.044,6.88),color=BLUE,fontsize=7.2,rotation=36)
ax.set_xlabel(r'measured signaling $\Delta_{\rm sig}$'); ax.set_ylabel(r'$S_4^{\rm op}$')
ax.set_xlim(0,.155); ax.set_ylim(5.9,7.05)
for s in ('top','right'): ax.spines[s].set_visible(False)
ax.grid(alpha=.25,lw=.4)
fig.tight_layout(); fig.savefig('fig_tradeoff.pdf'); plt.close(fig)
print("figures regenerated")
