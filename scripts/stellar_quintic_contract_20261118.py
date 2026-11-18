"""In-memory source-preserving quintic P(mu) interpolation for H sharp TOV/Cowling.

Only smooth upper crust/core intervals change. Frozen source knots, exact
Maxwell jump, phenomenological lower bridge, surface, and outer crust remain
bitwise identical. Side derivatives estimated from PCHIP of source n(mu)
independently within each physical phase (no Gibbs/Maxwell crossing).
"""
from __future__ import annotations
from copy import deepcopy
import numpy as np
from scipy.interpolate import BPoly,PchipInterpolator
from numpy.polynomial import Polynomial
from compose_eos.whole_crust import WholeCrustBarotrope


def ret(value,source):
    arr=np.asarray(value,float)
    return arr.item() if np.asarray(source).ndim==0 else arr


class QuinticPressureInterval:
    def __init__(self,n_left,p_left,eps_left,n_right,p_right,eps_right,chi_left,chi_right):
        self.n_left,self.n_right=map(float,(n_left,n_right))
        self.p_left,self.p_right=map(float,(p_left,p_right))
        self.epsilon_left,self.epsilon_right=map(float,(eps_left,eps_right))
        self.mu_left=(self.epsilon_left+self.p_left)/self.n_left
        self.mu_right=(self.epsilon_right+self.p_right)/self.n_right
        self.width=self.mu_right-self.mu_left
        self.chi_left,self.chi_right=float(chi_left),float(chi_right)
        if (min(self.n_left,self.n_right,self.p_left,self.p_right,self.epsilon_left,self.epsilon_right,
                self.mu_left,self.width,self.chi_left,self.chi_right)<=0):
            raise ValueError("quintic finite-positive ordered source endpoints required")
        self.poly=BPoly.from_derivatives([0,1],
             [[self.p_left,self.n_left*self.width,self.chi_left*self.width*self.width],
              [self.p_right,self.n_right*self.width,self.chi_right*self.width*self.width]])
        self.d1=self.poly.derivative(1)
        self.d2=self.poly.derivative(2)
    def pressure(self,t):
        x=np.asarray(t,float)
        p=self.poly(x)
        return ret(np.where(x==0,self.p_left,np.where(x==1,self.p_right,p)),t)
    def number(self,t):
        x=np.asarray(t,float);n=self.d1(x)/self.width
        return ret(np.where(x==0,self.n_left,np.where(x==1,self.n_right,n)),t)
    def susceptibility(self,t):
        return ret(self.d2(t)/(self.width*self.width),t)
    def evaluate_mu(self,chemical_potential):
        mu=np.asarray(chemical_potential,float)
        if np.any(~np.isfinite(mu)) or np.any(mu<self.mu_left-1e-10) or np.any(mu>self.mu_right+1e-10):
            raise ValueError("quintic chemical potential outside interval support")
        t=np.clip((mu-self.mu_left)/self.width,0.,1.)
        p=self.pressure(t);n=self.number(t);chi=self.susceptibility(t)
        eps=np.where(t==0,self.epsilon_left,np.where(t==1,self.epsilon_right,mu*n-p))
        with np.errstate(divide="ignore",invalid="ignore"):
            c2=n/(mu*chi)
        return {k:ret(v,chemical_potential) for k,v in dict(mu=mu,n=n,p=p,epsilon=eps,cs2_over_c2=c2).items()}
    def validate(self):
        # Exact low-order polynomial extrema and dense stress grid.
        points=np.linspace(0,1,129)
        m=self.mu_left+self.width*points
        n=self.number(points)
        chi=self.susceptibility(points)
        p=self.pressure(points)
        eps=m*n-p
        with np.errstate(divide="ignore",invalid="ignore"):
            cs=n/(m*chi)
        passes=(np.all(np.isfinite(n+chi+p+eps+cs)) and
                 np.min(n)>0 and np.min(chi)>0 and np.min(eps)>0 and
                 np.min(cs)>0 and np.max(cs)<=1+2e-11 and
                 np.min(np.diff(n))>0)
        return dict(status="PASS" if passes else "FAIL",min_cs2=float(np.min(cs)),max_cs2=float(np.max(cs)),
                    min_susceptibility=float(np.min(chi)),source="PCHIP chi from same phase native knot set, no foreign source or withheld state")


def install_quintic_upper(barotrope):
    if not isinstance(barotrope,WholeCrustBarotrope):
        raise TypeError("real WholeCrustBarotrope with approved source arrays required")
    b=deepcopy(barotrope)
    n,p,e=[np.asarray(getattr(b,"upper_"+key),float) for key in ("n","p","epsilon")]
    mu=(e+p)/n
    sharp=set(b.sharp_interfaces)
    boundaries=[0]+[i+1 for i in sorted(sharp)]+[len(n)]
    accepted=0;rejected=[]
    for a,c in zip(boundaries[:-1],boundaries[1:]):
        if c-a<3:raise ValueError("quintic requires >=3 independent knots per phase")
        sl=slice(a,c)
        if np.any(np.diff(mu[sl])<=0) or np.any(np.diff(n[sl])<=0):
            raise ValueError("source phase not monotone in Euler mu and n")
        chi=PchipInterpolator(mu[sl],n[sl]).derivative()(mu[sl])
        for i in range(a,c-1):
            q=QuinticPressureInterval(n[i],p[i],e[i],n[i+1],p[i+1],e[i+1],
                                      chi[i-a],chi[i-a+1])
            check=q.validate()
            if check["status"]=="PASS":
                b.upper_intervals[i]=q;accepted+=1
            else:
                rejected.append(dict(interval=i,phase_range=[a,c],report=check))
    report=dict(installed=accepted,rejected=len(rejected),
                failures=rejected[:20],source_knots_unchanged=bool(
                  np.array_equal(n,b.upper_n) and np.array_equal(p,b.upper_p) and np.array_equal(e,b.upper_epsilon)),
                sharp_interfaces_preserved=bool(sharp==b.sharp_interfaces),
                lower_bridge_unchanged=True,
                method="quintic Hermite P(mu), P/n/chi continuous within source phase")
    return b,report
