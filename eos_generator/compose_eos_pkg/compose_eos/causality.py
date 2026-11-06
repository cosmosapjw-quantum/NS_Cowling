"""Causality diagnostics and optional maximally-stiff causal completion.

A source EOS that is acausal is never silently clipped. A completion is a new,
model-dependent EOS in which the table is truncated at its first causal-limit
crossing and continued with a constant c_s^2 (default 1) branch.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import PchipInterpolator, CubicSpline
from scipy.optimize import brentq
from .eos_data import EOSData


def causality_profile(data: EOSData, *, n_min=0.0, max_log_step=0.25, tolerance=1e-6):
    n=np.asarray(data.nb,float); p=np.asarray(data.p,float); e=np.asarray(data.eps,float)
    x=np.log(n); lp=np.log(p); le=np.log(e)
    dln=np.diff(x); mask=np.ones(len(n),bool)
    mask[1:] &= dln<=max_log_step; mask[:-1] &= dln<=max_log_step
    if data.Iphase.shape==(len(n),):
        edges=np.where(np.diff(data.Iphase)!=0)[0]
        for i in edges: mask[max(0,i-1):min(len(n),i+3)]=False
    mask &= n>=n_min
    out={}
    for name in ('pchip','cubic_spline','centered_fd'):
        if name=='pchip': dp=PchipInterpolator(x,lp).derivative()(x); de=PchipInterpolator(x,le).derivative()(x)
        elif name=='cubic_spline': dp=CubicSpline(x,lp,bc_type='natural').derivative()(x); de=CubicSpline(x,le,bc_type='natural').derivative()(x)
        else: dp=np.gradient(lp,x,edge_order=2); de=np.gradient(le,x,edge_order=2)
        cs=(p/e)*dp/de; valid=mask & np.isfinite(cs) & (de>0)
        ids=np.where(valid & (cs>1+tolerance))[0]
        out[name]={'max_cs2':float(np.max(cs[valid])) if np.any(valid) else None,
          'min_cs2':float(np.min(cs[valid])) if np.any(valid) else None,
          'acausal_count':int(ids.size),'acausal_nb_range':[float(n[ids[0]]),float(n[ids[-1]])] if ids.size else None,
          'first_acausal_nb':float(n[ids[0]]) if ids.size else None,'n_valid':int(valid.sum()),
          'acausal':bool(ids.size)}
    robust=all(out[k]['acausal'] for k in out)
    onset=[out[k]['first_acausal_nb'] for k in out if out[k]['first_acausal_nb'] is not None]
    out['classification']={'robust_acausality':robust,'onset_spread_fm3':float(max(onset)-min(onset)) if onset else 0.0,
       'root_cause':'source_EOS_acausality' if robust else ('derivative_scheme_sensitive' if onset else 'no_acausality_detected')}
    return out


def causal_completion(data: EOSData, *, cs2=1.0, n_max=None, n_points=200):
    """Create a separate thermodynamically matched constant-sound-speed tail.

    The original EOS object is not modified. The completion is only a causal
    model alternative, not a correction uniquely implied by the source physics.
    """
    s=float(cs2)
    if not np.isfinite(s) or not 0<s<=1: raise ValueError('cs2 must be in (0,1]')
    n=np.asarray(data.nb,float); P=np.asarray(data.p,float); eps=np.asarray(data.eps,float)
    x=np.log(n); lp=np.log(P); le=np.log(eps)
    p_spline=PchipInterpolator(x,lp); e_spline=PchipInterpolator(x,le)
    dp=p_spline.derivative()(x); de=e_spline.derivative()(x)
    c=P/eps*dp/de
    brackets=[(i,i+1) for i in range(len(n)-1) if np.isfinite(c[i:i+2]).all() and c[i]<=1 and c[i+1]>1 and np.diff(x)[i]<=0.25]
    if not brackets: raise ValueError('no resolved first c_s^2=1 crossing found; causal completion not applied')
    i,j=brackets[0]
    def cs2_minus_one(z):
        return float(np.exp(p_spline(z)-e_spline(z))*p_spline.derivative()(z)/e_spline.derivative()(z)-1.0)
    root=brentq(cs2_minus_one,x[i],x[j],xtol=1e-13,rtol=1e-13)
    nt=float(np.exp(root)); pt=float(np.exp(p_spline(root))); et=float(np.exp(e_spline(root)))
    mu_t=(et+pt)/nt
    mu_source=float(PchipInterpolator(x,np.asarray(data.mub,float))(root))
    mu_match_rel=abs(mu_t-mu_source)/max(abs(mu_t),1.0)
    if mu_match_rel>5e-3: raise ValueError(f'first-law chemical potential mismatch at causal completion point: {mu_match_rel:.3e}')
    if n_max is None: n_max=float(n[-1])
    nmax=float(n_max)
    if nmax<=nt: raise ValueError('n_max must exceed causal-limit density')
    # Match mu_t=(epsilon_t+P_t)/n_t and integrate dP/dn=s(ε+P)/n.
    w=et+pt; tail=np.geomspace(nt,nmax,max(4,int(n_points)))
    q=tail/nt; p_tail=pt+s/(1+s)*w*(q**(1+s)-1); mu_tail=mu_t*q**s; e_tail=mu_tail*tail-p_tail
    keep=n<nt
    n_out=np.r_[n[keep],tail]; p_out=np.r_[P[keep],p_tail]; e_out=np.r_[eps[keep],e_tail]
    # Enforce the same Euler chemical potential on both pieces; source-table
    # roundoff/interpolation in Q3 must not create a numerical jump at n_t.
    mu_out=(e_out+p_out)/n_out
    if np.any(np.diff(n_out)<=0) or np.any(np.diff(p_out)<=0) or np.any(np.diff(e_out)<=0): raise ValueError('causal completion is not strictly monotone')
    mn=float(data.mn); F=e_out/n_out
    eos=EOSData(nb=n_out,mn=mn,mp=float(data.mp),Q1=p_out/n_out,Q2=np.zeros_like(n_out),
       Q3=mu_out/mn-1,Q4=np.zeros_like(n_out),Q5=np.zeros_like(n_out),Q6=F/mn-1,Q7=F/mn-1,
       q_extra=np.zeros_like(n_out),Iphase=np.zeros(len(n_out),dtype=int),
       Yn=np.zeros_like(n_out),Yp=np.zeros_like(n_out),Ye=np.zeros_like(n_out),Ymu=np.zeros_like(n_out),
       Aav=np.zeros_like(n_out),Zav=np.zeros_like(n_out),Yh=np.zeros_like(n_out),Il=int(data.Il),
       iT=int(data.iT),iYq=int(data.iYq),inb_start=1,temperature_mev=0,charge_fraction_grid=0,
       q_extra_count=0,composition_available=False,beta_equilibrated=True)
    meta={'schema':'compose-eos-causal-completion-v1','method':'constant_sound_speed_extension',
      'cs2_extension':s,'n_match_fm3':nt,'P_match_MeV_fm3':pt,'eps_match_MeV_fm3':et,
      'n_max_fm3':nmax,'source_composition_omitted':bool(data.composition_available),
      'source_max_cs2':float(np.nanmax(c)),'source_first_limiter_index':int(i),
      'chemical_potential_match_relative_residual':float(mu_match_rel),
      'sound_speed_continuity_at_join':bool(abs(s-1.0)<1e-12),
      'warning':'This is an explicit constant-cs alternative, not a unique repair of the microscopic source EOS. If cs2_extension differs from 1, c_s^2 has a finite kink at the first causal-limit point even though P, epsilon and mu_B remain continuous.'}
    return eos,meta
