"""Thermodynamically explicit Maxwell/Gibbs constructions and phase-aware TOV.

The legacy matcher remains a smooth core--crust interpolation. This module adds
separate, physically distinct sharp Maxwell and global-neutrality Gibbs paths.
Maxwell data carry an explicit pressure-continuous, density-discontinuous event;
Gibbs data require two-charge pressure surfaces P(mu_B,mu_Q), never a 1-D beta-
equilibrium barotrope. The l=2 tidal jump uses the pressure-corrected
Takatsy--Kovacs (2020) condition.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional
import json, os, tempfile
import numpy as np
from scipy.interpolate import PchipInterpolator, RectBivariateSpline
from scipy.optimize import brentq, root_scalar
from scipy.integrate import solve_ivp
from .eos_data import EOSData
from .io import read_compose, write_compose

MAXWELL_METHOD_ID = 'maxwell'

@dataclass
class MaxwellEOS:
    data: EOSData
    transition: dict
    metadata: dict = field(default_factory=dict)

@dataclass
class GibbsResult:
    eos: EOSData
    metadata: dict

@dataclass
class CausalCompletion:
    eos: EOSData
    metadata: dict


def _phase_arrays(data: EOSData):
    n=np.asarray(data.nb,float); mu=np.asarray(data.mub,float); p=np.asarray(data.p,float); eps=np.asarray(data.eps,float)
    if n.ndim!=1 or len(n)<5 or eps.shape!=n.shape or np.any(~np.isfinite(n+mu+p+eps)) or np.any(n<=0) or np.any(p<=0) or np.any(eps<=0):
        raise ValueError('phase EOS needs finite positive 1-D n, P, epsilon and mu arrays')
    if np.any(np.diff(n)<=0) or np.any(np.diff(mu)<=0) or np.any(np.diff(p)<=0) or np.any(np.diff(eps)<=0):
        raise ValueError('each phase must have strictly increasing n_b, mu_B, P and epsilon')
    # A sharp coexistence construction is only meaningful when each input is
    # already a stable, causal, cold barotrope with its first-law identity.
    euler=mu*n-p
    euler_rel=float(np.max(np.abs(euler-eps)/np.maximum(np.abs(eps),1.0)))
    if euler_rel>5e-3:
        raise ValueError(f'phase EOS violates epsilon=mu_B*n_B-P (max relative residual {euler_rel:.3e})')
    if not bool(data.beta_equilibrated):
        raise ValueError('1-D Maxwell input must be a beta-equilibrated, locally neutral phase barotrope')
    from .consistency import check_consistency
    rep=check_consistency(data,verbose=False)
    for key in ('C1','C2','C6','C10','C11'):
        if not bool(rep[key].get('pass_',False)):
            raise ValueError(f'phase {key} hard gate failed; Maxwell construction requires stable causal source phases')
    if data.composition_available and (not rep['C3']['pass_'] or not rep['C4']['pass_']):
        raise ValueError('phase composition fails local baryon/charge closure')
    return n,mu,p,eps,rep


def build_maxwell_eos(low_phase: EOSData, high_phase: EOSData, *, low_name='low', high_name='high',
                      mu_samples=800, crossing_index=0) -> MaxwellEOS:
    """Construct one exact local-neutrality Maxwell transition from two 1-D phases.

    Coexistence solves P_low(mu_B)=P_high(mu_B). The returned table contains
    two adjacent rows at the same P_t and mu_t but n_high>n_low and eps_high>eps_low.
    The two input EOS are understood as locally charge-neutral phase barotropes.
    """
    nl,mul,pl,el,repl=_phase_arrays(low_phase); nh,muh,ph,eh,reph=_phase_arrays(high_phase)
    lo=max(mul[0],muh[0]); hi=min(mul[-1],muh[-1])
    if not lo<hi: raise ValueError('the two phases have no common baryon-chemical-potential interval')
    Pl=PchipInterpolator(mul,pl,extrapolate=False); Ph=PchipInterpolator(muh,ph,extrapolate=False)
    nli=PchipInterpolator(mul,nl,extrapolate=False); nhi=PchipInterpolator(muh,nh,extrapolate=False)
    eli=PchipInterpolator(mul,el,extrapolate=False); ehi_sp=PchipInterpolator(muh,eh,extrapolate=False)
    grid=np.linspace(lo,hi,int(mu_samples)); delta=np.asarray(Pl(grid)-Ph(grid),float)
    brackets=[(grid[i],grid[i+1]) for i in range(len(grid)-1) if delta[i]==0 or delta[i]*delta[i+1]<0]
    if not brackets: raise ValueError('no P_low(mu_B)=P_high(mu_B) coexistence root in common interval')
    if crossing_index<0 or crossing_index>=len(brackets): raise ValueError(f'{len(brackets)} coexistence roots found; choose crossing_index explicitly')
    a,b=brackets[crossing_index]
    mut=float(a if delta[np.searchsorted(grid,a)]==0 else brentq(lambda x: float(Pl(x)-Ph(x)),a,b,xtol=1e-12))
    pt=float(0.5*(Pl(mut)+Ph(mut))); nlo=float(nli(mut)); nhi_t=float(nhi(mut))
    if not (pt>0 and mut>0 and nlo>0 and nhi_t>nlo):
        raise ValueError(f'coexistence root is not a physical density jump: P={pt}, mu={mut}, n_low={nlo}, n_high={nhi_t}')
    elo=float(mut*nlo-pt); epshi_t=float(mut*nhi_t-pt); de=epshi_t-elo
    euler_low_rel=abs(elo-float(eli(mut)))/max(abs(elo),1.0)
    euler_high_rel=abs(epshi_t-float(ehi_sp(mut)))/max(abs(epshi_t),1.0)
    if max(euler_low_rel,euler_high_rel)>5e-3:
        raise ValueError('interpolated coexistence point violates the phase first-law identity beyond 5e-3')
    if not (elo>0 and epshi_t>elo and de>0): raise ValueError('Maxwell energy-density jump is non-positive')
    # Keep the untouched tabulated branches and add exact interpolated endpoints.
    ml=mul<mut; mh=muh>mut
    if ml.sum()<3 or mh.sum()<3: raise ValueError('not enough tabulated points on both sides of the coexistence point')
    n_out=np.r_[nl[ml],nlo,nhi_t,nh[mh]]
    mu_out=np.r_[mul[ml],mut,mut,muh[mh]]
    p_out=np.r_[pl[ml],pt,pt,ph[mh]]
    phase_ids=np.r_[np.ones(ml.sum()+1,dtype=int),np.full(1+mh.sum(),2,dtype=int)]
    eps_out=mu_out*n_out-p_out
    if np.any(np.diff(n_out)<=0) or np.any(np.diff(p_out)<-1e-10) or np.any(np.diff(eps_out)<=0):
        raise ValueError('constructed Maxwell table is not monotone in n_b, P, and epsilon')
    mn=float(low_phase.mn); F=eps_out/n_out
    data=EOSData(nb=n_out,mn=mn,mp=float(low_phase.mp),
        Q1=p_out/n_out,Q2=np.zeros_like(n_out),Q3=mu_out/mn-1.0,
        Q4=np.zeros_like(n_out),Q5=np.zeros_like(n_out),Q6=F/mn-1.0,Q7=F/mn-1.0,
        q_extra=np.zeros_like(n_out),Iphase=phase_ids,
        Yn=np.zeros_like(n_out),Yp=np.zeros_like(n_out),Ye=np.zeros_like(n_out),Ymu=np.zeros_like(n_out),
        Aav=np.zeros_like(n_out),Zav=np.zeros_like(n_out),Yh=np.zeros_like(n_out),
        Il=1,iT=0,iYq=1,inb_start=1,temperature_mev=0.0,charge_fraction_grid=0.0,
        q_extra_count=0,composition_available=False,beta_equilibrated=True)
    il=int(ml.sum()); ih=il+1
    critical=0.5*elo+1.5*pt
    transition={'kind':'maxwell','method_id':MAXWELL_METHOD_ID,'pressure_MeV_fm3':pt,'muB_MeV':mut,
       'n_low_fm3':nlo,'n_high_fm3':nhi_t,'eps_low_MeV_fm3':elo,'eps_high_MeV_fm3':epshi_t,
       'delta_eps_MeV_fm3':de,'low_index':il,'high_index':ih,
       'seidov_critical_delta_eps_MeV_fm3':critical,
       'seidov_delta_eps_over_critical':de/critical,
       'seidov_small_core_stable':bool(de<critical),
       'tidal_jump_formula':'Takatsy-Kovacs-2020: Delta y=-4 pi r_t^3 Delta epsilon/(m_t+4 pi r_t^3 P_t)',
       'charge_neutrality':'local in each phase (Maxwell)'}
    # Store provenance and interpolation consistency at the coexistence point.
    source_e_low=float(eli(mut)); source_e_high=float(ehi_sp(mut))
    meta={'schema':'compose-eos-maxwell-v1','method_id':MAXWELL_METHOD_ID,'low_phase':low_name,'high_phase':high_name,
          'coexistence_root_count':len(brackets),'selected_root_index':int(crossing_index),
          'low_energy_identity_residual':elo-source_e_low,'high_energy_identity_residual':epshi_t-source_e_high,
          'low_source_hard_gates':{k:bool(repl[k].get('pass_',False)) for k in ('C1','C2','C6','C10','C11')},
          'high_source_hard_gates':{k:bool(reph[k].get('pass_',False)) for k in ('C1','C2','C6','C10','C11')},
          'transition':transition,'phase_ids':phase_ids.tolist(),
          'interpretation':'Exact sharp Maxwell interface between the supplied locally-neutral 1-D barotropes; not the legacy smooth maxwell bridge.'}
    return MaxwellEOS(data=data,transition=transition,metadata=meta)


def write_maxwell_eos(directory: str, eos: MaxwellEOS) -> None:
    """Write CompOSE arrays plus an explicit non-standard phase-event sidecar."""
    write_compose(directory,eos.data)
    target=Path(directory)/'eos.phase_transitions.json'; tmp=target.with_suffix('.json.tmp')
    payload={**eos.metadata,'transition':eos.transition,'phase_ids':eos.data.Iphase.astype(int).tolist()}
    tmp.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n',encoding='utf-8'); os.replace(tmp,target)


def read_maxwell_eos(directory: str) -> MaxwellEOS:
    data=read_compose(directory); path=Path(directory)/'eos.phase_transitions.json'
    if not path.is_file(): raise FileNotFoundError('EOS directory has no eos.phase_transitions.json sharp-interface metadata')
    meta=json.loads(path.read_text(encoding='utf-8'))
    if meta.get('schema')!='compose-eos-maxwell-v1' or meta.get('transition',{}).get('kind')!='maxwell' or meta.get('method_id', MAXWELL_METHOD_ID)!=MAXWELL_METHOD_ID:
        raise ValueError('unsupported or malformed Maxwell sidecar schema')
    ids=np.asarray(meta.get('phase_ids',[]),dtype=int)
    if ids.shape!=(data.n,): raise ValueError('phase_ids length does not match CompOSE table')
    data.Iphase=ids
    t=meta['transition']; il=int(t['low_index']); ih=int(t['high_index'])
    if ih!=il+1 or not (0<=il<data.n-1) or abs(data.p[il]-data.p[ih])>1e-8*max(1.0,abs(data.p[il])):
        raise ValueError('Maxwell sidecar indices do not identify a pressure-continuous adjacent density jump')
    comparisons=(('P_low',data.p[il],t['pressure_MeV_fm3']),('P_high',data.p[ih],t['pressure_MeV_fm3']),
      ('mu_low',data.mub[il],t['muB_MeV']),('mu_high',data.mub[ih],t['muB_MeV']),
      ('n_low',data.nb[il],t['n_low_fm3']),('n_high',data.nb[ih],t['n_high_fm3']),
      ('eps_low',data.eps[il],t['eps_low_MeV_fm3']),('eps_high',data.eps[ih],t['eps_high_MeV_fm3']))
    for label,actual,expected in comparisons:
        if abs(float(actual)-float(expected))>2e-8*max(abs(float(expected)),1.0):
            raise ValueError(f'Maxwell sidecar {label} mismatch: table={actual}, metadata={expected}')
    if not (data.nb[ih]>data.nb[il] and data.eps[ih]>data.eps[il]): raise ValueError('Maxwell jump is not upward in n and epsilon')
    if int(ids[il])!=1 or int(ids[ih])!=2: raise ValueError('Maxwell sidecar phase IDs do not bracket the density jump')
    return MaxwellEOS(data,t,meta)


def tidal_y_jump(y_inner: float, radius_km: float, mass_geometric_km: float,
                 pressure_geometric_km2: float, delta_epsilon_geometric_km2: float) -> float:
    """Apply the pressure-corrected l=2 tidal jump in G=c=1 units.

    ``mass_geometric_km`` is m(r_t) in km; pressure and energy density are km^-2.
    The returned value is y just outside the interface.
    """
    vals=np.asarray([y_inner,radius_km,mass_geometric_km,pressure_geometric_km2,delta_epsilon_geometric_km2],float)
    if np.any(~np.isfinite(vals)) or radius_km<=0 or mass_geometric_km<0 or delta_epsilon_geometric_km2<=0:
        raise ValueError('invalid tidal jump arguments')
    den=mass_geometric_km+4*np.pi*radius_km**3*pressure_geometric_km2
    if den<=0: raise ValueError('non-positive tidal-jump denominator')
    return float(y_inner-4*np.pi*radius_km**3*delta_epsilon_geometric_km2/den)


# Exact unit conversion: (MeV fm^-3) -> km^-2 in G=c=1.
_G=6.67430e-11; _C=299792458.0; _MEV=1.602176634e-13; _MSUN=1.988409870698051e30
_MEVFM3_TO_KM2=_G*_MEV/(1e-45*_C**4)*1e6
_MSUN_GEOM_KM=_G*_MSUN/_C**2/1000.0


def solve_tov_with_phase_transition(eos: MaxwellEOS, central_pressure_MeV_fm3: float,
                                    *, rtol=2e-8, max_step_km=0.12) -> dict:
    """Solve TOV plus l=2 Love equation with one sharp Maxwell interface.

    Integrates the high-density phase to P_t, applies the corrected y-jump,
    then integrates the low-density phase to the finite tabulated surface.
    All stellar structure equations are in G=c=1 geometric units.
    """
    d=eos.data; t=eos.transition; pt=float(t['pressure_MeV_fm3'])
    pc=float(central_pressure_MeV_fm3); pgrid=d.p; eps=d.eps
    if not np.isfinite(pc) or pc<=pgrid[0] or pc>pgrid[-1]: raise ValueError('central pressure outside tabulated EOS range')
    il=int(t['low_index']); ih=int(t['high_index'])
    def make_branch(start,stop):
        pp=np.asarray(pgrid[start:stop],float); ee=np.asarray(eps[start:stop],float)
        if len(pp)<4 or np.any(np.diff(pp)<=0) or np.any(np.diff(ee)<=0): raise ValueError('each smooth phase branch needs >=4 strictly increasing P/epsilon points')
        sp=PchipInterpolator(np.log(pp),np.log(ee),extrapolate=False)
        def state(P):
            lp=np.log(np.asarray(P,float)); ev=np.exp(sp(lp)); dedp=ev/np.asarray(P,float)*sp.derivative()(lp)
            if np.any(~np.isfinite(ev+dedp)) or np.any(dedp<=0): raise FloatingPointError('invalid d epsilon/dP in phase branch')
            return ev,dedp
        return (pp,ee,state)
    low=make_branch(0,il+1); high=make_branch(ih,d.n)
    active_high=pc>pt
    branch=high if active_high else low
    ec,de_dp=branch[2](pc); ec=float(ec); de_dp=float(de_dp); cs2c=1.0/de_dp
    if not (0<cs2c<=1+1e-5): raise ValueError(f'central sound speed is noncausal: {cs2c:g}')
    pcg=pc*_MEVFM3_TO_KM2; ecg=ec*_MEVFM3_TO_KM2
    r0=1e-4
    p0=pcg-(2*np.pi/3)*(ecg+pcg)*(ecg+3*pcg)*r0*r0
    m0=4*np.pi/3*ecg*r0**3
    y0=2-(4*np.pi/7)*(5*ecg+9*pcg+(ecg+pcg)/cs2c)*r0*r0
    if p0<=0 or (active_high and p0<=pt*_MEVFM3_TO_KM2): raise ValueError('central expansion step is too large or starts below the Maxwell interface')
    transitions=[]
    def integrate(beg_state, r_start, target_P, phase):
        pmin=float(phase[0][0]); pmax=float(phase[0][-1])
        if target_P<pmin*(1-1e-10) or target_P>pmax*(1+1e-10): raise ValueError('TOV event lies outside selected phase table')
        target_P_geom = target_P * _MEVFM3_TO_KM2
        def rhs(r,s):
            # The integrated pressure state is geometric (km^-2). Convert back
            # to CompOSE MeV/fm^3 before querying the phase interpolator.
            P_geom=max(float(s[0]),target_P_geom)
            P_mev=P_geom/_MEVFM3_TO_KM2
            eps_v,dedp_v=phase[2](P_mev)
            pg=P_geom; eg=float(eps_v)*_MEVFM3_TO_KM2
            m,y=s[1],s[2]; fac=1-2*m/r
            if fac<=0: raise FloatingPointError('trapped-surface condition reached')
            F=(1-4*np.pi*r*r*(eg-pg))/fac
            nu_prime=2*(m+4*np.pi*r**3*pg)/(r*(r-2*m))
            Q=4*np.pi*(5*eg+9*pg+(eg+pg)*float(dedp_v))/fac-6/(r*r*fac)-nu_prime**2
            dpdr=-(eg+pg)*(m+4*np.pi*r**3*pg)/(r*(r-2*m))
            dmdr=4*np.pi*r*r*eg
            dydr=-(y*y+y*F+r*r*Q)/r
            return [dpdr,dmdr,dydr]
        def event(r,s): return s[0]-target_P_geom
        event.terminal=True; event.direction=-1
        sol=solve_ivp(rhs,(r_start,40.0),beg_state,method='DOP853',rtol=rtol,atol=[1e-13,1e-12,1e-10],max_step=max_step_km,events=event)
        if not sol.success or not sol.t_events or len(sol.t_events[0])!=1: raise RuntimeError(f'TOV phase integration missed pressure event {target_P:g}: {sol.message}')
        return float(sol.t_events[0][0]),[float(v) for v in sol.y_events[0][0]]
    r=r0; state=[p0,m0,y0]
    if active_high:
        r,state=integrate(state,r,pt,high)
        de=float(t['delta_eps_MeV_fm3'])*_MEVFM3_TO_KM2
        pgeom=pt*_MEVFM3_TO_KM2
        y_before=state[2]; state[2]=tidal_y_jump(y_before,r,state[1],pgeom,de)
        transitions.append({'radius_km':r,'mass_geometric_km':float(state[1]),'pressure_geometric_km2':pgeom,
          'delta_epsilon_geometric_km2':de,'y_inside':y_before,'y_outside':state[2],'delta_y':state[2]-y_before,
          'delta_y_formula':'-4*pi*r_t^3*Delta_epsilon/(m_t+4*pi*r_t^3*P_t)'})
        r,state=integrate(state,r,float(pgrid[0]),low)
    else:
        r,state=integrate(state,r,float(pgrid[0]),low)
    R=r; Mgeom=state[1]; yR=state[2]; C=Mgeom/R
    if not (0<C<0.5): raise RuntimeError(f'unphysical compactness {C:g}')
    t1=(8/5)*C**5*(1-2*C)**2*(2+2*C*(yR-1)-yR)
    t2=2*C*(6-3*yR+3*C*(5*yR-8))
    t3=4*C**3*(13-11*yR+C*(3*yR-2)+2*C**2*(1+yR))
    t4=3*(1-2*C)**2*(2-yR+2*C*(yR-1))*np.log1p(-2*C)
    den=t2+t3+t4
    if abs(den)<1e-30: raise RuntimeError('Love-number denominator is singular')
    k2=float(t1/den); lam=float((2/3)*k2/C**5)
    return {'status':'SUCCESS','mass_Msun':Mgeom/_MSUN_GEOM_KM,'radius_km':R,'compactness':C,
            'k2':k2,'Lambda':lam,'central_pressure_MeV_fm3':pc,'surface_y':yR,
            'phase_interfaces':transitions,'tidal_jump_formula':'Takatsy-Kovacs-2020 Eq. (11)',
            'seidov_small_core_stable':bool(eos.transition['seidov_small_core_stable']),
            'radial_stability_scope':'Seidov small-core onset criterion only; no full radial-mode eigenfrequency calculation is performed.'}


@dataclass
class ChargeResolvedPhase:
    """Tabulated two-charge phase surface P(mu_B,mu_Q), in MeV fm^-3."""
    mu_b: np.ndarray
    mu_q: np.ndarray
    pressure: np.ndarray
    label: str='phase'
    _spline: object=field(init=False,repr=False)
    def __post_init__(self):
        self.mu_b=np.asarray(self.mu_b,float); self.mu_q=np.asarray(self.mu_q,float); self.pressure=np.asarray(self.pressure,float)
        if self.mu_b.ndim!=1 or self.mu_q.ndim!=1 or len(self.mu_b)<4 or len(self.mu_q)<4 or self.pressure.shape!=(len(self.mu_b),len(self.mu_q)):
            raise ValueError('two-charge surface needs >=4x4 pressure grid with shape (mu_b,mu_q)')
        if np.any(np.diff(self.mu_b)<=0) or np.any(np.diff(self.mu_q)<=0) or np.any(~np.isfinite(self.pressure)):
            raise ValueError('phase surface grids must be finite and strictly increasing')
        self._spline=RectBivariateSpline(self.mu_b,self.mu_q,self.pressure,kx=3,ky=3,s=0)
    @classmethod
    def from_npz(cls,path,label='phase'):
        z=np.load(path)
        return cls(z['mu_b'],z['mu_q'],z['pressure'],label)
    def evaluate(self,mu_b,mu_q):
        if not(self.mu_b[0]<=mu_b<=self.mu_b[-1] and self.mu_q[0]<=mu_q<=self.mu_q[-1]): raise ValueError('chemical potentials outside tabulated phase surface')
        P=float(self._spline.ev(mu_b,mu_q)); n=float(self._spline.ev(mu_b,mu_q,dx=1,dy=0)); q=float(self._spline.ev(mu_b,mu_q,dx=0,dy=1))
        hbb=float(self._spline.ev(mu_b,mu_q,dx=2,dy=0)); hbq=float(self._spline.ev(mu_b,mu_q,dx=1,dy=1)); hqq=float(self._spline.ev(mu_b,mu_q,dx=0,dy=2))
        hess=np.array([[hbb,hbq],[hbq,hqq]],float); eig=np.linalg.eigvalsh(hess)
        scale=max(float(np.max(np.abs(hess))),1e-12)
        if np.min(eig)<-1e-7*scale: raise ValueError(f'{self.label} pressure surface is not convex in (mu_B,mu_Q): Hessian eigenvalues {eig}')
        if not np.isfinite(P+n+q) or n<=0: raise ValueError(f'{self.label} has non-positive/non-finite baryon density')
        return P,n,q

@dataclass
class LeptonSector:
    pressure: Callable[[float],float]
    charge_density: Callable[[float],float]
    label: str='custom-lepton-sector'
    @staticmethod
    def electron_muon():
        hbarc=197.3269804; masses=(0.51099895,105.6583755)
        def p_of_mu_q(mu_q):
            mu_e=-float(mu_q); out=0.0
            for mass in masses:
                if mu_e>mass:
                    k=np.sqrt(mu_e**2-mass**2)
                    out+=(mu_e*k*(2*mu_e**2-5*mass**2)+3*mass**4*np.log((mu_e+k)/mass))/(24*np.pi**2*hbarc**3)
            return float(out)
        def q_of_mu_q(mu_q):
            mu_e=-float(mu_q); out=0.0
            for mass in masses:
                if mu_e>mass: out-=(mu_e**2-mass**2)**1.5/(3*np.pi**2*hbarc**3)
            return float(out)
        return LeptonSector(p_of_mu_q,q_of_mu_q,'zero-temperature electron+muon Fermi gases')


def gibbs_mixed_phase(hadronic: ChargeResolvedPhase, high_density: ChargeResolvedPhase,
                      mu_b_values, *, lepton_sector: Optional[LeptonSector]=None,
                      mu_q_bounds=None) -> GibbsResult:
    """Construct a globally neutral Gibbs mixed phase from two-charge surfaces.

    Requires P_i(mu_B,mu_Q), including charge-resolved baryon/charge densities.
    A one-dimensional beta-equilibrated CompOSE barotrope is insufficient by design.
    """
    if not isinstance(hadronic,ChargeResolvedPhase) or not isinstance(high_density,ChargeResolvedPhase):
        raise TypeError('Gibbs construction requires two ChargeResolvedPhase P(mu_B,mu_Q) surfaces; 1-D EOSData is insufficient')
    lb=max(hadronic.mu_b[0],high_density.mu_b[0]); ub=min(hadronic.mu_b[-1],high_density.mu_b[-1])
    mus=np.asarray(mu_b_values,float)
    if mus.ndim!=1 or len(mus)<4 or np.any(np.diff(mus)<=0) or mus[0]<lb or mus[-1]>ub: raise ValueError('mu_b_values must be an increasing grid inside both phase surfaces')
    qlo=max(hadronic.mu_q[0],high_density.mu_q[0]); qhi=min(hadronic.mu_q[-1],high_density.mu_q[-1])
    if mu_q_bounds is not None:qlo=max(qlo,float(mu_q_bounds[0]));qhi=min(qhi,float(mu_q_bounds[1]))
    if qlo>=qhi:raise ValueError('two phases have no common charge-chemical-potential interval')
    leptons=lepton_sector or LeptonSector.electron_muon()
    out=[]; pressure_equal_residuals=[]
    for mub in mus:
        qscan=np.linspace(qlo,qhi,160); vals=np.array([hadronic.evaluate(mub,q)[0]-high_density.evaluate(mub,q)[0] for q in qscan])
        br=[(qscan[i],qscan[i+1]) for i in range(len(qscan)-1) if vals[i]==0 or vals[i]*vals[i+1]<0]
        if len(br)!=1:raise ValueError(f'Gibbs pressure equality has {len(br)} roots at mu_B={mub:g}; a unique coexistence branch is required')
        qmu=float(brentq(lambda q:hadronic.evaluate(mub,q)[0]-high_density.evaluate(mub,q)[0],*br[0]))
        ph,nh,qh=hadronic.evaluate(mub,qmu); pq,nq,qq=high_density.evaluate(mub,qmu)
        pressure_equal_residuals.append(abs(ph-pq))
        ql=float(leptons.charge_density(qmu)); pl=float(leptons.pressure(qmu))
        denom=qq-qh
        if abs(denom)<1e-12:raise ValueError('phase charge densities are indistinguishable; mixed fraction is undefined')
        chi=float(-(qh+ql)/denom)
        if chi < -1e-8 or chi > 1+1e-8:raise ValueError(f'global neutrality requires mixed fraction chi={chi:g} outside [0,1] at mu_B={mub:g}')
        chi=float(np.clip(chi,0,1)); n=(1-chi)*nh+chi*nq; P=0.5*(ph+pq)+pl
        charge=(1-chi)*qh+chi*qq+ql; eps=mub*n-P
        if n<=0 or P<=0 or eps<=0 or abs(charge)>2e-7:raise ValueError('Gibbs point violates positivity or global charge neutrality')
        out.append((n,P,eps,mub,qmu,chi,charge))
    a=np.asarray(out); n,P,eps,mub,qmu,chi,qres=a.T
    if np.any(np.diff(n)<=0) or np.any(np.diff(P)<=0) or np.any(np.diff(eps)<=0):raise ValueError('Gibbs mixed phase is not monotone/stable on requested mu_B grid')
    mn=float(939.0); F=eps/n
    eos=EOSData(nb=n,mn=mn,mp=mn,Q1=P/n,Q2=np.zeros_like(n),Q3=mub/mn-1,
        Q4=qmu/mn,Q5=np.zeros_like(n),Q6=F/mn-1,Q7=F/mn-1,q_extra=np.zeros_like(n),
        Iphase=np.full(len(n),3,dtype=int),Yn=np.zeros_like(n),Yp=np.zeros_like(n),Ye=np.zeros_like(n),Ymu=np.zeros_like(n),
        Aav=np.zeros_like(n),Zav=np.zeros_like(n),Yh=np.zeros_like(n),Il=1,iT=0,iYq=1,inb_start=1,
        temperature_mev=0,charge_fraction_grid=0,q_extra_count=0,composition_available=False,beta_equilibrated=True)
    meta={'schema':'compose-eos-gibbs-v1','hadronic_phase':hadronic.label,'high_density_phase':high_density.label,
      'lepton_sector':leptons.label,'mu_b_MeV':mub.tolist(),'mu_q_MeV':qmu.tolist(),'mixed_fraction_chi':chi.tolist(),
      'global_charge_residual':qres.tolist(),'pressure_equality_residual_max_MeV_fm3':float(max(pressure_equal_residuals)),
      'construction':'two conserved charges; pressure equality plus global charge neutrality',
      'input_requirement':'P_H(mu_B,mu_Q), P_Q(mu_B,mu_Q), derivatives n_B=partial P/partial mu_B and q=partial P/partial mu_Q'}
    return GibbsResult(eos,meta)


def write_gibbs_eos(directory: str, result: GibbsResult):
    write_compose(directory,result.eos)
    p=Path(directory)/'eos.gibbs_construction.json'; tmp=p.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(result.metadata,indent=2,allow_nan=False)+'\n',encoding='utf-8'); os.replace(tmp,p)


def read_gibbs_eos(directory: str) -> GibbsResult:
    data=read_compose(directory); path=Path(directory)/'eos.gibbs_construction.json'
    if not path.is_file(): raise FileNotFoundError('EOS directory has no eos.gibbs_construction.json provenance sidecar')
    meta=json.loads(path.read_text(encoding='utf-8'))
    if meta.get('schema')!='compose-eos-gibbs-v1': raise ValueError('unsupported Gibbs sidecar schema')
    fields=('mu_b_MeV','mu_q_MeV','mixed_fraction_chi','global_charge_residual')
    if any(len(meta.get(k,[]))!=data.n for k in fields): raise ValueError('Gibbs chemical-potential/mixture arrays do not match EOS row count')
    if np.max(np.abs(np.asarray(meta['global_charge_residual'],float)))>2e-7: raise ValueError('Gibbs sidecar violates global charge neutrality')
    if np.any(np.asarray(meta['mixed_fraction_chi'],float)<-1e-8) or np.any(np.asarray(meta['mixed_fraction_chi'],float)>1+1e-8): raise ValueError('Gibbs mixed fraction lies outside [0,1]')
    return GibbsResult(data,meta)


def maxwell_tidal_jump_report(eos: MaxwellEOS, radius_km, mass_geometric_km):
    t=eos.transition; return {'transition':dict(t),'radius_km':float(radius_km),'mass_geometric_km':float(mass_geometric_km),
        'delta_y_for_y_inside_1':tidal_y_jump(1.0,radius_km,mass_geometric_km,
          t['pressure_MeV_fm3']*_MEVFM3_TO_KM2,t['delta_eps_MeV_fm3']*_MEVFM3_TO_KM2)-1.0}


def build_maxwell_css_eos(hadronic: EOSData, pressure_transition_MeV_fm3: float,
                          delta_epsilon_MeV_fm3: float, cs2_high: float, *,
                          n_max_fm3: float, n_points: int=240,
                          overlap_fraction: float=0.01, hadronic_name='hadronic',
                          high_phase_name='CSS-high-density') -> MaxwellEOS:
    """Create an exact Maxwell transition with a phenomenological CSS phase.

    The hadronic phase is taken from the supplied cold beta-equilibrated table
    below P_t. The high-density phase obeys dP/dε=cs2_high and is fixed by the
    specified energy-density jump Δε at equal P_t and μ_B. A tiny metastable
    overlap brackets the chemical-potential crossing and is removed from the
    final table. This is a generic CSS parameterization, not a microscopic QCD EOS.
    """
    from scipy.interpolate import PchipInterpolator
    n=np.asarray(hadronic.nb,float); p=np.asarray(hadronic.p,float); eps=np.asarray(hadronic.eps,float)
    pt=float(pressure_transition_MeV_fm3); de=float(delta_epsilon_MeV_fm3); cs=float(cs2_high); nmax=float(n_max_fm3)
    if not hadronic.beta_equilibrated: raise ValueError('CSS Maxwell builder requires a beta-equilibrated one-dimensional hadronic EOS')
    if not np.all(np.isfinite([pt,de,cs,nmax,overlap_fraction])) or not p[0]<pt<p[-1] or de<=0 or not 0<cs<=1:
        raise ValueError('invalid transition pressure, energy jump or high-phase sound speed')
    if n_points<40 or not 0<overlap_fraction<0.1: raise ValueError('n_points>=40 and overlap_fraction in (0,0.1) required')
    x=np.log(n); lp=PchipInterpolator(x,np.log(p)); le=PchipInterpolator(x,np.log(eps))
    lognt=brentq(lambda z: float(lp(z)-np.log(pt)),float(x[0]),float(x[-1]),xtol=1e-13)
    nt=float(np.exp(lognt)); et=float(np.exp(le(lognt))); mut=float((et+pt)/nt)
    mu_source=float(PchipInterpolator(x,np.asarray(hadronic.mub,float))(lognt))
    mu_residual=abs(mut-mu_source)/max(abs(mut),1.0)
    if mu_residual>5e-3: raise ValueError(f'hadronic phase violates mu_B=(epsilon+P)/n_B at P_t by {mu_residual:.3e}')
    # Limit the metastable overlap to keep the CSS pressure positive below P_t.
    epsq=et+de; nq=(epsq+pt)/mut
    alpha=min(float(overlap_fraction),0.2*pt/(cs*(epsq+pt)))
    if alpha<1e-7: raise ValueError('P_t is too small to resolve a positive-pressure CSS metastable overlap')
    ncut=nt*(1+alpha)
    if ncut>=n[-1] or nmax<=nq: raise ValueError('hadronic table does not extend beyond P_t overlap or n_max does not exceed n_Q')
    nlow0=n[n<ncut]
    if len(nlow0)<5: raise ValueError('insufficient hadronic points below the transition')
    plow0=p[n<ncut]; elow0=eps[n<ncut]
    nlow=np.r_[nlow0,ncut]; plow=np.r_[plow0,float(np.exp(lp(np.log(ncut))))]; elow=np.r_[elow0,float(np.exp(le(np.log(ncut))))]
    mulow=(elow+plow)/nlow
    q4sp=PchipInterpolator(x,np.asarray(hadronic.muq,float)); q5sp=PchipInterpolator(x,np.asarray(hadronic.mul,float))
    q4low=q4sp(np.log(nlow))/float(hadronic.mn); q5low=q5sp(np.log(nlow))/float(hadronic.mn)
    def make_phase(nn,pp,ee,mm,q4,q5):
        z=np.zeros_like(nn); F=ee/nn; mn=float(hadronic.mn)
        return EOSData(nb=nn,mn=mn,mp=float(hadronic.mp),Q1=pp/nn,Q2=z,Q3=mm/mn-1,
          Q4=q4,Q5=q5,Q6=F/mn-1,Q7=F/mn-1,q_extra=z,Iphase=np.zeros(len(nn),int),
          Yn=z,Yp=z,Ye=z,Ymu=z,Aav=z,Zav=z,Yh=z,Il=int(hadronic.Il),iT=0,iYq=1,inb_start=1,
          temperature_mev=0,charge_fraction_grid=0,q_extra_count=0,composition_available=False,beta_equilibrated=True)
    low=make_phase(nlow,plow,elow,mulow,q4low,q5low)
    nmin=nq/(1+alpha); npre=np.geomspace(nmin,nq,max(8,n_points//3)); npost=np.geomspace(nq,nmax,max(8,n_points-len(npre)+1))[1:]
    nh=np.r_[npre,npost]; ratio=nh/nq; muh=mut*ratio**cs
    ph=pt+cs/(1+cs)*(epsq+pt)*(ratio**(1+cs)-1); eh=muh*nh-ph
    if np.any(ph<=0) or np.any(eh<=0): raise ValueError('CSS high phase became non-positive')
    q4h=np.full(len(nh),float(q4sp(lognt))/float(hadronic.mn)); q5h=np.full(len(nh),float(q5sp(lognt))/float(hadronic.mn))
    high=make_phase(nh,ph,eh,muh,q4h,q5h)
    result=build_maxwell_eos(low,high,low_name=hadronic_name,high_name=high_phase_name)
    result.metadata.update({'high_phase_model':'constant_sound_speed','cs2_high':cs,
      'input_transition_pressure_MeV_fm3':pt,'input_delta_epsilon_MeV_fm3':de,
      'input_epsilon_hadronic_transition_MeV_fm3':et,'input_n_hadronic_transition_fm3':nt,
      'input_muB_transition_MeV':mut,'n_high_phase_transition_fm3':nq,
      'metastable_overlap_fraction':alpha,'hadronic_mu_match_relative_residual':mu_residual,
      'charge_potential_extension':'mu_Q and effective lepton chemical potential are held at the transition values as a one-fluid placeholder; no high-phase microscopic composition is claimed.',
      'interpretation':'Phenomenological CSS high-density phase. Its transition pressure, energy jump and sound speed are model parameters, not a microscopic quark-matter calculation.'})
    return result
