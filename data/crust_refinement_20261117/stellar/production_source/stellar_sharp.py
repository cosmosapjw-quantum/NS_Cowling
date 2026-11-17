"""One-sided sharp-interface TOV and idealized slow/rapid Cowling reference.

Metric (-,+,+,+); r,m in km, thermodynamics in MeV/fm^3, omega=2*pi*f/c.
The underlying piece laws and source rows remain unchanged. Bulk closure is
Gamma_1=Gamma_eq. Interface conversion is an explicit modeling choice.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from .stellar_background import (
    ContractEnthalpyEOS, StellarBackground, solve_background,
    MEV_FM3_TO_GEOM_KM2,
)
from .stellar_cowling import CowlingSettings, C_KM_S, scaled_rhs_radius, _nodes


class SharpContractEnthalpyEOS(ContractEnthalpyEOS):
    """Chemical-potential coordinate with two states at each coexistence point."""
    def __init__(self, barotrope):
        b = self.barotrope = barotrope
        self.a = b.surface.rest_energy_MeV
        self.h_max = float(np.log(b.upper_mu[-1]/self.a))
        self.density_domain, self.pressure_domain = b.density_domain, b.pressure_domain
        self.h_surface_anchor = float(b.surface.enthalpy(b.surface.n0))
        self.h_outer_end = float(np.log(b.outer.mu[-1]/self.a))
        self.h_bridge_end = float(np.log(b.bridge.mu_right/self.a))
        self._mu_search = b.upper_mu.copy()
        self.transitions = []
        for i in sorted(b.sharp_interfaces):
            mus = b.upper_mu[i:i+2]
            h = float(np.log(mus[0]/self.a))
            hs = np.log(mus/self.a)
            pad = 8*np.finfo(float).eps*max(1.,abs(h))
            self.transitions.append(dict(index=i,h=h,
                h_band=(float(min(hs)-pad),float(max(hs)+pad)),
                p=float(b.upper_p[i]),n_low=float(b.upper_n[i]),
                n_high=float(b.upper_n[i+1])))
            self._mu_search[i+1] = mus[0]  # search coordinate only
        self.transition_h = np.array([t["h"] for t in self.transitions])
        self.break_h = np.unique(np.log(np.r_[self.a,b.outer.mu,
                                              b.bridge.mu_right,self._mu_search]/self.a))

    def _transition(self,h):
        return next((t for t in self.transitions
                     if t["h_band"][0] <= h <= t["h_band"][1]),None)

    def state_at_enthalpy(self,h,side="high"):
        h=float(h)
        if side not in ("low","high"):
            raise ValueError("side must be low or high")
        if not np.isfinite(h) or h < 0 or h > self.h_max:
            raise ValueError("enthalpy outside actual contract support")
        if not self.transitions or h <= self.h_bridge_end:
            return ContractEnthalpyEOS.state_at_enthalpy(self,h)
        t=self._transition(h)
        if t is not None:
            s=self.barotrope._sharp_state(t["index"],t["p"],side)
        else:
            mu=self.a*np.exp(h)
            i=int(np.clip(np.searchsorted(self._mu_search,mu,side="right")-1,
                          0,len(self.barotrope.upper_intervals)-1))
            interval=self.barotrope.upper_intervals[i]
            if interval is None:
                raise ValueError("density gap selected without one-sided interface state")
            s=interval.evaluate_mu(float(np.clip(mu,interval.mu_left,interval.mu_right)))
        s["enthalpy"]=h
        return s

    def enthalpy_at_density(self,n):
        n=float(n)
        for t in self.transitions:
            if t["n_low"] < n < t["n_high"]:
                raise ValueError("density gap has no homogeneous stellar state")
            if n in (t["n_low"],t["n_high"]):
                return t["h"]
        return super().enthalpy_at_density(n)

    def side_for_segment(self,h,midpoint):
        t=self._transition(float(h))
        return "low" if t is not None and midpoint < t["h"] else "high"


class _PiecewiseSolution:
    """Dense outputs remain inside their own integrated phase segments."""
    def __init__(self,segments,jumps=()):
        self.segments=segments
        self.jumps=list(jumps)
        self.y=np.concatenate([s.y for s in segments],axis=1)
        self.nfev=sum(s.nfev for s in segments)
        self.success=all(s.success for s in segments)

    def sol(self,hs,side="high"):
        values=np.asarray(hs,dtype=float)
        out=[]
        for h in values.reshape(-1):
            eligible=[s for s in self.segments if min(s.t[0],s.t[-1]) <= h <= max(s.t[0],s.t[-1])]
            if not eligible:
                raise ValueError("dense query outside integrated segments")
            s=(max if side=="high" else min)(eligible,key=lambda x:.5*(x.t[0]+x.t[-1]))
            out.append(s.sol(float(h)))
        arr=np.asarray(out).T.reshape((2,)+values.shape)
        return arr


@dataclass
class SharpStellarBackground(StellarBackground):
    interfaces: list = field(default_factory=list)

    def state_at_enthalpy(self,h,side="high"):
        h=float(h)
        if h < self.surface_h or h > self.h_start:
            raise ValueError("background query outside integrated radial domain")
        r,m=map(float,self.solution.sol(h))
        s=self.eos.state_at_enthalpy(h,side)
        p=s["p"]*MEV_FM3_TO_GEOM_KM2
        e=s["epsilon"]*MEV_FM3_TO_GEOM_KM2
        phi=.5*np.log1p(-2*self.mass_geom_km/self.radius_km)-h+self.surface_h
        grad=(m+4*np.pi*r**3*p)/(r*(r-2*m))
        return dict(**s,r_km=r,m_geom_km=m,Phi=float(phi),
                    lambda_metric=float(-.5*np.log1p(-2*m/r)),
                    dPhi_dr_km=float(grad),p_geom_km2=p,epsilon_geom_km2=e)

    def summary(self):
        out=super().summary()
        out.update(backend="Python piece-contract sharp enthalpy DOP853 reference",
                   interfaces=self.interfaces)
        return out


def solve_sharp_background(eos,central_density_fm3,*,r0_km=.001,rtol=1e-9,
                           atol=1e-11,surface_h=0.,max_step_h=.005):
    controls=np.array([r0_km,rtol,atol,max_step_h,surface_h],float)
    if np.any(~np.isfinite(controls)) or min(controls[:4])<=0 or surface_h<0:
        raise ValueError("positive finite solver controls required")
    if not eos.transitions:
        smooth=solve_background(eos,central_density_fm3,r0_km=r0_km,rtol=rtol,
                                atol=atol,surface_h=surface_h,max_step_h=max_step_h)
        return SharpStellarBackground(**smooth.__dict__)
    hc=eos.enthalpy_at_density(central_density_fm3)
    if eos._transition(hc) is not None:
        raise ValueError("central density must lie inside a phase, not at its interface")
    st=eos.state_at_enthalpy(hc)
    pc,ec=st["p"]*MEV_FM3_TO_GEOM_KM2,st["epsilon"]*MEV_FM3_TO_GEOM_KM2
    h0=hc-2*np.pi*(ec/3+pc)*r0_km**2
    if not surface_h < h0 < hc:
        raise ValueError("center or surface controls leave no stellar interior")
    y=np.array([r0_km,4*np.pi*ec*r0_km**3/3])
    crossings=sorted([h for h in eos.transition_h if surface_h<h<h0],reverse=True)
    bounds=[h0,*crossings,float(surface_h)]
    segments=[];interfaces=[]
    for start,end in zip(bounds[:-1],bounds[1:]):
        midpoint=.5*(start+end)
        def rhs(h,y):
            r,m=y
            s=eos.state_at_enthalpy(float(np.clip(h,0.,hc)),eos.side_for_segment(h,midpoint))
            p,e=s["p"]*MEV_FM3_TO_GEOM_KM2,s["epsilon"]*MEV_FM3_TO_GEOM_KM2
            if r<=0 or r<=2*m:
                raise ValueError("nonpositive radius or trapped TOV integration state")
            dr=-r*(r-2*m)/(m+4*np.pi*r**3*p)
            return [dr,4*np.pi*r*r*e*dr]
        sol=solve_ivp(rhs,(start,end),y,method="DOP853",rtol=rtol,atol=atol,
                      max_step=max_step_h,dense_output=True)
        if not sol.success:
            raise RuntimeError("sharp TOV integration failed: "+sol.message)
        segments.append(sol);y=sol.y[:,-1].copy()
        if end in crossings:
            low=eos.state_at_enthalpy(end,"low");high=eos.state_at_enthalpy(end,"high")
            interfaces.append(dict(h=float(end),r_km=float(y[0]),m_geom_km=float(y[1]),
                r_jump_km=0.,m_jump_km=0.,p_MeV_fm3=low["p"],
                n_low_fm3=low["n"],n_high_fm3=high["n"],
                epsilon_low_MeV_fm3=low["epsilon"],epsilon_high_MeV_fm3=high["epsilon"]))
    r,m=map(float,y)
    if r<=2*m or m<=0 or r<=r0_km:
        raise RuntimeError("unphysical endpoint in sharp TOV reference")
    cfg=dict(r0_km=r0_km,rtol=rtol,atol=atol,surface_h=surface_h,max_step_h=max_step_h)
    return SharpStellarBackground(eos,float(central_density_fm3),hc,h0,float(surface_h),
                                  r,m,_PiecewiseSolution(segments),cfg,interfaces)


def transfer_interface(y,source,target,omega_geom,conversion):
    """Preserve DeltaP and w (slow) or v (rapid) with no surface layer.

    DeltaP=-(epsilon+P)e^-Lambda r^(ell-1) (A*v+B*w),
    A=omega^2 e^(Lambda-2Phi) r, B=Phi'. Rapid conversion also preserves
    xi^r-DeltaP/P0', where P0'=-(epsilon+P)B.
    """
    if conversion not in ("slow","rapid"):
        raise ValueError("conversion must be explicitly slow or rapid")
    if not np.isfinite(omega_geom) or omega_geom<=0:
        raise ValueError("positive finite frequency required")
    y=np.asarray(y,dtype=float)
    if y.shape!=(2,) or not np.all(np.isfinite(y)):
        raise ValueError("finite two-component mode state required")
    for key in ("r_km","Phi","lambda_metric","dPhi_dr_km","p"):
        if not np.isclose(source[key],target[key],rtol=4e-12,atol=1e-14):
            raise ValueError("interface states must share pressure and geometry")
    a=omega_geom**2*np.exp(source["lambda_metric"]-2*source["Phi"])*source["r_km"]
    b=source["dPhi_dr_km"]
    q=(source["epsilon"]+source["p"])/(target["epsilon"]+target["p"])
    if not np.isfinite(a+b+q) or min(a,b,q)<=0:
        raise ValueError("positive traction coefficients and enthalpy densities required")
    w,v=y
    if conversion=="slow":
        return np.array([w,q*v+(q-1)*b/a*w])
    return np.array([q*w+(q-1)*a/b*v,v])


def _integrate_mode_path(background,start,end,y,om,settings,conversion,dense):
    direction=np.sign(end-start)
    crossings=sorted([h for h in background.eos.transition_h if min(start,end)<h<max(start,end)],
                     reverse=direction<0)
    bounds=[start,*crossings,end];segments=[];jumps=[]
    for a,b in zip(bounds[:-1],bounds[1:]):
        midpoint=.5*(a+b)
        def rhs(h,yy):
            side=background.eos.side_for_segment(h,midpoint)
            st=background.state_at_enthalpy(h,side)
            return -scaled_rhs_radius(st,yy,om,settings.ell)/st["dPhi_dr_km"]
        sol=solve_ivp(rhs,(a,b),y,method="DOP853",rtol=settings.rtol,
                      atol=settings.atol,max_step=settings.max_step_h,dense_output=dense)
        if not sol.success:
            raise RuntimeError("sharp Cowling integration failed: "+sol.message)
        segments.append(sol);y=sol.y[:,-1].copy()
        if b in crossings:
            source_side,target_side=("high","low") if direction<0 else ("low","high")
            ss=background.state_at_enthalpy(b,source_side)
            ts=background.state_at_enthalpy(b,target_side)
            after=transfer_interface(y,ss,ts,om,conversion)
            aa=om**2*np.exp(ss["lambda_metric"]-2*ss["Phi"])*ss["r_km"]
            bb=ss["dPhi_dr_km"]
            q1=(ss["epsilon"]+ss["p"])*(aa*y[1]+bb*y[0])
            q2=(ts["epsilon"]+ts["p"])*(aa*after[1]+bb*after[0])
            jumps.append(dict(h=float(b),r_km=ss["r_km"],source_side=source_side,
                target_side=target_side,before=y.tolist(),after=after.tolist(),
                traction_relative_residual=float(abs(q1-q2)/max(abs(q1),abs(q2),1e-300)),
                kinematic_residual=float(after[0]-y[0] if conversion=="slow" else after[1]-y[1])))
            y=after
    return _PiecewiseSolution(segments,jumps)


def _integrations(background,frequency_hz,settings,conversion,dense=False):
    settings.validate(background)
    if conversion not in ("slow","rapid"):
        raise ValueError("conversion must be explicitly slow or rapid")
    if not np.isfinite(frequency_hz) or frequency_hz<=0:
        raise ValueError("frequency_hz must be positive and finite")
    if settings.surface_h<background.surface_h:
        raise ValueError("mode surface must lie in integrated background")
    om=2*np.pi*frequency_hz/C_KM_S
    hm=settings.match_h_fraction*background.h_c
    if background.eos._transition(hm) is not None:
        raise ValueError("choose a matching point away from a sharp interface")
    inside=_integrate_mode_path(background,background.h_start,hm,
                               [1.,-1./settings.ell],om,settings,conversion,dense)
    st=background.state_at_enthalpy(settings.surface_h)
    v=-st["dPhi_dr_km"]/(om**2*np.exp(st["lambda_metric"]-2*st["Phi"])*st["r_km"])
    outside=_integrate_mode_path(background,settings.surface_h,hm,[1.,v],
                                om,settings,conversion,dense)
    return inside,outside


def matching_residual(background,frequency_hz,settings=CowlingSettings(),*,conversion):
    a,b=_integrations(background,frequency_hz,settings,conversion)
    u,v=a.y[:,-1],b.y[:,-1]
    return float(np.linalg.det(np.column_stack((u,v)))/(np.linalg.norm(u)*np.linalg.norm(v)))


def scan_brackets(background,frequencies_hz,settings=CowlingSettings(),*,conversion):
    fs=np.asarray(frequencies_hz,dtype=float)
    if fs.ndim!=1 or len(fs)<2 or np.any(~np.isfinite(fs)) or np.any(np.diff(fs)<=0) or fs[0]<=0:
        raise ValueError("scan frequencies must be finite, increasing and positive")
    values=[matching_residual(background,f,settings,conversion=conversion) for f in fs]
    brackets=[(float(fs[i]),float(fs[i+1])) for i in range(len(fs)-1) if values[i]*values[i+1]<0]
    return brackets,[dict(frequency_hz=float(f),residual=float(v)) for f,v in zip(fs,values)]


def solve_mode(background,bracket_hz,settings=CowlingSettings(),*,conversion,profile_points=1200):
    lo,hi=map(float,bracket_hz)
    if not 0<lo<hi or profile_points<20:
        raise ValueError("positive ordered frequency bracket and >=20 profile points required")
    root,info=brentq(lambda f:matching_residual(background,f,settings,conversion=conversion),
                    lo,hi,xtol=settings.root_xtol_hz,rtol=4*np.finfo(float).eps,full_output=True)
    inside,outside=_integrations(background,root,settings,conversion,dense=True)
    hm=settings.match_h_fraction*background.h_c
    ui,uo=inside.y[:,-1],outside.y[:,-1]
    scale=np.dot(ui,uo)/np.dot(uo,uo)
    hgrid=np.r_[np.linspace(background.h_start,hm,profile_points//2),
                np.geomspace(settings.surface_h,hm,profile_points//2)[::-1][1:]]
    points=[(float(h),"high") for h in hgrid if background.eos._transition(h) is None]
    for h in background.eos.transition_h:
        if settings.surface_h<h<background.h_start:
            points.extend([(float(h),"high"),(float(h),"low")])
    points.sort(key=lambda x:(-x[0],0 if x[1]=="high" else 1))
    profile=[];ws=[];vs=[]
    for h,side in points:
        y=inside.sol(h,side) if h>=hm else outside.sol(h,side)*scale
        st=background.state_at_enthalpy(h,side);r=st["r_km"]
        W=float(r**(settings.ell+1)*y[0]);V=float(r**settings.ell*y[1])
        ws.append(W);vs.append(V)
        profile.append(dict(h=h,side=side,r_km=r,W=W,V_km_inverse=V))
    norm=max(abs(np.asarray(ws)))
    for row in profile:
        row["W"]/=norm;row["V_km_inverse"]/=norm
    st=background.state_at_enthalpy(settings.surface_h)
    ys=outside.sol(settings.surface_h)*scale;om=2*np.pi*root/C_KM_S
    terms=np.array([om**2*np.exp(st["lambda_metric"]-2*st["Phi"])*st["r_km"]*ys[1],
                    st["dPhi_dr_km"]*ys[0]])
    return dict(frequency_hz=float(root),ell=settings.ell,conversion=conversion,
                bracket_hz=[lo,hi],root_converged=bool(info.converged),root_iterations=int(info.iterations),
                matching_residual=matching_residual(background,root,settings,conversion=conversion),
                relative_join_mismatch=float(np.linalg.norm(ui-scale*uo)/np.linalg.norm(ui)),
                surface_BC_relative=float(abs(sum(terms))/sum(abs(terms))),
                W_sign_changes_including_interfaces=_nodes(ws),
                V_sign_changes_including_interfaces=_nodes(vs),
                interface_junctions=inside.jumps+outside.jumps,profile=profile,
                nfev=inside.nfev+outside.nfev)
