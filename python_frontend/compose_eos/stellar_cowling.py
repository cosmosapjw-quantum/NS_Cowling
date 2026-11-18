"""Exact-contract equilibrium-barotropic relativistic Cowling reference solver.

Geometric lengths are km, omega=2*pi*f/(c in km/s), Phi is half the temporal
metric potential. Gamma1=Gamma_eq; no frozen-composition or elastic claim.
Continuous EOS knots retain W,V. A sharp density jump is outside this solver.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

C_KM_S = 299792.458


@dataclass(frozen=True)
class CowlingSettings:
    ell: int = 2
    rtol: float = 1e-7
    atol: float = 1e-10
    surface_h: float = 1e-7
    match_h_fraction: float = 0.15
    max_step_h: float = 0.002
    root_xtol_hz: float = 1e-5

    def validate(self, background):
        if self.ell < 2:
            raise ValueError("This nonradial implementation requires ell >= 2")
        if not 0 < self.surface_h < self.match_h_fraction*background.h_c < background.h_start:
            raise ValueError("Require surface_h < match_h < h_start")
        if min(self.rtol,self.atol,self.max_step_h,self.root_xtol_hz) <= 0:
            raise ValueError("Tolerances and step cap must be positive")


def scaled_rhs_radius(state, y, omega_geom, ell=2):
    """Return d(w,v)/dr for W=r^(ell+1)w, V=r^ell v."""
    r=state["r_km"]; phi=state["Phi"]; gp=state["dPhi_dr_km"]
    el=np.exp(state["lambda_metric"]); cs=state["cs2_over_c2"]
    w,v=y
    if r <= 0 or cs <= 0:
        raise ValueError("Cowling RHS requires r>0 and positive sound speed")
    return np.array([(omega_geom**2*el*np.exp(-2*phi)*r*v+gp*w)/cs
                     -(ell+1)*(ell*el*v+w)/r,
                     (2*gp-ell/r)*v-el*w/r])


def _integrations(background,frequency_hz,settings,dense=False):
    settings.validate(background)
    if not np.isfinite(frequency_hz) or frequency_hz <= 0:
        raise ValueError("frequency_hz must be positive and finite")
    om=2*np.pi*frequency_hz/C_KM_S
    hm=settings.match_h_fraction*background.h_c
    def rhs(h,y):
        st=background.state_at_enthalpy(h)
        return -scaled_rhs_radius(st,y,om,settings.ell)/st["dPhi_dr_km"]
    kwargs=dict(method="DOP853",rtol=settings.rtol,atol=settings.atol,
                max_step=settings.max_step_h,dense_output=dense)
    inside=solve_ivp(rhs,(background.h_start,hm),
                     [1.,-1./settings.ell],**kwargs)
    st=background.state_at_enthalpy(settings.surface_h)
    v=-st["dPhi_dr_km"]/(om**2*np.exp(st["lambda_metric"]-2*st["Phi"])*st["r_km"])
    outside=solve_ivp(rhs,(settings.surface_h,hm),[1.,v],**kwargs)
    if not inside.success or not outside.success:
        raise RuntimeError("Cowling integration failed: "+inside.message+"; "+outside.message)
    return inside,outside


def matching_residual(background,frequency_hz,settings=CowlingSettings(), *,
                      backend="python", native_samples=8193, native_method="adams"):
    if backend == "fortran_dvode":
        from .stellar_cowling_native import NativeCowlingDVODE
        return NativeCowlingDVODE(background, settings, samples=native_samples,
                                 method=native_method).matching_residual(frequency_hz)
    if backend != "python":
        raise ValueError("unsupported Cowling backend")
    a,b=_integrations(background,frequency_hz,settings)
    u=a.y[:,-1]; v=b.y[:,-1]
    return float((u[0]*v[1]-u[1]*v[0])/(np.linalg.norm(u)*np.linalg.norm(v)))


def scan_brackets(background,frequencies_hz,settings=CowlingSettings(), *,
                  backend="python", native_samples=8193, native_method="adams"):
    if backend == "fortran_dvode":
        from .stellar_cowling_native import NativeCowlingDVODE
        return NativeCowlingDVODE(background, settings, samples=native_samples,
                                 method=native_method).scan_brackets(frequencies_hz)
    if backend != "python":
        raise ValueError("unsupported Cowling backend")
    fs=np.asarray(frequencies_hz,dtype=float)
    if fs.ndim != 1 or len(fs)<2 or np.any(np.diff(fs)<=0) or fs[0]<=0:
        raise ValueError("scan frequencies must be strictly increasing and positive")
    values=[matching_residual(background,f,settings) for f in fs]
    brackets=[(float(fs[i]),float(fs[i+1])) for i in range(len(fs)-1)
              if values[i]*values[i+1]<0]
    return brackets,[dict(frequency_hz=float(f),residual=float(v)) for f,v in zip(fs,values)]


def _nodes(a):
    a=np.asarray(a)
    valid=a[np.abs(a)>1e-9*np.max(abs(a))]
    return int(np.count_nonzero(valid[:-1]*valid[1:]<0))


def solve_mode(background,bracket_hz,settings=CowlingSettings(),profile_points=1200, *,
               backend="python", native_samples=8193, native_method="adams"):
    if backend == "fortran_dvode":
        from .stellar_cowling_native import NativeCowlingDVODE
        return NativeCowlingDVODE(background, settings, samples=native_samples,
                                 method=native_method).solve_mode(bracket_hz,profile_points)
    if backend != "python":
        raise ValueError("unsupported Cowling backend")
    lo,hi=map(float,bracket_hz)
    if not 0 < lo < hi:
        raise ValueError("positive ordered frequency bracket required")
    root,info=brentq(lambda f:matching_residual(background,f,settings),
                    lo,hi,xtol=settings.root_xtol_hz,rtol=4*np.finfo(float).eps,
                    full_output=True)
    inside,outside=_integrations(background,root,settings,dense=True)
    hm=settings.match_h_fraction*background.h_c
    ui=inside.y[:,-1]; uo=outside.y[:,-1]
    scale=np.dot(ui,uo)/np.dot(uo,uo)
    hi_grid=np.linspace(background.h_start,hm,profile_points//2)
    ho_grid=np.geomspace(settings.surface_h,hm,profile_points//2)[::-1][1:]
    hs=np.r_[hi_grid,ho_grid]
    ys=np.column_stack((inside.sol(hi_grid),outside.sol(ho_grid)*scale))
    states=[background.state_at_enthalpy(h) for h in hs]
    rs=np.array([s["r_km"] for s in states])
    W=rs**(settings.ell+1)*ys[0];V=rs**settings.ell*ys[1]
    norm=np.max(abs(W))
    profile=[dict(h=float(h),r_km=float(r),W=float(w/norm),V_km_inverse=float(v/norm))
             for h,r,w,v in zip(hs,rs,W,V)]
    residual=matching_residual(background,root,settings)
    reljoin=float(np.linalg.norm(ui-scale*uo)/np.linalg.norm(ui))
    st=states[-1];om=2*np.pi*root/C_KM_S
    terms=np.array([om**2*np.exp(st["lambda_metric"]-2*st["Phi"])*st["r_km"]*ys[1,-1],
                    st["dPhi_dr_km"]*ys[0,-1]])
    return dict(frequency_hz=float(root),ell=settings.ell,bracket_hz=[lo,hi],
                root_converged=bool(info.converged),root_iterations=int(info.iterations),
                matching_residual=residual,relative_join_mismatch=reljoin,
                surface_BC_relative=float(abs(sum(terms))/sum(abs(terms))),
                W_nodes=_nodes(W),V_nodes=_nodes(V),
                r_surface_mode_km=float(rs[-1]),profile=profile,
                nfev=inside.nfev+outside.nfev)
