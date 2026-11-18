"""Exact declared-contract TOV reference; geometry in km, thermodynamics MeV/fm^3.

Metric: ds²=-exp(2 Phi)c²dt²+(1-2m/r)^-1 dr²+r²dOmega².
m=G M/c²; h=int dP/(epsilon+P)=log(mu/mu_vac).
This reference does not replace the native Fortran/table backend.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from .whole_crust import WholeCrustBarotrope

C_SI = 299792458.
G_SI = 6.67430e-11
GM_SUN_NOMINAL = 1.3271244e20  # IAU nominal convention, matching native backend
M_SUN_KG = GM_SUN_NOMINAL/G_SI
MEV_FM3_TO_PA = 1.602176634e32
MEV_FM3_TO_GEOM_KM2 = G_SI*MEV_FM3_TO_PA/C_SI**4*1e6
MSUN_GEOM_KM = G_SI*M_SUN_KG/C_SI**2/1000.

class ContractEnthalpyEOS:
    """Evaluate the original piece laws directly in chemical potential.

    Sharp jumps are rejected here: this selected smooth baseline solver does
    not silently invent a discontinuity/interface prescription.
    """
    def __init__(self, barotrope):
        self.barotrope = barotrope
        if barotrope.sharp_interfaces:
            raise ValueError("smooth TOV reference requires a continuous contract")
        self.a = barotrope.surface.rest_energy_MeV
        self.h_max = float(np.log(barotrope.upper_mu[-1]/self.a))
        self.density_domain = barotrope.density_domain
        self.pressure_domain = barotrope.pressure_domain
        mu = np.r_[self.a, barotrope.outer.mu,
                   barotrope.bridge.mu_right, barotrope.upper_mu]
        self.break_h = np.unique(np.log(mu/self.a))
        self.h_surface_anchor = float(barotrope.surface.enthalpy(barotrope.surface.n0))
        self.h_outer_end = float(np.log(barotrope.outer.mu[-1]/self.a))
        self.h_bridge_end = float(np.log(barotrope.bridge.mu_right/self.a))

    @classmethod
    def from_path(cls, path):
        return cls(WholeCrustBarotrope.from_contract(json.loads(Path(path).read_text())))

    def state_at_enthalpy(self, h):
        h = float(h)
        if not np.isfinite(h) or h < 0 or h > self.h_max:
            raise ValueError("enthalpy outside actual contract support")
        b = self.barotrope
        mu = self.a*np.exp(h)
        if h <= self.h_surface_anchor:
            s = b.surface
            q = self.a*np.expm1(h)*(s.gamma-1)/s.gamma
            n = s.n0*(q/(s.p0/s.n0))**(1/(s.gamma-1))
            n = min(n,s.n0)  # only roundoff at the analytic endpoint
            state = s.evaluate(n)
        elif h <= self.h_outer_end:
            o = b.outer
            i = int(np.clip(np.searchsorted(o.mu,mu,side="right")-1,0,len(o.gamma)-1))
            gamma, q0 = o.gamma[i],o.p[i]/o.nb[i]
            delta = mu-o.mu[i]
            if abs(gamma-1) < 1e-10:
                x = delta/q0
            else:
                x = np.log1p(delta*(gamma-1)/(gamma*q0))/(gamma-1)
            n = o.nb[i]*np.exp(x)
            p = o.p[i]*np.exp(gamma*x)
            state = dict(n=float(n),p=float(p),epsilon=float(mu*n-p),
                         mu=float(mu),cs2_over_c2=float(gamma*p/(mu*n)))
        else:
            if h <= self.h_bridge_end:
                interval = b.bridge
            else:
                i = int(np.clip(np.searchsorted(b.upper_mu,mu,side="right")-1,
                                0,len(b.upper_intervals)-1))
                interval = b.upper_intervals[i]
            mu = float(np.clip(mu,interval.mu_left,interval.mu_right))
            state = interval.evaluate_mu(mu)
        state["enthalpy"] = h
        return state

    def enthalpy_at_density(self,n):
        n = float(n)
        if not np.isfinite(n) or n < 0 or n > self.density_domain[1]:
            raise ValueError("central density outside actual contract support")
        if n == 0: return 0.
        if n == self.density_domain[1]: return self.h_max
        return brentq(lambda h:self.state_at_enthalpy(h)["n"]-n,
                      0.,self.h_max,xtol=2e-15,rtol=8*np.finfo(float).eps)

@dataclass
class StellarBackground:
    eos: ContractEnthalpyEOS
    central_density_fm3: float
    h_c: float
    h_start: float
    surface_h: float
    radius_km: float
    mass_geom_km: float
    solution: object
    config: dict

    @property
    def mass_msun(self):
        return self.mass_geom_km/MSUN_GEOM_KM

    @property
    def break_h(self):
        return self.eos.break_h

    def state_at_enthalpy(self,h):
        h = float(h)
        if h < self.surface_h or h > self.h_start:
            raise ValueError("background query outside integrated radial domain")
        r,m = map(float,self.solution.sol(h))
        s = self.eos.state_at_enthalpy(h)
        p,e = s["p"]*MEV_FM3_TO_GEOM_KM2,s["epsilon"]*MEV_FM3_TO_GEOM_KM2
        phi = .5*np.log1p(-2*self.mass_geom_km/self.radius_km)-h+self.surface_h
        grad = (m+4*np.pi*r**3*p)/(r*(r-2*m))
        return dict(**s,r_km=r,m_geom_km=m,Phi=float(phi),
                    lambda_metric=float(-.5*np.log1p(-2*m/r)),
                    dPhi_dr_km=float(grad),p_geom_km2=p,epsilon_geom_km2=e)

    def state_at_radius(self,r_km):
        r0 = float(self.solution.y[0,0])
        r = float(r_km)
        if r < r0 or r > self.radius_km:
            raise ValueError("radius outside integrated radial domain")
        h = brentq(lambda x:float(self.solution.sol(x)[0])-r,
                   self.surface_h,self.h_start,xtol=2e-15)
        return self.state_at_enthalpy(h)

    def summary(self):
        return dict(central_density_fm3=self.central_density_fm3,h_c=self.h_c,
                    h_start=self.h_start,surface_h=self.surface_h,
                    central_pressure_MeV_fm3=self.eos.state_at_enthalpy(self.h_c)["p"],
                    radius_km=self.radius_km,mass_msun=self.mass_msun,
                    mass_geom_km=self.mass_geom_km,
                    compactness=self.mass_geom_km/self.radius_km,
                    nfev=self.solution.nfev,config=self.config,
                    backend="Python exact-contract enthalpy DOP853 reference")

def solve_background(eos, central_density_fm3, *, r0_km=.001,rtol=1e-9,
                     atol=1e-11,surface_h=0.,max_step_h=.005,
                     backend="python",native_samples=8193,native_method="adams"):
    """Integrate smooth beta-barotropic TOV within the declared source support.

    m0=4 pi epsilon_c r0³/3, hc-h0=2 pi (epsilon_c/3+Pc)r0².
    O(r0^4) and O(r0^5) omitted terms are controlled by center-radius convergence.
    """
    if backend == "fortran_dvode":
        from .stellar_cowling_tov_native import solve_background_native
        return solve_background_native(eos,central_density_fm3,r0_km=r0_km,
                   rtol=rtol,atol=atol,surface_h=surface_h,
                   max_step_h=max_step_h,native_samples=native_samples,
                   method=native_method)
    if backend != "python":
        raise ValueError("unsupported TOV backend")
    if not (np.isfinite(r0_km) and r0_km > 0 and rtol > 0 and atol > 0
            and max_step_h > 0 and surface_h >= 0):
        raise ValueError("positive finite solver controls required")
    hc = eos.enthalpy_at_density(central_density_fm3)
    s = eos.state_at_enthalpy(hc)
    p,e = s["p"]*MEV_FM3_TO_GEOM_KM2,s["epsilon"]*MEV_FM3_TO_GEOM_KM2
    h0 = hc-2*np.pi*(e/3+p)*r0_km**2
    if not surface_h < h0 < hc:
        raise ValueError("center or surface controls leave no stellar interior")
    m0 = 4*np.pi*e*r0_km**3/3
    def rhs(h,y):
        r,m = y
        state = eos.state_at_enthalpy(float(np.clip(h,0.,hc)))
        p = state["p"]*MEV_FM3_TO_GEOM_KM2
        e = state["epsilon"]*MEV_FM3_TO_GEOM_KM2
        if r <= 0 or r <= 2*m:
            raise ValueError("nonpositive radius or trapped TOV integration state")
        dr = -r*(r-2*m)/(m+4*np.pi*r**3*p)
        return [dr,4*np.pi*r*r*e*dr]
    sol = solve_ivp(rhs,(h0,float(surface_h)),[r0_km,m0],method="DOP853",
                    rtol=rtol,atol=atol,max_step=max_step_h,dense_output=True)
    if not sol.success:
        raise RuntimeError("TOV integration failed: "+sol.message)
    r,m = map(float,sol.y[:,-1])
    if r <= 2*m or m <= 0 or r <= r0_km:
        raise RuntimeError("unphysical endpoint in TOV reference")
    config = dict(r0_km=r0_km,rtol=rtol,atol=atol,surface_h=surface_h,
                  max_step_h=max_step_h)
    return StellarBackground(eos,float(central_density_fm3),hc,h0,float(surface_h),
                             r,m,sol,config)
