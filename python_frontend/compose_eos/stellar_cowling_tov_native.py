"""Native enthalpy-coordinate TOV background for smooth Cowling calculations.

Physics TOV integration is Fortran DVODE, while Python supplies the original
smooth EOS contract and central-state initialization. All coordinate conversions,
surface normalization, and profile diagnostics are explicitly recorded.
"""
from __future__ import annotations

from dataclasses import dataclass

import ctypes as ct
import numpy as np
from scipy.interpolate import PchipInterpolator

from .stellar_background import (G_SI, C_SI, M_SUN_KG, MSUN_GEOM_KM,
                                 MEV_FM3_TO_GEOM_KM2)
from .stellar_cowling_native import _LOCK, _native_library


@dataclass
class NativeTOVBackground:
    eos: object
    central_density_fm3: float
    h_c: float
    h_start: float
    surface_h: float
    radius_km: float
    mass_geom_km: float
    s: np.ndarray
    profile: np.ndarray
    thermo: np.ndarray
    nfev: int
    config: dict

    backend: str = "Fortran DVODE Adams TOV"

    def __post_init__(self):
        self._tov_pchip = PchipInterpolator(self.s, self.profile, axis=0,
                                            extrapolate=False)
        self._thermo_pchip = PchipInterpolator(self.s, self.thermo, axis=0,
                                               extrapolate=False)
        self.break_h = np.asarray(getattr(self.eos, "break_h", []), dtype=float)

    @property
    def mass_msun(self):
        return self.mass_geom_km/MSUN_GEOM_KM

    def _sample_native(self, s):
        """Get r,m,Pgeom,epsgeom,cs2 without calling the Python EOS."""
        query = np.asarray(s, dtype=float)
        if np.any(~np.isfinite(query)) or np.any(query < -1e-12) or np.any(query > self.s[-1]+1e-12):
            raise ValueError("background query outside native profile domain")
        r_m = self._tov_pchip(np.clip(query, self.s[0],self.s[-1]))
        therm = self._thermo_pchip(np.clip(query,self.s[0],self.s[-1]))
        return r_m,therm

    def coefficients_on_s(self, s):
        """Vectorized coefficient path directly into native Cowling integrator."""
        rm,t=self._sample_native(s)
        r,m=rm[...,0],rm[...,1]
        p=t[...,0]
        h=self.h_start-np.asarray(s)**2
        metric_factor=1.-2*m/r
        if np.any(metric_factor<=0):
            raise RuntimeError("trapped/nonphysical native TOV profile")
        phi=0.5*np.log1p(-2*self.mass_geom_km/self.radius_km)-h+self.surface_h
        lam=-0.5*np.log(metric_factor)
        gp=(m+4*np.pi*r**3*p)/(r*(r-2*m))
        return np.stack((r,phi,lam,gp,t[...,2]),axis=-1)

    def state_at_enthalpy(self, h):
        h=float(h)
        if not self.surface_h<=h<=self.h_start:
            raise ValueError("background query outside native TOV domain")
        s=np.sqrt(max(self.h_start-h,0.))
        r,m=self._tov_pchip(s)
        eos_state=self.eos.state_at_enthalpy(h)
        p=eos_state["p"]*MEV_FM3_TO_GEOM_KM2
        e=eos_state["epsilon"]*MEV_FM3_TO_GEOM_KM2
        lam=-0.5*np.log1p(-2*m/r)
        phi=0.5*np.log1p(-2*self.mass_geom_km/self.radius_km)-h+self.surface_h
        gp=(m+4*np.pi*r**3*p)/(r*(r-2*m))
        return dict(**eos_state,r_km=float(r),m_geom_km=float(m),
                    Phi=float(phi),lambda_metric=float(lam),
                    dPhi_dr_km=float(gp),p_geom_km2=float(p),
                    epsilon_geom_km2=float(e))

    def summary(self):
        return dict(central_density_fm3=self.central_density_fm3,
                    h_c=self.h_c,h_start=self.h_start,
                    surface_h=self.surface_h,
                    radius_km=self.radius_km,mass_msun=self.mass_msun,
                    mass_geom_km=self.mass_geom_km,
                    compactness=self.mass_geom_km/self.radius_km,
                    nfev=self.nfev,config=self.config,backend=self.backend)


def solve_background_native(eos, central_density_fm3, *, r0_km=0.001,
                            rtol=1e-9, atol=1e-11,surface_h=0.,
                            max_step_h=.005,native_samples=8193,
                            method="adams"):
    """TOV numerical evolution exclusively with Fortran DVODE.

    Source EOS evaluation & tabulation remain on the exact Python contract.
    The native profile uses s=sqrt(h_start-h) to regularize center sampling.
    """
    if method not in ("adams","bdf"):
        raise ValueError("native TOV DVODE method must be adams or bdf")
    if not isinstance(native_samples,int) or native_samples<33:
        raise ValueError("native_samples must be an integer >= 33")
    if (not np.all(np.isfinite([central_density_fm3,r0_km,rtol,atol,
                               surface_h,max_step_h]))
        or central_density_fm3<=0 or r0_km<=0 or rtol<=0
        or atol<=0 or max_step_h<=0 or surface_h<0):
        raise ValueError("positive finite native TOV controls required")
    if (bool(getattr(eos,"sharp_interfaces",False)) or
        bool(getattr(getattr(eos,"barotrope",None),"sharp_interfaces",False))):
        raise ValueError("native smooth TOV does not admit Maxwell jump interfaces")
    hc=float(eos.enthalpy_at_density(central_density_fm3))
    center=eos.state_at_enthalpy(hc)
    pc=center["p"]*MEV_FM3_TO_GEOM_KM2
    ec=center["epsilon"]*MEV_FM3_TO_GEOM_KM2
    hstart=hc-2*np.pi*(ec/3+pc)*r0_km*r0_km
    if not surface_h<hstart<hc:
        raise ValueError("center/surface TOV controls leave no stellar interior")
    r0=float(r0_km)
    m0=float(4*np.pi*ec*r0**3/3)
    smax=float(np.sqrt(hstart-surface_h))
    sgrid=np.linspace(0,smax,native_samples,dtype=np.float64)
    if surface_h>0.:
        hg=np.geomspace(surface_h,hstart,max(1025,native_samples//2))
        sgrid=np.unique(np.r_[sgrid,np.sqrt(np.maximum(hstart-hg,0.))])
    break_h=np.asarray(getattr(eos,"break_h",[]),dtype=float)
    break_h=break_h[(break_h>surface_h)&(break_h<hstart)]
    if len(break_h):
        break_s=np.sqrt(hstart-break_h)
        d=smax/(native_samples-1)
        # Explicit two-sided EOS derivatives at each piece-join.
        widths=np.minimum(0.01*d,np.maximum(1e-10,2e-12/np.maximum(break_s,1e-4)))
        x=np.r_[break_s-widths,break_s,break_s+widths]
        x=x[(x>0)&(x<smax)]
        sgrid=np.unique(np.r_[sgrid,x])
    hvalues=hstart-sgrid*sgrid
    hvalues[0]=hstart
    hvalues[-1]=surface_h
    thermo=np.empty((len(sgrid),3),dtype=np.float64)
    for i,h in enumerate(hvalues):
        state=eos.state_at_enthalpy(float(h))
        thermo[i,:]=[state["p"]*MEV_FM3_TO_GEOM_KM2,
                     state["epsilon"]*MEV_FM3_TO_GEOM_KM2,
                     state["cs2_over_c2"]]
    if not np.all(np.isfinite(thermo)) or np.any(thermo[:,1]<0):
        raise ValueError("native smooth TOV source has nonfinite/invalid thermodynamics")
    thermo_spline=PchipInterpolator(sgrid,thermo[:,:2],axis=0,
                                   extrapolate=False)
    slopes=thermo_spline(sgrid,nu=1)
    n=len(sgrid)
    result=np.empty((2,n),dtype=np.float64,order="F")
    nfev=ct.c_int()
    status=ct.c_int()
    mf={"adams":10,"bdf":22}[method]
    max_ds=min(smax/200.,max_step_h/(2*smax))
    with _LOCK:
        dll=_native_library()
        dll.tov_native_profile(n,np.ascontiguousarray(sgrid),
            np.ascontiguousarray(thermo[:,0]),
            np.ascontiguousarray(thermo[:,1]),
            np.ascontiguousarray(slopes[:,0]),np.ascontiguousarray(slopes[:,1]),
            r0,m0,rtol,atol,min(atol*1.e-3,1.e-15),max_ds,mf,
            result,ct.byref(nfev),ct.byref(status))
    if status.value:
        raise RuntimeError("Fortran TOV DVODE integration rejected: status "+
                           str(status.value))
    radius,mass=result[:,-1]
    if (radius<=2*mass or mass<=0 or radius<=r0
        or not np.isfinite(radius) or not np.isfinite(mass)):
        raise RuntimeError("nonphysical native TOV endpoint")
    config=dict(r0_km=r0_km,rtol=rtol,atol=atol,surface_h=surface_h,
                max_step_h=max_step_h,native_samples=native_samples,
                actual_samples=n,method=method,
                native_profile_coordinate="s=sqrt(h_start-h)",
                original_eos_contract_evaluated_in_python=True)
    return NativeTOVBackground(eos,central_density_fm3,hc,hstart,surface_h,
                               float(radius),float(mass),sgrid,result.T.copy(),
                               thermo,int(nfev.value),config,
                               backend="Fortran DVODE "+method.capitalize()+" TOV")
