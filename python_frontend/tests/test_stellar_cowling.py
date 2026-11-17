"""Physical controls for the dedicated Cowling reference."""
import numpy as np
import pytest
from compose_eos.stellar_cowling import (
    CowlingSettings,scaled_rhs_radius,solve_mode,scan_brackets,C_KM_S)


class UniformDensityControl:
    """Exact Schwarzschild interior, incompressible perturbations.

    This control is deliberately acausal: only its weak-gravity Kelvin/Cowling
    limit omega²=ell GM/R³ is used, never as a physical EOS.
    """
    def __init__(self,compactness=1e-4,R=10.,r0=.01):
        self.C=compactness;self.R=R
        self.s=np.sqrt(1-2*self.C)
        self.h_c=np.log(2*self.s/(3*self.s-1))
        self.h_start=self.enthalpy(r0)
    def enthalpy(self,r):
        return np.log(2*self.s/(3*self.s-np.sqrt(1-2*self.C*(r/self.R)**2)))
    def state_at_enthalpy(self,h):
        z=self.s*(3-2*np.exp(-h))
        r=self.R*np.sqrt((1-z*z)/(2*self.C))
        ep=(3*self.s-z)/2
        return dict(r_km=r,Phi=np.log(ep),lambda_metric=-np.log(z),
                    dPhi_dr_km=self.C*r/(self.R**2*z*ep),cs2_over_c2=np.inf)


def test_weak_gravity_incompressible_cowling_frequency():
    bg=UniformDensityControl()
    f0=C_KM_S*np.sqrt(2*bg.C/bg.R**2)/(2*np.pi)
    settings=CowlingSettings(surface_h=1e-10,rtol=1e-9,atol=1e-12,
                             max_step_h=1e-6)
    result=solve_mode(bg,[.99*f0,1.01*f0],settings,200)
    assert abs(result["frequency_hz"]/f0-1)<5e-4
    assert result["W_nodes"]==0
    assert abs(result["matching_residual"])<1e-7
    assert result["relative_join_mismatch"]<1e-7


def test_scaled_rhs_recovers_unscaled_covariant_equations():
    s=dict(r_km=7.,Phi=-.2,lambda_metric=.1,dPhi_dr_km=.03,cs2_over_c2=.4)
    r=s["r_km"]; w,v=.7,-.3;ell=2;k=.06
    dp=scaled_rhs_radius(s,[w,v],k,ell)
    W=r**3*w;V=r**2*v
    Wp=3*r*r*w+r**3*dp[0]
    Vp=2*r*v+r*r*dp[1]
    Q=k*k*r*r*np.exp(.5)*V+.03*W
    assert Wp==pytest.approx(Q/.4-6*np.exp(.1)*V,rel=2e-14)
    assert Vp==pytest.approx(.06*V-np.exp(.1)*W/r**2,rel=2e-14)


def test_invalid_scans_and_surface_domains_rejected():
    bg=UniformDensityControl()
    with pytest.raises(ValueError,match="strictly increasing"):
        scan_brackets(bg,[1.,1.])
    with pytest.raises(ValueError,match="surface_h"):
        CowlingSettings(surface_h=1.).validate(bg)
    with pytest.raises(ValueError,match="ell"):
        CowlingSettings(ell=0).validate(bg)
