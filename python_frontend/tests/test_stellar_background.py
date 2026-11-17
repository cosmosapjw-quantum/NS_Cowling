"""Scientific controls for the selected exact-contract TOV reference."""
from pathlib import Path
import numpy as np
import pytest
from compose_eos.stellar_background import (ContractEnthalpyEOS,solve_background,
    MEV_FM3_TO_GEOM_KM2,MSUN_GEOM_KM)

ROOT=Path(__file__).resolve().parents[2]
@pytest.fixture(scope="module")
def eos():
    return ContractEnthalpyEOS.from_path(ROOT/"data/core_support_20261117/primary_v1/REPRESENTATION.json")

def test_enthalpy_consumer_reuses_original_piece_laws(eos):
    # Include all regions and both sides of each representation boundary.
    hs=np.unique(np.r_[np.geomspace(1e-7,eos.h_max,300),
                       eos.break_h[1:]])
    for h in hs:
        s=eos.state_at_enthalpy(h)
        q=eos.barotrope.state_at_pressure(min(s["p"],eos.pressure_domain[1]))
        for key in ("n","p","epsilon","mu","cs2_over_c2"):
            # At an exact join cs2 has valid distinct one-sided limits; an
            # inverse rounded pressure can select the other side. Compare
            # derivatives only at interior points; all states also at joins.
            if key == "cs2_over_c2" and np.min(abs(eos.break_h-h)) < 1e-14:
                continue
            assert np.isclose(s[key],q[key],rtol=2e-7,atol=1e-20),(h,key,s[key],q[key])

def test_first_law_in_enthalpy(eos):
    for lo,hi in zip(eos.break_h[:-1],eos.break_h[1:]):
        if hi-lo < 1e-9: continue
        h=(hi+lo)/2;dh=(hi-lo)*1e-3
        a,b,c=(eos.state_at_enthalpy(v) for v in (h-dh,h,h+dh))
        deriv=(c["p"]-a["p"])/(2*dh)
        assert np.isclose(deriv,b["epsilon"]+b["p"],rtol=3e-5)

def test_closed_support_and_analytic_vacuum(eos):
    z=eos.state_at_enthalpy(0.)
    assert z["n"]==z["p"]==z["epsilon"]==0.
    assert eos.enthalpy_at_density(1.)==eos.h_max
    for n in (-1.,np.nextafter(1.,np.inf),np.nan):
        with pytest.raises(ValueError): eos.enthalpy_at_density(n)
    for h in (-1e-10,np.nextafter(eos.h_max,np.inf),np.inf):
        with pytest.raises(ValueError): eos.state_at_enthalpy(h)
    with pytest.raises(ValueError): solve_background(eos,1.001)
    with pytest.raises(ValueError): solve_background(eos,.5,r0_km=-1)

class IncompressibleControl:
    # Exact Schwarzschild interior control, deliberately not causal matter.
    def __init__(self,epsilon=300.):
        self.epsilon=epsilon
    def enthalpy_at_density(self,n): return float(n)
    def state_at_enthalpy(self,h):
        return dict(p=self.epsilon*np.expm1(h),epsilon=self.epsilon)

def test_tov_against_analytic_schwarzschild_interior():
    e=IncompressibleControl()
    hc=.25
    ratio=np.expm1(hc)
    compactness=.5*(1-((1+ratio)/(1+3*ratio))**2)
    radius=np.sqrt(3*compactness/(4*np.pi*e.epsilon*MEV_FM3_TO_GEOM_KM2))
    b=solve_background(e,hc,r0_km=.001,rtol=1e-11,atol=1e-13,max_step_h=.001)
    assert np.isclose(b.radius_km,radius,rtol=1e-8)
    assert np.isclose(b.mass_geom_km,compactness*radius,rtol=1e-8)

def test_surface_normalization_center_and_differential_equations(eos):
    b=solve_background(eos,.54,rtol=2e-10,atol=1e-12)
    surf=b.state_at_enthalpy(0.)
    assert np.isclose(np.exp(2*surf["Phi"]),1-2*b.mass_geom_km/b.radius_km,rtol=1e-14)
    assert surf["p"]==surf["epsilon"]==0.
    for frac in (.1,.4,.8):
        h=frac*b.h_start;dh=1e-7
        l,c,r=(b.state_at_enthalpy(v) for v in (h-dh,h,h+dh))
        dr=(r["r_km"]-l["r_km"])/(2*dh)
        dm=(r["m_geom_km"]-l["m_geom_km"])/(2*dh)
        assert np.isclose(dr,-1/c["dPhi_dr_km"],rtol=3e-6)
        assert np.isclose(dm/dr,4*np.pi*c["r_km"]**2*c["epsilon_geom_km2"],rtol=3e-6)
        assert np.isclose(b.state_at_radius(c["r_km"])["enthalpy"],h,rtol=1e-12)
