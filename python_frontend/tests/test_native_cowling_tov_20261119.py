"""Native Fortran enthalpy TOV profile, then native Cowling modes."""
import json
from pathlib import Path

import numpy as np
import pytest

from compose_eos.stellar_background import ContractEnthalpyEOS, solve_background
from compose_eos.stellar_cowling import CowlingSettings, solve_mode
from compose_eos.stellar_cowling_tov_native import solve_background_native

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"data/core_support_20261117/primary_v1/REPRESENTATION.json"
CANONICAL=ROOT/"data/stellar_consumers_20261117/tov/CANONICAL_1P4.json"


@pytest.fixture(scope="module")
def eos():
    return ContractEnthalpyEOS.from_path(SOURCE)


@pytest.fixture(scope="module")
def comparison(eos):
    nc=json.loads(CANONICAL.read_text())["central_density_fm3"]
    kw=dict(surface_h=1e-4, rtol=1e-10,atol=1e-12,max_step_h=.002)
    return solve_background(eos,nc,**kw), solve_background_native(
        eos,nc,native_samples=8193,**kw)


def test_tov_native_matches_exact_smooth_star(comparison):
    old,new=comparison
    assert abs(new.mass_msun-old.mass_msun)/old.mass_msun<1e-5
    assert abs(new.radius_km-old.radius_km)/old.radius_km<1e-5
    assert new.mass_geom_km>0 and 2*new.mass_geom_km < new.radius_km
    assert new.backend=="Fortran DVODE Adams TOV"
    for f in (0.1,0.3,0.5,0.8):
        h=old.h_start*f
        o=old.state_at_enthalpy(h)
        n=new.state_at_enthalpy(h)
        for field in ("r_km","m_geom_km","Phi","lambda_metric","dPhi_dr_km"):
            assert np.isclose(o[field],n[field],rtol=5e-5,atol=2e-7),(field,h,o[field],n[field])


def test_native_tov_cowling_f_and_p1(comparison):
    old,new=comparison
    settings=CowlingSettings(surface_h=1e-4,rtol=1e-9,atol=1e-12)
    for label,br,nodes,reference in (("f",[2250,2500],0,2401.381874932175),
                                     ("p1",[6500,6750],1,6721.957056883109)):
        result=solve_mode(new,br,settings,240,backend="fortran_dvode",
                          native_samples=8193)
        assert abs(result["frequency_hz"]-reference)<0.03,(label,result["frequency_hz"])
        assert result["W_nodes"]==nodes
        assert result["V_nodes"]==nodes


def test_native_tov_rejects_unphysical_inputs(eos):
    for n in (-1.,0.,2.):
        with pytest.raises(ValueError):
            solve_background_native(eos,n,surface_h=1e-4)
    with pytest.raises(ValueError):
        solve_background_native(eos,.5,native_samples=10)


def test_public_background_backend_opt_in(eos):
    x=solve_background(eos,.54,surface_h=1e-4,backend="fortran_dvode",
                       native_samples=4097)
    assert x.mass_msun>0
    assert x.backend=="Fortran DVODE Adams TOV"


class UniformDensityExactTOV:
    """Acausal analytic Schwarzschild interior as a numerical limit test."""
    def __init__(self, epsilon=300.):
        self.epsilon=epsilon

    def enthalpy_at_density(self,n):
        return float(n)

    def state_at_enthalpy(self,h):
        return dict(p=self.epsilon*np.expm1(h),epsilon=self.epsilon,
                    cs2_over_c2=1.e100)


def test_native_tov_exact_schwarzschild_interior():
    from compose_eos.stellar_background import MEV_FM3_TO_GEOM_KM2
    eos=UniformDensityExactTOV()
    hc=0.25
    ratio=np.expm1(hc)
    C=0.5*(1-((1+ratio)/(1+3*ratio))**2)
    radius=np.sqrt(3*C/(4*np.pi*eos.epsilon*MEV_FM3_TO_GEOM_KM2))
    native=solve_background_native(eos,hc,native_samples=4097,
                                   surface_h=0.,rtol=1e-11,
                                   atol=1e-13,max_step_h=.001)
    assert abs(native.radius_km/radius-1)<2e-8
    assert abs(native.mass_geom_km/(C*radius)-1)<2e-8


@pytest.mark.parametrize("density",[0.35,0.54,0.7])
def test_native_tov_mass_radius_multi_density(eos,density):
    opts=dict(surface_h=1e-4,rtol=1e-10,atol=1e-12,max_step_h=.002)
    py=solve_background(eos,density,**opts)
    native=solve_background_native(eos,density,native_samples=4097,**opts)
    assert abs(native.mass_msun-py.mass_msun)<1.e-5
    assert abs(native.radius_km-py.radius_km)<1.e-5


def test_native_tov_extension_high_density():
    path=ROOT/"data/stellar_extension_20261117/eos/fine/REPRESENTATION.json"
    e=ContractEnthalpyEOS.from_path(path)
    opts=dict(surface_h=1e-4,rtol=1e-10,atol=1e-12,max_step_h=.002)
    a=solve_background(e,1.15,**opts)
    b=solve_background_native(e,1.15,native_samples=4097,**opts)
    assert abs(a.mass_msun-b.mass_msun)<2e-5
    assert abs(a.radius_km-b.radius_km)<1e-5


def test_native_tov_supports_vacuum_surface(eos):
    result=solve_background_native(eos,.54,native_samples=4097,surface_h=0.)
    assert result.state_at_enthalpy(0.)["p"]==0.
    assert result.state_at_enthalpy(0.)["epsilon"]==0.


def test_native_tov_rejects_sharp_EOS_before_sampling():
    class FakeSharp:
        sharp_interfaces=True
    class FakeEOS:
        barotrope=FakeSharp()
    with pytest.raises(ValueError,match="Maxwell"):
        solve_background_native(FakeEOS(),.5)


def test_native_tov_enthalpy_ode_and_surface_metric(comparison):
    _, star = comparison
    from compose_eos.stellar_background import MEV_FM3_TO_GEOM_KM2
    surface=star.state_at_enthalpy(star.surface_h)
    expected=1-2*star.mass_geom_km/star.radius_km
    assert np.isclose(np.exp(2*surface["Phi"]),expected,rtol=3e-14)
    for frac in (.22,.43,.71):
        h=frac*star.h_start
        s=np.sqrt(star.h_start-h)
        r,m=star._tov_pchip(s)
        r_s,m_s=star._tov_pchip.derivative()(s)
        st=star.eos.state_at_enthalpy(h)
        p=st["p"]*MEV_FM3_TO_GEOM_KM2
        e=st["epsilon"]*MEV_FM3_TO_GEOM_KM2
        expected_dr=-r*(r-2*m)/(m+4*np.pi*r**3*p)
        observed_dr=-r_s/(2*s)
        expected_dm=4*np.pi*r*r*e*expected_dr
        observed_dm=-m_s/(2*s)
        assert abs(observed_dr/expected_dr-1)<2e-4,(frac,"r",observed_dr,expected_dr)
        assert abs(observed_dm/expected_dm-1)<3e-4,(frac,"m",observed_dm,expected_dm)


def test_native_cowling_two_backgrounds_no_cross_contamination(eos):
    from compose_eos.stellar_cowling_native import NativeCowlingDVODE
    kw=dict(surface_h=1e-4,rtol=1e-10,atol=1e-12,max_step_h=.002,
            native_samples=4097)
    a=solve_background_native(eos,.35,**kw)
    b=solve_background_native(eos,.7,**kw)
    settings=CowlingSettings(surface_h=1e-4,rtol=1e-9,atol=1e-12)
    x=NativeCowlingDVODE(a,settings,samples=4097)
    y=NativeCowlingDVODE(b,settings,samples=4097)
    f=2040.
    first=x.matching_residual(f)
    _=y.matching_residual(2520.)
    second=x.matching_residual(f)
    assert first==second
