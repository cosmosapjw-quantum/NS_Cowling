"""Physical junction conditions and an independent pressure-coordinate TOV control."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from compose_eos.whole_crust import SurfacePolytrope,PiecewisePolytrope,HermiteInterval,WholeCrustBarotrope
from compose_eos.stellar_background import solve_background,MEV_FM3_TO_GEOM_KM2
from compose_eos.stellar_cowling import matching_residual as smooth_residual,CowlingSettings
from compose_eos.stellar_sharp import (
    SharpContractEnthalpyEOS,solve_sharp_background,transfer_interface,matching_residual,
)

def contract(sharp=True):
    surface=SurfacePolytrope(1.,.1,10.)
    outer=PiecewisePolytrope([1.,2.],[.1,.4],10.)
    bridge=HermiteInterval(2.,.4,20.2,3.,.9,30.6)
    n=np.array([3.,4.,5.,6.]) if sharp else np.array([3.,4.,5.])
    p=np.array([.9,1.6,1.6,2.304]) if sharp else .1*n*n
    e=np.array([30.6,41.2,51.9,62.664]) if sharp else 9.9*n+p
    return dict(surface=surface.contract(),outer=outer.contract(),lower_bridge=bridge.contract(),
                upper=dict(n=n.tolist(),p=p.tolist(),epsilon=e.tolist(),sharp_interfaces=[1] if sharp else []))

@pytest.mark.parametrize("conversion",["slow","rapid"])
def test_transfer_preserves_physical_conditions_and_round_trip(conversion):
    state=dict(r_km=11.,Phi=-.2,lambda_metric=.3,dPhi_dr_km=.02,p=3.)
    source=dict(state,epsilon=170.)
    target=dict(state,epsilon=120.)
    y=np.array([.7,-.2]);om=.04
    z=transfer_interface(y,source,target,om,conversion)
    a=om**2*np.exp(state["lambda_metric"]-2*state["Phi"])*state["r_km"]
    b=state["dPhi_dr_km"]
    traction=lambda q,s:(s["epsilon"]+s["p"])*(a*q[1]+b*q[0])
    assert traction(y,source)==pytest.approx(traction(z,target),rel=3e-15,abs=1e-14)
    if conversion=="slow":
        assert z[0]==y[0]
    else:
        # xi^r-DeltaP/P' = -r^(l-1)e^-Lambda (A/B) v is continuous.
        assert z[1]==y[1]
    np.testing.assert_allclose(transfer_interface(z,target,source,om,conversion),y,rtol=2e-15,atol=2e-15)
    np.testing.assert_array_equal(transfer_interface(y,source,source,om,conversion),y)

def test_one_sided_contract_does_not_fill_density_gap():
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(contract()))
    h=eos.transition_h[0]
    lo=eos.state_at_enthalpy(h,side="low");hi=eos.state_at_enthalpy(h,side="high")
    assert lo["n"]==4. and hi["n"]==5.
    assert lo["p"]==hi["p"]==1.6
    assert lo["mu"]==pytest.approx(hi["mu"],rel=2e-15)
    assert lo["cs2_over_c2"]>0 and hi["cs2_over_c2"]>0
    with pytest.raises(ValueError,match="density gap"):
        eos.enthalpy_at_density(4.5)

def test_smooth_limit_reuses_validated_background_and_mode_equations():
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(contract(False)))
    cfg=dict(rtol=3e-10,atol=1e-12,max_step_h=.002)
    a=solve_sharp_background(eos,4.6,**cfg)
    b=solve_background(eos,4.6,**cfg)
    assert a.mass_msun==pytest.approx(b.mass_msun,rel=2e-10)
    assert a.radius_km==pytest.approx(b.radius_km,rel=2e-10)
    settings=CowlingSettings(rtol=1e-7,atol=1e-10)
    for frequency in [100.,300.,900.]:
        assert matching_residual(a,frequency,settings,conversion="slow")==pytest.approx(
            smooth_residual(b,frequency,settings),abs=2e-8)

def test_tov_mass_radius_agree_with_independent_pressure_coordinate():
    bar=WholeCrustBarotrope.from_contract(contract())
    eos=SharpContractEnthalpyEOS(bar)
    pstop=1e-8
    hstop=bar.state_at_pressure(pstop)["enthalpy"]
    bg=solve_sharp_background(eos,5.6,r0_km=.001,rtol=2e-11,atol=1e-13,
                              max_step_h=.001,surface_h=hstop)
    ec=eos.state_at_enthalpy(bg.h_c)["epsilon"]*MEV_FM3_TO_GEOM_KM2
    y=[.001,4*np.pi*ec*.001**3/3]
    pstart=eos.state_at_enthalpy(bg.h_start)["p"]
    for start,end,side in [(pstart,1.6,"high"),(1.6,pstop,"low")]:
        def rhs(p,y):
            radius,m=y
            state=bar.state_at_pressure(max(p,pstop),transition_side=side)
            ep=state["epsilon"]*MEV_FM3_TO_GEOM_KM2
            pp=p*MEV_FM3_TO_GEOM_KM2
            dr=-radius*(radius-2*m)/((ep+pp)*(m+4*np.pi*radius**3*pp))*MEV_FM3_TO_GEOM_KM2
            return [dr,4*np.pi*radius**2*ep*dr]
        sol=solve_ivp(rhs,(start,end),y,method="DOP853",rtol=2e-11,atol=1e-13,max_step=.01)
        assert sol.success
        y=sol.y[:,-1]
    assert bg.radius_km==pytest.approx(y[0],rel=2e-7)
    assert bg.mass_geom_km==pytest.approx(y[1],rel=2e-7)
    assert len(bg.interfaces)==1
    assert bg.interfaces[0]["r_jump_km"]==0 and bg.interfaces[0]["m_jump_km"]==0

def test_invalid_conversion_and_zero_frequency_rejected():
    s=dict(r_km=10.,Phi=-.2,lambda_metric=.3,dPhi_dr_km=.02,p=3.,epsilon=100.)
    with pytest.raises(ValueError):transfer_interface([1.,2.],s,s,.03,"unspecified")
    with pytest.raises(ValueError):transfer_interface([1.,2.],s,s,0.,"slow")

def test_center_taylor_start_must_not_swallow_a_sharp_interface():
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(contract()))
    with pytest.raises(ValueError,match="smaller r0"):
        solve_sharp_background(eos,5.+1e-8,r0_km=.001)
    resolved=solve_sharp_background(eos,5.+1e-8,r0_km=.0001)
    assert len(resolved.interfaces)==1
    assert resolved.interfaces[0]["r_km"]>.0001
