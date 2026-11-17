import pytest
from compose_eos.apr_uniform import nuclear_state, beta_state, M_N, M_P

@pytest.mark.parametrize("branch",["low","high"])
@pytest.mark.parametrize("nn,np_",[(.1,.03),(.3,.04),(.7,.15)])
def test_apr_chemical_potentials_are_energy_derivatives(branch,nn,np_):
    s=nuclear_state(nn,np_,branch)
    h=1e-6
    dn=(nuclear_state(nn+h,np_,branch)["epsilon"]-nuclear_state(nn-h,np_,branch)["epsilon"])/(2*h)
    dp=(nuclear_state(nn,np_+h,branch)["epsilon"]-nuclear_state(nn,np_-h,branch)["epsilon"])/(2*h)
    assert s["mu_n"]==pytest.approx(dn,rel=2e-9)
    assert s["mu_p"]==pytest.approx(dp,rel=2e-9)

def test_apr_symmetric_saturation_reference():
    s=nuclear_state(.08,.08)
    binding=s["epsilon"]/.16-(M_N+M_P)/2
    assert binding==pytest.approx(-16,abs=.15)
    assert abs(s["pressure"]) < .1

@pytest.mark.parametrize("branch",["low","high"])
@pytest.mark.parametrize("muons",[False,True])
def test_apr_beta_charge_and_thermodynamic_derivative(branch,muons):
    n=.4
    s=beta_state(n,branch,muons)
    assert s["Y_n"]+s["Y_p"]==pytest.approx(1)
    assert s["Y_p"]==pytest.approx(s["Y_e"]+s["Y_mu"],abs=2e-13)
    assert abs(s["euler_residual_mev_fm3"])<1e-8
    h=1e-5
    derivative=(beta_state(n+h,branch,muons)["epsilon"]-beta_state(n-h,branch,muons)["epsilon"])/(2*h)
    assert derivative==pytest.approx(s["mu_n"],rel=2e-8)


@pytest.mark.parametrize("muons",[False,True])
def test_maxwell_equal_pressure_chemical_potential_and_tangent(muons):
    from compose_eos.apr_uniform import maxwell_coexistence, equilibrium_state
    import numpy as np
    phase=maxwell_coexistence(muons)
    left,right=phase["left"],phase["right"]
    assert .19 < left["n_b"] < right["n_b"] < .26
    assert abs(phase["pressure_residual_mev_fm3"]) < 1e-9
    assert abs(phase["chemical_residual_mev"]) < 1e-9
    assert phase["multistart_max_density_difference_fm3"] < 1e-10
    for n in np.linspace(left["n_b"],right["n_b"],21)[1:-1]:
        mixed=equilibrium_state(float(n),phase,muons)
        assert mixed["branch"]=="mixed"
        assert mixed["mu_e"] is None and mixed["mu_p"] is None
        assert mixed["Y_n"]+mixed["Y_p"]==pytest.approx(1,abs=2e-14)
        assert mixed["Y_p"]==pytest.approx(mixed["Y_e"]+mixed["Y_mu"],abs=2e-13)
        assert abs(mixed["euler_residual_mev_fm3"]) < 1e-9
        for branch in ("low","high"):
            assert mixed["epsilon"] <= beta_state(float(n),branch,muons)["epsilon"]+1e-10
        h=1e-6
        derivative=(equilibrium_state(n+h,phase,muons)["epsilon"]
                    -equilibrium_state(n-h,phase,muons)["epsilon"])/(2*h)
        assert derivative==pytest.approx(phase["mu_b"],rel=1e-9)


def test_maxwell_composition_is_density_weighted_and_preserves_endpoints():
    from compose_eos.apr_uniform import maxwell_coexistence, equilibrium_state
    phase=maxwell_coexistence()
    left,right=phase["left"],phase["right"]
    n=(left["n_b"]+right["n_b"])/2
    mixed=equilibrium_state(n,phase)
    assert mixed["high_phase_volume_fraction"]==pytest.approx(.5)
    assert mixed["Y_p"]==pytest.approx(
        (left["n_b"]*left["Y_p"]+right["n_b"]*right["Y_p"])/(2*n))
    assert equilibrium_state(left["n_b"],phase)["branch"]=="low"
    assert equilibrium_state(right["n_b"],phase)["branch"]=="high"
    with pytest.raises(ValueError,match="species"):
        equilibrium_state(n,phase,False)


@pytest.mark.parametrize("density",[float("nan"),float("inf"),0.,-.1,2.1])
def test_maxwell_refuses_unsupported_density(density):
    from compose_eos.apr_uniform import equilibrium_state
    with pytest.raises(ValueError,match="density"):
        equilibrium_state(density)


@pytest.mark.parametrize("nn,np_",[(1e308,1e308),(1e-300,0.),(2.,.1)])
def test_nuclear_state_refuses_unsupported_total_density(nn,np_):
    with pytest.raises(ValueError,match="density"):
        nuclear_state(nn,np_)


def test_generator_phase_endpoints_remain_inside_requested_window():
    import runpy
    from pathlib import Path
    import numpy as np
    module=runpy.run_path(str(Path(__file__).resolve().parents[2]/"scripts/build_apr_uniform.py"))
    phase={"left":{"n_b":.20},"right":{"n_b":.24}}
    select=module["maxwell_native_density_grid"]
    grid=np.array([.1,.18,.21,.22,.25,.5])
    assert select(grid,phase)==pytest.approx([.1,.18,.20,.24,.25,.5])
    assert select(np.array([.5,.6,.7]),phase)==pytest.approx([.5,.6,.7])
    assert select(np.array([.21,.22,.25,.5]),phase)==pytest.approx([.24,.25,.5])
