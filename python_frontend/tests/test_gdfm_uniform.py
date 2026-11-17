"""Independent author NMP and thermodynamic derivative checks of GDFM."""
import pytest

from compose_eos.gdfm_uniform import (
    MASS, MODEL_IDS, baryons, beta_state, coupling, parameters, symmetry_energy,
)


@pytest.mark.parametrize("model", MODEL_IDS)
def test_author_archived_nuclear_matter_parameters(model):
    p = parameters()[model]
    n = p["rho0_fm3"]
    h = n*1e-3
    state = baryons(model,n,.5)
    incompressibility = 9*(baryons(model,n+h,.5)["pressure_baryons"]
                          -baryons(model,n-h,.5)["pressure_baryons"])/(2*h)
    slope = 3*n*(symmetry_energy(model,n+h)-symmetry_energy(model,n-h))/(2*h)
    target = p["NMP"]
    assert state["epsilon_baryons"]/n-MASS == pytest.approx(target[1], abs=1e-4)
    assert incompressibility == pytest.approx(target[2], abs=1e-3)
    assert symmetry_energy(model,n) == pytest.approx(target[5], abs=1e-4)
    assert slope == pytest.approx(target[6], abs=1e-3)
    # The author file stores Gamma_rho = g_rho/2, resolved independently by J.
    assert abs(symmetry_energy(model,n,rho_vertex_scale=.5)-target[5]) > 5


@pytest.mark.parametrize("model", MODEL_IDS)
@pytest.mark.parametrize("density", [.03,.16,.8])
def test_chemical_potentials_are_energy_partial_derivatives(model,density):
    nn, np_ = .8*density,.2*density
    step = density*1e-5

    def energy(n,p):
        return baryons(model,n+p,p/(n+p))["epsilon_baryons"]
    state = baryons(model,density,.2)
    mun = (energy(nn+step,np_)-energy(nn-step,np_))/(2*step)
    mup = (energy(nn,np_+step)-energy(nn,np_-step))/(2*step)
    assert state["mu_n"] == pytest.approx(mun, rel=3e-8, abs=2e-5)
    assert state["mu_p"] == pytest.approx(mup, rel=3e-8, abs=2e-5)


@pytest.mark.parametrize("model", MODEL_IDS)
def test_beta_composition_and_euler_identity(model):
    for density in [.001,.03,.16,.8]:
        state = beta_state(model,density)
        assert state["Yp"] == pytest.approx(state["Ye"]+state["Ymu"],abs=2e-10)
        assert state["Yn"]+state["Yp"] == pytest.approx(1.0,abs=1e-15)
        assert abs(state["euler_residual"]) < 2e-8
        assert abs(state["scalar_residual_fm3"]) < 1e-11
    assert beta_state(model,.001)["Ymu"] == 0
    assert beta_state(model,.16)["Ymu"] > 0


def test_exact_coupling_derivative():
    model = parameters()["SPG_M1"]
    n, h = .22, 1e-6
    for coeff in model["couplings"].values():
        _, derivative = coupling(coeff,n,model["rho0_fm3"])
        numerical = (coupling(coeff,n+h,model["rho0_fm3"])[0]
                     -coupling(coeff,n-h,model["rho0_fm3"])[0])/(2*h)
        assert derivative == pytest.approx(numerical,rel=2e-9)
