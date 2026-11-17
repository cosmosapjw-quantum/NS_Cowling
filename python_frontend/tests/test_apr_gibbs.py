import numpy as np
import pytest
from compose_eos.apr_gibbs import (gibbs_coexistence,gibbs_equilibrium_state,
    gibbs_state_at_volume_fraction)
from compose_eos.apr_uniform import beta_state,equilibrium_state,maxwell_coexistence
from compose_eos.cold_leptons import leptons


@pytest.fixture(scope="module",params=[False,True])
def coexistence(request):
    return request.param,gibbs_coexistence(request.param)


def test_gibbs_two_charge_endpoints_match_neutral_beta_branches(coexistence):
    muons,phase=coexistence
    assert phase["left"]["n_b"]<phase["right"]["n_b"]
    assert phase["multistart_max_density_or_fraction_difference"]<1e-10
    assert phase["sampled_tangent_plane"]["minimum_excess_mev_fm3"]>=0
    for endpoint,branch in ((phase["left"],"low"),(phase["right"],"high")):
        pure=beta_state(endpoint["n_b"],branch,muons)
        for key in ("epsilon","pressure","mu_n","mu_p","mu_e","Y_p","Y_e","Y_mu"):
            assert endpoint[key]==pytest.approx(pure[key],rel=2e-10,abs=2e-10)
        assert max(abs(v) for v in endpoint["gibbs_residuals"].values())<1e-9


def test_gibbs_density_charge_species_and_thermodynamics(coexistence):
    muons,phase=coexistence
    path=phase["volume_fraction_path"]
    assert np.all(np.diff([row["n_b"] for row in path])>0)
    assert np.all(np.diff([row["pressure"] for row in path])>0)
    for row in path:
        chi=row["high_phase_volume_fraction"]
        n=(1-chi)*row["n_low"]+chi*row["n_high"]
        yp=((1-chi)*row["n_low"]*row["Yp_low"]+chi*row["n_high"]*row["Yp_high"])/n
        assert row["n_b"]==pytest.approx(n)
        assert row["Y_p"]==pytest.approx(yp)
        assert row["Y_n"]+row["Y_p"]==pytest.approx(1)
        assert row["Y_p"]==pytest.approx(row["Y_e"]+row["Y_mu"],abs=1e-12)
        assert abs(row["euler_residual_mev_fm3"])<1e-9
        assert max(abs(v) for v in row["gibbs_residuals"].values())<1e-9
        assert row["mu_n"]-row["mu_p"]==pytest.approx(row["mu_e"],abs=1e-10)


def test_gibbs_local_phases_have_opposite_charge_and_global_leptons(coexistence):
    muons,phase=coexistence
    row=gibbs_state_at_volume_fraction(.5,muons)
    lep=leptons(row["mu_e"],muons)
    nl=lep["n_e"]+lep["n_mu"]
    charge_low=row["n_low"]*row["Yp_low"]-nl
    charge_high=row["n_high"]*row["Yp_high"]-nl
    assert charge_low*charge_high<0
    assert .5*(charge_low+charge_high)==pytest.approx(0,abs=1e-12)
    low,high=row["nuclear_phase_states"].values()
    assert row["epsilon"]==pytest.approx(.5*(low["epsilon"]+high["epsilon"])+lep["epsilon_leptons"])
    assert row["pressure"]==pytest.approx(.5*(low["pressure"]+high["pressure"])+lep["pressure_leptons"])


def test_gibbs_energy_envelope_derivative_and_pressure(coexistence):
    muons,phase=coexistence
    for fraction in (.2,.5,.8):
        n=phase["left"]["n_b"]+fraction*(phase["right"]["n_b"]-phase["left"]["n_b"])
        state=gibbs_equilibrium_state(n,phase,muons)
        assert state["n_b"]==pytest.approx(n,abs=1e-12)
        for h in (1e-5,3e-6):
            left=gibbs_equilibrium_state(n-h,phase,muons)
            right=gibbs_equilibrium_state(n+h,phase,muons)
            derivative=(right["epsilon"]-left["epsilon"])/(2*h)
            assert derivative==pytest.approx(state["mu_n"],rel=2e-9)
            assert n*derivative-state["epsilon"]==pytest.approx(state["pressure"],abs=3e-7)
            assert (right["pressure"]-left["pressure"])/(right["epsilon"]-left["epsilon"])>0


def test_gibbs_energy_does_not_exceed_locally_neutral_maxwell(coexistence):
    muons,phase=coexistence
    maxwell=maxwell_coexistence(muons)
    for n in np.linspace(phase["left"]["n_b"],phase["right"]["n_b"],15):
        state=gibbs_equilibrium_state(float(n),phase,muons)
        alternative=equilibrium_state(float(n),maxwell,muons)
        assert state["epsilon"]<=alternative["epsilon"]+1e-9
    with pytest.raises(ValueError,match="species"):
        gibbs_equilibrium_state(.22,phase,not muons)


@pytest.mark.parametrize("chi",[-.1,1.1,float("nan"),float("inf")])
def test_gibbs_refuses_invalid_phase_fraction(chi):
    with pytest.raises(ValueError,match="volume fraction"):
        gibbs_state_at_volume_fraction(chi)


@pytest.mark.parametrize("density",[0.,2.1,float("nan"),float("inf")])
def test_gibbs_refuses_invalid_density(density):
    with pytest.raises(ValueError,match="density"):
        gibbs_equilibrium_state(density)
