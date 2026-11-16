"""Independent derivatives, published saturation properties and composition contracts."""
import math
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

from compose_eos.skyrme_uniform import (
    PARAMETERS, HBARC, MN, MP, ME, MMU, nuclear_state, lepton_state, beta_equilibrium,
)

@pytest.mark.parametrize("model", PARAMETERS)
@pytest.mark.parametrize("nn,np_", [(0.08,0.08),(0.15,0.01),(0.4,0.05),(0.9,0.1)])
def test_analytic_chemical_potentials_against_energy_gradient(model,nn,np_):
    h=1e-6
    state=nuclear_state(model,nn,np_)
    dn=(nuclear_state(model,nn+h,np_)["energy_density_mev_fm3"]-
        nuclear_state(model,nn-h,np_)["energy_density_mev_fm3"])/(2*h)
    dp=(nuclear_state(model,nn,np_+h)["energy_density_mev_fm3"]-
        nuclear_state(model,nn,np_-h)["energy_density_mev_fm3"])/(2*h)
    assert state["mu_n_mev"]==pytest.approx(dn,rel=2e-9)
    assert state["mu_p_mev"]==pytest.approx(dp,rel=2e-9)

@pytest.mark.parametrize("model,nsat,energy", [
    ("RG_KDE0v1",.165,-16.23),("RG_SLY4",.160,-15.97)])
def test_published_symmetric_matter_saturation(model,nsat,energy):
    fit=minimize_scalar(lambda n:nuclear_state(model,n/2,n/2)["binding_energy_density_mev_fm3"]/n,
                        bounds=(.12,.19),method="bounded")
    assert fit.success
    # Published rounded coefficients/kinetic constant have finite precision.
    assert fit.x==pytest.approx(nsat,abs=0.001)
    assert fit.fun==pytest.approx(energy,abs=.03)

@pytest.mark.parametrize("mu,mass",[(.511,ME),(1.,ME),(100.,ME),(106.,MMU),(300.,MMU)])
def test_relativistic_lepton_energy_against_quadrature(mu,mass):
    k=math.sqrt((mu-mass)*(mu+mass))
    expected=quad(lambda p:p*p*math.sqrt(p*p+mass*mass),0,k,
                  epsabs=1e-12,epsrel=2e-13)[0]/(math.pi**2*HBARC**3)
    state=lepton_state(mu,mass)
    assert state["energy_density_mev_fm3"]==pytest.approx(expected,rel=3e-12,abs=1e-30)

@pytest.mark.parametrize("model", PARAMETERS)
@pytest.mark.parametrize("density",[.08,.16,.5,1.6])
def test_beta_equilibrium_charge_baryon_and_first_law(model,density):
    state=beta_equilibrium(model,density)
    assert state["Y_n"]+state["Y_p"]==pytest.approx(1,abs=2e-15)
    assert state["Y_p"]==pytest.approx(state["Y_e"]+state["Y_mu"],abs=2e-13)
    assert state["euler_residual_mev_fm3"]==pytest.approx(0,abs=2e-9)
    h=density*1e-5
    plus=beta_equilibrium(model,density+h)["energy_density_mev_fm3"]
    minus=beta_equilibrium(model,density-h)["energy_density_mev_fm3"]
    assert (plus-minus)/(2*h)==pytest.approx(state["mu_n_mev"],rel=2e-8)

@pytest.mark.parametrize("model",PARAMETERS)
def test_optional_muons_change_high_density_composition(model):
    npe=beta_equilibrium(model,.5,include_muons=False)
    npemu=beta_equilibrium(model,.5)
    assert npe["Y_mu"]==0
    assert npemu["Y_mu"]>0
    assert npemu["Y_p"]>npe["Y_p"]

@pytest.mark.parametrize("density",[0,-1,math.nan,math.inf])
def test_invalid_density_refused(density):
    with pytest.raises(ValueError):
        beta_equilibrium("RG_SLY4",density)

def test_kde_density_exponent_uses_published_value():
    assert PARAMETERS["RG_KDE0v1"].alpha==.1673
    assert PARAMETERS["RG_KDE0v1"].alpha != 1/6

@pytest.mark.parametrize("model",PARAMETERS)
def test_generated_table_retains_composition_and_acausal_rows(model):
    from compose_eos.skyrme_uniform import compose_table, sound_speed_ratio
    table, rows=compose_table(model,[.08,.16,.5,1.6])
    assert table.n==4
    assert np.all(table.Yn_present & table.Yp_present & table.Ye_present & table.Ymu_present)
    assert table.composition_available and table.beta_equilibrated
    assert np.allclose(table.Yp,table.Ye+table.Ymu,atol=2e-13,rtol=0)
    assert sound_speed_ratio(model,1.6)>1
    # Mathematical generation includes this row and reports domain separately.
    assert rows[-1]["n_b_fm3"]==1.6

@pytest.mark.parametrize("grid",[[.1],[.1,.1],[.2,.1],[.1,math.nan]])
def test_invalid_generation_grid_refused(grid):
    from compose_eos.skyrme_uniform import compose_table
    with pytest.raises(ValueError):
        compose_table("RG_SLY4",grid)

@pytest.mark.parametrize("model",PARAMETERS)
def test_sound_speed_derivative_refines(model):
    from compose_eos.skyrme_uniform import sound_speed_ratio
    fine=sound_speed_ratio(model,.5,relative_step=1e-4)
    coarse=sound_speed_ratio(model,.5,relative_step=2e-4)
    assert fine==pytest.approx(coarse,rel=1e-7)
