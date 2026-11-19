"""Cold RMF one-density decomposition and independent author benchmarks."""
from pathlib import Path
import numpy as np
import pytest
from compose_eos.rmf_onepoint_components import (
    nonlinear_component_budget, read_shen_eos4_t00_exact_row,
)
from compose_eos.gppva_rmf import beta_state,baryons

ROOT=Path(__file__).resolve().parents[2]
EOS4=ROOT/"data/cowling_original_rmf_onepoint_20261119/external_reference/shen_EOS4_original_eos4.t00.gz"

@pytest.mark.parametrize("model,n",[
    ("FSU2",0.19160198),("TM1e",0.1904367),
    ("FSU2",0.1505),("TM1e",0.145)])
def test_local_component_pressure_energy_and_field_residuals(model,n):
    s=beta_state(model,n)
    b=nonlinear_component_budget(model,n,s["Y_p"],
          s["Y_e"]*n,s["Y_mu"]*n)
    assert b["source_type"]=="LOCAL_RMF_FIELDS_NOT_ORIGINAL_PROVIDER_FIELDS"
    assert abs(b["epsilon_recomposition_residual_MeV_fm3"])<1e-8
    assert abs(b["pressure_recomposition_residual_MeV_fm3"])<1e-8
    assert max(abs(x) for x in b["field_equation_residuals_fm3"].values())<1e-6
    assert abs(b["charge_residual_fm3"])<1e-10
    assert abs(b["mu_n_recomposition_residual_MeV"])<1e-9
    assert abs(b["mu_p_recomposition_residual_MeV"])<1e-9
    assert "energy_omega_rho_mixed" in b["epsilon_components_MeV_fm3"]
    assert "pressure_omega_rho_mixed" in b["pressure_components_MeV_fm3"]
    assert b["mstar_MeV"]>0


def test_explicit_fixed_source_fraction_is_not_automatically_beta_equilibrium():
    model="TM1e";n=.1904367;yp=.1
    b=nonlinear_component_budget(model,n,yp,n*yp,0.)
    assert b["charge_residual_fm3"]==pytest.approx(0,abs=1e-14)
    assert b["epsilon_recomposition_residual_MeV_fm3"]==pytest.approx(0,abs=1e-9)
    assert b["pressure_recomposition_residual_MeV_fm3"]==pytest.approx(0,abs=1e-9)
    assert abs(b["beta_residual_MeV"])>1.


def test_original_shen_tm1e_author_uniform_mstar_p_energy_mu():
    r=read_shen_eos4_t00_exact_row(EOS4,log10_rho_gcm3=14.5,Yp=.1)
    assert r["raw_numeric_source_rows"]==7150
    assert r["temperature_mev"]==0.
    assert r["n_B_fm3"]==pytest.approx(.1904367,abs=0)
    assert r["Mstar_MeV"]==pytest.approx(523.8616,abs=0)
    assert r["pressure_baryons_MeV_fm3"]==pytest.approx(4.639332,abs=0)
    assert r["free_energy_per_baryon_relative938_MeV"]==pytest.approx(8.674307,abs=0)
    assert r["free_baryon_neutron_fraction"]==pytest.approx(.9,abs=0)
    assert r["heavy_nucleus_fraction"]==0.
    b=nonlinear_component_budget("TM1e",r["n_B_fm3"],r["Yp"],0.,0.)
    assert abs(b["mstar_MeV"]-r["Mstar_MeV"])<.02
    assert abs(b["pressure_baryons_MeV_fm3"]-r["pressure_baryons_MeV_fm3"])<.001
    assert abs(b["epsilon_baryons_MeV_fm3"]/r["n_B_fm3"]-938.-r["free_energy_per_baryon_relative938_MeV"])<.01
    assert abs(b["mu_n_MeV"]-938.-r["mu_n_relative938_MeV"])<.05
    assert abs(b["mu_p_MeV"]-938.-r["mu_p_relative938_MeV"])<.05


def test_original_FSU2_literature_saturation_mstar_and_pure_neutron_P():
    n=.1505
    snm=baryons("FSU2",n,.5)
    pnm=baryons("FSU2",n,0.)
    assert abs(snm["effective_mass_n"]/939.-.593)<.002
    assert abs(pnm["pressure_baryons"]-5.81)<.03


def test_original_shen_reference_exact_grid_rejects_interpolation():
    with pytest.raises(ValueError,match="exact"):
        read_shen_eos4_t00_exact_row(EOS4,14.53,.1)
    with pytest.raises(ValueError,match="exact"):
        read_shen_eos4_t00_exact_row(EOS4,14.5,.103)
