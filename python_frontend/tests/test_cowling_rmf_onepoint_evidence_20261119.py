"""Exact locked provider one-density and independent TM1e author evidence."""
import hashlib,json
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"data/cowling_original_rmf_onepoint_20261119"
PROVIDER=ROOT/"data/cowling_maxwell_native_20261119/provider_original_unified"

def load(p):
    return json.loads(Path(p).read_text())
def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def test_exact_original_gppva_onepoint_source_authority_and_phase():
    prior=load(ROOT/"data/cowling_rmf_reconciliation_20261119/EXACT_DENSITY_AUDIT.json")["models"]
    one=load(DATA/"RESULT.json")["models"]
    assert set(one)=={"FSU2","TM1e","DDME2"}
    for name,rec in one.items():
        p=PROVIDER/name/"raw/eos.thermo"
        expected=prior[name]["source_identity"]["original_eos_thermo_sha256"]
        assert sha(p)==expected==rec["source_authority"]["original_GPPVA_source_thermo_sha256"]
        assert rec["selected_density_is_original_provider_exact_knot"]
        assert rec["no_source_interpolation"] is True
        assert .19 < rec["source"]["original_source_n_B_fm3"] < .21
        assert not rec["source_authority"]["original_GPPVA_eos_micro_available"]
        assert not rec["source_authority"]["original_GPPVA_software_generator_obtained"]


def test_non_equivalent_Pais_GPPVA_bulk_is_not_hidden_behind_correct_local_fields():
    one=load(DATA/"RESULT.json")["models"]
    for name in ("FSU2","TM1e"):
        z=one[name]
        delta=z["source_minus_local"]
        assert abs(delta["provider_minus_local_pressure_MeV_fm3"])>.1
        assert abs(delta["provider_minus_local_pressure_MeV_fm3"])<.15
        assert abs(z["source"]["source_muEuler_minus_muQ3_MeV"])>.5
        assert abs(z["source"]["original_electron_mu_Q_minus_Fermi_MeV"])>1.
        assert abs(z["source_minus_local_fixed_provider_Yp_pressure_MeV_fm3"]-
                   delta["provider_minus_local_pressure_MeV_fm3"])<.001
        for key in ("local_mean_field_components_at_local_beta_equilibrium",
                    "local_mean_field_components_at_original_printed_Yp_Ye_Ymu"):
            b=z[key]
            assert abs(b["epsilon_recomposition_residual_MeV_fm3"])<1e-8
            assert abs(b["pressure_recomposition_residual_MeV_fm3"])<1e-8
            assert abs(b["epsilon_explicit_source_potential_residual_MeV_fm3"])<1e-8
            assert abs(b["pressure_explicit_source_potential_residual_MeV_fm3"])<1e-8
            assert max(abs(v) for v in b["field_equation_residuals_fm3"].values())<1e-8
            assert b["source_type"]=="LOCAL_RMF_FIELDS_NOT_ORIGINAL_PROVIDER_FIELDS"
            energy=b["epsilon_components_MeV_fm3"]
            pressure=b["pressure_components_MeV_fm3"]
            assert energy["energy_omega_quartic"]==pytest.approx(
                   3*pressure["pressure_omega_quartic"])
            assert energy["energy_omega_rho_mixed"]==pytest.approx(
                   3*pressure["pressure_omega_rho_mixed"])
            assert pressure["pressure_sigma_mass"]==pytest.approx(
                   -energy["energy_sigma_mass"])
    dd=one["DDME2"]
    assert abs(dd["source_minus_local"]["provider_minus_local_pressure_MeV_fm3"])<.001
    assert abs(dd["source"]["source_muEuler_minus_muQ3_MeV"])<1e-4


def test_original_shen_author_exact_row_reproduced_without_leptons():
    doc=load(DATA/"RESULT.json")
    compare=doc["independent_author_TM1e"]
    author=compare["original_author_source"]
    assert compare["status"]=="INDEPENDENT_SHEN_TM1e_AUTHOR_BULK_EQUATION_VERIFIED"
    assert compare["accepted_by_preregistered_tolerance"]
    assert compare["same_GPPVA_generator"] is False
    assert author["leptons_not_included_in_original_author_EOS4"]
    assert author["Mstar_MeV"]==523.8616
    assert author["pressure_baryons_MeV_fm3"]==4.639332
    delta=compare["local_minus_shen_table"]
    assert compare["derived_sigma_status"].startswith("DERIVED_FROM_PUBLISHED_MSTAR")
    assert abs(compare["author_implied_sigma_MeV"]+41.2944989)<1e-5
    assert abs(compare["local_sigma_in_author_sign_MeV"]-compare["author_implied_sigma_MeV"])<.00025
    assert abs(delta["Mstar_MeV"])<.002
    assert abs(delta["pressure_baryons_MeV_fm3"])<.0001
    assert abs(delta["E_per_baryon_relative938_MeV"])<.005
    assert abs(delta["mu_n_relative938_MeV"])<.01
    assert abs(delta["mu_p_relative938_MeV"])<.02


def test_bounded_FSU2_reference_supports_only_saturation_benchmark():
    m=load(DATA/"RESULT.json")["literature_FSU2"]
    assert "NOT_PAIS_GPPVA_ONEPOINT" in m["source_type"]
    assert "THREE_DIGIT" in m["original_field_status"]
    assert abs(m["local_scalar_phi_MeV"]-m["original_author_implied_scalar_phi_from_rounded_mstar_MeV"])<.006
    assert abs(m["local_FSU2_SNM_mstar_over_mN"]-.593)<.001
    assert abs(m["local_FSU2_PNM_pressure_MeV_fm3"]-5.81)<.02


def test_no_original_provider_fields_fabricated_inside_source_row():
    allcases=load(DATA/"RESULT.json")["models"]
    for name,z in allcases.items():
        assert "original_phi" not in z["source"]
        assert "original_omega0" not in z["source"]
        assert "original_rho" not in z["source"]
        assert "original_component_pressure" not in z["source"]
        assert z["physical_gate"]=="NORMALIZED_P_EPS_MU_Y_ONLY_NO_PAIS_MEAN_FIELDS"
