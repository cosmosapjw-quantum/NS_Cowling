"""Observed original-source evidence for stage-05, not invented model equivalence."""
import hashlib,json
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"data/cowling_rmf_reconciliation_20261119"
PREV=ROOT/"data/cowling_eos_source_calibration_20261119"
M03=ROOT/"data/cowling_maxwell_native_20261119"

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def test_full_original_archive_source_hashes_are_locked():
    rows=read(DATA/"EXACT_DENSITY_AUDIT.json")["models"]
    for name in ("FSU2","TM1e","DDME2"):
        base=M03/"provider_original_unified"/name/"raw/eos.thermo"
        assert sha(base)==rows[name]["source_identity"]["original_eos_thermo_sha256"]
        assert rows[name]["original_compose_archive_authentic_record"][
            "same_bytes_as_compose_download"] is True


def test_close_original_composition_does_not_conceal_thermo_mismatch():
    rows=read(DATA/"EXACT_DENSITY_AUDIT.json")["models"]
    for name in ("FSU2","TM1e"):
        q=rows[name]["original_beta_equilibrium_fraction_comparison"]
        assert q["max_abs_composition_fraction_difference"]<8e-6
        assert q["pressure_max_abs_rel"]>.03
        assert q["bulk_pressure_source_identical"] is False
        assert rows[name]["provider_thermodynamic_derivative"][
            "median_abs_rel_Pminus_n2dF"]>.005
        assert rows[name]["independent_local_thermodynamic_derivative"][
            "median_abs_rel_Pminus_n2dF"]<2e-7
        sat=rows[name]["symmetric_matter_saturation"]
        assert abs(sat["K_MeV"]-sat["published_saturation"]["K_MeV"])<.5
    dd=rows["DDME2"]["original_beta_equilibrium_fraction_comparison"]
    assert dd["pressure_max_abs_rel"]<1e-4
    assert rows["DDME2"]["provider_thermodynamic_derivative"][
         "median_abs_rel_Pminus_n2dF"]<1e-5


def test_original_charge_Q4_electron_relation_fails_without_model_fitting():
    rows=read(DATA/"charge_chemical_potentials/RESULT.json")["models"]
    for name in ("FSU2","TM1e"):
        q=rows[name]
        assert q["within_rounding_rows"]==0
        assert q["max_abs_mu_e_discrepancy_MeV"]>.9
        assert max(abs(z["source_lepton_mu_Q5_MeV"]) for z in q["rows"])==0
    dd=rows["DDME2"]
    assert dd["max_abs_mu_e_discrepancy_MeV"]<.01
    assert dd["within_rounding_rows"]>=60


def test_original_Yp_fixed_independent_RMF_pressure_still_disagrees():
    x=read(DATA/"fixed_composition/RESULT.json")["models"]
    for name in ("FSU2","TM1e"):
        c=x[name]
        assert c["max_composition_error_to_source_pressure_gap_ratio"]<.001
        assert c["min_source_pressure_abs_error_MeVfm3"]>.05
        assert c["max_proton_fraction_perturbation_effect_MeVfm3"]<.0003


def test_unfitted_literature_precision_variation_cannot_close_TM1e_gap():
    z=read(DATA/"EXACT_DENSITY_AUDIT.json")["models"]["TM1e"]
    control=z["parameter_table_precision_hypothesis"]
    assert control["no_parameter_fitting"]
    for baseline,refined in zip(control["baseline_pressure_error"],
                                control["high_precision_pressure_error"]):
        assert abs(refined)>0.99*abs(baseline)


@pytest.mark.parametrize("model,inner_len",[
    ("FSU2",52),("TM1e",86),("DDME2",71)])
def test_paired_stellar_crust_completely_identical_and_core_replace_only(
        model,inner_len):
    z=read(DATA/"paired_core/RESULT.json")["cases"][model]["css"]
    assert z["status"]=="PAIRED_THREE_LAYER_RUN"
    assert z["exactly_identical_original_inner_samples"]==inner_len
    assert z["crust_frozen_against_provider"] is True
    actual=z["star"]
    assert actual["status"]=="FIXED_1P4_NATIVE_TOV_F_P1_BOTH_CONVERSIONS_VERIFIED"
    assert abs(actual["native_TOV"]["M_Msun"]-1.4)<3e-6
    for conv in ("slow","rapid"):
        r=actual["native_cowling"][conv]
        assert r["status"]=="F_P1_COMPLETE"
        assert set(r["accepted_modes"])=={"f","p1"}
        assert all(x["max_interface_traction"]<1e-8
                   for x in r["accepted_modes"].values())


def test_different_homogeneous_core_changes_observables_in_model_specific_way():
    cases=read(DATA/"paired_core/RESULT.json")["cases"]
    old=read(PREV/"phase_local_fixed_mass/RESULT.json")["cases"]
    deltas={}
    for name in ("FSU2","TM1e","DDME2"):
        new=cases[name]["css"]["star"]["native_TOV"]["R_km"]
        original=old[name]["native_TOV"]["R_km"]
        deltas[name]=new-original
    assert .085<deltas["FSU2"]<.10
    assert .07<deltas["TM1e"]<.08
    assert abs(deltas["DDME2"])<.001
    assert cases["FSU2"]["maxwell"]["status"]=="PAIRED_THREE_LAYER_RUN"
    for model in ("TM1e","DDME2"):
        assert cases[model]["maxwell"]["status"]=="NATIVE_OR_STELLAR_REJECTED"


def test_machine_readable_summary_never_claims_original_fsu2_maxwell():
    x=read(DATA/"SOURCE_RECONCILIATION_OUTCOME.json")
    assert x["source_equivalent_homogeneous_core_for_FSU2_TM1e"] is False
    for name in ("FSU2","TM1e","DDME2"):
        z=x["models"][name]
        assert z["generated_RMF_plus_original_inner_crust_Maxwell_test"][
            "original_provider_phase6_Maxwell_not_verified"] is True
    assert x["source_immutable"] is True


def test_original_provider_eos_mr_not_independent_unified_stellar_reference():
    mr=read(DATA/"provider_mr/RESULT.json")
    assert mr["status"]=="FULL_UNIFIED_MR_BYTES_MATCH_CRUST_ONLY_ALL_THREE"
    for name,v in mr["models"].items():
        assert v["exact_same_eos_mr_bytes_as_crust_only"]
        assert v["original_thermo_fields_differ"]
        assert v["unified_source_number_of_density_knots"]>v["crust_only_number_of_density_knots"]
        assert v["unified_eos_mr_sha256"]==v["crust_only_eos_mr_sha256"]
        original=(M03/"provider_original_unified"/name/"raw/eos.mr")
        crust=(ROOT/"data/maxwell_three_families_20261117/intake/references"/
               ("GPPVA_"+name)/"eos.mr")
        assert sha(original)==sha(crust)
        assert sha(original)==v["unified_eos_mr_sha256"]
    overall=read(DATA/"SOURCE_RECONCILIATION_OUTCOME.json")
    assert overall[
        "provider_eos_mr_is_reused_crust_only_not_independent_unified_TOV_reference"]


def test_same_density_pressure_gap_is_algebraically_partitioned_without_fitting():
    a=read(DATA/"EXACT_DENSITY_AUDIT.json")["models"]
    for name in ("FSU2","TM1e","DDME2"):
        q=a[name]["pressure_difference_decomposition"]
        assert q["max_abs_balance_error_MeV_fm3"]<1e-8
        assert q["rms_model_source_pressure_difference_MeV_fm3"]>0
        assert q["rms_energy_term_over_observed"]>0
        assert q["rms_provider_identity_over_observed"]>0
        assert q["source_eos_unchanged"]
        if name in ("FSU2","TM1e"):
            assert q["rms_provider_identity_over_observed"]>.25
            assert q["rms_model_identity_over_observed"]<.001


def test_stage05_summary_keeps_both_pressure_mismatch_mechanisms():
    report=read(DATA/"SOURCE_RECONCILIATION_OUTCOME.json")["models"]
    for name in ("FSU2","TM1e","DDME2"):
        comp=report[name]["pressure_gap_energy_and_identity_decomposition"]
        assert comp["max_abs_balance_error_MeV_fm3"]<1e-8
        assert comp["rms_energy_term_over_observed"]>0
        assert comp["rms_provider_identity_over_observed"]>0
        assert comp["source_eos_unchanged"] is True
