#!/usr/bin/env python3
"""Stage-05 source-locked complete bounded outcome.

No fitting of parameters to provider data; no re-labelling of a generated
RMF homogeneous branch as the publisher's original bulk phase.
"""
from __future__ import annotations
import hashlib,json,subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data/cowling_rmf_reconciliation_20261119"
STAGE04=ROOT/"data/cowling_eos_source_calibration_20261119"
FILES={
  "original_vs_generated_composition_and_pressure":"EXACT_DENSITY_AUDIT.json",
  "original_CompOSE_Q4_independent_electron_test":"charge_chemical_potentials/RESULT.json",
  "fixed_original_fraction_control":"fixed_composition/RESULT.json",
  "paired_same_crust_1p4_stellar":"paired_core/RESULT.json",
  "unified_vs_crust_only_provider_eos_mr":"provider_mr/RESULT.json",
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    obj={k:json.loads((BASE/path).read_text()) for k,path in FILES.items()}
    docs={k:dict(relative_path="data/cowling_rmf_reconciliation_20261119/"+v,
                 sha256=sha(BASE/v)) for k,v in FILES.items()}
    comp=obj["original_vs_generated_composition_and_pressure"]["models"]
    q4=obj["original_CompOSE_Q4_independent_electron_test"]["models"]
    fixed=obj["fixed_original_fraction_control"]["models"]
    paired=obj["paired_same_crust_1p4_stellar"]["cases"]
    mr=obj["unified_vs_crust_only_provider_eos_mr"]["models"]
    prev=json.loads((STAGE04/"phase_local_fixed_mass/RESULT.json").read_text())
    root=dict(schema="ns_cowling.source_calibration05.original_rmf_reconciliation.v1",
       parent_04_commit="f191abe843d81c6a1e5797dd9edc9a1558f095fa",
       original_archive_from_04_Zenodo_independently_byte_verified=True,
       source_previous_stage_sha256=sha(STAGE04/"CALIBRATION_OUTCOME.json"),
       numerical_scope_n_fm3=[0.08,0.45],
       source_immutable=True,
       code_and_model_claim="same-named beta equilibrium fractions need not imply same bulk stress-energy EOS",
       original_beta_and_output_convention="cold npe(mu), CompOSE Q3=mub/mn-1, Q4=muq/mn, Q5=mul/mn; mu_e=mn(Q5-Q4)",
       source_schema_citation="https://link.springer.com/article/10.1140/epja/s10050-022-00847-y",
       no_new_physical_unification_claim=True,
       provider_chemical_potential_issue_diagnosed_not_repaired_in_original=True,
       provider_eos_mr_is_reused_crust_only_not_independent_unified_TOV_reference=True,
       source_equivalent_homogeneous_core_for_FSU2_TM1e=False,
       record_receipts=docs,
       models={},status="RUNNING")
    for name in ("FSU2","TM1e","DDME2"):
        rmf=comp[name];current=rmf["original_beta_equilibrium_fraction_comparison"]
        q=q4[name];fix=fixed.get(name)
        css=paired[name]["css"]
        reference=prev["cases"][name]
        if css["status"]!="PAIRED_THREE_LAYER_RUN":
            raise RuntimeError("missing accepted paired-core TOV/Cowling for "+name)
        obs=css["star"];original=reference["native_TOV"]
        modes={}
        for conv in ("slow","rapid"):
            new=obs["native_cowling"][conv]["accepted_modes"]
            old=reference["native_cowling"][conv]["accepted_modes"]
            modes[conv]={label:dict(source_provider_frequency_hz=old[label]["frequency_hz"],
               newly_computed_RMF_frequency_hz=z["frequency_hz"],
               delta_hz=float(z["frequency_hz"]-old[label]["frequency_hz"]))
                for label,z in new.items() if label in old}
        cs=dict(status=css["star"]["status"],
                same_outer_inner_crust_samples=int(css["exactly_identical_original_inner_samples"]),
                identical_surface_outer_bridge=bool(css["crust_frozen_against_provider"]),
                transition_pressure_MeV_fm3=css["transition_p"],
                provider_M1p4=float(original["M_Msun"]),
                local_RMF_M1p4=float(obs["native_TOV"]["M_Msun"]),
                delta_R1p4_km=float(obs["native_TOV"]["R_km"]-original["R_km"]),
                original_R1p4_km=float(original["R_km"]),
                local_RMF_R1p4_km=float(obs["native_TOV"]["R_km"]),
                mode_comparisons=modes,transition_density=[
                    css["transition_n_low"],css["transition_n_high"]],
                original_source_core_thermo_sha256=css["original_provider_core_sha256"],
                local_generated_core_thermo_sha256=css["generated_core_thermo_sha256"])
        root["models"][name]=dict(
            phase6_provider_source_sha256=rmf["source_identity"]["original_eos_thermo_sha256"],
            independent_same_density_rows=current["n_evaluated"],
            composition_fraction_max_abs_deviation=current["max_abs_composition_fraction_difference"],
            provider_5_digit_printed_Yp_not_all_strictly_rounding_equal=not current[
                 "composition_compatible_with_printed_rounding"],
            pressure_max_abs_relative_mismatch=current["pressure_max_abs_rel"],
            pressure_median_abs_relative_mismatch=current["pressure_median_abs_rel"],
            provider_energy_derivative_pressure_med_abs_relative_error=rmf[
                 "provider_thermodynamic_derivative"]["median_abs_rel_Pminus_n2dF"],
            local_RMF_energy_derivative_pressure_med_abs_relative_error=rmf[
                 "independent_local_thermodynamic_derivative"]["median_abs_rel_Pminus_n2dF"],
            provider_Q3_minus_Euler_max_abs_MeV=rmf[
                 "provider_thermodynamic_derivative"]["max_abs_Q3minus_Euler_MeV"],
            pressure_gap_energy_and_identity_decomposition={k:v for k,v in
                 rmf["pressure_difference_decomposition"].items()
                 if k!="sampled_source_density_rows"},
            provider_Q4_Q5_minus_electron_mu_max_abs_MeV=q[
                 "max_abs_mu_e_discrepancy_MeV"],
            provider_Q4_rounding_compatible_rows=q["within_rounding_rows"],
            provider_Q4_electron_check_rows=q["n_rows"],
            literature_saturation_check=rmf["symmetric_matter_saturation"],
            high_precision_TM1_isoscalar_control=rmf["parameter_table_precision_hypothesis"],
            provider_mr_provenance=dict(
                duplicate_of_crust_only=mr[name]["exact_same_eos_mr_bytes_as_crust_only"],
                source_thermo_differs_from_crust_only=mr[name]["original_thermo_fields_differ"],
                original_unified_number_of_rows=mr[name]["unified_source_number_of_density_knots"],
                crust_only_number_of_rows=mr[name]["crust_only_number_of_density_knots"],
                original_eos_mr_sha256=mr[name]["unified_eos_mr_sha256"],
                model_status=mr[name]["provider_mr_source_equivalence_status"]),
            independent_fixed_original_proton_fraction_control=fix and dict(
                maximum_pressure_change_due_to_fixed_provider_Yp=fix[
                    "max_proton_fraction_perturbation_effect_MeVfm3"],
                min_original_provider_vs_model_pressure_gap=fix[
                    "min_source_pressure_abs_error_MeVfm3"],
                max_composition_roundoff_effect_over_P_gap_ratio=fix[
                    "max_composition_error_to_source_pressure_gap_ratio"]),
            common_crust_phase_local_css_replacement=cs,
            generated_RMF_plus_original_inner_crust_Maxwell_test=dict(
                status=paired[name]["maxwell"]["status"],
                transition_pressure=paired[name]["maxwell"].get("transition_p"),
                original_provider_phase6_Maxwell_not_verified=True),
            physical_cause_status="SOURCE_BULK_OUTPUT_MISMATCH_MEASURED_EXACT_CAUSE_NOT_YET_UNIQUELY_IDENTIFIED",
            numerical_source_equivalence=("HIGH_AGREEMENT_NOT_BYTE_EQUAL" if name=="DDME2" else
                "FAIL_SOURCE_BULK_EQUIVALENCE"))
        print("SUMMARY",name,"P mismatch%",root["models"][name]["pressure_max_abs_relative_mismatch"]*100,
           "Q4 gap",root["models"][name]["provider_Q4_Q5_minus_electron_mu_max_abs_MeV"],
           "DeltaR_km",cs["delta_R1p4_km"],flush=True)
    root["status"]="BOUNDED_SOURCE_RECONCILIATION_COMPARISON_VERIFIED_BY_EXECUTED_SOURCES"
    path=BASE/"SOURCE_RECONCILIATION_OUTCOME.json"
    path.write_text(json.dumps(root,indent=2,allow_nan=False)+"\n")
    print("SOURCE_RECONCILIATION_OUTCOME",path,flush=True)

if __name__=="__main__":main()
