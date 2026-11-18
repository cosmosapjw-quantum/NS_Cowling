"""Deterministic machine-readable delivery and independent evidence identities."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"data/all_origin_maxwell_css_20261118"
DOC=ROOT/"docs/all_origin_maxwell_css_20261118"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
c=read(D/"CAMPAIGN.json")
src=read(D/"SOURCE_CATALOG.json")
a=read(D/"INDEPENDENT_AUDIT.json")
s=read(D/"RESULTS_SUMMARY.json")
full=read(D/"FULL_SUITE_RECEIPT.json")
st=read(D/"whole_stars/CAMPAIGN.json")
root=read(D/"HERMITE_ROOT_CROSSCHECK.json")
assert c["status"]=="COMPLETED_WITH_EXPLICIT_ORIGIN_AND_SUPPORT_GATES" and len(c["cases"])==635
assert len(src["sources"])==95 and s["number_crust_sources"]==70
assert a["status"]=="PASS_SCOPED" and a["checked"]>=600 and not a["failed"]
assert full["actual_exit_code"]==0 and full["source_bytes_unchanged"]
assert st["both_conditional_whole"]==5 and root["eligible_dual_pair_count"]==7
assert s["finite_dual_success"]==7 and s["true_original_full_whole_success"]==0
assert s["matching_methods"]==["maxwell","css"]
assert s["input_interpolation_methods"]==["cubic","quintic"]
assert s["interpolation_matrix_2x2_total"]==28
assert s["actual_native_maxwell_interpolated_source_accepted"]==14
assert s["actual_native_css_interpolated_source_accepted"]==14
for name in ("ALL_CANDIDATE_PAIRS.csv","INDEPENDENT_CRUST_INVENTORY.csv","DUAL_FINITE_SUCCESS.csv",
             "CONDITIONAL_WHOLE_TOV.csv","MATCHING_2X2_VERIFIED.csv"):
    assert (D/name).exists()
assert len(list(DOC.glob("*.md")))>=2
files=sorted(set([*list(D.rglob("*.json")),*list(D.rglob("*.csv")),
                  *list(D.rglob("eos.*")),*list(DOC.glob("*.md")),
                  *list((ROOT/"scripts").glob("*20261118.py")),
                  ROOT/"python_frontend/tests/test_all_origin_maxwell_css_20261118.py"]))
files=[f for f in files if f.is_file() and f!=D/"DELIVERY_RECEIPT.json"]
receipt=dict(schema="ns_cowling.origin_isolated_maxwell_css_deliverable.v1",
  repository="cosmosapjw-quantum/NS_Cowling",
  source_head="bc09c87010fccb280c455d765e6239d3910fb0db",
  target_branch="research/all-origin-maxwell-css-20261118",
  scope="cold EOS original/corrected/reconstructed never cross-mixed; original outer+inner required for unconditional whole NS",
  advertised_outputs=["ALL_CANDIDATE_PAIRS.csv","INDEPENDENT_CRUST_INVENTORY.csv",
   "DUAL_FINITE_SUCCESS.csv","CONDITIONAL_WHOLE_TOV.csv","MATCHING_2X2_VERIFIED.csv","HERMITED_INPUT_2X2/NATIVE_MAXWELL/CAMPAIGN.json",
   "RESULTS_SUMMARY.json","HERMITE_ROOT_CROSSCHECK.json",
   "docs/all_origin_maxwell_css_20261118/RESULTS_KO.md",
   "docs/all_origin_maxwell_css_20261118/CORE_INDEPENDENT_CRUST_AND_FUNCTIONAL_CONTROLS_KO.md"],
  metrics=dict(source_entries=95,crust_entries=70,core_entries=25,
    same_origin_pairings=635,finite_dual_css_maxwell=7,
    native_interpolated_maxwell_accepted=14,
    native_interpolated_css_accepted=14,
    interpolation_x_matching_2_by_2=28,
    unconstrained_source_native_whole_verified=0,
    conditional_reconstructed_whole_TOV=5,independent_audit_checks=a["checked"],
    full_suite_exit_code=full["actual_exit_code"],
    full_suite_source_files=full["source_snapshot_count"]),
  runtime=dict(max_workers=c["max_workers"],seconds=c["elapsed_seconds"],
    average_candidates_per_second=len(c["cases"])/c["elapsed_seconds"],
    repeated_exact_input_hash_reused_legacy_Maxwell=195,
    optimization="cached EOS reads per source; vectorized bounded P(mu) root screening; 6 parallel workers; final SSOT file published atomically; independent native root gates"),
  provenance=dict(provider_original="source physical data kept byte identical",
    energy_primary_corrected="never paired to original",
    native_reconstructed="never paired to original or corrected",
    nuda_wheel_sha256=src["nuda_evidence"]["sha256"],
    raw_user_ZIP_never_republished=True),
  scientific_limitations=a["explicit_limitations"],
  exact_file_sha256={str(f.relative_to(ROOT)):sha(f) for f in files})
(D/"DELIVERY_RECEIPT.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
print("DELIVERY_PASS",len(files),"file identities",receipt["metrics"],"rate",receipt["runtime"]["average_candidates_per_second"])
