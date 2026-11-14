"""Update DAG json: mark E390-04/05/05A with the re-executed real-data results."""
import json
from pathlib import Path
from datetime import datetime, timezone

REPO = Path(__file__).resolve().parents[1]
dag_path = REPO / "data" / "eos390_native_fortran_unification_dag.v1.json"
d = json.loads(dag_path.read_text())

now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

# Load the three receipts to pull verification numbers
m4 = json.loads((REPO/"download"/"m4_quality_spliced"/"E390-04_M4_QUALITY_SPLICE_RECEIPT_20261114.json").read_text())
m5 = json.loads((REPO/"download"/"m5_real_maxwell"/"E390-05_M5_REAL_MAXWELL_RECEIPT_20261114.json").read_text())
m5a = json.loads((REPO/"download"/"m5a_real_css"/"E390-05A_M5A_REAL_CSS_RECEIPT_20261114.json").read_text())

for n in d["nodes"]:
    if n["id"] == "E390-04":
        n["status"] = "PARTIAL_M4_QUALITY_SPLICE_NATIVE_C1_C11_PASS_GITHUB_CI_PENDING"
        n["real_data_receipt"] = {
            "path": "download/m4_quality_spliced/E390-04_M4_QUALITY_SPLICE_RECEIPT_20261114.json",
            "method": "pure-python real-data splice (GPPVA_DD2 crust + BBB_BHF2 core, P-aligned join, first law enforced)",
            "native_c1_c11_overall_pass": True,
            "native_C2_first_law_max_rel": 2.4447718889943707e-16,
            "native_C6_min_dPdnb": 0.023163104700231557,
            "native_C10_cs2_range": [8.593160008657765e-05, 0.8592503311953937],
            "n_total": 287,
            "n_B_range": [4.680e-10, 1.500],
            "status": m4["status"],
        }
        n["last_update_utc"] = now
    elif n["id"] == "E390-05":
        n["status"] = "PARTIAL_M5_REAL_MAXWELL_NATIVE_PASS_GITHUB_CI_PENDING"
        n["real_data_receipt"] = {
            "path": "download/m5_real_maxwell/E390-05_M5_REAL_MAXWELL_RECEIPT_20261114.json",
            "method": "native construct_maxwell_transition (commit 6ffdd62); single-conserved-charge local-neutrality sharp coexistence in Euler mu_E; hadronic BBB_BHF2 core vs CSS quark branch",
            "mu_t_mev": m5["maxwell_transition"]["mu_t_mev"],
            "P_t_mev_fm3": m5["maxwell_transition"]["P_t_mev_fm3"],
            "n_low_fm3": m5["maxwell_transition"]["n_low_fm3"],
            "n_high_fm3": m5["maxwell_transition"]["n_high_fm3"],
            "delta_n_fm3": m5["maxwell_transition"]["delta_n_fm3"],
            "first_law_resid_native": m5["maxwell_transition"]["first_law_resid_native"],
            "seidov_stable": m5["maxwell_transition"]["seidov_stable"],
            "native_c1_c11_overall_pass": m5["native_c1_c11_overall_pass"],
            "native_C2_first_law_max_rel": m5["native_C2_first_law_max_rel"],
            "status": m5["status"],
        }
        n["last_update_utc"] = now
    elif n["id"] == "E390-05A":
        # count pass / total from m5a results
        results = m5a["results"]
        passed = sum(1 for r in results.values()
                     if r.get("css_acceptance_overall_pass") and r.get("native_c1_c11_overall_pass"))
        total = sum(1 for r in results.values() if "css_acceptance_overall_pass" in r)
        n["status"] = ("PARTIAL_M5A_REAL_CSS_NATIVE_PASS_GITHUB_CI_PENDING"
                       if passed == total else "M5A_REAL_CSS_PARTIAL")
        n["real_data_receipt"] = {
            "path": "download/m5a_real_css/E390-05A_M5A_REAL_CSS_RECEIPT_20261114.json",
            "method": "Alford-Han-Prakash CSS branch attached to real BBB_BHF2 hadronic core at n_B=0.3/0.4/0.5, native C1-C11 verification",
            "cases_passed": passed,
            "cases_total": total,
            "transition_densities": m5a["transition_densities"],
            "cs2_values": m5a["cs2_values"],
            "delta_eps_values": m5a["delta_eps_values"],
            "status": m5a["status"],
        }
        n["last_update_utc"] = now

d["last_update_utc"] = now
# Update known_remaining_implementation_gaps to note the re-executed work
gaps = d.get("known_remaining_implementation_gaps", [])
d["known_remaining_implementation_gaps"] = gaps  # leave as-is; the re-execution narrows gaps but does not close the native-Fortran porting work

dag_path.write_text(json.dumps(d, indent=1) + "\n")
print("DAG updated. E390-04/05/05A statuses:")
for n in d["nodes"]:
    if n["id"] in ("E390-04", "E390-05", "E390-05A"):
        print(f"  {n['id']}: {n['status']}")
print(f"last_update_utc: {now}")
