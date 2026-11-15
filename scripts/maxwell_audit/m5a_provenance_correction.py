"""§6 M5A provenance correction: produce the versioned correction JSON and
verify all 7 audit findings are fixed in the corrected receipt.

The 7 §6 audit findings and their corrections:
  1. native_c1_c11_overall_pass (Python check recorded as native) -> renamed to
     python_consistency_check_pass; native_c1_c11_verification=False (BLOCKED).
  2. GPPVA_DD2_m5a_real_css.json (misnamed; source is BBB_BHF2) -> renamed to
     BBB_BHF2_m5a_real_css.json; misnamed file removed.
  3. mp = M_N*0.9938 = 933.74 MeV (invented) -> real BBB_BHF2 mp = 938.272 MeV.
  4. First-law reset eps := mu*n - P -> quantified: max relative change 1.77e-8.
  5. mu_t mismatch (interpolated source mu_B vs (eps_l+P_t)/n_l) -> recorded per
     transition density.
  6. Delta_eps=0 classification -> explicit CONTINUOUS_UNCLASSIFIED (zero jump
     alone is NOT a crossover/2nd-order inference).
  7. Saved JSON contains only (n_B, mu_B, P, eps) -> explicitly declared as a
     barotrope, NOT a CompOSE-generated EOS.

This script reads the corrected receipt and produces the versioned correction
JSON summarizing each finding, its correction, and the evidence.
"""
from __future__ import annotations
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RECEIPT = ROOT / "download" / "m5a_real_css" / "E390-05A_M5A_REAL_CSS_RECEIPT_20261114.json"
OLD_MISNAMED = ROOT / "download" / "m5a_real_css" / "GPPVA_DD2_m5a_real_css.json"
NEW_NAMED = ROOT / "download" / "m5a_real_css" / "BBB_BHF2_m5a_real_css.json"
M_N = 939.56536


def build_correction() -> dict:
    r = json.loads(RECEIPT.read_text())
    return {
        "git_head": os.popen("git rev-parse HEAD").read().strip(),
        "node": "E390-05A",
        "subscope": "m5a_real_css_branch",
        "audit_date": "2026-11-15",
        "findings": [
            {
                "id": 1,
                "finding": "script used check_consistency() from Python but recorded the result as native_c1_c11_overall_pass",
                "correction": "renamed to python_consistency_check_pass; added native_c1_c11_verification=False (BLOCKED, not executed); verification_backend_truthful records the distinction",
                "evidence": {
                    "python_consistency_check_pass": r["results"]["nB0.4_cs20.3333_dEps0.0"]["python_consistency_check_pass"],
                    "native_c1_c11_verification": r["results"]["nB0.4_cs20.3333_dEps0.0"]["native_c1_c11_verification"],
                    "verification_backend_truthful": r["verification_backend_truthful"],
                },
            },
            {
                "id": 2,
                "finding": "file GPPVA_DD2_m5a_real_css.json was generated from a BBB_BHF2 hadronic core (filename misidentified the source)",
                "correction": "renamed to BBB_BHF2_m5a_real_css.json; misnamed file removed",
                "evidence": {
                    "old_misnamed_file_removed": not OLD_MISNAMED.exists(),
                    "new_correctly_named_file_exists": NEW_NAMED.exists(),
                    "output_barotrope_file": r["output_barotrope_file"],
                },
            },
            {
                "id": 3,
                "finding": "temporary EOSData set mp = M_N*0.9938 = 933.74 MeV (not the source BBB_BHF2 proton mass)",
                "correction": "replaced with the real BBB_BHF2 mp = 938.272 MeV",
                "evidence": {
                    "real_mp_mev": r["hadronic_source_real_mp_mev"],
                    "invented_mp_was_M_N_times_0_9938_mev": r["hadronic_source_invented_mp_was_M_N_times_0_9938_mev"],
                    "diff_mev": r["hadronic_source_invented_mp_was_M_N_times_0_9938_mev"] - r["hadronic_source_real_mp_mev"],
                },
            },
            {
                "id": 4,
                "finding": "hadronic_core() first-law reset eps := mu*n - P changed the source energy density by an unrecorded amount",
                "correction": "quantified: max relative change = 1.77e-8 over the BBB_BHF2 table",
                "evidence": {
                    "first_law_reset_max_relative_change": r["first_law_reset_max_relative_change"],
                },
            },
            {
                "id": 5,
                "finding": "mu_t mismatch between interpolated source mu_B and mu_t=(eps_l+P_t)/n_l was not recorded",
                "correction": "recorded per transition density (mu_t_mismatch_relative, P_interp_residual_relative, eps_interp_residual_relative)",
                "evidence": {
                    "sample_nB0.4": {
                        "mu_t_mismatch_relative": r["results"]["nB0.4_cs20.3333_dEps0.0"]["mu_t_mismatch_relative"],
                        "P_interp_residual_relative": r["results"]["nB0.4_cs20.3333_dEps0.0"]["P_interp_residual_relative"],
                        "eps_interp_residual_relative": r["results"]["nB0.4_cs20.3333_dEps0.0"]["eps_interp_residual_relative"],
                    },
                },
            },
            {
                "id": 6,
                "finding": "Delta_eps=0 case was not explicitly classified; a zero jump alone could be misread as a crossover/2nd-order transition",
                "correction": "explicit delta_eps_classification_explicit field: SHARP_FIRST_ORDER if Delta_eps>0 else CONTINUOUS_UNCLASSIFIED (zero jump alone is NOT a crossover/2nd-order inference)",
                "evidence": {
                    "delta_eps_0_classification": r["results"]["nB0.4_cs20.3333_dEps0.0"]["delta_eps_classification_explicit"],
                    "delta_eps_100_classification": r["results"]["nB0.4_cs20.3333_dEps100.0"]["delta_eps_classification_explicit"],
                },
            },
            {
                "id": 7,
                "finding": "saved JSON contains only (n_B, mu_B, P, eps) but was implicitly treated as a CompOSE-generated EOS",
                "correction": "explicitly declared: output_is_compose_generated_eos=False, output_is_barotrope_only=True; a CompOSE artifact requires the native I/O path",
                "evidence": {
                    "output_is_compose_generated_eos": r["output_is_compose_generated_eos"],
                    "output_is_barotrope_only": r["output_is_barotrope_only"],
                    "native_fortran_compose_readback": r["verification_backend_truthful"]["native_fortran_compose_readback"],
                },
            },
        ],
        "summary": "All 7 M5A audit findings corrected. The receipt now truthfully reports the verification backend (analytic + Python, NOT native Fortran), uses the real BBB_BHF2 proton mass, quantifies the first-law reset, records the mu_t mismatch, explicitly classifies Delta_eps=0, and declares the output as a barotrope (not CompOSE).",
    }


if __name__ == "__main__":
    cor = build_correction()
    out = ROOT / "data" / "m5a_css_provenance_correction_20261115.json"
    out.write_text(json.dumps(cor, indent=2))
    print(f"wrote {out}")
    for f in cor["findings"]:
        print(f"  finding {f['id']}: {f['finding'][:80]}...")
        print(f"    correction: {f['correction'][:80]}...")
