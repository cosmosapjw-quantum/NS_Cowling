"""Record the 2026-11-15/16 continuation evidence in both DAGs.

Each touched node receives a new ``status``; the previous status is appended
to ``status_history`` and the evidence/remaining blockers are stored under
``update_20261116``.  Untouched nodes are left as they were.
Idempotent: re-running does not duplicate history entries.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECEIPTS = "data/receipts/core_crust_audit_20261115/"
DOC = "docs/CORE_CRUST_REAUDIT_CONTINUATION_20261115.md"

CLOUD = {
    "CLOUD-00": ("COMPLETE_ISOLATED_WORKTREE_DROPBOX_RECONCILED", {
        "evidence": [RECEIPTS + "DROPBOX_RECONCILIATION.json",
                     "work performed in an isolated worktree at the exact remote HEAD"],
        "finding": "Dropbox checkout holds only older committed versions or empty files; no unique work"}),
    "EOS-03": ("PARTIAL_SPG_MASS_AND_CHEMICAL_POTENTIAL_CORRECTED", {
        "evidence": ["SPG nucleon mass 938.9 MeV (Goegelein 2008; CompOSE 317/319 headers)",
                     "SPG mu_{n,p} now equal d eps/d n_i (test_rmf_chemical_potential_consistency)"],
        "remaining": ["SPG_M2/M4 parameter sets are placeholders of M1 (Scurto 2024 curves only)",
                      "RG cores: provider Q5 is a density-independent mu_l offset (1.2996 / 0.56698 MeV); "
                      "convention unresolved, beta contract not overridden",
                      "APR core acausal at high density (max c_s^2 1.4709)"]}),
    "EOS-04": ("PARTIAL_INNER_CRUST_V5_ADMITTED_OUTER_CRUST_SPLICE_MISSING", {
        "evidence": [RECEIPTS + "CLD_COLUMN_DEFECT_EVIDENCE.json",
                     "download/reconstruction_v5/regeneration_index.json",
                     RECEIPTS + "CLD_GAS_EQUILIBRIUM_v5.json"],
        "finding": ("former CLD writer put n_B in the Q1 slot (P = n_B^2); corrected writer derives "
                    "P and mu_B from eps(n_B), trims optimizer-box rows, refuses non-convex, "
                    "acausal and gas-nonequilibrium tables"),
        "remaining": ["outer-crust splice (CLD tables start at n_B >= 1e-3 fm^-3)",
                      "SPG_M3 rejected (interior box-bound minima)",
                      "composition columns not emitted"]}),
    "EOS-06": ("CSS_JOIN_RUN_NATIVE_8_ADMITTED_SURVIVORS", {
        "evidence": [RECEIPTS + "CENSUS_a05d0c9.json", "scripts/match_report_admission.py"],
        "remaining": ["per-task source-backed CSS parameter manifest (cs2 is solved per pair, not a source parameter)"]}),
    "COW-02": ("PHASE_SPLIT_ADAPTER_IMPLEMENTED_PROFILE_EXPORT_SCHEMA_PENDING", {
        "evidence": ["python_frontend/compose_cowling/css_adapter.py",
                     "python_frontend/tests/test_css_cowling_adapter.py"],
        "remaining": ["W/V/r/q/P/eps/h/chi/nu/lambda/nu' per-side profile export with schema and SHA receipt"]}),
    "COW-07": ("INTEGRATED_ON_CLD_COMPOSITE_EOS_SOURCE_ADMISSION_BLOCKED", {
        "evidence": [RECEIPTS + "CSS_SURVIVOR_STARS_a05d0c9.json"],
        "finding": ("phase-segmented TOV and slow/rapid Cowling roots computed for every admitted "
                    "survivor; radii are truncation radii (no outer crust); roots unlabelled"),
        "remaining": ["source-admitted microscopic crust with outer crust and composition",
                      "mode identity (COW-06)"]}),
}

E390 = {
    "E390-03": ("RESOLVED_20261114_NATIVE_DDRMF_CLD_TF_GAS_COMPOSITION_PASS_COMPOSE_CLOSURE_PASS"
                "__20261116_PYTHON_CLD_COMPOSE_WRITER_DEFECT_FIXED_V5", {
        "evidence": [RECEIPTS + "CLD_COLUMN_DEFECT_EVIDENCE.json", DOC]}),
    "E390-05": ("REMEDIATION_20261116_MAXWELL_MASS_GATE_ADJUDICATED_CONSTITUENT_MASS_DIFFERENCE", {
        "evidence": [DOC + " (section 'B re-adjudicated on v3 crusts')", RECEIPTS + "CENSUS_a05d0c9.json"],
        "finding": ("all 25 evaluated Maxwell attempts stop at the 1e-12 mass gate: DD-ME2/TW crusts use "
                    "939/939, SPG 938.9/938.9, DD2 939.56536/938.27203 versus provider cores "
                    "939.56533/938.272 -- Hamiltonian-level constituent-mass differences, not normalization"),
        "remaining": ["explicit labelled decision on mass-constant equivalence (DD2 vs cores, 3e-5 MeV)"]}),
    "E390-05A": ("CI_PASS", {"evidence": ["GitHub Actions runs 607-609 on PR #15 (fortran-eos-ci)"]}),
    "E390-06": ("NATIVE_TOV_LOVE_PASS_PHASE_SEGMENTED_COWLING_ADAPTER_ON_GENERATED_EOS", {
        "evidence": ["python_frontend/compose_cowling/css_adapter.py",
                     RECEIPTS + "CSS_SURVIVOR_STARS_a05d0c9.json"],
        "remaining": ["Cowling end-to-end tests still skip without raw M0 inputs"]}),
}


def apply(path: Path, updates: dict, head: str) -> list[str]:
    dag = json.loads(path.read_text(encoding="utf-8"))
    nodes = dag["nodes"] if isinstance(dag["nodes"], list) else list(dag["nodes"].values())
    changed = []
    for node in nodes:
        if node.get("id") not in updates:
            continue
        status, extra = updates[node["id"]]
        if node.get("status") != status:
            node.setdefault("status_history", []).append(node.get("status"))
            node["status"] = status
        node["update_20261116"] = {"git_head_before_update": head, **extra}
        changed.append(node["id"])
    path.write_text(json.dumps(dag, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return changed


def main() -> int:
    # This writer describes the earlier continuation, not the later re-audit.
    # Preserve the current authoritative state instead of replaying old labels.
    for filename in ("ns_cowling_cloud_agent_dag.v1.json",
                     "eos390_native_fortran_unification_dag.v1.json"):
        current = json.loads((ROOT / "data" / filename).read_text(encoding="utf-8"))
        if "adversarial_reaudit_20261116" in current:
            print("Historical DAG updater superseded by the 2026-11-16 re-audit; no changes.")
            return 0
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()
    print(apply(ROOT / "data" / "ns_cowling_cloud_agent_dag.v1.json", CLOUD, head))
    print(apply(ROOT / "data" / "eos390_native_fortran_unification_dag.v1.json", E390, head))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
