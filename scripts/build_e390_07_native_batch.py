"""E390-07 (M7): batch orchestration native pipeline integration receipt.

The M7 batch orchestrator stubs (run_first_law_constrained_hermite_pressure_bridge_task,
run_css_task, run_smooth_task) are replaced with a native Fortran pipeline:
  read crust -> internal; read core -> internal; match-core METHOD n_bridge
  drop_composition=true; write internal -> CompOSE; verify readback.

This script records the integration status. The native match-core step
succeeds (combined EOS produced), but the readback verification step
(write -> read) currently returns 0 rows for drop_composition=true tables,
which is the remaining M7 acceptance gap.
"""
from __future__ import annotations
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DAG_PATH = ROOT / "data" / "eos390_native_fortran_unification_dag.v1.json"
RECEIPT_DIR = ROOT / "data" / "receipts" / "eos390_m7"


def main() -> int:
    RECEIPT_DIR.mkdir(parents=True, exist_ok=True)
    sha = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    receipt = {
        "schema_version": "1.0.1",
        "node": "E390-07",
        "gate": "M7",
        "title": "Implement native manifest-driven 390 batch orchestration",
        "status": "PARTIAL_NATIVE_PIPELINE_MATCH_CORE_PASS_READBACK_REMAINS",
        "built_date": "2026-11-14",
        "outputs": {
            "pair_manifest": "DELIVERED: m7_manifest_generator.py produces the 390-task manifest (13 crust x 5 core x 6 methods) with task_id, provenance, output_folder, status slots",
            "fortran_batch_command": "PARTIAL: the native pipeline (read -> match-core -> write -> verify) is wired into m7_batch_orchestrator.py via _native_pipeline; match-core succeeds (combined EOS, e.g. 92 rows) but the write->verify readback returns 0 rows for drop_composition=true tables",
            "per_pair_receipts": "DELIVERED: the orchestrator emits per-task receipts with PASS/FAIL/BLOCKED status, hashes, and native evidence; the remaining gap is the readback row count for composition-dropped tables",
        },
        "delivered": [
            "m7_manifest_generator.py: 390-task manifest (13 crust x 5 core x 6 smooth methods)",
            "m7_batch_orchestrator.py: native Fortran pipeline wired into all 4 method stubs (FLCHPB, css, cubic_F/P, quintic_F/P)",
            "per-task fail-closed receipts with native evidence (execution_backend, combined_n, match_report)",
        ],
        "remaining": [
            "readback verification: write -> read returns 0 rows when drop_composition=true (the compo file is omitted, so the verify step needs to accept composition-absent tables)",
            "full 390-task batch run once the readback gap is closed",
            "maxwell method remains a separate M5 workflow (not in the M7 smooth matrix)",
        ],
        "physics_contract": {
            "pipeline": "read crust CompOSE -> internal; read core CompOSE -> internal; match-core METHOD n_bridge drop_composition; write internal -> CompOSE; verify readback",
            "methods": "first_law_constrained_hermite_pressure_bridge, css, cubic_F, cubic_P, quintic_F, quintic_P (all native match-core smooth methods)",
            "drop_composition": "true (the M7 smooth matrix is thermodynamic-only; the native composition bridge is the E390-03 TF-gas subscope)",
            "maxwell_excluded": "maxwell is the sharp first-order transition (E390-05), not a smooth EOS390 v2 method",
        },
        "non_claims": [
            "The readback gap is a verify-step acceptance for composition-dropped tables, not a physics gap; match-core produces a complete combined thermodynamic EOS.",
            "No Python numerical fallback; the pipeline is native Fortran end-to-end.",
        ],
        "provenance": {
            "parent_commit": sha,
            "orchestrator": "scripts/m7/m7_batch_orchestrator.py",
            "manifest_generator": "scripts/m7/m7_manifest_generator.py",
        },
    }
    (RECEIPT_DIR / "E390-07_NATIVE_BATCH_RECEIPT.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=False, ensure_ascii=False) + "\n",
        encoding="utf-8")

    dag = json.loads(DAG_PATH.read_text(encoding="utf-8"))
    for n in dag["nodes"]:
        if n["id"] != "E390-07":
            continue
        n["status"] = "PARTIAL_NATIVE_PIPELINE_MATCH_CORE_PASS_READBACK_REMAINS"
        n["outputs"] = [
            "pair manifest (DELIVERED: 390-task manifest, 13 crust x 5 core x 6 methods)",
            "Fortran batch command (PARTIAL: native pipeline wired; match-core succeeds, readback gap for composition-dropped tables remains)",
            "per-pair receipts (DELIVERED: fail-closed per-task receipts with native evidence)",
        ]
        n["receipt"] = "data/receipts/eos390_m7/E390-07_NATIVE_BATCH_RECEIPT.json"
        n["remaining_acceptance"] = [
            "readback verification for composition-dropped (drop_composition=true) tables",
            "full 390-task batch run",
        ]
        break
    if "metadata" in dag and "parent_commit" in dag["metadata"]:
        dag["metadata"]["parent_commit"] = sha
    DAG_PATH.write_text(json.dumps(dag, indent=2, sort_keys=False, ensure_ascii=False) + "\n",
                        encoding="utf-8")
    print(f"Wrote E390-07 receipt and updated DAG: PARTIAL_NATIVE_PIPELINE_MATCH_CORE_PASS_READBACK_REMAINS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
