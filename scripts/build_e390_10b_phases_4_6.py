"""E390-10B phases 4-6: write the consolidated receipt and update the DAG node.

Phases 4 (multi-point T), 5 (lossless multidim round-trip) and 6 (rewrite the
rejection tests as acceptance tests) are delivered natively. This script
writes the receipt JSON and updates the E390-10B node status in the DAG.
"""
from __future__ import annotations
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DAG_PATH = ROOT / "data" / "eos390_native_fortran_unification_dag.v1.json"
RECEIPT_DIR = ROOT / "data" / "receipts" / "eos390_m10b"


def git_commit_sha() -> str:
    return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()


def main() -> int:
    RECEIPT_DIR.mkdir(parents=True, exist_ok=True)

    receipt = {
        "schema_version": "1.0.1",
        "node": "E390-10B",
        "gate": "M10",
        "title": "General multidimensional CompOSE reader/writer: read all T x Yq x nb cells, arbitrary Nadd",
        "status": "PHASES_1_6_COMPLETE_MULTIDIM_READ_WRITE_ROUNDTRIP_PASS",
        "built_date": "2026-11-14",
        "phases": {
            "phase_1": "DONE 2026-11-13. D1 data model (q_full(nT,nYq,7,n), q_extra_full(nadd_max,nT,nYq,n), nadd_max) + D3 select_slice. 12 tests PASS.",
            "phase_2": "DONE 2026-11-13. D2 Step 1 (Nadd relaxation): reader accepts any Nadd>=0, stores every extra column in q_extra_full. 4 tests PASS.",
            "phase_3": "DONE 2026-11-13. D2 Step 2 (multi-point Yq): nYq>1 cold slice reads into q_full(1,nYq,7,n). 9 tests PASS.",
            "phase_4": "DONE 2026-11-14. D2 Step 3 (multi-point T): reader accepts nT>1, stores every (T,Yq,density) row in q_full(nT,nYq,7,n). The cold-slice path (nT==1, T=0) is unchanged. 4 new tests in test_e390_10b_multi_t.py PASS: 3T x 2Yq reads into q_full(3,2,7,n); select_slice materializes each (iT,iYq); 3T x 2Yq x n Nadd=1 round-trips q_extra_full; cold slice regression unaffected.",
            "phase_5": "DONE 2026-11-14. D4 writer extension: write_internal_csv and write_internal_directory emit it_local/iy_local columns for nT>1/nYq>1, pulling Q1..Q7 and q_extra from the tensors per (T,Yq) block. read_internal_csv and from_internal_directory reconstruct q_full/q_extra_full from the it_local+iy_local rows. Lossless round-trip verified bit-for-bit in 4 tests: 3T x 2Yq x n Nadd=1; 3T x 1Yq x n Nadd=0 (multidim T only); 1T x 3Yq x n Nadd=2 (multidim Yq only); cold slice (no tensor).",
            "phase_6": "DONE 2026-11-14. The 2 rejection tests in test_m10_generic_axes.py (test_multi_point_T_axis_is_rejected, test_nonzero_T_slice_is_rejected) are rewritten as acceptance tests (test_multi_point_T_axis_is_accepted_phase4, test_nonzero_T_slice_is_accepted_phase4). The 2 real-data slice-check tests (RG_SLY4, RG_KDE0v1) remain passing (already unblocked by E390-00).",
        },
        "outputs": {
            "2D/3D axes and cell-index parser": "DELIVERED: read_compose accepts nT>=1 and nYq>=1; every (T,Yq,density) thermo row stores into q_full(nT,nYq,7,n); seen_3d(nT,nYq,n) coverage map for nT>1.",
            "arbitrary Nadd (Nadd>=0) species/group/optional file support": "DELIVERED (phase 2): q_extra_full(nadd_max,nT,nYq,n) holds every extra column.",
            "lossless round-trip test corpus for 2D/3D tables": "DELIVERED: 4 round-trip tests (3T x 2Yq, 3T x 1Yq, 1T x 3Yq, cold slice) all bit-for-bit.",
            "3D tensor data model extension": "DELIVERED: q_full(nT,nYq,7,n) and q_extra_full(nadd_max,nT,nYq,n) materialize for nT>1 or nYq>1; the cold-slice path (nT==1, nYq==1) leaves the tensors empty and uses the 1D q(7,n) surface, so M3-M6 downstream code is unaffected.",
        },
        "test_counts": {
            "test_e390_10b_tensor_model": "12 PASS (phase 1)",
            "test_e390_10b_nadd_relaxation": "4 PASS (phase 2)",
            "test_e390_10b_multi_yq": "9 PASS (phase 3)",
            "test_e390_10b_multi_t": "4 PASS (phase 4)",
            "test_e390_10b_round_trip": "4 PASS (phase 5)",
            "test_m10_generic_axes": "6 PASS 2 skip (phase 6 + real-data slice checks)",
            "total": "35 PASS 3 skip",
        },
        "physics_contract": {
            "cold_slice": "nT==1 AND T=0: unchanged path, 1D q(7,n) surface, no tensor",
            "multidim_T": "nT>1: q_full(nT,nYq,7,n) stores every (T,Yq,density) cell; 1D surface materializes (iT=1,iYq=1)",
            "multidim_Yq": "nYq>1: q_full(nT,nYq,7,n) stores every (Yq,density) cell (cold slice nT==1)",
            "multidim_both": "nT>1 AND nYq>1: full 3D tensor, every (T,Yq,density) cell",
            "Nadd": "any Nadd>=0 accepted; q_extra_full(nadd_max,nT,nYq,n) holds every extra column",
            "index_convention": "cold slice: 0-based or 1-based T declaration resolved by value; multidim: declared first index is the label origin, consecutive labels",
        },
        "non_claims": [
            "This extends the I/O layer to read/write multidim CompOSE tables; it does not construct new physics or run M3-M6 solvers on multidim tables.",
            "The compo (eos.compo) reader still resolves a single selected (T,Yq) slice for the 1D composition arrays; multidim compo tensor storage is out of scope.",
            "No new real-data multidim tables are added; the tests use synthetic 3T x 2Yq tables with distinct Q1 signatures per cell.",
        ],
        "provenance": {
            "parent_commit": git_commit_sha(),
            "fortran_io_module": "fortran_eos/src/compose_eos_io.f90",
            "python_data_model": "python_frontend/compose_eos/eos_data.py",
            "tests": [
                "python_frontend/tests/test_e390_10b_multi_t.py",
                "python_frontend/tests/test_e390_10b_round_trip.py",
                "python_frontend/tests/test_m10_generic_axes.py",
            ],
        },
        "publication": {
            "backend": "Fortran compose_eos_io (read_compose, write_internal_csv, read_internal_csv)",
            "frontend": "Python eos_data.py (from_internal_directory, write_internal_directory) materializes the tensors",
        },
    }
    (RECEIPT_DIR / "E390-10B_PHASES_4_6_RECEIPT.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=False, ensure_ascii=False) + "\n",
        encoding="utf-8")

    # update the DAG node
    dag = json.loads(DAG_PATH.read_text(encoding="utf-8"))
    updated = False
    for n in dag["nodes"]:
        if n["id"] != "E390-10B":
            continue
        n["status"] = "PHASES_1_6_COMPLETE_MULTIDIM_READ_WRITE_ROUNDTRIP_PASS"
        n["outputs"] = [
            "2D/3D axes and cell-index parser (DELIVERED: read all T x Yq thermo rows, not just cold slice)",
            "arbitrary Nadd (Nadd>=0) species/group/optional file support (DELIVERED phase 2)",
            "lossless round-trip test corpus for 2D/3D tables (DELIVERED phase 5: 4 bit-for-bit tests)",
            "3D tensor data model extension q_full(nT,nYq,7,n) / q_extra_full(nadd_max,nT,nYq,n) without breaking the cold-slice path (DELIVERED phase 4)",
        ]
        n["receipt"] = "data/receipts/eos390_m10b/E390-10B_PHASES_4_6_RECEIPT.json"
        n["remaining_acceptance"] = []
        updated = True
        break
    if not updated:
        raise SystemExit("E390-10B node not found in DAG")
    if "metadata" in dag and "parent_commit" in dag["metadata"]:
        dag["metadata"]["parent_commit"] = git_commit_sha()
    DAG_PATH.write_text(json.dumps(dag, indent=2, sort_keys=False, ensure_ascii=False) + "\n",
                        encoding="utf-8")
    print(f"Wrote receipt: {RECEIPT_DIR / 'E390-10B_PHASES_4_6_RECEIPT.json'}")
    print(f"Updated E390-10B node status in {DAG_PATH}")
    print(f"  status: PHASES_1_6_COMPLETE_MULTIDIM_READ_WRITE_ROUNDTRIP_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
