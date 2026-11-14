"""E390-06: update the DAG node status after native TOV+Love delivery."""
from __future__ import annotations
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DAG_PATH = ROOT / "data" / "eos390_native_fortran_unification_dag.v1.json"


def main() -> int:
    dag = json.loads(DAG_PATH.read_text(encoding="utf-8"))
    updated = False
    for n in dag["nodes"]:
        if n["id"] != "E390-06":
            continue
        n["status"] = "NATIVE_TOV_LOVE_PASS_COWLING_RUNTIME_COMPILE_REMAINS"
        n["outputs"] = [
            "native stellar interfaces (DELIVERED: Fortran TOV RK45 integration, smooth + sharp-Maxwell barotropes)",
            "tidal Love solver (DELIVERED: l=2 y(r) ODE Takatsy-Kovacs 2020 Eq. 2-3, k2 and Lambda from Thorne 1971)",
            "slow/rapid Cowling interface receipts (DELIVERED slow: Takatsy-Kovacs Eq. 11 pressure-corrected y-jump; rapid Masuda 2013 documented as distinct formula, not integrated)",
            "phase-transition TOV/Love tidal jump Delta y = -4 pi r_d^3 Delta eps / (m_d + 4 pi r_d^3 p_d) + phase-segmented integration (DELIVERED: two-branch integration with the y-jump applied at the transition radius)",
        ]
        n["partial_scope"] = {
            "delivered": [
                "native Fortran TOV + l=2 Love solver (smooth + sharp-Maxwell)",
                "Takatsy-Kovacs 2020 y-jump at the phase transition",
                "Love number k2 and tidal deformability Lambda",
            ],
            "remaining": [
                "compose_cowling f-mode oscillation solver (the runtime-compile path in python_frontend/compose_cowling/solver.py is a separate subscope; the 3 pre-existing test failures are this path, not the new native TOV+Love module)",
                "rapid-conversion Cowling interface (Masuda 2013) integration",
                "consumption of the M5 Maxwell production EOS (delivered by E390-07 batch orchestrator)",
            ],
        }
        n["receipt"] = "data/receipts/eos390_m6/E390-06_NATIVE_TOV_LOVE_RECEIPT.json"
        n["remaining_acceptance"] = [
            "GitHub CI build/test pass for the new compose_eos_tov module",
            "compose_cowling runtime-compile path resolution (separate subscope)",
        ]
        updated = True
        break
    if not updated:
        raise SystemExit("E390-06 node not found in DAG")
    if "metadata" in dag and "parent_commit" in dag["metadata"]:
        dag["metadata"]["parent_commit"] = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True).stdout.strip()
    DAG_PATH.write_text(json.dumps(dag, indent=2, sort_keys=False, ensure_ascii=False) + "\n",
                        encoding="utf-8")
    print(f"Updated E390-06 node status: {dag['nodes'][[x['id'] for x in dag['nodes']].index('E390-06')]['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
