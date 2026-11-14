"""E390-03: update the DAG node status and emit the consolidated native readback
receipt after the CLD composition reconstruction and CompOSE group-abundance
closure were delivered natively.

This script updates only bookkeeping fields (status, outputs, receipt pointers,
remaining_acceptance) in the DAG JSON; it does NOT change the contract, the
dependency graph, the gate, the legacy_method_semantics, or any other node.
"""
from __future__ import annotations
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DAG_PATH = ROOT / "data" / "eos390_native_fortran_unification_dag.v1.json"


def git_commit_sha() -> str:
    return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()


def main() -> int:
    dag = json.loads(DAG_PATH.read_text(encoding="utf-8"))
    updated = False
    for n in dag["nodes"]:
        if n["id"] != "E390-03":
            continue
        # status: CLD composition reconstruction + CompOSE closure delivered natively;
        # the only remaining subscope is the full Thomas-Fermi Wigner-Seitz gas
        # composition (currently the pure-neutron-gas approximation) and the
        # per-model absolute Abar/Zbar (blocked by the CLD continuous-phase model).
        n["status"] = (
            "PARTIAL_NATIVE_DDRMF_HOMOGENEOUS_AND_CLD_COMPOSITION_PASS_"
            "COMPOSE_CLOSURE_PASS_TF_GAS_ABSOLUTE_ABAR_ZBAR_REMAIN"
        )
        # outputs: mark the three delivered subscopes, leave the TF-gas one as remaining
        n["outputs"] = [
            "native DD-RMF homogeneous symmetry energy, beta-equilibrium and "
            "electron/muon branch diagnostics (DELIVERED, b3d80cc)",
            "native CLD cell composition reconstruction for GPPVA_DD2/DDME2/TW over "
            "the inner-crust density grid (DELIVERED: 3 models x 20 points = 60 rows, "
            "5 physics invariants PASS; Scurto 2024 CLD on the Fortran DD-RMF solver; "
            "Lattimer surface tension, RPW Coulomb lattice, Nelder-Mead 3-var)",
            "row-level CompOSE-standard group abundances with baryon/charge closure "
            "(DELIVERED: C3 Yn+Yp+Yh=1 PASS, C4 Yp+(Zbar/Abar)Yh=Ye+Ymu PASS at 1e-9; "
            "Zbar/Abar=rho_Ip/rho_I exact; Abar/Zbar absolute = NaN, CLD continuous-"
            "phase model does not fix cluster nucleon count)",
            "per-model reconstructed composition and readback receipts (DELIVERED: "
            "data/receipts/eos390_m3/cld/ and data/receipts/eos390_m3/compose_closure/; "
            "full CompOSE eos.thermo/eos.compo table assembly is E390-07 scope)",
            "full Thomas-Fermi Wigner-Seitz gas composition beyond the pure-neutron "
            "approximation (REMAINING: requires model-specific TF functionals, "
            "separate from the CLD continuous-phase cluster model)",
        ]
        n["partial_scope"] = {
            "delivered": [
                "native DD-RMF homogeneous solver (b3d80cc)",
                "native CLD cell reconstruction (this commit)",
                "CompOSE group-abundance closure C3/C4 (this commit)",
            ],
            "remaining": [
                "full Thomas-Fermi Wigner-Seitz gas composition",
                "absolute Abar/Zbar (blocked by CLD continuous-phase model)",
            ],
        }
        n["receipt"] = {
            "ddrmf_homogeneous": "data/receipts/eos390_m3/E390-03_DDRMF_NATIVE_SOLVER_20261114.json",
            "cld_composition": "data/receipts/eos390_m3/cld/E390-03_CLD_NATIVE_SOLVER_receipt.json",
            "compose_closure": "data/receipts/eos390_m3/compose_closure/E390-03_COMPOSE_CLOSURE_receipt.json",
        }
        n["remaining_acceptance"] = [
            "GitHub CI build/test pass for the new compose_eos_cld module",
            "full Thomas-Fermi Wigner-Seitz gas composition (separate subscope)",
        ]
        updated = True
        break

    if not updated:
        raise SystemExit("E390-03 node not found in DAG")

    # bump the provenance parent_commit reference in the DAG metadata if present
    if "metadata" in dag and "parent_commit" in dag["metadata"]:
        dag["metadata"]["parent_commit"] = git_commit_sha()

    DAG_PATH.write_text(json.dumps(dag, indent=2, sort_keys=False, ensure_ascii=False) + "\n",
                        encoding="utf-8")
    print(f"Updated E390-03 node status in {DAG_PATH}")
    print(f"  status: {dag['nodes'][5]['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
