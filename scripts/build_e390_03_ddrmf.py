#!/usr/bin/env python3
"""Build a reproducible native DD-RMF beta-equilibrium/symmetry diagnostic grid.

This is deliberately a diagnostic grid, not a CompOSE EOS table and not a CLD
or Thomas-Fermi crust reconstruction. The Fortran executable owns all numerical
physics; Python writes JSON provenance and hashes only.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos.native import evaluate_ddrmf, executable  # noqa: E402

MODEL_IDS = ("GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW")
N_POINTS = 32
RHO_MIN = 1.0e-4
RHO_MAX = 1.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: build_e390_03_ddrmf.py OUTPUT_DIRECTORY")
    out = Path(sys.argv[1]).expanduser().resolve()
    if out.exists():
        raise FileExistsError(f"refusing to overwrite existing attempt directory: {out}")
    out.mkdir(parents=True)

    executable_path = executable()
    density_grid = [math.exp(math.log(RHO_MIN) + i * (math.log(RHO_MAX) - math.log(RHO_MIN)) / (N_POINTS - 1))
                    for i in range(N_POINTS)]
    registry = ROOT / "data" / "eos390_rmf_parameter_sets.v1.json"
    source_module = ROOT / "fortran_eos" / "src" / "compose_eos_rmf.f90"
    python_frontend = ROOT / "python_frontend" / "compose_eos" / "rmf_solver.py"

    model_receipts = []
    total_rows = 0
    for model_id in MODEL_IDS:
        rows = []
        for rho_b in density_grid:
            result = evaluate_ddrmf(model_id, rho_b)
            if result.get("status") != "ok":
                raise RuntimeError(f"native solver failed for {model_id} at {rho_b}: {result}")
            for field, tolerance in (("beta_residual_mev", 2e-8),
                                     ("charge_residual_fm3", 2e-12),
                                     ("pressure_identity_residual_mev_fm3", 2e-7),
                                     ("scalar_field_residual_mev", 2e-8)):
                if abs(float(result[field])) > tolerance:
                    raise RuntimeError(f"{model_id} {field}={result[field]} exceeds {tolerance}")
            rows.append(result)
        data_path = out / f"{model_id}_ddrmf_beta_symmetry_grid.v1.json"
        write_json(data_path, {
            "schema_version": "1.0.0",
            "model_id": model_id,
            "output_type": "native_homogeneous_diagnostic_grid_not_compose_eos",
            "density_support_fm3": [RHO_MIN, RHO_MAX],
            "density_point_count": len(rows),
            "symmetry_energy_definitions": {
                "second_derivative": "five-point central finite difference of E/A at alpha=0, h=0.05, S=0.5*d2(E/A)/d(alpha)^2",
                "finite_difference": "E/A(alpha=1)-E/A(alpha=0), with source-specific linear rest-mass term subtracted in each state",
            },
            "rows": rows,
        })
        model_receipt = {
            "model_id": model_id,
            "status": "NATIVE_DDRMF_BETA_SYMMETRY_GRID_PASS",
            "row_count": len(rows),
            "output_file": data_path.name,
            "output_sha256": sha256(data_path),
            "muon_inactive_rows": sum(r["muon_branch"] == "INACTIVE_BELOW_THRESHOLD" for r in rows),
            "muon_active_rows": sum(r["muon_branch"] == "ACTIVE" for r in rows),
            "max_abs_beta_residual_mev": max(abs(r["beta_residual_mev"]) for r in rows),
            "max_abs_charge_residual_fm3": max(abs(r["charge_residual_fm3"]) for r in rows),
            "max_abs_pressure_identity_residual_mev_fm3": max(abs(r["pressure_identity_residual_mev_fm3"]) for r in rows),
            "max_abs_scalar_field_residual_mev": max(abs(r["scalar_field_residual_mev"]) for r in rows),
        }
        model_receipts.append(model_receipt)
        total_rows += len(rows)

    receipt = {
        "schema_version": "1.0.0",
        "node": "E390-03",
        "subscope": "native_DD_RMF_homogeneous_symmetry_and_beta_equilibrium",
        "built_date": "2026-11-14",
        "status": "NATIVE_DDRMF_DIAGNOSTIC_GRID_PASS",
        "models": list(MODEL_IDS),
        "density_support_fm3": [RHO_MIN, RHO_MAX],
        "density_point_count_per_model": N_POINTS,
        "total_output_rows": total_rows,
        "physics_contract": {
            "matter": "cold homogeneous nucleonic matter, T=0",
            "weak_equilibrium": "mu_n-mu_p=mu_e; mu_e=mu_mu only when the muon branch is populated",
            "charge_neutrality": "n_p=n_e+n_mu",
            "leptons": "ideal zero-temperature relativistic Fermi gases with physical electron and muon masses",
            "density_dependence": "Typel-Wolter rational sigma/omega couplings and exponential rho coupling",
            "rearrangement": "included with the literature-positive sign in nucleon chemical potentials and pressure",
            "rho_isospin_convention": "registry Gamma_rho uses the original no-1/2 tau_3 convention; equivalent Scurto convention is g_rho=2*Gamma_rho with explicit 1/2",
            "symmetry_energy": "derivative and pure-neutron finite-difference definitions are both reported, not conflated",
        },
        "provenance": {
            "parameter_registry": "data/eos390_rmf_parameter_sets.v1.json",
            "parameter_registry_sha256": sha256(registry),
            "native_fortran_source": "fortran_eos/src/compose_eos_rmf.f90",
            "native_fortran_source_sha256": sha256(source_module),
            "python_frontend_source": "python_frontend/compose_eos/rmf_solver.py",
            "python_frontend_source_sha256": sha256(python_frontend),
            "native_executable_sha256": sha256(executable_path),
            "solver_backend": "native_fortran",
        },
        "model_receipts": model_receipts,
        "non_claims": [
            "This is not a reconstructed CompOSE EOS table.",
            "No crust construction, CLD minimization, Thomas-Fermi Wigner-Seitz calculation or cluster abundance reconstruction is performed here.",
            "The provider composition semantics for SPG/GPPVA Iphase=7 remain unresolved and are not used as solver targets.",
            "The GPPVA_TW crust provider conflict remains separate; TW is evaluated only as a homogeneous DD-RMF functional.",
            "No interpolation, extrapolation, smoothing or modification of source EOS tables is performed.",
        ],
    }
    receipt_path = out / "E390-03_DDRMF_NATIVE_SOLVER_RECEIPT.v1.json"
    write_json(receipt_path, receipt)
    manifest = {
        "schema_version": "1.0.0",
        "status": "ALL_DDRMF_DIAGNOSTIC_GRIDS_NATIVE_FORTRAN_VERIFIED",
        "model_count": len(MODEL_IDS),
        "density_point_count_per_model": N_POINTS,
        "total_output_rows": total_rows,
        "model_files": [
            {"path": m["output_file"], "sha256": m["output_sha256"], "rows": m["row_count"]}
            for m in model_receipts
        ],
        "receipt": {"path": receipt_path.name, "sha256": sha256(receipt_path)},
    }
    write_json(out / "ddrmf_grid_manifest.v1.json", manifest)
    print(json.dumps(manifest, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
