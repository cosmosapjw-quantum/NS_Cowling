"""E390-03 CLD native Fortran reconstruction verification.

Drives the native compose_eos_fortran `cld-reconstruct` subcommand over the
inner-crust density grid for the three GPPVA DD-RMF models (DD2, DDME2, TW),
then verifies the CLD cell physics invariants the E390-03 contract requires:

  1. baryon conservation: rho_b == f*rho_I + (1-f)*rho_II
  2. cell charge neutrality: rho_e == f*rho_Ip + (1-f)*rho_IIp  (rho_IIp=0)
  3. beta-equilibrium at the cluster: |mu_n - mu_p - mu_e| small
  4. physical bounds: 0 < f < 1, R_d > 0, rho_I > rho_II >= 0, yp_cell in [0,1]
  5. thermodynamic monotonicity: E_binding_per_baryon monotonic (non-decreasing)
     across the inner-crust density ramp for each model

Python owns only path/manifest validation, native invocation, receipt
generation and presentation. The Fortran backend owns all numerical physics.
"""
from __future__ import annotations
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos.native import executable  # noqa: E402

MODEL_IDS = ("GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW")
N_POINTS = 20
RHO_MIN = 1.0e-3
RHO_MAX = 8.0e-2
SIGMA_0 = 1.0   # MeV/fm^2 Lattimer surface tension scale
B_EXP = 1.0     # Lattimer exponent

TOL_BARYON = 1.0e-9
TOL_CHARGE = 1.0e-9
TOL_BETA = 1.0e-3   # MeV; cluster beta-equilibrium residual tolerance


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")


def run_native_cld(model_id: str) -> list[dict]:
    exe = executable()
    cmd = [str(exe), "cld-reconstruct", model_id, str(SIGMA_0), str(B_EXP),
           str(RHO_MIN), str(RHO_MAX), str(N_POINTS)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"native cld-reconstruct failed for {model_id} (rc={result.returncode})\n"
            f"stderr:\n{result.stderr}")
    return json.loads(result.stdout)


def verify_point(p: dict, model_id: str) -> list[str]:
    """Return list of invariant violations for one CLD point."""
    issues = []
    rho_b = p["rho_b_fm3"]; f = p["f"]; r_d = p["R_d_fm"]
    rho_i = p["rho_I_fm3"]; rho_ip = p["rho_Ip_fm3"]
    rho_ii = p["rho_II_fm3"]; rho_iip = p["rho_IIp_fm3"]
    rho_e = p["rho_e_fm3"]
    mu_n = p["mu_n_MeV"]; mu_p = p["mu_p_MeV"]; mu_e = p["mu_e_MeV"]

    # 1. baryon conservation
    baryon_resid = abs(rho_b - (f*rho_i + (1.0 - f)*rho_ii))
    if baryon_resid > TOL_BARYON * max(1.0, rho_b):
        issues.append(f"baryon conservation residual {baryon_resid:.3e} exceeds tol")

    # 2. cell charge neutrality (rho_iip must be 0)
    charge_resid = abs(rho_e - (f*rho_ip + (1.0 - f)*rho_iip))
    if charge_resid > TOL_CHARGE * max(1.0, rho_e):
        issues.append(f"charge neutrality residual {charge_resid:.3e} exceeds tol")
    if abs(rho_iip) > 1.0e-12:
        issues.append(f"gas proton density rho_IIp={rho_iip:.3e} must be 0 (pure neutron gas)")

    # 3. beta-equilibrium at the cluster
    beta_resid = abs(mu_n - mu_p - mu_e)
    if beta_resid > TOL_BETA:
        issues.append(f"beta-equilibrium residual {beta_resid:.3e} MeV exceeds tol {TOL_BETA}")

    # 4. physical bounds
    if not (0.0 < f < 1.0):
        issues.append(f"volume fraction f={f} out of (0,1)")
    if r_d <= 0.0:
        issues.append(f"cluster radius R_d={r_d} non-positive")
    if not (rho_i > rho_ii >= 0.0):
        issues.append(f"density ordering violated: rho_I={rho_i} rho_II={rho_ii}")
    yp = p["yp_cell"]
    if not (0.0 <= yp <= 1.0):
        issues.append(f"cell proton fraction yp={yp} out of [0,1]")
    if p["status"] != "ok":
        issues.append(f"point status={p['status']}")
    if not p["converged"]:
        issues.append("point not converged")
    return issues


def verify_model(points: list[dict], model_id: str) -> dict:
    all_issues = []
    for idx, p in enumerate(points):
        iss = verify_point(p, model_id)
        if iss:
            all_issues.append({"point": idx, "rho_b": p["rho_b_fm3"], "issues": iss})

    # 5. thermodynamic monotonicity of E_binding_per_baryon
    ebind = [p["E_binding_per_baryon_MeV"] for p in points]
    monotonic_violations = 0
    for i in range(1, len(ebind)):
        if ebind[i] < ebind[i-1] - 1.0e-6:
            monotonic_violations += 1

    return {
        "model_id": model_id,
        "n_points": len(points),
        "rho_min": RHO_MIN, "rho_max": RHO_MAX, "n_points_requested": N_POINTS,
        "sigma_0": SIGMA_0, "b": B_EXP,
        "invariant_checks": {
            "baryon_conservation_tol": TOL_BARYON,
            "charge_neutrality_tol": TOL_CHARGE,
            "beta_equilibrium_tol_mev": TOL_BETA,
            "point_violations": all_issues,
            "monotonic_e_binding_violations": monotonic_violations,
        },
        "all_invariants_pass": (len(all_issues) == 0 and monotonic_violations == 0),
        "first_point": points[0] if points else None,
        "last_point": points[-1] if points else None,
    }


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: build_e390_03_cld_native.py OUTPUT_DIRECTORY")
    out = Path(sys.argv[1]).expanduser().resolve()
    if out.exists():
        raise FileExistsError(f"refusing to overwrite existing attempt directory: {out}")
    out.mkdir(parents=True)

    exe = executable()
    exe_sha = sha256(exe)
    src_module = ROOT / "fortran_eos" / "src" / "compose_eos_cld.f90"
    rmf_module = ROOT / "fortran_eos" / "src" / "compose_eos_rmf.f90"
    python_cld_legacy = ROOT / "python_frontend" / "compose_eos" / "cld_cell.py"
    registry = ROOT / "data" / "eos390_rmf_parameter_sets.v1.json"

    model_receipts = []
    total_points = 0
    overall_pass = True
    for model_id in MODEL_IDS:
        points = run_native_cld(model_id)
        total_points += len(points)
        verification = verify_model(points, model_id)
        model_receipts.append(verification)
        if not verification["all_invariants_pass"]:
            overall_pass = False
        # save per-model CLD table
        write_json(out / f"{model_id}_cld_grid.json", points)
        write_json(out / f"{model_id}_cld_verification.json", verification)

    receipt = {
        "schema_version": "1.0.1",
        "node": "E390-03",
        "subscope": "native_CLD_Thomas_Fermi_composition_reconstruction",
        "status": "NATIVE_CLD_INVARIANT_PASS" if overall_pass else "NATIVE_CLD_INVARIANT_FAIL",
        "built_date": "2026-11-14",
        "models": list(MODEL_IDS),
        "model_receipts": model_receipts,
        "total_output_rows": total_points,
        "density_point_count_per_model": N_POINTS,
        "density_support_fm3": [RHO_MIN, RHO_MAX],
        "cld_parameters": {
            "sigma_0_MeV_fm2": SIGMA_0,
            "b_lattimer_exponent": B_EXP,
            "lattice_function": "Ravenhall-Pethick-Wilson dilute-sphere Phi(u)=0.5*sqrt(2-3u)",
            "surface_tension": "Lattimer empirical sigma=sigma_0*(yp*(1-yp))^b",
            "minimizer": "coarse grid scan + Nelder-Mead 3-variable (f, log R_d, rho_I)",
            "gas_model": "pure neutron gas (rho_IIp=0), standard low-density inner-crust limit",
        },
        "physics_contract": {
            "cld_cell": "Scurto 2024 PRD 109 103015 eq. 21-25 on the validated DD-RMF homogeneous solver",
            "cluster_beta_equilibrium": "mu_n(I) - mu_p(I) = mu_e(cell), bisection on rho_Ip in [0, rho_I]",
            "baryon_conservation": "rho_B = f*rho_I + (1-f)*rho_II (rho_II from conservation)",
            "cell_charge_neutrality": "rho_e = f*rho_Ip + (1-f)*rho_IIp with rho_IIp=0",
            "homogeneous_solver": "Fortran solve_ddrmf_homogeneous (Typel-Wolter DD-RMF, rearrangement included)",
            "rest_mass_convention": "E_cell uses eps_b (binding, rest-mass-free), matching Python binding_energy",
        },
        "invariants_checked": [
            "baryon conservation: rho_b == f*rho_I + (1-f)*rho_II",
            "cell charge neutrality: rho_e == f*rho_Ip + (1-f)*rho_IIp",
            "gas is pure neutron: rho_IIp == 0",
            "beta-equilibrium at cluster: |mu_n - mu_p - mu_e| < 1e-3 MeV",
            "physical bounds: 0<f<1, R_d>0, rho_I>rho_II>=0, yp_cell in [0,1]",
            "thermodynamic monotonicity: E_binding_per_baryon non-decreasing",
        ],
        "non_claims": [
            "This is not a reconstructed CompOSE EOS table.",
            "No core-crust matching, TOV, tidal or Cowling output is generated.",
            "The CLD gas phase is the standard low-density pure-neutron approximation (rho_IIp=0); "
            "full Thomas-Fermi Wigner-Seitz gas composition is a separate subscope.",
            "The Lattimer surface-tension empirical form sigma=sigma_0*(yp*(1-yp))^b is a parameterisation, "
            "not a microscopic derivation; sigma_0 and b are inputs, not fitted here.",
            "Python cld_cell.py uses a different (richer) RMF solver and is NOT a numerical reference for "
            "the native Fortran CLD; the Fortran DD-RMF homogeneous solver is the canonical backend.",
        ],
        "provenance": {
            "native_executable_sha256": exe_sha,
            "fortran_cld_module_sha256": sha256(src_module),
            "fortran_rmf_module_sha256": sha256(rmf_module),
            "python_cld_legacy_sha256": sha256(python_cld_legacy),
            "rmf_parameter_registry_sha256": sha256(registry),
            "parent_commit": subprocess.run(
                ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                capture_output=True, text=True).stdout.strip(),
        },
        "publication": {
            "backend": "Fortran compose_eos_cld module + compose_eos_rmf homogeneous solver",
            "frontend": "Python receipt/verification only; no numerical fallback",
        },
    }
    write_json(out / "E390-03_CLD_NATIVE_SOLVER_receipt.json", receipt)

    # human-readable summary
    summary_lines = [f"E390-03 CLD native reconstruction: {len(MODEL_IDS)} models x {N_POINTS} points = {total_points} rows"]
    for mr in model_receipts:
        v = mr["invariant_checks"]
        summary_lines.append(
            f"  {mr['model_id']}: pass={mr['all_invariants_pass']} "
            f"point_violations={len(v['point_violations'])} "
            f"monotonic_violations={v['monotonic_e_binding_violations']}")
    (out / "SUMMARY.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("\n".join(summary_lines))
    return 0 if overall_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
