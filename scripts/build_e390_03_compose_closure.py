"""E390-03 frontend: emit CompOSE-standard group abundances from the native
CLD reconstruction and verify the baryon/charge closure invariants.

The native Fortran compose_eos_cld module owns all CLD cell physics (cluster
density, gas density, volume fraction, beta-equilibrium).  This module only
maps the native output to the CompOSE eos.compo group-abundance tuple and
verifies the C3/C4 closure invariants defined in the E390-02A receipt:

  C3 (baryon closure):  Yn + Yp + Yh == 1   (constructed from baryon conservation)
  C4 (charge closure):  Yp + (Zbar/Abar)*Yh == Ye + Ymu
                        (with Ye = rho_e/rho_B, Ymu = 0 in the cold no-muon inner
                         crust, and the cluster providing Zbar/Abar = rho_Ip/rho_I)

CLD cell -> CompOSE group mapping (one heavy group = the cluster phase):
  Yh  = f * rho_I / rho_B          (heavy-group baryon fraction)
  Yn  = (1-f) * rho_II / rho_B     (free-neutron fraction; gas is pure neutron)
  Yp  = (1-f) * rho_IIp / rho_B    (free-proton fraction; = 0 for pure-neutron gas)
  Ye  = rho_e / rho_B              (electron fraction)
  Ymu = 0                          (no muon in the cold inner crust)
  Abar = rho_I / n_b_cluster       (cluster average mass number; see note)
  Zbar = rho_Ip / n_b_cluster      (cluster average charge number)
  where n_b_cluster = rho_I / Abar so that Abar = rho_I * V_cluster / N_b,cluster.

Because the CLD model treats the cluster as a continuous liquid-drop phase
(not a discrete nucleus with a fixed nucleon count), Abar is not a direct CLD
output.  The closure C3 (Yn+Yp+Yh=1) is constructed from baryon conservation
and holds exactly by construction, independent of Abar.  The C4 charge
closure uses Zbar/Abar = rho_Ip/rho_I (the cluster proton fraction), which is
a direct CLD output and is independent of the cluster nucleon count.

This module therefore:
  1. reads the native CLD grid for each model,
  2. emits the CompOSE group-abundance tuple per density point,
  3. verifies C3 and C4 closure to machine precision,
  4. writes the per-model composition table + closure receipt.

Python owns mapping, closure verification and receipt emission only; the
native Fortran backend owns all numerical physics.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

MODEL_IDS = ("GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW")

C3_TOL = 1.0e-9   # Yn + Yp + Yh == 1, constructed from baryon conservation
C4_TOL = 1.0e-9   # charge closure


def compose_group_tuple(p: dict) -> dict:
    """Map one native CLD point to the CompOSE eos.compo group-abundance tuple."""
    rho_b = p["rho_b_fm3"]
    f = p["f"]
    rho_i = p["rho_I_fm3"]
    rho_ip = p["rho_Ip_fm3"]
    rho_ii = p["rho_II_fm3"]
    rho_iip = p["rho_IIp_fm3"]
    rho_e = p["rho_e_fm3"]

    if rho_b <= 0.0:
        raise ValueError(f"non-positive rho_b={rho_b}")

    # heavy group (cluster) baryon fraction
    yh = f * rho_i / rho_b
    # free neutron fraction (gas is pure neutron: rho_iip = 0)
    yn = (1.0 - f) * rho_ii / rho_b
    # free proton fraction
    yp = (1.0 - f) * rho_iip / rho_b
    # electron fraction
    ye = rho_e / rho_b
    # muon fraction (cold inner crust: no muons)
    ymu = 0.0

    # cluster proton fraction (charge-to-baryon ratio of the heavy group)
    yp_cluster = rho_ip / rho_i if rho_i > 0.0 else 0.0
    # Abar / Zbar: the CLD model treats the cluster as a continuous phase, so
    # the absolute nucleon count is not a CLD output.  We record Abar as the
    # cluster baryon density normalised to a unit number density, which makes
    # Zbar/Abar = rho_Ip/rho_I exact and the C4 closure exact by construction.
    # Abar itself is left as a placeholder (NaN) since the CLD model does not
    # fix the cluster nucleon count; the C3/C4 closures do not require it.
    return {
        "rho_b_fm3": rho_b,
        "Iphase": 1,            # single inner-crust phase
        "Yn": yn,
        "Yp": yp,
        "Yh": yh,
        "Ye": ye,
        "Ymu": ymu,
        "Abar": float("nan"),   # CLD continuous-phase model: cluster nucleon count not fixed
        "Zbar": float("nan"),   # ditto; Zbar/Abar = rho_Ip/rho_I is the exact charge ratio
        "cluster_proton_fraction": yp_cluster,
        "cluster_baryon_density_fm3": rho_i,
        "gas_baryon_density_fm3": rho_ii,
        "volume_fraction": f,
        "provenance": "native_fortran_cld_reconstruct",
    }


def verify_closure(rows: list[dict]) -> dict:
    """Verify C3 (baryon) and C4 (charge) closure for the composition table."""
    c3_violations = []
    c4_violations = []
    for idx, r in enumerate(rows):
        # C3: Yn + Yp + Yh == 1
        c3 = r["Yn"] + r["Yp"] + r["Yh"]
        if abs(c3 - 1.0) > C3_TOL:
            c3_violations.append({"point": idx, "rho_b": r["rho_b_fm3"],
                                  "sum": c3, "residual": abs(c3 - 1.0)})
        # C4: Yp + (Zbar/Abar)*Yh == Ye + Ymu
        # Zbar/Abar = cluster_proton_fraction (exact CLD output)
        z_over_a = r["cluster_proton_fraction"]
        c4_lhs = r["Yp"] + z_over_a * r["Yh"]
        c4_rhs = r["Ye"] + r["Ymu"]
        if abs(c4_lhs - c4_rhs) > C4_TOL:
            c4_violations.append({"point": idx, "rho_b": r["rho_b_fm3"],
                                  "lhs": c4_lhs, "rhs": c4_rhs,
                                  "residual": abs(c4_lhs - c4_rhs)})
    return {
        "C3_baryon_closure": {
            "invariant": "Yn + Yp + Yh == 1",
            "tolerance": C3_TOL,
            "n_points": len(rows),
            "n_violations": len(c3_violations),
            "violations": c3_violations,
            "pass": len(c3_violations) == 0,
        },
        "C4_charge_closure": {
            "invariant": "Yp + (Zbar/Abar)*Yh == Ye + Ymu",
            "tolerance": C4_TOL,
            "n_points": len(rows),
            "n_violations": len(c4_violations),
            "violations": c4_violations,
            "pass": len(c4_violations) == 0,
            "note": "Zbar/Abar = cluster_proton_fraction = rho_Ip/rho_I (exact CLD output); "
                    "Abar/Zbar absolute values are NaN because the CLD continuous-phase model "
                    "does not fix the cluster nucleon count, but the charge ratio is exact.",
        },
    }


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: build_e390_03_compose_closure.py OUTPUT_DIRECTORY")
    out = Path(sys.argv[1]).expanduser().resolve()
    if out.exists():
        raise FileExistsError(f"refusing to overwrite existing attempt directory: {out}")
    out.mkdir(parents=True)

    cld_dir = ROOT / "data" / "receipts" / "eos390_m3" / "cld"
    model_receipts = []
    overall_pass = True
    total_rows = 0
    for model_id in MODEL_IDS:
        cld_path = cld_dir / f"{model_id}_cld_grid.json"
        if not cld_path.is_file():
            raise FileNotFoundError(f"missing native CLD grid: {cld_path}")
        cld_points = json.loads(cld_path.read_text(encoding="utf-8"))
        rows = [compose_group_tuple(p) for p in cld_points]
        total_rows += len(rows)
        closure = verify_closure(rows)
        model_pass = closure["C3_baryon_closure"]["pass"] and closure["C4_charge_closure"]["pass"]
        if not model_pass:
            overall_pass = False
        write_json(out / f"{model_id}_compose_compo.json", rows)
        write_json(out / f"{model_id}_closure_verification.json", closure)
        model_receipts.append({
            "model_id": model_id,
            "n_points": len(rows),
            "C3_pass": closure["C3_baryon_closure"]["pass"],
            "C4_pass": closure["C4_charge_closure"]["pass"],
            "C3_violations": closure["C3_baryon_closure"]["n_violations"],
            "C4_violations": closure["C4_charge_closure"]["n_violations"],
        })

    receipt = {
        "schema_version": "1.0.1",
        "node": "E390-03",
        "subscope": "compose_group_abundances_closure",
        "status": "COMPOSE_CLOSURE_PASS" if overall_pass else "COMPOSE_CLOSURE_FAIL",
        "built_date": "2026-11-14",
        "models": list(MODEL_IDS),
        "model_receipts": model_receipts,
        "total_output_rows": total_rows,
        "closure_invariants": {
            "C3_baryon": "Yn + Yp + Yh == 1 (constructed from baryon conservation)",
            "C4_charge": "Yp + (Zbar/Abar)*Yh == Ye + Ymu (Zbar/Abar = rho_Ip/rho_I exact)",
            "C3_tolerance": C3_TOL,
            "C4_tolerance": C4_TOL,
        },
        "mapping": {
            "Yh": "f * rho_I / rho_B  (heavy-group baryon fraction = cluster)",
            "Yn": "(1-f) * rho_II / rho_B  (free neutron, gas is pure neutron)",
            "Yp": "(1-f) * rho_IIp / rho_B  (= 0 for pure-neutron gas)",
            "Ye": "rho_e / rho_B  (electron fraction)",
            "Ymu": "0  (cold inner crust, no muons)",
            "Abar": "NaN (CLD continuous-phase model does not fix cluster nucleon count)",
            "Zbar": "NaN (ditto; Zbar/Abar = rho_Ip/rho_I is the exact charge ratio)",
        },
        "non_claims": [
            "Abar/Zbar absolute values are not produced; the CLD continuous-phase model "
            "treats the cluster as a liquid-drop phase without a fixed nucleon count. "
            "The C3/C4 closures use only the ratio Zbar/Abar = rho_Ip/rho_I, which is exact.",
            "This is a composition table, not a complete CompOSE eos.thermo/eos.compo pair; "
            "the thermodynamic table comes from the native DD-RMF + CLD reconstruction.",
            "No core-crust matching, TOV, tidal or Cowling output is generated.",
        ],
        "provenance": {
            "native_cld_grids": str(cld_dir),
            "parent_commit": subprocess_run(["git", "-C", str(ROOT), "rev-parse", "HEAD"]),
        },
    }
    write_json(out / "E390-03_COMPOSE_CLOSURE_receipt.json", receipt)
    print(f"E390-03 CompOSE group-abundance closure: {len(MODEL_IDS)} models, {total_rows} rows")
    for mr in model_receipts:
        print(f"  {mr['model_id']}: C3={mr['C3_pass']} C4={mr['C4_pass']} "
              f"(violations: C3={mr['C3_violations']} C4={mr['C4_violations']})")
    return 0 if overall_pass else 1


def subprocess_run(cmd):
    import subprocess
    return subprocess.run(cmd, capture_output=True, text=True).stdout.strip()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
