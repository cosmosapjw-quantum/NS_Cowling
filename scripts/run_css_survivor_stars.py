"""Phase-segmented TOV and Cowling scans for CSS census survivors.

For every survivor of a census manifest: build the phase-segmented EOS with
compose_cowling.css_adapter, integrate an M(P_c) sequence, bisect for the
1.4 M_sun star on the rising branch, and scan real-frequency l=2 Cowling
roots for slow and rapid interface conversion.  Roots are returned
unlabelled (no f/p mode identity is asserted) together with node counts.

Scope: radii are truncation radii at the CLD crust table's minimum pressure
(no outer crust); composition is absent.

Usage: python scripts/run_css_survivor_stars.py CENSUS_DIR OUT_JSON
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_cowling.css_adapter import css_match_maxwell_eos  # noqa: E402
from compose_cowling.maxwell import (  # noqa: E402
    MEV_FM3_TO_J_M3, scan_maxwell_cowling_modes, solve_maxwell_tov)

TARGET_MASS = 1.4
FREQ_GRID_HZ = np.linspace(800.0, 4000.0, 25)


def star_row(eos, pc_mev):
    s = solve_maxwell_tov(eos, pc_mev * MEV_FM3_TO_J_M3)
    return s, {"p_c_MeV_fm3": pc_mev, "mass_solar": s.mass_solar,
               "truncation_radius_km": s.radius_m / 1e3,
               "interface_radius_km": s.inner.r_end_m / 1e3}


def _plain(v):
    if isinstance(v, (bool, str)) or v is None:
        return v
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return float(v)
    return str(v)


def main(census_dir: Path, out_json: Path) -> int:
    manifest = json.loads((census_dir / "matching_manifest.json").read_text())
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()
    out = {"schema": "ns-cowling.css-survivor-stars.v1", "git_head": head,
           "census_manifest_sha256": hashlib.sha256(
               (census_dir / "matching_manifest.json").read_bytes()).hexdigest(),
           "census_git_head": manifest.get("git_head"),
           "scope": ("phase-segmented reference TOV/Cowling; radius = truncation radius at the "
                     "CLD crust table minimum pressure (no outer crust); no composition; "
                     "Cowling roots unlabelled"),
           "frequency_grid_hz": FREQ_GRID_HZ.tolist(), "survivors": []}
    for s in manifest["survivors"]:
        pair = s["pair_id"]
        t0 = time.time()
        rec = {"pair_id": pair}
        try:
            eos, meta = css_match_maxwell_eos(census_dir / pair / "compose", s["match_report"])
            rec["adapter"] = meta
            rec["interface"] = {k: _plain(v) for k, v in eos.interface_thermo.items()}
            p_lo = eos.transition_pressure_pa / MEV_FM3_TO_J_M3
            p_hi = eos.high_pressure_max_pa / MEV_FM3_TO_J_M3
            grid = np.geomspace(max(1.5 * p_lo, 5.0), 0.95 * p_hi, 14)
            seq = []
            for pc in grid:
                try:
                    seq.append(star_row(eos, float(pc))[1])
                except (ValueError, RuntimeError, ArithmeticError) as exc:
                    seq.append({"p_c_MeV_fm3": float(pc), "error": f"{type(exc).__name__}: {exc}"})
            rec["sequence"] = seq
            ok = [r for r in seq if "mass_solar" in r]
            masses = [r["mass_solar"] for r in ok]
            imax = int(np.argmax(masses)) if masses else None
            rec["max_mass_on_grid_solar"] = masses[imax] if masses else None
            rising = ok[: imax + 1] if imax is not None else []
            bracket = next(((a, b) for a, b in zip(rising, rising[1:])
                            if a["mass_solar"] <= TARGET_MASS <= b["mass_solar"]), None)
            if bracket is None:
                rec["canonical_star"] = {"status": "1.4 Msun not on the rising branch of the grid"}
            else:
                lo, hi = bracket[0]["p_c_MeV_fm3"], bracket[1]["p_c_MeV_fm3"]
                star, row = None, None
                for _ in range(40):
                    mid = float(np.sqrt(lo * hi))
                    star, row = star_row(eos, mid)
                    if row["mass_solar"] < TARGET_MASS:
                        lo = mid
                    else:
                        hi = mid
                    if hi / lo - 1.0 < 1e-7:
                        break
                rec["canonical_star"] = row
                for conv in ("slow", "rapid"):
                    roots = scan_maxwell_cowling_modes(star, FREQ_GRID_HZ, conversion=conv,
                                                       max_step_m=200.0)
                    rec[f"cowling_{conv}"] = [
                        {k: _plain(v) for k, v in r.items() if not isinstance(v, np.ndarray)}
                        for r in roots]
            rec["status"] = "OK"
        except Exception as exc:  # recorded per survivor, never silently dropped
            rec["status"] = "ERROR"
            rec["error"] = f"{type(exc).__name__}: {exc}"
        rec["elapsed_s"] = round(time.time() - t0, 1)
        out["survivors"].append(rec)
        print(f"{pair:28s} {rec['status']} Mmax={rec.get('max_mass_on_grid_solar')} "
              f"star={rec.get('canonical_star', {}).get('truncation_radius_km')} "
              f"slow={len(rec.get('cowling_slow', []))} rapid={len(rec.get('cowling_rapid', []))} "
              f"{rec['elapsed_s']}s", flush=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(out, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1]), Path(sys.argv[2])))
