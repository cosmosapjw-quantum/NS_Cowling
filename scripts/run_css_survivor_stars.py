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
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import scipy

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
    """Serialize numeric evidence without turning nested diagnostics into reprs."""
    if v is None or isinstance(v, str):
        return v
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        if not np.isfinite(v):
            raise ValueError("non-finite numerical evidence")
        return float(v)
    if isinstance(v, dict):
        return {str(k): _plain(value) for k, value in v.items()}
    if isinstance(v, (list, tuple, np.ndarray)):
        return [_plain(value) for value in v]
    raise TypeError(f"unsupported evidence type: {type(v).__name__}")


def _validated_survivor_dir(census_dir: Path, survivor: dict) -> Path:
    """Bind the stellar input to the admitted census bytes before solving."""
    pair = survivor.get("pair_id")
    if not isinstance(pair, str) or not pair or Path(pair).name != pair or pair in (".", ".."):
        raise ValueError("invalid census pair_id")
    if survivor.get("method") != "css" or survivor.get("admission", {}).get("status") != "ADMITTED":
        raise ValueError("stellar scan requires an ADMITTED css survivor")
    directory = (census_dir / pair / "compose").resolve()
    if not directory.is_relative_to(census_dir.resolve()):
        raise ValueError("survivor input leaves census directory")
    hashes = survivor.get("output_sha256", {})
    required = {"eos.nb", "eos.t", "eos.yq", "eos.thermo"}
    actual = {p.name for p in directory.glob("eos.*") if p.is_file()}
    if not isinstance(hashes, dict) or not required.issubset(hashes) or actual != set(hashes):
        raise ValueError("survivor output sha256 inventory is incomplete or changed")
    for name, expected in hashes.items():
        if Path(name).name != name or name in (".", ".."):
            raise ValueError("invalid survivor output filename")
        file_path = directory / name
        if not file_path.resolve().is_relative_to(directory):
            raise ValueError("survivor output leaves compose directory")
        actual_hash = hashlib.sha256(file_path.read_bytes()).hexdigest()
        if actual_hash != expected:
            raise ValueError(f"survivor output sha256 mismatch: {name}")
    return directory


def _first_rising_bracket(sequence: list[dict], target_mass: float):
    """Bracket only within the first connected increasing sampled segment.

    A failed integration or sampled turning point terminates that segment.
    This is a grid selection rule, not a radial-stability certificate.
    """
    for a, b in zip(sequence, sequence[1:]):
        if "mass_solar" not in a or "mass_solar" not in b:
            return None
        if b["p_c_MeV_fm3"] <= a["p_c_MeV_fm3"] or b["mass_solar"] <= a["mass_solar"]:
            return None
        if a["mass_solar"] <= target_mass <= b["mass_solar"]:
            return a, b
    return None


def main(census_dir: Path, out_json: Path) -> int:
    if out_json.exists():
        raise FileExistsError(f"refusing to overwrite stellar receipt: {out_json}")
    manifest = json.loads((census_dir / "matching_manifest.json").read_text())
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()
    status = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--untracked-files=no"],
                            capture_output=True, text=True, check=True)
    sources = ["scripts/run_css_survivor_stars.py", "python_frontend/compose_cowling/css_adapter.py",
               "python_frontend/compose_cowling/css_barotrope.py", "python_frontend/compose_cowling/maxwell.py",
               "python_frontend/compose_cowling/direct_barotrope.py"]
    out = {"schema": "ns-cowling.css-survivor-stars.v1", "git_head": head,
           "git_tracked_tree_dirty": bool(status.stdout.strip()),
           "executed_source_sha256": {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sources},
           "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
           "census_manifest_sha256": hashlib.sha256(
               (census_dir / "matching_manifest.json").read_bytes()).hexdigest(),
           "census_git_head": manifest.get("git_head"),
           "scope": ("phase-segmented reference TOV/Cowling; radius = truncation radius at the "
                     "CLD crust table minimum pressure (no outer crust); no composition; "
                     "Cowling roots unlabelled; finite-pressure Delta P=0 boundary and "
                     "vacuum-normalized lapse at cutoff are an artificial truncation model; "
                     "frequencies are not complete-star observables or stability certification"),
           "frequency_grid_hz": FREQ_GRID_HZ.tolist(), "survivors": []}
    for s in manifest["survivors"]:
        pair = s["pair_id"]
        t0 = time.time()
        rec = {"pair_id": pair}
        try:
            compose_dir = _validated_survivor_dir(census_dir, s)
            rec["input_sha256"] = s["output_sha256"]
            eos, meta = css_match_maxwell_eos(compose_dir, s["match_report"])
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
            rec["mass_grid_scope"] = "sampled M(Pc); no maximum-mass or radial-stability certification"
            bracket = _first_rising_bracket(seq, TARGET_MASS)
            if bracket is None:
                rec["canonical_star"] = {"status": "NOT_FOUND_ON_GRID",
                                         "reason": "1.4 Msun not bracketed on sampled rising branch"}
                rec["status"] = "NO_CANONICAL_STAR_ON_GRID"
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
                    rec[f"cowling_scan_{conv}"] = {
                        "status": "ROOTS_FOUND_ON_GRID" if roots else "NO_ROOT_FOUND_ON_GRID",
                        "roots_found": len(roots),
                        "search_method": "sign changes and exact endpoint zeros",
                        "frequency_min_hz": float(FREQ_GRID_HZ[0]),
                        "frequency_max_hz": float(FREQ_GRID_HZ[-1]),
                        "spectral_completeness": "not_certified",
                    }
                    rec[f"cowling_{conv}"] = [
                        {**{k: _plain(v) for k, v in r.items() if k != "shot"},
                         "accepted_steps": r["shot"].accepted_steps}
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
    with out_json.open("x") as handle:
        handle.write(json.dumps(out, indent=2, allow_nan=False) + "\n")
    return int(any(r["status"] == "ERROR" for r in out["survivors"]))


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1]), Path(sys.argv[2])))
