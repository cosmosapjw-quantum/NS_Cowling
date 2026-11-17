"""Regenerate CLD inner-crust tables with the corrected CompOSE writer.

Output: download/reconstruction_v6_proxy/<MODEL>/ (v2 kept rows on the optimizer
box bound; v3 trims them at the table ends and refuses interior ones; the historical
download/reconstruction/ tree is left untouched as evidence of the former
column-layout defect).  For each model the original density grid is reused,
the CLD minimization is rerun with the parameters that reproduce the stored
rows (sigma_0 = 1 MeV fm^-2, b = 1), and the table is emitted only if it is an
admissible cold proxy barotrope. The interface functional is exploratory,
not the source-calibrated Scurto functional; every receipt states this scope.  Otherwise ``REJECTED.json`` records why and no
CompOSE file is written.  SPG_M2/SPG_M4 parameter sets are placeholders equal
to SPG_M1 and are not regenerated (``NOT_GENERATED.json``).

Usage: python scripts/regenerate_cld_crusts_v2.py [OUT_ROOT]
OUT_ROOT must be absent or empty.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos.reconstruct_cld import (  # noqa: E402
    CLDThermodynamicsError, emit_compose_set, reconstruct_inner_crust)
from compose_eos.rmf_solver import load_models  # noqa: E402

SOURCE_ROOT = ROOT / "download" / "reconstruction"
PARAMS = ROOT / "data" / "eos390_rmf_parameter_sets.v1.json"
MODELS = ["SPG_M1", "SPG_M2", "SPG_M3", "SPG_M4", "SPG_M5",
          "GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW"]
PLACEHOLDERS = {"SPG_M2": "SPG_M1", "SPG_M4": "SPG_M1"}
SIGMA_0, B_EXP = 1.0, 1.0


def original_grid(model: str) -> list[float]:
    """Density grid of the historical table (eos.nb, else the misplaced Q1 slot)."""
    d = SOURCE_ROOT / model
    if (d / "eos.nb").is_file():
        return [float(x) for x in (d / "eos.nb").read_text().split()[2:]]
    rows = [l.split() for l in (d / "eos.thermo").read_text().splitlines()[1:] if l.strip()]
    return [float(r[3]) for r in rows]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(out_root: Path) -> int:
    if out_root.exists() and any(out_root.iterdir()):
        raise FileExistsError(f"output root is not empty: {out_root}")
    out_root.mkdir(parents=True, exist_ok=True)
    models = load_models(PARAMS)
    params_json = json.loads(PARAMS.read_text())
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()
    source_paths = ["python_frontend/compose_eos/cld_cell.py",
                    "python_frontend/compose_eos/reconstruct_cld.py",
                    "python_frontend/compose_eos/rmf_solver.py",
                    "scripts/regenerate_cld_crusts_v2.py"]
    source_hashes = {p: sha256(ROOT / p) for p in source_paths}
    tracked_dirty = bool(subprocess.run(["git", "-C", str(ROOT), "diff", "HEAD", "--name-only"],
                                       capture_output=True, text=True).stdout.strip())
    index = {"schema": "ns-cowling.cld-regeneration-index.v1", "git_head": head,
             "sigma_0": SIGMA_0, "b": B_EXP, "parameter_sets_sha256": sha256(PARAMS),
             "scientific_scope": "exploratory_proxy_CLD", "source_calibrated": False,
             "source_sha256": source_hashes, "git_tracked_dirty": tracked_dirty, "models": {}}
    for name in MODELS:
        out = out_root / name
        out.mkdir()
        if name in PLACEHOLDERS:
            rec = {"status": "NOT_GENERATED",
                   "reason": f"parameter set is a placeholder identical to {PLACEHOLDERS[name]} "
                             "(Scurto 2024 publishes M1..M5 couplings only as curves)"}
            (out / "NOT_GENERATED.json").write_text(json.dumps(rec, indent=2) + "\n")
            index["models"][name] = rec
            continue
        model = models[name]
        grid = original_grid(name)
        masses_recorded = (params_json.get("models", params_json).get(name, {})
                           .get("nucleon_masses_MeV"))
        provenance = {"sigma_0": SIGMA_0, "b": B_EXP, "density_grid_source":
                      str((SOURCE_ROOT / name).relative_to(ROOT)), "git_head": head,
                      "nucleon_masses_recorded_in_parameter_set": masses_recorded,
                      "hamiltonian_masses_used_MeV": {"neutron": model.mn, "proton": model.mp},
                      "source_sha256": source_hashes, "git_tracked_dirty": tracked_dirty}
        cells = reconstruct_inner_crust(model, SIGMA_0, B_EXP, densities=grid)
        try:
            receipt = emit_compose_set(cells, out, name, mn=model.mn, mp=model.mp,
                                       provenance=provenance)
        except CLDThermodynamicsError as exc:
            rec = {"status": "REJECTED", "reason": str(exc), "provenance": provenance,
                   "cells": [{"n_B": c.rho_B, "E_cell": c.E_cell, "y_p": c.y_p_cell,
                              "f": c.f, "R_d": c.R_d, "converged": c.converged} for c in cells]}
            (out / "REJECTED.json").write_text(json.dumps(rec, indent=2) + "\n")
            index["models"][name] = {"status": "REJECTED", "reason": str(exc)}
            print(f"{name:12s} REJECTED  {exc}")
            continue
        hashes = {f: sha256(out / f) for f in ("eos.nb", "eos.t", "eos.yq", "eos.thermo")}
        index["models"][name] = {"status": "EMITTED", "sha256": hashes,
                                 "dropped_box_bound_rows": receipt["dropped_box_bound_rows"],
                                 "max_euler_relative_residual": receipt["max_euler_relative_residual"]}
        print(f"{name:12s} EMITTED   rows={receipt['rows']} dropped_box_bound_rows={receipt['dropped_box_bound_rows']}")
    (out_root / "regeneration_index.json").write_text(json.dumps(index, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "download" / "reconstruction_v6_proxy"
    raise SystemExit(main(target))
