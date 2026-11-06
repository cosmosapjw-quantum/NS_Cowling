#!/usr/bin/env python3
"""Portable, read-only adversarial audit of a directory of cold CompOSE EOS tables.

Example:
  python scripts/audit_eos_archive.py /path/to/eos_tables --json audit.json --csv audit.csv

Only directories containing eos.nb, eos.t, eos.yq and eos.thermo are considered.
The program never edits source tables.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PKG = HERE.parent / "compose_eos_pkg"
sys.path.insert(0, str(PKG))
from compose_eos.io import read_compose
from compose_eos.consistency import check_consistency


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def portable_root_label(root: Path) -> str:
    """Return a stable archive-root label without leaking a machine-local path."""
    return root.name or root.anchor or "."


def relative_source_path(path: Path, root: Path) -> str:
    """Represent source locations relative to the supplied audit root."""
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.name


def audit_one(directory: Path, root: Path) -> dict:
    data = read_compose(str(directory))
    r = check_consistency(data, verbose=False)
    return {
        "model": directory.name,
        "path": directory.relative_to(root).as_posix(),
        "n": data.n,
        "nb_min_fm-3": float(data.nb[0]),
        "nb_max_fm-3": float(data.nb[-1]),
        "T_MeV": float(data.temperature_mev),
        "Yq_grid": float(data.charge_fraction_grid),
        "Il": int(data.Il),
        "composition_available": bool(data.composition_available),
        "phase_ids": sorted(int(x) for x in np.unique(data.Iphase)) if data.composition_available else [],
        "Nadd": int(data.q_extra_count or 0),
        "checks": {
            k: {kk: vv for kk, vv in v.items() if kk not in {"cs2", "derivative_mask", "dPdnb", "dmubdnb", "depsdnb"}}
            for k, v in r.items() if k.startswith("C") or k == "overall" or k == "grid"
        },
        "C10_min_cs2": r["C10"]["min_cs2"],
        "C10_max_cs2": r["C10"]["max_cs2"],
        "C10_pass": r["C10"]["pass_"],
        "C2_max_rel": r["C2"]["max_rel"],
        "C3_max_abs": r["C3"]["max_abs"],
        "C4_max_abs": r["C4"]["max_abs"],
        "overall": r["overall"]["status"],
        "failed_checks": r["overall"]["failed_checks"],
        "file_sha256": {name: sha256(directory / name) for name in ("eos.nb", "eos.t", "eos.yq", "eos.thermo") + (("eos.compo",) if (directory / "eos.compo").is_file() else ())},
    }


def _json_safe(value):
    """Convert NaN/Inf to JSON null recursively; never emit nonstandard tokens."""
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input_root", type=Path)
    ap.add_argument("--json", type=Path, default=Path("eos_audit.json"))
    ap.add_argument("--csv", type=Path, default=Path("eos_audit.csv"))
    args = ap.parse_args()
    root = args.input_root.resolve()
    if not root.is_dir():
        ap.error(f"input_root is not a directory: {root}")
    dirs = sorted({p.parent for name in ("eos.nb",) for p in root.rglob(name)})
    if not dirs:
        ap.error(f"no eos.nb files found below {root}")
    results, errors = [], []
    for d in dirs:
        try:
            results.append(audit_one(d, root))
        except Exception as exc:
            errors.append({"model": d.name, "path": relative_source_path(d, root), "error": f"{type(exc).__name__}: {exc}"})
    payload = {
        "schema": "compose-eos-adversarial-audit-v3",
        "input_root": portable_root_label(root),
        "input_path_policy": "all table paths are relative to input_root; no machine-local absolute path is serialized",
        "tables_found": len(dirs),
        "audited": len(results),
        "errors": errors,
        "results": results,
        "limitations": [
            "C7-C9 and C11 are finite-grid derivative diagnostics, not analytic proofs.",
            "C10 is evaluated only on smooth-branch stencils; phase boundaries and large density gaps are masked.",
            "The audit does not establish model provenance or reproduce the original EOS generator without original raw inputs.",
        ],
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(_json_safe(payload), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    fields = ["model", "n", "nb_min_fm-3", "nb_max_fm-3", "composition_available", "C2_max_rel", "C3_max_abs", "C4_max_abs", "C10_min_cs2", "C10_max_cs2", "C10_pass", "overall", "failed_checks"]
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
        for row in results:
            w.writerow({k: json.dumps(row[k]) if isinstance(row[k], (list, dict)) else row.get(k) for k in fields})
    print(f"audited={len(results)} errors={len(errors)} JSON={args.json} CSV={args.csv}")
    for row in results:
        print(f"{row['model']}: {row['overall']} C2={row['C2_max_rel']:.3e} "
              f"cs2=[{row['C10_min_cs2']:.5g},{row['C10_max_cs2']:.5g}] "
              f"failed={row['failed_checks']}")
    return 0 if results and not errors else 2

if __name__ == "__main__":
    raise SystemExit(main())
