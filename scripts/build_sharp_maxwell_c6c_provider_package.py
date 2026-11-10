#!/usr/bin/env python3
"""C6C: build compact provider-clarification evidence tables (read-only).

Reads the C6B row-audit CSVs and the pinned raw eos.compo blobs (git objects at
the source authority commit) and writes small sample/aggregate tables.
It performs no scientific reinterpretation, changes no production code, and
never writes into the C6B evidence directory.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

AUTHORITY = "a09ec30d74b363652d2d1aefe990557c468679dd"
SOURCES = {
    "SPG_M1": {
        "path": "eos_generator/eos_data/SPG_M1/eos.compo",
        "blob": "ab79da517271bea1f7de37294e87e019de611699",
        "rows": [1, 217, 218, 291, 365],
    },
    "GPPVA_FSU2": {
        "path": "eos_generator/eos_data/GPPVA_FSU2/eos.compo",
        "blob": "fa1ad2b809eb1c95825e2ee6eb2650cf3f248693",
        "rows": [1, 217, 218, 243, 269],
    },
}
C6B_DIR_DEFAULT = "artifacts/sharp_maxwell/C6B_COMPOSITION_SEMANTICS_20261110T100954Z"
SAMPLE_COLS = [
    "row", "iT", "inb", "iYq", "Iphase", "eos_thermo_extra_geometry", "Yn", "Yp", "Ye", "Ymu",
    "group_id", "Abar", "Zbar", "YI_raw", "Yh_C3_implied", "Yh_C4_implied",
    "r_C3_direct", "r_C4_direct", "Zbar_over_Abar",
]


def git_show(repo: Path, rev_path: str) -> bytes:
    return subprocess.run(["git", "-C", str(repo), "show", rev_path], check=True,
                          capture_output=True).stdout


def num(x: str):
    return float(x) if x not in ("", "NA") else None


def rng(vals):
    vals = [v for v in vals if v is not None]
    return [min(vals), max(vals)] if vals else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    ap.add_argument("--c6b-dir", default=C6B_DIR_DEFAULT)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    repo = Path(a.repo).resolve()
    out = Path(a.out_dir)
    if "C6B_COMPOSITION" in str(out):
        raise SystemExit("refusing to write into the C6B evidence directory")
    out.mkdir(parents=True, exist_ok=True)

    report = {"schema": "sharp_maxwell_c6c_attachment_tables.v1",
              "source_authority_commit": AUTHORITY, "sources": {}}
    for sid, cfg in SOURCES.items():
        raw = git_show(repo, f"{AUTHORITY}:{cfg['path']}")
        blob = subprocess.run(["git", "-C", str(repo), "rev-parse", f"{AUTHORITY}:{cfg['path']}"],
                              check=True, capture_output=True, text=True).stdout.strip()
        if blob != cfg["blob"]:
            raise SystemExit(f"{sid}: blob mismatch {blob}")
        lines = raw.decode().splitlines()
        with open(repo / a.c6b_dir / f"row_audit_{sid}.csv", newline="") as fh:
            rows = list(csv.DictReader(fh))
        by_row = {int(r["row"]): r for r in rows}
        if len(rows) != len(lines):
            raise SystemExit(f"{sid}: row count mismatch")

        # representative rows: raw text + audited fields
        sample = []
        for n in cfg["rows"]:
            r = by_row[n]
            sample.append({"raw_line": lines[n - 1].strip(), **{c: r[c] for c in SAMPLE_COLS}})
        with open(out / f"sample_rows_{sid}.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["raw_line"] + SAMPLE_COLS)
            w.writeheader()
            w.writerows(sample)

        # aggregates per Iphase
        agg = {}
        for ip in sorted({r["Iphase"] for r in rows}):
            g = [r for r in rows if r["Iphase"] == ip]
            c3 = [num(r["Yh_C3_implied"]) for r in g]
            c4 = [num(r["Yh_C4_implied"]) for r in g]
            ye = [num(r["Ye"]) for r in g]
            yi = [num(r["YI_raw"]) for r in g]
            ab = [num(r["Abar"]) for r in g]
            zb = [num(r["Zbar"]) for r in g]
            za = [num(r["Zbar_over_Abar"]) for r in g]
            gap = [abs(num(r["C3_vs_C4_implied_gap"])) for r in g]
            agg[f"Iphase={ip}"] = {
                "n_rows": len(g),
                "row_range": [min(int(r["row"]) for r in g), max(int(r["row"]) for r in g)],
                "n_Yn_and_Yp_zero": sum(1 for r in g if num(r["Yn"]) == 0 and num(r["Yp"]) == 0),
                "n_YI_equals_Ye_plus_Ymu_1e-9": sum(
                    1 for r in g if abs(num(r["YI_raw"]) - num(r["Ye"]) - num(r["Ymu"])) < 1e-9),
                "n_Zbar_over_Abar_matches_Ye_1e-3_rel": sum(
                    1 for r in g if abs(num(r["Yq_over_ZoverA"]) - 1.0) < 1e-3),
                "Abar_range": rng(ab), "Zbar_range": rng(zb), "Zbar_over_Abar_range": rng(za),
                "Ye_range": rng(ye), "YI_raw_range": rng(yi),
                "Yh_C3_implied_range": rng(c3), "Yh_C4_implied_range": rng(c4),
                "max_abs_C3_vs_C4_implied_gap": max(gap),
                "n_C3_C4_implied_consistent_1e-3": sum(1 for x in gap if x < 1e-3),
                "unaccounted_baryon_fraction_range_if_C3_implied_Yh": (
                    None if rng(c3) is None else [1 - rng(c3)[1], 1 - rng(c3)[0]]),
                "max_abs_r_C3_direct": max(abs(num(r["r_C3_direct"])) for r in g),
                "max_abs_r_C4_direct": max(abs(num(r["r_C4_direct"])) for r in g),
                "geometry_codes": sorted({r["eos_thermo_extra_geometry"] for r in g}),
            }
        report["sources"][sid] = {
            "raw_path_at_authority_commit": cfg["path"],
            "git_blob_sha1": blob,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "bytes": len(raw), "rows": len(lines),
            "aggregates": agg,
        }
    (out / "aggregate_statistics.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v["sha256"] for k, v in report["sources"].items()}))


if __name__ == "__main__":
    main()
