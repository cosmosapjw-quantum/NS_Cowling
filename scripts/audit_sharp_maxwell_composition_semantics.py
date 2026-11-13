#!/usr/bin/env python3
"""Read-only C6B composition-semantics audit for SPG_M1 and GPPVA_FSU2.

Independent of the C5 reader/C6 helper parser: parses eos.compo from raw text with
its own grammar, validates grid keys against eos.thermo/eos.nb, and tabulates the
candidate Y_h interpretations row by row.  Every candidate is a *diagnostic*; none
is an approved mapping.  No raw file is written; no phase-dependent branch is used
to compute any quantity (phase is only a grouping label for aggregation).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data" / "sharp_maxwell_65_cohort_manifest.v1.json"
DEFAULT_SOURCE_ROOT = ROOT / "eos_generator" / "eos_data"
SOURCES = ("SPG_M1", "GPPVA_FSU2")
AUX_FILES = ("eos.nb", "eos.thermo", "eos.compo", "eos.t", "eos.yq")
PARTICLES = {10: "Yn", 11: "Yp", 0: "Ye", 1: "Ymu"}
TOL = 2.0e-5  # native default composition tolerance, used only to label residuals


def git_blob_sha1(data: bytes) -> str:
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    data = path.read_bytes()
    return {"bytes": len(data), "git_blob_sha1": git_blob_sha1(data), "sha256": hashlib.sha256(data).hexdigest()}


def parse_compo_line(line: str, line_no: int) -> dict[str, Any]:
    """Strict positional grammar: iT inb iYq Iphase Npairs (id Y)*Npairs Nquad (id A Z Y)*Nquad."""
    tok = line.split()
    if len(tok) < 6:
        raise ValueError(f"line {line_no}: too few fields")
    rec: dict[str, Any] = {"iT": int(tok[0]), "inb": int(tok[1]), "iYq": int(tok[2]),
                           "Iphase": int(tok[3]), "Npairs": int(tok[4])}
    pos = 5
    pairs = []
    for _ in range(rec["Npairs"]):
        pairs.append((int(tok[pos]), float(tok[pos + 1])))
        pos += 2
    rec["pairs"] = pairs
    rec["Nquad"] = int(tok[pos])
    pos += 1
    quads = []
    for _ in range(rec["Nquad"]):
        quads.append((int(tok[pos]), float(tok[pos + 1]), float(tok[pos + 2]), float(tok[pos + 3])))
        pos += 4
    rec["quads"] = quads
    if pos != len(tok):
        raise ValueError(f"line {line_no}: {len(tok) - pos} trailing fields")
    rec["raw_field_count"] = len(tok)
    return rec


def candidate_row(rec: dict[str, Any]) -> dict[str, Any]:
    """Explicit candidate interpretations for one parsed row (no phase branches)."""
    ids = [p[0] for p in rec["pairs"]]
    frac = {PARTICLES.get(i, f"id{i}"): v for i, v in rec["pairs"]}
    yn, yp, ye, ym = (frac.get(k, 0.0) for k in ("Yn", "Yp", "Ye", "Ymu"))
    gid, a, z, yi = rec["quads"][0]
    out = {
        "Yn": yn, "Yp": yp, "Ye": ye, "Ymu": ym, "group_id": gid, "Abar": a, "Zbar": z, "YI_raw": yi,
        "particle_ids": " ".join(map(str, ids)),
        "Yq_charge": ye + ym,
        "Yh_direct": yi,
        "Yh_standard_Abar_YI": a * yi,
        "Yh_C3_implied": (a / z) * (ye + ym - yp) if z > 0 else None,
        "Yh_C4_implied": 1.0 - yn - yp,
        "r_C3_direct": yp + (z / a) * yi - ye - ym if a > 0 else None,
        "r_C4_direct": yn + yp + yi - 1.0,
        "r_C3_standard": yp + z * yi - ye - ym,
        "r_C4_standard": yn + yp + a * yi - 1.0,
    }
    # diagnostics (not candidate mappings): relations between raw fields
    out["YI_minus_Yq"] = yi - (ye + ym)
    out["Zbar_over_Abar"] = z / a if a > 0 else None
    out["Yq_over_ZoverA"] = (ye + ym) / (z / a) if a > 0 and z > 0 else None  # baryon fraction carried by cell charge
    out["C3_vs_C4_implied_gap"] = (out["Yh_C3_implied"] - out["Yh_C4_implied"]) if out["Yh_C3_implied"] is not None else None
    out["direct_C3_pass"] = abs(out["r_C3_direct"]) <= TOL if out["r_C3_direct"] is not None else None
    out["direct_C4_pass"] = abs(out["r_C4_direct"]) <= TOL
    out["standard_C3_pass"] = abs(out["r_C3_standard"]) <= TOL
    out["standard_C4_pass"] = abs(out["r_C4_standard"]) <= TOL
    return out


def _numbers(path: Path) -> list[list[str]]:
    return [ln.split() for ln in path.read_text().splitlines() if ln.strip()]


def parse_grid_keys(source_dir: Path) -> dict[str, Any]:
    nb_tok = _numbers(source_dir / "eos.nb")
    nb_count = int(nb_tok[1][0])
    nb_vals = [float(t) for row in nb_tok[2:] for t in row]
    thermo = _numbers(source_dir / "eos.thermo")
    th_keys = [(int(r[0]), int(r[1]), int(r[2])) for r in thermo[1:]]
    return {"nb_declared": nb_count, "nb_found": len(nb_vals), "nb_strictly_increasing":
            all(b > a for a, b in zip(nb_vals, nb_vals[1:])), "thermo_keys": th_keys,
            "t_lines": [" ".join(r) for r in _numbers(source_dir / "eos.t")],
            "yq_lines": [" ".join(r) for r in _numbers(source_dir / "eos.yq")],
            "thermo_header": " ".join(thermo[0])}


def audit_source(source_id: str, source_dir: Path, expected_blobs: dict[str, str]) -> dict[str, Any]:
    identity = {f: file_identity(source_dir / f) for f in AUX_FILES}
    mismatches = {f: (identity[f]["git_blob_sha1"], sha) for f, sha in expected_blobs.items()
                  if identity[f]["git_blob_sha1"] != sha}
    if mismatches:
        raise SystemExit(f"{source_id}: pinned identity mismatch {mismatches}")
    lines = [ln for ln in (source_dir / "eos.compo").read_text().splitlines() if ln.strip()]
    recs = [parse_compo_line(ln, i + 1) for i, ln in enumerate(lines)]
    grid = parse_grid_keys(source_dir)
    keys = [(r["iT"], r["inb"], r["iYq"]) for r in recs]
    checks = {
        "rows": len(recs),
        "compo_keys_equal_thermo_keys_in_order": keys == grid["thermo_keys"],
        "inb_is_1_to_N_consecutive": [k[1] for k in keys] == list(range(1, len(keys) + 1)),
        "rows_equal_nb_count": len(recs) == grid["nb_declared"] == grid["nb_found"],
        "nb_strictly_increasing": grid["nb_strictly_increasing"],
        "keys_unique": len(set(keys)) == len(keys),
        "iT_values": sorted({k[0] for k in keys}), "iYq_values": sorted({k[2] for k in keys}),
        "Npairs_values": sorted({r["Npairs"] for r in recs}), "Nquad_values": sorted({r["Nquad"] for r in recs}),
        "particle_id_orders": sorted({" ".join(str(p[0]) for p in r["pairs"]) for r in recs}),
        "group_ids": sorted({q[0] for r in recs for q in r["quads"]}),
        "field_counts": sorted({r["raw_field_count"] for r in recs}),
        "eos_t_lines": grid["t_lines"], "eos_yq_lines": grid["yq_lines"], "thermo_header": grid["thermo_header"],
    }
    thermo = _numbers(source_dir / "eos.thermo")
    nadd = int(thermo[0][2]) if len(thermo[0]) > 2 else 0
    rows = []
    for r, tr in zip(recs, thermo[1:]):
        c = candidate_row(r)
        extra = tr[11] if nadd >= 1 and len(tr) > 11 else ""
        rows.append({"source_id": source_id, "row": r["inb"], "iT": r["iT"], "inb": r["inb"], "iYq": r["iYq"],
                     "Iphase": r["Iphase"], "eos_thermo_extra_geometry": extra, "Npairs": r["Npairs"], "Nquad": r["Nquad"],
                     **c, "source_compo_blob_sha1": identity["eos.compo"]["git_blob_sha1"],
                     "parser": "audit_sharp_maxwell_composition_semantics.py:parse_compo_line"})
    return {"identity": identity, "checks": checks, "rows": rows, "nadd": nadd}


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    def agg(sub: list[dict[str, Any]]) -> dict[str, Any]:
        def mx(k: str) -> float | None:
            v = [abs(r[k]) for r in sub if r[k] is not None]
            return max(v) if v else None
        return {"n": len(sub),
                "max_abs_r_C3_direct": mx("r_C3_direct"), "max_abs_r_C4_direct": mx("r_C4_direct"),
                "max_abs_r_C3_standard": mx("r_C3_standard"), "max_abs_r_C4_standard": mx("r_C4_standard"),
                "n_direct_C3_pass": sum(1 for r in sub if r["direct_C3_pass"]),
                "n_direct_C4_pass": sum(1 for r in sub if r["direct_C4_pass"]),
                "n_standard_C3_pass": sum(1 for r in sub if r["standard_C3_pass"]),
                "n_standard_C4_pass": sum(1 for r in sub if r["standard_C4_pass"]),
                "max_abs_YI_minus_Yq": mx("YI_minus_Yq"),
                "n_YI_equals_Yq_within_1e-9": sum(1 for r in sub if abs(r["YI_minus_Yq"]) <= 1e-9),
                "max_abs_C3_vs_C4_implied_gap": mx("C3_vs_C4_implied_gap"),
                "n_C3_C4_implied_consistent_within_1e-3": sum(1 for r in sub if r["C3_vs_C4_implied_gap"] is not None and abs(r["C3_vs_C4_implied_gap"]) <= 1e-3),
                "n_Yn_and_Yp_zero": sum(1 for r in sub if r["Yn"] == 0.0 and r["Yp"] == 0.0),
                "Yh_C3_implied_range": [min(r["Yh_C3_implied"] for r in sub if r["Yh_C3_implied"] is not None),
                                         max(r["Yh_C3_implied"] for r in sub if r["Yh_C3_implied"] is not None)],
                "Yq_over_ZoverA_range": [min(r["Yq_over_ZoverA"] for r in sub if r["Yq_over_ZoverA"] is not None),
                                          max(r["Yq_over_ZoverA"] for r in sub if r["Yq_over_ZoverA"] is not None)]}
    out: dict[str, Any] = {}
    for key in sorted({(r["Iphase"], r["eos_thermo_extra_geometry"]) for r in rows}, key=str):
        sub = [r for r in rows if (r["Iphase"], r["eos_thermo_extra_geometry"]) == key]
        out[f"Iphase={key[0]}|extra={key[1] or 'NA'}"] = agg(sub)
    out["ALL"] = agg(rows)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    cat = json.loads(args.manifest.read_text())["source_catalog"]["crust_models"]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    report: dict[str, Any] = {"schema": "sharp_maxwell_c6b_composition_audit.v1", "sources": {}}
    for sid in SOURCES:
        expected = {f: v["git_blob_sha"] for f, v in cat[sid]["files"].items()}
        res = audit_source(sid, args.source_root / sid, expected)
        rows = res["rows"]
        (args.out_dir / f"row_audit_{sid}.json").write_text(json.dumps(rows, indent=1, sort_keys=True) + "\n")
        with (args.out_dir / f"row_audit_{sid}.csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        report["sources"][sid] = {"identity": res["identity"], "schema_checks": res["checks"], "nadd": res["nadd"],
                                  "closure_summary": summarize(rows)}
    (args.out_dir / "closure_summary.json").write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print(json.dumps({s: v["closure_summary"]["ALL"]["n"] for s, v in report["sources"].items()}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
