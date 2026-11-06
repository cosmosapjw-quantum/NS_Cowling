#!/usr/bin/env python3
"""Read-only cumulative acoustic-time diagnostic for an existing background.

This tool does not solve TOV or Cowling equations and does not alter the input
background.  It applies the same declared center-hold + trapezoid quadrature
used by the structure diagnostic to expose tau(r) as a cumulative series.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

C_SI = 299_792_458.0
DIAGNOSTIC_ID = "cumulative_acoustic_time_center_hold_v1"
REQUIRED = {"r_m", "chi", "nu", "lambda"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_background(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None or not REQUIRED.issubset(reader.fieldnames):
            raise ValueError("background_missing_required_columns")
        rows = []
        for i, row in enumerate(reader):
            try:
                values = {key: float(row[key]) for key in REQUIRED}
            except (TypeError, ValueError, KeyError) as exc:
                raise ValueError(f"background_non_numeric_row:{i}") from exc
            if not all(math.isfinite(v) for v in values.values()):
                raise ValueError(f"background_nonfinite_row:{i}")
            rows.append(values)
    if not rows:
        raise ValueError("background_empty")
    if rows[0]["r_m"] <= 0.0:
        raise ValueError("background_first_radius_must_be_positive")
    for i, row in enumerate(rows):
        if not (0.0 < row["chi"] <= 1.0):
            raise ValueError(f"background_chi_out_of_range:{i}")
        if i and row["r_m"] <= rows[i-1]["r_m"]:
            raise ValueError(f"background_radius_not_strictly_increasing:{i}")
    return rows


def cumulative_acoustic_time(rows: list[dict[str, float]]) -> list[dict[str, float]]:
    integrand = [math.exp(r["lambda"] - r["nu"]) / (C_SI * math.sqrt(r["chi"])) for r in rows]
    tau = rows[0]["r_m"] * integrand[0]
    out = [{"r_m": rows[0]["r_m"], "cumulative_acoustic_time_s": tau,
            "integrand_s_per_m": integrand[0]}]
    for i in range(1, len(rows)):
        dr = rows[i]["r_m"] - rows[i-1]["r_m"]
        tau += 0.5 * dr * (integrand[i-1] + integrand[i])
        out.append({"r_m": rows[i]["r_m"], "cumulative_acoustic_time_s": tau,
                    "integrand_s_per_m": integrand[i]})
    return out


def write_output(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["q", "r_m", "cumulative_acoustic_time_s", "integrand_s_per_m"])
        writer.writeheader()
        for i, row in enumerate(rows):
            writer.writerow({"q": i, "r_m": f"{row['r_m']:.17e}",
                             "cumulative_acoustic_time_s": f"{row['cumulative_acoustic_time_s']:.17e}",
                             "integrand_s_per_m": f"{row['integrand_s_per_m']:.17e}"})


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--background", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--metadata", required=True, type=Path)
    args = p.parse_args(argv)
    try:
        background = args.background.resolve(strict=True)
        rows = read_background(background)
        result = cumulative_acoustic_time(rows)
        write_output(args.output, result)
        meta = {
            "schema_version": 1,
            "diagnostic_id": DIAGNOSTIC_ID,
            "input": str(background),
            "input_sha256": sha256(background),
            "rows": len(result),
            "center_closure": "hold_first_accepted_integrand_over_0_to_r0",
            "quadrature": "trapezoid_on_stored_background",
            "terminal_radius_m": result[-1]["r_m"],
            "terminal_acoustic_time_s": result[-1]["cumulative_acoustic_time_s"],
            "output_sha256": sha256(args.output),
            "solver_modified": False,
        }
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    except (OSError, ValueError, OverflowError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
