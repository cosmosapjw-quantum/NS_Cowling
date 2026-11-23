"""Run the native layered TOV solver on indexed, unpacked CompOSE tables."""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--tables", type=Path, required=True, help="directory containing cell_ID/eos.*")
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--causality", type=Path, help="public audit keyed on census_cell_id")
    parser.add_argument("--private-only", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.index.open(newline="") as handle:
        selected = [row for row in csv.DictReader(handle)
                    if all(row[key] == "True" for key in
                           ("qualified", "has_BPS", "has_synthetic_bridge", "has_genuine_NV"))
                    and (row["private"] == "True") == args.private_only]
    if not selected:
        parser.error("no qualified inputs in requested visibility class")
    audited = {}
    if args.causality:
        with args.causality.open(newline="") as handle:
            audited = {row["census_cell_id"]: row for row in csv.DictReader(handle)}
        if any(row["census_cell_id"] not in audited for row in selected):
            parser.error("causal interval index missing for selected cell")
    native = args.binary.resolve(strict=True)
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    fields = ("M_msun", "R_km", "compactness", "k2", "Lambda", "y_surface",
              "central_pressure_MeV_fm3", "central_eps_MeV_fm3", "n_steps", "transition_count",
              "surface_pressure_MeV_fm3", "surface_epsilon_MeV_fm3", "surface_jump_scale")
    meta = ("census_cell_id", "core_id", "core_model", "method", "request", "max_dlogp",
            "first_bad_P_low_MeV_fm3", "first_bad_P_high_MeV_fm3", "causal_interval_status", "status")

    def run_case(row):
        identity = row["census_cell_id"]
        table = args.tables / identity
        if not (table / "eos.thermo").is_file():
            raise FileNotFoundError(table / "eos.thermo")
        audit = audited.get(identity, {})
        left, right = (audit.get("first_bad_P_" + bound + "_MeV_fm3", "") for bound in ("low", "high"))
        audit_options = (["--first-failure", left, right] if left else ["--audited-clean"]) if audit else []
        output = []
        for step in (".04", ".02", ".01"):
            trials = [(f"Pc{pc}", ["layered", str(table), pc, step]) for pc in ("25", "50", "100")]
            trials.append(("M1.4", ["layered-target", str(table), "1.4", "5", "300", step]))
            for request, cli in trials:
                proc = subprocess.run([str(native), "tov-love", *cli, *audit_options], env=env, capture_output=True,
                                      text=True, timeout=120, check=False)
                result = json.loads(proc.stdout)
                if proc.returncode and result.get("status") not in {"fail", "failed", "error"}:
                    raise RuntimeError(f"malformed failure from {identity}: {proc.stdout} {proc.stderr}")
                good = proc.returncode == 0 and result.get("converged") is True
                record = dict(census_cell_id=identity, core_id=row["core_id"], core_model=audit.get("core_model", ""),
                              method=row["method"], request=request, max_dlogp=step,
                              first_bad_P_low_MeV_fm3=left, first_bad_P_high_MeV_fm3=right)
                if good:
                    record.update({key: result[key] for key in fields})
                record["status"] = "ok" if good else result.get("reason", result.get("message", str(result)))
                record["causal_interval_status"] = result["causal_interval_status"] if good else "no_numerical_solution"
                output.append(record)
        return output

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.workers) as workers, args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=meta + fields, lineterminator="\n")
        writer.writeheader()
        for records in workers.map(run_case, selected):
            writer.writerows(records)


if __name__ == "__main__":
    main()
