#!/usr/bin/env python3
"""Execute an explicit native HF experiment plan without Python EOS arithmetic."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos.nv_hf_frontend import run_native_cell, sha256


def save(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=300)
    args = parser.parse_args()
    if args.workers not in range(1, 5) or not 0 < args.timeout <= 3600:
        raise ValueError("workers must be1..4 and timeout in(0,3600]")
    plan = json.loads(args.plan.read_text())
    if plan.get("schema") != "ns_cowling.nv_hf_experiment_plan.v1":
        raise ValueError("Unknown HF plan schema")
    jobs = plan["jobs"]
    ids = [job["id"] for job in jobs]
    if not jobs or len(ids) != len(set(ids)):
        raise ValueError("Empty or duplicate HF job identities")
    for job in jobs:
        if set(job) != {"id", "purpose", "config"} or not job["id"] or not all(c.isalnum() or c in "_-" for c in job["id"]):
            raise ValueError("Invalid job schema or identifier")
    args.out.mkdir(parents=True, exist_ok=False)
    files = [args.plan.resolve(), args.executable.resolve(), args.source.resolve(),
             Path(__file__).resolve(), ROOT / "python_frontend/compose_eos/nv_hf_frontend.py"]
    for relative in ("fortran_eos/src/compose_eos_nv_hf.f90", "fortran_eos/src/compose_eos_nv_dme.f90",
                     "fortran_eos/src/compose_eos_nv_original.f90", "fortran_eos/app/nv_hf_probe.f90"):
        files.append(ROOT / relative)
    before = {str(p): sha256(p) for p in files}
    settings = dict(schema="ns_cowling.nv_hf_source_campaign.v1", plan=plan,
                    source_head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                    source_status=subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True),
                    identities_before=before, workers=args.workers, timeout=args.timeout,
                    solver_threads_per_process=1)
    save(args.out / "SETTINGS.json", settings)
    started = time.perf_counter()
    results = []
    with (args.out / "RESULTS.jsonl").open("x") as journal, ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending = {pool.submit(run_native_cell, args.executable, args.source, args.out / "cells" / (job["id"] + ".txt"),
                               job["config"], timeout=args.timeout): job for job in jobs}
        for future in as_completed(pending):
            job = pending[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {"status": "FRONTEND_EXECUTION_FAILURE", "error": str(exc), "config": job["config"]}
            result.update(id=job["id"], purpose=job["purpose"])
            journal.write(json.dumps(result, allow_nan=False) + "\n")
            journal.flush()
            results.append(result)
            print(json.dumps({"id": result["id"], "status": result["status"], "completed": len(results), "total": len(jobs)}), flush=True)
    after = {str(p): sha256(p) for p in files}
    fields = ["id", "purpose", "status", "neutrons", "protons", "nb_fm3", "mesh", "seed_kind", "iterations",
              "energy_per_baryon_MeV", "pressure_MeV_fm3", "beta_residual_MeV", "density_residual",
              "energy_residual_MeV", "particle_error", "dilation_pressure_error_MeV_fm3", "wall_seconds"]
    with (args.out / "CELLS.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        for result in sorted(results, key=lambda x: x["id"]):
            native = result.get("native", {}).get("fields", {})
            writer.writerow({key: result.get(key, native.get(key, "")) for key in fields})
    unchanged = before == after and all(r.get("inputs_unchanged") for r in results)
    summary = dict(schema=settings["schema"], execution="COMPLETED", input_identities_unchanged=unchanged,
                   jobs=len(jobs), unique_jobs=len({r["id"] for r in results}),
                   statuses=dict(Counter(r["status"] for r in results)),
                   wall_seconds=time.perf_counter() - started, identities_after=after,
                   source_eos_admitted=False,
                   interpretation="Fixed-composition cell experiments; separate composition search, grid convergence and thermodynamic admission are required.")
    save(args.out / "SUMMARY.json", summary)
    print(json.dumps(summary, allow_nan=False), flush=True)
    if not unchanged:
        raise SystemExit(4)


if __name__ == "__main__":
    main()
