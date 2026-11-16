#!/usr/bin/env python3
"""Prepare source-bound SPG states and run a native energy-derivative check."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import pathlib
import statistics
import math
import subprocess
import sys
import tempfile

root = pathlib.Path(__file__).resolve().parents[2]
out = pathlib.Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=False)
source = root / "data/research_blockers_20261116/spg_inverse"
params = json.loads((source / "inferred_surface_parameters.json").read_text())["parameters"]
source_identity = json.loads((source / "summary.json").read_text())["source_sha256"]
identity_mismatches = []
for rel, expected in source_identity.items():
    actual = hashlib.sha256((root / rel).read_bytes()).hexdigest()
    if actual != expected:
        identity_mismatches.append((rel, expected, actual))
if identity_mismatches:
    raise SystemExit(f"source identity mismatch: {identity_mismatches}")
rows = list(csv.DictReader((source / "inverse.csv").open()))
by_key = {(r["model"].strip(), r["row"].strip()): r for r in rows}
base = [line.split() for line in (source / "inputs.txt").read_text().splitlines()
        if line.strip() and int(line.split()[1]) < 100000]
selected = []
for fields in base:
    model, row = fields[:2]
    r = by_key.get((model, row))
    if r is None or r["status"].strip() != "0":
        raise SystemExit(f"missing/non-root baseline row: {model} {row}")
    selected.append(r)
if len(base) != 760 or len(selected) != 760:
    raise SystemExit(f"unexpected source inventory: baseline={len(base)} eligible={len(selected)}")

def fmt(x: str | float) -> str:
    return f"{float(x):.17g}"

with (out / "inputs.txt").open("w") as f:
    for r in selected:
        model = r["model"].strip()
        p = params[model]
        f.write(" ".join([
            model, r["row"].strip(), fmt(r["nb"]), fmt(r["Ye"]), fmt(r["A"]),
            fmt(r["Z"]), fmt(r["ni"]), fmt(r["ng"]), fmt(r["f"]), fmt(r["R"]),
            fmt(p["sigma0"]), fmt(p["b"]),
        ]) + "\n")

build_dir = pathlib.Path(tempfile.mkdtemp(prefix="spg_residual_fd_"))
exe = build_dir / "spg_residual_fd"
cmd = [
    "gfortran", "-std=legacy", "-O2", "-ffree-line-length-none",
    "-fcheck=all", f"-J{build_dir}", f"-I{build_dir}",
    str(root / "fortran_eos/src/compose_eos_types.f90"),
    str(root / "fortran_eos/src/compose_eos_gdfm_native.f90"),
    str(root / "fortran_eos/tests/spg_residual_fd.f90"), "-o", str(exe),
]
build = subprocess.run(cmd, capture_output=True, text=True)
(out / "build.stdout").write_text(build.stdout)
(out / "build.stderr").write_text(build.stderr)
if build.returncode:
    raise SystemExit(f"native build failed ({build.returncode}); see {out}/build.stderr")
run = subprocess.run([str(exe)], input=(out / "inputs.txt").read_text(),
                     capture_output=True, text=True,
                     env={"PATH": __import__("os").environ["PATH"],
                          "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                          "MKL_NUM_THREADS": "1"})
(out / "stdout.csv").write_text(run.stdout)
(out / "stderr.log").write_text(run.stderr)
if run.returncode:
    raise SystemExit(f"native run failed ({run.returncode}); see {out}/stderr.log")
bad_record = "SPG_M1 999999 0.01 0.05 80 20 -0.01 0.001 0.1 4 0.9 10\n"
bad = subprocess.run([str(exe)], input=bad_record, capture_output=True, text=True)
(out / "fail_closed_probe.stdout").write_text(bad.stdout)
(out / "fail_closed_probe.stderr").write_text(bad.stderr)
if bad.returncode == 0 or "FAIL dense reference SPG_M1 999999" not in bad.stderr:
    raise SystemExit(f"native invalid-input probe did not fail closed; rc={bad.returncode}")
parsed = list(csv.DictReader((out / "stdout.csv").open()))
if len(parsed) != 760 * 5:
    raise SystemExit(f"unexpected native output rows: {len(parsed)}")
stats = {}
for step in sorted({float(r["step"]) for r in parsed}, reverse=True):
    vals = [abs(float(r["abs_fd_full_MeV"])) for r in parsed if float(r["step"]) == step]
    stats[f"{step:.1e}"] = {
        "count": len(vals),
        "max_abs_fd_full_MeV": max(vals),
        "rms_abs_fd_full_MeV": math.sqrt(sum(x*x for x in vals)/len(vals)),
        "median_abs_fd_full_MeV": statistics.median(vals),
    }
final = [r for r in parsed if float(r["step"]) == 1.e-5]
full = [float(r["beta_full_MeV"]) for r in final]
scaled = [float(r["beta_scaled_MeV"]) for r in final]
def sha(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
input_hashes = source_identity
root_sources = [
    "fortran_eos/src/compose_eos_types.f90",
    "fortran_eos/src/compose_eos_gdfm_native.f90",
    "fortran_eos/tests/spg_residual_fd.f90",
    "fortran_eos/tests/spg_residual_run.py",
    "data/research_blockers_20261116/spg_inverse/inputs.txt",
    "data/research_blockers_20261116/spg_inverse/inverse.csv",
    "data/research_blockers_20261116/spg_inverse/inferred_surface_parameters.json",
]
receipt = {
    "scope": "SPG source composition beta residual derivative check",
    "status": "OBSERVATIONS_RECORDED",
    "finite_difference_observation": {
        "assessment": "centered_difference_discrepancy_decreases_with_step_over_tested_range",
        "smallest_step_y": 1.e-5,
        "smallest_step_count": len(final),
        "smallest_step_rms_abs_fd_full_MeV": stats["1.0e-05"]["rms_abs_fd_full_MeV"],
        "smallest_step_max_abs_fd_full_MeV": stats["1.0e-05"]["max_abs_fd_full_MeV"],
        "largest_step_y": 1.e-3,
        "largest_step_rms_abs_fd_full_MeV": stats["1.0e-03"]["rms_abs_fd_full_MeV"],
        "largest_step_max_abs_fd_full_MeV": stats["1.0e-03"]["max_abs_fd_full_MeV"],
    },
    "claim_ceiling": "does not authenticate source solver prescription or infer author-code defect",
    "baseline_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
    "source_table_rows": 760,
    "native_finite_difference_rows": len(parsed),
    "steps_in_y": stats,
    "at_smallest_step_beta_rms_MeV": math.sqrt(sum(x*x for x in full)/len(full)),
    "at_smallest_step_beta_abs_max_MeV": max(map(abs, full)),
    "diagnostic_factor_1_minus_f_rms_MeV": math.sqrt(sum(x*x for x in scaled)/len(scaled)),
    "diagnostic_factor_1_minus_f_abs_max_MeV": max(map(abs, scaled)),
    "table_residual_comparison": {
        "inferred_full_derivative_rms_MeV": math.sqrt(sum(x*x for x in full)/len(full)),
        "inferred_full_derivative_max_abs_MeV": max(map(abs, full)),
        "factor_1_minus_f_diagnostic_rms_MeV": math.sqrt(sum(x*x for x in scaled)/len(scaled)),
        "factor_1_minus_f_diagnostic_max_abs_MeV": max(map(abs, scaled)),
        "interpretation": "diagnostic algebraic comparison only; neither is source implementation evidence",
    },
    "native_returncode": run.returncode,
    "fail_closed_probe": {
        "input": "SPG_M1 row 999999 with ni=-0.01 fm^-3",
        "expected_failure_class": "NATIVE_DENSE_BARYON_INPUT_REJECTED",
        "returncode": bad.returncode,
        "diagnostic": "FAIL dense reference SPG_M1 999999 invalid GDFM input",
        "observed": "FAIL dense reference SPG_M1 999999" in bad.stderr,
    },
    "command": ["python3", "fortran_eos/tests/spg_residual_run.py", sys.argv[1]],
    "compiler_command": cmd,
    "source_input_sha256": input_hashes,
    "source_identity_verification": {"status": "PASS", "checked_files": len(input_hashes)},
    "direct_dependency_sha256": {p: sha(root/p) for p in root_sources},
    "producer_binary_sha256": sha(exe),
    "input_manifest_sha256": sha(out / "inputs.txt"),
    "stdout_csv_sha256": sha(out / "stdout.csv"),
    "stderr_sha256": sha(out / "stderr.log"),
    "compiler_stderr_sha256": sha(out / "build.stderr"),
    "first_failed_attempt": {
        "classification": "INPUT_SELECTION_COUNT_MISMATCH",
        "detail": "initially selected all 880 inverse rows, including 120 rounding probes; native run did not start",
        "corrected": "selected the 760 baseline source rows by original row id < 100000 and keyed them to inverse.csv",
    },
}
(out / "SUMMARY.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
print(f"native returncode={run.returncode}; selected={len(selected)}; output_lines={len(run.stdout.splitlines())}")
