#!/usr/bin/env python3
"""Stage E selected causal-tail equilibrium Cowling modes; reuse D consumers.

This runner does not extend the EOS, change D tolerances, or infer radial stability.
Its source-grid comparison is a measured representation sensitivity, not an error
bound. Every attempted configuration and rejected residual remains in the output.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import asdict, replace
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
import traceback
import numpy as np
import scipy
from compose_eos.stellar_background import ContractEnthalpyEOS, solve_background
from compose_eos.stellar_cowling import CowlingSettings, matching_residual, scan_brackets, solve_mode

ROOT = Path(__file__).resolve().parents[1]
BASE_COMMIT = "300e5a0532dd137cb09738568fca0f35cfcc9282"
RESIDUAL_LIMIT = 2e-6
REUSED_PATHS = [
    "python_frontend/compose_eos/stellar_background.py",
    "python_frontend/compose_eos/stellar_cowling.py",
    "python_frontend/tests/test_stellar_background.py",
    "python_frontend/tests/test_stellar_cowling.py",
    "data/stellar_consumers_20261117/cowling/VALIDATION.json",
    "docs/stellar_consumers_20261117/FORMULATION_KO.md",
    "docs/stellar_consumers_20261117/INDEPENDENT_REVIEW_KO.md",
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def nodes(values):
    values = np.asarray(values, dtype=float)
    kept = values[np.abs(values) > 1e-9*np.max(np.abs(values))]
    return int(np.count_nonzero(kept[:-1]*kept[1:] < 0))


def source_identity():
    rows = []
    for name in REUSED_PATHS:
        old = subprocess.check_output(["git", "show", BASE_COMMIT+":"+name], cwd=ROOT)
        current = digest(ROOT/name)
        old_digest = hashlib.sha256(old).hexdigest()
        if current != old_digest:
            raise RuntimeError("D reuse identity changed: "+name)
        rows.append(dict(path=name, sha256=current, baseline_sha256=old_digest,
                         relation="BYTE_IDENTITY", unchanged=True))
    return rows


def background_profile(bg):
    central = bg.eos.state_at_enthalpy(bg.h_c)
    phi_c = .5*np.log1p(-2*bg.mass_geom_km/bg.radius_km)-bg.h_c
    row = dict(location="analytic_center_limit", h=bg.h_c, r_km=0., m_geom_km=0.,
               n_fm3=central["n"], pressure_MeV_fm3=central["p"],
               energy_MeV_fm3=central["epsilon"], cs2_over_c2=central["cs2_over_c2"],
               Phi=phi_c, lambda_metric=0., dPhi_dr_km=0.)
    rows = [row]
    for h in np.linspace(bg.h_start, 0., 800):
        state = bg.state_at_enthalpy(h)
        rows.append(dict(location="vacuum_surface" if h == 0. else "integrated",
                         h=float(h), r_km=state["r_km"], m_geom_km=state["m_geom_km"],
                         n_fm3=state["n"], pressure_MeV_fm3=state["p"],
                         energy_MeV_fm3=state["epsilon"], cs2_over_c2=state["cs2_over_c2"],
                         Phi=state["Phi"], lambda_metric=state["lambda_metric"],
                         dPhi_dr_km=state["dPhi_dr_km"]))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", default="data/stellar_extension_20261117/tov/MODE_SELECTION.json")
    parser.add_argument("--contract")
    parser.add_argument("--comparison-contract")
    parser.add_argument("--models", help="Comma separated label:density_fm3 pairs")
    parser.add_argument("--output", default="data/stellar_extension_20261117/cowling")
    args = parser.parse_args()
    out = ROOT/args.output
    out.mkdir(parents=True, exist_ok=True)
    start_all = time.monotonic()
    selection_path = ROOT/args.selection
    selection = json.loads(selection_path.read_text()) if selection_path.exists() else {}
    contract_name = args.contract or selection.get("primary_contract")
    if not contract_name:
        raise ValueError("An explicit causal extension contract is required")
    contract_path = ROOT/contract_name
    comparison_name = args.comparison_contract or selection.get("comparison_contract")
    if args.models:
        models = [dict(label=item.split(":", 1)[0], central_density_fm3=float(item.split(":", 1)[1]),
                       reason="explicit selected Stage E background") for item in args.models.split(",")]
    else:
        models = selection["models"]
    if not 1 <= len(models) <= 2:
        raise ValueError("Stage E runner is bounded to one or two selected backgrounds")
    for model in models:
        label = model["label"]
        if not label or any(c not in "abcdefghijklmnopqrstuvwxyz0123456789_" for c in label):
            raise ValueError("Model labels must be lowercase ASCII letters, digits, or underscore")
    eos = ContractEnthalpyEOS.from_path(contract_path)
    comparison_eos = ContractEnthalpyEOS.from_path(ROOT/comparison_name) if comparison_name else None
    identity = source_identity()
    provenance = dict(
        schema="ns-cowling.extension-modes-provenance.v1", baseline_commit=BASE_COMMIT,
        run_head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        command=shlex.join([sys.executable, *sys.argv]), cwd=str(ROOT),
        python=sys.version, numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform(),
        thread_environment={name:os.environ.get(name) for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]},
        pythonpath=os.environ.get("PYTHONPATH"),
        runner=dict(path=relative(Path(__file__)), sha256=digest(Path(__file__))),
        primary_contract=dict(path=relative(contract_path), sha256=digest(contract_path), density_domain=list(eos.density_domain)),
        comparison_contract=(dict(path=comparison_name, sha256=digest(ROOT/comparison_name),
                                  density_domain=list(comparison_eos.density_domain)) if comparison_name else None),
        selection=(dict(path=relative(selection_path), sha256=digest(selection_path), content=selection) if selection_path.exists() else None),
        reused_D_controls=identity,
        reused_control_status="Previously executed D analytic and regression controls; source/tests/receipt byte identity verified, not rerun",
        residual_acceptance=RESIDUAL_LIMIT,
        selected_models=models)
    write_json(out/"PROVENANCE.json", provenance)
    base = CowlingSettings(rtol=1e-8, atol=1e-11, surface_h=1e-9, max_step_h=.002)
    tight = replace(base, rtol=1e-9, atol=1e-12)
    scan_settings = replace(base, rtol=1e-5, atol=1e-8, surface_h=1e-7, max_step_h=.005)
    configurations = [
        ("baseline", base, {}),
        ("ode_tight", tight, {}),
        ("surface_h1e_7", base.__class__(**{**asdict(base), "surface_h":1e-7}), {}),
        ("center_0p1m", base, dict(r0_km=.0001)),
    ]
    results = []
    comparison_skips = []
    for model in models:
        label = model["label"]
        nc = float(model["central_density_fm3"])
        directory = out/label
        directory.mkdir(exist_ok=True)
        bg = solve_background(eos, nc, rtol=1e-10, atol=1e-12)
        write_csv(directory/"background_profile.csv", background_profile(bg))
        brackets, scan = scan_brackets(bg, np.arange(250., 12001., 250.), scan_settings)
        write_json(directory/"SCAN.json", dict(settings=asdict(scan_settings), scan=scan, brackets=brackets))
        print(label, "nc", nc, "M", bg.mass_msun, "R", bg.radius_km, "brackets", brackets, flush=True)
        if len(brackets) < 2:
            raise RuntimeError("Two real mode brackets not found for "+label)
        todo = [(name, settings, changes, eos, "fine") for name, settings, changes in configurations]
        if comparison_eos is not None:
            if nc <= comparison_eos.density_domain[1]:
                todo.append(("source_grid_coarse", tight, {}, comparison_eos, "coarse"))
            else:
                comparison_skips.append(dict(model=label, reason="central density exceeds comparison contract's causal support", central_density_fm3=nc, comparison_upper_density_fm3=comparison_eos.density_domain[1]))
        model_rows = []
        for config_label, settings, changes, selected_eos, source_grid in todo:
            start = time.monotonic()
            kw = dict(rtol=1e-10, atol=1e-12)
            kw.update(changes)
            selected_bg = bg if selected_eos is eos and not changes else solve_background(selected_eos, nc, **kw)
            row = dict(model=label, label=config_label, source_grid=source_grid,
                       settings=asdict(settings), background=selected_bg.summary(), modes=[])
            for j, bracket in enumerate(brackets[:2]):
                name = "f" if j == 0 else "p1"
                try:
                    signs = [matching_residual(selected_bg, f, settings) for f in bracket]
                    mode = solve_mode(selected_bg, bracket, settings)
                    profile = mode.pop("profile")
                    profile_path = directory/(config_label+"_"+name+"_profile.csv")
                    write_csv(profile_path, profile)
                    profile_readback = list(csv.DictReader(profile_path.open()))
                    w_nodes = nodes([p["W"] for p in profile_readback])
                    v_nodes = nodes([p["V_km_inverse"] for p in profile_readback])
                    center = profile_readback[0]
                    center_ratio = -settings.ell*float(center["r_km"])*float(center["V_km_inverse"])/float(center["W"])
                    residual_ok = abs(mode["matching_residual"]) <= RESIDUAL_LIMIT and mode["relative_join_mismatch"] <= RESIDUAL_LIMIT
                    mode.update(name=name, bracket_endpoint_residuals=signs,
                                residual_accepted=bool(residual_ok),
                                root_and_nodes_accepted=bool(mode["root_converged"] and signs[0]*signs[1] < 0 and mode["W_nodes"] == j and mode["V_nodes"] == j),
                                csv_W_nodes=w_nodes, csv_V_nodes=v_nodes,
                                csv_center_regularity_ratio=center_ratio,
                                profile_path=relative(profile_path), profile_sha256=digest(profile_path))
                    mode["accepted"] = bool(mode["residual_accepted"] and mode["root_and_nodes_accepted"] and w_nodes == j and v_nodes == j)
                except Exception as error:
                    mode = dict(name=name, accepted=False, status="ERROR", error=repr(error), traceback=traceback.format_exc())
                row["modes"].append(mode)
            row["wall_seconds"] = time.monotonic()-start
            row["accepted"] = all(mode["accepted"] for mode in row["modes"])
            model_rows.append(row)
            results.append(row)
            write_json(directory/(config_label+".json"), row)
            write_json(out/"PROGRESS.json", dict(status="RUNNING", completed_configurations=results, comparison_skips=comparison_skips))
            print(label, config_label, [(m.get("frequency_hz"),m.get("W_nodes"),m.get("matching_residual"),m["accepted"]) for m in row["modes"]], "seconds", row["wall_seconds"], flush=True)
        reference = model_rows[0]["modes"]
        fine_tight = next(row for row in model_rows if row["label"] == "ode_tight")["modes"]
        for row in model_rows:
            reference_modes = fine_tight if row["source_grid"] == "coarse" else reference
            row["frequency_shifts_hz"] = [m.get("frequency_hz",0)-r.get("frequency_hz",0) if "frequency_hz" in m and "frequency_hz" in r else None for m,r in zip(row["modes"],reference_modes)]
            row["shift_reference"] = "fine_ode_tight" if row["source_grid"] == "coarse" else "fine_baseline"
            write_json(directory/(row["label"]+".json"), row)
    rejected = [dict(model=row["model"],configuration=row["label"],mode=mode) for row in results for mode in row["modes"] if not mode["accepted"]]
    all_primary = all(row["accepted"] for row in results if row["label"] == "baseline")
    output = dict(schema="ns-cowling.extension-modes.v1", stage="E", status="SELECTED_PRIMARY_MODES_ACCEPTED" if all_primary else "PRIMARY_MODE_REJECTION",
                  closure="equilibrium barotropic perfect fluid: Gamma1=Gamma_eq",
                  convention="(-+++), geometric km, Phi half temporal metric, frequency at infinity", ell=2,
                  center="regular series at finite r0; center-radius sensitivity tested",
                  surface="TOV h=0 vacuum surface; mode DeltaP=0 at positive hcut, surface sensitivity tested",
                  eos_interface="continuous W,V at continuous EOS knots; no sharp Maxwell jump",
                  residual_acceptance=RESIDUAL_LIMIT, models=models, configurations=results,
                  rejected_attempts=rejected, source_grid_comparison_skips=comparison_skips,
                  wall_seconds=time.monotonic()-start_all,
                  uncertainty_interpretation="Observed control and source-grid sensitivities; no rigorous global error bound",
                  scope="equilibrium l=2 Cowling f/p1 only; no radial-stability inference, frozen response, g-mode, elasticity, superfluidity, or full-GR modes",
                  provenance_path=relative(out/"PROVENANCE.json"))
    write_json(out/"SUMMARY.json", output)
    write_json(out/"PROGRESS.json", dict(status="COMPLETE", summary_path=relative(out/"SUMMARY.json"), rejected_attempt_count=len(rejected)))
    print("Saved",out/"SUMMARY.json", "status",output["status"], "rejected",len(rejected),flush=True)
    return 0 if all_primary else 2


if __name__ == "__main__":
    raise SystemExit(main())
