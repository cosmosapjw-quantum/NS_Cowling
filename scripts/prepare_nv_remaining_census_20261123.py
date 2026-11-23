#!/usr/bin/env python3
"""Enumerate unavailable-source admissions; never invoke numerical matching.

Python performs catalog, provenance, process scheduling and safe projection only.
The existing Cartesian runner remains byte-unchanged. Its native-call entry point
is denied in every worker as an additional invariant for this specific census.
All full catalogs, settings, identity manifests and runner output stay outside
the repository. Only an allowlisted, value-free model/status projection is public.
"""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import run_eos_all_pairs_20261123 as runner

ORIGINAL_INIT_WORKER = runner.init_worker
EXPECTED_MATCHER_SHA256 = "d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588"
CRUST_STATES = {
    "NV1973_DME_HF_printed": "QUARANTINED_TABLE_I_INTERNAL_INCONSISTENCY",
    "NV1973_DME_D03x100_diagnostic": "NOT_ADMITTED_DIAGNOSTIC_CELLS_NOT_GROUND_STATE_EOS",
    "BPS_UCRL11980_1965": "BLOCKED_MISSING_APPENDIX_D",
}
SOURCE_PROOFS = {
    "NV1973_DME_HF_printed": [
        "data/eos_nv_dme_20261123/source_contract/nv1972_table1_coefficients.json",
        "data/eos_nv_dme_20261123/source_contract/coefficient_consistency_audit.txt",
        "docs/NV1972_DME_FUNCTIONAL_20261123.md",
    ],
    "NV1973_DME_D03x100_diagnostic": [
        "data/eos_nv_dme_20261123/source_contract/FUNCTIONAL_VERIFICATION.json",
        "docs/NV1972_DME_FUNCTIONAL_20261123.md",
    ],
    "BPS_UCRL11980_1965": ["docs/EOS_ORIGINAL_PAPERS_20261123.md"],
}


def digest(path):
    return runner.digest(path)


def save(path, value):
    runner.save(path, value)


def deny_native_call(action, *args):
    """A call here is a census bug: record it privately and prevent execution."""
    events = Path(os.environ["NS_COWLING_SOURCE_GATE_GUARD_DIR"])
    events.mkdir(parents=True, exist_ok=True)
    event = events / f"attempt_{os.getpid()}_{time.time_ns()}.json"
    save(event, {"action": str(action), "native_executed": False})
    raise RuntimeError("Native matching is forbidden in this unavailable-source census")


def init_source_gate_worker(settings):
    ORIGINAL_INIT_WORKER(settings)
    runner.invoke = deny_native_call


def public_name(value):
    if not isinstance(value, str) or not value or len(value) > 200:
        raise ValueError("Invalid public model name")
    if any(token in value for token in ("/", "\\", "://", "\n", "\r", "\x00")):
        raise ValueError("Model name contains path or control characters")
    return value


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--original-catalog", type=Path, required=True)
    ap.add_argument("--native-exe", type=Path, required=True)
    ap.add_argument("--runtime-root", type=Path, required=True)
    ap.add_argument("--public-out", type=Path,
                    default=ROOT / "data/eos_nv_dme_20261123/remaining_census")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    if args.workers != 8:
        ap.error("This bounded census records exactly eight source-gate workers")
    runtime = args.runtime_root.resolve()
    public = args.public_out.resolve()
    if runtime.is_relative_to(ROOT) or runtime.is_relative_to(public) or public.is_relative_to(runtime):
        raise ValueError("Full runtime metadata must remain outside repository and public tree")
    expected_public = (ROOT / "data/eos_nv_dme_20261123/remaining_census").resolve()
    if public != expected_public:
        raise ValueError("This task owns only the declared new remaining_census directory")
    if runtime.exists() or public.exists():
        raise FileExistsError("Use fresh census runtime and public directories; existing results are immutable")
    original_sha = digest(args.original_catalog)
    binary_sha = digest(args.native_exe)
    script_sha = digest(__file__)
    runner_sha = digest(runner.__file__)
    if binary_sha != EXPECTED_MATCHER_SHA256:
        raise ValueError("Pinned numerical matcher identity differs")
    original = json.loads(args.original_catalog.read_text())
    cores = sorted((copy.deepcopy(s) for s in original["sources"] if s["role"] == "core"),
                   key=lambda s:s["source_id"])
    if len(cores) != 32 or len({s["source_id"] for s in cores}) != 32:
        raise ValueError("Expected the unchanged 32-core source set")
    # Different source representations may share a physical model name.
    # Opaque public labels preserve all 32 identities without disclosing paths.
    core_public_ids = {s["source_id"]:f"core_{i:03d}" for i,s in enumerate(cores,1)}
    availability = []
    for source in cores:
        if source.get("available") is not True:
            raise ValueError("An inherited core is unavailable before this census")
        if set(source.get("full_selected_sha256", {})) != {"eos_data.csv", "eos_data.nml"}:
            raise ValueError("Core lacks complete selected-source byte identities")
        if not source.get("matching_views"):
            raise ValueError("Core has no prescribed matching views")
        availability.append({
            "core_source_id": core_public_ids[source["source_id"]],
            "core_model": public_name(source["model"]),
            "private": bool(source.get("private")),
            "available": True,
            "declared_source_files_present_and_hash_verified": True,
            "declared_views_validated": True,
        })
    # These are explicit unavailable records, not fabricated numerical sources.
    crusts = [{
        "source_id": "crust:" + model, "role": "crust", "model": model,
        "available": False, "status": state, "private": False,
        "matching_views": [], "source_proof": SOURCE_PROOFS[model],
        "ground_state_eos_admitted": False,
    } for model, state in CRUST_STATES.items()]
    for model, paths in SOURCE_PROOFS.items():
        for path in paths:
            if not (ROOT/path).is_file():
                raise FileNotFoundError("Required source proof missing for " + model)
    for crust in crusts:
        if any(k in crust for k in ("internal_path", "compose_path", "full_selected_rows")):
            raise AssertionError("Unavailable crust must not contain a fabricated selected EOS")
    sources = cores + crusts
    runner.validate_catalog(sources)
    cells = runner.enumerate_cells(sources)
    if len(cells) != 672 or len(runner.METHODS) != 7:
        raise ValueError("Cartesian census dimensions differ from 32 x 3 x 7")
    runtime.mkdir(parents=True)
    catalog_path = runtime / "SOURCE_CATALOG_FULL.json"
    catalog = {
        "schema": "ns_cowling.nv_remaining_unavailable_source_catalog.v1",
        "source_catalog_sha256": original_sha,
        "native_executable_sha256": binary_sha,
        "methods": list(runner.METHODS), "sources": sources,
        "counts": {"core": 32, "unavailable_crust": 3, "cells": 672},
        "purpose": "admission census only; unavailable inputs never reach matcher",
    }
    save(catalog_path, catalog)
    save(runtime/"PUBLIC_SOURCE_ID_MAP.json", core_public_ids)
    campaign = runtime / "campaign"
    private_artifacts = runtime / "private_artifacts"
    guard_events = runtime / "native_call_guard_events"
    guard_events.mkdir()
    os.environ["NS_COWLING_SOURCE_GATE_GUARD_DIR"] = str(guard_events)
    source_proof_hashes = {
        path: digest(ROOT/path) for paths in SOURCE_PROOFS.values() for path in paths
    }
    save(runtime / "PREPARATION.json", {
        "original_catalog_path": str(args.original_catalog.resolve()),
        "original_catalog_sha256": original_sha, "cores_available_and_validated": 32,
        "matcher_path": str(args.native_exe.resolve()), "matcher_sha256": binary_sha,
        "runner_sha256": runner_sha, "preparation_script_sha256": script_sha,
        "source_proof_hashes": source_proof_hashes,
        "raw_catalog_and_numerical_paths_public": False,
    })
    saved_argv = sys.argv
    saved_init, saved_invoke = runner.init_worker, runner.invoke
    command = [str(Path(runner.__file__).resolve()), "--catalog", str(catalog_path),
               "--native-exe", str(args.native_exe.resolve()), "--out", str(campaign),
               "--private-artifacts", str(private_artifacts), "--work", str(runtime/"work"),
               "--workers", str(args.workers)]
    save(runtime / "RUNNER_COMMAND.json", command)
    try:
        runner.init_worker = init_source_gate_worker
        runner.invoke = deny_native_call
        sys.argv = command
        # Runner stdout includes private settings, so keep it in the runtime tree.
        with (runtime/"runner_stdout.log").open("w") as stream:
            old_stdout = sys.stdout
            try:
                sys.stdout = stream
                runner.main()
            finally:
                sys.stdout = old_stdout
    finally:
        runner.init_worker, runner.invoke = saved_init, saved_invoke
        sys.argv = saved_argv
    raw_summary = json.loads((campaign/"SUMMARY.json").read_text())
    raw_results = [json.loads(line) for line in (campaign/"RESULTS.jsonl").read_text().splitlines() if line]
    if len(raw_results) != 672 or len({r["cell_id"] for r in raw_results}) != 672:
        raise AssertionError("Missing or duplicate source-gate records")
    if any(r["status"] != "BLOCKED_SOURCE" or r["native_candidate_attempts"] != 0 for r in raw_results):
        raise AssertionError("An unavailable source escaped the input gate")
    if list(guard_events.iterdir()):
        raise AssertionError("Native invocation was attempted; no public success projection may be made")
    if raw_summary["status"] != "COMPLETED" or not raw_summary["input_and_binary_hashes_unchanged"]:
        raise AssertionError("Runtime input identity verification failed")
    if digest(args.original_catalog) != original_sha or digest(__file__) != script_sha:
        raise AssertionError("Preparation inputs changed during source-gate census")
    if digest(runner.__file__) != runner_sha or digest(args.native_exe) != binary_sha:
        raise AssertionError("Runner or numerical binary changed during census")
    projected = []
    for row in sorted(raw_results, key=lambda r: r["cell_id"]):
        model = public_name(row["crust_model"])
        if model not in CRUST_STATES:
            raise AssertionError("Unexpected crust in result")
        expected_states = {"crust:" + model: CRUST_STATES[model]}
        if row.get("source_states") != expected_states:
            raise AssertionError("Unexpected unavailable source state")
        projected.append({
            "cell_id": row["cell_id"], "core_source_id": core_public_ids[row["core_source_id"]],
            "core_model": public_name(row["core_model"]),
            "crust_model": model, "method": row["method"], "private": bool(row["private"]),
            "status": "BLOCKED_SOURCE", "source_status": CRUST_STATES[model],
            "native_candidate_attempts": 0, "native_matching_executed": False,
            "historical_output_reused": False, "numerical_matching_result": False,
        })
    expected_triples = {(core_public_ids[c["source_id"]], s["model"], m) for c in cores for s in crusts for m in runner.METHODS}
    actual_triples = {(r["core_source_id"], r["crust_model"], r["method"]) for r in projected}
    if actual_triples != expected_triples:
        raise AssertionError("Public projection lost Cartesian coverage")
    projection_keys = set(projected[0])
    if any(set(r) != projection_keys for r in projected):
        raise AssertionError("Projection contains non-allowlisted fields")
    for row in projected:
        for value in row.values():
            if isinstance(value, str) and any(s in value for s in ("/tmp/", "/home/", "://", "internal_path")):
                raise AssertionError("Private metadata leaked into public cell projection")
    public.mkdir(parents=True)
    save(public/"ALL_CELLS.json", projected)
    with (public/"ALL_CELLS.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(projected[0]))
        writer.writeheader();writer.writerows(projected)
    save(public/"CORE_AVAILABILITY.json", sorted(availability, key=lambda x:x["core_source_id"]))
    public_summary = {
        "schema": "ns_cowling.nv_remaining_source_admission_census.v1",
        "status": "SOURCE_ADMISSION_CENSUS_COMPLETED",
        "scope": "quarantine and missing-input closure; no numerical matching was run",
        "core_sources": 32, "unique_core_model_names":len({s["model"] for s in cores}),
        "core_source_ids_are_opaque_public_labels":True,
        "unavailable_crust_sources": 3, "matching_methods": 7,
        "source_pairs": 96, "requested_cells": 672, "recorded_cells": 672,
        "statuses": dict(Counter(r["status"] for r in projected)),
        "source_status_counts": dict(Counter(r["source_status"] for r in projected)),
        "by_method": {m:{"records":96,"blocked_source":96,"native_candidate_attempts":0} for m in runner.METHODS},
        "native_candidate_attempts": 0, "native_matcher_invocations": 0,
        "native_invocation_guard_events": 0, "native_invocation_guard_installed_in_each_worker": True,
        "numerical_matching_runs": 0, "surviving_new_matched_eos": 0,
        "source_gate_worker_processes": 8, "core_source_files_validated": True,
        "input_and_binary_hashes_unchanged": True,
        "native_binary_sha256": binary_sha, "existing_runner_sha256": runner_sha,
        "preparation_script_sha256": script_sha, "original_catalog_sha256": original_sha,
        "runtime_catalog_sha256": digest(catalog_path),
        "public_projection_fields": list(projected[0]),
        "private_numerical_paths_or_endpoints_published": False,
        "historical_1344_results_modified": False, "historical_650_public_results_modified": False,
        "scientific_no_go_claim": False,
        "negative_result_scope": "The three crust source contracts are currently unavailable; these records contain no physical matching verdict.",
        "source_proof_files": SOURCE_PROOFS,
        "diagnostic_hf_receipt_binding": "Later fixed-cell receipts may establish numerical behavior; they do not admit a ground-state EOS without its own source contract.",
    }
    save(public/"SUMMARY.json", public_summary)
    save(public/"PUBLICATION_MANIFEST.json", {
        "schema":"ns_cowling.value_free_census_projection.v1",
        "projection_audited":True,
        "numerical_source_rows_or_endpoints_in_projection":False,
        "runtime_settings_and_input_identities_kept_private":True,
        "files":{name:digest(public/name) for name in ("ALL_CELLS.json","ALL_CELLS.csv","CORE_AVAILABILITY.json","SUMMARY.json")},
    })
    print(json.dumps({"status":public_summary["status"],"records":672,"native_candidate_attempts":0,
                      "native_matcher_invocations":0,"workers":8,"source_status_counts":public_summary["source_status_counts"],
                      "public_manifest_sha256":digest(public/"PUBLICATION_MANIFEST.json")},sort_keys=True))


if __name__ == "__main__":
    main()
