#!/usr/bin/env python3
"""Fresh Cartesian EOS campaign: Python schedules; Fortran computes physics.

One deterministic certified output is retained for each source-pair/method.
Every requested triple has an outcome, including unavailable sources. Inputs
are prescribed source phase views, never intervals cut at failed physics
checks. Historical matches and Python root/sign prefilters are not used.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import traceback
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

# Set before importing frontend helpers; each native child uses one BLAS thread.
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT / "python_frontend")]
from match_core_crust import (
    verify_internal_readback, verify_internal_metadata_readback,
    _consistency_report_errors,
)
from match_report_admission import admit_match_report, classify_failure

METHODS = (
    "first_law_constrained_hermite_pressure_bridge",
    "cubic_F", "quintic_F", "cubic_P", "quintic_P", "css", "maxwell",
)
SCHEMA = "ns_cowling.fresh_cartesian_native_campaign.v1"
SETTINGS = {}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n")
    temporary.replace(path)


def resolved(path):
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def enumerate_cells(sources):
    crusts = sorted((s for s in sources if s["role"] == "crust"),
                    key=lambda s: s["source_id"])
    cores = sorted((s for s in sources if s["role"] == "core"),
                   key=lambda s: s["source_id"])
    ids = [s["source_id"] for s in sources]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate source identity")
    return [
        (f"cell_{i:06d}", crust, core, method)
        for i, (crust, core, method) in enumerate(
            itertools.product(crusts, cores, METHODS), start=1)
    ]


def source_available(source):
    return source.get("available") is True


def allowed_views(source, method):
    return [
        v for v in source.get("matching_views", [])
        if method in v.get("methods", []) and int(v.get("rows", 0)) >= 3
    ]


def validate_catalog(sources):
    for source in sources:
        if not source_available(source):
            continue
        for name, expected in source.get("full_selected_sha256", {}).items():
            if digest(resolved(source["internal_path"]) / name) != expected:
                raise ValueError("Prepared source byte identity mismatch: " + source["source_id"])
        n = source.get("full_selected_rows")
        if type(n) is not int or n < 3:
            raise ValueError("invalid selected source row count: " + source["source_id"])
        for view in source.get("matching_views", []):
            a, b, rows = view.get("start"), view.get("stop"), view.get("rows")
            if not all(type(x) is int for x in (a, b, rows)):
                raise ValueError("noninteger source view bounds")
            if not (0 <= a < b <= n and rows == b-a):
                raise ValueError("source view escapes declared selected source")
            if not set(view.get("methods", [])) <= set(METHODS):
                raise ValueError("unknown source view method permission")


def input_identities(sources, exe, catalog_path):
    files = {str(Path(exe).resolve()), str(Path(catalog_path).resolve()),
             str(Path(__file__).resolve()),
             str(ROOT / "scripts/match_core_crust.py"),
             str(ROOT / "scripts/match_report_admission.py")}
    for source in sources:
        paths = [source.get("internal_path")]
        paths += [v.get("internal_path") for v in source.get("matching_views", [])]
        for raw in paths:
            if not raw:
                continue
            directory = resolved(raw)
            for name in ("eos_data.csv", "eos_data.nml"):
                p = directory / name
                if not p.is_file():
                    raise FileNotFoundError(p)
                files.add(str(p))
    return {p: digest(p) for p in sorted(files)}


def validate_private_location(private, public):
    private, public = Path(private).resolve(), Path(public).resolve()
    if (private.is_relative_to(ROOT) or private.is_relative_to(public)
            or public.is_relative_to(private)):
        raise ValueError("Private artifacts must be outside the repository and public output tree")


def init_worker(settings):
    global SETTINGS
    SETTINGS = settings


def invoke(action, *args):
    command = [SETTINGS["exe"], action, *(str(v) for v in args)]
    start = time.perf_counter()
    try:
        result = subprocess.run(command, capture_output=True, text=True,
                                timeout=SETTINGS["timeout"])
    except subprocess.TimeoutExpired as exc:
        return {
            "ok": False, "action": action, "elapsed_seconds": time.perf_counter()-start,
            "returncode": None, "payload": {},
            "failure_category": "runtime_timeout",
            "message": f"Native operation exceeded {SETTINGS['timeout']} s",
        }
    payload = None
    for line in reversed(result.stdout.splitlines()):
        try:
            item = json.loads(line)
        except (ValueError, TypeError):
            continue
        if isinstance(item, dict):
            payload = item
            break
    valid_payload = isinstance(payload, dict)
    okay = (result.returncode == 0 and valid_payload
            and payload.get("status") == "ok"
            and payload.get("operation") == action)
    message = (payload or {}).get("message")
    if not okay and not message:
        message = (result.stderr.strip() or result.stdout.strip()
                   or "Missing or invalid structured native receipt")
    return {
        "ok": okay, "action": action, "returncode": result.returncode,
        "elapsed_seconds": time.perf_counter()-start, "payload": payload or {},
        "message": message, "stderr": result.stderr[-4000:],
        "failure_category": (
            "runtime_signal" if result.returncode < 0 else
            "protocol_failure" if not valid_payload else None),
    }


def failure(stage, call=None, message=None, category=None, details=None):
    call = call or {}
    payload = call.get("payload", {})
    message = message or call.get("message") or "Unspecified failed gate"
    category = category or call.get("failure_category")
    if category is None:
        if stage == "matching" and payload.get("window_search_exhausted") is True:
            category = "no_certified_bridge_on_supplied_endpoint_grid"
        else:
            category = classify_failure(stage, message, payload.get("category"))
    return {
        "status": "REJECTED", "stage": stage, "failure_category": category,
        "message": message, "native_call": call, "details": details,
    }


def check_phase_transport(report, matched_csv, crust_csv, core_csv):
    """Categorical row/phase transport only; no thermodynamic reconstruction."""
    def labels(path):
        with Path(path).open(newline="") as f:
            return [int(row["Iphase"]) for row in csv.DictReader(f)]
    output, low, high = map(labels, (matched_csv, crust_csv, core_csv))
    left = int(report["i1"]) + 1
    right = int(report["i2"])
    high_start = len(output) - (len(high)-right)
    errors = []
    if output != report.get("output_phase_ids"):
        errors.append("native phase receipt differs from serialized output")
    if output[:left] != low[:left] or output[high_start:] != high[right:]:
        errors.append("retained source phase labels differ")
    if any(v != report["bridge_phase_id"] for v in output[left:high_start]):
        errors.append("bridge phase label differs")
    return {"status": "PASS" if not errors else "FAIL", "errors": errors}


def interface_admitted(receipt, method):
    checks = receipt.get("checks")
    return (
        receipt.get("schema") == "ns-cowling.match-report-admission.v1"
        and receipt.get("method") == method
        and receipt.get("status") == "ADMITTED"
        and receipt.get("failed_checks") == []
        and isinstance(checks, dict) and bool(checks)
        and all(isinstance(item, dict) and item.get("pass") is True
                for item in checks.values())
    )


def try_candidate(directory, crust_view, core_view, method):
    low = resolved(crust_view["internal_path"])
    high = resolved(core_view["internal_path"])
    internal = directory / "internal"
    compose = directory / "compose"
    readback = directory / "readback"
    for p in (internal, compose, readback):
        p.mkdir(parents=True)
    action = "maxwell-transition" if method == "maxwell" else "match-core"
    args = ([low, high, internal, "true"] if method == "maxwell" else
            [low, high, internal, method, SETTINGS["n_bridge"], "true"])
    call = invoke(action, *args)
    if not call["ok"]:
        return failure("matching", call)
    report = call["payload"]
    save(directory / "NATIVE_MATCH.json", report)
    if type(report.get("n")) is not int or report["n"] < 3:
        return failure("protocol", call, "Native output row count missing/invalid",
                       "protocol_failure")
    n = report["n"]
    if method != "maxwell":
        if report.get("method") != method:
            return failure("protocol", call, "Native method identity mismatch",
                           "protocol_failure")
        phase = check_phase_transport(
            report, internal / "eos_data.csv", low / "eos_data.csv",
            high / "eos_data.csv")
        save(directory / "PHASE_TRANSPORT.json", phase)
        if phase["status"] != "PASS":
            return failure("admission", call, "Phase transport failed",
                           details=phase)
        first_law = report.get("first_law_relative_residual")
        if not isinstance(first_law, (int, float)) or not 0 <= first_law <= 5e-9:
            return failure("admission", call, "First-law endpoint residual failed",
                           details={"residual": first_law, "limit": 5e-9})
        if report.get("dm_mass_gauge_meV") != 0.0:
            return failure("admission", call, "Unexpected physical energy shift")
    consistency_path = directory / "CONSISTENCY_INTERNAL.json"
    check = invoke("verify-csv", internal, consistency_path)
    if not check["ok"]:
        return failure("consistency", check)
    consistency = json.loads(consistency_path.read_text())
    errors = _consistency_report_errors(consistency, n)
    if check["payload"].get("overall") != "PASS" or errors:
        return failure("consistency", check,
                       "Assembled native consistency failed: "
                       + str(consistency.get("failed_checks", "")),
                       details={"errors": errors, "receipt": consistency})

    write = invoke("write", internal, compose)
    if not write["ok"]:
        return failure("write", write)
    read = invoke("read", compose, readback)
    if not read["ok"]:
        return failure("readback", read)
    rows = verify_internal_readback(
        internal / "eos_data.csv", readback / "eos_data.csv",
        expected_n=n, actual_n=read["payload"].get("n", -1))
    metadata = verify_internal_metadata_readback(
        internal / "eos_data.nml", readback / "eos_data.nml")
    save(directory / "READBACK_ROWS.json", rows)
    save(directory / "READBACK_METADATA.json", metadata)
    if rows["status"] != "PASS" or metadata["status"] != "PASS":
        return failure("readback", read, "Serialized rows/metadata changed",
                       details={"rows": rows, "metadata": metadata})
    written_check = invoke("verify", compose, directory / "CONSISTENCY_COMPOSE.json")
    if not written_check["ok"]:
        return failure("consistency", written_check)
    written_report = json.loads((directory / "CONSISTENCY_COMPOSE.json").read_text())
    errors = _consistency_report_errors(written_report, n)
    if written_check["payload"].get("overall") != "PASS" or errors:
        return failure("consistency", written_check,
                       "Written CompOSE consistency failed",
                       details={"errors": errors, "receipt": written_report})
    if method in ("css", "maxwell"):
        admission = admit_match_report(
            method, report, readback / "eos_data.csv", readback / "eos_data.nml")
        save(directory / "INTERFACE_ADMISSION.json", admission)
        if not interface_admitted(admission, method):
            return failure("admission", message="Serialized interface admission failed",
                           details=admission)
    else:
        admission = {
            "status": "PASS", "method": method,
            "continuous_bridge_certificate": "native analytic interval gate",
            "first_law_limit": 5e-9,
            "serialized_phase_and_thermodynamic_identity": "verified",
            "scope": "bridge certificate plus finite-table source verifier; no global-minimum claim",
        }
        save(directory / "INTERFACE_ADMISSION.json", admission)
    return {
        "status": "SURVIVED", "native_report": report,
        "consistency": written_report, "native_seconds": call["elapsed_seconds"],
        "verification": {
            "native_internal": "PASS", "compose_readback_all_columns": "PASS",
            "compose_readback_all_metadata": "PASS", "native_compose": "PASS",
            "interface_admission": "PASS",
        },
    }


def archive_candidate(directory, path):
    """Lossless deterministic archive; full readback duplicate is omitted."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as f:
        with gzip.GzipFile(fileobj=f, mode="wb", mtime=0, filename="") as z:
            with tarfile.open(fileobj=z, mode="w") as tar:
                for p in sorted(directory.rglob("*")):
                    if not p.is_file() or "readback" in p.relative_to(directory).parts:
                        continue
                    info = tar.gettarinfo(str(p), arcname=str(p.relative_to(directory)))
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mtime = 0
                    with p.open("rb") as source:
                        tar.addfile(info, source)
    return {"path": str(path), "sha256": digest(path), "bytes": path.stat().st_size}


def compact_native(report):
    omit = {"output_phase_ids"}
    return {k: v for k, v in report.items() if k not in omit}


def run_cell(cell):
    cell_id, crust, core, method = cell
    start = time.perf_counter()
    base = {
        "cell_id": cell_id, "crust_source_id": crust["source_id"],
        "core_source_id": core["source_id"], "method": method,
        "crust_model": crust.get("model"), "core_model": core.get("model"),
        "private": bool(crust.get("private") or core.get("private")),
        "status": "RUNNING", "native_candidate_attempts": 0,
        "historical_output_reused": False, "python_physical_prefilter": False,
    }
    if not source_available(crust) or not source_available(core):
        base.update(status="BLOCKED_SOURCE", failure_category="data_provenance",
                    source_states={s["source_id"]: s.get("status", s.get("state"))
                                   for s in (crust, core) if not source_available(s)})
        base["elapsed_seconds"] = time.perf_counter()-start
        return base
    low_views, high_views = allowed_views(crust, method), allowed_views(core, method)
    if not low_views or not high_views:
        base.update(status="BLOCKED_SOURCE_SUPPORT", failure_category="input_contract",
                    message="No prescribed source phase view with >=3 knots for this method")
        base["elapsed_seconds"] = time.perf_counter()-start
        return base
    trials = []
    try:
        with tempfile.TemporaryDirectory(prefix=cell_id+"_", dir=SETTINGS["work"]) as temp:
            for index, (low, high) in enumerate(itertools.product(low_views, high_views)):
                candidate = Path(temp) / f"candidate_{index:04d}"
                try:
                    trial = try_candidate(candidate, low, high, method)
                except Exception as exc:
                    trial = failure("runtime", message=f"{type(exc).__name__}: {exc}",
                                    category="python_orchestration_failure",
                                    details={"traceback": traceback.format_exc()})
                receipt_path = candidate / "NATIVE_MATCH.json"
                if trial["status"] != "SURVIVED" and receipt_path.is_file():
                    first_match = json.loads(receipt_path.read_text())
                    trial["initial_matching_receipt"] = compact_native(first_match)
                    if method != "maxwell":
                        trial["underlying_failure_category"] = trial["failure_category"]
                        trial["failure_category"] = "post_bridge_gate_failure_search_incomplete"
                        trial["search_scope"] = "later endpoint pairs in this phase view remain untested"
                trial["crust_view_id"] = low["id"]
                trial["core_view_id"] = high["id"]
                base["native_candidate_attempts"] += 1
                if trial["status"] == "SURVIVED":
                    # Ordinary bridge trimming at the interface is intended.
                    # The outer crust origin and highest selected core row must
                    # be present to call the retained endpoint chain complete.
                    crust_rows = int(crust["full_selected_rows"])
                    core_rows = int(core["full_selected_rows"])
                    low_start = int(low.get("start", low.get("start_zero_based", -1)))
                    high_stop = int(high.get("stop", high.get("stop_exclusive", -1)))
                    covers_endpoints = low_start == 0 and high_stop == core_rows
                    coverage = ("SELECTED_SOURCE_ENDPOINT_CHAIN" if covers_endpoints
                                else "FINITE_PHASE_COMPONENT_ONLY")
                    provenance = {
                        **base, "status": "SURVIVED", "source_coverage": coverage,
                        "crust_view": low, "core_view": high,
                        "source_catalog_sha256": SETTINGS["catalog_sha256"],
                        "native_binary_sha256": SETTINGS["binary_sha256"],
                        "n_bridge": SETTINGS["n_bridge"],
                        "crust_selected_rows": crust_rows, "core_selected_rows": core_rows,
                        "prior_failed_phase_pairs": len(trials),
                        "surface_to_star_admitted": False,
                        "global_phase_ground_state_proven": False,
                        "local_charge_neutrality": "assumed_from_cold_beta_source_contract",
                    }
                    save(candidate / "PROVENANCE.json", provenance)
                    archive_root = (SETTINGS["private_artifacts"] if base["private"]
                                    else str(Path(SETTINGS["out"]) / "artifacts"))
                    artifact = archive_candidate(candidate, Path(archive_root)/(cell_id+".tar.gz"))
                    base.update(
                        status="SURVIVED", source_coverage=coverage,
                        crust_view_id=low["id"], core_view_id=high["id"],
                        verification=trial["verification"],
                        native_report=compact_native(trial["native_report"]),
                        artifact=artifact,
                        source_view_bounds={"crust": [low_start, low.get("stop")],
                                            "core": [high.get("start"), high_stop]},
                    )
                    break
                trials.append(trial)
                # Failed output rows are not candidate EOS deliverables.
                if candidate.exists():
                    shutil.rmtree(candidate)
            else:
                base["status"] = "NO_ADMITTED_MATCH"
                base["failure_categories"] = dict(Counter(
                    t["failure_category"] for t in trials))
            if trials:
                trace_root = (SETTINGS["private_artifacts"] if base["private"] else SETTINGS["out"])
                trace_path = Path(trace_root) / "failure_traces" / (cell_id+".json.gz")
                trace_path.parent.mkdir(parents=True, exist_ok=True)
                with trace_path.open("wb") as raw:
                    with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as z:
                        z.write((json.dumps(trials, ensure_ascii=False, allow_nan=False)+"\n").encode())
                base["failure_trace"] = {"path": str(trace_path), "sha256": digest(trace_path)}
                base["prior_or_total_failure_categories"] = dict(Counter(
                    t["failure_category"] for t in trials))
    except Exception as exc:
        base.update(status="RUNTIME_FAILURE", failure_category="python_orchestration_failure",
                    message=f"{type(exc).__name__}: {exc}", traceback=traceback.format_exc())
    base["elapsed_seconds"] = time.perf_counter()-start
    return base


def summarize(cells, results, elapsed, settings, unchanged):
    by_method = {}
    for method in METHODS:
        group = [r for r in results if r["method"] == method]
        by_method[method] = {
            "requested": len(group),
            "statuses": dict(Counter(r["status"] for r in group)),
            "survivor_coverage": dict(Counter(
                r["source_coverage"] for r in group if r["status"] == "SURVIVED")),
            "native_phase_pair_attempts": sum(r["native_candidate_attempts"] for r in group),
            "elapsed_cell_seconds": sum(r["elapsed_seconds"] for r in group),
        }
    return {
        "schema": SCHEMA, "status": "COMPLETED" if unchanged else "INPUT_IDENTITY_CHANGED",
        "requested_cells": len(cells), "recorded_cells": len(results),
        "source_pairs": len({(c[1]["source_id"], c[2]["source_id"]) for c in cells}),
        "statuses": dict(Counter(r["status"] for r in results)),
        "survivors": sum(r["status"] == "SURVIVED" for r in results),
        "by_method": by_method, "wall_seconds": elapsed,
        "input_and_binary_hashes_unchanged": unchanged, "settings": settings,
        "selection_rule": "first independently admitted output per prescribed source-pair/method",
        "coverage_limit": "finite supplied source domains; surface/TOV/global-ground-state not certified",
        "negative_result_scope": (
            "native failures on every prescribed phase-pair tested; a later "
            "assembled-table gate failure is not an impossibility theorem "
            "for every continuum interpolation or alternate representation"),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--catalog", required=True, type=Path)
    p.add_argument("--native-exe", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--private-artifacts", required=True, type=Path)
    p.add_argument("--work", type=Path, default=Path("/tmp/ns_eos_reaudit_20261123/campaign_work"))
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--n-bridge", type=int, default=80)
    p.add_argument("--timeout", type=float, default=180.0)
    p.add_argument("--pilot-pairs", type=int, default=0,
                   help="Stratified pair sample, all seven methods; never labelled full census.")
    p.add_argument("--resume", action="store_true")
    args = p.parse_args()
    if args.workers < 1 or args.n_bridge < 8 or args.timeout <= 0:
        p.error("positive workers/timeout and n-bridge>=8 required")
    catalog = json.loads(args.catalog.read_text())
    sources = catalog["sources"]
    validate_catalog(sources)
    cells = enumerate_cells(sources)
    full_requested = len(cells)
    if args.pilot_pairs:
        pairs = [cells[i:i+len(METHODS)] for i in range(0, len(cells), len(METHODS))]
        count = min(args.pilot_pairs, len(pairs))
        indices = sorted(set(round(i*(len(pairs)-1)/max(1, count-1)) for i in range(count)))
        cells = [cell for i in indices for cell in pairs[i]]
    validate_private_location(args.private_artifacts, args.out)
    args.out.mkdir(parents=True, exist_ok=True)
    args.work.mkdir(parents=True, exist_ok=True)
    args.private_artifacts.mkdir(parents=True, exist_ok=True)
    identities = input_identities(sources, args.native_exe, args.catalog)
    settings = {
        "exe": str(args.native_exe.resolve()), "out": str(args.out.resolve()),
        "private_artifacts": str(args.private_artifacts.resolve()),
        "work": str(args.work.resolve()), "n_bridge": args.n_bridge,
        "timeout": args.timeout, "workers": args.workers,
        "catalog_sha256": digest(args.catalog), "binary_sha256": digest(args.native_exe),
        "runner_sha256": digest(__file__), "full_requested_cells": full_requested,
        "pilot": bool(args.pilot_pairs), "all_methods": list(METHODS),
        "input_hash_manifest": "INPUT_IDENTITIES.json",
    }
    previous = []
    journal = args.out / "RESULTS.jsonl"
    if journal.exists():
        if not args.resume:
            raise FileExistsError("Output contains a campaign; use an empty directory or explicit --resume")
        if json.loads((args.out/"SETTINGS.json").read_text()) != settings:
            raise ValueError("Resume settings/source/binary/runner identities differ")
        if json.loads((args.out/"INPUT_IDENTITIES.json").read_text()) != identities:
            raise ValueError("Resume source byte identities differ")
        previous = [json.loads(line) for line in journal.read_text().splitlines() if line.strip()]
        if len({r["cell_id"] for r in previous}) != len(previous):
            raise ValueError("Duplicate cells in resume journal")
        for result in previous:
            for key in ("artifact", "failure_trace"):
                if key in result and digest(result[key]["path"]) != result[key]["sha256"]:
                    raise ValueError("Resume artifact identity mismatch")
    else:
        save(args.out/"SETTINGS.json", settings)
        save(args.out/"INPUT_IDENTITIES.json", identities)
    done = {r["cell_id"] for r in previous}
    expected = {c[0] for c in cells}
    if not done <= expected:
        raise ValueError("Resume journal contains out-of-scope cells")
    pending = [c for c in cells if c[0] not in done]
    start = time.perf_counter()
    results = list(previous)
    print(json.dumps({"event": "START", "requested": len(cells),
                      "pending": len(pending), "workers": args.workers,
                      "pilot": settings["pilot"]}), flush=True)
    with journal.open("a") as stream:
        with ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker,
                                 initargs=(settings,)) as pool:
            futures = {pool.submit(run_cell, c): c[0] for c in pending}
            for future in as_completed(futures):
                result = future.result()
                stream.write(json.dumps(result, ensure_ascii=False, allow_nan=False)+"\n")
                stream.flush()
                results.append(result)
                if len(results) % 100 == 0 or len(results) == len(cells):
                    print(json.dumps({
                        "event": "PROGRESS", "completed": len(results), "total": len(cells),
                        "survivors": sum(r["status"] == "SURVIVED" for r in results),
                        "elapsed_seconds": round(time.perf_counter()-start, 3),
                    }), flush=True)
    if {r["cell_id"] for r in results} != expected or len(results) != len(expected):
        raise RuntimeError("Cartesian accounting mismatch")
    unchanged = input_identities(sources, args.native_exe, args.catalog) == identities
    summary = summarize(cells, results, time.perf_counter()-start, settings, unchanged)
    results.sort(key=lambda r: r["cell_id"])
    survivors = [r for r in results if r["status"] == "SURVIVED"]
    if not unchanged:
        summary["candidate_survivors_quarantined"] = len(survivors)
        summary["survivors"] = 0
        save(args.out/"SUMMARY.json", summary)
        save(args.out/"QUARANTINED_RESULTS.json", results)
        raise RuntimeError("Input or executable changed; all outputs quarantined")
    save(args.out/"SUMMARY.json", summary)
    save(args.out/"SURVIVORS.json", survivors)
    fields = ["cell_id", "crust_source_id", "core_source_id", "method", "status",
              "source_coverage", "crust_view_id", "core_view_id", "private",
              "native_candidate_attempts", "elapsed_seconds"]
    for name, subset in (("ALL_CELLS.csv", results), ("SURVIVORS.csv", survivors)):
        with (args.out/name).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(subset)
    print(json.dumps({"event": "COMPLETE", **summary}, ensure_ascii=False), flush=True)
    if not unchanged:
        raise RuntimeError("Input or executable changed during campaign")


if __name__ == "__main__":
    main()
