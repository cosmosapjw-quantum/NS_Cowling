#!/usr/bin/env python3
"""Execute a new, bounded campaign of sealed native-generated EOS pairs.

This runner does not modify the historical one-shot intake or its CLAIM. The
default cohort has 30 method cells and reuses its one existing admitted CSS
artifact after source/output identity and fresh serialized verification.
Physics remains in the native backend; smooth construction and independent
CSS/Maxwell admission retain different evidence statuses.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "python_frontend"), str(ROOT / "scripts")]
import compose_eos.api as eos_api
import match_core_crust as pipeline
from compose_eos import read_compose, write_compose
from compose_eos.native import NativeExecutionError, executable, run_native
from match_report_admission import admit_match_report, classify_failure

CRUSTS = tuple(f"SPG_M{i}" for i in range(1, 6))
CORES = ("GM_GM1", "RG_SLY4")
METHODS = ("quintic_P", "css", "maxwell")
SPG_ROOT = "data/research_matching_20261117/spg/native_129_repaired"
CORE_ROOT = "data/research_matching_20261117/integration/real_pairs_v2/core"
REUSE_EXECUTION = "data/current_generator_intake_20261117/runtime/run_001/MATCH_EXECUTION.json"
REUSE_CELL = "SPG_M1__GM_GM1__css"
REQUIRED = {"eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo", "SUPPORT.json"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def _safe(root: Path, relative: str) -> Path:
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts or not relative:
        raise ValueError(f"expected repository-relative path: {relative!r}")
    resolved = (root / path).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError(f"path escapes repository: {relative!r}")
    return resolved


def _json(path: Path) -> dict:
    record = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(record, dict):
        raise ValueError(f"JSON object required: {path}")
    return record


def _segments(directory: Path) -> tuple[list[dict], list[dict]]:
    tokens = (directory / "eos.nb").read_text().split()
    first, last = int(tokens[0]), int(tokens[1])
    grid = [float(x.replace("D", "E").replace("d", "e")) for x in tokens[2:]]
    if (not grid or len(grid) != last - first + 1
            or any(not math.isfinite(x) or x <= 0 for x in grid)
            or any(b <= a for a, b in zip(grid, grid[1:]))):
        raise ValueError("invalid sealed density grid")
    phases = {}
    for line in (directory / "eos.compo").read_text().splitlines():
        if line.strip():
            row = line.split()
            index = int(row[1]) - first
            if index in phases or index < 0 or index >= len(grid):
                raise ValueError("duplicate/out-of-range phase row")
            phases[index] = int(row[3])
    if set(phases) != set(range(len(grid))):
        raise ValueError("missing sealed phase rows")
    segments = []
    for index in range(len(grid)):
        if segments and segments[-1]["iphase"] == phases[index]:
            segments[-1]["row_end"] = index
        else:
            segments.append(dict(iphase=phases[index], row_start=index, row_end=index))
    return ([dict(row_start=0, row_end=len(grid) - 1, n_b_min=grid[0], n_b_max=grid[-1])], segments)


def source_snapshot(root: Path, family: str, model: str) -> dict:
    """Validate the actual outer/native receipt, then snapshot its bound bytes."""
    if family not in ("crust", "core") or model not in (CRUSTS if family == "crust" else CORES):
        raise ValueError("source outside the explicit campaign cohort")
    parent = f"{SPG_ROOT if family == 'crust' else CORE_ROOT}/{model}"
    outer_path = f"{parent}/{'RECEIPT.json' if family == 'crust' else 'DISPATCH.json'}"
    outer = _json(_safe(root, outer_path))
    native_path = f"{parent}/generated/RESULT.json"
    native = _json(_safe(root, native_path))
    if (outer.get("status") != "REQUESTED_ROWS_NATIVE_VERIFIED" or outer.get("returncode") != 0
            or (family == "core" and outer.get("child_returncode") != 0)):
        raise ValueError("outer receipt does not certify requested rows")
    if (native.get("status") != "REQUESTED_ROWS_NATIVE_VERIFIED"
            or native.get("backend") != "fortran" or native.get("native_readback") is not True
            or native.get("model") != model or native.get("family") != ("spg" if family == "crust" else "core")
            or outer["request"]["model"] != model):
        raise ValueError("native/outer model identity or readback invalid")
    request_path = f"{parent}/request.nml"
    if sha(_safe(root, request_path)) != outer["request_sha256"]:
        raise ValueError("request hash mismatch")
    inventory = outer["output_sha256"] if family == "crust" else outer["evidence"]["output_sha256"]
    if not isinstance(inventory, dict):
        raise ValueError("outer output inventory is not an object")
    sealed_paths = {(f"{parent}/{name}" if family == "crust" else f"{parent}/generated/{name}"): digest
                    for name, digest in inventory.items()}
    for path, digest in sealed_paths.items():
        if not isinstance(digest, str) or len(digest) != 64 or sha(_safe(root, path)) != digest:
            raise ValueError(f"sealed artifact missing/changed: {path}")
    prefix = f"{parent}/generated/"
    artifacts = {name[len(prefix):]: digest for name, digest in sealed_paths.items() if name.startswith(prefix)}
    if not REQUIRED.issubset(artifacts):
        raise ValueError("required CompOSE/SUPPORT files missing from outer inventory")
    support_path = f"{parent}/generated/SUPPORT.json"
    support = _json(_safe(root, support_path))
    if any(support.get(key) != "UNAPPROVED" for key in ("cross_phase_interpolation", "cross_segment_interpolation")):
        raise ValueError("SUPPORT restrictions differ from this campaign contract")
    density, phases = _segments(_safe(root, f"{parent}/generated"))
    if native["rows"] != density[0]["row_end"] + 1:
        raise ValueError("native row count differs from sealed table")
    if len(phases) != 1:
        raise ValueError("source SUPPORT requires a single contiguous phase branch")
    receipts = {outer_path: sha(_safe(root, outer_path)), native_path: sha(_safe(root, native_path))}
    if family == "crust":
        run_path = f"{SPG_ROOT}/RUN_RECEIPT.json"
        run = _json(_safe(root, run_path))
        if run["models"].get(model) != "REQUESTED_ROWS_NATIVE_VERIFIED":
            raise ValueError("SPG run receipt does not verify this model")
        receipts[run_path] = sha(_safe(root, run_path))
        generation_identity = {"scope": "historical_RUN_RECEIPT", "binaries_sha256": {
            path: digest for path, digest in run["files"].items()
            if path.endswith(("/compose_eos_generate", "/compose_eos_fortran"))}}
        if len(generation_identity["binaries_sha256"]) != 2:
            raise ValueError("SPG historical binary inventory incomplete")
        surface = outer.get("surface_input", {})
        dependencies = {surface["path"]: surface["sha256"]} if surface.get("path") else {}
    else:
        generation_identity = {"scope": "historical_dispatch", "status": "NOT_INDEPENDENTLY_OBSERVED",
                               "command": outer["command"]}
        dependencies = {}
        for path, digest in outer["source_sha256"].items():
            for marker in ("/python_frontend/", "/fortran_eos/"):
                if marker in path:
                    dependencies[marker[1:] + path.split(marker, 1)[1]] = digest
                    break
    drift = []
    for path, expected in dependencies.items():
        actual = sha(_safe(root, path)) if _safe(root, path).is_file() else None
        if actual != expected:
            drift.append(dict(path=path, expected_sha256=expected, observed_sha256=actual))
    return dict(id=f"{family}:{model}", family=family, model=model,
                variant=outer["request"].get("variant"), request_contract=outer["request"].get("contract", {}),
                compose_path=f"{parent}/generated", request_path=request_path,
                request_sha256=outer["request_sha256"], receipt_paths_sha256=receipts,
                sealed_paths_sha256=sealed_paths, artifact_sha256=artifacts,
                support_path=support_path, support_sha256=artifacts["SUPPORT.json"],
                density_segments=density, phase_segments=phases,
                generation_executable_identity=generation_identity,
                historical_dependencies_sha256=dependencies, dependency_differences=drift,
                dependency_status="HISTORICAL_GENERATOR_DEPENDENCY_VERSION_DIFFERENCE" if drift else "MATCH")


def validate_source(root: Path, source: dict) -> None:
    if source_snapshot(root, source["family"], source["model"]) != source:
        raise ValueError("source changed after campaign snapshot")


def reuse_identity(root: Path, execution_path: Path, cell: dict, sources: dict) -> dict:
    """Check immutable old inputs/settings/output before fresh readback admission."""
    execution = _json(execution_path)
    if (cell["id"] != REUSE_CELL or execution["cell_id"] != cell["id"]
            or execution["n_bridge"] != cell["n_bridge"] or execution.get("admitted") is not True):
        raise ValueError("historical reuse cell/settings/admission mismatch")
    manifest_path = _safe(root, execution["input_manifest_path"])
    if sha(manifest_path) != execution["input_manifest_sha256"]:
        raise ValueError("historical input manifest hash mismatch")
    historical = {source["id"]: source for source in _json(manifest_path)["sources"]}
    result = execution["result"]
    if (result.get("survived") is not True or result["method"] != cell["method"]
            or result["admission"]["status"] != "ADMITTED"
            or result["match_report"]["n_bridge"] != cell["n_bridge"]):
        raise ValueError("historical result/method admission mismatch")
    for side in ("crust", "core"):
        source = sources[cell[f"{side}_source_id"]]
        old = historical[source["id"]]
        for field in ("compose_path", "request_sha256", "artifact_sha256", "support_sha256"):
            if old[field] != source[field]:
                raise ValueError(f"reuse {side} {field} identity mismatch")
        expected = {name: digest for name, digest in source["artifact_sha256"].items() if name in pipeline.REQUIRED_SOURCE_FILES or name == "eos.compo"}
        if result["source_sha256"][side] != expected:
            raise ValueError(f"reuse {side} CompOSE identity mismatch")
    directory = execution_path.parent / cell["id"]
    if not set(pipeline.REQUIRED_SOURCE_FILES).issubset(result["output_sha256"]):
        raise ValueError("historical output inventory incomplete")
    for name, digest in result["output_sha256"].items():
        if sha(_safe(directory / "compose", name)) != digest:
            raise ValueError(f"historical output hash mismatch: {name}")
    return dict(execution=execution, directory=directory)


@contextmanager
def _observe_smooth_kernel(state: dict, work: Path):
    """Observe the native boundary only; sequential runner restores it in finally."""
    original = eos_api.run_native
    def observed(action, *args):
        try:
            payload = original(action, *args)
        except NativeExecutionError as exc:
            # An operation-tagged backend failure proves process execution.
            # Executable/build failures carry no operation and stay unattempted.
            if action == "match-core" and exc.operation == action:
                state["invoked"] = True
            raise
        if action == "match-core":
            state["invoked"] = True
            save(work / "native_match_report.json", payload)
            for label, source in zip(("crust_internal", "core_internal", "matched_internal"), args[:3]):
                shutil.copytree(source, work / label)
        return payload
    eos_api.run_native = observed
    try:
        yield
    finally:
        eos_api.run_native = original


def _verify_serialized(directory: Path, work: Path, expected_n: int, *, expected_internal: Path | None = None) -> dict:
    back = work / "readback_internal"
    back.mkdir(parents=True, exist_ok=False)
    read_report = run_native("read", directory, back)
    if read_report.get("status") != "ok" or read_report.get("n") != expected_n:
        raise ValueError("serialized native readback failed or row count differs")
    result = dict(readback_report=read_report)
    if expected_internal is not None:
        numeric = pipeline.verify_internal_readback(expected_internal / "eos_data.csv", back / "eos_data.csv", expected_n=expected_n, actual_n=read_report["n"])
        metadata = pipeline.verify_internal_metadata_readback(expected_internal / "eos_data.nml", back / "eos_data.nml")
        result.update(readback_validation=numeric, metadata_validation=metadata)
        save(work / "readback_validation.json", result)
        if numeric["status"] != "PASS" or metadata["status"] != "PASS":
            raise ValueError("serialized numeric/metadata readback mismatch")
    verification_path = work / "consistency_verification.json"
    report = run_native("verify", directory, verification_path)
    consistency = _json(verification_path)
    errors = pipeline._consistency_report_errors(consistency, expected_n)
    result.update(consistency_native=report, consistency_report=consistency)
    if report.get("status") != "ok" or report.get("overall") != "PASS" or errors:
        raise ValueError(f"serialized consistency verification failed: {errors or report}")
    return result


def run_smooth(crust: Path, core: Path, destination: Path, cell: dict) -> dict:
    work = destination / "work"
    work.mkdir(parents=True, exist_ok=False)
    state = {"invoked": False}
    stage = "read_sources"
    try:
        cr, co = read_compose(crust), read_compose(core)
        stage = "matching"
        with _observe_smooth_kernel(state, work):
            matched, info = eos_api.match_crust_core(cr, co, cr_name=cell["crust"], co_name=cell["core"],
                                                   method=cell["method"], n_bridge=cell["n_bridge"], drop_composition=True)
        save(work / "api_match_report.json", info)
        stage = "write"
        write_compose(destination / "compose", matched)
        stage = "consistency"
        verified = _verify_serialized(destination / "compose", work, matched.n, expected_internal=work / "matched_internal")
        return dict(status="CONSTRUCTED_VERIFIED", matching_kernel_invoked=state["invoked"],
                    method=cell["method"], report=info, rows=matched.n, verification=verified,
                    output_sha256=pipeline._source_compose_hashes(destination / "compose"),
                    evidence_scope="native_smooth_construction_and_serialized_verification; no independent method-specific admission")
    except Exception as exc:
        return dict(status="REJECTED" if state["invoked"] else "NOT_EVALUATED", stage=stage,
                    matching_kernel_invoked=state["invoked"], error=f"{type(exc).__name__}: {exc}",
                    failure_category=classify_failure(stage, str(exc), getattr(exc, "category", None)))


def census(cells: list[dict]) -> dict:
    return dict(requested=len(cells), eligible=sum(c["eligible"] for c in cells),
                new_kernel_attempts=sum(c["attempted"] for c in cells),
                new_admitted=sum(c["attempted"] and c["admitted"] for c in cells),
                new_constructed_verified=sum(c["attempted"] and c["status"] == "CONSTRUCTED_VERIFIED" for c in cells),
                reused_admitted=sum(c["reused"] and c["admitted"] for c in cells),
                outcomes=dict(Counter(c["status"] for c in cells)),
                by_method={method: dict(Counter(c["status"] for c in cells if c["method"] == method)) for method in METHODS})


def execute(root: Path, output: Path, *, crusts=CRUSTS, cores=CORES, methods=METHODS, n_bridge=40,
            reuse_execution: Path | None = None) -> dict:
    root, output = root.resolve(), output.resolve()
    if not output.is_relative_to(root):
        raise ValueError("campaign output must be inside the repository for stable artifact references")
    if isinstance(n_bridge, bool) or not isinstance(n_bridge, int) or n_bridge < 8:
        raise ValueError("n_bridge must be an integer >= 8")
    for values, allowed in ((crusts, CRUSTS), (cores, CORES), (methods, METHODS)):
        if not values or len(values) != len(set(values)) or not set(values).issubset(allowed):
            raise ValueError("empty, duplicate, or out-of-cohort selection")
    output.mkdir(parents=True, exist_ok=False)
    sources, issues = {}, {}
    for family, models in (("crust", crusts), ("core", cores)):
        for model in models:
            sid = f"{family}:{model}"
            try:
                sources[sid] = source_snapshot(root, family, model)
            except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
                issues[sid] = f"{type(exc).__name__}: {exc}"
    cells = [dict(id=f"{crust}__{core}__{method}", crust=crust, core=core, method=method,
                  crust_source_id=f"crust:{crust}", core_source_id=f"core:{core}", n_bridge=n_bridge,
                  requested=True, eligible=not any(x in issues for x in (f"crust:{crust}", f"core:{core}")),
                  attempted=False, admitted=False, reused=False, status="PENDING")
             for crust in crusts for core in cores for method in methods]
    deps = sorted((root / "fortran_eos/src").glob("*.f90")) + sorted((root / "fortran_eos/app").glob("*.f90")) + sorted((root / "python_frontend/compose_eos").glob("*.py"))
    deps += [root / "scripts" / name for name in ("generated_pair_campaign.py", "match_core_crust.py", "match_report_admission.py")]
    initial = {str(p.relative_to(root)): sha(p) for p in deps if p.is_file()}
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, check=True).stdout.strip()
    receipt = dict(schema="ns-cowling.generated-pair-campaign.v1", base_commit=head,
                   started_utc=datetime.now(timezone.utc).isoformat(), scope="finite inner-crust+core; no surface/stellar/Cowling certification",
                   matching_physics="native_Fortran", composition_policy="explicit_thermo_only",
                   workers=1, threads=1, sources=list(sources.values()), source_issues=issues,
                   source_sha256_before=initial, cells=cells, status="RUNNING")
    save(output / "MANIFEST.json", receipt)
    save(output / "RECEIPT.json", receipt)
    for cell in cells:
        destination = output / "pairs" / cell["id"]
        if not cell["eligible"]:
            cell.update(status="HOLD_SOURCE", error={sid: issues[sid] for sid in (cell["crust_source_id"], cell["core_source_id"]) if sid in issues})
        else:
            stage = "source_revalidation"
            try:
                cr, co = (sources[cell[f"{side}_source_id"]] for side in ("crust", "core"))
                validate_source(root, cr)
                validate_source(root, co)
                stage = "executable"
                exe = executable()
                if not os.access(exe, os.X_OK):
                    raise PermissionError(f"native executable is not executable: {exe}")
                cell["verification_executable_sha256"] = sha(exe)
                if reuse_execution is not None and cell["id"] == REUSE_CELL:
                    stage = "reuse_identity"
                    prior = reuse_identity(root, reuse_execution, cell, sources)
                    result = prior["execution"]["result"]
                    stage = "reuse_readback"
                    verified = _verify_serialized(prior["directory"] / "compose", destination / "work", result["combined_n"])
                    admission = admit_match_report(cell["method"], result["match_report"], destination / "work/readback_internal/eos_data.csv", destination / "work/readback_internal/eos_data.nml")
                    save(destination / "work/admission.json", admission)
                    if admission["status"] != "ADMITTED":
                        raise ValueError("reused output fails current independent admission")
                    cell.update(status="REUSED_ADMITTED", reused=True, admitted=True,
                                compose_path=str((prior["directory"] / "compose").relative_to(root)),
                                reuse_execution_path=str(reuse_execution.relative_to(root)), reuse_execution_sha256=sha(reuse_execution),
                                historical_matching_executable_sha256=result["fortran_executable_sha256"],
                                output_sha256=result["output_sha256"], verification=verified, admission=admission,
                                match_report=result["match_report"])
                elif cell["method"] == "quintic_P":
                    stage = "smooth_pipeline"
                    result = run_smooth(_safe(root, cr["compose_path"]), _safe(root, co["compose_path"]), destination, cell)
                    cell.update(status=result["status"], attempted=bool(result["matching_kernel_invoked"]), result=result)
                    if result["status"] == "CONSTRUCTED_VERIFIED":
                        cell["compose_path"] = str((destination / "compose").relative_to(root))
                else:
                    stage = "pipeline"
                    old_n_bridge = pipeline.N_BRIDGE
                    try:
                        pipeline.N_BRIDGE = n_bridge
                        result = pipeline.run_pipeline(cell["method"], _safe(root, cr["compose_path"]), _safe(root, co["compose_path"]), destination)
                    finally:
                        pipeline.N_BRIDGE = old_n_bridge
                    cell.update(status="ADMITTED" if result.get("survived") else "REJECTED" if result.get("matching_kernel_invoked") else "NOT_EVALUATED",
                                attempted=bool(result.get("matching_kernel_invoked")), admitted=bool(result.get("survived")), result=result)
                    if cell["admitted"]:
                        cell["compose_path"] = str((destination / "compose").relative_to(root))
            except Exception as exc:
                cell.update(status="HOLD_REUSE" if stage.startswith("reuse") else "NOT_EVALUATED", stage=stage,
                            error=f"{type(exc).__name__}: {exc}", failure_category="data_provenance" if stage in ("source_revalidation", "reuse_identity") else classify_failure(stage, str(exc), getattr(exc, "category", None)))
        save(destination / "PAIR.json", cell)
        receipt["census"] = census(cells)
        save(output / "RECEIPT.json", receipt)
        print(cell["id"], cell["status"], flush=True)
    final = {str(p.relative_to(root)): sha(p) for p in deps if p.is_file()}
    receipt.update(finished_utc=datetime.now(timezone.utc).isoformat(), source_sha256_after=final,
                   sources_unchanged=initial == final, status="COMPLETE" if initial == final else "IDENTITY_CHANGED")
    save(output / "RECEIPT.json", receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--crusts", nargs="+", choices=CRUSTS, default=CRUSTS)
    parser.add_argument("--cores", nargs="+", choices=CORES, default=CORES)
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=METHODS)
    parser.add_argument("--n-bridge", type=int, default=40)
    parser.add_argument("--reuse-execution", type=Path, default=ROOT / REUSE_EXECUTION)
    parser.add_argument("--no-reuse", action="store_true")
    args = parser.parse_args()
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                      COMPOSE_EOS_FORTRAN_EXE=str(args.build.resolve() / "compose_eos_fortran"))
    result = execute(ROOT, args.out.resolve(), crusts=args.crusts, cores=args.cores, methods=args.methods,
                     n_bridge=args.n_bridge, reuse_execution=None if args.no_reuse else args.reuse_execution.resolve())
    print(json.dumps(result["census"], indent=2))
    return 0 if result["status"] == "COMPLETE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
