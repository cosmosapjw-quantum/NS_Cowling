"""M7: fail-closed orchestrator for the historical EOS390 method matrix.

Python owns manifest/path validation, native invocation, receipt generation and
presentation only. It is not a numerical fallback. A task is PASS only when a
native Fortran backend proves complete crust-core output and successful
CompOSE read-back. Python oracles and placeholders cannot generate production
PASS receipts.
"""
from __future__ import annotations
import json
import os
import sys
import time
import traceback
import hashlib
from pathlib import Path

FIRST_LAW_HERMITE_PRESSURE_BRIDGE = "first_law_constrained_hermite_pressure_bridge"

REPO_ROOT = Path(__file__).resolve().parents[2]
PYTHON_FRONTEND = REPO_ROOT / "python_frontend"
if str(PYTHON_FRONTEND) not in sys.path:
    sys.path.insert(0, str(PYTHON_FRONTEND))


class BackendUnavailable(RuntimeError):
    """The requested production kernel is not yet available in native Fortran."""


NATIVE_PASS_EVIDENCE = (
    "execution_backend", "complete_combined_eos", "compose_readback_pass",
    "fortran_executable_sha256", "input_source_sha256", "physics_contract_id",
)


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def run_first_law_constrained_hermite_pressure_bridge_task(crust_dir: str, core_dir: str, output_dir: str) -> dict:
    """Run the historical smooth NPPHB parameterization under its new method ID.

    The native Fortran smooth Hermite kernel is not ported/validated yet, so the
    390-matrix task remains blocked rather than calling a Python numerical oracle.
    """
    raise BackendUnavailable(
        "first_law_constrained_hermite_pressure_bridge requires the native Fortran M2/M10 smooth Hermite kernel"
    )


def run_maxwell_task(crust_dir: str, core_dir: str, output_dir: str) -> dict:
    """The canonical 'maxwell' method is the true sharp local-neutrality transition.

    This separate physical transition path requires common P(mu_B) coexistence,
    equal P_t and mu_t, a positive density/energy jump, and composition-aware
    native M5 integration. It must not call the historical NPPHB kernel.
    """
    raise BackendUnavailable(
        "maxwell is reserved for sharp Maxwell coexistence; native Fortran M5 composition-aware transition integration is pending"
    )


def run_css_task(crust_dir: str, core_dir: str, output_dir: str) -> dict:
    """Require explicit CSS parameters and a native full-table join.

    The former Python backend hard-coded Delta_epsilon=0.5*epsilon_low and
    c_s^2/c^2=1/3, then wrote only the original core EOS plus an analytic sidecar.
    This is not a complete combined crust-core table. No parameter values are
    inferred from Alford, Han & Prakash, PRD 88, 083013 (2013), arXiv:1302.4732.
    """
    raise BackendUnavailable(
        "CSS is blocked pending source-backed per-task parameters and a native Fortran "
        "complete crust-core branch join with CompOSE read-back"
    )


def run_smooth_task(method: str, crust_dir: str, core_dir: str, output_dir: str) -> dict:
    """Do not report a copied crust table as a successful interpolation."""
    if method not in ("cubic_F", "cubic_P", "quintic_F", "quintic_P"):
        raise ValueError(f"unsupported smooth method: {method}")
    raise BackendUnavailable(
        f"{method} requires the native Fortran M2/M10 interpolation kernel; "
        "Python placeholder copies are forbidden"
    )


def run_task(task: dict) -> dict:
    """Run one task without overwriting prior evidence.

    PASS requires native Fortran provenance, a complete combined EOS,
    successful CompOSE read-back, executable/input hashes and a physics
    contract identifier. A function that merely returns normally is not PASS.
    """
    out_dir = Path(task["output_folder"])
    if out_dir.exists():
        raise FileExistsError(
            f"refusing to overwrite existing task output directory: {out_dir}; "
            "choose a new versioned output root"
        )
    out_dir.mkdir(parents=True, exist_ok=False)
    receipt = {
        "task_id": task["task_id"],
        "crust_model": task["crust_model"],
        "core_model": task["core_model"],
        "method": task["method"],
        "method_kind": task["method_kind"],
        "method_reference": task["method_reference"],
        "acceptance_node": task["acceptance_node"],
        "output_folder": str(out_dir),
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "status": "RUNNING",
    }
    if task["status"] == "BLOCKED_MISSING_SOURCE":
        receipt["status"] = "BLOCKED"
        receipt["blocked_reason"] = task["blocked_reason"]
        receipt["completed_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        return receipt
    try:
        crust_dir = task["crust_eos_path"]
        core_dir = task["core_eos_path"]
        if not Path(crust_dir).is_dir():
            raise FileNotFoundError(f"crust EOS source missing: {crust_dir}")
        if not Path(core_dir).is_dir():
            raise FileNotFoundError(f"core EOS source missing: {core_dir}")
        method = task["method"]
        if method == FIRST_LAW_HERMITE_PRESSURE_BRIDGE:
            result = run_first_law_constrained_hermite_pressure_bridge_task(crust_dir, core_dir, str(out_dir))
        elif method == "maxwell":
            raise BackendUnavailable(
                "maxwell is the sharp first-order transition method and is not part of the six-method smooth EOS390 v2 matrix; "
                "use the separate M5 sharp-Maxwell workflow"
            )
        elif method == "css":
            result = run_css_task(crust_dir, core_dir, str(out_dir))
        elif method in ("cubic_F", "cubic_P", "quintic_F", "quintic_P"):
            result = run_smooth_task(method, crust_dir, core_dir, str(out_dir))
        else:
            raise ValueError(f"unknown method: {method}")

        missing = [
            key for key in NATIVE_PASS_EVIDENCE
            if key not in result or result[key] in (None, "", False, {})
        ]
        if missing:
            receipt["status"] = "BLOCKED"
            receipt["blocked_reason"] = (
                "backend lacks mandatory native-production evidence: " + ", ".join(missing)
            )
            receipt["backend_result"] = result
        elif result["execution_backend"] != "native_fortran":
            receipt["status"] = "BLOCKED"
            receipt["blocked_reason"] = "PASS requires execution_backend=native_fortran"
            receipt["backend_result"] = result
        elif result["complete_combined_eos"] is not True or result["compose_readback_pass"] is not True:
            receipt["status"] = "BLOCKED"
            receipt["blocked_reason"] = (
                "native backend did not prove a complete combined EOS and successful CompOSE read-back"
            )
            receipt["backend_result"] = result
        else:
            output_files = result.get("output_files", [])
            if not output_files:
                receipt["status"] = "BLOCKED"
                receipt["blocked_reason"] = "native backend returned no output files"
            else:
                hashes = {}
                for filename in output_files:
                    file_path = out_dir / filename
                    if not file_path.is_file():
                        raise FileNotFoundError(f"backend output is missing: {file_path}")
                    hashes[filename] = sha256_file(str(file_path))
                receipt["status"] = "PASS"
                receipt["output_files"] = output_files
                receipt["file_hashes"] = hashes
                receipt["result"] = result
    except BackendUnavailable as exc:
        receipt["status"] = "BLOCKED"
        receipt["blocked_reason"] = str(exc)
    except Exception as exc:
        receipt["status"] = "FAIL"
        receipt["error_type"] = type(exc).__name__
        receipt["error_message"] = str(exc)[:500]
        receipt["traceback"] = traceback.format_exc()[-1500:]
    receipt["completed_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    return receipt


def run_batch(manifest_path: str, *, max_pending: int = None) -> dict:
    """Run all tasks in the manifest. Returns the batch summary.
    max_pending: if set, cap the number of PENDING tasks actually dispatched
    (useful in the sandbox to avoid long runs; None = run all)."""
    with open(manifest_path, 'r', encoding='utf-8') as f:
        manifest = json.load(f)
    if manifest.get("schema") != "ns-cowling.eos390-m7-batch-manifest" or manifest.get("schema_version") != "2.0.0":
        raise ValueError(
            "M7 accepts only the v2 six-smooth-method manifest. Historical v1 manifests used "
            "method=maxwell for NPPHB and are immutable; regenerate a v2 manifest with "
            "first_law_constrained_hermite_pressure_bridge. The sharp method ID maxwell "
            "belongs to the separate M5 transition workflow."
        )
    if "maxwell" in manifest.get("matrix", {}).get("methods", []):
        raise ValueError("method=maxwell is forbidden in the M7 smooth EOS390 manifest")
    tasks = manifest['tasks']
    receipts = []
    counts = {'PASS': 0, 'FAIL': 0, 'BLOCKED': 0, 'PARTIAL': 0, 'PENDING': 0}
    pending_dispatched = 0
    for i, task in enumerate(tasks):
        if task['status'] == 'PENDING' and max_pending is not None and pending_dispatched >= max_pending:
            # Leave remaining PENDING tasks as PENDING (not dispatched)
            receipt = {
                'task_id': task['task_id'],
                'crust_model': task['crust_model'],
                'core_model': task['core_model'],
                'method': task['method'],
                'status': 'NOT_DISPATCHED_CAP',
                'reason': f'batch capped at {max_pending} pending; {len(tasks)-i} tasks remain',
            }
            receipts.append(receipt)
            counts['PENDING'] += 1
            continue
        receipt = run_task(task)
        receipts.append(receipt)
        counts[receipt['status']] = counts.get(receipt['status'], 0) + 1
        if task['status'] == 'PENDING':
            pending_dispatched += 1
        # Write the per-task receipt into the task's output folder
        if receipt.get('output_folder'):
            rpath = Path(receipt['output_folder']) / 'task_receipt.json'
            rpath.write_text(json.dumps(receipt, indent=2, allow_nan=False,
                                         default=str) + '\n', encoding='utf-8')
        else:
            # NOT_DISPATCHED_CAP tasks have no folder yet; write to a cap-area
            cap_dir = Path(manifest_path).expanduser().resolve().parent / 'outputs' / '_cap_receipts'
            cap_dir.mkdir(parents=True, exist_ok=True)
            (cap_dir / f'{receipt["task_id"]}_receipt.json').write_text(
                json.dumps(receipt, indent=2, default=str) + '\n', encoding='utf-8')
    summary = {
        'schema': 'ns-cowling.eos390-m7-batch-summary',
        'schema_version': '1.0.0',
        'manifest_path': manifest_path,
        'total_tasks': len(tasks),
        'receipts_count': len(receipts),
        'status_counts': counts,
        'cardinality_check': len(tasks) == 390,
        'uniqueness_check': len(set(t['task_id'] for t in tasks)) == 390,
        'no_silent_skip_check': len(receipts) == len(tasks),
        'fail_closed_check': all(r['status'] in ('PASS','FAIL','BLOCKED','PARTIAL','NOT_DISPATCHED_CAP','PENDING') for r in receipts),
        'pending_attempted': pending_dispatched,
        'pending_dispatched': pending_dispatched,
        'max_pending_cap': max_pending,
    }
    return summary, receipts


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="exact versioned 390-key manifest")
    ap.add_argument("--summary", required=True, help="new versioned summary JSON path")
    ap.add_argument("--max-pending", type=int, default=None,
                    help="optional cap on attempted tasks; omitted means all manifest rows")
    args = ap.parse_args()
    summary, receipts = run_batch(args.manifest, max_pending=args.max_pending)
    output = Path(args.summary).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        json.dump({"summary": summary, "receipts": receipts}, stream,
                  indent=2, allow_nan=False, default=str)
        stream.write("\\n")
    print(json.dumps(summary, indent=2))
    print(f"\\nbatch summary written: {output}")
