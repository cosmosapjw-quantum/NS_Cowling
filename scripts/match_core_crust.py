"""Run native core--crust matching and verify the serialized CompOSE readback.

Runs the native Fortran Maxwell and CSS methods for each configured crust/core
pair. A pair survives only if the native match, CompOSE write, CompOSE readback,
reported point counts, and every numerical/metadata column in the internal
readback CSV all agree.

Physics:
  - maxwell: sharp first-order local-neutrality coexistence, with common
    pressure and baryon chemical potential and a finite density/energy jump.
  - css: one constant-sound-speed CSS bridge. At crust--CSS, P and mu_B are
    continuous while n_B and epsilon jump; at CSS--core, n_B, P, epsilon and
    mu_B are continuous. This is not a smooth density-continuous join over
    the complete crust-plus-bridge domain.

Both methods are thermodynamic-only in this pipeline. Composition is explicitly
dropped; no composition-preservation claim is made.
"""
from __future__ import annotations
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos.native import executable, run_native  # noqa: E402

CRUST_MODELS = [
    "SPG_M1", "SPG_M2", "SPG_M3", "SPG_M4", "SPG_M5",
    "GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW",
]
CORE_MODELS = ["APR_APR", "BBB_BHF2", "GM_GM1", "RG_SLY4", "RG_KDE0v1"]
CRUST_ROOT = ROOT / "download" / "reconstruction"
CORE_ROOT = ROOT / "eos_generator" / "eos_data"
DEFAULT_RUN_STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
DEFAULT_OUT_ROOT = ROOT / "download" / "matched_runs" / DEFAULT_RUN_STAMP
OUT_ROOT = Path(os.environ.get("NS_COWLING_MATCH_OUT_ROOT", str(DEFAULT_OUT_ROOT))).expanduser()
N_BRIDGE = 8


def sha256_file(path: Path | str) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# Integer-valued metadata columns in the native Fortran internal CSV.
_INTEGER_COLUMNS = {
    "Iphase", "row_status", "yn_present", "yp_present", "ye_present",
    "ymu_present", "abar_present", "zbar_present", "yh_present",
    "c3_applicable", "c4_applicable", "iy_local", "it_local",
}

_METADATA_INTEGER_KEYS = {
    "il", "iT", "iYq", "inb_start", "nadd", "nadd_max", "naxis_t", "naxis_yq",
}
_METADATA_BOOLEAN_KEYS = {
    "composition_file_present", "composition_available", "beta_equilibrated",
}
_METADATA_FLOAT_KEYS = {
    "mn", "mp", "temperature_mev", "charge_fraction_grid",
}
_METADATA_REQUIRED_KEYS = (
    _METADATA_INTEGER_KEYS | _METADATA_BOOLEAN_KEYS | _METADATA_FLOAT_KEYS
)


def verify_internal_readback(
    expected_csv: Path | str,
    actual_csv: Path | str,
    *,
    expected_n: int,
    actual_n: int,
    rtol: float = 5e-13,
    atol: float = 1e-14,
) -> dict:
    """Compare the full native internal EOS table before/after CompOSE I/O.

    This is intentionally stricter than checking a successful process exit or
    an unchanged row count. Every numeric column emitted by write_internal_csv
    is compared row-by-row, including Q1..Q7, q_extra slots, phase labels,
    composition values, row status, presence masks and C3/C4 applicability.
    Integer metadata is exact; floating-point columns use the declared
    mixed absolute/relative tolerance.
    """
    expected_path = Path(expected_csv)
    actual_path = Path(actual_csv)
    errors: list[str] = []
    report = {
        "schema": "native-compose-internal-readback-v1",
        "status": "FAIL",
        "expected_n_reported": int(expected_n),
        "actual_n_reported": int(actual_n),
        "expected_path": str(expected_path),
        "actual_path": str(actual_path),
        "rtol": float(rtol),
        "atol": float(atol),
        "columns": [],
        "expected_rows": None,
        "actual_rows": None,
        "max_abs_residuals": {},
        "max_rel_residuals": {},
        "errors": errors,
        "metadata_status": "NOT_CHECKED",
        "metadata_max_abs_residuals": {},
        "metadata_max_rel_residuals": {},
    }

    if expected_n < 0 or actual_n < 0:
        errors.append("native readback row counts must be nonnegative")
    if int(expected_n) != int(actual_n):
        errors.append(
            f"native readback row count mismatch: expected {expected_n}, got {actual_n}"
        )
    if not expected_path.is_file():
        errors.append(f"expected internal CSV is missing: {expected_path}")
    if not actual_path.is_file():
        errors.append(f"readback internal CSV is missing: {actual_path}")
    if errors:
        return report

    def load_table(path: Path):
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            header = reader.fieldnames
            if not header:
                raise ValueError(f"{path}: missing CSV header")
            if len(set(header)) != len(header):
                raise ValueError(f"{path}: duplicate CSV column names")
            rows = list(reader)
            for line_number, row in enumerate(rows, start=2):
                if None in row:
                    raise ValueError(
                        f"{path}:{line_number}: row has more values than the header"
                    )
                if any(value is None for value in row.values()):
                    raise ValueError(
                        f"{path}:{line_number}: row has fewer values than the header"
                    )
            return list(header), rows

    try:
        expected_header, expected_rows = load_table(expected_path)
        actual_header, actual_rows = load_table(actual_path)
    except (OSError, UnicodeError, csv.Error, ValueError) as exc:
        errors.append(f"internal CSV parse failure: {type(exc).__name__}: {exc}")
        return report

    report["expected_rows"] = len(expected_rows)
    report["actual_rows"] = len(actual_rows)
    report["columns"] = expected_header

    if len(expected_rows) != int(expected_n):
        errors.append(
            f"expected internal CSV row count mismatch: native match reported "
            f"{expected_n}, file contains {len(expected_rows)}"
        )
    if len(actual_rows) != int(actual_n):
        errors.append(
            f"readback internal CSV row count mismatch: native read reported "
            f"{actual_n}, file contains {len(actual_rows)}"
        )
    expected_columns = set(expected_header)
    actual_columns = set(actual_header)
    if expected_columns != actual_columns:
        missing = sorted(expected_columns - actual_columns)
        extra = sorted(actual_columns - expected_columns)
        errors.append(
            f"internal CSV column mismatch: missing={missing}, extra={extra}"
        )
    if len(expected_rows) != len(actual_rows):
        errors.append(
            f"internal CSV row count mismatch: expected {len(expected_rows)}, "
            f"got {len(actual_rows)}"
        )
        return report
    if expected_columns != actual_columns:
        return report

    # Validate the serialized schema needed to preserve the EOS contract.
    required = {"nb", *(f"Q{i}" for i in range(1, 8)), "Iphase"}
    missing_required = sorted(required - expected_columns)
    if missing_required:
        errors.append(f"internal CSV lacks required EOS columns: {missing_required}")
        return report

    max_abs = {name: 0.0 for name in expected_header}
    max_rel = {name: 0.0 for name in expected_header}
    for row_index, (expected_row, actual_row) in enumerate(
        zip(expected_rows, actual_rows), start=1
    ):
        for name in expected_header:
            try:
                expected_value = float(expected_row[name])
                actual_value = float(actual_row[name])
            except (TypeError, ValueError):
                errors.append(
                    f"row {row_index} column {name}: non-numeric value "
                    f"(expected={expected_row.get(name)!r}, actual={actual_row.get(name)!r})"
                )
                continue
            if not (math.isfinite(expected_value) and math.isfinite(actual_value)):
                errors.append(
                    f"row {row_index} column {name}: non-finite value "
                    f"(expected={expected_value!r}, actual={actual_value!r})"
                )
                continue
            difference = abs(actual_value - expected_value)
            scale = max(abs(expected_value), abs(actual_value), 1e-30)
            max_abs[name] = max(max_abs[name], difference)
            max_rel[name] = max(max_rel[name], difference / scale)

            if name in _INTEGER_COLUMNS:
                if not expected_value.is_integer() or not actual_value.is_integer():
                    errors.append(
                        f"row {row_index} column {name}: integer metadata is not integral"
                    )
                elif int(expected_value) != int(actual_value):
                    errors.append(
                        f"row {row_index} column {name}: exact metadata mismatch "
                        f"(expected={int(expected_value)}, actual={int(actual_value)})"
                    )
            elif difference > atol + rtol * abs(expected_value):
                errors.append(
                    f"row {row_index} column {name}: numerical mismatch "
                    f"(expected={expected_value:.17g}, actual={actual_value:.17g}, "
                    f"abs={difference:.3e}, rel={difference/scale:.3e})"
                )

    report["max_abs_residuals"] = max_abs
    report["max_rel_residuals"] = max_rel
    report["status"] = "PASS" if not errors else "FAIL"
    return report


def _parse_fortran_metadata(path: Path) -> dict[str, bool | int | float]:
    """Parse the scalar fields emitted by the native write_metadata routine."""
    text = path.read_text(encoding="utf-8")
    values = {}
    for match in re.finditer(
        r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*([^,\s/]+)",
        text,
        flags=re.MULTILINE,
    ):
        key, raw = match.group(1), match.group(2).strip()
        if key in _METADATA_BOOLEAN_KEYS:
            normalized = raw.strip(".").upper()
            if normalized in {"T", "TRUE"}:
                value = True
            elif normalized in {"F", "FALSE"}:
                value = False
            else:
                raise ValueError(f"{path}: invalid Fortran logical {key}={raw!r}")
        elif key in _METADATA_INTEGER_KEYS:
            value_float = float(raw.replace("D", "E").replace("d", "e"))
            if not value_float.is_integer():
                raise ValueError(f"{path}: non-integral metadata {key}={raw!r}")
            value = int(value_float)
        else:
            value = float(raw.replace("D", "E").replace("d", "e"))
            if not math.isfinite(value):
                raise ValueError(f"{path}: non-finite metadata {key}={raw!r}")
        if key in values:
            raise ValueError(f"{path}: duplicate metadata key {key}")
        values[key] = value
    missing = sorted(_METADATA_REQUIRED_KEYS - set(values))
    if missing:
        raise ValueError(f"{path}: missing required metadata keys {missing}")
    return values


def verify_internal_metadata_readback(
    expected_metadata: Path | str,
    actual_metadata: Path | str,
    *,
    rtol: float = 5e-13,
    atol: float = 1e-14,
) -> dict:
    """Compare all scalar EOS metadata written before and after CompOSE I/O."""
    expected_path = Path(expected_metadata)
    actual_path = Path(actual_metadata)
    report = {
        "status": "FAIL",
        "expected_path": str(expected_path),
        "actual_path": str(actual_path),
        "rtol": float(rtol),
        "atol": float(atol),
        "max_abs_residuals": {},
        "max_rel_residuals": {},
        "expected_values": {},
        "actual_values": {},
        "errors": [],
    }
    if not expected_path.is_file():
        report["errors"].append(f"expected internal metadata is missing: {expected_path}")
    if not actual_path.is_file():
        report["errors"].append(f"readback internal metadata is missing: {actual_path}")
    if report["errors"]:
        return report
    try:
        expected = _parse_fortran_metadata(expected_path)
        actual = _parse_fortran_metadata(actual_path)
    except (OSError, UnicodeError, ValueError) as exc:
        report["errors"].append(f"internal metadata parse failure: {type(exc).__name__}: {exc}")
        return report
    report["expected_values"] = expected
    report["actual_values"] = actual
    if set(expected) != set(actual):
        report["errors"].append(
            f"internal metadata key mismatch: missing={sorted(set(expected)-set(actual))}, "
            f"extra={sorted(set(actual)-set(expected))}"
        )
        return report
    max_abs = {}
    max_rel = {}
    for key in sorted(expected):
        ev, av = expected[key], actual[key]
        if key in _METADATA_BOOLEAN_KEYS or key in _METADATA_INTEGER_KEYS:
            if type(ev) is not type(av) or ev != av:
                report["errors"].append(
                    f"metadata {key}: exact mismatch (expected={ev!r}, actual={av!r})"
                )
        else:
            difference = abs(float(av) - float(ev))
            scale = max(abs(float(ev)), abs(float(av)), 1e-30)
            max_abs[key] = difference
            max_rel[key] = difference / scale
            if difference > atol + rtol * abs(float(ev)):
                report["errors"].append(
                    f"metadata {key}: numerical mismatch (expected={ev:.17g}, actual={av:.17g}, "
                    f"abs={difference:.3e}, rel={difference/scale:.3e})"
                )
    report["max_abs_residuals"] = max_abs
    report["max_rel_residuals"] = max_rel
    report["status"] = "PASS" if not report["errors"] else "FAIL"
    return report


def _source_compose_hashes(directory: Path) -> dict[str, str]:
    """Return hashes of the CompOSE source files actually present."""
    hashes = {}
    for name in ("eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo"):
        path = directory / name
        if path.is_file():
            hashes[name] = sha256_file(path)
    return hashes


def run_pipeline(method: str, crust_dir: Path, core_dir: Path, out_dir: Path) -> dict:
    """Run one (crust, core, method) match. Returns a result dict.

    method in {"maxwell", "css"}.
    """
    work = out_dir / "work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    (out_dir / "compose").mkdir(parents=True, exist_ok=True)

    crust_internal = work / "crust_internal"
    core_internal = work / "core_internal"
    matched_internal = work / "matched_internal"
    crust_internal.mkdir(exist_ok=True)
    core_internal.mkdir(exist_ok=True)
    matched_internal.mkdir(exist_ok=True)

    exe = executable()
    exe_sha = sha256_file(exe)

    try:
        # 1. read crust -> internal
        r1 = run_native("read", str(crust_dir), str(crust_internal))
        if r1.get("status") != "ok":
            return {"survived": False, "method": method, "stage": "read_crust",
                    "error": r1.get("message", "?")}
        # 2. read core -> internal
        r2 = run_native("read", str(core_dir), str(core_internal))
        if r2.get("status") != "ok":
            return {"survived": False, "method": method, "stage": "read_core",
                    "error": r2.get("message", "?")}

        if method == "maxwell":
            # maxwell-transition: low=crust, high=core, local_neutrality=true
            r3 = run_native("maxwell-transition", str(crust_internal), str(core_internal),
                            str(matched_internal), "true")
        else:  # css
            r3 = run_native("match-core", str(crust_internal), str(core_internal),
                            str(matched_internal), "css", str(N_BRIDGE), "true")
        if r3.get("status") != "ok":
            return {"survived": False, "method": method, "stage": method,
                    "error": r3.get("message", "?"),
                    "report": r3}

        combined_n = int(r3.get("n", 0))
        if combined_n < 4:
            return {"survived": False, "method": method, "stage": "match",
                    "error": f"combined n={combined_n} too small", "report": r3}

        # 4. write internal -> CompOSE output
        r4 = run_native("write", str(matched_internal), str(out_dir / "compose"))
        if r4.get("status") != "ok":
            return {"survived": False, "method": method, "stage": "write",
                    "error": r4.get("message", "?")}

        # 5. readback verification
        readback_internal = work / "readback_internal"
        readback_internal.mkdir(exist_ok=True)
        try:
            r5 = run_native("read", str(out_dir / "compose"), str(readback_internal))
        except Exception as exc:
            return {
                "survived": False, "method": method, "stage": "readback",
                "error": f"native CompOSE readback failed: {type(exc).__name__}: {exc}",
                "combined_n": combined_n,
                "readback_n": 0,
                "readback_status": "ERROR",
                "readback_validation": {
                    "schema": "native-compose-internal-readback-v1",
                    "status": "FAIL",
                    "errors": [f"native CompOSE readback failed: {type(exc).__name__}: {exc}"],
                },
                "match_report": r3,
                "write_report": r4,
            }

        readback_n = int(r5.get("n", 0)) if r5.get("status") == "ok" else 0
        expected_csv = matched_internal / "eos_data.csv"
        actual_csv = readback_internal / "eos_data.csv"
        expected_nml = matched_internal / "eos_data.nml"
        actual_nml = readback_internal / "eos_data.nml"
        readback_validation = verify_internal_readback(
            expected_csv, actual_csv,
            expected_n=combined_n, actual_n=readback_n,
        )
        metadata_validation = verify_internal_metadata_readback(
            expected_nml, actual_nml,
        )
        readback_validation["metadata_status"] = metadata_validation["status"]
        readback_validation["metadata"] = metadata_validation
        readback_validation["metadata_max_abs_residuals"] = metadata_validation["max_abs_residuals"]
        readback_validation["metadata_max_rel_residuals"] = metadata_validation["max_rel_residuals"]
        readback_validation["errors"].extend(metadata_validation["errors"])
        if metadata_validation["status"] != "PASS":
            readback_validation["status"] = "FAIL"
        if r5.get("status") != "ok":
            readback_validation["status"] = "FAIL"
            readback_validation["errors"].append(
                f"native readback returned status={r5.get('status')!r}"
            )
        if readback_validation["status"] != "PASS":
            return {
                "survived": False,
                "method": method,
                "stage": "readback",
                "error": "; ".join(readback_validation["errors"][:8]),
                "combined_n": combined_n,
                "readback_n": readback_n,
                "readback_status": r5.get("status", "?"),
                "readback_validation": readback_validation,
                "fortran_executable_sha256": exe_sha,
                "match_report": r3,
                "write_report": r4,
            }

        # 6. Verify the actual serialized CompOSE artefact with the native
        # C1/C2/C3/C4/C5/C6/C10/C11 checker. C7-C9 remain diagnostics only.
        consistency_path = work / "consistency_verification.json"
        try:
            r6 = run_native("verify", str(out_dir / "compose"), str(consistency_path))
        except Exception as exc:
            return {
                "survived": False, "method": method, "stage": "consistency",
                "error": f"native consistency verification failed: {type(exc).__name__}: {exc}",
                "combined_n": combined_n,
                "readback_n": readback_n,
                "readback_status": r5.get("status", "?"),
                "readback_validation": readback_validation,
                "match_report": r3,
                "write_report": r4,
            }
        if r6.get("status") != "ok" or r6.get("overall") != "PASS":
            consistency_report = None
            try:
                consistency_report = json.loads(consistency_path.read_text(encoding="utf-8"))
            except (OSError, UnicodeError, json.JSONDecodeError):
                pass
            return {
                "survived": False,
                "method": method,
                "stage": "consistency",
                "error": f"native consistency overall={r6.get('overall', 'UNKNOWN')}; "
                         f"failed_checks={r6.get('failed_checks', 'unavailable')}",
                "combined_n": combined_n,
                "readback_n": readback_n,
                "readback_status": r5.get("status", "?"),
                "readback_validation": readback_validation,
                "consistency_status": r6.get("overall", "UNKNOWN"),
                "consistency_report": consistency_report,
                "match_report": r3,
                "write_report": r4,
            }
        try:
            consistency_report = json.loads(consistency_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            return {
                "survived": False, "method": method, "stage": "consistency",
                "error": f"native verifier did not produce valid JSON: {type(exc).__name__}: {exc}",
                "combined_n": combined_n,
                "readback_n": readback_n,
                "readback_status": r5.get("status", "?"),
                "readback_validation": readback_validation,
                "consistency_status": "INVALID_REPORT",
                "match_report": r3,
                "write_report": r4,
            }

        output_files = sorted(p.name for p in (out_dir / "compose").iterdir() if p.is_file())
        output_hashes = _source_compose_hashes(out_dir / "compose")
        source_hashes = {
            "crust": _source_compose_hashes(crust_dir),
            "core": _source_compose_hashes(core_dir),
        }

        return {
            "survived": True,
            "method": method,
            "stage": "complete",
            "combined_n": combined_n,
            "readback_n": readback_n,
            "readback_status": r5.get("status", "?"),
            "readback_validation": readback_validation,
            "consistency_status": "PASS",
            "consistency_report": consistency_report,
            "source_sha256": source_hashes,
            "output_sha256": output_hashes,
            "output_files": output_files,
            "fortran_executable_sha256": exe_sha,
            "match_report": r3,
            "write_report": r4,
        }
    except Exception as exc:
        return {"survived": False, "method": method, "stage": "exception",
                "error": f"{type(exc).__name__}: {exc}"}


def main() -> int:
    # Never delete a previous match census. A fresh timestamped run root is
    # the default; an explicit override is allowed only when its directory is
    # absent or empty. This makes reruns and interrupted jobs fail closed.
    if OUT_ROOT.exists():
        if not OUT_ROOT.is_dir():
            raise NotADirectoryError(f"matching output root is not a directory: {OUT_ROOT}")
        if any(OUT_ROOT.iterdir()):
            raise FileExistsError(
                f"matching output root is not empty: {OUT_ROOT}; "
                "choose a new NS_COWLING_MATCH_OUT_ROOT"
            )
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    results = []
    survivors = []
    for crust in CRUST_MODELS:
        crust_dir = CRUST_ROOT / crust
        if not crust_dir.is_dir():
            continue
        for core in CORE_MODELS:
            core_dir = CORE_ROOT / core
            if not core_dir.is_dir():
                continue
            for method in ("maxwell", "css"):
                pair_id = f"{crust}__{core}__{method}"
                out_dir = OUT_ROOT / pair_id
                out_dir.mkdir(parents=True, exist_ok=True)
                t0 = time.time()
                res = run_pipeline(method, crust_dir, core_dir, out_dir)
                res["crust_model"] = crust
                res["core_model"] = core
                res["pair_id"] = pair_id
                res["elapsed_s"] = round(time.time() - t0, 2)
                results.append(res)
                status = "SURVIVED" if res["survived"] else f"FAILED@{res.get('stage')}"
                print(f"  [{status:20s}] {crust:14s} + {core:12s} / {method:7s}  "
                      f"(n={res.get('combined_n','-')})  {res.get('elapsed_s','-')}s  "
                      f"{res.get('error','')[:60]}")
                if res["survived"]:
                    survivors.append({
                        "pair_id": pair_id,
                        "crust_model": crust,
                        "core_model": core,
                        "method": method,
                        "combined_n": res["combined_n"],
                        "readback_n": res["readback_n"],
                        "readback_status": res["readback_status"],
                        "readback_validation": res["readback_validation"],
                        "consistency_status": res["consistency_status"],
                        "consistency_report": res["consistency_report"],
                        "source_sha256": res["source_sha256"],
                        "output_sha256": res["output_sha256"],
                        "fortran_executable_sha256": res["fortran_executable_sha256"],
                        "output_dir": str(out_dir / "compose"),
                        "match_report": res["match_report"],
                    })
                else:
                    # clean up failed attempt dir
                    shutil.rmtree(out_dir, ignore_errors=True)

    # write manifest
    n_maxwell = sum(1 for s in survivors if s["method"] == "maxwell")
    n_css = sum(1 for s in survivors if s["method"] == "css")
    built_at = datetime.now(timezone.utc)
    manifest = {
        "schema": "ns-cowling.core-crust-matching.v2",
        "built_at_utc": built_at.isoformat(),
        "built_date_utc": built_at.date().isoformat(),
        "output_root": str(OUT_ROOT),
        "crust_models": CRUST_MODELS,
        "core_models": CORE_MODELS,
        "methods": ["maxwell", "css"],
        "method_contracts": {
            "maxwell": {
                "construction_kind": "sharp_first_order_local_neutrality_coexistence",
                "interfaces": {
                    "crust_core": "common pressure and baryon chemical potential; finite density and energy-density jump"
                },
                "composition": "dropped",
            },
            "css": {
                "construction_kind": "css_single_speed_with_crust_density_jump",
                "interfaces": {
                    "crust_css": "common pressure and baryon chemical potential; finite positive density and energy-density jump with Delta_epsilon=mu_t Delta_n",
                    "css_core": "continuous n_B, pressure, energy density, and baryon chemical potential",
                },
                "sound_speed": "one constant cs2 across the CSS bridge",
                "composition": "dropped",
                "smooth_over_entire_crust_plus_bridge": False,
            },
        },
        "total_attempts": len(results),
        "total_survivors": len(survivors),
        "survivors_maxwell": n_maxwell,
        "survivors_css": n_css,
        "survivors": survivors,
        "all_results": results,
    }
    manifest_path = OUT_ROOT / "matching_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
                             encoding="utf-8")
    print(f"\nMatching complete: {len(survivors)}/{len(results)} survived "
          f"(maxwell={n_maxwell}, css={n_css})")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
