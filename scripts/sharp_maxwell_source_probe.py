#!/usr/bin/env python3
"""Read-only native reader/verifier probe for the sharp-Maxwell source cohort.

The default vertical slice uses SPG_M1 (Nadd=0) and GPPVA_FSU2 (Nadd=1).
It preserves source identity, checks row-wise parsing, exercises native
read -> write -> read, and runs the native C1-C11 consistency verifier.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data" / "sharp_maxwell_65_cohort_manifest.v1.json"
DEFAULT_SOURCE_ROOT = ROOT / "eos_generator" / "eos_data"
DEFAULT_SOURCE_IDS = ("SPG_M1", "GPPVA_FSU2")
COMPOSE_ARRAYS = (
    "nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "q_extra",
    "Iphase", "Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh",
)


def _json_safe(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if hasattr(value, "item"):
        try:
            return _json_safe(value.item())
        except (TypeError, ValueError):
            pass
    return value


def _fingerprint(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "bytes": len(raw),
        "git_blob_sha1": hashlib.sha1(f"blob {len(raw)}\0".encode("ascii") + raw).hexdigest(),
        "sha256": hashlib.sha256(raw).hexdigest(),
    }


def _meaningful_lines(path: Path) -> list[str]:
    return [
        line.strip() for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def _read_raw_density(path: Path) -> tuple[int, list[float]]:
    tokens = " ".join(_meaningful_lines(path)).split()
    if len(tokens) < 3:
        raise ValueError(f"short eos.nb: {path}")
    il, count = int(tokens[0]), int(tokens[1])
    values = [float(token) for token in tokens[2:]]
    if len(values) != count:
        raise ValueError(f"eos.nb declares {count} knots but contains {len(values)}: {path}")
    return il, values


def _read_raw_thermo(path: Path) -> tuple[list[int], list[list[float]], list[float], tuple[float, float, int]]:
    lines = _meaningful_lines(path)
    if not lines:
        raise ValueError(f"empty eos.thermo: {path}")
    header = lines[0].split()
    if len(header) != 3:
        raise ValueError(f"invalid eos.thermo header: {path}")
    masses = (float(header[0]), float(header[1]), int(header[2]))
    nadds, q_rows, extras = [], [], []
    for line_no, line in enumerate(lines[1:], 2):
        tokens = line.split()
        if len(tokens) < 11:
            raise ValueError(f"short eos.thermo row at {path}:{line_no}")
        q_rows.append([float(x) for x in tokens[3:10]])
        nadd = int(tokens[10])
        if nadd not in (0, 1) or len(tokens) != 11 + nadd:
            raise ValueError(f"unsupported Nadd/row width at {path}:{line_no}")
        nadds.append(nadd)
        extras.append(float(tokens[11]) if nadd == 1 else 0.0)
    return nadds, q_rows, extras, masses


def _read_raw_composition(path: Path) -> tuple[list[int], dict[str, list[float]]]:
    """Decode the selected CompOSE eos.compo fields without reinterpretation."""
    phase = []
    names = ("Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh")
    fields = {name: [] for name in names}
    species = {10: "Yn", 11: "Yp", 0: "Ye", 1: "Ymu"}
    for line_no, line in enumerate(_meaningful_lines(path), 1):
        tokens = line.split()
        if len(tokens) < 6:
            raise ValueError(f"short eos.compo row at {path}:{line_no}")
        phase.append(int(tokens[3]))
        row = {name: 0.0 for name in names}
        npairs = int(tokens[4])
        pos = 5
        for _ in range(npairs):
            if pos + 1 >= len(tokens):
                raise ValueError(f"truncated particle pair at {path}:{line_no}")
            particle_id, fraction = int(tokens[pos]), float(tokens[pos + 1])
            if particle_id not in species:
                raise ValueError(f"unsupported particle id {particle_id} at {path}:{line_no}")
            row[species[particle_id]] = fraction
            pos += 2
        if pos >= len(tokens):
            raise ValueError(f"missing Nquad at {path}:{line_no}")
        nquad = int(tokens[pos])
        pos += 1
        if nquad not in (0, 1):
            raise ValueError(f"unsupported Nquad={nquad} at {path}:{line_no}")
        for _ in range(nquad):
            if pos + 3 >= len(tokens):
                raise ValueError(f"truncated nuclear group at {path}:{line_no}")
            group_id = int(tokens[pos])
            abar, zbar, yi = map(float, tokens[pos + 1:pos + 4])
            if group_id != 1:
                raise ValueError(f"unsupported nuclear-group id {group_id} at {path}:{line_no}")
            # CompOSE Reference Manual Eq. (4.22): the stored group Y_I is
            # already the combined group fraction; do not multiply by Abar.
            row["Aav"], row["Zav"], row["Yh"] = abar, zbar, yi
            pos += 4
        if pos != len(tokens):
            raise ValueError(f"trailing eos.compo fields at {path}:{line_no}")
        for name in names:
            fields[name].append(row[name])
    return phase, fields


def _histogram(values: Any) -> dict[str, int]:
    counts = Counter(int(value) for value in values)
    return {str(key): counts[key] for key in sorted(counts)}


def _extra_histogram(values: Any) -> dict[str, int]:
    converted = []
    for item in values:
        number = float(item)
        nearest = round(number)
        if not math.isfinite(number) or number != nearest:
            raise ValueError(f"q_extra contains a non-integer-valued provider code: {number!r}")
        converted.append(int(nearest))
    return _histogram(converted)


def _phase_boundaries(values: Any) -> list[dict[str, int]]:
    phases = [int(value) for value in values]
    return [
        {"left_knot": i + 1, "right_knot": i + 2, "from": left, "to": right}
        for i, (left, right) in enumerate(zip(phases, phases[1:])) if left != right
    ]


def _expected_boundary(entry: dict[str, Any]) -> list[dict[str, int]]:
    boundary = entry.get("Iphase_boundary")
    if not boundary:
        return []
    return [{
        "left_knot": int(boundary["left_knot"]),
        "right_knot": int(boundary["right_knot"]),
        "from": int(boundary["from"]),
        "to": int(boundary["to"]),
    }]


def _finite_arrays(data: Any) -> dict[str, bool]:
    import numpy as np
    return {
        name: bool(
            len(getattr(data, name)) == data.n
            and np.all(np.isfinite(getattr(data, name)))
        )
        for name in COMPOSE_ARRAYS if name != "Iphase"
    } | {"Iphase": len(data.Iphase) == data.n}


def _rowwise_equal(a: Any, b: Any, *, integer: bool = False, atol: float = 5e-13) -> bool:
    if len(a) != len(b):
        return False
    if integer:
        return all(int(x) == int(y) for x, y in zip(a, b))
    return all(
        math.isclose(float(x), float(y), rel_tol=0.0, abs_tol=atol)
        for x, y in zip(a, b)
    )


def _git_revision() -> str | None:
    try:
        p = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True
        )
        return p.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return os.environ.get("GITHUB_SHA")


def _tool_version(command: list[str]) -> str | None:
    try:
        p = subprocess.run(command, capture_output=True, text=True, check=True)
        return (p.stdout or p.stderr).splitlines()[0].strip()
    except (OSError, subprocess.CalledProcessError, IndexError):
        return None


def _inspect_source(source_id: str, entry: dict[str, Any], source_root: Path) -> dict[str, Any]:
    source_dir = source_root / Path(entry["raw_directory"]).name
    result: dict[str, Any] = {
        "source_id": source_id,
        "source_directory": str(source_dir),
        "source_identity_status": "NOT_RUN",
        "raw_source_parse_status": "NOT_RUN",
        "reader_status": "NOT_RUN",
        "structural_status": "NOT_RUN",
        "compose_roundtrip_status": "NOT_RUN",
        "verifier_execution_status": "NOT_RUN",
        "scientific_status": "NOT_EVALUATED",
        "root_search": "NOT_RUN",
        "tov_cowling": "NOT_RUN",
    }

    before: dict[str, dict[str, Any]] = {}
    try:
        for filename, metadata in entry.get("files", {}).items():
            if metadata is None:
                continue
            path = source_dir / filename
            if not path.is_file():
                raise FileNotFoundError(f"missing source file: {path}")
            fingerprint = _fingerprint(path)
            expected = metadata.get("git_blob_sha")
            if expected and fingerprint["git_blob_sha1"] != expected:
                raise ValueError(
                    f"Git blob identity mismatch for {path}: expected {expected}, "
                    f"observed {fingerprint['git_blob_sha1']}"
                )
            before[filename] = fingerprint
        result["source_identity"] = {"files_before": before}
        result["source_identity_status"] = "PASS"
    except Exception as exc:
        result["source_identity_status"] = "ERROR"
        result["source_identity_error"] = f"{type(exc).__name__}: {exc}"
        return result

    raw_ok = False
    try:
        il, nb = _read_raw_density(source_dir / "eos.nb")
        nadds, q_rows, extras, thermo_header = _read_raw_thermo(source_dir / "eos.thermo")
        compo_path = source_dir / "eos.compo"
        if entry.get("eos_compo_present", True) and compo_path.is_file():
            raw_phase, raw_composition = _read_raw_composition(compo_path)
        else:
            raw_phase, raw_composition = [], {}
        result["raw_source_parse"] = {
            "il": il, "nb": nb, "nadd": nadds, "q": q_rows,
            "q_extra": extras, "thermo_header": thermo_header,
            "iphase": raw_phase, "composition": raw_composition,
        }
        result["raw_source_parse_status"] = "PASS"
        raw_ok = True
    except Exception as exc:
        result["raw_source_parse_status"] = "ERROR"
        result["raw_source_parse_error"] = f"{type(exc).__name__}: {exc}"

    sys.path.insert(0, str(ROOT / "python_frontend"))
    try:
        from compose_eos.api import read_compose, write_compose, check_consistency
        data = read_compose(source_dir)
        result["reader_status"] = "PASS"
        result["reader"] = {
            "n": int(data.n),
            "temperature_mev": float(data.temperature_mev),
            "mn_mev": float(data.mn),
            "mp_mev": float(data.mp),
            "nadd": int(data.q_extra_count),
            "composition_file_present": bool(data.composition_file_present),
            "composition_available": bool(data.composition_available),
            "beta_equilibrated": bool(data.beta_equilibrated),
        }
    except Exception as exc:
        result["reader_status"] = "ERROR"
        result["reader_error"] = f"{type(exc).__name__}: {exc}"
        result["source_identity"]["files_after"] = {
            name: _fingerprint(source_dir / name) for name in before
        }
        result["source_identity"]["unchanged_during_probe"] = (
            result["source_identity"]["files_before"] == result["source_identity"]["files_after"]
        )
        result["structural_status"] = "NOT_RUN"
        return result

    try:
        finite = _finite_arrays(data)
        expected_masses = entry["source_masses_mev"]
        expected_nadd_values = [int(x) for x in entry.get("Nadd_values", [])]
        expected_nadd = expected_nadd_values[0] if len(set(expected_nadd_values)) == 1 else None
        expected_phase_hist = {str(k): int(v) for k, v in entry.get("Iphase_histogram", {}).items()}
        expected_extra_hist = {
            str(k): int(v) for k, v in entry.get("additional_column_value_counts", {}).items()
        }
        observed_phase_hist = _histogram(data.Iphase)
        observed_extra_hist = _extra_histogram(data.q_extra) if data.q_extra_count else {}
        observed_boundaries = _phase_boundaries(data.Iphase)

        raw = result.get("raw_source_parse", {})
        if raw_ok:
            rowwise = {
                "nb": _rowwise_equal(raw["nb"], data.nb),
                "Q1_to_Q7": len(raw["q"]) == data.n and all(
                    _rowwise_equal([row[k] for row in raw["q"]], getattr(data, f"Q{k+1}"))
                    for k in range(7)
                ),
                "q_extra": _rowwise_equal(raw["q_extra"], data.q_extra, atol=0.0),
                "Iphase": _rowwise_equal(raw["iphase"], data.Iphase, integer=True),
            }
            for name, values in raw["composition"].items():
                rowwise[name] = _rowwise_equal(values, getattr(data, name))
            nadd_consistent = (
                len(raw["nadd"]) == data.n and len(set(raw["nadd"])) == 1
                and expected_nadd is not None and raw["nadd"][0] == expected_nadd
            )
            il_match = int(raw["il"]) == int(data.Il)
            thermo_header_matches = (
                math.isclose(raw["thermo_header"][0], float(expected_masses["m_n"]), rel_tol=0.0, abs_tol=1e-8)
                and math.isclose(raw["thermo_header"][1], float(expected_masses["m_p"]), rel_tol=0.0, abs_tol=1e-8)
                and int(raw["thermo_header"][2]) == int(data.Il)
            )
        else:
            rowwise = {}
            nadd_consistent = il_match = thermo_header_matches = False

        roundtrip_checks: dict[str, bool] = {}
        roundtrip_error = None
        try:
            with tempfile.TemporaryDirectory(prefix="sharp-maxwell-compose-roundtrip-") as td:
                roundtrip_dir = Path(td) / "roundtrip"
                write_compose(roundtrip_dir, data)
                reread = read_compose(roundtrip_dir)
                for name in COMPOSE_ARRAYS:
                    roundtrip_checks[name] = _rowwise_equal(
                        getattr(data, name), getattr(reread, name), integer=(name == "Iphase"),
                        atol=5e-13 if name != "q_extra" else 0.0
                    )
                roundtrip_checks["Nadd"] = int(data.q_extra_count) == int(reread.q_extra_count)
                roundtrip_checks["composition_file_present"] = bool(reread.composition_file_present)
                roundtrip_checks["composition_available"] = bool(reread.composition_available)
        except Exception as exc:
            roundtrip_error = f"{type(exc).__name__}: {exc}"
        result["compose_roundtrip_checks"] = roundtrip_checks
        result["compose_roundtrip_status"] = (
            "PASS" if roundtrip_checks and all(roundtrip_checks.values()) else
            "ERROR" if roundtrip_error else "FAIL"
        )
        if roundtrip_error:
            result["compose_roundtrip_error"] = roundtrip_error

        checks = {
            "raw_source_parser_succeeded": raw_ok,
            "row_count_matches_manifest": int(data.n) == int(entry["rows"]),
            "il_matches_grid_header": il_match,
            "thermo_masses_match_manifest": thermo_header_matches,
            "nb_preserved_rowwise": rowwise.get("nb", False),
            "Q1_to_Q7_preserved_rowwise": rowwise.get("Q1_to_Q7", False),
            "q_extra_preserved_rowwise": rowwise.get("q_extra", False),
            "Nadd_constant_and_matches_manifest": nadd_consistent and int(data.q_extra_count) == expected_nadd,
            "temperature_is_zero": abs(float(data.temperature_mev)) <= 1e-12,
            "native_density_positive_and_strictly_increasing": (
                all(float(x) > 0.0 for x in data.nb)
                and all(float(b) > float(a) for a, b in zip(data.nb, data.nb[1:]))
            ),
            "all_native_arrays_finite_and_aligned": all(finite.values()),
            "composition_file_presence_matches_manifest": (
                bool(data.composition_file_present) == bool(entry.get("eos_compo_present", False))
            ),
            "composition_available": bool(data.composition_available),
            "composition_rows_match_manifest": (
                entry.get("eos_compo_rows") is None or int(entry["eos_compo_rows"]) == int(data.n)
            ),
            "composition_fields_preserved_rowwise": all(
                rowwise.get(name, False) for name in ("Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh")
            ),
            "Iphase_histogram_preserved": observed_phase_hist == expected_phase_hist,
            "Iphase_preserved_rowwise": rowwise.get("Iphase", False),
            "Iphase_boundaries_preserved": observed_boundaries == _expected_boundary(entry),
            "q_extra_histogram_preserved": observed_extra_hist == expected_extra_hist,
            "native_compose_write_read_roundtrip": (
                bool(roundtrip_checks) and all(roundtrip_checks.values())
            ),
        }
        result["structural_checks"] = checks
        result["array_finiteness"] = finite
        result["preserved_metadata"] = {
            "observed_iphase_histogram": observed_phase_hist,
            "expected_iphase_histogram": expected_phase_hist,
            "observed_iphase_boundaries": observed_boundaries,
            "expected_iphase_boundaries": _expected_boundary(entry),
            "observed_q_extra_histogram": observed_extra_hist,
            "expected_q_extra_histogram": expected_extra_hist,
            "raw_to_native_rowwise": rowwise,
            "roundtrip": roundtrip_checks,
        }
        result["structural_status"] = "PASS" if all(checks.values()) else "FAIL"
    except Exception as exc:
        result["structural_status"] = "ERROR"
        result["structural_error"] = f"{type(exc).__name__}: {exc}"

    try:
        report = check_consistency(data)
        result["verifier_execution_status"] = "PASS"
        result["native_verifier"] = report
        overall = report.get("overall", {})
        result["scientific_status"] = "PASS" if overall.get("pass_") is True else "FAIL"
        result["scientific_failed_checks"] = list(overall.get("failed_checks", []))
        result["scientific_scope_note"] = (
            "Native C1-C6/C10-C11 knot-level gate; C7-C9 are finite-difference diagnostics. "
            "This is not whole-Hermite-interval admission."
        )
    except Exception as exc:
        result["verifier_execution_status"] = "ERROR"
        result["verifier_error"] = f"{type(exc).__name__}: {exc}"

    after = {name: _fingerprint(source_dir / name) for name in before}
    unchanged = before == after
    result["source_identity"]["files_after"] = after
    result["source_identity"]["unchanged_during_probe"] = unchanged
    result.setdefault("structural_checks", {})["raw_source_bytes_unchanged"] = unchanged
    if not unchanged:
        result["structural_status"] = "FAIL"
    elif result.get("structural_status") == "PASS" and all(result["structural_checks"].values()):
        result["structural_status"] = "PASS"
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--source-id", action="append", dest="source_ids")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--strict-physics", action="store_true")
    args = parser.parse_args(argv)

    manifest_path = args.manifest.expanduser().resolve()
    source_root = args.source_root.expanduser().resolve()
    selected = args.source_ids or list(DEFAULT_SOURCE_IDS)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    catalog = manifest["source_catalog"]["crust_models"]
    receipt: dict[str, Any] = {
        "schema": "sharp_maxwell_c5_native_source_probe_receipt.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "research_loop_stage": "C5",
        "scope": "read-only native reader, Compose round-trip, structural preservation and C1-C11 vertical slice",
        "manifest": {
            "path": str(manifest_path),
            "manifest_id": manifest.get("manifest_id"),
            "identity": _fingerprint(manifest_path),
            "git_revision": _git_revision(),
        },
        "runner": {"path": str(Path(__file__).resolve()), "sha256": _fingerprint(Path(__file__).resolve())["sha256"]},
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "gfortran": _tool_version(["gfortran", "--version"]),
        },
        "selected_source_ids": selected,
        "sources": [],
        "runner_status": "PASS",
        "scientific_admission": "NOT_DECIDED_BY_C5_PROBE",
        "root_search": "NOT_RUN",
        "tov_cowling": "NOT_RUN",
        "scope_limitations": [
            "No source repair, normalization, projection, interpolation, matching or raw-source rewrite.",
            "Iphase and Nadd provider semantics remain unresolved beyond format-level preservation.",
            "C7-C9 are diagnostics; no whole-Hermite-interval extrema or full physical admission.",
            "No density reconstruction, rest-mass transformation, pair compatibility or Maxwell root search.",
            "No TOV, stable-branch, Seidov or Cowling calculation.",
        ],
    }

    for source_id in selected:
        if source_id not in catalog:
            item = {
                "source_id": source_id,
                "source_identity_status": "NOT_RUN",
                "reader_status": "NOT_RUN",
                "structural_status": "FAIL",
                "compose_roundtrip_status": "NOT_RUN",
                "verifier_execution_status": "NOT_RUN",
                "scientific_status": "NOT_EVALUATED",
                "error": "unknown source ID in manifest",
            }
        else:
            try:
                item = _inspect_source(source_id, catalog[source_id], source_root)
            except Exception as exc:
                item = {
                    "source_id": source_id,
                    "reader_status": "ERROR",
                    "structural_status": "NOT_RUN",
                    "compose_roundtrip_status": "NOT_RUN",
                    "verifier_execution_status": "NOT_RUN",
                    "scientific_status": "NOT_EVALUATED",
                    "root_search": "NOT_RUN",
                    "tov_cowling": "NOT_RUN",
                    "error": f"{type(exc).__name__}: {exc}",
                }
        receipt["sources"].append(item)
        required = (
            item.get("source_identity_status") == "PASS"
            and item.get("reader_status") == "PASS"
            and item.get("structural_status") == "PASS"
            and item.get("compose_roundtrip_status") == "PASS"
            and item.get("verifier_execution_status") == "PASS"
        )
        if not required:
            receipt["runner_status"] = "FAIL"

    if receipt["runner_status"] == "PASS":
        try:
            sys.path.insert(0, str(ROOT / "python_frontend"))
            from compose_eos.native import executable
            exe = executable()
            receipt["native_executable"] = {
                "path": str(exe),
                "sha256": _fingerprint(exe)["sha256"],
            }
        except Exception as exc:
            receipt["native_executable_error"] = f"{type(exc).__name__}: {exc}"
            receipt["runner_status"] = "FAIL"

    rendered = json.dumps(_json_safe(receipt), indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        output = args.output.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")

    if receipt["runner_status"] != "PASS":
        return 2
    if args.strict_physics and any(x.get("scientific_status") == "FAIL" for x in receipt["sources"]):
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
