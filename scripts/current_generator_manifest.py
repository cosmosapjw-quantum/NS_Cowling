"""Explicit, non-executing intake of current generator matching sources.

This module inventories a 13 x 5 x 6 request surface. It does not promote
historical matches or run the generator or matcher.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

BASE_COMMIT = "23b1f2099dc15470ac63a100e6b3aae0c397a36b"
SCHEMA_VERSION = "ns-cowling.current-generator-intake.v1"
CRUSTS = ("SPG_M1", "SPG_M2", "SPG_M3", "SPG_M4", "SPG_M5",
          "GPPVA_DD2", "GPPVA_DDME2", "GPPVA_FSU2", "GPPVA_FSU2H",
          "GPPVA_FSU2R", "GPPVA_NL3wrL55", "GPPVA_TM1e", "GPPVA_TW")
CORES = ("APR_APR", "BBB_BHF2", "GM_GM1", "RG_KDE0v1", "RG_SLY4")
METHODS = ("cubic_F", "cubic_P", "quintic_F", "quintic_P", "css", "maxwell")
SOURCE_FIELDS = ("id", "model", "family", "variant", "representation", "compose_path",
                 "request_path", "request_sha256", "source_dependencies_sha256",
                 "receipt_paths_sha256", "artifact_sha256", "support_path", "support_sha256",
                 "density_segments", "phase_segments", "composition_policy",
                 "generation_executable_identity")
CELL_FIELDS = ("id", "crust_source_id", "core_source_id", "method", "n_bridge",
               "requested", "eligible", "attempted", "admitted", "status",
               "reason_code", "evidence_refs", "reuse_mode")
COMPOSE_FILES = ("eos.nb", "eos.t", "eos.yq", "eos.thermo")
_PINNED = {
    "crust:SPG_M1": "data/research_matching_20261117/spg/native_129_repaired/SPG_M1",
    "core:GM_GM1": "data/research_matching_20261117/integration/real_pairs_v2/core/GM_GM1",
}
_SEALED_RECEIPTS = {
    "data/research_matching_20261117/spg/native_129_repaired/RUN_RECEIPT.json": "955dec02748f124a72fab3d0010e379c702e1475aacec25ca3cfba7d3b451ea4",
    "data/research_matching_20261117/spg/native_129_repaired/SPG_M1/RECEIPT.json": "5c1d5f70033697417930f40aaf038ab16750b66d939cb7daf5381f46e913e20b",
    "data/research_matching_20261117/spg/native_129_repaired/SPG_M1/generated/RESULT.json": "315b41e499da5350710c3418a4f21d97251c5e6cfe0f9a8b784905553f8a2bff",
    "data/research_matching_20261117/integration/real_pairs_v2/core/GM_GM1/DISPATCH.json": "f7d9774a26524479f0d6076b94af2c156abc8cb12912e8b66e0787c04def9447",
    "data/research_matching_20261117/integration/real_pairs_v2/core/GM_GM1/generated/RESULT.json": "0108fb20c83d13ec20eecfee57e71d6ec96b0a24e126cab3ab28955643e9752c",
}


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def _safe(root: Path, path: str) -> Path:
    p = Path(path)
    if p.is_absolute() or ".." in p.parts or not path:
        raise ValueError(f"path is not repository relative: {path}")
    return root / p


def _observed_segments(compose_dir: Path) -> tuple[list[dict], list[dict]]:
    """Read row metadata only; these ranges assert no phase-boundary law."""
    nb_lines = [line.split() for line in (compose_dir / "eos.nb").read_text().splitlines()
                if line.strip()]
    n = int(nb_lines[1][0])
    densities = [float(row[0]) for row in nb_lines[2:]]
    if n < 1 or len(densities) != n or any(not value > 0 for value in densities):
        raise ValueError("invalid sealed eos.nb row inventory")
    phase_by_row = {}
    for line in (compose_dir / "eos.compo").read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split()
        index = int(fields[1]) - 1
        if index < 0 or index >= n or index in phase_by_row:
            raise ValueError("invalid or duplicate CompOSE phase row")
        phase_by_row[index] = int(fields[3])
    if set(phase_by_row) != set(range(n)):
        raise ValueError("missing CompOSE phase rows")
    phases = []
    for index in range(n):
        label = phase_by_row[index]
        if phases and phases[-1]["iphase"] == label:
            phases[-1]["row_end"] = index
        else:
            phases.append(dict(iphase=label, row_start=index, row_end=index))
    density = [dict(row_start=0, row_end=n - 1,
                    n_b_min=min(densities), n_b_max=max(densities))]
    return density, phases


def _empty_source(family: str, model: str) -> dict:
    sid = f"{family}:{model}"
    source = dict.fromkeys(SOURCE_FIELDS)
    source.update(id=sid, model=model, family=family, variant=None,
                  representation=None, compose_path=None, request_path=None,
                  request_sha256=None, source_dependencies_sha256={},
                  receipt_paths_sha256={}, artifact_sha256={}, support_path=None,
                  support_sha256=None, density_segments=[], phase_segments=[],
                  composition_policy=None, generation_executable_identity=None)
    return source


def _source(root: Path, family: str, model: str) -> dict:
    source = _empty_source(family, model)
    sid = source["id"]
    if sid not in _PINNED:
        return source
    parent = _PINNED[sid]
    receipt_name = "RECEIPT.json" if family == "crust" else "DISPATCH.json"
    receipt_path = f"{parent}/{receipt_name}"
    receipt = json.loads(_safe(root, receipt_path).read_text())
    native = json.loads(_safe(root, f"{parent}/generated/RESULT.json").read_text())
    density_segments, phase_segments = _observed_segments(_safe(root, f"{parent}/generated"))
    if density_segments[0]["row_end"] + 1 != native["rows"]:
        raise ValueError("native RESULT row count differs from sealed CompOSE table")
    output = receipt["output_sha256"] if family == "crust" else receipt["evidence"]["output_sha256"]
    # The outer seal's output inventory is authoritative, including optional
    # diagnostics that were present when the receipt was written.
    artifacts = {name.removeprefix("generated/"): value for name, value in output.items()
                 if name.startswith("generated/")} if family == "crust" else dict(output)
    deps = {}
    if family == "core":
        deps = {("python_frontend/" if "/python_frontend/" in p else "fortran_eos/")
                + p.split("/python_frontend/" if "/python_frontend/" in p else "/fortran_eos/", 1)[1]: h
                for p, h in receipt["source_sha256"].items()}
    else:
        surface = receipt.get("surface_input", {})
        if surface.get("path") and surface.get("sha256"):
            deps[surface["path"]] = surface["sha256"]
    if family == "crust":
        run_receipt_path = "data/research_matching_20261117/spg/native_129_repaired/RUN_RECEIPT.json"
        if _sha(_safe(root, run_receipt_path)) != _SEALED_RECEIPTS[run_receipt_path]:
            raise ValueError("SPG RUN_RECEIPT byte identity mismatch")
        run_receipt = json.loads(_safe(root, run_receipt_path).read_text())
        binary_hashes = run_receipt["files"]
        generation_identity = {
            "path": "build/spg_finish_repaired/compose_eos_generate",
            "sha256": next(value for path, value in binary_hashes.items()
                           if path.endswith("/build/spg_finish_repaired/compose_eos_generate")),
            "readback_path": "build/spg_finish_repaired/compose_eos_fortran",
            "readback_sha256": next(value for path, value in binary_hashes.items()
                                    if path.endswith("/build/spg_finish_repaired/compose_eos_fortran")),
            "identity_scope": "historical_RUN_RECEIPT",
            "receipt_path": run_receipt_path,
            "receipt_sha256": _SEALED_RECEIPTS[run_receipt_path],
        }
    else:
        generation_identity = {
            "status": "NOT_INDEPENDENTLY_OBSERVED",
            "path": None, "sha256": None,
            "launcher_path": "fortran_eos/run-generator.sh",
            "launcher_sha256": deps.get("fortran_eos/run-generator.sh"),
        }
    source.update(variant=receipt["request"].get("variant"),
                  representation=receipt["request"].get("contract", {}).get("representation"),
                  compose_path=f"{parent}/generated", request_path=f"{parent}/request.nml",
                  request_sha256=receipt["request_sha256"],
                  source_dependencies_sha256=deps,
                  receipt_paths_sha256={receipt_path: _sha(_safe(root, receipt_path)),
                                        f"{parent}/generated/RESULT.json": _sha(_safe(root, f"{parent}/generated/RESULT.json"))},
                  artifact_sha256=artifacts,
                  support_path=f"{parent}/generated/SUPPORT.json",
                  support_sha256=artifacts["SUPPORT.json"],
                  density_segments=density_segments,
                  phase_segments=phase_segments,
                  composition_policy=receipt["request"].get("contract", {}).get("composition"),
                  generation_executable_identity=generation_identity)
    return source


def source_dependency_mismatches(root: Path, source: dict) -> list[dict]:
    mismatches = []
    for path, expected in source["source_dependencies_sha256"].items():
        actual = _sha(_safe(root, path)) if _safe(root, path).is_file() else None
        if actual != expected:
            mismatches.append(dict(path=path, expected_sha256=expected,
                                   observed_sha256=actual))
    return mismatches


def build_manifest(root: Path) -> dict:
    """Build all coverage requests; only one cell is eligible for fresh work."""
    sources = []
    issues = {}
    for family, models in (("crust", CRUSTS), ("core", CORES)):
        for model in models:
            try:
                sources.append(_source(root, family, model))
            except (OSError, KeyError, ValueError, TypeError, StopIteration, IndexError,
                    json.JSONDecodeError) as exc:
                sources.append(_empty_source(family, model))
                code = ("HOLD_SOURCE_INTAKE_MISSING" if isinstance(exc, FileNotFoundError)
                        else "HOLD_SOURCE_INTAKE_MALFORMED" if isinstance(exc, (KeyError, TypeError, IndexError, json.JSONDecodeError))
                        else "BYTE_IDENTITY_MISMATCH")
                issues[f"{family}:{model}"] = dict(source_id=f"{family}:{model}",
                                                     reason_code=code,
                                                     detail=f"{type(exc).__name__}: {exc}")
    pinned = {source["id"]: source for source in sources}
    for sid in _PINNED:
        if sid in issues:
            continue
        try:
            validate_source(root, pinned[sid])
        except (OSError, KeyError, ValueError, TypeError, StopIteration, IndexError,
                json.JSONDecodeError) as exc:
            code = ("HOLD_SOURCE_INTAKE_MISSING" if isinstance(exc, FileNotFoundError)
                    else "HOLD_SOURCE_INTAKE_MALFORMED" if isinstance(exc, (KeyError, TypeError, IndexError, json.JSONDecodeError))
                    else "BYTE_IDENTITY_MISMATCH")
            issues[sid] = dict(source_id=sid, reason_code=code,
                               detail=f"{type(exc).__name__}: {exc}")
    drift = any(source_dependency_mismatches(root, pinned[sid]) for sid in _PINNED if sid not in issues)
    sealed_hold = bool(issues)
    cells = []
    for crust in CRUSTS:
        for core in CORES:
            for method in METHODS:
                minimum = (crust, core, method) == ("SPG_M1", "GM_GM1", "css")
                implemented = method in ("css", "maxwell")
                cells.append(dict(id=f"{crust}__{core}__{method}",
                                  crust_source_id=f"crust:{crust}", core_source_id=f"core:{core}",
                                  method=method, n_bridge=40 if minimum else None,
                                  requested=True, eligible=minimum and not sealed_hold,
                                  attempted=False, admitted=False,
                                  status="READY" if minimum and not sealed_hold else "HOLD",
                                  reason_code=(next(iter(issues.values()))["reason_code"] if sealed_hold else
                                               "HISTORICAL_GENERATOR_DEPENDENCY_VERSION_DIFFERENCE" if drift
                                               else "FRESH_MINIMUM_PENDING") if minimum else
                                  ("HOLD_SOURCE_INTAKE_NOT_PINNED" if implemented else
                                   "HOLD_METHOD_INTAKE_NOT_IMPLEMENTED"),
                                  evidence_refs=list(issues.values()) if sealed_hold else [],
                                  reuse_mode="none"))
    return dict(schema_version=SCHEMA_VERSION, base_commit=BASE_COMMIT,
                sources=sources, cells=cells)


def validate_source(root: Path, source: dict) -> None:
    """Validate all bytes bound by the actual outer seal and native RESULT."""
    if set(source) != set(SOURCE_FIELDS) or source["id"] not in _PINNED:
        raise ValueError("source is not a pinned current-generator source")
    parent = _PINNED[source["id"]]
    if source["compose_path"] != f"{parent}/generated" or source["request_path"] != f"{parent}/request.nml":
        raise ValueError("source paths differ from pinned intake")
    receipt_name = "RECEIPT.json" if source["family"] == "crust" else "DISPATCH.json"
    rp = f"{parent}/{receipt_name}"
    if set(source["receipt_paths_sha256"]) != {rp, f"{parent}/generated/RESULT.json"}:
        raise ValueError("outer/native receipt set differs")
    for path, digest in source["receipt_paths_sha256"].items():
        if digest != _SEALED_RECEIPTS[path] or _sha(_safe(root, path)) != digest:
            raise ValueError(f"receipt hash mismatch: {path}")
    receipt = json.loads(_safe(root, rp).read_text())
    native = json.loads(_safe(root, f"{parent}/generated/RESULT.json").read_text())
    if receipt["status"] != "REQUESTED_ROWS_NATIVE_VERIFIED" or receipt["returncode"] != 0:
        raise ValueError("outer seal does not verify requested rows")
    if source["family"] == "core" and receipt["child_returncode"] != 0:
        raise ValueError("dispatcher child failed")
    if (native["status"] != "REQUESTED_ROWS_NATIVE_VERIFIED" or native["backend"] != "fortran"
            or native["native_readback"] is not True or native["model"] != source["model"]
            or native["family"] != ("spg" if source["family"] == "crust" else "core")
            or receipt["request"]["model"] != source["model"]):
        raise ValueError("native RESULT identity or readback invalid")
    if source["request_sha256"] != receipt["request_sha256"] or _sha(_safe(root, source["request_path"])) != source["request_sha256"]:
        raise ValueError("request hash mismatch")
    output = receipt["output_sha256"] if source["family"] == "crust" else receipt["evidence"]["output_sha256"]
    sealed = {k.removeprefix("generated/"): v for k, v in output.items() if k.startswith("generated/")} if source["family"] == "crust" else output
    if source["artifact_sha256"] != sealed or not set(COMPOSE_FILES).issubset(sealed):
        raise ValueError("artifact inventory differs from outer seal")
    for name, digest in sealed.items():
        if _sha(_safe(root, f"{parent}/generated/{name}")) != digest:
            raise ValueError(f"sealed artifact missing or changed: {name}")
    if source["support_path"] != f"{parent}/generated/SUPPORT.json" or source["support_sha256"] != sealed.get("SUPPORT.json"):
        raise ValueError("SUPPORT identity mismatch")
    if "eos.compo" not in sealed:
        raise ValueError("phase source is not sealed")
    density_segments, phase_segments = _observed_segments(_safe(root, source["compose_path"]))
    if (source["density_segments"] != density_segments or source["phase_segments"] != phase_segments
            or density_segments[0]["row_end"] + 1 != native["rows"]):
        raise ValueError("sealed density or phase row segments differ")
    support = json.loads(_safe(root, source["support_path"]).read_text())
    if support.get("cross_phase_interpolation") != "UNAPPROVED" or support.get("cross_segment_interpolation") != "UNAPPROVED":
        raise ValueError("SUPPORT restrictions changed")
    if source["family"] == "core":
        expected_dependencies = {("python_frontend/" if "/python_frontend/" in p else "fortran_eos/")
                                 + p.split("/python_frontend/" if "/python_frontend/" in p else "/fortran_eos/", 1)[1]: h
                                 for p, h in receipt["source_sha256"].items()}
    else:
        surface = receipt.get("surface_input", {})
        expected_dependencies = {surface["path"]: surface["sha256"]} if surface.get("path") else {}
    if source["source_dependencies_sha256"] != expected_dependencies:
        raise ValueError("historical dependency pins differ from receipt")
    if source["generation_executable_identity"] != _source(root, source["family"], source["model"])["generation_executable_identity"]:
        raise ValueError("generation executable identity differs from sealed receipt")


def validate_manifest(manifest: dict) -> None:
    if set(manifest) != {"schema_version", "base_commit", "sources", "cells"} or manifest["schema_version"] != SCHEMA_VERSION or manifest["base_commit"] != BASE_COMMIT:
        raise ValueError("manifest schema/base mismatch")
    sources, cells = manifest["sources"], manifest["cells"]
    expected_sources = {f"crust:{model}" for model in CRUSTS} | {f"core:{model}" for model in CORES}
    if (len(sources) != 18 or {s["id"] for s in sources} != expected_sources
            or any(set(s) != set(SOURCE_FIELDS) for s in sources)):
        raise ValueError("source grid invalid")
    expected = {f"{crust}__{core}__{method}" for crust in CRUSTS for core in CORES for method in METHODS}
    if len(cells) != 390 or {c["id"] for c in cells} != expected or any(set(c) != set(CELL_FIELDS) for c in cells):
        raise ValueError("390 cell grid invalid")
    if any(cell["requested"] is not True for cell in cells):
        raise ValueError("all 390 coverage cells must be requested")
    for cell in cells:
        crust, core, method = cell["id"].split("__")
        if (cell["crust_source_id"] != f"crust:{crust}" or cell["core_source_id"] != f"core:{core}"
                or cell["method"] != method or cell["admitted"] or cell["attempted"]
                or (cell["reuse_mode"] != "none" and (cell["eligible"] or cell["admitted"]))):
            raise ValueError("cell identity or unsupported promotion")
        if method not in ("css", "maxwell") and (cell["status"] != "HOLD" or cell["reason_code"] != "HOLD_METHOD_INTAKE_NOT_IMPLEMENTED"):
            raise ValueError("smooth method intake must remain on HOLD")
        if cell["eligible"] and cell["id"] != "SPG_M1__GM_GM1__css":
            raise ValueError("only minimum fresh cell may be eligible")
        if cell["n_bridge"] != (40 if cell["id"] == "SPG_M1__GM_GM1__css" else None):
            raise ValueError("bridge execution parameter exceeds single-cell cap")
