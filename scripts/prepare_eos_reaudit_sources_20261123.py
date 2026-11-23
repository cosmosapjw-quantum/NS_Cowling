#!/usr/bin/env python3
"""Freeze exact source views for the 2026-11-23 native all-pair campaign.

Python only verifies provenance, dispatches native readers/importers and copies
selected rows. No EOS interpolation, reconstruction, causal screening, energy
shift or numerical rejection prefilter is implemented here. Phase components
are prescribed contiguous Iphase runs; failed physical intervals never create
new components. Private input and all private table derivatives stay outside Git.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos.api import read_compose, _thermo_only_view
from compose_eos.eos_data import EOSData
from compose_eos.maxwell_v2 import _segment_eosdata

METHODS = ["maxwell", "css", "cubic_F", "cubic_P", "quintic_F", "quintic_P",
           "first_law_constrained_hermite_pressure_bridge"]
OLD_CATALOG = ROOT / "data/all_origin_maxwell_css_20261118/SOURCE_CATALOG.json"
OLD_CATALOG_SHA = "d8a2fe3ad429b4dcde2bcec7f1e33c7caa2bd31dc5bf1cd72dcbdfafb842873f"
DEFAULT_NUDA = Path("/tmp/ns_eos_reaudit_20261123/generator/NUDA_NATIVE_IMPORT_LAYOUTS.json")
NUDA_MANIFEST_SHA = "320958a27f9f5f9a86b8347c92d6f8eb80ecc1a973e1a6709c40d2a09ea61320"
PRIVATE_ZIP_SHA = "ca42416dfdab5a352139e14ce8a1b86cb6e2a516fcc9c0e07f48350ccd09271c"
PROVIDER_ROOT = ROOT
CONTROLS = ROOT / "data/cowling_maxwell_native_20261119/provider_original_unified/ORIGINAL_CONTROLS.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sha(path: Path, expected: str) -> str:
    actual = sha(path)
    if actual != expected:
        raise ValueError(f"SOURCE_HASH_MISMATCH: {path}: {actual} != {expected}")
    return actual


def ref(path: Path) -> str:
    p = path.resolve()
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def absolute(path: str | Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def save(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    tmp.replace(path)


def hashes(directory: Path, names: tuple[str, ...] | None = None) -> dict:
    files = [directory / n for n in names] if names else sorted(directory.glob("eos.*"))
    return {p.name: sha(p) for p in files if p.is_file()}


def verify_old(source: dict) -> dict:
    verified = {}
    for key, path_key in (("source_table_sha256", "source_path"),
                          ("internal_sha256", "internal_path"),
                          ("compose_sha256", "compose_path")):
        expected = source.get(key)
        if not isinstance(expected, dict):
            continue
        path = absolute(source[path_key])
        verified[path_key] = {name: check_sha(path / name, digest)
                              for name, digest in expected.items()}
    if source.get("source_receipt_path") and source.get("source_receipt_sha256"):
        verified["receipt"] = check_sha(absolute(source["source_receipt_path"]),
                                         source["source_receipt_sha256"])
    return verified


def base_record(source: dict) -> dict:
    return dict(source_id=source["source_id"], role=source["role"], model=source["model"],
                origin=source["source_id"].split(":", 1)[0],
                provenance_class=source.get("provenance_class"),
                available=False, status=source.get("status", "UNPREPARED"),
                private=source["source_id"].startswith("user_zip:"),
                compose_path=None, internal_path=None, source_sha256={},
                selection=source.get("selection") or {}, components=[], matching_views=[],
                whole_selected_matching_allowed=False,
                fullcoveragepolicy="WHOLE_SELECTED_RETAINED_CHAIN requires all retained source rows and inherited boundaries; a finite component is not whole-source admission",
                source_catalog_entry=source,
                legacy_diagnostic_components=source.get("components", []),
                old_status=source.get("status"),
                full_outer_inner=source.get("full_outer_inner", False),
                outer_present=source.get("outer_present", False),
                inner_present=source.get("inner_present", False))


def phase_runs(data: EOSData) -> list[tuple[int, int, int]]:
    # Copy source-declared labels. No physical-value predicate enters the split.
    labels = [int(x) for x in data.Iphase]
    runs = []
    start = 0
    for stop in range(1, len(labels) + 1):
        if stop == len(labels) or labels[stop] != labels[start]:
            runs.append((labels[start], start, stop))
            start = stop
    return runs


def freeze_view(record: dict, data: EOSData, destination: Path) -> dict:
    if data.n == 0:
        raise ValueError("SOURCE_EMPTY_SELECTED_DOMAIN")
    destination.mkdir(parents=True, exist_ok=True)
    # Only the abundance inventory is explicitly omitted. Q1-Q7, all already
    # present auxiliary quantities, labels, density knots and source axes survive.
    view = _thermo_only_view(data)
    full = destination / "full_selected" / "internal"
    view.write_internal_directory(full)
    restored = EOSData.from_internal_directory(full)
    import numpy as np
    for field in ("nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "q_extra", "Iphase", "q_extra_full"):
        if not np.array_equal(getattr(view, field), getattr(restored, field)):
            raise ValueError(f"SOURCE_ROW_MARSHAL_CHANGED: {field}")
    density = [float(x) for x in view.nb]
    axis_file = destination / "SOURCE_DENSITIES.json"
    save(axis_file, {"source_id": record["source_id"], "nb_fm3": density})
    record.update(available=True, status="SOURCE_SELECTED_VIEW_PREPARED", internal_path=ref(full),
                  full_selected_internal_path=ref(full),
                  full_selected_rows=view.n, selected_density_fm3=[density[0], density[-1]],
                  source_densities=None if record["private"] else density,
                  source_densities_path=ref(axis_file),
                  full_selected_sha256=hashes(full, ("eos_data.csv", "eos_data.nml")),
                  composition_policy="EXPLICIT_THERMO_ONLY_ABUNDANCE_OMISSION",
                  q_extra_count=int(view.q_extra_count), nadd_max=int(view.nadd_max),
                  numerical_transformation="NONE: exact row selection and lossless serialization only")
    matching = view
    matching_full = full
    auxiliary_projection = None
    if view.q_extra_count or view.nadd_max:
        # Explicit output-schema projection, not a numerical EOS transform.
        # Original columns remain in full_selected/internal and this sidecar.
        sidecar = destination / "AUXILIARY_PROJECTION.json"
        save(sidecar, dict(schema="ns_cowling.exact_auxiliary_omission.v1",
            source_id=record["source_id"], source_internal_path=ref(full),
            source_internal_sha256=record["full_selected_sha256"],
            source_q_extra_count=int(view.q_extra_count), source_nadd_max=int(view.nadd_max),
            q_extra=view.q_extra.tolist(), q_extra_full=view.q_extra_full.tolist(),
            q_extra_full_shape=list(view.q_extra_full.shape),
            semantics="auxiliary_not_consumed_by_thermo_matcher; source meaning not inferred",
            operation="OMIT_AUXILIARY_COLUMNS_FROM_MATCHING_VIEW_ONLY"))
        matching = view.copy()
        matching.q_extra_count = 0
        matching.nadd_max = 0
        matching.q_extra = np.zeros_like(view.q_extra)
        matching.q_extra_full = np.zeros(0, dtype=float)
        matching_full = destination / "matching_projection" / "internal"
        matching.write_internal_directory(matching_full)
        reloaded = EOSData.from_internal_directory(matching_full)
        for field in ("nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Iphase", "t_axis", "yq_axis", "q_full"):
            if not np.array_equal(getattr(view, field), getattr(reloaded, field)):
                raise ValueError(f"AUXILIARY_OMISSION_CHANGED_PRIMARY_FIELD: {field}")
        for field in ("mn", "mp", "inb_start", "temperature_mev", "charge_fraction_grid", "nT", "nYq"):
            if getattr(view, field) != getattr(reloaded, field):
                raise ValueError(f"AUXILIARY_OMISSION_CHANGED_SOURCE_METADATA: {field}")
        auxiliary_projection = dict(
            kind="EXPLICIT_THERMO_ONLY_AUXILIARY_OMISSION",
            source_nadd=int(view.q_extra_count), matching_nadd=0,
            full_source_internal_path=ref(full), auxiliary_path=ref(sidecar),
            auxiliary_sha256=sha(sidecar), source_columns_preserved_in_full_view=True,
            primary_columns_phase_axes_bitwise_equal=True,
            semantics="auxiliary_not_consumed_by_thermo_matcher")
        record["auxiliary_policy"] = auxiliary_projection
        record["matching_internal_path"] = ref(matching_full)
    record["matching_nadd"] = int(matching.q_extra_count)
    runs = phase_runs(view)
    whole = len(runs) == 1
    record["whole_selected_matching_allowed"] = whole
    record["whole_selected_matching_basis"] = ("one source-declared phase; native physical gates remain mandatory"
        if whole else "cross-phase interpolation is not authorized; source Iphase runs are the prescribed matching components")
    components = []
    for number, (label, start, stop) in enumerate(runs):
        cid = f"phase_{number:03d}"
        directory = matching_full if whole else destination / "components" / cid / "internal"
        if not whole:
            _segment_eosdata(matching, (label, start, stop)).write_internal_directory(directory)
        component = dict(id=cid, phase_id=label, start=start, stop=stop,
                         start_zero_based=start, stop_exclusive=stop,
                         rows=stop-start, internal_path=ref(directory),
                         source_original_start=record["selection"].get("source_start_zero_based", 0) + start,
                         source_original_stop=record["selection"].get("source_start_zero_based", 0) + stop,
                         density_fm3=[density[start], density[stop-1]],
                         source_support_scope="WHOLE_SELECTED_RETAINED_CHAIN" if whole else "FINITE_COMPONENT_ONLY",
                         component_basis="SOURCE_DECLARED_CONTIGUOUS_IPHASE_RUN",
                         physical_admission="NOT_TESTED_BY_PREPARATION", candidate_eligible=stop-start >= 3,
                         internal_sha256=hashes(directory, ("eos_data.csv", "eos_data.nml")))
        components.append(component)
    record["components"] = components
    record["matching_views"] = [dict(id="full_selected" if whole else c["id"],
        internal_path=c["internal_path"], start=c["start"], stop=c["stop"], rows=c["rows"],
        phase_id=c["phase_id"], methods=METHODS if c["candidate_eligible"] else [],
        source_support_scope=c["source_support_scope"], component_basis=c["component_basis"],
        auxiliary_projection=auxiliary_projection)
        for c in components]
    record["source_phase_count"] = len(runs)
    record["source_identity_sha256"] = hashlib.sha256(json.dumps({
        "source_id": record["source_id"], "source_sha256": record["source_sha256"],
        "selection": record["selection"]}, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return record


def prepare_snapshot(source: dict, destination: Path) -> tuple[dict, EOSData]:
    record = base_record(source)
    record["verified_historical_hashes"] = verify_old(source)
    record["source_sha256"] = source.get("source_table_sha256", {})
    record["source_path"] = source.get("source_path")
    path = absolute(source["compose_path"])
    record["input_compose_path"] = ref(path)
    record["input_compose_sha256"] = hashes(path)
    data = read_compose(path)
    selection = dict(source.get("selection") or {})
    if source.get("provenance_class") == "new_public":
        start, stop = (int(selection[k]) for k in ("source_start_zero_based", "source_stop_exclusive"))
        record["input_rows_before_selection"] = data.n
        if stop > data.n or start < 0 or stop <= start:
            raise ValueError("SOURCE_SELECTION_OUT_OF_RANGE")
        data = _segment_eosdata(data, (0, start, stop))
        selection["application"] = "APPLIED_ONCE_TO_HISTORICAL_FULL_PROJECTION"
        aux = path.parent / "SOURCE_AUXILIARY.json"
        check_sha(aux, selection["auxiliary_sha256"])
        destination.mkdir(parents=True, exist_ok=True)
        shutil.copy2(aux, destination / "SOURCE_AUXILIARY.json")
        record["auxiliary_policy"] = {
            "kind": "EXISTING_SOURCE_AUTHORIZED_Q1_Q7_PROJECTION",
            "source_nadd": selection["nadd_original"], "matching_nadd": int(data.q_extra_count),
            "matching_thermo_omits_original_auxiliary_columns": True,
            "exact_full_raw_and_auxiliary_sidecar_preserved": True,
            "auxiliary_path": ref(destination / "SOURCE_AUXILIARY.json"),
            "auxiliary_sha256": sha(destination / "SOURCE_AUXILIARY.json"),
            "projection_reapplied": False}
    elif selection:
        start, stop = (int(selection[k]) for k in ("source_start_zero_based", "source_stop_exclusive"))
        if data.n != stop - start:
            raise ValueError("SOURCE_PREPARED_SELECTION_ROWCOUNT_MISMATCH")
        selection["application"] = "ALREADY_SELECTED_SNAPSHOT_NOT_SLICED_AGAIN"
    else:
        selection = dict(source_start_zero_based=0, source_stop_exclusive=data.n,
                         excluded_lower_rows=0, excluded_upper_rows=0,
                         basis="entire finest-resolution source snapshot; no duplicate resolution variants",
                         application="ENTIRE_INPUT")
    record["selection"] = selection
    record["original_source_rows"] = source.get("original_source_rows") or record.get("input_rows_before_selection", data.n)
    return record, data


def raw_import(source: dict, spec: dict, raw: Path, destination: Path, importer: Path,
               density_bounds: tuple[str, str] | None = None) -> tuple[dict, EOSData]:
    record = base_record(source)
    raw_sha = check_sha(raw, spec["sha256"])
    projection = destination / "native_projection"
    command = [str(importer), spec["layout"], str(raw), str(projection),
               str(spec["mn"]), str(spec["mp"])]
    if density_bounds:
        command += list(density_bounds)
    destination.mkdir(parents=True, exist_ok=True)
    # Reuse only exact immutable bytes produced by the same importer invocation.
    run_record = destination / "IMPORT_INVOCATION.json"
    identity = dict(argv=command, executable_sha256=sha(importer), raw_sha256=raw_sha)
    if projection.exists():
        previous = json.loads(run_record.read_text())
        if previous["identity"] != identity or previous["output_sha256"] != hashes(projection, tuple(previous["output_sha256"])):
            raise ValueError("SOURCE_IMPORT_CACHE_IDENTITY_MISMATCH")
    else:
        result = subprocess.run(command, text=True, capture_output=True, check=False, timeout=120)
        (destination / "IMPORT.stdout").write_text(result.stdout)
        (destination / "IMPORT.stderr").write_text(result.stderr)
        if result.returncode:
            raise RuntimeError(f"NATIVE_IMPORT_FAILED({result.returncode}): {result.stderr.strip() or result.stdout.strip()}")
        save(run_record, dict(identity=identity, output_sha256=hashes(projection, tuple(p.name for p in projection.iterdir() if p.is_file()))))
    receipt = json.loads((projection / "IMPORT.json").read_text())
    if receipt["status"] != "NATIVE_RAW_PROJECTION_READBACK_VERIFIED":
        raise ValueError("SOURCE_IMPORT_NOT_READBACK_VERIFIED")
    data = EOSData.from_internal_directory(projection)
    record.update(source_path=ref(raw), source_sha256={raw.name: raw_sha},
        original_source_rows=int(receipt["raw_rows"]),
        selection=dict(source_start_zero_based=int(receipt["source_row_start_zero_based"]),
            source_stop_exclusive=int(receipt["source_row_stop_exclusive"]),
            excluded_lower_rows=int(receipt["excluded_lower_rows"]),
            excluded_upper_rows=int(receipt["excluded_upper_rows"]),
            basis="explicit native density selection" if density_bounds else "entire typed raw table",
            application="APPLIED_ONCE_INSIDE_NATIVE_IMPORTER",
            selected_density_fm3=receipt["selected_density_domain_fm3"]),
        compose_path=ref(projection), native_import_receipt=receipt,
        native_import_receipt_path=ref(projection / "IMPORT.json"),
        native_importer_sha256=sha(importer), source_layout=spec,
        phase_policy="RAW_TABLE_DECLARES_NO_PHASE_LABEL: single unsplit source branch; no inference from failed physical intervals")
    return record, data


def extra_sources() -> list[dict]:
    controls = json.loads(CONTROLS.read_text())
    extra = []
    generated = ROOT / "data/cowling_maxwell_native_20261119/same_rmf/generated_core"
    for model in ("DDME2", "FSU2", "TM1e"):
        full = "GPPVA_" + model
        receipt_path = generated / (full + ".json")
        receipt = json.loads(receipt_path.read_text())
        extra.append(dict(source_id=f"native_reconstructed:{full}:core", model=full, role="core",
            provenance_class="historical_homogeneous_generator_snapshot", status="SOURCE_SNAPSHOT_AVAILABLE",
            compose_path=ref(generated / full), source_path=ref(generated / full),
            source_table_sha256={"eos.thermo": receipt["output_eos_thermo_sha256"]},
            source_receipt_path=ref(receipt_path), source_receipt_sha256=sha(receipt_path),
            historical_generation_receipt=receipt, original_source_rows=receipt["n_points"],
            historical_evaluator=receipt["evaluator"], new_numerical_generation_this_turn=False))
        item = controls["models"][model]
        base = Path("data/cowling_maxwell_native_20261119/provider_original_unified") / model
        raw = PROVIDER_ROOT / base / "raw"
        projection = item.get("source_projection")
        used = PROVIDER_ROOT / projection["projection_path"] if projection else raw
        layout = item["original_layout"]
        start = layout["outer_rows"] + layout["inner_rows"]
        extra.append(dict(source_id=f"public_original:{full}:core", model=full, role="core",
            provenance_class="provider_original_homogeneous_core", status="SOURCE_SNAPSHOT_AVAILABLE",
            compose_path=str(used), source_path=str(raw),
            source_table_sha256=item["source_file_hashes"],
            source_receipt_path=ref(CONTROLS), source_receipt_sha256=sha(CONTROLS),
            original_source_rows=projection["original_rows"] if projection else layout["total_rows"],
            selection=dict(source_start_zero_based=start, source_stop_exclusive=layout["total_rows"],
                excluded_lower_rows=start, excluded_upper_rows=projection["original_rows"]-layout["total_rows"] if projection else 0,
                basis="provider original homogeneous phase; recorded unchanged monotone prefix where applicable"),
            provider_controls={k:item.get(k) for k in ("official_source_url", "provider_eos_zip_sha256", "provider_checksum_recorded",
                "provider_zip_checksum_matches", "independent_thermo_download_sha256", "independent_thermo_matches_zip_member", "source_projection")},
            apply_selection_once=True, projection_sha256=projection.get("projected_eos_thermo_sha256") if projection else None))
    return extra


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-exe", type=Path, required=True)
    parser.add_argument("--importer-exe", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "data/eos_matching_reaudit_20261123")
    parser.add_argument("--nuda-manifest", type=Path, default=DEFAULT_NUDA)
    parser.add_argument("--private-zip", type=Path, default=Path("/tmp/ns_eos_reaudit_20261123/private_intake/source.zip"))
    parser.add_argument("--private-output", type=Path, default=Path("/tmp/ns_eos_reaudit_20261123/private_intake/prepared"))
    args = parser.parse_args()
    native = args.native_exe.resolve()
    importer = (args.importer_exe or native.with_name("compose_eos_import_cold_table")).resolve()
    native_sha, importer_sha = sha(native), sha(importer)
    os.environ["COMPOSE_EOS_FORTRAN_EXE"] = str(native)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    check_sha(OLD_CATALOG, OLD_CATALOG_SHA)
    check_sha(args.nuda_manifest, NUDA_MANIFEST_SHA)
    old = json.loads(OLD_CATALOG.read_text())
    nuda = json.loads(args.nuda_manifest.read_text())
    wheel = Path(nuda["wheel"])
    check_sha(wheel, nuda["wheel_sha256"])
    typed = {s["source_id"]: s for s in nuda["records"]}
    source_list = old["sources"] + extra_sources()
    if len(source_list) != 101 or len({s["source_id"] for s in source_list}) != 101:
        raise ValueError("SOURCE_CENSUS_IDENTITY_OR_COUNT_MISMATCH")
    output = args.output.resolve()
    if args.private_output.resolve().is_relative_to(ROOT):
        raise ValueError("PRIVATE_OUTPUT_MUST_BE_OUTSIDE_GIT_WORKTREE")
    private_sha = check_sha(args.private_zip, PRIVATE_ZIP_SHA)
    private_raw = args.private_zip.parent
    records = []
    with zipfile.ZipFile(wheel) as archive:
        for source in source_list:
            sid = source["source_id"]
            folder = re.sub(r"[^A-Za-z0-9_.-]+", "__", sid)
            private = sid.startswith("user_zip:")
            destination = (args.private_output if private else output / "sources") / folder
            record = base_record(source)
            try:
                if sid in typed:
                    spec = typed[sid]
                    raw = destination / "raw" / Path(spec["path_in_wheel"]).name
                    raw.parent.mkdir(parents=True, exist_ok=True)
                    blob = archive.read(spec["path_in_wheel"])
                    if hashlib.sha256(blob).hexdigest() != spec["sha256"]:
                        raise ValueError("NUDA_WHEEL_MEMBER_HASH_MISMATCH")
                    raw.write_bytes(blob)
                    importer_spec = dict(spec, mn=spec["rest_mass_offset_mev"], mp=spec["proton_header_mass_mev"])
                    record, data = raw_import(source, importer_spec, raw, destination, importer)
                    record["source_wheel_sha256"] = nuda["wheel_sha256"]
                elif private:
                    core = source["role"] == "core"
                    name = "EoS.d" if core else "RHF-TF.data"
                    spec = dict(layout="user_core_energy_density" if core else "user_crust_energy_density",
                                sha256=source["raw_sha256"], mn=939.5654133, mp=938.2720813,
                                energy_gauge="source total energy density; zero added rest-mass offset")
                    record, data = raw_import(source, spec, private_raw / name, destination, importer,
                                             (".079", "1.7976931348623157e308") if core else ("0", ".06757"))
                    record["original_zip_sha256"] = private_sha
                    record["publication_policy"] = "HASHES_AND_METADATA_ONLY; source and derived tables private backup only"
                elif not source.get("compose_path"):
                    record["status"] = source["status"]
                    if source.get("nuda_path"):
                        member = source["nuda_path"]
                        blob = archive.read(member)
                        actual = hashlib.sha256(blob).hexdigest()
                        if actual != source["raw_sha256"]:
                            raise ValueError("BLOCKED_SOURCE_HASH_MISMATCH")
                        record["source_sha256"] = {Path(member).name: actual}
                    records.append(record)
                    print(json.dumps({"source_id":sid, "available":False, "status":record["status"]}), flush=True)
                    continue
                elif source.get("apply_selection_once"):
                    record = base_record(source)
                    record["verified_historical_hashes"] = verify_old(source)
                    path = absolute(source["compose_path"])
                    if source.get("projection_sha256"):
                        check_sha(path / "eos.thermo", source["projection_sha256"])
                    data = read_compose(path)
                    selection = dict(source["selection"])
                    start, stop = selection["source_start_zero_based"], selection["source_stop_exclusive"]
                    if data.n != stop:
                        raise ValueError("PROVIDER_AUTHORIZED_FULL_INPUT_ROWCOUNT_MISMATCH")
                    data = _segment_eosdata(data, (0, start, stop))
                    selection["application"] = "APPLIED_ONCE_TO_PROVIDER_OR_AUTHORIZED_PREFIX"
                    record.update(selection=selection, source_sha256=source["source_table_sha256"],
                        source_path=source["source_path"], input_compose_path=ref(path),
                        input_compose_sha256=hashes(path), original_source_rows=source["original_source_rows"],
                        provider_controls=source["provider_controls"])
                else:
                    record, data = prepare_snapshot(source, destination)
                record = freeze_view(record, data, destination)
            except Exception as exc:
                record.update(available=False, status="SOURCE_PREPARATION_FAILED",
                              failure={"type":type(exc).__name__, "message":str(exc)})
            records.append(record)
            print(json.dumps({"source_id":sid, "available":record["available"], "status":record["status"],
                "rows":record.get("full_selected_rows"), "components":len(record["components"]), "failure":record.get("failure")}), flush=True)
    check_sha(native, native_sha)
    check_sha(importer, importer_sha)
    counts = dict(total=len(records), roles=dict(Counter(s["role"] for s in records)),
                  available=sum(s["available"] for s in records),
                  available_roles=dict(Counter(s["role"] for s in records if s["available"])),
                  statuses=dict(Counter(s["status"] for s in records)),
                  phase_components=sum(len(s["components"]) for s in records),
                  cartesian_cells=31*70*len(METHODS))
    catalog = dict(schema="ns_cowling.eos_reaudit_source_catalog.v1", source_catalog=ref(OLD_CATALOG),
                   source_catalog_sha256=OLD_CATALOG_SHA, methods=METHODS,
                   source_head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                   native_executable=str(native), native_executable_sha256=native_sha,
                   native_importer=str(importer), native_importer_sha256=importer_sha,
                   nuda_manifest_sha256=NUDA_MANIFEST_SHA, nuda_wheel_sha256=nuda["wheel_sha256"],
                   preparation_policy="No numerical screening or physical-interval trimming; all 101 source records retained",
                   counts=counts, sources=records)
    save(output / "SOURCE_CATALOG.json", catalog)
    print(json.dumps({"catalog":ref(output / "SOURCE_CATALOG.json"), "counts":counts}), flush=True)
    if counts["statuses"].get("SOURCE_PREPARATION_FAILED"):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
