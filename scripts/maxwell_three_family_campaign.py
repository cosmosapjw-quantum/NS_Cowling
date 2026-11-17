#!/usr/bin/env python3
"""Execute all same-family registered core/crust Maxwell combinations.

The source diagnostic is native, independent of a matching partner, and never
repairs a source. Every connected admissible component (at least three knots)
is retained. Production roots, interpolation and physical gates are Fortran.
This runner records limited support; an interval match is not a complete-star
EOS or a proof that excluded source intervals cannot contain a transition.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "python_frontend"), str(ROOT / "scripts")]
from compose_eos.api import read_compose, write_compose, diagnose_maxwell_phase, _thermo_only_view
import compose_eos.api as eos_api
from compose_eos.eos_data import EOSData
from compose_eos.maxwell_v2 import construct_maxwell_transition_v2, _segment_eosdata
from compose_eos.native import executable, run_native
from match_report_admission import admit_match_report
from match_core_crust import verify_internal_readback, verify_internal_metadata_readback
from reconstruct_public_thermo_20261117 import reconstruct_data

FAMILIES = ("native_generated", "public_original", "public_reconstructed")
BASE = ROOT / "data/maxwell_three_families_20261117"
CORE_EXTENDED = ("GM_GM1", "RG_KDE0v1", "RG_SLY4")
FILES = ("eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path):
    return str(Path(path).resolve().relative_to(ROOT))


def plain(value):
    if is_dataclass(value):
        return plain(asdict(value))
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, slice):
        return {"start": value.start, "stop": value.stop, "step": value.step}
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(plain(value), indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def load(path):
    return json.loads(Path(path).read_text())


def safe_path(relative):
    p = (ROOT / relative).resolve()
    if not p.is_relative_to(ROOT) or Path(relative).is_absolute():
        raise ValueError("repository-relative source path required")
    return p


def verify_snapshot(entry):
    directory = safe_path(entry["path"])
    for name, expected in entry["table_sha256"].items():
        if sha(directory / name) != expected:
            raise ValueError(f"source changed after inventory: {entry['model']}/{name}")
    if entry.get("receipt_path") and sha(safe_path(entry["receipt_path"])) != entry["receipt_sha256"]:
        raise ValueError("source receipt changed after inventory")
    return directory


def extended_core_entry(old):
    """Prefer a newly evaluated larger native domain, only after sealed success."""
    model = old["model"]
    parent = BASE / "native_core_extended_v2" / model
    path = parent / "generated"
    if model not in CORE_EXTENDED or not path.is_dir():
        return dict(old)
    receipt = load(parent / "DISPATCH.json")
    if receipt.get("status") != "REQUESTED_ROWS_NATIVE_VERIFIED" or receipt.get("returncode") != 0:
        raise ValueError("extended native core has no successful dispatch receipt")
    for name, expected in receipt["evidence"]["output_sha256"].items():
        if sha(path / name) != expected:
            raise ValueError("extended native core receipt hash mismatch")
    return {**old, "path": rel(path), "receipt_path": rel(parent / "DISPATCH.json"),
            "receipt_sha256": sha(parent / "DISPATCH.json"),
            "table_sha256": {name: sha(path/name) for name in FILES},
            "request": receipt["request"], "replaces_finite_source_path": old["path"],
            "selection_reason": "partner-independent homogeneous functional evaluation on 0.02-1.0 fm^-3; physical gates applied afterwards"}


def role_slice(data, entry, role, family):
    """Exclude attached provider crust from core; select inner crust explicitly."""
    start, stop = 0, data.n
    basis = "entire native generated finite domain"
    if family != "native_generated":
        if role == "core":
            selection = entry["homogeneous_core_selection"]
            start, stop = selection["start_zero_based"], selection["stop_exclusive"]
            basis = selection["evidence"]
        else:
            # These thirteen provider archives explicitly label their computed
            # inner crust 7 and their attached BSk22 lower crust 0.
            selected = np.flatnonzero(data.Iphase == 7)
            if not len(selected) or np.any(np.diff(selected) != 1):
                raise ValueError("registered public inner-crust phase7 is not one contiguous supplied range")
            start, stop = int(selected[0]), int(selected[-1])+1
            basis = "provider phase7 inner crust; attached phase0 lower crust remains preserved in immutable original"
    result = _segment_eosdata(data, (0, start, stop))
    # A slice's beta-summary is recomputed from its actual Q5 values, exactly
    # as the native selected-slice reader does; this does not change any Q5.
    result.beta_equilibrated = bool(np.max(np.abs(result.Q5)) <= 1e-8)
    return result, {"source_start_zero_based": start, "source_stop_exclusive": stop,
                    "excluded_lower_rows": start, "excluded_upper_rows": data.n-stop,
                    "basis": basis, "selected_density_fm3": [result.nb[0], result.nb[-1]],
                    "selected_beta_flag_basis": "native Q5 threshold on exact selected rows; no physical Q5 changes"}


def valid_components(intervals, n):
    """Maximal connected original intervals, never bridge a rejected interval."""
    if len(intervals) != n-1:
        raise ValueError("native diagnostic did not cover every source interval")
    components, start = [], None
    for i, item in enumerate(intervals + [{"valid": False}]):
        if i < n-1 and (item["left_row_zero_based"] != i or item["right_row_zero_based"] != i+1):
            raise ValueError("native diagnostic interval indexing mismatch")
        if item["valid"] and start is None:
            start = i
        if not item["valid"] and start is not None:
            if i-start+1 >= 3:
                components.append((start, i+1))
            start = None
    return components


def auxiliary_projection(data, entry, selection, target):
    """Preserve provider payload; carry documented geometry into phase support.

    No thermodynamic Q1-Q7 column is modified. Unknown quantities are retained
    without assigning guessed physical meanings. Nadd omission is explicit.
    """
    result = data.copy()
    count = result.q_extra_count
    receipt = {"input_nadd": count, "thermodynamic_columns_modified": False}
    if not count:
        return result, receipt
    if count != 1 or result.q_extra.shape != result.nb.shape:
        raise ValueError("unhandled auxiliary layout; require lossless indexed payload projection")
    semantics = entry.get("extra_column_semantics", "UNKNOWN_PROVIDER_AUXILIARY")
    payload = {"source_path": entry["path"], "source_table_sha256": entry["table_sha256"],
               "semantics": semantics, "primary_evidence": entry.get("extra_column_primary_evidence"),
               "rows": [{"source_row_zero_based": selection["source_start_zero_based"]+i,
                         "n_b_fm3": float(n), "source_iphase": int(phase), "nadd": count,
                         "payload": [float(extra)]} for i, (n, phase, extra) in
                        enumerate(zip(result.nb, result.Iphase, result.q_extra))]}
    if entry["model"].startswith("GPPVA_") and entry.get("family") == "public_original":
        if not entry.get("extra_column_primary_evidence"):
            raise ValueError("GPPVA geometry projection requires its actual provider PDF evidence")
        geometry = np.rint(result.q_extra).astype(int)
        if np.any(np.abs(result.q_extra-geometry) > 1e-12) or np.any((geometry < 0) | (geometry > 6)):
            raise ValueError("provider geometry payload is not its documented integer code")
        pairs = sorted(set(zip(result.Iphase.tolist(), geometry.tolist())))
        labels = {pair: i+1 for i, pair in enumerate(pairs)}
        result.Iphase = np.array([labels[pair] for pair in zip(result.Iphase.tolist(), geometry.tolist())], dtype=np.int64)
        receipt["working_phase_dictionary"] = [{"working_iphase": value, "source_iphase": key[0],
                                                "provider_geometry": key[1]} for key, value in labels.items()]
        receipt["geometry_boundary_policy"] = "no interpolation across a change in documented geometry; no inferred first-order jump"
    save(target / "SOURCE_AUXILIARY.json", payload)
    receipt.update(payload_path=rel(target / "SOURCE_AUXILIARY.json"),
                   payload_sha256=sha(target / "SOURCE_AUXILIARY.json"),
                   matching_transport="explicit Nadd=0 thermodynamic-only projection; source payload retained separately")
    result.q_extra[:] = 0.0
    result.q_extra_count = result.nadd_max = 0
    result.q_extra_full = np.zeros(0)
    return result, receipt


def prepare_one(entry, family, role, output):
    directory = verify_snapshot(entry)
    original = read_compose(directory)
    selected, selection = role_slice(original, entry, role, family)
    target = output / family / entry["model"]
    target.mkdir(parents=True, exist_ok=False)
    record = {"id": f"{family}:{entry['model']}", "model": entry["model"],
              "family": family, "role": role, "source_path": entry["path"],
              "source_table_sha256": entry["table_sha256"], "selection": selection,
              "original_source_rows": original.n, "source_receipt_path": entry.get("receipt_path"),
              "source_receipt_sha256": entry.get("receipt_sha256"),
              "model_variant": entry.get("request", {}).get("variant"),
              "source_contract": entry.get("request", {}).get("contract", {}),
              "microscopic_source_limitations": entry.get("scientific_scope", entry.get("support", {})),
              "state": "PREPARING", "components": []}
    selected, projection = auxiliary_projection(selected, entry, selection, target)
    record["auxiliary_projection"] = projection
    if family == "public_reconstructed":
        selected, reconstruction = reconstruct_data(selected, raw_source_dir=directory, model=entry["model"])
        save(target / "RECONSTRUCTION.json", reconstruction)
        record["reconstruction_receipt"] = rel(target / "RECONSTRUCTION.json")
        record["reconstruction_receipt_sha256"] = sha(target / "RECONSTRUCTION.json")
    else:
        record["thermodynamic_reconstruction"] = "NONE"
    # Abundances are not transported by the sharp one-charge kernel. All Q1-Q7,
    # phase labels, source support and provenance remain unchanged here.
    prepared = _thermo_only_view(selected)
    record["composition_policy"] = "explicitly omitted from thermo-only matching; immutable source composition retained"
    record["nadd"] = prepared.q_extra_count
    prepared.write_internal_directory(target / "internal")
    record["internal_path"] = rel(target / "internal")
    record["internal_sha256"] = {name: sha(target / "internal" / name) for name in ("eos_data.csv", "eos_data.nml")}
    export = prepared.copy()
    if np.any(export.Iphase != 0):
        export.composition_file_present = True
        export.row_status[:] = 1
    write_compose(target / "compose", export)
    record["compose_path"] = rel(target / "compose")
    record["compose_sha256"] = {p.name: sha(p) for p in (target / "compose").glob("eos.*")}
    diagnostic = diagnose_maxwell_phase(prepared)
    save(target / "NATIVE_INTERVAL_DIAGNOSTIC.json", diagnostic)
    record["diagnostic_path"] = rel(target / "NATIVE_INTERVAL_DIAGNOSTIC.json")
    record["interval_failures_by_category"] = dict(Counter(x["category"] for x in diagnostic["intervals"] if not x["valid"]))
    ranges = valid_components(diagnostic["intervals"], prepared.n)
    for i, (a, b) in enumerate(ranges):
        record["components"].append({"id": f"component_{i:03d}", "start_zero_based": a, "stop_exclusive": b,
                                     "original_source_rows_zero_based": [selection["source_start_zero_based"]+a,
                                                                          selection["source_start_zero_based"]+b-1],
                                     "rows": b-a, "phase_id": int(prepared.Iphase[a]),
                                     "density_fm3": [prepared.nb[a], prepared.nb[b-1]],
                                     "mu_E_MeV": [float((prepared.eps[a]+prepared.p[a])/prepared.nb[a]),
                                                  float((prepared.eps[b-1]+prepared.p[b-1])/prepared.nb[b-1])]})
    record["state"] = "FINITE_COMPONENTS_AVAILABLE" if ranges else "NO_ADMISSIBLE_THREE_KNOT_COMPONENT"
    record["whole_selected_domain_admissible"] = not record["interval_failures_by_category"]
    record["matching_binary_sha256"] = sha(executable())
    save(target / "SOURCE.json", record)
    print(json.dumps({"prepared": record["id"], "components": len(ranges), "failures": record["interval_failures_by_category"]}), flush=True)
    return plain(record)


def prepare(inventory_path, output):
    inventory = load(inventory_path)
    output.mkdir(parents=True, exist_ok=False)
    entries = {"native_generated": inventory["native"], "public_original": inventory["public_original"],
               "public_reconstructed": inventory["public_original"]}
    records = []
    for family in FAMILIES:
        for initial in entries[family]:
            entry = extended_core_entry(initial) if family == "native_generated" else initial
            role = "core" if entry["model"] in inventory["core_ids"] else "crust"
            try:
                record = prepare_one(entry, family, role, output)
            except Exception as exc:
                record = {"id": f"{family}:{entry['model']}", "family": family, "role": role, "model": entry["model"],
                          "source_path": entry["path"], "state": "PREPARATION_FAILED", "components": [],
                          "exception_type": type(exc).__name__, "error": str(exc),
                          "native_category": getattr(exc, "category", None)}
                save(output / family / entry["model"] / "SOURCE.json", record)
                print(json.dumps({"prepared": record["id"], "error": record["error"]}), flush=True)
            records.append(record)
    manifest = {"schema": "ns-cowling.maxwell-three-family-preparation.v1", "created_at_utc": datetime.now(timezone.utc).isoformat(),
                "inventory_path": rel(inventory_path), "inventory_sha256": sha(inventory_path),
                "core_ids": inventory["core_ids"], "crust_ids": inventory["crust_ids"],
                "families": list(FAMILIES), "source_records": records,
                "rule": "partner-independent maximal native-admissible contiguous components; no changed source values to manufacture roots",
                "matching_executable_sha256": sha(executable())}
    save(output / "MANIFEST.json", manifest)
    return manifest


def read_prepared(record):
    path = safe_path(record["internal_path"])
    for name, expected in record["internal_sha256"].items():
        if sha(path/name) != expected:
            raise ValueError("prepared source changed after snapshot")
    original = safe_path(record["source_path"])
    for name, expected in record["source_table_sha256"].items():
        if sha(original/name) != expected:
            raise ValueError("original source changed after preparation")
    return EOSData.from_internal_directory(path, source_dir=str(original))


def check_same_family(crust, core):
    if crust["family"] != core["family"] or crust["family"] not in FAMILIES:
        raise ValueError("cross-family matching is forbidden")
    if crust["role"] != "crust" or core["role"] != "core":
        raise ValueError("crust/core roles reversed or unspecified")


@contextmanager
def observe_kernel(state):
    original = eos_api.run_native
    def observed(action, *args):
        if action == "maxwell-transition":
            state["kernel_invoked"] = True
        return original(action, *args)
    eos_api.run_native = observed
    try:
        yield
    finally:
        eos_api.run_native = original


def run_candidate(crust, core, low, high, lc, hc, destination):
    check_same_family(crust, core)
    destination.mkdir(parents=True, exist_ok=False)
    result = {"crust_component": lc, "core_component": hc, "status": "STARTED", "kernel_invoked": False}
    try:
        ldata = _segment_eosdata(low, (lc["phase_id"], lc["start_zero_based"], lc["stop_exclusive"]))
        hdata = _segment_eosdata(high, (hc["phase_id"], hc["start_zero_based"], hc["stop_exclusive"]))
        # Local recomputation of metadata only; Q5 values remain untouched.
        for d in (ldata, hdata):
            d.beta_equilibrated = bool(np.max(np.abs(d.Q5)) <= 1e-8)
        with observe_kernel(result):
            matched, report = construct_maxwell_transition_v2(ldata, hdata,
                assume_local_charge_neutrality=True, return_report=True)
        native = report.native_report
        result["v2_report"] = plain(report)
        matched.write_internal_directory(destination / "expected_internal")
        write_compose(destination / "compose", matched)
        back = read_compose(destination / "compose")
        back.write_internal_directory(destination / "readback_internal")
        numerical = verify_internal_readback(destination / "expected_internal/eos_data.csv",
            destination / "readback_internal/eos_data.csv", expected_n=matched.n, actual_n=back.n)
        metadata = verify_internal_metadata_readback(destination / "expected_internal/eos_data.nml",
            destination / "readback_internal/eos_data.nml")
        result["readback"] = {"numerical": numerical, "metadata": metadata}
        if numerical["status"] != "PASS" or metadata["status"] != "PASS":
            raise ValueError("serialized matched EOS failed lossless readback")
        result["admission"] = admit_match_report("maxwell", native,
            destination / "readback_internal/eos_data.csv", destination / "readback_internal/eos_data.nml")
        if result["admission"]["status"] != "ADMITTED":
            raise ValueError("independent serialized-interface admission rejected")
        # The v1 API already checked its internal output; independently invoke
        # the verifier on the actual delivered CompOSE bytes as well.
        verification = run_native("verify", destination / "compose", destination / "CONSISTENCY.json")
        result["serialized_native_verification"] = verification
        if verification.get("overall") != "PASS":
            raise ValueError("native serialized consistency verification rejected")
        result["status"] = "GENERATED_VERIFIED_FINITE_DOMAIN"
        result["output_path"] = rel(destination / "compose")
        result["output_sha256"] = {p.name: sha(p) for p in (destination / "compose").glob("eos.*")}
        result["density_fm3"] = [matched.nb[0], matched.nb[-1]]
        result["transition"] = {k: native[k] for k in ("mu_transition", "pressure_transition", "n_low", "n_high", "epsilon_low", "epsilon_high", "density_jump", "latent_energy_density", "maxwell_identity_residual")}
        result["scope"] = "two admitted original connected components joined by one native sharp Maxwell transition; finite inner-crust/core support; no whole-star or source-microphysics equivalence claim"
        result["source_exclusions"] = {"crust_before_component": lc["start_zero_based"], "core_after_component": high.n-hc["stop_exclusive"],
                                       "excluded_intervals_not_interpolated": True}
    except Exception as exc:
        result["status"] = "REJECTED"
        result["exception_type"] = type(exc).__name__
        result["error"] = str(exc)
        result["category"] = getattr(exc, "category", None) or "EXECUTION_OR_VALIDATION_FAILURE"
        if hasattr(exc, "report"):
            result["v2_report"] = plain(exc.report)
    save(destination / "RESULT.json", result)
    return plain(result)


def execute(prepared_manifest, output):
    manifest = load(prepared_manifest)
    if manifest["matching_executable_sha256"] != sha(executable()):
        raise ValueError("native matching binary differs from interval preparation; explicitly regenerate source diagnostics")
    output.mkdir(parents=True, exist_ok=False)
    records = {(r["family"], r["model"]): r for r in manifest["source_records"]}
    if set(manifest["families"]) != set(FAMILIES):
        raise ValueError("family universe differs from fixed three-family contract")
    prepared = {r["id"]: read_prepared(r) for r in records.values() if r.get("internal_path")}
    cells = []
    start = time.monotonic()
    for family, crust_id, core_id in itertools.product(FAMILIES, manifest["crust_ids"], manifest["core_ids"]):
        crust, core = records[(family, crust_id)], records[(family, core_id)]
        check_same_family(crust, core)
        cell_id = f"{crust_id}__{core_id}"
        destination = output / family / cell_id
        destination.mkdir(parents=True, exist_ok=False)
        cell = {"family": family, "id": cell_id, "crust": crust_id, "core": core_id,
                "crust_source_id": crust["id"], "core_source_id": core["id"],
                "candidate_results": [], "source_states": {"crust": crust["state"], "core": core["state"]},
                "source_interval_failures": {"crust": crust.get("interval_failures_by_category"), "core": core.get("interval_failures_by_category")}}
        for i, (lc, hc) in enumerate(itertools.product(crust["components"], core["components"])):
            result = run_candidate(crust, core, prepared[crust["id"]], prepared[core["id"]], lc, hc,
                                   destination / f"candidate_{i:03d}")
            cell["candidate_results"].append(result)
        successful = [r for r in cell["candidate_results"] if r["status"] == "GENERATED_VERIFIED_FINITE_DOMAIN"]
        cell["generated_count"] = len(successful)
        cell["kernel_invocations"] = sum(r["kernel_invoked"] for r in cell["candidate_results"])
        cell["status"] = ("GENERATED_VERIFIED_FINITE_DOMAIN" if successful else "NO_MATCH_ON_ADMITTED_COMPONENTS"
                          if cell["candidate_results"] else "SOURCE_HAS_NO_ADMITTED_COMPONENT_PAIR")
        cell["root_search_domain_is_complete_source"] = bool(crust.get("whole_selected_domain_admissible") and core.get("whole_selected_domain_admissible"))
        save(destination / "CELL.json", cell)
        cells.append(cell)
        print(json.dumps({"family": family, "pair": cell_id, "status": cell["status"], "generated": len(successful),
                          "attempts": cell["kernel_invocations"]}), flush=True)
    summary = {family: {"model_pairs": sum(c["family"] == family for c in cells),
                         "pairs_with_output": sum(c["family"] == family and c["generated_count"] > 0 for c in cells),
                         "generated_tables": sum(c["generated_count"] for c in cells if c["family"] == family),
                         "native_invocations": sum(c["kernel_invocations"] for c in cells if c["family"] == family),
                         "pair_statuses": dict(Counter(c["status"] for c in cells if c["family"] == family)),
                         "candidate_rejections": dict(Counter(r["category"] for c in cells if c["family"] == family
                                                              for r in c["candidate_results"] if r["status"] == "REJECTED"))}
               for family in FAMILIES}
    receipt = {"schema": "ns-cowling.maxwell-three-family-campaign.v1", "created_at_utc": datetime.now(timezone.utc).isoformat(),
               "prepared_manifest": rel(prepared_manifest), "prepared_manifest_sha256": sha(prepared_manifest),
               "matching_executable_sha256": sha(executable()), "campaign_script_sha256": sha(__file__),
               "source_preparation_snapshot_unchanged": True, "same_family_only": True,
               "git_base_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
               "wall_seconds": time.monotonic()-start, "summary": summary, "cells": cells}
    # Check all input bytes again after real production. Successful outputs do
    # not excuse a source mutation during execution.
    for r in records.values():
        if r.get("internal_path"):
            read_prepared(r)
    save(output / "CAMPAIGN.json", receipt)
    return receipt


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--inventory", type=Path, default=BASE / "intake/SOURCE_INVENTORY.json")
    p.add_argument("--output", type=Path, default=BASE / "prepared")
    e = sub.add_parser("execute")
    e.add_argument("--manifest", type=Path, default=BASE / "prepared/MANIFEST.json")
    e.add_argument("--output", type=Path, default=BASE / "campaign")
    args = parser.parse_args(argv)
    result = prepare(args.inventory, args.output) if args.command == "prepare" else execute(args.manifest, args.output)
    print(json.dumps(result.get("summary", {"prepared_sources": len(result.get("source_records", []))})), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
