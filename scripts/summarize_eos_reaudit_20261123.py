#!/usr/bin/env python3
"""Summarize immutable EOS campaign receipts; no numerical EOS construction."""
from __future__ import annotations
import argparse, csv, json
from collections import Counter, defaultdict
from pathlib import Path

METHODS = ("first_law_constrained_hermite_pressure_bridge", "cubic_F", "quintic_F",
           "cubic_P", "quintic_P", "css", "maxwell")

def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n")

def write_csv(path, fields, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    args = parser.parse_args()
    root = args.campaign
    summary = json.loads((root/"SUMMARY.json").read_text())
    if summary["status"] != "COMPLETED" or not summary["input_and_binary_hashes_unchanged"]:
        raise ValueError("Only a completed immutable campaign may be summarized")
    results = [json.loads(line) for line in (root/"RESULTS.jsonl").read_text().splitlines()]
    if len(results) != summary["requested_cells"] or len({r["cell_id"] for r in results}) != len(results):
        raise ValueError("Incomplete or duplicate Cartesian receipts")
    survivors = [r for r in results if r["status"] == "SURVIVED"]
    by_pair = defaultdict(dict)
    failed_trials = Counter()
    per_source = defaultdict(list)
    for result in results:
        per_source[("crust", result["crust_source_id"])].append(result)
        per_source[("core", result["core_source_id"])].append(result)
        failed_trials.update(result.get("prior_or_total_failure_categories", {}))
        if result["status"] == "SURVIVED":
            by_pair[(result["crust_source_id"], result["core_source_id"])][result["method"]] = result
    all_seven = []
    for (crust, core), methods in sorted(by_pair.items()):
        if set(methods) != set(METHODS):
            continue
        maxwell = methods["maxwell"]
        native = maxwell["native_report"]
        all_seven.append({
            "crust_source_id": crust, "core_source_id": core,
            "crust_model": maxwell["crust_model"], "core_model": maxwell["core_model"],
            "maxwell_cell_id": maxwell["cell_id"], "maxwell_source_coverage": maxwell["source_coverage"],
            "P_t_MeV_fm3": native["pressure_transition"],
            "mu_t_MeV": native["mu_transition"], "n_low_fm3": native["n_low"],
            "n_high_fm3": native["n_high"], "delta_epsilon_MeV_fm3": native["latent_energy_density"],
            "private": maxwell["private"],
            **{m+"_cell_id": methods[m]["cell_id"] for m in METHODS if m != "maxwell"},
        })
    method_rows = []
    for method in METHODS:
        group = [r for r in results if r["method"] == method]
        admitted = [r for r in group if r["status"] == "SURVIVED"]
        method_rows.append({
            "method": method, "requested": len(group),
            "available_source_cells": sum(r["status"] != "BLOCKED_SOURCE" for r in group),
            "survived": len(admitted),
            "selected_source_endpoint_chain": sum(r["source_coverage"] == "SELECTED_SOURCE_ENDPOINT_CHAIN" for r in admitted),
            "finite_phase_component_only": sum(r["source_coverage"] == "FINITE_PHASE_COMPONENT_ONLY" for r in admitted),
            "not_admitted": sum(r["status"] == "NO_ADMITTED_MATCH" for r in group),
            "blocked_source": sum(r["status"] == "BLOCKED_SOURCE" for r in group),
            "native_phase_pair_attempts": sum(r["native_candidate_attempts"] for r in group),
        })
    source_rows = []
    for (role, identifier), group in sorted(per_source.items()):
        admitted = [r for r in group if r["status"] == "SURVIVED"]
        source_rows.append({
            "role": role, "source_id": identifier, "requested_cells": len(group),
            "survived_cells": len(admitted),
            "maxwell_survived_cells": sum(r["method"] == "maxwell" for r in admitted),
            "contains_private_pair": any(r["private"] for r in group),
            **{m: sum(r["method"] == m for r in admitted) for m in METHODS},
        })
    analysis = {
        "schema": "ns_cowling.reaudit_census_summary.v1",
        "requested_cells": len(results), "available_source_cells": len(results)-summary["statuses"]["BLOCKED_SOURCE"],
        "surviving_cells": len(survivors), "surviving_source_pairs": len(by_pair),
        "all_seven_method_pairs": len(all_seven),
        "all_six_bridge_method_pairs": sum(set(METHODS)-{"maxwell"} <= set(m) for m in by_pair.values()),
        "survivor_coverage": dict(Counter(r["source_coverage"] for r in survivors)),
        "public_survivors": sum(not r["private"] for r in survivors),
        "private_survivors": sum(r["private"] for r in survivors),
        "native_phase_pair_attempts": sum(r["native_candidate_attempts"] for r in results),
        "failed_native_or_admission_trials": dict(failed_trials),
        "post_bridge_incomplete_trials": failed_trials["post_bridge_gate_failure_search_incomplete"],
        "runtime_or_protocol_failures": sum(
            v for k,v in failed_trials.items()
            if k in ("runtime_timeout", "runtime_signal", "protocol_failure", "python_orchestration_failure")),
        "important_counting_rule": "trial failure counts include multiple prescribed phase-pair attempts and cannot be added to cell counts",
        "source_identity_rule": "original/corrected/reconstructed variants remain distinct records; counts do not mean independent nuclear functionals",
        "by_method": method_rows,
    }
    with (root/"SURVIVORS.csv").open(newline="") as source_csv:
        reader = csv.DictReader(source_csv)
        fields = list(reader.fieldnames or [])
        flat_survivors = list(reader)
    if len(flat_survivors) != len(survivors) or {r["cell_id"] for r in flat_survivors} != {r["cell_id"] for r in survivors}:
        raise ValueError("CSV survivor inventory differs from journal")
    write_csv(root/"ENDPOINT_CHAIN_SURVIVORS.csv", fields,
              (r for r in flat_survivors if r["source_coverage"] == "SELECTED_SOURCE_ENDPOINT_CHAIN"))
    write_csv(root/"FINITE_COMPONENT_SURVIVORS.csv", fields,
              (r for r in flat_survivors if r["source_coverage"] == "FINITE_PHASE_COMPONENT_ONLY"))
    write_json(root/"ANALYSIS.json", analysis)
    write_csv(root/"METHOD_COUNTS.csv", list(method_rows[0]), method_rows)
    write_csv(root/"SOURCE_OUTCOME_COUNTS.csv", list(source_rows[0]), source_rows)
    if all_seven:
        write_csv(root/"ALL_SEVEN_SURVIVOR_PAIRS.csv", list(all_seven[0]), all_seven)
    print(json.dumps(analysis, ensure_ascii=False))
if __name__ == "__main__":
    main()
