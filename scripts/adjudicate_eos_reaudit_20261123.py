#!/usr/bin/env python3
"""Attach a fixed-prefix native certificate without rewriting campaign receipts."""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def write_fresh(path, value):
    with path.open("x") as stream:
        stream.write(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--campaign", type=Path, required=True)
    p.add_argument("--certificate", type=Path, required=True)
    args = p.parse_args()
    root = args.campaign.resolve(strict=True)
    journal_raw = (root / "RESULTS.jsonl").read_bytes()
    summary_raw = (root / "SUMMARY.json").read_bytes()
    cert_raw = args.certificate.read_bytes()
    summary, cert = json.loads(summary_raw), json.loads(cert_raw)
    records = [json.loads(line) for line in journal_raw.splitlines()]
    if summary["status"] != "COMPLETED" or summary["input_and_binary_hashes_unchanged"] is not True:
        raise ValueError("completed immutable campaign required")
    if len(records) != summary["requested_cells"] or len({r["cell_id"] for r in records}) != len(records):
        raise ValueError("missing or duplicate campaign cells")
    if dict(Counter(r["status"] for r in records)) != summary["statuses"]:
        raise ValueError("completion census differs from journal")
    if cert["schema"] != "ns_cowling.fixed_retained_prefix_certificate.v1" or cert["status"] != "CERTIFIED":
        raise ValueError("unsupported or incomplete native certificate")
    if cert["classification"] != "FIXED_RETAINED_PREFIX_C11_REJECTION" or any(cert["mutations"].values()):
        raise ValueError("certificate is not an immutable fixed-prefix resolution")
    if cert["files"]["campaign_results"]["sha256"] != sha(journal_raw):
        raise ValueError("certificate targets a different journal")
    if cert["files"]["campaign_summary"]["sha256"] != sha(summary_raw):
        raise ValueError("certificate targets a different completion receipt")
    if cert["files"]["native_executable"]["sha256"] != summary["settings"]["binary_sha256"]:
        raise ValueError("certificate targets a different numerical binary")
    if cert["files"]["source_catalog"]["sha256"] != summary["settings"]["catalog_sha256"]:
        raise ValueError("certificate targets a different source census")
    pending_key = "post_bridge_gate_failure_search_incomplete"
    pending = {r["cell_id"]: r for r in records if r.get("prior_or_total_failure_categories", {}).get(pending_key, 0)}
    cells = cert["cells"]
    resolved = {r["cell_id"]: r for r in cells}
    if len(resolved) != len(cells) or set(resolved) != set(pending):
        raise ValueError("certificate does not cover exactly the pending cells")
    trials = 0
    for cell, r in pending.items():
        c = resolved[cell]
        if r["status"] != "NO_ADMITTED_MATCH" or c["original_status"] != r["status"]:
            raise ValueError("resolution cannot change the admitted survivor set")
        if r["crust_source_id"] != cert["source_id"] or r["core_source_id"] != c["core_source_id"] or r["method"] != c["method"]:
            raise ValueError("certificate source/method identity differs")
        if r["failure_trace"]["sha256"] != c["failure_trace_sha256"]:
            raise ValueError("failure trace identity differs")
        count = r["prior_or_total_failure_categories"][pending_key]
        if count != len(c["affected_trials"]):
            raise ValueError("certificate trial accounting differs")
        if any(t["native_failed_checks"] != "C11" or t["initial_i1_zero_based"] < 2 for t in c["affected_trials"]):
            raise ValueError("certificate contains a different failed gate or source support")
        trials += count
    if len(pending) != cert["affected_cells"] or trials != cert["affected_trials"]:
        raise ValueError("certificate accounting differs")
    classifications = []
    for r in sorted(records, key=lambda item: item["cell_id"]):
        classification = r["status"]
        reason = ""
        if r["status"] == "SURVIVED":
            classification = "SURVIVED_" + r["source_coverage"]
        elif r["cell_id"] in resolved:
            classification = "FIXED_RETAINED_PREFIX_C11_REJECTION"
            reason = sha(cert_raw)
        elif r.get("prior_or_total_failure_categories", {}).get("implementation_scope"):
            classification = "UNSUPPORTED_MULTIPLE_COEXISTENCE_ROOTS"
        classifications.append({k:r[k] for k in ("cell_id", "crust_source_id", "core_source_id", "method", "private")} |
                               dict(original_status=r["status"], final_classification=classification, certificate_sha256=reason))
    if (root / "RESULTS.jsonl").read_bytes() != journal_raw or (root / "SUMMARY.json").read_bytes() != summary_raw:
        raise ValueError("campaign changed during adjudication")
    decision = dict(
        schema="ns_cowling.campaign_adjudication.v1", status="RESOLVED_BY_FIXED_NATIVE_CERTIFICATE",
        journal_sha256=sha(journal_raw), summary_sha256=sha(summary_raw),
        certificate_sha256=sha(cert_raw), certificate_path=str(args.certificate),
        resolved_trial_count=trials, resolved_cell_count=len(resolved),
        remaining_post_bridge_search_incomplete_trials=0,
        original_journal_modified=False, original_counters_rewritten=False,
        survivor_set_changed=False, numerical_thresholds_changed=False,
        proof_scope=cert["scope"], cells=cells,
    )
    census = dict(
        schema="ns_cowling.final_adjudicated_census.v1", status="COMPLETED_WITH_EXPLICIT_SOURCE_AND_IMPLEMENTATION_LIMITS",
        requested_cells=len(records), original_status_counts=summary["statuses"],
        final_classification_counts=dict(Counter(r["final_classification"] for r in classifications)),
        source_catalog_sha256=summary["settings"]["catalog_sha256"],
        native_binary_sha256=summary["settings"]["binary_sha256"],
        immutable_journal_sha256=sha(journal_raw), certificate_sha256=sha(cert_raw),
        original_post_bridge_incomplete_trials=trials, resolved_post_bridge_trials=trials,
        resolved_post_bridge_cells=len(resolved), remaining_post_bridge_search_incomplete_trials=0,
        unresolved_implementation_scope_cells=sum(r["final_classification"] == "UNSUPPORTED_MULTIPLE_COEXISTENCE_ROOTS" for r in classifications),
        runtime_or_protocol_failures=json.loads((root / "ANALYSIS.json").read_bytes())["runtime_or_protocol_failures"],
        counting_scope="source identity x method cells, not semantically distinct EOS functions",
        negative_scope="current supplied views, finite endpoint grid and frozen admission; no universal physical impossibility claim",
    )
    write_fresh(root / "FINAL_ADJUDICATION.json", decision)
    write_fresh(root / "FINAL_CENSUS.json", census)
    with (root / "FINAL_CELL_CLASSIFICATION.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, list(classifications[0]))
        writer.writeheader()
        writer.writerows(classifications)
    print(json.dumps(census))


if __name__ == "__main__":
    main()
