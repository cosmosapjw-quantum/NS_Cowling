#!/usr/bin/env python3
"""Frontend comparison of C19/C20 serialized science; no stellar calculation."""
from __future__ import annotations
import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path

SOURCE_SHA = "299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
FREQUENCY_ATOL_HZ = .005  # inherited C13 PCHIP migration comparator
BACKGROUND_RTOL = 2e-7  # inherited C15 independent native TOV comparator
PHYSICAL_GATE = 1e-5
ADMITTED = {"LOCAL_NODE_CONTINUATION_ADMITTED", "OVERLAP_AND_NODE_CONTINUATION_ADMITTED"}
OVERLAP = {"absolute_radial_W_overlap", "high_phase_overlap", "low_phase_overlap"}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def audit(baseline, current):
    errors = []
    metrics = dict(max_root_frequency_difference_hz=0.,
                   max_scan_determinant_difference=0.,
                   max_overlap_difference=0., compared_scan_points=0,
                   compared_candidates=0, compared_branches=0)
    backgrounds = {k: 0. for k in
        ("mass_msun", "radius_km", "maxwell_interface_r_km", "compactness")}
    artifacts, counts = {}, {}

    def check(condition, label):
        if not condition:
            errors.append(label)

    def equal(a, b, label):
        check(a == b, label + ": exact contract mismatch")

    def finite(obj, path):
        if isinstance(obj, float):
            check(math.isfinite(obj), path + ": nonfinite JSON number")
        elif isinstance(obj, dict):
            for k, v in obj.items():
                finite(v, path + "." + k)
        elif isinstance(obj, list):
            for i, v in enumerate(obj):
                finite(v, path + "[" + str(i) + "]")

    def freq(a, b, path):
        equal(a is None, b is None, path + ".presence")
        if a is not None and b is not None:
            delta = abs(a-b)
            metrics["max_root_frequency_difference_hz"] = max(
                metrics["max_root_frequency_difference_hz"], delta)
            check(delta <= FREQUENCY_ATOL_HZ, path + ": frequency migration limit")

    def overlap_differences(a, b):
        if isinstance(a, dict) and isinstance(b, dict):
            for k in a.keys() & b.keys():
                if k in OVERLAP:
                    metrics["max_overlap_difference"] = max(
                        metrics["max_overlap_difference"], abs(a[k]-b[k]))
                else:
                    overlap_differences(a[k], b[k])
        elif isinstance(a, list) and isinstance(b, list):
            for x, y in zip(a, b):
                overlap_differences(x, y)

    def load(root, model):
        path = root / (model + ".json")
        raw = path.read_bytes()
        doc = json.loads(raw)
        finite(doc, str(path))
        compressed = (root / doc["original_scan_payload_path"]).read_bytes()
        packed = gzip.decompress(compressed)
        equal(digest(compressed), doc["original_scan_payload_gzip_sha256"], str(path)+".gzip")
        equal(digest(packed), doc["original_scan_payload_uncompressed_sha256"], str(path)+".scan")
        scans = json.loads(packed)
        finite(scans, str(path)+".scans")
        artifacts[str(path)] = dict(json_sha256=digest(raw), scan_gzip_sha256=digest(compressed))
        for row in doc["densities"]:
            key = format(row["density_factor"], ".5f")
            check(key in scans, str(path)+": missing scan density "+key)
            if key in scans:
                buf = json.dumps(scans[key], sort_keys=True, separators=(",", ":"),
                                 allow_nan=False).encode()
                equal(digest(buf), row["native_scan_points_sha256"], str(path)+"."+key+".scan")
        return doc, scans

    for model in ("M2", "M3"):
        old, old_scan = load(baseline, model)
        new, new_scan = load(current, model)
        equal(new["original_source_archive_sha256"], SOURCE_SHA, model+".source")
        for key in old:
            if key not in {"densities", "original_scan_payload_gzip_sha256",
                           "original_scan_payload_uncompressed_sha256"}:
                equal(old[key], new.get(key), model+"."+key)
        equal(len(old["densities"]), 20, model+".baseline_count")
        equal(len(new["densities"]), 20, model+".current_count")
        equal(set(old_scan), set(new_scan), model+".scan_density_keys")
        counts[model] = {"f": Counter(), "p1": Counter()}
        for a, b in zip(old["densities"], new["densities"]):
            label = model + "." + str(a["density_factor"])
            for key in ("density_factor", "n_c_fm3", "source_restart_count",
                        "scan_backend", "scan_grid_hz", "scan_point_failures",
                        "scan_sign_change_candidate_count", "background_backend",
                        "relative_to_transition_high_density", "local_branch_classification",
                        "mode_order_valid"):
                equal(a.get(key), b.get(key), label+"."+key)
            equal(b["scan_point_failures"], 0, label+".scan_failures")
            for key in backgrounds:
                delta = abs(a[key]-b[key])/abs(a[key])
                backgrounds[key] = max(backgrounds[key], delta)
                check(delta <= BACKGROUND_RTOL, label+"."+key+": migration limit")
            ca, cb = a["candidate_eigenvalues"], b["candidate_eigenvalues"]
            equal(len(ca), len(cb), label+".candidate_count")
            for i, (x, y) in enumerate(zip(ca, cb)):
                cpath = label+".candidate"+str(i)
                metrics["compared_candidates"] += 1
                for key in ("status", "bracket_hz", "W_nodes", "V_nodes",
                            "original_strong_knot_restarts", "phasewise_nodal_backend",
                            "source_eos_rhs_callback_calls"):
                    equal(x.get(key), y.get(key), cpath+"."+key)
                freq(x.get("frequency_hz"), y.get("frequency_hz"), cpath)
                if y["status"] == "PHYSICAL_ROOT_ADMITTED":
                    check(abs(y["normalized_signed_determinant"]) < PHYSICAL_GATE,
                          cpath+": determinant gate")
                    check(0 <= y["relative_join_mismatch"] < PHYSICAL_GATE,
                          cpath+": independent full-state join gate")
                    equal(y["source_eos_rhs_callback_calls"], 0, cpath+".callbacks")
            for mode in ("f", "p1"):
                x, y = a["local_branches"][mode], b["local_branches"][mode]
                bpath = label+"."+mode
                metrics["compared_branches"] += 1
                counts[model][mode][y["status"]] += 1
                for key in ("status", "selector_status", "candidate_index",
                            "W_nodes", "V_nodes", "high_low_shape_min"):
                    equal(x.get(key), y.get(key), bpath+"."+key)
                freq(x.get("frequency_hz"), y.get("frequency_hz"), bpath)
                if y["status"] in ADMITTED:
                    equal(y["W_nodes"], 0 if mode=="f" else 1, bpath+".W_identity")
                    equal(y["V_nodes"], 0 if mode=="f" else 1, bpath+".V_identity")
                    check(y["relative_join_mismatch"] < PHYSICAL_GATE, bpath+".join")
                    check(y["relative_frequency_step"] <= .12, bpath+".step")
                    if y["status"] == "OVERLAP_AND_NODE_CONTINUATION_ADMITTED":
                        d = y["cross_density_shape_overlap"]
                        check(d["high_phase_overlap"] >= .95, bpath+".high_phase")
                        check(d["low_phase_overlap"] >= .95, bpath+".low_phase")
                overlap_differences(x, y)
            key = format(a["density_factor"], ".5f")
            left, right = old_scan[key], new_scan[key]
            equal(len(left), 111, label+".baseline_scan_size")
            equal(len(right), 111, label+".current_scan_size")
            for i, (x, y) in enumerate(zip(left, right)):
                spath = label+".scan"+str(i)
                metrics["compared_scan_points"] += 1
                for k in ("frequency_hz", "source_restarts", "status"):
                    equal(x[k], y[k], spath+"."+k)
                equal(y["frequency_hz"], new["scan_grid_hz"][i], spath+".grid")
                check(x["determinant"]*y["determinant"] > 0 or
                      x["determinant"] == y["determinant"] == 0, spath+".sign")
                metrics["max_scan_determinant_difference"] = max(
                    metrics["max_scan_determinant_difference"],
                    abs(x["determinant"]-y["determinant"]))
        counts[model] = {k: dict(v) for k,v in counts[model].items()}
    check(metrics["compared_scan_points"] == 4440, "full original scan coverage")
    check(metrics["compared_branches"] == 80, "full original branch coverage")
    metrics["max_background_relative_difference"] = backgrounds
    return dict(schema="ns_cowling.c20.migration_comparison.v1",
                status="PASS_SCOPED_C19_C20_MIGRATION" if not errors else "FAIL_C19_C20_MIGRATION",
                limits=dict(root_frequency_atol_hz=FREQUENCY_ATOL_HZ,
                            background_rtol=BACKGROUND_RTOL, unchanged_physical_gate=PHYSICAL_GATE),
                baseline=str(baseline), current=str(current),
                assertion_tested="source-locked 40-star migration equivalence; finite original scan/branch coverage",
                metrics=metrics, branch_counts=counts, artifacts=artifacts, errors=errors,
                interpretation="Finite migration check, not an absolute eigenfrequency error bound. Adaptive RHS counts may differ.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("output exists; select a new evidence file")
    result = audit(args.baseline, args.current)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print(json.dumps({k: result[k] for k in ("status", "metrics", "branch_counts", "errors")},
                     indent=2, sort_keys=True, allow_nan=False))
    raise SystemExit(0 if not result["errors"] else 1)


if __name__ == "__main__":
    main()
