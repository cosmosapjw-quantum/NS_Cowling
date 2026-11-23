#!/usr/bin/env python3
"""Compare two native matching journals and their archived products without recomputing physics."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import tarfile


def read_json(path):
    return json.loads(path.read_bytes())


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def journal(path):
    result = {}
    for line in (path / "RESULTS.jsonl").read_text().splitlines():
        row = json.loads(line)
        key = row["cell_id"]
        require(key not in result, "duplicate case")
        result[key] = row
    return result


def checked_file(row, field):
    file = Path(row[field]["path"])
    require(digest(file) == row[field]["sha256"], f"{field} hash mismatch")
    return file


def compare_artifacts(left, right, baseline_binary, candidate_binary):
    with tarfile.open(checked_file(left, "artifact")) as a, tarfile.open(checked_file(right, "artifact")) as b:
        names = set(a.getnames())
        require(names == set(b.getnames()) and len(names) == len(a.getmembers()) == len(b.getmembers()),
                "archive member sets differ or contain duplicates")
        require({"internal/eos_data.csv", "internal/eos_data.nml", "compose/eos.nb", "compose/eos.thermo",
                 "NATIVE_MATCH.json", "PROVENANCE.json"} <= names, "incomplete surviving artifact")
        for name in names:
            old, new = a.extractfile(name).read(), b.extractfile(name).read()
            if name == "PROVENANCE.json":
                old, new = json.loads(old), json.loads(new)
                require(old.pop("native_binary_sha256") == baseline_binary, "unrecognized baseline binary")
                require(new.pop("native_binary_sha256") == candidate_binary, "unrecognized candidate binary")
            elif name in {"READBACK_ROWS.json", "READBACK_METADATA.json"}:
                old, new = json.loads(old), json.loads(new)
                # Only the per-process temporary directory can differ.
                for receipt in (old, new):
                    for field in ("actual_path", "expected_path"):
                        receipt[field] = Path(receipt[field]).parts[-2:]
            require(old == new, f"archive contents differ in {name}")
    return sum(n.startswith("compose/") for n in names)


def trace(row):
    with gzip.open(checked_file(row, "failure_trace")) as stream:
        trials = json.load(stream)
    for trial in trials:
        if trial.get("native_call"):
            trial["native_call"].pop("elapsed_seconds", None)
    return trials


def compare_lane(baseline_dir, candidate_dir):
    summary_left, summary_right = [read_json(p / "SUMMARY.json") for p in (baseline_dir, candidate_dir)]
    settings = [s["settings"] for s in (summary_left, summary_right)]
    for s in (summary_left, summary_right):
        require(s["status"] == "COMPLETED" and s["input_and_binary_hashes_unchanged"], "incomplete or changed input")
    for key in ("catalog_sha256", "runner_sha256", "all_methods", "n_bridge", "full_requested_cells", "pilot"):
        require(settings[0][key] == settings[1][key], f"settings mismatch for {key}")
    require(not settings[0]["pilot"], "pilot cannot serve as a full regression")
    binaries = [s["binary_sha256"] for s in settings]
    # Each arm hashes the exact input files both before and after running. The
    # executable itself is the only allowed difference in those manifests.
    identities = [read_json(p / "INPUT_IDENTITIES.json") for p in (baseline_dir, candidate_dir)]
    for i, manifest in enumerate(identities):
        require(manifest.pop(settings[i]["exe"]) == binaries[i], "executable hash mismatch")
    require(identities[0] == identities[1], "source bytes or input paths changed")
    a, b = journal(baseline_dir), journal(candidate_dir)
    require(a.keys() == b.keys() and len(a) == settings[0]["full_requested_cells"], "case coverage mismatch")
    counts = Counter()
    ignore = {"elapsed_seconds", "artifact", "failure_trace"}
    for key in sorted(a):
        old, new = a[key], b[key]
        label = f"{baseline_dir.name}/{key}"
        for rec in (old, new):
            require(rec["status"] in ("SURVIVED", "NO_ADMITTED_MATCH") and rec["native_candidate_attempts"] > 0,
                    f"{label}: no physical native attempt")
        require({k: v for k, v in old.items() if k not in ignore} ==
                {k: v for k, v in new.items() if k not in ignore}, f"{label}: native result changed")
        counts["cases"] += 1
        counts["source_view_attempts"] += old["native_candidate_attempts"]
        require(("artifact" in old) == ("artifact" in new), f"{label}: artifact availability changed")
        require(("failure_trace" in old) == ("failure_trace" in new), f"{label}: trace availability changed")
        if old["status"] == "SURVIVED":
            try:
                counts["identical_compose_members"] += compare_artifacts(old, new, *binaries)
            except ValueError as exc:
                raise ValueError(f"{label}: {exc}") from exc
            counts["identical_survivors"] += 1
        if "failure_trace" in old:
            lhs, rhs = trace(old), trace(new)
            require(lhs == rhs, f"{label}: failed attempt changed")
            counts["identical_failed_attempts"] += len(lhs)
    return {"baseline_binary": binaries[0], "candidate_binary": binaries[1],
            "catalog_sha256": settings[0]["catalog_sha256"], "counts": dict(counts)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("lanes", nargs="+")
    args = parser.parse_args()
    require(len(set(args.lanes)) == len(args.lanes), "duplicate lane")
    reports = {lane: compare_lane(args.baseline / lane, args.candidate / lane) for lane in args.lanes}
    total = Counter()
    for result in reports.values():
        total.update(result["counts"])
    args.output.write_text(json.dumps({"lanes": reports, "totals": dict(total)}, indent=2) + "\n")


if __name__ == "__main__":
    main()
