#!/usr/bin/env python3
"""Read-only aggregation of frozen v2 post-bridge failures; no EOS arithmetic."""
from pathlib import Path
import argparse, collections, gzip, hashlib, json

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = args.root / "data/eos_matching_reaudit_20261123/campaign_v2"
    journal = base / "RESULTS.jsonl"
    assert sha(journal) == "c7b30f71ad74df5f9cfe6882fa84f3ee66c2f24fd7d764ad4a651c88aabbf27c"
    needle = "post_bridge_gate_failure_search_incomplete"
    hits, trials = [], []
    for line in journal.read_text().splitlines():
        r = json.loads(line)
        if needle not in json.dumps(r):
            continue
        hits.append(r)
        trace = Path(r["failure_trace"]["path"])
        if not trace.exists():
            trace = base / "failure_traces" / (r["cell_id"] + ".json.gz")
        assert sha(trace) == r["failure_trace"]["sha256"]
        for t in json.loads(gzip.decompress(trace.read_bytes())):
            if t.get("failure_category") != needle:
                continue
            rec = t.get("details", {}).get("receipt", {})
            match = t.get("initial_matching_receipt", {})
            native = t.get("native_call", {})
            trials.append(dict(cell_id=r["cell_id"], core=r["core_source_id"],
                crust=r["crust_source_id"], method=r["method"],
                crust_view=t.get("crust_view_id"), core_view=t.get("core_view_id"),
                stage=t.get("stage"), native_payload=native.get("payload"),
                failed_checks=rec.get("failed_checks"),
                min_dmu=rec.get("min_dmu_b_dnb"),
                derivative_points=rec.get("derivative_points"),
                masked_points=rec.get("derivative_masked_points"),
                invalid_points=rec.get("invalid_derivative_points"),
                i1=match.get("i1"), i2=match.get("i2"), n1=match.get("n1"),
                pairs_total=match.get("window_pairs_total"),
                pairs_examined=match.get("window_pairs_examined"),
                bridge_rejected=match.get("window_candidates_rejected"),
                bridge_candidates=match.get("window_candidates_constructed")))
    result = dict(cell_count=len(hits), trial_count=len(trials), trials=trials)
    blob = (json.dumps(result, indent=2) + "\n").encode()
    expected = "1c10e880056670e79d7eb1737d1908988eaa1fa6d97c60db88194e450f1b476b"
    assert hashlib.sha256(blob).hexdigest() == expected
    if args.output.exists():
        assert args.output.read_bytes() == blob
    else:
        with args.output.open("xb") as stream:
            stream.write(blob)
    print(json.dumps(dict(cells=len(hits), trials=len(trials), output=str(args.output),
        sha256=expected, failed_checks=dict(collections.Counter(t["failed_checks"] for t in trials)),
        methods=dict(collections.Counter(t["method"] for t in trials)),
        crusts=dict(collections.Counter(t["crust"] for t in trials)))))

if __name__ == "__main__":
    main()
