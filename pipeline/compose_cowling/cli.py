from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
from collections import Counter
from .adapter import build_adapter_csv, sha256
from .intake import prepare_from_matched
from .direct_barotrope import audit_directory

def verify_spectrum(release: Path) -> dict:
    summary = release / "data" / "spectrum_summary_234.csv"
    rows = list(csv.DictReader(summary.open(newline="", encoding="utf-8")))
    if len(rows) != 234:
        raise ValueError(f"expected 234 spectrum rows, found {len(rows)}")
    names = [r["name"] for r in rows]
    if len(set(names)) != 234:
        raise ValueError("spectrum EOS keys are not unique")
    cores = Counter(r["core"] for r in rows)
    if cores != Counter({"BBB_BHF2":78,"RG_SLY4":78,"GM_GM1":78}):
        raise ValueError(f"unexpected core cohort: {dict(cores)}")
    for r in rows:
        if r["mode_status"] != "VALID" or r["projected_native_cell_audit"] != "PASS":
            raise ValueError(f"{r['name']}: mode or native-cell audit is not valid")
        nodes = tuple(int(r[k]) for k in ("f_inner_W_nodes","f_inner_V_nodes","p1_inner_W_nodes","p1_inner_V_nodes"))
        if nodes != (0,0,1,1):
            raise ValueError(f"{r['name']}: wrong f/p1 node identity: {nodes}")
        if float(r["f_hz"]) <= 0 or float(r["p1_hz"]) <= 0:
            raise ValueError(f"{r['name']}: nonpositive eigenfrequency")
    return {"status":"PASS","rows":len(rows),"core_counts":dict(cores),"recomputed_modes":0,
            "meaning":"Existing completed spectrum reused; this command does not integrate any mode."}

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Reproducible joined-EOS to relativistic-Cowling workflow")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare", help="read successful matched CompOSE tables with the bundled EOS generator")
    p.add_argument("--matched-root", type=Path, required=True)
    p.add_argument("--classifier", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--generator-package", type=Path, default=Path(__file__).resolve().parents[2]/"eos_generator"/"compose_eos_pkg")
    a = sub.add_parser("build-adapters", help="build the optional legacy q-linear Hermite solver representation")
    a.add_argument("--prepared-dir", type=Path, required=True)
    a.add_argument("--output-dir", type=Path, required=True)
    a.add_argument("--no-pav-repair", action="store_true")
    d = sub.add_parser("audit-direct", help="audit original P-epsilon EOS knots with the default no-rewrite barotrope")
    d.add_argument("--source-dir", type=Path, default=Path(__file__).resolve().parents[2]/"data"/"eos_original_normalized")
    d.add_argument("--output-dir", type=Path, required=True)
    v = sub.add_parser("verify-existing", help="verify and reuse saved f/p1 results without rerunning")
    v.add_argument("--release", type=Path, default=Path(__file__).resolve().parents[2]/"oscillation_release_v1")
    args = parser.parse_args(argv)
    if args.command == "prepare":
        result = prepare_from_matched(args.matched_root,args.classifier,args.output_dir,args.generator_package)
        print(json.dumps({"counts":result["counts"],"prepared":len(result["entries"])},indent=2))
    elif args.command == "build-adapters":
        args.output_dir.mkdir(parents=True,exist_ok=True)
        receipts = []
        for source in sorted(args.prepared_dir.glob("*.csv")):
            if source.name in {"spectrum_summary_234.csv"}:
                continue
            target = args.output_dir/source.name
            receipts.append(build_adapter_csv(source,target,enforce_monotone_q=not args.no_pav_repair))
        manifest = {"schema":"compose-cowling-adapter-manifest-v1","n_adapters":len(receipts),"receipts":receipts}
        (args.output_dir/"adapter_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
        print(json.dumps({"n_adapters":len(receipts),"pooled_knots":sum(x["metrics"]["pav_pooled_knots"] for x in receipts)},indent=2))
    elif args.command == "audit-direct":
        manifest = audit_directory(args.source_dir,args.output_dir)
        print(json.dumps({"n_tables":manifest["n_tables"],"all_admissible":manifest["all_admissible"],
                          "source_tables_modified":manifest["source_tables_modified"]},indent=2))
    else:
        print(json.dumps(verify_spectrum(args.release),indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
