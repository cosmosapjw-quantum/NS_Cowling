#!/usr/bin/env python3
"""Measure serial/process scheduling with identical native EOS work and verify bytes."""
from __future__ import annotations
import argparse, hashlib, json, tarfile, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import run_eos_all_pairs_20261123 as runner

def eos_hashes(result):
    if result["status"] != "SURVIVED":
        return None
    with tarfile.open(result["artifact"]["path"], "r:gz") as archive:
        return {name:hashlib.sha256(archive.extractfile(name).read()).hexdigest()
                for name in ("internal/eos_data.csv", "internal/eos_data.nml")}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--native-exe", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args=parser.parse_args()
    sources=json.loads(args.catalog.read_text())["sources"]
    runner.validate_catalog(sources)
    all_cells=runner.enumerate_cells(sources)
    pairs=[all_cells[i:i+7] for i in range(0,len(all_cells),7)]
    indices=sorted(set(round(i*(len(pairs)-1)/13) for i in range(14)))
    selected={c[0]:c for i in indices for c in pairs[i]}
    survivors=json.loads((args.campaign/"SURVIVORS.json").read_text())
    sharp=sorted((r for r in survivors if r["method"]=="maxwell"),
                 key=lambda r:r["cell_id"])
    extra_pairs={(r["crust_source_id"],r["core_source_id"]) for r in (sharp[0],sharp[-1])}
    for cell in all_cells:
        if (cell[1]["source_id"],cell[2]["source_id"]) in extra_pairs:
            selected[cell[0]]=cell
    cells=[selected[k] for k in sorted(selected)]
    if args.out.exists():
        raise FileExistsError("Benchmark output must be fresh")
    args.out.mkdir(parents=True)
    before=runner.input_identities(sources,args.native_exe,args.catalog)
    runs={}
    for workers in (1,8):
        public=args.out/f"workers_{workers}"
        private=args.out/f"private_workers_{workers}"
        work=args.out/f"work_{workers}"
        for path in (public,private,work):path.mkdir()
        runner.validate_private_location(private,public)
        settings=dict(exe=str(args.native_exe.resolve()),out=str(public.resolve()),
            private_artifacts=str(private.resolve()),work=str(work.resolve()),
            n_bridge=80,timeout=180,workers=workers,
            catalog_sha256=runner.digest(args.catalog),binary_sha256=runner.digest(args.native_exe))
        start=time.perf_counter()
        with ProcessPoolExecutor(max_workers=workers,initializer=runner.init_worker,
                                 initargs=(settings,)) as pool:
            results=list(pool.map(runner.run_cell,cells))
        elapsed=time.perf_counter()-start
        normalized={}
        for result in results:
            normalized[result["cell_id"]]={
                "status":result["status"],"method":result["method"],
                "native_candidate_attempts":result["native_candidate_attempts"],
                "crust_view_id":result.get("crust_view_id"),"core_view_id":result.get("core_view_id"),
                "coverage":result.get("source_coverage"),
                "native_report":result.get("native_report"),
                "failure_categories":result.get("failure_categories"),
                "eos_file_hashes":eos_hashes(result),
            }
        runner.save(public/"RESULTS.json",results)
        runner.save(public/"NORMALIZED.json",normalized)
        runs[workers]=dict(seconds=elapsed,normalized=normalized)
        print(json.dumps({"workers":workers,"seconds":elapsed,"cells":len(cells)}),flush=True)
    unchanged=runner.input_identities(sources,args.native_exe,args.catalog)==before
    exact=runs[1]["normalized"]==runs[8]["normalized"]
    receipt=dict(schema="ns_cowling.native_parallel_equivalence.v1",
        cells=len(cells),pair_count=len(cells)//7,methods=list(runner.METHODS),
        sample_rule="14 evenly spaced source pairs plus first/last genuinely surviving Maxwell pairs",
        one_worker_seconds=runs[1]["seconds"],eight_worker_seconds=runs[8]["seconds"],
        measured_speedup=runs[1]["seconds"]/runs[8]["seconds"],
        exact_status_endpoint_and_eos_csv_nml_byte_equivalence=exact,
        input_and_binary_hashes_unchanged=unchanged,
        success_cells=sum(r["status"]=="SURVIVED" for r in runs[1]["normalized"].values()),
        all_eos_hashes={k:v["eos_file_hashes"] for k,v in runs[1]["normalized"].items()},
        limitation="one paired benchmark on this host; not a universal speedup or equal-archive-byte claim")
    runner.save(args.out/"PARALLEL_EQUIVALENCE.json",receipt)
    print(json.dumps({k:v for k,v in receipt.items() if k!="all_eos_hashes"}),flush=True)
    if not exact or not unchanged:
        raise RuntimeError("Parallel numerical equivalence or identity failed")
if __name__=="__main__":main()
