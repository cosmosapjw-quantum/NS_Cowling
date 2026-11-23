#!/usr/bin/env python3
"""Publish complete public receipts and only metadata of private combinations."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
from publish_eos_primary_crust_projection_20261123 import campaign_view, sha, save

def publish_failure_traces(directory):
    records=[json.loads(line) for line in (directory/"RESULTS.jsonl").read_text().splitlines()]
    trace_dir=directory/"failure_traces"
    trace_dir.mkdir()
    trace_manifest={}
    for r in records:
        if r["private"] or "failure_trace" not in r:
            continue
        original=Path(r["failure_trace"]["path"])
        if sha(original)!=r["failure_trace"]["sha256"]:
            raise ValueError("Public failure trace identity changed")
        target=trace_dir/(r["cell_id"]+".json.gz")
        os.link(original,target)
        if sha(target)!=r["failure_trace"]["sha256"]:
            raise ValueError("Public failure trace copy changed")
        r["failure_trace"]["path"]=str(target.resolve())
        trace_manifest[target.name]=dict(sha256=sha(target),bytes=target.stat().st_size)
    (directory/"RESULTS.jsonl").write_text("".join(json.dumps(r,ensure_ascii=False,allow_nan=False)+"\n" for r in records))
    save(directory/"SURVIVORS.json",[r for r in records if r["status"]=="SURVIVED"])
    save(directory/"FAILURE_TRACE_INDEX.json",trace_manifest)
    return dict(public_outputs={p.name:sha(p) for p in directory.iterdir() if p.is_file()},
                public_failure_traces=len(trace_manifest))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--catalog",required=True,type=Path)
    p.add_argument("--native-campaign",required=True,type=Path)
    p.add_argument("--adjudicated",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    if a.out.exists():
        raise FileExistsError(a.out)
    a.out.mkdir(parents=True)
    catalog=json.loads(a.catalog.read_text())
    sources=[]
    for s in catalog["sources"]:
        if not s["private"]:
            sources.append(s)
        else:
            sources.append({k:s[k] for k in ("source_id","role","model","private")}|
                dict(available=False,status="PRIVATE_SOURCE_NUMERICAL_CONTENT_WITHHELD",
                     publication="METADATA_ONLY",execution_source_available=True,matching_views=[]))
    full=[json.loads(line) for line in (a.native_campaign/"RESULTS.jsonl").read_text().splitlines()]
    public_cells=sum(not r["private"] for r in full)
    private_cells=sum(r["private"] for r in full)
    if len(full)!=catalog["counts"]["total_cells"] or private_cells!=catalog["counts"]["private_cells"]:
        raise ValueError("Source catalog and executed cell census differ")
    view={k:v for k,v in catalog.items() if k not in ("sources","counts")}
    view.update(schema="ns_cowling.original_paper_source_public_projection.v1",sources=sources,
        full_catalog_sha256=sha(a.catalog),full_execution_counts=catalog["counts"],
        public_projection_counts=dict(public_sources=sum(not s["private"] for s in sources),
            private_metadata_sources=sum(s["private"] for s in sources),
            public_replay_available_cells=public_cells,privately_executed_cells_withheld=private_cells))
    save(a.out/"SOURCE_CATALOG.json",view)
    native=campaign_view(a.native_campaign.resolve(),a.out/"native_campaign",False)
    qualified=campaign_view(a.adjudicated.resolve(),a.out/"qualified_campaign",True)
    native.update(publish_failure_traces(a.out/"native_campaign"))
    qualified.update(publish_failure_traces(a.out/"qualified_campaign"))
    receipt=dict(schema="ns_cowling.original_paper_public_projection.v1",
        source_catalog_full_sha256=sha(a.catalog),source_catalog_public_sha256=sha(a.out/"SOURCE_CATALOG.json"),
        native=native,qualified=qualified,original_full_receipts_modified=False,
        private_policy="Private source tables, endpoints, numerical reports and candidate files withheld; identities, statuses and aggregate counts only")
    save(a.out/"PUBLIC_PROJECTION.json",receipt)
    print(json.dumps(dict(out=str(a.out),qualified=qualified),ensure_ascii=False))

if __name__=="__main__":
    main()
