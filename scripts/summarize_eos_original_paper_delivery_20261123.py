#!/usr/bin/env python3
"""Summarize provenance/coverage and exact serialized-table duplicates; no EOS physics."""
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path
import tarfile
from adjudicate_eos_primary_crust_20261123 import sha,save,write_csv,need
from extract_eos_reaudit_candidate_20261123 import validated_members

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--campaign",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    need(not a.out.exists(),"Summary target already exists")
    a.out.mkdir(parents=True)
    records=[json.loads(l) for l in (a.campaign/"RESULTS.jsonl").read_text().splitlines()]
    survivors=[r for r in records if r["status"]=="SURVIVED"]
    public=[r for r in survivors if not r["private"]]
    groups=defaultdict(list);metadata=[];strict=[];model_rows=[]
    for r in public:
        path=Path(r["artifact"]["path"])
        need(sha(path)==r["artifact"]["sha256"],"Qualified archive changed")
        with tarfile.open(path,"r:gz") as tar:
            members={m.name:m for m,_ in validated_members(tar) if m.isfile()}
            content={name:tar.extractfile(members[name]).read()
                     for name in ("internal/eos_data.csv","internal/eos_data.nml")}
        hashes={name:hashlib.sha256(value).hexdigest() for name,value in content.items()}
        identity=hashlib.sha256(json.dumps(hashes,sort_keys=True,separators=(",",":")).encode()).hexdigest()
        groups[identity].append(r["cell_id"])
        complete=(r["retained_source_extent"]["all_declared_crust_phases_present"]
                  and r["source_coverage"]=="SELECTED_SOURCE_ENDPOINT_CHAIN")
        replay=r["crust_model"]=="NV1973_author_fit"
        row=dict(cell_id=r["cell_id"],crust_model=r["crust_model"],core_model=r["core_model"],
            method=r["method"],source_coverage=r["source_coverage"],
            all_declared_crust_phases_present=r["retained_source_extent"]["all_declared_crust_phases_present"],
            retained_original_crust_rows=r["retained_source_extent"]["retained_original_rows"],
            original_crust_rows=r["retained_source_extent"]["original_source_rows"],
            semantic_duplicate_of="NV1974_author_fit" if replay else "",
            exact_serialized_native_table_sha256=identity,
            phase_complete_endpoint_chain=complete)
        metadata.append(row)
        if complete and not replay:
            strict.append(row)
    for model in sorted({r["crust_model"] for r in records}):
        group=[r for r in records if r["crust_model"]==model]
        pos=[r for r in group if r["status"]=="SURVIVED"]
        complete=[r for r in pos if r["retained_source_extent"]["all_declared_crust_phases_present"]
                  and r["source_coverage"]=="SELECTED_SOURCE_ENDPOINT_CHAIN"]
        model_rows.append(dict(model=model,requested=len(group),
            native_survivors=sum(r["native_status"]=="SURVIVED" for r in group),
            source_rejected=sum(r["status"]=="REJECTED_RETAINED_SOURCE_THERMODYNAMICS" for r in group),
            representation_candidates=len(pos),public=sum(not r["private"] for r in pos),
            private=sum(r["private"] for r in pos),
            all_crust_phases_present=sum(r["retained_source_extent"]["all_declared_crust_phases_present"] for r in pos),
            phase_complete_endpoint_chain=len(complete),phase_complete_public=sum(not r["private"] for r in complete)))
    write_csv(a.out/"PUBLIC_CANDIDATES_WITH_COVERAGE.csv",list(metadata[0]),metadata)
    write_csv(a.out/"PHASE_COMPLETE_NEW_PUBLIC_CANDIDATES.csv",list(metadata[0]),strict)
    write_csv(a.out/"MODEL_COUNTS.csv",list(model_rows[0]),model_rows)
    save(a.out/"EXACT_TABLE_IDENTITY_GROUPS.json",dict(
        identity_definition="SHA256 of canonical JSON containing exact internal/eos_data.csv and eos_data.nml SHA256; same means byte identity, different does not prove physical inequivalence",
        public_candidates=len(public),exact_table_groups=len(groups),
        duplicate_groups={key:ids for key,ids in groups.items() if len(ids)>1},
        all_groups=dict(groups)))
    summary=dict(schema="ns_cowling.original_paper_delivery_summary.v1",requested=len(records),
        statuses=dict(Counter(r["status"] for r in records)),qualified=len(survivors),
        public=len(public),private=len(survivors)-len(public),
        nv1973_semantic_replay=sum(r["crust_model"]=="NV1973_author_fit" for r in survivors),
        public_exact_table_groups=len(groups),phase_complete_new_public_candidates=len(strict),
        phase_complete_new_public_exact_table_groups=len({r["exact_serialized_native_table_sha256"] for r in strict}),
        strict_selection="all declared crust phases retained AND selected source endpoint chain AND not the NV1973 semantic replay; metadata coverage filter, not an additional microphysics or stellar proof",
        source_inputs={name:sha(a.campaign/name) for name in ("RESULTS.jsonl","SURVIVORS.json","SUMMARY.json")},
        models=model_rows)
    save(a.out/"SUMMARY.json",summary)
    print(json.dumps(summary,ensure_ascii=False))

if __name__=="__main__":
    main()
