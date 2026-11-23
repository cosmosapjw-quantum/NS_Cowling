#!/usr/bin/env python3
"""Audit actual retained crust prefixes in Fortran and qualify native candidates.

Python handles archives, exact density-row identity and receipts only. Physical
unit conversion, uncertainty propagation and cold-first-law bounds are Fortran.
A positive finite-table bound is necessary, not a continuous first-law proof.
Native journals and original numerical archive members remain unchanged.
"""
from __future__ import annotations
import argparse
from collections import Counter
import copy
import csv
import io
import json
from pathlib import Path
import re
import subprocess
import time

from adjudicate_eos_primary_crust_20261123 import (
    ROOT, INPUTS, need, sha, save, locate, qualify_archive, write_csv)
from extract_eos_reaudit_candidate_20261123 import validated_members
from prepare_eos_original_paper_sources_20261123 import LIMITATIONS
import tarfile

GOOD="FINITE_TABLE_NECESSARY_BOUNDS_SATISFIED"
BAD="REFUTED_BY_FINITE_TABLE"

def rows_and_mass(csv_text,nml_text):
    match=re.search(r"^\s*mn\s*=\s*([^\s,/]+)",nml_text,re.I|re.M)
    need(match is not None,"Missing source mass convention")
    rows=list(csv.DictReader(io.StringIO(csv_text)))
    need(bool(rows),"Empty physical table")
    return rows,match[1].replace("D","E").replace("d","e")

def write_input(path,rows,mn):
    # Copy serialization fields verbatim. No energy/pressure conversion here.
    with path.open("x") as f:
        f.write("mn "+mn+"\n")
        for r in rows:
            f.write(" ".join(r[k] for k in ("nb","Q1","Q7","Iphase"))+"\n")

def native_audit(exe,directory,rows,mn):
    directory.mkdir(parents=True)
    inp=directory/"prefix.tsv"
    output=directory/"AUDIT.json"
    write_input(inp,rows,mn)
    run=subprocess.run([str(exe),str(inp),str(output)],text=True,capture_output=True,timeout=30)
    (directory/"execution.log").write_text(run.stdout+run.stderr)
    need(run.returncode==0 and output.is_file(),"Source auditor execution failed: "+str(directory))
    receipt=json.loads(output.read_text())
    need(receipt["status"] in (GOOD,BAD),"Source auditor did not produce a finite decision")
    need(receipt["rows"]==len(rows),"Auditor row count mismatch")
    return receipt,dict(input_sha256=sha(inp),audit_sha256=sha(output),audit_path=str(output))

def candidate_prefix(record,source,campaign_identities,input_ledger):
    archive=locate(record["artifact"]["path"])
    need(sha(archive)==record["artifact"]["sha256"],"Candidate hash differs")
    need(archive.stat().st_size==record["artifact"]["bytes"],"Candidate byte count differs")
    with tarfile.open(archive,"r:gz") as tar:
        members={m.name:m for m,_ in validated_members(tar) if m.isfile()}
        required=("internal/eos_data.csv","internal/eos_data.nml")
        need(all(k in members and members[k].size<20000000 for k in required),"Missing archive table")
        text={k:tar.extractfile(members[k]).read().decode() for k in required}
    rows,mn=rows_and_mass(text[required[0]],text[required[1]])
    report=record["native_report"]
    if record["method"]=="maxwell":
        count=report["low_phase_rows"]
        endpoint=report["n_low"]
    else:
        count=report["i1"]+1
        endpoint=report["n1"]
    need(type(count) is int and 2<=count<=len(rows),"Retained prefix has insufficient or invalid rows")
    prefix=rows[:count]
    need(float(prefix[-1]["nb"])==endpoint,"Native crust endpoint and archived rows differ")
    need(all(float(r["nb"])<=endpoint for r in prefix),"Archived crust prefix exceeds interface")
    view=next(v for v in source["matching_views"] if v["id"]==record["crust_view_id"])
    srcdir=locate(view["internal_path"])
    for name in ("eos_data.csv","eos_data.nml"):
        path=str((srcdir/name).resolve())
        need(path in campaign_identities,"Component absent from frozen campaign identities")
        need(sha(path)==campaign_identities[path],"Matching component changed since campaign")
        input_ledger[path]=campaign_identities[path]
    src=list(csv.DictReader(io.StringIO((srcdir/"eos_data.csv").read_text())))
    if record["method"]=="maxwell":
        # Native Maxwell trims by its chemical-potential tolerance. The
        # reported count defines the exact retained prefix; a new density
        # threshold here would create a conflicting endpoint policy.
        originals=src[:count-1]
        need([float(r["nb"]) for r in prefix[:-1]]==[float(r["nb"]) for r in originals],
             "Retained Maxwell source density knots differ")
    else:
        need([float(r["nb"]) for r in prefix]==[float(r["nb"]) for r in src[:count]],
             "Retained bridge source density knots differ")
    return prefix,mn

def retained_extent(record,source,prefix):
    view=next(v for v in source["matching_views"] if v["id"]==record["crust_view_id"])
    original_rows=len(prefix)-(1 if record["method"]=="maxwell" else 0)
    if record["method"]=="maxwell":
        phases={view["phase_id"]}
    else:
        phases={int(r["Iphase"]) for r in prefix}
    source_phases={c["phase_id"] for c in source["components"]}
    return dict(retained_original_rows=original_rows,
        original_source_rows=source["full_selected_rows"],
        original_row_range_zero_based=[view["start"],view["start"]+original_rows],
        retained_original_phase_ids=sorted(phases),source_phase_count=source["source_phase_count"],
        all_declared_crust_phases_present=phases==source_phases,
        endpoint_chain_label_is_not_full_density_domain_admission=True)

def record_view_starts_at_zero(record):
    return record["source_view_bounds"]["crust"][0]==0

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--campaign",required=True,type=Path)
    p.add_argument("--catalog",required=True,type=Path)
    p.add_argument("--auditor",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--private-artifacts",required=True,type=Path)
    p.add_argument("--audit-root",required=True,type=Path)
    a=p.parse_args()
    campaign=a.campaign.resolve();out=a.out.resolve();private=a.private_artifacts.resolve()
    auditroot=a.audit_root.resolve();exe=a.auditor.resolve()
    need(not private.is_relative_to(ROOT) and not private.is_relative_to(out) and not out.is_relative_to(private),
         "Private archives must stay outside Git and public outputs")
    for target in (out,private,auditroot):
        need(not target.exists(),"Refusing to overwrite "+str(target))
    before={str(campaign/f):sha(campaign/f) for f in INPUTS}
    before.update({str(a.catalog.resolve()):sha(a.catalog),str(exe):sha(exe),
        str(Path(__file__).resolve()):sha(__file__)})
    catalog=json.loads(a.catalog.read_text())
    sources={s["source_id"]:s for s in catalog["sources"]}
    initial=json.loads((campaign/"SUMMARY.json").read_text())
    need(initial["status"]=="COMPLETED" and initial["input_and_binary_hashes_unchanged"],
         "Native campaign incomplete or input identities changed")
    settings=json.loads((campaign/"SETTINGS.json").read_text())
    need(sha(a.catalog)==settings["catalog_sha256"]==initial["settings"]["catalog_sha256"],
         "Adjudication catalog is not the frozen campaign catalog")
    need(settings==initial["settings"],"Frozen campaign settings differ from completion summary")
    campaign_identities=json.loads((campaign/"INPUT_IDENTITIES.json").read_text())
    records=[json.loads(l) for l in (campaign/"RESULTS.jsonl").read_text().splitlines()]
    need(len(records)==1344 and len({r["cell_id"] for r in records})==1344,"Cell census differs")
    need(Counter(r["status"] for r in records)==Counter(initial["statuses"]),"Status census differs")
    started=time.perf_counter()
    out.mkdir(parents=True);(out/"artifacts").mkdir();private.mkdir(parents=True)
    auditroot.mkdir(parents=True)
    source_audits={}
    for s in sources.values():
        if s["role"]!="crust":
            continue
        receipt=locate(s["generation_receipt_path"])
        need(sha(receipt)==s["generation_receipt_sha256"],"Generation receipt changed")
        dispatch=json.loads(receipt.read_text())
        generated=locate(s["compose_path"])
        for name,digest in dispatch["evidence"]["output_sha256"].items():
            need(sha(generated/name)==digest,"Generated member changed: "+name)
        directory=locate(s["internal_path"])
        for name,digest in s["full_selected_sha256"].items():
            need(sha(directory/name)==digest,"Frozen source changed")
            before[str((directory/name).resolve())]=digest
        rows,mn=rows_and_mass((directory/"eos_data.csv").read_text(),(directory/"eos_data.nml").read_text())
        result,identity=native_audit(exe,auditroot/"sources"/s["model"],rows,mn)
        source_audits[s["model"]]=dict(status=result["status"],rows=result["rows"],
            pressure_decrease_count=result["pressure_decrease_count"],
            firstlaw_violation_count=result["firstlaw_violation_count"],
            first_bad_interval_one_based=result["first_bad_interval_one_based"],**identity)
    adjudicated=[];audit_metadata=[];overturned=0
    for original in records:
        r=copy.deepcopy(original);r["native_status"]=r["status"]
        r["source_accuracy_qualification"]=LIMITATIONS[r["crust_model"]]
        if original["status"]=="SURVIVED":
            source=sources[r["crust_source_id"]]
            prefix,mn=candidate_prefix(r,source,campaign_identities,before)
            r["retained_source_extent"]=retained_extent(r,source,prefix)
            result,identity=native_audit(exe,auditroot/("private" if r["private"] else "public")/r["cell_id"],prefix,mn)
            r["source_admission"]=dict(status=result["status"],retained_rows=len(prefix),**identity)
            # This public metadata ledger contains no private endpoints or values.
            audit_metadata.append(dict(cell_id=r["cell_id"],private=r["private"],
                status=result["status"],audit_sha256=identity["audit_sha256"],
                input_sha256=identity["input_sha256"]))
            if result["status"]==BAD:
                r["status"]="REJECTED_RETAINED_SOURCE_THERMODYNAMICS"
                overturned+=1
                whole_first_interval=(r["crust_model"]=="BBP1971_original_published"
                    and source_audits[r["crust_model"]]["first_bad_interval_one_based"]==1
                    and len(source["matching_views"])==1
                    and source["matching_views"][0]["start"]==0
                    and record_view_starts_at_zero(r))
                r["all_allowed_source_prefixes_refuted"]=whole_first_interval
                r["source_rejection_scope"]=("all_allowed_prefixes_of_frozen_single_view_source"
                    if whole_first_interval else "saved_native_candidate_only; alternate_candidates_not_researched")
                if whole_first_interval:
                    r["retention_proof"]=dict(first_refuted_source_interval_one_based=1,
                        minimum_original_rows=2 if r["method"]=="maxwell" else 3,
                        frozen_kernel_sha256=settings["binary_sha256"],
                        excludes_original_microphysical_no_go=True)
            else:
                scope=dict(scope="native_interface_and_declared_finite_thermodynamic_representation",
                    source_interval_test="necessary_bounds_only_not_continuous_sufficiency",
                    retained_prefix_audit_sha256=identity["audit_sha256"],
                    complete_star_validated=False,global_microphysical_minimum_proven=False)
                r["admission_scope"]=scope
                target=(private if r["private"] else out/"artifacts")/(r["cell_id"]+".tar.gz")
                r["artifact"]=qualify_archive(original,target,r["source_accuracy_qualification"],scope)
        adjudicated.append(r)
    survivors=[r for r in adjudicated if r["status"]=="SURVIVED"]
    summary=copy.deepcopy(initial)
    summary.update(schema="ns_cowling.original_paper_adjudicated_campaign.v1",
        statuses=dict(Counter(r["status"] for r in adjudicated)),survivors=len(survivors),
        native_survivors=initial["survivors"],overturned_native_survivors=overturned,
        native_campaign=str(campaign),adjudication_seconds=time.perf_counter()-started,
        selection_rule="frozen native matching gates plus Fortran audit of actual retained crust rows",
        scientific_scope="representation-specific finite-source candidates; source accuracy, finite-search limits and original HF non-reproduction remain explicit",
        negative_result_scope=dict(native_failure_scope=initial.get("negative_result_scope"),
            source_overturn_scope="saved native candidate refuted; all-prefix rejection only with first-interval/retention certificate; otherwise alternate candidates remain unsearched"),
        all_allowed_source_prefixes_refuted_cells=sum(r.get("all_allowed_source_prefixes_refuted",False) for r in adjudicated),
        candidate_only_source_rejections=sum(r["status"]=="REJECTED_RETAINED_SOURCE_THERMODYNAMICS" and not r.get("all_allowed_source_prefixes_refuted",False) for r in adjudicated),
        nv1973_semantic_duplicate_of="NV1974_author_fit",
        semantic_duplicate_survivors=sum(r["crust_model"]=="NV1973_author_fit" for r in survivors))
    summary["by_method"]={}
    for method in initial["by_method"]:
        group=[r for r in adjudicated if r["method"]==method]
        summary["by_method"][method]=dict(requested=len(group),statuses=dict(Counter(r["status"] for r in group)),
            survivor_coverage=dict(Counter(r["source_coverage"] for r in group if r["status"]=="SURVIVED")))
    after={path:sha(path) for path in before}
    need(before==after,"Adjudication input changed")
    save(out/"SUMMARY.json",summary)
    save(out/"SETTINGS.json",json.loads((campaign/"SETTINGS.json").read_text()))
    save(out/"INPUT_IDENTITIES.json",json.loads((campaign/"INPUT_IDENTITIES.json").read_text()))
    (out/"RESULTS.jsonl").write_text("".join(json.dumps(r,ensure_ascii=False,allow_nan=False)+"\n" for r in adjudicated))
    save(out/"SURVIVORS.json",survivors)
    save(out/"SOURCE_ADJUDICATION.json",dict(schema="ns_cowling.original_paper_source_admission.v1",
        auditor_sha256=sha(exe),source_audits=source_audits,retained_prefix_audit_metadata=audit_metadata,
        inputs_before=before,inputs_after=after,kernel_journal_modified=False,source_rows_modified=False,
        positive_result_is_only_necessary_condition=True,original_microphysical_no_go_asserted=False))
    counts=[]
    for model in sorted({r["crust_model"] for r in adjudicated}):
        for method in initial["by_method"]:
            group=[r for r in adjudicated if r["crust_model"]==model and r["method"]==method]
            pos=[r for r in group if r["status"]=="SURVIVED"]
            counts.append(dict(model=model,method=method,requested=len(group),
                native_survivors=sum(r["native_status"]=="SURVIVED" for r in group),
                source_rejected=sum(r["status"]=="REJECTED_RETAINED_SOURCE_THERMODYNAMICS" for r in group),
                representation_survivors=len(pos),public_survivors=sum(not r["private"] for r in pos),
                private_survivors=sum(r["private"] for r in pos),
                endpoint_chain=sum(r["source_coverage"]=="SELECTED_SOURCE_ENDPOINT_CHAIN" for r in pos),
                finite_component=sum(r["source_coverage"]=="FINITE_PHASE_COMPONENT_ONLY" for r in pos),
                all_crust_phases_present=sum(r["retained_source_extent"]["all_declared_crust_phases_present"] for r in pos),
                semantic_duplicate_of="NV1974_author_fit" if model=="NV1973_author_fit" else ""))
    write_csv(out/"MODEL_METHOD_COUNTS.csv",list(counts[0]),counts)
    fields=["cell_id","crust_model","core_model","crust_source_id","core_source_id","method",
            "private","source_coverage","crust_view_id","core_view_id","source_accuracy_qualification"]
    write_csv(out/"SURVIVOR_CATALOG.csv",fields,survivors)
    print(json.dumps(dict(out=str(out),statuses=summary["statuses"],survivors=len(survivors),
        overturned_native_survivors=overturned,public=sum(not r["private"] for r in survivors),
        private=sum(r["private"] for r in survivors),source_audits=source_audits),ensure_ascii=False))

if __name__=="__main__":
    main()
