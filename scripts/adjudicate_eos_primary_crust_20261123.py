#!/usr/bin/env python3
"""Apply a frozen native source-first-law certificate to actual matched rows.

The native campaign is immutable. This creates a separate, fully enumerated
delivery view. Python performs receipts, archive identity and serialization
checks only. The thermodynamic impossibility bound is computed in Fortran by
compose_eos_bbp_primary and independently reproduced by the native audit probe.
No source values, interpolation choices or physical tolerances are changed.
"""
from __future__ import annotations
import argparse
from collections import Counter
import copy
import csv
import hashlib
import gzip
import io
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import tarfile
import time

ROOT=Path(__file__).resolve().parents[1]
RAW_MODEL="BBP1971_table5"
RAW_SOURCE_SHA="5a2711d3230c4aebb5589818c1f9843582b58118f97fff4b505fc755e58a72f8"
QUALIFICATIONS={
 "BPS1971_sequence":"restricted12nuclei; outer_crust_only; retained_phase_extent_is_candidate_specific",
 "BBP1971_pressure_integral":"derived_pressure_barotrope; absolute_energy_anchor_and_near_drip_neutron_gas_interpretation_limited",
 "NV1974_author_fit":"author_energy_fit; near_drip_mu_below_neutron_rest_energy; originalHF_not_recomputed",
}
INPUTS=("SUMMARY.json","RESULTS.jsonl","SURVIVORS.json","SETTINGS.json","INPUT_IDENTITIES.json")

def sha(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as f:
  for block in iter(lambda:f.read(1<<20),b""): h.update(block)
 return h.hexdigest()

def need(ok,message):
 if not ok: raise ValueError(message)

def save(path,value):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def locate(path):
 p=Path(path);return p if p.is_absolute() else ROOT/p

def physical_prefix(csv_text,nml_text):
 # Independent serialization identity check, not an EOS reconstruction.
 match=re.search(r"^\s*mn\s*=\s*([^\s,/]+)",nml_text,re.I|re.M)
 need(match is not None,"Missing source mass convention")
 mn=float(match[1].replace("D","E").replace("d","e"))
 rows=list(itertools.islice(csv.DictReader(io.StringIO(csv_text)),2))
 need(len(rows)==2,"Matched output does not retain two source rows")
 prefix=[]
 for r in rows:
  n=float(r["nb"]);p=n*float(r["Q1"]);eps=n*mn*(1+float(r["Q7"]))
  need(all(math.isfinite(x) for x in (n,p,eps)) and min(n,p,eps)>0,"Invalid prefix serialization")
  prefix.append((n,p,eps))
 return prefix

def candidate_prefix(record):
 archive=locate(record["artifact"]["path"])
 need(sha(archive)==record["artifact"]["sha256"],"Candidate archive hash differs")
 need(archive.stat().st_size==record["artifact"]["bytes"],"Candidate byte count differs")
 with tarfile.open(archive,"r:gz") as tar:
  content={}
  for name in ("internal/eos_data.csv","internal/eos_data.nml"):
   members=[m for m in tar.getmembers() if m.name==name]
   need(len(members)==1 and members[0].isfile() and members[0].size<20000000,
        "Invalid or duplicate archive evidence member")
   content[name]=tar.extractfile(members[0]).read().decode()
 return physical_prefix(content["internal/eos_data.csv"],content["internal/eos_data.nml"])

def qualify_archive(original,target,qualification,admission_scope):
 from extract_eos_reaudit_candidate_20261123 import validated_members
 archive=locate(original["artifact"]["path"])
 need(sha(archive)==original["artifact"]["sha256"],"Original numerical archive changed")
 with tarfile.open(archive,"r:gz") as old:
  payload={}
  for member,_ in validated_members(old):
   if member.isdir(): continue
   need(member.size<20000000,"Unexpectedly large candidate member")
   payload[member.name]=old.extractfile(member).read()
 need("SOURCE_ADMISSION.json" not in payload,"Candidate already has source qualification")
 original_hashes={name:hashlib.sha256(content).hexdigest() for name,content in payload.items()}
 sidecar=dict(schema="ns_cowling.qualified_matched_representation.v1",
   cell_id=original["cell_id"],crust_source_id=original["crust_source_id"],
   core_source_id=original["core_source_id"],private=original["private"],
   native_status=original["status"],final_status="REPRESENTATION_SPECIFIC_CANDIDATE",
   admission_scope=admission_scope,source_accuracy_qualification=qualification,
   original_microphysics_reproduced=False,complete_star_validated=False,
   native_original_archive_sha256=original["artifact"]["sha256"],
   original_members_sha256=original_hashes,original_member_bytes_unchanged=True)
 payload["SOURCE_ADMISSION.json"]=(json.dumps(sidecar,indent=2,ensure_ascii=False,allow_nan=False)+"\n").encode()
 with target.open("xb") as raw:
  with gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0,compresslevel=6) as zipped:
   with tarfile.open(fileobj=zipped,mode="w",format=tarfile.PAX_FORMAT) as new:
    for name,content in sorted(payload.items()):
     info=tarfile.TarInfo(name);info.size=len(content);info.mode=0o644
     info.uid=info.gid=0;info.uname=info.gname="";info.mtime=0
     new.addfile(info,io.BytesIO(content))
 with tarfile.open(target,"r:gz") as check:
  back={m.name:hashlib.sha256(check.extractfile(m).read()).hexdigest()
        for m,_ in validated_members(check) if m.isfile()}
 expected={name:hashlib.sha256(content).hexdigest() for name,content in payload.items()}
 need(back==expected,"Qualified archive changed original bytes or sidecar")
 need(sha(archive)==original["artifact"]["sha256"],"Original archive changed during wrapping")
 return dict(path=str(target),sha256=sha(target),bytes=target.stat().st_size,
             native_original_sha256=original["artifact"]["sha256"],
             qualification_sha256=expected["SOURCE_ADMISSION.json"],
             original_members_unchanged=True)

def write_csv(path,fields,records):
 with Path(path).open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore");w.writeheader();w.writerows(records)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--campaign",required=True,type=Path)
 p.add_argument("--catalog",required=True,type=Path)
 p.add_argument("--out",required=True,type=Path)
 p.add_argument("--private-artifacts",required=True,type=Path)
 a=p.parse_args();campaign=a.campaign.resolve();out=a.out.resolve()
 private=a.private_artifacts.resolve()
 need(not private.is_relative_to(ROOT) and not private.is_relative_to(out) and not out.is_relative_to(private),
      "Private qualified artifacts must stay outside Git and public outputs")
 need(not private.exists(),"Refusing to overwrite private qualified archives")
 need(not out.exists(),"Refusing to overwrite adjudicated results")
 before={str(campaign/f):sha(campaign/f) for f in INPUTS}
 before[str(a.catalog.resolve())]=sha(a.catalog)
 before[str(Path(__file__).resolve())]=sha(Path(__file__).resolve())
 catalog=json.loads(a.catalog.read_text())
 sources={s["source_id"]:s for s in catalog["sources"]}
 raw_source=next(s for s in sources.values() if s["model"]==RAW_MODEL and s["role"]=="crust")
 source_file=locate(raw_source["source_contract"]["source_path"])
 need(sha(source_file)==RAW_SOURCE_SHA,"Primary BBP source identity differs from proved source")
 dispatch_path=locate(raw_source["generation_receipt_path"])
 need(sha(dispatch_path)==raw_source["generation_receipt_sha256"],"Native generation receipt changed")
 dispatch=json.loads(dispatch_path.read_text())
 diag_path=locate(raw_source["compose_path"])/"primary_source_diagnostics.json"
 need(sha(diag_path)==dispatch["evidence"]["output_sha256"]["primary_source_diagnostics.json"],
      "Native source certificate is not sealed")
 certificate=json.loads(diag_path.read_text())
 need(certificate["model"]==RAW_MODEL and certificate["stable_firstlaw_within_printed_rounding"]==
      "REFUTED_FOR_SELECTED_SOURCE" and certificate["first_incompatible_interval_one_based"]==1,
      "Native source certificate does not establish the required first interval")
 need(certificate["first_interval_delta_e_bounds_mev"][0]>
      certificate["first_interval_pressure_integral_bounds_mev"][1],
      "Native firstlaw intervals do not prove incompatibility")
 source_dir=locate(raw_source["internal_path"])
 source_prefix=physical_prefix((source_dir/"eos_data.csv").read_text(),(source_dir/"eos_data.nml").read_text())
 for path in (source_file,dispatch_path,diag_path,source_dir/"eos_data.csv",source_dir/"eos_data.nml"):
  before[str(path.resolve())]=sha(path)
 initial=json.loads((campaign/"SUMMARY.json").read_text())
 need(initial["status"]=="COMPLETED" and initial["input_and_binary_hashes_unchanged"],
      "Native campaign is not complete and immutable")
 records=[json.loads(l) for l in (campaign/"RESULTS.jsonl").read_text().splitlines()]
 need(len(records)==896 and len({r["cell_id"] for r in records})==896,"Wrong or duplicate cell census")
 need(initial["recorded_cells"]==len(records),"Native census mismatch")
 need(Counter(r["status"] for r in records)==Counter(initial["statuses"]),"Native status census mismatch")
 started=time.perf_counter();adjudicated=[];prefix_evidence=[]
 out.mkdir(parents=True);(out/"artifacts").mkdir();private.mkdir(parents=True)
 for original in records:
  record=copy.deepcopy(original)
  record["native_status"]=original["status"]
  if record["crust_model"]==RAW_MODEL:
   need(record["crust_view_id"]=="full_selected" if original["status"]=="SURVIVED" else True,
        "Unexpected raw BBP phase projection")
   if original["status"]=="SURVIVED":
    need(record["source_view_bounds"]["crust"][0]==0,"Raw BBP prefix was removed")
    actual=candidate_prefix(record);errors=[]
    for observed,expected in zip(actual,source_prefix):
     need(observed[0]==expected[0],"Original source density knot was not retained")
     row_errors=[abs(x-y)/max(abs(x),abs(y)) for x,y in zip(observed,expected)]
     need(max(row_errors)<=64*2.220446049250313e-16,"Retained BBP energy/pressure prefix differs")
     errors.append(row_errors)
    prefix_evidence.append(dict(cell_id=record["cell_id"],private=record["private"],
                                archive_sha256=record["artifact"]["sha256"],
                                max_relative_serialization_difference=max(map(max,errors))))
   record["status"]="REJECTED_SOURCE_FIRST_LAW"
   record["source_admission"]={"certificate":str(diag_path.relative_to(ROOT)),
       "source_sha256":RAW_SOURCE_SHA,"retained_required_prefix_rows":[1,2],
       "scope":"stable cold firstlaw impossible inside the specified printing boxes; not a no-go for BBP microphysics"}
  elif original["status"]=="SURVIVED":
   record["source_accuracy_qualification"]=QUALIFICATIONS[record["crust_model"]]
   record["admission_scope"]="interface_and_declared_thermodynamic_representation; not_full_microphysical_or_stellar_admission"
   target=(private if record["private"] else out/"artifacts")/(record["cell_id"]+".tar.gz")
   record["artifact"]=qualify_archive(original,target,record["source_accuracy_qualification"],
                                     record["admission_scope"])
  adjudicated.append(record)
 survivors=[r for r in adjudicated if r["status"]=="SURVIVED"]
 settings=json.loads((campaign/"SETTINGS.json").read_text())
 summary=copy.deepcopy(initial)
 summary.update(schema="ns_cowling.primary_crust_adjudicated_campaign.v1",
   statuses=dict(Counter(r["status"] for r in adjudicated)),survivors=len(survivors),
   native_survivors=initial["survivors"],native_campaign=str(campaign),
   native_campaign_inputs_sha256={f:before[str(campaign/f)] for f in INPUTS},
   adjudication_seconds=time.perf_counter()-started,
   overturned_native_survivors=len(prefix_evidence),
   primary_target_cells=672,derived_sensitivity_cells=224,
   selection_rule="native interface/readback admission plus immutable original-prefix source certificate",
   scientific_scope="representation-specific matched candidates; source precision and neutron-gas interpretation remain qualified")
 by_method={}
 for method in initial["by_method"]:
  group=[r for r in adjudicated if r["method"]==method]
  pos=[r for r in group if r["status"]=="SURVIVED"]
  by_method[method]=dict(requested=len(group),statuses=dict(Counter(r["status"] for r in group)),
                        survivor_coverage=dict(Counter(r["source_coverage"] for r in pos)))
 summary["by_method"]=by_method
 before_after={path:sha(path) for path in before}
 need(before==before_after,"Adjudication input changed")
 save(out/"SUMMARY.json",summary)
 save(out/"SETTINGS.json",settings)
 old_inputs=json.loads((campaign/"INPUT_IDENTITIES.json").read_text())
 save(out/"INPUT_IDENTITIES.json",old_inputs)
 (out/"RESULTS.jsonl").write_text("".join(json.dumps(r,ensure_ascii=False,allow_nan=False)+"\n" for r in adjudicated))
 save(out/"SURVIVORS.json",survivors)
 save(out/"SOURCE_ADJUDICATION.json",dict(schema="ns_cowling.primary_source_adjudication.v1",
      status="VERIFIED_RETAINED_PREFIX_SOURCE_REJECTION",native_certificate=certificate,
      source_certificate_sha256=sha(diag_path),inputs_before=before,inputs_after=before_after,
      retained_prefix_checks=prefix_evidence,source_rows_modified=False,kernel_journal_modified=False,
      original_microphysical_no_go_asserted=False))
 method_rows=[]
 for model in sorted({r["crust_model"] for r in adjudicated}):
  for method in initial["by_method"]:
   group=[r for r in adjudicated if r["crust_model"]==model and r["method"]==method]
   pos=[r for r in group if r["status"]=="SURVIVED"]
   method_rows.append(dict(model=model,method=method,requested=len(group),
       native_survivors=sum(r["native_status"]=="SURVIVED" for r in group),
       source_rejected=sum(r["status"]=="REJECTED_SOURCE_FIRST_LAW" for r in group),
       representation_survivors=len(pos),public_survivors=sum(not r["private"] for r in pos),
       private_survivors=sum(r["private"] for r in pos),
       endpoint_chain=sum(r["source_coverage"]=="SELECTED_SOURCE_ENDPOINT_CHAIN" for r in pos),
       finite_component=sum(r["source_coverage"]=="FINITE_PHASE_COMPONENT_ONLY" for r in pos)))
 write_csv(out/"MODEL_METHOD_COUNTS.csv",list(method_rows[0]),method_rows)
 fields=["cell_id","crust_model","core_model","crust_source_id","core_source_id","method",
         "private","source_coverage","crust_view_id","core_view_id","source_accuracy_qualification"]
 write_csv(out/"SURVIVOR_CATALOG.csv",fields,survivors)
 save(out/"INNER_SOURCE_QUALIFICATIONS.json",{
     "scope":"scientific accuracy warning, not a modification of raw kernel outcomes",
     "BPS1971_sequence":"original published candidate sequence restricted, outer-only",
     "BBP1971_pressure_integral":"first printed energy anchor gives mu=938.876804478 MeV below modern neutron rest energy; full neutron-gas microphysics not established",
     "NV1974_author_fit":"near-drip author fit mu below neutron rest energy, original HF table discrepancy and derivative-fit noise",
     "candidate_ids_requiring_inner_source_interpretation":[r["cell_id"] for r in survivors if r["crust_model"]!="BPS1971_sequence"]})
 print(json.dumps({"out":str(out),"statuses":summary["statuses"],"survivors":len(survivors),
                   "overturned_native_survivors":len(prefix_evidence),
                   "public":sum(not r["private"] for r in survivors),
                   "private":sum(r["private"] for r in survivors)},ensure_ascii=False))

if __name__=="__main__": main()
