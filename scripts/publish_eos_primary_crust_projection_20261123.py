#!/usr/bin/env python3
"""Create public receipts with private numerical results strictly withheld."""
from __future__ import annotations
import argparse
from collections import Counter
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil

META=("cell_id","crust_source_id","core_source_id","crust_model","core_model",
      "method","private","status","native_status","source_coverage",
      "crust_view_id","core_view_id","source_accuracy_qualification","admission_scope")

def sha(p):
 return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def record_view(r):
 if not r["private"]: return copy.deepcopy(r)
 return {**{k:r[k] for k in META if k in r},
         "publication":"PRIVATE_METADATA_ONLY_NUMERICAL_RESULTS_WITHHELD"}

def campaign_view(source,out,include_archives):
 out.mkdir(parents=True)
 originals=[json.loads(l) for l in (source/"RESULTS.jsonl").read_text().splitlines()]
 records=[record_view(r) for r in originals]
 if include_archives: (out/"artifacts").mkdir()
 for original,record in zip(originals,records):
  if original["private"]:
   assert set(record)<=set(META)|{"publication"}
   assert not any(k in record for k in ("native_report","verification","artifact","failure_trace","source_view_bounds"))
  elif record["status"]=="SURVIVED" and include_archives:
   source_archive=Path(original["artifact"]["path"])
   assert sha(source_archive)==original["artifact"]["sha256"]
   target=out/"artifacts"/source_archive.name
   os.link(source_archive,target)
   record["artifact"]["path"]=str(target.resolve())
  elif not include_archives:
   record["publication_scope"]="raw native kernel receipt; use qualified campaign for deliverable EOSs"
 summary=json.loads((source/"SUMMARY.json").read_text())
 settings=json.loads((source/"SETTINGS.json").read_text())
 settings["private_artifacts"]="PRIVATE_LOCAL_BACKUP_NOT_PUBLISHED"
 summary["settings"]=settings
 summary["publication_projection"]={
   "full_journal_sha256":sha(source/"RESULTS.jsonl"),
   "full_survivors_sha256":sha(source/"SURVIVORS.json"),
   "full_settings_sha256":sha(source/"SETTINGS.json"),
   "private_numerical_content":"WITHHELD",
   "public_full_records":sum(not r["private"] for r in records),
   "private_metadata_records":sum(r["private"] for r in records)}
 assert Counter(r["status"] for r in records)==Counter(summary["statuses"])
 (out/"RESULTS.jsonl").write_text("".join(json.dumps(r,ensure_ascii=False,allow_nan=False)+"\n" for r in records))
 save(out/"SURVIVORS.json",[r for r in records if r["status"]=="SURVIVED"])
 save(out/"SUMMARY.json",summary);save(out/"SETTINGS.json",settings)
 # This is only path/SHA identity metadata, never private EOS numerical rows.
 shutil.copyfile(source/"INPUT_IDENTITIES.json",out/"INPUT_IDENTITIES.json")
 for name in ("SOURCE_ADJUDICATION.json","MODEL_METHOD_COUNTS.csv","SURVIVOR_CATALOG.csv",
              "INNER_SOURCE_QUALIFICATIONS.json"):
  if (source/name).is_file(): shutil.copyfile(source/name,out/name)
 return {"full_inputs":{n:sha(source/n) for n in ("RESULTS.jsonl","SURVIVORS.json","SUMMARY.json","SETTINGS.json")},
         "public_outputs":{p.name:sha(p) for p in out.iterdir() if p.is_file()},
         "public_records":sum(not r["private"] for r in records),
         "private_metadata_records":sum(r["private"] for r in records),
         "public_survivors":sum(not r["private"] and r["status"]=="SURVIVED" for r in records),
         "private_survivors_metadata":sum(r["private"] and r["status"]=="SURVIVED" for r in records)}

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--catalog",required=True,type=Path)
 p.add_argument("--native-campaign",required=True,type=Path)
 p.add_argument("--adjudicated",required=True,type=Path)
 p.add_argument("--out",required=True,type=Path)
 a=p.parse_args();assert not a.out.exists(),"public projection already exists"
 a.out.mkdir(parents=True)
 catalog=json.loads(a.catalog.read_text())
 sources=[]
 for s in catalog["sources"]:
  if not s["private"]: sources.append(s)
  else:
   sources.append({k:s[k] for k in ("source_id","role","model","private")}|
     {"available":False,"status":"PRIVATE_SOURCE_NUMERICAL_CONTENT_WITHHELD",
      "publication":"METADATA_ONLY","execution_source_available":True,"matching_views":[]})
 public_catalog={k:v for k,v in catalog.items() if k not in ("sources","counts")}
 public_catalog.update(schema="ns_cowling.primary_crust_source_public_projection.v1",
     sources=sources,full_catalog_sha256=sha(a.catalog),
     full_execution_counts=catalog["counts"],
     public_projection_counts={"public_sources":sum(not s["private"] for s in sources),
                               "private_metadata_sources":sum(s["private"] for s in sources),
                               "public_replay_available_cells":868,
                               "privately_executed_cells_withheld":28})
 save(a.out/"SOURCE_CATALOG.json",public_catalog)
 native=campaign_view(a.native_campaign.resolve(),a.out/"native_campaign",False)
 qualified=campaign_view(a.adjudicated.resolve(),a.out/"qualified_campaign",True)
 receipt={"schema":"ns_cowling.primary_crust_public_projection.v1",
          "source_catalog_full_sha256":sha(a.catalog),
          "source_catalog_public_sha256":sha(a.out/"SOURCE_CATALOG.json"),
          "native":native,"qualified":qualified,
          "private_policy":"Private source tables, endpoint values, native numerical reports and candidate files remain in local backup; only identities/status/count metadata is published",
          "original_full_receipts_modified":False}
 save(a.out/"PUBLIC_PROJECTION.json",receipt)
 print(json.dumps({"out":str(a.out),"qualified":qualified},ensure_ascii=False))

if __name__=="__main__": main()
