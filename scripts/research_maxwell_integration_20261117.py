#!/usr/bin/env python3
"""Re-evaluate source-pair Maxwell feasibility after reference reencoding repair."""
from __future__ import annotations
import argparse, collections, hashlib, json, os, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.generation_dispatch import CORE_IDS,SPG_IDS
from match_core_crust import run_pipeline

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--build",type=Path,required=True)
 p.add_argument("--spg-root",type=Path,required=True)
 p.add_argument("--core-root",type=Path,required=True)
 p.add_argument("--out",type=Path,required=True)
 a=p.parse_args()
 out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
 os.environ.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
     COMPOSE_EOS_FORTRAN_ROOT=str(ROOT),
     COMPOSE_EOS_FORTRAN_EXE=str(a.build.resolve()/"compose_eos_fortran"))
 deps=[ROOT/x for x in ("fortran_eos/src/compose_eos_types.f90",
       "fortran_eos/src/compose_eos_phase_transition.f90",
       "fortran_eos/app/compose_eos_fortran.f90","python_frontend/compose_eos/api.py",
       "scripts/match_core_crust.py","scripts/match_report_admission.py",
       "scripts/research_maxwell_integration_20261117.py")]
 before={str(x.relative_to(ROOT)):sha(x)for x in deps}
 r={"started_unix":time.time(),"scope":"conditional native source pairs; feasibility within provided support",
    "physical_energy_shift_mev":0.0,"source_sha256_before":before,
    "binary_sha256":sha(a.build.resolve()/"compose_eos_fortran"),"pairs":[]}
 save(out/"RECEIPT.json",r)
 for crust in SPG_IDS:
  for core in CORE_IDS:
   dest=out/(crust+"__"+core)
   result=run_pipeline("maxwell",a.spg_root.resolve()/crust/"generated",
                       a.core_root.resolve()/core/"generated",dest)
   dest.mkdir(parents=True,exist_ok=True)
   save(dest/"PAIR.json",result)
   r["pairs"].append({"crust":crust,"core":core,"survived":result.get("survived",False),
       "stage":result.get("stage"),"failure_category":result.get("failure_category"),
       "native_category":result.get("native_category"),"error":result.get("error"),
       "receipt":str((dest/"PAIR.json").relative_to(out))})
   save(out/"RECEIPT.json",r)
   print(crust,core,result.get("survived"),result.get("native_category"),result.get("error"),flush=True)
 after={str(x.relative_to(ROOT)):sha(x)for x in deps}
 r.update(finished_unix=time.time(),source_sha256_after=after,sources_unchanged=before==after,
          admitted=sum(x["survived"]for x in r["pairs"]),
          categories=dict(collections.Counter(x["failure_category"]or"admitted"for x in r["pairs"])),
          execution_status="COMPLETE"if len(r["pairs"])==25 and before==after else"IDENTITY_CHANGED")
 save(out/"RECEIPT.json",r)
 return 0 if r["execution_status"]=="COMPLETE"else 1
if __name__=="__main__":raise SystemExit(main())
