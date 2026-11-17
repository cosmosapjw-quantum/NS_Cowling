#!/usr/bin/env python3
"""Run bounded native SPG campaign with a versioned conditional surface input."""
import concurrent.futures, hashlib, json, os, subprocess, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.generation_dispatch import generation_request
from compose_eos.native_generation import surface_input,write_request,verify_native_generation
from compose_eos.native_model_contracts import SPG_INFERRED_VARIANT
BUILD=Path(sys.argv[1]).resolve();OUT=Path(sys.argv[2]).resolve()
OUT.mkdir(parents=True,exist_ok=True)
SOURCE=Path(sys.argv[3]).resolve() if len(sys.argv)>3 else ROOT
ENV={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','COMPOSE_EOS_FORTRAN_EXE':str(BUILD/'compose_eos_fortran')}
os.environ.update({k:ENV[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','COMPOSE_EOS_FORTRAN_EXE']})
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
deps=[SOURCE/'fortran_eos/src'/n for n in ['compose_eos_types.f90','compose_eos_optimize_native.f90','compose_eos_gdfm_native.f90','compose_eos_spg_native.f90']]
deps += [SOURCE/'fortran_eos/app/compose_eos_generate.f90',BUILD/'compose_eos_generate',BUILD/'compose_eos_fortran',ROOT/'fortran_eos/tests/spg_finish_run.py']
manifest={'baseline':'478c7b5de1514ac512f91156f36551478f36f688','compiled_source_root':str(SOURCE),'scope':'conditional inferred surface; 129 native stationary rows per model; no original-calibration or full-crust equivalence claim','workers':2,'threads':1,'files':{str(p):sha(p) for p in deps},'started_unix':time.time()}
(OUT/'INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
def one(i):
 model=f'SPG_M{i}';out=OUT/model
 out.mkdir(exist_ok=False)
 request=generation_request(model,variant=SPG_INFERRED_VARIANT,points=129)
 nml=write_request(ROOT,out/'generated',request)
 command=[str(BUILD/'compose_eos_generate'),str(nml),str(out/'generated')]
 receipt={'request':request,'surface_input':surface_input(ROOT,model,SPG_INFERRED_VARIANT),'command':command,'request_sha256':sha(nml),'started_unix':time.time()}
 with (out/'execution.log').open('w') as log:
  run=subprocess.run(command,cwd=ROOT,env=ENV,stdout=log,stderr=subprocess.STDOUT)
 receipt.update(returncode=run.returncode,finished_unix=time.time())
 if run.returncode==0:
  try:receipt['evidence']=verify_native_generation(out/'generated',request);receipt['status']='REQUESTED_ROWS_NATIVE_VERIFIED'
  except Exception as exc:receipt.update(status='EVIDENCE_FAILED',error=str(exc))
 else:receipt['status']='NATIVE_EXECUTION_FAILED'
 receipt['output_sha256']={str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()}
 (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
 print(model,receipt['status'],round(receipt['finished_unix']-receipt['started_unix'],2),flush=True)
 return model,receipt['status']
with concurrent.futures.ThreadPoolExecutor(max_workers=2)as pool: results=dict(pool.map(one,range(1,6)))
manifest.update(finished_unix=time.time(),models=results,physical_dependencies_unchanged=all(sha(p)==manifest['files'][str(p)]for p in deps[:4]))
(OUT/'RUN_RECEIPT.json').write_text(json.dumps(manifest,indent=2)+'\n')
raise SystemExit(0 if all(s=='REQUESTED_ROWS_NATIVE_VERIFIED'for s in results.values())else 1)
