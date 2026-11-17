#!/usr/bin/env python3
"""Recover native readback after unrelated concurrent source-set changes, preserving failures."""
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.native_generation import verify_native_generation,source_files
BASE=ROOT/'data/research_blockers_20261116/tf/ambiguous_N128'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 result={'status':'READBACK_RECOVERY','original_failed_dispatches_preserved':True,'relation':'BYTE_IDENTITY for linked physical source files; added unlinked source files do not invalidate native output','models':{}}
 needed=['fortran_eos/src/'+x+'.f90' for x in ('compose_eos_types','compose_eos_io','compose_eos_optimize_native','compose_eos_core_native','compose_eos_gdfm_native','compose_eos_lower_crust_native','compose_eos_tf_models','compose_eos_spg_native','compose_eos_spatial_tf')]
 needed+=['fortran_eos/app/compose_eos_generate.f90','fortran_eos/vendor/slsqp/slsqp_optmz.f','fortran_eos/build.sh','fortran_eos/run-generator.sh','python_frontend/compose_eos/generation_dispatch.py','python_frontend/compose_eos/native_generation.py']
 for model in ('FSU2R','TM1e'):
  folder=BASE/model/'output';dispatch=folder/'DISPATCH.json';old=json.loads(dispatch.read_text());recorded=old['source_sha256']
  checks={p:{'recorded':recorded.get(str(ROOT/p)),'current':sha(ROOT/p)} for p in needed}
  mismatch=[p for p,v in checks.items() if v['recorded']!=v['current']]
  current={str(p):sha(p) for p in source_files(ROOT)}
  changes={'added_source_paths':[p for p in current if p not in recorded],'changed_source_paths':[p for p,h in current.items() if p in recorded and recorded[p]!=h],'removed_source_paths':[p for p in recorded if p.startswith(str(ROOT/'fortran_eos')) and p not in current]}
  if mismatch or old['child_returncode']!=0:raise RuntimeError((model,mismatch,old['child_returncode']))
  readback=verify_native_generation(folder/'generated',old['request'])
  result['models'][model]={'status':'NATIVE_RESULT_VALID_READBACK_RECOVERED','original_dispatch_status':old['status'],'original_dispatch_sha256':sha(dispatch),'native_child_returncode':old['child_returncode'],'linked_source_identity':checks,'mismatches':mismatch,'source_set_changes':changes,'native_readback':readback,'binary_sha256':sha(ROOT/'build/blocker_tf_native/compose_eos_generate')}
 (BASE/'READBACK_RECOVERY.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps({m:v['status'] for m,v in result['models'].items()}))
if __name__=='__main__':main()
