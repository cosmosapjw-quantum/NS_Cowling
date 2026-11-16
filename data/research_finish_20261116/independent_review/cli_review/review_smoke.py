#!/usr/bin/env python3
"""Independently read native public-CLI smoke products and check their receipts."""
from pathlib import Path
import sys,json,gzip,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos import generation_dispatch as dispatch
from compose_eos.api import read_compose
BASE=ROOT/'data/research_finish_20261116/cli_smoke';OUT=Path(__file__).parent
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
summary=read(BASE/'SUMMARY.json');results=[]
catalog=dispatch.generator_catalog()
assert set(catalog)==set(dispatch.MODEL_IDS) and len(catalog)==18
assert len(dispatch.CORE_IDS)==5 and len(dispatch.SPG_IDS)==5 and len(dispatch.TF_IDS)==8
for case in ['core3','apr3','spg17']:
 p=BASE/case;receipt=read(p/'DISPATCH.json');request=receipt['request']
 evidence=dispatch.verify_generation(p/'generated',request)
 native=read_compose(Path(evidence['native_directory']))
 assert native.n==summary['cases'][case]['native_rows']
 for key in ['Yn_present','Yp_present','Ye_present','Ymu_present','C3_applicable','C4_applicable']:assert getattr(native,key).all()
 assert np.all(np.diff(native.nb)>0)
 sources=summary['cases'][case]
 assert sources['initial_source_sha256']==sources['after_source_sha256']
 source_matches={name:sha(ROOT/name)==value for name,value in sources['after_source_sha256'].items()}
 assert all(source_matches.values())
 errors={}
 if case!='spg17':
  record=read(p/'generated'/request['model']/'RECEIPT.json');rows=record['rows']
  mapping={'eps':'epsilon','p':'pressure','mub':'mu_n','Yn':'Y_n','Yp':'Y_p','Ye':'Y_e','Ymu':'Y_mu'}
 else:
  directory=p/'generated'/request['model']/request['variant']
  rows=json.load(gzip.open(directory/'all_attempts.json.gz','rt'))['fine']
  rows=[{**r['best'],**r['best']['composition']} for r in rows]
  mapping={'eps':'energy','p':'pressure_MeV_fm3','mub':'mu_B_MeV','Yn':'Yn','Yp':'Yp','Ye':'Ye','Ymu':'Ymu'}
 for field,key in mapping.items():
  errors[field]=float(max(abs(getattr(native,field)-np.asarray([r[key] for r in rows]))))
 assert max(errors.values())<1e-10
 results.append(dict(case=case,model=request['model'],requested_rows=request['points'],native_rows=native.n,source_matches=source_matches,errors=errors))
result=dict(scope='18-entry inventory and actual retained GM1/APR/SPG native smoke; no TF CLI solve or whole-model calibration certificate',
 catalog={m:{'family':c['family'],'region':c['region']} for m,c in catalog.items()},cases=results,
 native_executable_sha256=sha(ROOT/'build/compose_eos_fortran/compose_eos_fortran'),
 dispatcher_sha256=sha(Path(dispatch.__file__)),runner_sha256=sha(Path(__file__)))
assert result['native_executable_sha256']==summary['native_executable_sha256']
(OUT/'CLI_SMOKE_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
