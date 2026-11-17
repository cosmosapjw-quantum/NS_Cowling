#!/usr/bin/env python3
"""Inspect all retained SPG drip critical directions, including inward gas."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.candidate_gdfm_cld import reduced_cell,SurfaceParameters
base=ROOT/'data/research_finish_20261116/spg/grid257';checks=[]
for path in sorted(base.glob('*/*/RECEIPT.json')):
 r=json.loads(path.read_text());model=r['model'];p=SurfaceParameters(**r['surface_parameters'])
 d=r['drip'];n=d['n_B'];z=np.array([1.,d['coordinates'][2],d['coordinates'][3]])
 g0=reduced_cell(model,n,z,surface=p)[1];spectra=[]
 for ht in [2e-5,2e-6]:
  columns=[]
  for j,h in enumerate([ht,2e-7,2e-7]):
   a=z.copy();b=z.copy();b[j]-=h
   if j:a[j]+=h
   columns.append((reduced_cell(model,n,a,surface=p)[1]-reduced_cell(model,n,b,surface=p)[1])/(a[j]-b[j]))
  H=np.array(columns).T;scale=np.array([1.,.1,.1])
  spectra.append(float(np.linalg.eigvalsh((H+H.T)*.5*scale[:,None]*scale[None,:])[0]))
 checks.append(dict(model=model,label=r['label'],n_B=n,boundary_gradient=g0.tolist(),
  min_scaled_full_hessian_at_inward_steps=spectra))
assert len(checks)==25 and min(min(r['min_scaled_full_hessian_at_inward_steps']) for r in checks)>0
result=dict(scope='25 finite drip states; one-sided gas insertion curvature included. The neutron-gas Fermi derivative is nonanalytic exactly at zero density, so these finite-step Hessians are diagnostics, not a uniform endpoint theorem.',
 cases=checks,min_scaled_eigenvalue=min(min(r['min_scaled_full_hessian_at_inward_steps']) for r in checks),
 runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('SPG_DRIP_CURVATURE_ORACLE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))
