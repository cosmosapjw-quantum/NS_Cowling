#!/usr/bin/env python3
"""Replay dense captured inputs against final production and independent MP residuals."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.spatial_tf import RadialGrid,E2
from compose_eos.spatial_tf_equilibrium import screened_electrons
mp.mp.dps=70
module=ROOT/'python_frontend/compose_eos/spatial_tf_equilibrium.py'
initialhash=hashlib.sha256(module.read_bytes()).hexdigest()
cases=[]
for path in sorted((ROOT/'data/research_finish_20261116/electron_repair').glob('GPPVA_*_slab_failure_input.json')):
 r=json.loads(path.read_text());N=r['points'];d=r['dimension'];R=mp.mpf(r['radius'])
 grid=RadialGrid(N,d)
 initial=None if r['initial_v'] is None else (np.asarray(r['initial_v']),r['initial_mu'])
 ne,v,mu,error=screened_electrons(grid,r['protons'],r['radius'],initial)
 wm=[(mp.mpf(i+1)/N)**d-(mp.mpf(i)/N)**d for i in range(N)]
 c=[d*(mp.mpf(i+1)/N)**(d-1)*N/R**2 for i in range(N-1)]
 vm=list(map(mp.mpf,v));nm=list(map(mp.mpf,ne));pm=list(map(mp.mpf,r['protons']))
 F=[c[i]*(vm[i]-vm[i+1]) for i in range(N-1)]
 res=[((F[i] if i<N-1 else 0)-(F[i-1] if i else 0))/wm[i]-4*mp.pi*mp.mpf(E2)*(pm[i]-nm[i]) for i in range(N)]
 rec=dict(case=path.name,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),production_PDE=error,
  independent_70digit_PDE=float(max(abs(q) for q in res)),
  relative_neutrality=float(abs(mp.fsum(wm[i]*(nm[i]-pm[i]) for i in range(N)))/mp.fsum(wm[i]*pm[i] for i in range(N))),
  mean_potential_MeV=grid.mean(v),output_dtype=str(v.dtype))
 assert rec['independent_70digit_PDE']<2e-10
 assert rec['relative_neutrality']<1e-10
 cases.append(rec)
result=dict(source_sha256=initialhash,source_unchanged=initialhash==hashlib.sha256(module.read_bytes()).hexdigest(),
 cases=cases,unique_input_count=len({r['input_sha256'] for r in cases}),
 max_independent_PDE=max(r['independent_70digit_PDE'] for r in cases),
 max_relative_neutrality=max(r['relative_neutrality'] for r in cases),
 runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('FINAL_DENSE_ELECTRON_ORACLE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2));assert result['source_unchanged']
