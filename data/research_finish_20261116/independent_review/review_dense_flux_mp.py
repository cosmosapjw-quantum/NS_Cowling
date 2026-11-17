#!/usr/bin/env python3
"""Independent high-precision FV action on the six captured dense failures."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.spatial_tf import RadialGrid,E2
mp.mp.dps=70
OUT=Path(__file__).parent
cases=[]
for path in sorted((ROOT/'data/research_finish_20261116/electron_repair').glob('GPPVA_*_slab_failure_input.json')):
 r=json.loads(path.read_text());N=r['points'];d=r['dimension'];R=mp.mpf(r['radius'])
 v=list(map(mp.mpf,r['v']));np_=list(map(mp.mpf,r['protons']));ne=list(map(mp.mpf,r['ne']))
 # Exact geometric FV operator, rather than supplied floating dense matrix.
 w=[(mp.mpf(i+1)/N)**d-(mp.mpf(i)/N)**d for i in range(N)]
 c=[d*(mp.mpf(i+1)/N)**(d-1)*N/R**2 for i in range(N-1)]
 F=[c[i]*(v[i]-v[i+1]) for i in range(N-1)]
 lv=[((F[i] if i<N-1 else 0)-(F[i-1] if i else 0))/w[i] for i in range(N)]
 reference=np.array(list(map(float,lv)));vf=np.asarray(r['v'])
 grid=RadialGrid(N,d);L=grid.minus_laplacian/r['radius']**2
 dense=L@vf
 flux=-np.diag(L,1)*grid.weights[:-1]*(vf[:-1]-vf[1:])
 stable=(np.r_[flux,0.]-np.r_[0.,flux])/grid.weights
 charge=np.asarray(r['protons'])-np.asarray(r['ne'])
 rec=dict(case=path.name,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
  dense_max_operator_error=float(max(abs(dense-reference))),
  flux_max_operator_error=float(max(abs(stable-reference))),
  stored_iteration_error=r['error'],dense_PDE_max=float(max(abs(dense-4*np.pi*E2*charge))),
  flux_PDE_max=float(max(abs(stable-4*np.pi*E2*charge))),
  exact_geometry_70digit_PDE_max=float(max(abs(lv[i]-4*mp.pi*mp.mpf(E2)*(np_[i]-ne[i])) for i in range(N))),
  flux_weighted_conservation=float(abs(grid.weights@stable)))
 cases.append(rec)
assert len(cases)==6
result=dict(method='70-digit exact rational radial geometry FV operator on retained potential and densities, independently compared with float dense and face-flux actions',cases=cases,
 max_dense_operator_error=max(r['dense_max_operator_error'] for r in cases),
 max_flux_operator_error=max(r['flux_max_operator_error'] for r in cases),
 runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(OUT/'DENSE_FLUX_ORACLE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
