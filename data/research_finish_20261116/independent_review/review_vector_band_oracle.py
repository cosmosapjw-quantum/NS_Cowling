#!/usr/bin/env python3
"""Read-only oracle for a future coupled-vector band solve; production untouched."""
from pathlib import Path
import sys,gzip,json,time,hashlib
import numpy as np
from scipy.linalg import solve_banded
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.spatial_tf import RadialGrid,HBARC,HC3
from compose_eos import gppva_rmf as rmf
cases=[]
for model in ['GPPVA_FSU2R','GPPVA_FSU2H','GPPVA_TM1e']:
 p=ROOT/f'data/research_completion_20261116/gppva_tf/final/{model}/selected_states.json.gz'
 old=json.loads(gzip.decompress(p.read_bytes()));state=min(old,key=lambda r:abs(r['n_b']-.05))
 for points in [192,384]:
  g=RadialGrid(points,state['dimension']);profile=state['profile']
  nn,np_,w,r=[np.interp(g.x,profile['x'],profile[k]) for k in ['nn','np','omega','rho']]
  lap=HBARC**2*g.minus_laplacian/state['radius_fm']**2
  v=rmf.vector_potential(model,w,r);c=rmf.couplings(model,nn+np_)
  jw=lap+np.diag(np.broadcast_to(v['dww'],nn.shape));jr=lap+np.diag(np.broadcast_to(v['drr'],nn.shape))
  cross=np.broadcast_to(v['dwr'],nn.shape);A=np.block([[jw,np.diag(cross)],[np.diag(cross),jr]])
  rhs=-np.r_[lap@w+v['dw']-c['gw']*(nn+np_)*HC3,lap@r+v['dr']-c['gr']*(np_-nn)*HC3/2]
  # Interleaved (omega0,rho0,omega1,rho1,...) ordering has |row-column|<=2.
  band=np.zeros((5,2*points))
  band[2,0::2]=np.diag(jw);band[2,1::2]=np.diag(jr)
  band[1,1::2]=cross;band[3,0::2]=cross
  band[0,2::2]=np.diag(lap,1);band[0,3::2]=np.diag(lap,1)
  band[4,:-2:2]=np.diag(lap,-1);band[4,1:-2:2]=np.diag(lap,-1)
  order=np.array([j for i in range(points) for j in [i,i+points]])
  timed={'dense':[],'band':[]}
  for _ in range(5):
   start=time.perf_counter();dense=np.linalg.solve(A,rhs);timed['dense'].append(time.perf_counter()-start)
   start=time.perf_counter();interleaved=solve_banded((2,2),band,rhs[order]);timed['band'].append(time.perf_counter()-start)
  alternative=np.empty_like(dense);alternative[order]=interleaved
  difference=float(max(abs(dense-alternative)))
  residual=float(max(abs(A@alternative-rhs)))
  assert np.allclose(dense,alternative,rtol=1e-10,atol=1e-13)
  cases.append(dict(model=model,points=points,source_profile_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
    max_solution_difference=difference,max_linear_residual=residual,
    dense_median_seconds=float(np.median(timed['dense'])),band_median_seconds=float(np.median(timed['band'])),
    median_speed_ratio=float(np.median(timed['dense'])/np.median(timed['band']))))
result=dict(scope='Read-only algebraic oracle on six actual-profile-derived Jacobians. It is not integrated, no complete cell rerun, no scientific output replacement.',
 cases=cases,runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('VECTOR_BAND_ORACLE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
