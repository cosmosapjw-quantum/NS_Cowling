"""Independent optimized-branch thermodynamic and mesh check for DDME2."""
import hashlib,json,pathlib
import numpy as np
from compose_eos.spatial_tf import solve_cell,HBARC
HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[2]
n=.03
coarse=solve_cell('DDME2',n,points=32,dimension=3,maxiter=1000,ftol=1e-12)
fine=solve_cell('DDME2',n,points=64,dimension=3,maxiter=1400,ftol=1e-13)
initial=dict(nn=fine['profile']['nn'],np=fine['profile']['np'],ne=fine['profile']['ne'],radius=fine['radius_fm'])
near=[]
for h in (3e-5,1e-5):
 rows=[]
 for sign in (-1,1):
  ni=n+sign*h
  init={k:(np.asarray(v)*ni/n).tolist() for k,v in initial.items() if k!='radius'}
  init['radius']=initial['radius']
  rows.append(solve_cell('DDME2',ni,points=64,dimension=3,maxiter=1400,ftol=1e-13,initial=init))
 pfd=n*n*(rows[1]['energy_per_baryon']-rows[0]['energy_per_baryon'])/(2*h)
 near.append(dict(h=h,pressure_fd=pfd,pressure_euler=fine['pressure'],relative_difference=(pfd-fine['pressure'])/fine['pressure'],states=rows))
source=dict(paper='Grill2014 arXiv1404.2753 Table4',n_b=n,epsilon_MeV_fm3=1.43578*HBARC/10,pressure_MeV_fm3=.40065*HBARC/1000,policy='rounded published table diagnostic, no exactmatching gate')
summary=dict(coarse_status=coarse['numerical_status'],fine_status=fine['numerical_status'],
 coarse_pressure=coarse['pressure'],fine_pressure=fine['pressure'],
 mesh_pressure_relative_change=(fine['pressure']-coarse['pressure'])/fine['pressure'],
 coarse_energy_per_baryon=coarse['energy_per_baryon'],fine_energy_per_baryon=fine['energy_per_baryon'],
 mesh_energy_change_MeV=fine['energy_per_baryon']-coarse['energy_per_baryon'],
 fine_kkt_MeV=fine['max_chemical_KKT_residual_MeV'],
 fine_radius_stationarity_MeV=fine['radius_stationarity_MeV'],
 source_pressure_relative_difference=(fine['pressure']-source['pressure_MeV_fm3'])/source['pressure_MeV_fm3'],
 source_energy_relative_difference=(fine['epsilon']-source['epsilon_MeV_fm3'])/source['epsilon_MeV_fm3'],
 pressure_derivative_relative_differences=[r['relative_difference'] for r in near])
result=dict(scope='One independently optimized spherical DDME2 branch; no geometry-global or full-density claim',summary=summary,coarse=coarse,fine=fine,near=near,source_diagnostic=source,
 code_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'python_frontend/compose_eos/spatial_tf.py',ROOT/'python_frontend/compose_eos/gppva_rmf.py']})
(HERE/'optimized_branch_review.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(summary,indent=2))
