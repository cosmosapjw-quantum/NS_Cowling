"""Interpolate saved physical trial profiles to radial grids for numerical checks."""
from pathlib import Path
import sys,json,importlib.util
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'python_frontend'))
OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('compose_eos._electron_flux_meshes',OUT/'flux_spatial_tf_equilibrium.py')
eq=importlib.util.module_from_spec(spec);spec.loader.exec_module(eq)
results=[]
for p in sorted(OUT.glob('*slab_failure_input.json')):
 c=json.loads(p.read_text());oldgrid=eq.RadialGrid(c['points'],c['dimension'])
 for n in (192,384,768):
  grid=eq.RadialGrid(n,c['dimension']);protons=np.interp(grid.x,oldgrid.x,c['protons'])
  initial=None if c['initial_v'] is None else (np.interp(grid.x,oldgrid.x,c['initial_v']),c['initial_mu'])
  try:
   ne,v,mu,res=eq.screened_electrons(grid,protons,c['radius'],initial=initial)
   r=dict(case=p.name,points=n,status='PASS',pde=res,charge=grid.mean(ne-protons),mu=mu)
  except Exception as e:r=dict(case=p.name,points=n,status='FAIL',error=str(e))
  results.append(r);print(r,flush=True)
(OUT/'FLUX_MESH_CHECKS.json').write_text(json.dumps(results,indent=2))
