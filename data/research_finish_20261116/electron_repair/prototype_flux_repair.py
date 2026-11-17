"""Evaluate actual slab failures with conservative flux differences only."""
from pathlib import Path
import sys,json,inspect,importlib.util
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'python_frontend'))
OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('compose_eos._electron_flux_prototype',OUT/'neutral_pinned_spatial_tf_equilibrium.py')
eq=importlib.util.module_from_spec(spec);spec.loader.exec_module(eq)
s=inspect.getsource(eq.screened_electrons)
s=s.replace('return ne,d,lap@v-coupling*(np_-ne)', '''flux=-np.diag(lap,1)*weights[:-1]*np.diff(-v)
        action=np.r_[flux,0.]-np.r_[0.,flux]
        return ne,d,action/weights-coupling*(np_-ne)''')
exec(s,eq.__dict__)
results=[]
for p in sorted(OUT.glob('*failure_input.json')):
 c=json.loads(p.read_text());grid=eq.RadialGrid(c['points'],c['dimension'])
 try:
  ne,v,mu,res=eq.screened_electrons(grid,np.asarray(c['protons']),c['radius'],initial=(np.asarray(c['initial_v']),c['initial_mu']) if c['initial_v'] is not None else None)
  r=dict(case=p.name,status='PASS',pde=res,charge=grid.mean(ne-np.asarray(c['protons'])),mu=mu)
 except Exception as e:r=dict(case=p.name,status='FAIL',error=str(e))
 results.append(r);print(r,flush=True)
(OUT/'FLUX_PROTOTYPE.json').write_text(json.dumps(results,indent=2))
