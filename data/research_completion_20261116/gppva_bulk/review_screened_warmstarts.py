"""Adversarial stale electron and pure gauge initial guesses; no source edits."""
import json,hashlib,sys
from pathlib import Path
import numpy as np
from compose_eos.spatial_tf import RadialGrid,ELECTRON_MASS
from compose_eos.spatial_tf_equilibrium import screened_electrons
rows=[]
for d in [1,2,3]:
 g=RadialGrid(32,d);p=.003*(1+.8*np.cos(np.pi*g.x));R=35.
 base=screened_electrons(g,p,R)
 guesses=[('stale_empty_mu',(np.zeros(32),ELECTRON_MASS)),('gauge_shift',(base[1]+1234.,base[2]-1234.))]
 try:guesses.append(('stale_low_density',screened_electrons(g,p*1e-12,R)[1:3]))
 except Exception as exc:rows.append(dict(dimension=d,initial='low_density_cold_solve',error=repr(exc)))
 for name,initial in guesses:
  try:
   got=screened_electrons(g,p,R,initial)
   row=dict(dimension=d,initial=name,ne_difference=float(max(abs(got[0]-base[0]))),V_difference=float(max(abs(got[1]-base[1]))),mu_difference=abs(got[2]-base[2]),charge=abs(g.mean(p-got[0])),pde=got[3],gauge=abs(g.mean(got[1])))
   assert row['ne_difference']<1e-10 and row['V_difference']<1e-8 and row['mu_difference']<1e-8
  except Exception as exc:row=dict(dimension=d,initial=name,error=repr(exc))
  rows.append(row)
 p0=screened_electrons(g,np.zeros(32),R,base[1:3])
 assert np.all(p0[0]==0) and np.all(p0[1]==0) and p0[2]==ELECTRON_MASS
out={'scope':'Stale chemical potential, near-empty cached electrons, gauge-shift warm starts and exact zero-charge branch in all three geometries','source_sha256':hashlib.sha256(Path('python_frontend/compose_eos/spatial_tf_equilibrium.py').read_bytes()).hexdigest(),'rows':rows}
Path(sys.argv[1]).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps(rows,indent=2))
assert not any('error' in row for row in rows)
