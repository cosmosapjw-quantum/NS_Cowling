"""Read-only reevaluation of every selected spatial TF campaign profile."""
import gzip,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from compose_eos.spatial_tf import CellFunctional,RadialGrid
from compose_eos.spatial_tf_equilibrium import screened_electrons
start=time.monotonic();folder=Path(sys.argv[1]);rows=[];files={}
source_paths=['python_frontend/compose_eos/spatial_tf.py','python_frontend/compose_eos/spatial_tf_equilibrium.py','python_frontend/compose_eos/gppva_rmf.py']
source_hashes={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in source_paths}
for path in sorted(folder.glob('*/selected_states.json.gz')):
 files[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
 states=json.loads(gzip.decompress(path.read_bytes()))
 for old in states:
  row=dict(model_id=old['model_id'],n_b=old['n_b'],dimension=old['dimension'],points=old['points'],status=old['numerical_status'])
  try:
   g=RadialGrid(old['points'],old['dimension']);n=old['n_b'];R=old['radius_fm']
   nn=np.asarray(old['profile']['nn']);p=np.asarray(old['profile']['np'])
   e,v,muq,eres=screened_electrons(g,p,R)
   s=CellFunctional(old['model_id'],g,n).evaluate(nn,p,e,R)
   an=nn>max(1e-10,n*1e-7);ap=p>max(1e-10,n*1e-7)
   # Independent neutron-reservoir estimate rather than the optimizer's joint fit.
   mub=float(g.weights[an]@s['mu_n'][an]/sum(g.weights[an]))
   kn=np.where(an,abs(s['mu_n']-mub),np.maximum(mub-s['mu_n'],0.))
   kp=np.where(ap,abs(s['mu_p']+muq-mub),np.maximum(mub-s['mu_p']-muq,0.))
   def energy(scale):
    a,b=nn*scale,p*scale;z,_,_,_=screened_electrons(g,b,R)
    return CellFunctional(old['model_id'],g,n*scale).evaluate(a,b,z,R)['epsilon']
   h=1e-5;p_density=(energy(1+h)-energy(1-h))/(2*h)-s['epsilon'];p_euler=n*mub-s['epsilon']
   row.update(energy_difference=s['epsilon']-old['epsilon'],max_KKT_MeV=float(max(max(kn),max(kp))),electron_chemical_error=float(max(abs(s['mu_e']-muq))),electron_pde=eres,charge=abs(g.mean(p-e)),baryon_error=abs(g.mean(nn+p)-n),mu_B_difference=mub-old['mu_B'],pressure_derivative_difference=p_density-p_euler,radius_stationarity_MeV=abs(s['d_epsilon_d_log_radius']/n),max_local_baryon_density=float(max(nn+p)))
  except Exception as exc:row['error']=repr(exc)
  rows.append(row)
summary=dict(models=len(files),selected_profiles=len(rows),runtime_errors=sum('error'in r for r in rows))
for key in ['energy_difference','max_KKT_MeV','electron_chemical_error','electron_pde','charge','baryon_error','mu_B_difference','pressure_derivative_difference','radius_stationarity_MeV']:
 summary['max_abs_'+key]=max(abs(r.get(key,0.)) for r in rows)
result=dict(scope='Every selected finite-grid profile independently reevaluated, including density-derived pressure; no proof of global minimum or unsampled convergence',source_sha256=source_hashes,input_sha256=files,summary=summary,rows=rows,runtime_seconds=time.monotonic()-start)
Path(sys.argv[2]).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(summary,indent=2))
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in source_hashes.items())
assert summary['models']>=int(sys.argv[3] if len(sys.argv)>3 else 8) and summary['runtime_errors']==0
assert summary['max_abs_energy_difference']<1e-8
assert summary['max_abs_max_KKT_MeV']<.02
assert summary['max_abs_electron_chemical_error']<1e-7
assert summary['max_abs_charge']<2e-11
assert summary['max_abs_pressure_derivative_difference']<2e-5
