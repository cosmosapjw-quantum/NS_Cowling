"""Independently recompute four screened-electron/KKT repaired TF profiles."""
import gzip,hashlib,json,sys
from pathlib import Path
import numpy as np
from compose_eos.spatial_tf import CellFunctional,RadialGrid
from compose_eos.spatial_tf_equilibrium import screened_electrons

paths=['data/research_completion_20261116/gppva_bulk/repair_inputs/reduced_GPPVA_FSU2_0.03.json.gz',
       'data/research_completion_20261116/gppva_bulk/repair_inputs/reduced_GPPVA_DD2_0.03.json.gz',
       'data/research_completion_20261116/gppva_bulk/repair_inputs/reduced_GPPVA_TW_0.07.json.gz',
       'data/research_completion_20261116/gppva_bulk/repair_inputs/nucleon_polished_1e-13.json.gz']
rows=[]
for path in paths:
 p=Path(path);old=json.loads(gzip.decompress(p.read_bytes()));g=RadialGrid(old['points'],old['dimension'])
 n=old['n_b'];R=old['radius_fm'];nn=np.asarray(old['profile']['nn']);np_=np.asarray(old['profile']['np'])
 ne,v,muq,eres=screened_electrons(g,np_,R)
 f=CellFunctional(old['model_id'],g,n);s=f.evaluate(nn,np_,ne,R)
 # Use the neutron reservoir alone for a genuinely different multiplier estimator.
 active_n=nn>max(1e-10,n*1e-7);active_p=np_>max(1e-10,n*1e-7)
 mub=float(g.weights[active_n]@s['mu_n'][active_n]/sum(g.weights[active_n]))
 residuals=[];species=[]
 for label,mask,density,chemical in [('n',active_n,nn,s['mu_n']),('p',active_p,np_,s['mu_p']+muq)]:
  kkt=np.where(mask,abs(chemical-mub),np.maximum(mub-chemical,0.))
  i=int(np.argmax(kkt));residuals.extend(kkt)
  species.append(dict(species=label,max_KKT_MeV=float(kkt[i]),worst_cell=i,worst_cell_density=float(density[i]),density_weighted_stationarity=float(g.mean(density*(chemical-mub))/n)))
 def energy(scale):
  a,b=scale*nn,scale*np_
  e,_,_,_=screened_electrons(g,b,R)
  return CellFunctional(old['model_id'],g,n*scale).evaluate(a,b,e,R)['epsilon']
 h=1e-5
 p_energy=(energy(1+h)-energy(1-h))/(2*h)-s['epsilon']
 p_euler=n*mub-s['epsilon']
 row=dict(path=path,input_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),model_id=old['model_id'],n_b=n,points=old['points'],dimension=old['dimension'],energy_recompute_difference=s['epsilon']-old['epsilon'],electron_profile_difference=float(max(abs(ne-np.asarray(old['profile']['ne'])))),electron_chemical_residual=float(max(abs(s['mu_e']-muq))),neutrality=abs(g.mean(np_-ne)),baryon_error=abs(g.mean(nn+np_)-n),independent_mu_B=mub,stored_mu_B=old['mu_B'],independent_pressure=p_euler,energy_derivative_pressure=p_energy,pressure_derivative_difference=p_energy-p_euler,species=species,max_KKT_MeV=float(max(residuals)),radius_stationarity_MeV=abs(s['d_epsilon_d_log_radius']/n),old_status=old['numerical_status'])
 rows.append(row)
 assert row['max_KKT_MeV']<.02
 assert row['electron_chemical_residual']<1e-7
 assert abs(row['energy_recompute_difference'])<1e-8
 assert abs(row['pressure_derivative_difference'])<2e-5
sources=['python_frontend/compose_eos/spatial_tf.py','python_frontend/compose_eos/spatial_tf_equilibrium.py','python_frontend/compose_eos/gppva_rmf.py']
result=dict(scope='Four saved red/green cases independently reevaluated; neutron-only chemical multiplier, full constrained density derivative; finite local stationarity, not global ground-state proof',source_sha256={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in sources},rows=rows)
Path(sys.argv[1]).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps([{k:r[k] for k in ['model_id','n_b','max_KKT_MeV','electron_chemical_residual','pressure_derivative_difference','radius_stationarity_MeV']} for r in rows],indent=2))
