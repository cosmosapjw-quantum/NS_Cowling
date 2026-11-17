from pathlib import Path
import sys,json,hashlib,subprocess,datetime
import numpy as np
from scipy.optimize import brentq
f=Path('/home/cosmosapjw/ns_cowling_wt/new_model_maxwell_20261117')
g=Path('/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117')
sys.path.insert(0,str(f/'python_frontend'))
from compose_eos.rmf_candidate import load_parameters,evaluate_states,sha256,HBARC
out=g/'data/crust_refinement_20261117/source_physics/nmp_definition_diagnostic';out.mkdir(parents=True,exist_ok=True)
exe=f/'data/new_model_maxwell_20261117/uniform_validation/build/rmf_candidate_driver'
source=f/'data/new_model_maxwell_20261117/sources/PARAMETER_INPUTS.json'
records=[]
for base in ['DDMEX_ORIGINAL_TABLE_RHO0','DDMEX_PROVIDER_ROUNDED_N0']:
 for mass in [938.5,939.0]:
  p=load_parameters(source,base)
  p.update(id=base+'_M'+str(mass),mn=mass,mp=mass)
  folder=out/p['id'];folder.mkdir(exist_ok=True)
  counter=[0];cache={}
  def state(n,yp=.5):
   key=(float(n),float(yp))
   if key not in cache:
    counter[0]+=1
    cache[key]=evaluate_states(exe,p,[n],mode='baryons',yp=yp,work_dir=folder/f'eval_{counter[0]:03d}')
   return cache[key]
  def pressure(n):return float(state(n)['pressure'][0])
  ns=brentq(pressure,.14,.16,xtol=2e-13)
  def epa(n,yp):return float(state(n,yp)['epsilon'][0])/n
  def quad(n,delta):
   return (epa(n,(1-delta)/2)+epa(n,(1+delta)/2)-2*epa(n,.5))/(2*delta*delta)
  def diff(n):return epa(n,0)-epa(n,.5)
  steps=[]
  for h in [1e-5,5e-6]:
   for delta in [.01,.005]:
    steps.append({'density_step':h,'asymmetry_step':delta,'J_quadratic':quad(ns,delta),
     'L_quadratic':3*ns*(quad(ns+h,delta)-quad(ns-h,delta))/(2*h),
     'J_PNM_minus_SNM':diff(ns),'L_PNM_minus_SNM':3*ns*(diff(ns+h)-diff(ns-h))/(2*h),
     'K':9*(pressure(ns+h)-pressure(ns-h))/(2*h)})
  rec={'parameter_variant':base,'diagnostic_equal_mass_MeV':mass,'n0_reference':p['n0'],'n_sat':ns,
   'binding_at_saturation':epa(ns,.5)-mass,'steps':steps,'state_calls':counter[0],
   'scope':'Independent diagnostic variants; NOT replacement of frozen F/G matching parameters; no parameter fitting',
   'source_parameter_sha256':sha256(source),'executable_sha256':sha256(exe)}
  records.append(rec)
  (folder/'SUMMARY.json').write_text(json.dumps(rec,indent=2)+'\n')
  print(json.dumps(rec))
(out/'NMP_DIAGNOSTIC.json').write_text(json.dumps({'created_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'interpretation':'Controlled mass and symmetry-energy definition test; original author definition requires primary documentation.','records':records},indent=2)+'\n')
