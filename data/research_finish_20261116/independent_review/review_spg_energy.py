#!/usr/bin/env python3
"""Independent reconstruction of the new explicit-surface energy and derivatives."""
import sys,json,hashlib,math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from compose_eos import candidate_gdfm_cld as cld
from compose_eos.gdfm_uniform import baryons,MASS
from compose_eos.cold_leptons import free_fermion,ELECTRON_MASS
from reoptimize_spg_surfaces import prescriptions
OUT=Path(__file__).parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
sources={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'python_frontend/compose_eos/candidate_gdfm_cld.py',ROOT/'python_frontend/compose_eos/spg_crust_generator.py']}
def energy(model,n,z,p):
    # Direct Eq. 21-25 reconstruction; do not call production CLD, phi, sigma or gradient.
    t,ni,y,v=map(float,z)
    f=t*n/ni;ng=n*(1-t)/(1-f);yg=v*y;ne=f*ni*y+(1-f)*ng*yg
    phi=(2-3*f**(1/3)+f)/5
    sigma=p.sigma0_MeV_fm2*(16+p.b)/(y**(-3)+p.b+(1-y)**(-3))
    contrast=ni*y-ng*yg
    radius=(3*sigma/(4*math.pi*1.43996448*phi*contrast**2))**(1/3)
    bulk=f*baryons(model,ni,y)['epsilon_baryons']+(1-f)*baryons(model,ng,yg)['epsilon_baryons']
    eps=bulk+free_fermion(ne,ELECTRON_MASS)[0]+3*f*sigma/radius+2*math.pi*1.43996448*f*phi*radius**2*contrast**2
    return eps/n-MASS
checks=[]
for model in [f'SPG_M{i}' for i in range(1,6)]:
 for label,p in prescriptions(model).items():
  for n,z in [(0.03,np.array([.3,.14,.1,.02])),(.008,np.array([.7,.15,.25,.005]))]:
   e,g,d=cld.reduced_two_phase_cell(model,n,z,surface=p)
   independent=energy(model,n,z,p)
   numeric=[];changes=[]
   for j,h in enumerate([1e-4,2e-5,2e-5,1e-4]):
    estimates=[]
    for hh in [h,h/2]:
     zs=[]
     for factor in [-2,-1,1,2]:
      zz=z.copy();zz[j]+=factor*hh;zs.append(energy(model,n,zz,p))
     estimates.append((zs[0]-8*zs[1]+8*zs[2]-zs[3])/(12*hh))
    numeric.append(estimates[1]);changes.append(abs(estimates[1]-estimates[0]))
   checks.append({'model':model,'label':label,'n_B':n,'z':z.tolist(),'absolute_energy_error_MeV_per_baryon':abs(e-independent),
     'analytic_gradient':g.tolist(),'independent_gradient':numeric,'absolute_gradient_error':abs(g-numeric).tolist(),
     'half_step_derivative_change':changes,'surface_parameters':d['surface_parameters']})
result={'scope':'50 interior coordinate cases across all25 explicit surfaces; independent direct energy reconstruction and fourth-order central finite differences, not a calibration or global-minimum proof',
 'units':'n_B fm^-3, energy MeV per baryon, coordinate derivatives in corresponding reduced coordinates',
 'source_sha256':sources,'source_unchanged':all(sha(ROOT/p)==v for p,v in sources.items()),
 'max_energy_absolute_error':max(c['absolute_energy_error_MeV_per_baryon'] for c in checks),
 'max_gradient_absolute_error':max(max(c['absolute_gradient_error']) for c in checks),
 'max_derivative_halfstep_change':max(max(c['half_step_derivative_change']) for c in checks),
 'cases':checks,'runner_sha256':sha(Path(__file__))}
(OUT/'SPG_ENERGY_ORACLE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))
assert result['source_unchanged']
assert result['max_energy_absolute_error']<5e-12
assert result['max_gradient_absolute_error']<2e-6
