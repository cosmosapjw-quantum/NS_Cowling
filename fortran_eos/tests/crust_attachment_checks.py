#!/usr/bin/env python3
"""Bounded native tests. Python source pressure is an independent implementation oracle."""
import hashlib,json,os,pathlib,subprocess,sys,time
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
import numpy as np
from scipy.integrate import quad
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.bsk22_analytic import state as source_state
OUT=ROOT/'data/research_blockers_20261116/lower'
OUT.mkdir(parents=True,exist_ok=True)
BUILD=OUT/'build';BUILD.mkdir(exist_ok=True)
EXE=BUILD/'attachment'
SOURCES=[ROOT/p for p in ['fortran_eos/src/compose_eos_types.f90','fortran_eos/src/compose_eos_lower_crust_native.f90','fortran_eos/src/compose_eos_crust_attachment.f90','fortran_eos/tests/crust_attachment_driver.f90']]
cmd=['gfortran','-std=f2018','-O2','-ffree-line-length-none','-fimplicit-none','-fcheck=all','-J'+str(BUILD),'-I'+str(BUILD),*[str(p) for p in SOURCES],'-o',str(EXE)]
p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True)
(OUT/'build.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
n0=2.69403e-4

def run(n,anchor=n0,aerr=0.,r=0.,tol=1e-10,ext=False,check=True):
 p=subprocess.run([str(EXE),'state',str(n),str(anchor),str(aerr),str(r),str(tol),'extended' if ext else 'nominal'],capture_output=True,text=True,cwd=ROOT)
 if not check:return p
 assert p.returncode==0,p.stdout+p.stderr
 return np.fromstring(p.stdout,sep=' ')

rows=[];max_oracle=0;max_firstlaw=0;max_pressure_derivative=0;max_sound_difference=0
source_anchor=source_state(n0,thermodynamic_mode='source_fit')
e0=source_anchor['epsilon']/n0
for n in np.geomspace(n0,.002,41):
 row=run(n);rows.append(row.tolist());orig=source_state(n,thermodynamic_mode='source_fit')
 integ,err=quad(lambda x:source_state(n0*np.exp(x),thermodynamic_mode='source_fit')['pressure']/(n0*np.exp(x)),0.,np.log(n/n0),epsabs=2e-12,epsrel=2e-13)
 assert abs(row[2]-orig['pressure'])<2e-16
 assert abs(row[4]-orig['epsilon'])<2e-14
 max_oracle=max(max_oracle,abs(row[1]/n-e0-integ))
 assert abs(row[1]+row[2]-n*row[3])<3e-15
 assert max(abs(row[19:21]))<3e-16
 assert row[8]>0 and 0<row[9]<1
 if n>n0*1.001 and n<.00199:
  h=n*1e-5;lo=run(n-h);hi=run(n+h)
  max_firstlaw=max(max_firstlaw,abs((hi[1]-lo[1])/(2*h)-row[3]))
  max_sound_difference=max(max_sound_difference,abs((hi[2]-lo[2])/(hi[1]-lo[1])-row[9]))
  max_pressure_derivative=max(max_pressure_derivative,abs(n*n*((hi[1]/(n+h))-(lo[1]/(n-h)))/(2*h)-row[2]))
assert max_oracle<4e-11,max_oracle
assert max_firstlaw<3e-7,max_firstlaw
assert max_pressure_derivative<3e-11,max_pressure_derivative
assert max_sound_difference<2e-11,max_sound_difference
at_anchor=run(n0);assert at_anchor[1]==at_anchor[4] and at_anchor[18]==0
sample=run(.002,aerr=.01,r=.04)
f=.04/.96
assert abs(sample[13]-(.01+f*sample[12]))<1e-14
assert abs(sample[14]-(sample[13]+f*sample[2]/sample[0]))<1e-14
# A pressure envelope alone cannot bound a derivative: numerical cs2 is checked
# separately above and no sound-speed uncertainty is invented in the API.
refine=[]
for tol in (1e-6,1e-8,1e-10,1e-12):
 q=run(.002,tol=tol);refine.append({'epsabs_MeV':tol,'energy_per_baryon_MeV':q[1]/q[0],'reported_quadrature_estimate_MeV':q[7]})
assert abs(refine[-1]['energy_per_baryon_MeV']-refine[-2]['energy_per_baryon_MeV'])<2e-10
assert run(.00201,check=False).returncode!=0
extended=run(.01,ext=True);assert 0<extended[9]<1
for kwargs in ({'r':1.},{'r':-.1},{'r':float('nan')},{'aerr':-.1},{'tol':0.},{'anchor':.003}):
 assert run(.001,check=False,**kwargs).returncode!=0,kwargs
seam=subprocess.run([str(EXE),'seam'],capture_output=True,text=True,cwd=ROOT)
assert seam.returncode==0,seam.stdout+seam.stderr
np.savetxt(OUT/'reference_integral_rows.csv',np.array(rows),delimiter=',',header='n,epsilon,P,mu,source_epsilon,source_P,source_thermo_residual,quadrature_estimate,gamma,cs2,composition_mu_residual,anchor_e,integral,energy_model_bound,mu_model_bound,epsilon_model_bound,epsilon_quadrature_estimate,mu_quadrature_estimate,delta_source_e,closure_baryon,closure_charge')
receipt={'status':'PASS','scope':'source-boundary diagnostics and source-energy-anchored analytic pressure integral; not microscopic ETFSI','base':'a49aaf9f82ddfda484e0f965aefdfe73aa7918f2','compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],'command':cmd,'rows':len(rows),'max_independent_quadrature_difference_MeV':max_oracle,'max_firstlaw_residual_MeV':max_firstlaw,'max_pressure_derivative_residual_MeV_fm3':max_pressure_derivative,'max_sound_speed_dimensionless_difference':max_sound_difference,'max_source_energy_correction_MeV':max(abs(row[18]) for row in rows),'max_composition_mu_residual_MeV':max(abs(row[10]) for row in rows),'conditional_example_at_002':{'anchor_absolute_error_MeV':.01,'pressure_same_n_relative_envelope':.04,'e_bound_MeV':sample[13],'mu_bound_MeV':sample[14],'epsilon_bound_MeV_fm3':sample[15],'warning':'conditional illustration; literature P(rho) typical error does not certify a uniform P(n) envelope'},'quadrature_refinement':refine,'seam_fixture':seam.stdout.strip(),'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCES}}
(OUT/'TEST_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
