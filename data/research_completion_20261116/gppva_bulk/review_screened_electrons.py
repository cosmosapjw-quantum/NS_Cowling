"""Independent directional, constrained-envelope and Poisson review.

Run from repository root with PYTHONPATH=python_frontend. This script makes
no source mutations. It checks all eight models in all three radial geometries.
"""
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
from compose_eos.gppva_rmf import MODELS
from compose_eos.spatial_tf import RadialGrid,CellFunctional,initial_profile,fermions,E2,ELECTRON_MASS,HBARC
from compose_eos.spatial_tf_equilibrium import ReducedCellFunctional,screened_electrons

start=time.monotonic();rows=[];seed=817204
rng=np.random.default_rng(seed)
paths=['python_frontend/compose_eos/spatial_tf_equilibrium.py','python_frontend/compose_eos/spatial_tf.py','python_frontend/compose_eos/gppva_rmf.py']
hashes={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}
for d in [1,2,3]:
 g=RadialGrid(24,d);n=.03;R=23.
 # Smooth positive neutron/proton profiles, including nonuniform electrons.
 np_=n*(.055+.035*np.cos(np.pi*g.x));nn=n*(1+.3*np.cos(np.pi*g.x))-np_
 nn*= (n-g.mean(np_))/g.mean(nn)
 ne,v,mu,error=screened_electrons(g,np_,R)
 lap=g.minus_laplacian/R**2
 # Independently solve augmented Poisson/gauge equations by dense algebra.
 aug=np.block([[lap,np.ones((g.points,1))],[g.weights[None,:],np.zeros((1,1))]])
 direct=np.linalg.solve(aug,np.r_[4*np.pi*E2*(np_-ne),0.])[:-1]
 electron_mu=fermions(ne,np.full_like(ne,ELECTRON_MASS))['ef']-v
 e_c=.5*g.mean((np_-ne)*v)
 field_c=g.mean(v*(lap@v))/(8*np.pi*E2)
 # A charge-preserving electron perturbation must raise the eliminated energy.
 perturb=np.cos(2*np.pi*g.x);perturb-=g.mean(perturb)
 perturbed=ne+.03*min(ne)*perturb/max(abs(perturb))
 def electron_energy(e):
  pot=g.poisson(np_-e,R)
  return g.mean(fermions(e,np.full_like(e,ELECTRON_MASS))['epsilon'])+.5*g.mean((np_-e)*pot)
 e_increase=electron_energy(perturbed)-electron_energy(ne)
 electron=dict(dimension=d,poisson_residual=float(max(abs(lap@v-4*np.pi*E2*(np_-ne)))),neutrality=abs(g.mean(np_-ne)),gauge=abs(g.mean(v)),dense_poisson_difference=float(max(abs(v-direct))),electron_chemical_spread=float(max(abs(electron_mu-mu))),coulomb_charge_field_difference=abs(e_c-field_c),electron_energy_increase=e_increase,weighted_laplacian_symmetry=float(np.max(np.abs(g.weights[:,None]*lap-(g.weights[:,None]*lap).T))),constant_nullspace=float(max(abs(lap@np.ones(g.points)))))
 for model in MODELS:
  f=ReducedCellFunctional(model,g,n)
  x=np.r_[nn*g.sqrtw/n,np_*g.sqrtw/n,np.log(R)]
  value,grad=f.objective(x)
  base=f.evaluate(nn,np_,ne,R)
  def values(z):
   val=f.objective(z)[0]
   vguess,muguess=f._electron_guess
   q=muguess+vguess
   electrons=np.where(q>ELECTRON_MASS,(np.sqrt(np.maximum(q*q-ELECTRON_MASS**2,0.))/HBARC)**3/(3*np.pi**2),0.)
   protons=z[g.points:-1]*n/g.sqrtw
   # Constrained Lagrangian removes the known finite neutrality residual only
   # for derivative diagnosis; raw energy differences remain separately reported.
   correction=mu*g.mean(protons-electrons)/n
   return val,val+correction
  checks=[]
  # Separate neutron, proton (including changed total charge), mixed and radius directions.
  for kind in ['neutron','proton','mixed','radius']:
   direction=np.zeros_like(x)
   if kind=='radius':direction[-1]=1.
   elif kind=='mixed':direction[:-1]=rng.normal(size=len(x)-1)
   elif kind=='neutron':direction[:g.points]=rng.normal(size=g.points)
   else:direction[g.points:-1]=rng.normal(size=g.points)
   direction/=np.linalg.norm(direction)
   h=1e-5
   plus,cp=values(x+h*direction);minus,cm=values(x-h*direction)
   coarse=(plus-minus)/(2*h);corrected_coarse=(cp-cm)/(2*h)
   plus,cp=values(x+h/2*direction);minus,cm=values(x-h/2*direction)
   fine=(plus-minus)/h;corrected_fine=(cp-cm)/h
   finite=(4*fine-coarse)/3;corrected_finite=(4*corrected_fine-corrected_coarse)/3;analytic=float(grad@direction)
   a=f.unpack(x+h*direction);b=f.unpack(x-h*direction)
   dne=(a[2]-b[2])/(2*h)
   dnn=n*direction[:g.points]/g.sqrtw;dnp=n*direction[g.points:-1]/g.sqrtw
   # Chain rule through full functional, including electron response and total neutrality.
   chain=(g.mean((base['mu_n']-f.model['mn'])*dnn+(base['mu_p']-f.model['mn'])*dnp+base['mu_e']*dne)+base['d_epsilon_d_log_radius']*direction[-1])/n
   checks.append(dict(direction=kind,analytic=analytic,coarse_finite_difference=coarse,fine_finite_difference=fine,finite_difference=finite,full_chain_rule=float(chain),fd_error=abs(analytic-finite),chain_error=abs(analytic-chain),neutrality_corrected_fd_error=abs(analytic-corrected_finite),neutrality_corrected_chain_error=abs(analytic-chain-mu*(g.mean(dnp)-g.mean(dne))/n),charge_response_error=abs(g.mean(dne)-g.mean(dnp))))
  rows.append(dict(model=model,dimension=d,electron=electron,derivatives=checks))
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items()),'Code changed during review; rerun.'
summary=dict(max_neutrality_corrected_fd_error=max(c['neutrality_corrected_fd_error'] for r in rows for c in r['derivatives']),max_neutrality_corrected_chain_error=max(c['neutrality_corrected_chain_error'] for r in rows for c in r['derivatives']),max_fd_error=max(c['fd_error'] for r in rows for c in r['derivatives']),max_chain_error=max(c['chain_error'] for r in rows for c in r['derivatives']),max_charge_response_error=max(c['charge_response_error'] for r in rows for c in r['derivatives']),max_electron_chemical_error=max(r['electron']['electron_chemical_spread'] for r in rows),max_dense_poisson_difference=max(r['electron']['dense_poisson_difference'] for r in rows),min_electron_energy_increase=min(r['electron']['electron_energy_increase'] for r in rows))
result=dict(scope='All8 x 3 geometries, four reduced/full-functional directions each; finite sampled implementation review, not global minimum proof',seed=seed,source_sha256=hashes,summary=summary,rows=rows,runtime_seconds=time.monotonic()-start)
Path(sys.argv[1]).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(summary,indent=2))
assert summary['max_neutrality_corrected_fd_error']<2e-6
assert summary['max_neutrality_corrected_chain_error']<2e-7
assert summary['max_fd_error']<2e-5
assert summary['max_chain_error']<2e-6
assert summary['max_charge_response_error']<2e-9
assert summary['max_electron_chemical_error']<2e-7
assert summary['max_dense_poisson_difference']<2e-8
assert summary['min_electron_energy_increase']>0
