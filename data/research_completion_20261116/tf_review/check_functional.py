"""Independent numerical review of finite-range TF functional, not a production gate.
Run: PYTHONPATH=python_frontend python3 data/research_completion_20261116/tf_review/check_functional.py
"""
import hashlib,json,math,pathlib
import numpy as np
from scipy.integrate import quad
from compose_eos import gppva_rmf as rmf
from compose_eos.spatial_tf import CellFunctional,RadialGrid,fermions,HBARC,E2
from compose_eos.cold_leptons import free_fermion,ELECTRON_MASS
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent
records=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
hashes={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'python_frontend/compose_eos/spatial_tf.py',ROOT/'python_frontend/compose_eos/gppva_rmf.py']}
kin=[]
for n in (0.,1e-16,1e-10,1e-6,.003,.08):
 for m in (ELECTRON_MASS,500.,939.):
  k=HBARC*(3*math.pi**2*n)**(1/3)
  factor=1/(math.pi**2*HBARC**3)
  ke=quad(lambda t:t*t*math.hypot(t,m),0,k,epsabs=1e-14,epsrel=5e-13)[0]*factor
  ns=quad(lambda t:t*t*m/math.hypot(t,m),0,k,epsabs=1e-14,epsrel=5e-13)[0]*factor
  ds=quad(lambda t:t**4/math.hypot(t,m)**3,0,k,epsabs=1e-14,epsrel=5e-13)[0]*factor
  got=fermions(n,m)
  errors={q:abs(float(got[q])-v)/max(abs(v),1e-300) for q,v in [('epsilon',ke),('scalar',ns),('ds_dm',ds)]}
  kin.append(dict(n=n,m=m,relative_errors=errors))
for name in rmf.MODELS:
 for dimension in (1,2,3):
  g=RadialGrid(24,dimension);x=g.x
  nn=.035+.008*np.cos(math.pi*x)
  np_=.0016+.0007*np.cos(math.pi*x)
  ne=.0016+.00015*np.cos(math.pi*x);ne+=g.mean(np_-ne)
  radius=28.;f=CellFunctional(name,g,g.mean(nn+np_))
  s=f.evaluate(nn,np_,ne,radius)
  direction=.3+np.sin(1.7*math.pi*x)
  derivatives=[]
  for i,key in enumerate(('mu_n','mu_p','mu_e')):
   h=2e-7;args=[nn.copy(),np_.copy(),ne.copy()]
   args[i]+=h*direction;plus=f.evaluate(*args,radius)['epsilon']
   args[i]-=2*h*direction;minus=f.evaluate(*args,radius)['epsilon']
   fd=(plus-minus)/(2*h);analytic=g.mean(direction*s[key])
   derivatives.append(dict(species=key,finite_difference=fd,analytic=analytic,absolute_error=abs(fd-analytic)))
  h=2e-4
  erp=f.evaluate(nn,np_,ne,radius*math.exp(h))['epsilon']
  erm=f.evaluate(nn,np_,ne,radius*math.exp(-h))['epsilon']
  rfd=(erp-erm)/(2*h);ra=s['d_epsilon_d_log_radius']
  charge=np_-ne;v=s['coulomb'];lap=g.minus_laplacian/radius**2
  sym=g.weights[:,None]*g.minus_laplacian
  null=np.max(abs(g.minus_laplacian@np.ones(g.points)))
  poisson=np.max(abs(lap@v-4*math.pi*E2*(charge-g.mean(charge))))
  ec_field=g.mean(v*(lap@v))/(8*math.pi*E2)
  n=.04;y=.025;b=rmf.baryons(name,n,y)
  un=np.full(g.points,n*(1-y));up=np.full(g.points,n*y)
  u=CellFunctional(name,g,n).evaluate(un,up,up,radius)
  eref=b['epsilon_baryons']+free_fermion(n*y,ELECTRON_MASS)[0]
  records.append(dict(model=name,dimension=dimension,density_derivatives=derivatives,
   log_radius_fd=rfd,log_radius_analytic=ra,log_radius_absolute_error=abs(rfd-ra),
   poisson_residual=poisson,weighted_self_adjoint_error=float(np.max(abs(sym-sym.T))),
   laplacian_constant_residual=float(null),mean_potential=g.mean(v),
   coulomb_energy=s['coulomb_energy'],coulomb_field_energy=ec_field,
   coulomb_energy_identity_error=abs(ec_field-s['coulomb_energy']),
   uniform_energy_error=abs(u['epsilon']-eref),uniform_mu_n_error=float(np.max(abs(u['mu_n']-b['mu_n']))),
   uniform_mu_p_error=float(np.max(abs(u['mu_p']-b['mu_p']))),
   scalar_pde_residual=s['scalar_pde_residual'],vector_pde_residual=s['vector_pde_residual']))
summary=dict(
 kinetic_max_relative_error=max(max(r['relative_errors'].values()) for r in kin),
 density_derivative_max_absolute_error=max(d['absolute_error'] for r in records for d in r['density_derivatives']),
 log_radius_derivative_max_absolute_error=max(r['log_radius_absolute_error'] for r in records),
 uniform_energy_max_absolute_error=max(r['uniform_energy_error'] for r in records),
 uniform_mu_max_absolute_error=max(max(r['uniform_mu_n_error'],r['uniform_mu_p_error']) for r in records),
 poisson_max_absolute_residual=max(r['poisson_residual'] for r in records),
 coulomb_energy_identity_max_absolute_error=max(r['coulomb_energy_identity_error'] for r in records),
 weighted_self_adjoint_max_error=max(r['weighted_self_adjoint_error'] for r in records),
 laplacian_constant_max_residual=max(r['laplacian_constant_residual'] for r in records),
 scalar_pde_max_residual=max(r['scalar_pde_residual'] for r in records),
 vector_pde_max_residual=max(r['vector_pde_residual'] for r in records))
conditions=dict(
 kinetic_quadrature=summary['kinetic_max_relative_error']<2e-7,
 density_derivatives=summary['density_derivative_max_absolute_error']<2e-4,
 log_radius_derivative=summary['log_radius_derivative_max_absolute_error']<2e-7,
 uniform_reduction=summary['uniform_energy_max_absolute_error']<1e-8 and summary['uniform_mu_max_absolute_error']<1e-7,
 poisson=summary['poisson_max_absolute_residual']<1e-10,
 coulomb_positive_and_identity=all(r['coulomb_energy']>=0 for r in records) and summary['coulomb_energy_identity_max_absolute_error']<1e-12,
 weighted_operator=summary['weighted_self_adjoint_max_error']<1e-10 and summary['laplacian_constant_max_residual']<1e-10)
receipt=dict(status='PASS_SCOPED' if all(conditions.values()) else 'FAIL',scope='Analytic kinetic integrals, all8 model spatial energy derivatives and uniform reduction on24 point grids in3 dimensions; no equilibrium/convergence/globalminimum claim',source_hashes=hashes,summary=summary,conditions=conditions,kinetic_cases=kin,spatial_cases=records)
(OUT/'functional_review.json').write_text(json.dumps(receipt,indent=2,default=lambda v:v.item())+'\n')
print(json.dumps({k:receipt[k] for k in ('status','scope','summary','conditions')},indent=2,default=lambda v:v.item()))
assert all(conditions.values())
