"""Screened-electron reduction of the same spatial TF energy functional.

Electrons obey sqrt(kFe(r)^2+me^2)-V(r)=constant exactly at each outer
iteration. This removes weakly weighted electron variables from the density
minimizer without changing the physical model or relaxing KKT tolerances.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import minimize
from .spatial_tf import CellFunctional,RadialGrid,initial_profile,_band_solve,fermions,E2,HBARC,ELECTRON_MASS,MUON_MASS


def screened_electrons(grid,protons,radius,initial=None):
    """Solve local zero-T electron equilibrium, global neutrality and Poisson."""
    np_=np.asarray(protons,float);average=grid.mean(np_)
    if np.any(np_<0):raise ValueError('nonnegative proton profile required')
    if average==0:return np.zeros_like(np_),np.zeros_like(np_),ELECTRON_MASS,0.
    lap=grid.minus_laplacian/radius**2
    if initial is None:
        mu=math.hypot(HBARC*(3*np.pi**2*average)**(1/3),ELECTRON_MASS)
        v=grid.poisson(np_-average,radius)
    else:v,mu=initial[0].copy(),float(initial[1])
    def terms(v,mu):
        chemical=mu+v;k=np.sqrt(np.maximum(chemical**2-ELECTRON_MASS**2,0.))
        k=np.where(chemical>ELECTRON_MASS,k,0.)
        ne=(k/HBARC)**3/(3*np.pi**2)
        derivative=np.where(chemical>ELECTRON_MASS,chemical*k/(np.pi**2*HBARC**3),0.)
        residual=lap@v-4*np.pi*E2*(np_-ne)
        return ne,derivative,residual
    error=math.inf
    for _ in range(45):
        ne,d,res=terms(v,mu);error=float(np.max(abs(res)))
        if error<2e-12 and abs(grid.mean(v))<1e-11:break
        jac=lap+np.diag(4*np.pi*E2*d)
        z=_band_solve(jac,-res)
        # jac*1 = 4*pi*e²*d, so Schur elimination of the gauge is exact.
        dmu=grid.mean(z)+grid.mean(v);dv=z-dmu
        before=np.linalg.norm(res)
        for scale in (1.,.5,.25,.125,.0625,.03125):
            tv,tm=v+scale*dv,mu+scale*dmu
            if np.linalg.norm(terms(tv,tm)[2])<before:
                v,mu=tv,tm;break
        else:raise RuntimeError('screened-electron Newton line search failed')
    if error>2e-10:raise RuntimeError(f'screened electron PDE residual {error}')
    ne,_,res=terms(v,mu)
    if abs(grid.mean(np_-ne))>2e-11:raise RuntimeError('screened electron neutrality failure')
    return ne,v,mu,float(np.max(abs(res)))


class ReducedCellFunctional(CellFunctional):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);self._electron_guess=None

    def unpack(self,x):
        nn,np_=np.split(x[:-1],2);nn=nn*self.nbar/self.grid.sqrtw;np_=np_*self.nbar/self.grid.sqrtw
        radius=float(np.exp(x[-1]))
        ne,v,mu,res=screened_electrons(self.grid,np_,radius,self._electron_guess)
        self._electron_guess=(v.copy(),mu)
        return nn,np_,ne,radius,mu,res

    def objective(self,x):
        if self._x is not None and np.array_equal(x,self._x):return self._result
        nn,np_,ne,radius,mu,res=self.unpack(x);s=self.evaluate(nn,np_,ne,radius)
        m=self.model['mn'];w=self.grid.sqrtw
        value=(s['epsilon']-m*s['mean_baryon_density'])/self.nbar
        grad=np.r_[w*(s['mu_n']-m),w*(s['mu_p']+mu-m),s['d_epsilon_d_log_radius']/self.nbar]
        self._x=x.copy();self._result=(value,grad);return self._result


def solve_screened_cell(model_id,nbar,*,points=48,dimension=3,seed='droplet',radius=35.,
               radius_bounds=(5.,180.),maxiter=900,initial=None,ftol=1e-13):
    grid=RadialGrid(points,dimension);f=ReducedCellFunctional(model_id,grid,nbar);w=grid.sqrtw
    if initial is None:nn,np_,_=initial_profile(grid,nbar,seed,radius)
    else:
        nn,np_=(np.asarray(initial[k],float) for k in ('nn','np'))
        radius=float(initial.get('radius',radius))
    x0=np.r_[nn*w/nbar,np_*w/nbar,np.log(radius)]
    constraint=np.r_[w,w,0.]
    upper=np.r_[.4*w/nbar,.4*w/nbar]
    bounds=list(zip(np.zeros(2*points),upper))+[(np.log(radius_bounds[0]),np.log(radius_bounds[1]))]
    opt=minimize(f.objective,x0,jac=True,method='SLSQP',bounds=bounds,
        constraints={'type':'eq','fun':lambda x:constraint@x-1,'jac':lambda x:constraint},
        options={'ftol':ftol,'maxiter':maxiter,'disp':False})
    nn,np_,ne,radius,muq,eres=f.unpack(opt.x);s=f.evaluate(nn,np_,ne,radius)
    active=[nn>max(1e-10,nbar*1e-7),np_>max(1e-10,nbar*1e-7)]
    coeff=np.r_[grid.weights[active[0]],grid.weights[active[1]]]
    chemical=np.r_[s['mu_n'][active[0]],(s['mu_p']+muq)[active[1]]]
    mub=float(coeff@chemical/sum(coeff))
    kkt=[]
    for mask,mu in zip(active,[s['mu_n'],s['mu_p']+muq]):
        kkt.extend(abs(mu[mask]-mub));kkt.extend(np.maximum(mub-mu[~mask],0.))
    kkt.extend(abs(s['mu_e']-muq))
    bound=bool(min(abs(np.log(radius/radius_bounds[0])),abs(np.log(radius/radius_bounds[1])))<1e-5)
    density_bound=bool(max(max(nn),max(np_))>.4*(1-1e-7))
    contrast=float((max(nn+np_)-min(nn+np_))/nbar)
    result={k:v for k,v in s.items() if np.isscalar(v)}
    result.update(model_id=model_id,n_b=nbar,points=points,dimension=dimension,seed=seed,
        radius_fm=radius,optimizer_success=bool(opt.success),optimizer_message=str(opt.message),iterations=opt.nit,
        evaluations=opt.nfev,energy_per_baryon=s['epsilon']/nbar,
        mu_B=mub,mu_Q_lepton=muq,pressure=float(nbar*mub-s['epsilon']),
        Y_n=grid.mean(nn)/nbar,Y_p=grid.mean(np_)/nbar,Y_e=grid.mean(ne)/nbar,Y_mu=0.,
        baryon_constraint_error=abs(s['mean_baryon_density']-nbar),charge_constraint_error=abs(s['mean_charge_density']),
        max_chemical_KKT_residual_MeV=float(max(kkt)),electron_pde_residual=eres,
        radius_stationarity_MeV=abs(s['d_epsilon_d_log_radius']/nbar),radius_at_search_bound=bound,
        density_at_search_bound=density_bound,baryon_density_contrast=contrast,
        muon_threshold_crossed=s['max_local_electron_Fermi_energy']>MUON_MASS,
        electron_algorithm='SCREENED_ELECTRON_ENVELOPE_REDUCTION',
        profile={k:v.tolist() for k,v in dict(x=grid.x,weights=grid.weights,nn=nn,np=np_,ne=ne,
            phi=s['phi'],omega=s['omega'],rho=s['rho'],coulomb=s['coulomb'],
            mu_n=s['mu_n'],mu_p=s['mu_p'],mu_e=s['mu_e']).items()})
    result['numerical_status']='STATIONARY' if (opt.success and max(kkt)<.02 and result['radius_stationarity_MeV']<2e-4 and not bound and not density_bound) else 'REQUIRES_REFINEMENT'
    if contrast<1e-4 and max(kkt)<.02 and not density_bound:result['numerical_status']='UNIFORM_STATIONARY'
    return result
