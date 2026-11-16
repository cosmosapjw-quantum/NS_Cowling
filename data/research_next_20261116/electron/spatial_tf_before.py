"""Cold radial finite-range RMF Thomas--Fermi cells in MeV and fm.

The constrained variables are local neutron, proton and electron densities.
Meson and Coulomb fields are solved spatially, with natural Neumann boundaries.
Geometries d=1,2,3 describe slabs, cylinders and spheres; dense and hollow seeds
are compared. This module does not substitute a fitted surface tension for TF.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import minimize
from . import gppva_rmf as rmf
from .cold_leptons import HBARC, ELECTRON_MASS, MUON_MASS

E2 = 1.43996448
HC3 = HBARC**3


def fermions(n, mass):
    """Vectorized spin-two zero-T kinetic energy, scalar density and derivatives."""
    n,mass=np.broadcast_arrays(np.asarray(n,float),np.asarray(mass,float))
    if np.any(n<0) or np.any(mass<=0) or not np.all(np.isfinite(n+mass)):
        raise ValueError('finite n>=0 and mass>0 required')
    k=HBARC*np.cbrt(3*np.pi**2*n);x=k/mass; ef=np.hypot(k,mass)
    small=x<.08
    with np.errstate(divide='ignore',invalid='ignore'):
        log=np.arcsinh(x)
        scalar=mass*(k*ef-mass**2*log)/(2*np.pi**2*HC3)
        ds=(k*ef+2*mass**2*k/ef-3*mass**2*log)/(2*np.pi**2*HC3)
        energy=(k*ef*(2*k*k+mass*mass)-mass**4*log)/(8*np.pi**2*HC3)
    scalar=np.where(small,n*(1-3*x*x/10+9*x**4/56-5*x**6/48+105*x**8/1408),scalar)
    ds=np.where(small,n/mass*(3*x*x/5-9*x**4/14+5*x**6/8-105*x**8/176),ds)
    energy=np.where(small,n*mass*(1+3*x*x/10-3*x**4/56+x**6/48-15*x**8/1408),energy)
    return dict(epsilon=energy,scalar=scalar,ds_dm=ds,ef=ef)


@dataclass
class RadialGrid:
    points: int=48
    dimension: int=3

    def __post_init__(self):
        if self.points<8 or self.dimension not in (1,2,3):
            raise ValueError('at least8 radial cells and dimension1,2,3 required')
        e=np.linspace(0.,1.,self.points+1); self.x=(e[:-1]+e[1:])/2
        self.weights=np.diff(e**self.dimension)
        self.sqrtw=np.sqrt(self.weights)
        a=np.zeros((self.points,self.points))
        for i in range(self.points-1):
            flux=self.dimension*e[i+1]**(self.dimension-1)/(self.x[i+1]-self.x[i])
            a[i,i]+=flux/self.weights[i];a[i,i+1]-=flux/self.weights[i]
            a[i+1,i+1]+=flux/self.weights[i+1];a[i+1,i]-=flux/self.weights[i+1]
        self.minus_laplacian=a
        sym=self.sqrtw[:,None]*a/self.sqrtw[None,:]
        vals,vecs=np.linalg.eigh(sym)
        self.poisson_inverse=(vecs[:,1:]/vals[1:])@vecs[:,1:].T

    def mean(self,v): return float(self.weights@v)

    def poisson(self,charge,radius):
        charge=np.asarray(charge,float)-self.mean(charge)
        v=(self.poisson_inverse@(self.sqrtw*charge))/self.sqrtw
        v*=4*np.pi*E2*radius**2
        return v-self.mean(v)


def _band_solve(matrix,rhs):
    n=len(rhs); band=np.zeros((3,n))
    band[0,1:]=np.diag(matrix,1);band[1]=np.diag(matrix);band[2,:-1]=np.diag(matrix,-1)
    return solve_banded((1,1),band,rhs,check_finite=False)


def _coupled_band_solve(jw,jr,cross,rw,rr):
    """Solve coupled nearest-neighbor fields in interleaved (w_i,r_i) order.

    Local cross derivatives occupy bands +/-1; same-field spatial neighbors
    occupy bands +/-2. No symmetry or positive definiteness is assumed.
    """
    n=len(rw)
    band=np.zeros((5,2*n))
    band[2,0::2]=np.diag(jw);band[2,1::2]=np.diag(jr)
    band[1,1::2]=cross;band[3,0::2]=cross
    band[0,2::2]=np.diag(jw,1);band[0,3::2]=np.diag(jr,1)
    band[4,:-2:2]=np.diag(jw,-1);band[4,1:-2:2]=np.diag(jr,-1)
    rhs=np.empty(2*n);rhs[0::2]=rw;rhs[1::2]=rr
    delta=solve_banded((2,2),band,rhs,check_finite=False)
    return delta[0::2],delta[1::2]


class CellFunctional:
    """Eliminate mean fields, retaining the analytic density/radius gradient."""
    def __init__(self,model_id:str,grid:RadialGrid,mean_density:float):
        self.model_id=model_id;self.model=rmf.model(model_id);self.grid=grid
        self.nbar=float(mean_density)
        if not math.isfinite(self.nbar) or not 0<self.nbar<.3:
            raise ValueError('mean density must lie in (0,.3) fm^-3')
        self._fields=None;self._x=None;self._result=None

    def fields(self,nn,np_,radius):
        n=nn+np_;c=rmf.couplings(self.model_id,n);g=self.grid
        lap=HBARC**2*g.minus_laplacian/radius**2
        if self._fields is None:
            phi=c['gs']*n*HC3/self.model['ms']**2
            w=c['gw']*n*HC3/self.model['mw']**2
            r=c['gr']*(np_-nn)*HC3/(2*self.model['mr']**2)
        else:phi,w,r=(v.copy() for v in self._fields)
        phi=np.clip(phi,0.,.95*np.minimum(self.model['mn'],self.model['mp'])/c['gs'])
        def scalar_res(v):
            mn=self.model['mn']-c['gs']*v;mp=self.model['mp']-c['gs']*v
            kn,kp=fermions(nn,mn),fermions(np_,mp)
            pot=rmf.scalar_potential(self.model_id,v)
            res=lap@v+pot['dU']-c['gs']*(kn['scalar']+kp['scalar'])*HC3
            jac=lap+np.diag(pot['ddU']+c['gs']**2*(kn['ds_dm']+kp['ds_dm'])*HC3)
            return res,jac
        scalar_error=math.inf
        for _ in range(45):
            res,jac=scalar_res(phi);scalar_error=np.max(abs(res))/HC3
            if scalar_error<2e-10:break
            step=_band_solve(jac,-res);before=np.linalg.norm(res)
            for scale in (.1**0,.5,.25,.125,.0625,.03125,.015625,.0078125):
                trial=phi+scale*step
                if np.any(trial<0) or np.any(c['gs']*trial>=.999*min(self.model['mn'],self.model['mp'])):continue
                if np.linalg.norm(scalar_res(trial)[0])<before:
                    phi=trial;break
            else:raise RuntimeError('scalar field Newton line search failed')
        if scalar_error>2e-8:raise RuntimeError(f'scalar PDE residual {scalar_error}')
        n3=np_-nn
        def vector_res(w,r):
            p=rmf.vector_potential(self.model_id,w,r)
            resw=lap@w+p['dw']-c['gw']*n*HC3
            resr=lap@r+p['dr']-c['gr']*n3*HC3/2
            return resw,resr,p
        vector_error=math.inf
        for _ in range(30):
            rw,rr,p=vector_res(w,r);vector_error=max(np.max(abs(rw)),np.max(abs(rr)))/HC3
            if vector_error<2e-10:break
            jw=lap+np.diag(np.broadcast_to(p['dww'],n.shape));jr=lap+np.diag(np.broadcast_to(p['drr'],n.shape))
            cross=np.diag(np.broadcast_to(p['dwr'],n.shape))
            if not np.any(cross):dw=_band_solve(jw,-rw);dr=_band_solve(jr,-rr)
            else:
                dw,dr=_coupled_band_solve(jw,jr,np.diag(cross),-rw,-rr)
            before=np.linalg.norm(np.r_[rw,rr])
            for scale in (1.,.5,.25,.125,.0625,.03125):
                tw,tr=w+scale*dw,r+scale*dr
                a,b,_=vector_res(tw,tr)
                if np.linalg.norm(np.r_[a,b])<before:w,r=tw,tr;break
            else:raise RuntimeError('vector field Newton line search failed')
        if vector_error>2e-8:raise RuntimeError(f'vector PDE residual {vector_error}')
        self._fields=(phi.copy(),w.copy(),r.copy())
        return phi,w,r,c,scalar_error,vector_error

    def evaluate(self,nn,np_,ne,radius):
        if not math.isfinite(radius) or radius<=0:raise ValueError('positive finite cell radius required')
        nn,np_,ne=(np.asarray(v,float) for v in (nn,np_,ne))
        if any(v.shape!=(self.grid.points,) for v in (nn,np_,ne)):raise ValueError('density shape mismatch')
        if np.any(nn<0) or np.any(np_<0) or np.any(ne<0):raise ValueError('nonnegative densities required')
        phi,w,r,c,serr,verr=self.fields(nn,np_,radius);g=self.grid;n=nn+np_;n3=np_-nn
        kn=fermions(nn,self.model['mn']-c['gs']*phi)
        kp=fermions(np_,self.model['mp']-c['gs']*phi)
        ke=fermions(ne,np.full_like(ne,ELECTRON_MASS))
        ns=kn['scalar']+kp['scalar'];u=rmf.scalar_potential(self.model_id,phi);vp=rmf.vector_potential(self.model_id,w,r)
        v=g.poisson(np_-ne,radius)
        grad=lambda field:g.mean(field*(g.minus_laplacian@field))/(2*HBARC*radius**2)
        gs,gw,gr=grad(phi),grad(w),grad(r)
        ec=.5*g.mean((np_-ne)*v)
        energy=g.mean(kn['epsilon']+kp['epsilon']+ke['epsilon']+u['U']/HC3+c['gw']*w*n+c['gr']*r*n3/2-vp['W']/HC3)+gs-gw-gr+ec
        rearr=c['dgw']*w*n+c['dgr']*r*n3/2-c['dgs']*phi*ns
        mun=kn['ef']+c['gw']*w-c['gr']*r/2+rearr
        mup=kp['ef']+c['gw']*w+c['gr']*r/2+rearr+v
        mue=ke['ef']-v
        dlogr=2*(-gs+gw+gr+ec)
        return dict(epsilon=energy,mu_n=mun,mu_p=mup,mu_e=mue,
            phi=phi,omega=w,rho=r,coulomb=v,rearrangement=rearr,
            kinetic_n=kn['epsilon'],kinetic_p=kp['epsilon'],kinetic_e=ke['epsilon'],
            gradient_sigma=gs,gradient_omega=gw,gradient_rho=gr,coulomb_energy=ec,
            d_epsilon_d_log_radius=dlogr,scalar_pde_residual=serr,vector_pde_residual=verr,
            mean_baryon_density=g.mean(n),mean_charge_density=g.mean(np_-ne),
            max_local_electron_Fermi_energy=float(np.max(ke['ef'])))

    def unpack(self,x):
        nn,np_,ne=np.split(x[:-1],3)
        return nn*self.nbar/self.grid.sqrtw,np_*self.nbar/self.grid.sqrtw,ne*self.nbar/self.grid.sqrtw,float(np.exp(x[-1]))

    def objective(self,x):
        if self._x is not None and np.array_equal(x,self._x):return self._result
        nn,np_,ne,radius=self.unpack(x);state=self.evaluate(nn,np_,ne,radius)
        mref=self.model['mn'];w=self.grid.sqrtw
        objective=(state['epsilon']-mref*state['mean_baryon_density'])/self.nbar
        gradient=np.r_[w*(state['mu_n']-mref),w*(state['mu_p']-mref),w*state['mu_e'],state['d_epsilon_d_log_radius']/self.nbar]
        self._x=x.copy();self._result=(objective,gradient)
        return self._result


def initial_profile(grid,nbar,kind='droplet',radius=35.):
    if kind=='uniform':
        yp=.035;nn=np.full(grid.points,nbar*(1-yp));np_=np.full(grid.points,nbar*yp)
    else:
        f=min(.65,max((2/grid.points)**grid.dimension,.45*nbar/.13));x0=f**(1/grid.dimension)
        shape=1/(1+np.exp(np.clip((grid.x-x0)*radius/.7,-60,60)))
        if kind=='hole':shape=1-shape
        elif kind!='droplet':raise ValueError('seed is uniform/droplet/hole')
        np_=.10*nbar*shape/max(grid.mean(shape),1e-12)
        nn=nbar-np_
        nn=.55*nbar+(.45*nbar-grid.mean(np_))*shape/max(grid.mean(shape),1e-12)
        nn*= (nbar-grid.mean(np_))/grid.mean(nn)
    ne=np.full(grid.points,grid.mean(np_))
    return nn,np_,ne


def solve_cell(model_id,nbar,*,points=48,dimension=3,seed='droplet',radius=35.,
               radius_bounds=(5.,100.),maxiter=500,initial=None,ftol=2e-10):
    """Constrained local npe equilibrium and radius minimization in a WS cell."""
    grid=RadialGrid(points,dimension);func=CellFunctional(model_id,grid,nbar);w=grid.sqrtw
    if initial is None:nn,np_,ne=initial_profile(grid,nbar,seed,radius)
    else:
        nn,np_,ne=(np.asarray(initial[k],float) for k in ('nn','np','ne'))
        radius=float(initial.get('radius',radius))
    x0=np.r_[nn*w/nbar,np_*w/nbar,ne*w/nbar,np.log(radius)]
    constraint=np.zeros((2,3*points+1));constraint[0,:points]=w;constraint[0,points:2*points]=w
    constraint[1,points:2*points]=w;constraint[1,2*points:3*points]=-w
    target=np.array([1.,0.])
    # Nuclear support is an explicit search domain, not clipping after solving.
    upper=np.r_[.40*w/nbar,.40*w/nbar,.10*w/nbar]
    bounds=list(zip(np.zeros(3*points),upper))+[(np.log(radius_bounds[0]),np.log(radius_bounds[1]))]
    opt=minimize(func.objective,x0,jac=True,method='SLSQP',bounds=bounds,
        constraints={'type':'eq','fun':lambda x:constraint@x-target,'jac':lambda x:constraint},
        options={'ftol':ftol,'maxiter':maxiter,'disp':False})
    nn,np_,ne,radius=func.unpack(opt.x);state=func.evaluate(nn,np_,ne,radius)
    active=[nn>max(1e-10,nbar*1e-7),np_>max(1e-10,nbar*1e-7),ne>max(1e-10,nbar*1e-7)]
    design=[];observed=[];weights=[]
    for mask,mu,coeff in zip(active,[state['mu_n'],state['mu_p'],state['mu_e']],[(1,0),(1,-1),(0,1)]):
        design.extend([coeff]*int(sum(mask)));observed.extend(mu[mask]);weights.extend(grid.sqrtw[mask])
    design=np.asarray(design);observed=np.asarray(observed);weights=np.asarray(weights)
    mub,muq=np.linalg.lstsq(design*weights[:,None],observed*weights,rcond=None)[0]
    kkt=[]
    for mask,mu,eq in zip(active,[state['mu_n'],state['mu_p'],state['mu_e']],[mub,mub-muq,muq]):
        kkt.extend(abs(mu[mask]-eq));kkt.extend(np.maximum(eq-mu[~mask],0.))
    radius_at_bound=min(abs(np.log(radius/radius_bounds[0])),abs(np.log(radius/radius_bounds[1])))<1e-5
    contrast=float((max(nn+np_)-min(nn+np_))/nbar)
    result={k:v for k,v in state.items() if np.isscalar(v)}
    result.update(model_id=model_id,n_b=nbar,points=points,dimension=dimension,seed=seed,
        radius_fm=radius,optimizer_success=bool(opt.success),optimizer_message=str(opt.message),iterations=opt.nit,
        evaluations=opt.nfev,energy_per_baryon=state['epsilon']/nbar,
        mu_B=float(mub),mu_Q_lepton=float(muq),pressure=float(nbar*mub-state['epsilon']),
        Y_n=grid.mean(nn)/nbar,Y_p=grid.mean(np_)/nbar,Y_e=grid.mean(ne)/nbar,Y_mu=0.,
        baryon_constraint_error=abs(state['mean_baryon_density']-nbar),
        charge_constraint_error=abs(state['mean_charge_density']),
        max_chemical_KKT_residual_MeV=float(max(kkt)),
        radius_stationarity_MeV=abs(state['d_epsilon_d_log_radius']/nbar),
        radius_at_search_bound=bool(radius_at_bound),baryon_density_contrast=contrast,
        muon_threshold_crossed=state['max_local_electron_Fermi_energy']>MUON_MASS,
        profile={k:v.tolist() for k,v in dict(x=grid.x,weights=grid.weights,nn=nn,np=np_,ne=ne,
            phi=state['phi'],omega=state['omega'],rho=state['rho'],coulomb=state['coulomb'],
            mu_n=state['mu_n'],mu_p=state['mu_p'],mu_e=state['mu_e']).items()})
    result['numerical_status']='STATIONARY' if (opt.success and max(kkt)<.02 and result['radius_stationarity_MeV']<2e-4 and not radius_at_bound) else 'REQUIRES_REFINEMENT'
    if contrast<1e-4 and max(kkt)<.02:result['numerical_status']='UNIFORM_STATIONARY'
    return result
