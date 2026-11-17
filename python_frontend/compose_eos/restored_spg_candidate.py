"""Explicit species-rest-energy-restored SPG candidate, restricted zero-p gas.

This changes the beta-equilibrium functional; it is not an energy-reference
gauge transformation and is not the originally calibrated SPG model. Bulk
kinetic masses, couplings and the explicit candidate surface remain unchanged.
The restricted gas branch is admitted only when proton insertion is costly.
"""
from dataclasses import asdict, dataclass
from copy import deepcopy
from functools import lru_cache
import math
import numpy as np
from scipy.optimize import minimize,least_squares
from .candidate_gdfm_cld import reduced_cell,MASS,default_surface_parameters,SurfaceParameters

MN=939.56542052
MP=938.27208816

def _local_minimum_check(gradient,z,lo,hi):
    """Refuse artificial boundaries and negative curvature on the free subspace.

    At t=1, only strict complementarity (g_t < -1e-6) permits the
    (ni,yi) tangent subspace. At a marginal boundary the inward t direction
    is checked with a one-sided Hessian column. No global minimum is proved.
    """
    z=np.asarray(z,float)
    active=[j for j in (1,2) if min(z[j]-lo[j],hi[j]-z[j])<1e-6*(hi[j]-lo[j])]
    if z[0]-lo[0]<1e-6*(hi[0]-lo[0]):active.append(0)
    if active:raise RuntimeError("artificial search bound active: "+str(active))
    physical_boundary=z[0]>1-1e-8
    strict_boundary=physical_boundary and gradient(z)[0]<-1e-6
    free=[1,2] if strict_boundary else [0,1,2]
    steps=np.array([2e-6,2e-7,2e-7]);columns=[]
    for j in free:
        a=z.copy();b=z.copy()
        if j==0 and physical_boundary:
            h=min(steps[j],.25*(z[j]-lo[j]));b[j]-=h
            columns.append((gradient(a)[free]-gradient(b)[free])/h)
        else:
            h=min(steps[j],.25*(z[j]-lo[j]),.25*(hi[j]-z[j]))
            a[j]+=h;b[j]-=h
            columns.append((gradient(a)[free]-gradient(b)[free])/(2*h))
    matrix=np.asarray(columns).T;scale=np.array([1.,.1,.1])[free]
    scaled=(matrix+matrix.T)/2*scale[:,None]*scale[None,:]
    if not np.isfinite(scaled).all():raise RuntimeError("nonfinite local curvature")
    mineig=float(np.linalg.eigvalsh(scaled)[0])
    if mineig< -1e-6:raise RuntimeError("negative local curvature: "+str(mineig))
    return dict(minimum_scaled_hessian_eigenvalue=mineig,
                curvature_free_coordinates=free,strict_zero_gas_complementarity=bool(strict_boundary),
                artificial_search_bound_active=False)


@dataclass(frozen=True,init=False)
class RestoredSPGCandidate:
    model: str
    surface: SurfaceParameters
    mn: float
    mp: float

    def __init__(self,model,*,surface=None,mn=MN,mp=MP):
        surface=default_surface_parameters(model) if surface is None else surface
        if not isinstance(surface,SurfaceParameters) or surface.model_id!=model:
            raise ValueError("surface/bulk model mismatch")
        if not all(math.isfinite(v) and v>0 for v in (mn,mp)):
            raise ValueError("finite positive rest energies required")
        object.__setattr__(self,"model",model)
        object.__setattr__(self,"surface",surface)
        object.__setattr__(self,"mn",float(mn))
        object.__setattr__(self,"mp",float(mp))

    @property
    def dn(self):return self.mn-MASS

    @property
    def dp(self):return self.mp-MASS

    @property
    def difference(self):return self.dp-self.dn

    def state(self,n):
        """Return an independent state; cached nested records stay private."""
        return deepcopy(self._cached_state(float(n)))

    def cell(self,n,z):
        e,g,d=reduced_cell(self.model,n,z,surface=self.surface)
        t,ni,y=z;diff=self.difference
        return e+self.dn+diff*t*y,g+np.array([diff*y,0.,diff*t]),d

    @lru_cache(maxsize=20000)
    def _cached_state(self,n):
        n=float(n)
        if not math.isfinite(n) or not 1e-5<=n<=.03:
            raise ValueError("declared candidate support is 1e-5<=n<=.03; proton-gas KKT also required")
        lo=np.array([1e-7,max(.035,n*(1+1e-6)),.001]);hi=np.array([1.,.25,.499])
        fun=lambda z:self.cell(n,z)[:2]
        records=[]
        for initial in ([.9,.16,.3],[.5,.15,.2],[.2,.14,.12]):
            opt=minimize(fun,initial,jac=True,method="L-BFGS-B",bounds=list(zip(lo,hi)),
                options={"ftol":1e-14,"gtol":1e-10,"maxiter":500,"maxls":50})
            z=opt.x.copy();boundary=z[0]>1-1e-6 and fun(z)[1][0]<0
            if boundary:
                fit=least_squares(lambda q:fun([1.,*q])[1][1:],z[1:],bounds=(lo[1:]+1e-12,hi[1:]-1e-12),
                    jac="3-point",xtol=1e-12,ftol=1e-12,gtol=1e-12)
                z=np.r_[1.,fit.x]
            else:
                fit=least_squares(lambda q:fun(q)[1],np.clip(z,lo+1e-10,hi-1e-10),
                    bounds=(lo+1e-12,hi-1e-12),jac="3-point",xtol=1e-12,ftol=1e-12,gtol=1e-12)
                z=fit.x
            e,g,d=self.cell(n,z)
            projected=g.copy()
            if z[0]>1-1e-8:projected[0]=max(g[0],0.)
            mu=d["state"]["mu_n_gas"]+self.dn+z[0]*g[0]
            eps=n*(MASS+e)
            margin=d["state"]["gas_proton_insertion_margin_MeV"]+self.difference
            records.append(dict(n_B=n,coordinates=z.tolist(),energy_per_baryon=MASS+e,
                epsilon=eps,pressure=n*mu-eps,mu_B=mu,Y_e=float(z[0]*z[2]),
                Y_n_free=float(1-z[0]),gradient=g.tolist(),
                projected_gradient_max=float(max(abs(projected))),
                gas_proton_insertion_margin_MeV=margin,
                muon_margin_MeV=d["state"]["muon_insertion_margin_MeV"],
                initial_minimize_success=bool(opt.success),polish_success=bool(fit.success)))
        winner=min(records,key=lambda r:r["energy_per_baryon"])
        if winner["projected_gradient_max"]>1e-6 or winner["gas_proton_insertion_margin_MeV"]<0 or winner["muon_margin_MeV"]<0 or not winner["polish_success"]:
            raise RuntimeError("unresolved restricted restored candidate: "+str(winner))
        stability=_local_minimum_check(lambda z:fun(z)[1],winner["coordinates"],lo,hi)
        return dict(winner,**stability,all_three_attempts=records,surface_parameters=asdict(self.surface),
                    functional_variant="species_rest_energy_restored_before_beta_minimization",
                    rest_energies_MeV={"mn":self.mn,"mp":self.mp,"kinetic_common_M":MASS})
