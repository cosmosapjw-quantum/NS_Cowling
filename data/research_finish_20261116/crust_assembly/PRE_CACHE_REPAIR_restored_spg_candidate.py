"""Explicit species-rest-energy-restored SPG candidate, restricted zero-p gas.

This changes the beta-equilibrium functional; it is not an energy-reference
gauge transformation and is not the originally calibrated SPG model. Bulk
kinetic masses, couplings and the explicit candidate surface remain unchanged.
The restricted gas branch is admitted only when proton insertion is costly.
"""
from dataclasses import asdict
from functools import lru_cache
import math
import numpy as np
from scipy.optimize import minimize,least_squares
from .candidate_gdfm_cld import reduced_cell,MASS,default_surface_parameters

MN=939.56542052
MP=938.27208816

class RestoredSPGCandidate:
    def __init__(self,model,*,surface=None,mn=MN,mp=MP):
        self.model=model
        self.surface=default_surface_parameters(model) if surface is None else surface
        if self.surface.model_id!=model:raise ValueError("surface/bulk model mismatch")
        if not all(math.isfinite(v) and v>0 for v in (mn,mp)):
            raise ValueError("finite positive rest energies required")
        self.mn=float(mn);self.mp=float(mp)
        self.dn=self.mn-MASS;self.dp=self.mp-MASS;self.difference=self.dp-self.dn

    def cell(self,n,z):
        e,g,d=reduced_cell(self.model,n,z,surface=self.surface)
        t,ni,y=z;diff=self.difference
        return e+self.dn+diff*t*y,g+np.array([diff*y,0.,diff*t]),d

    @lru_cache(maxsize=20000)
    def state(self,n):
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
        return dict(winner,all_three_attempts=records,surface_parameters=asdict(self.surface),
                    functional_variant="species_rest_energy_restored_before_beta_minimization",
                    rest_energies_MeV={"mn":self.mn,"mp":self.mp,"kinetic_common_M":MASS})
