"""Explicit candidate vacuum-droplet calibration, not Scurto's original fit.

The chosen finite-nucleus functional is the vacuum spherical limit of the
interface energy, with no shell, pairing, curvature, Coulomb exchange or
electron screening. Its applicability to the original mass-fit code is not
established. It is useful as a labeled constructive calibration experiment.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from .gdfm_uniform import baryons, MASS
from .scurto_interface import surface_sigma

AME_SHA256 = "2167f57a2a98331e4649b2dd2b658a9006ed4fba1975729ebfe52a42b4b9218a"
E2 = 1.43996448  # e^2/(4*pi*epsilon0), MeV fm; explicit candidate constant
M_NEUTRON = 939.56536
M_PROTON = 938.27203

def ame_candidate(path: str | Path) -> dict[str, np.ndarray]:
    raw = Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != AME_SHA256:
        raise ValueError("AME2016 input identity mismatch")
    rows = []
    for line_no, line in enumerate(raw.decode("ascii").splitlines()[39:], 40):
        n, z, a = int(line[4:9]), int(line[9:14]), int(line[14:19])
        if n+z != a:
            raise ValueError("invalid AME isotope")
        if n>=8 and z>=8 and "#" not in line[28:42]:
            binding = float(line[54:65])*a/1000
            rows.append((n,z,a,binding,line_no))
    if len(rows) != 2408:
        raise ValueError("unexpected AME candidate count")
    values=np.array(rows)
    return dict(zip(("N","Z","A","binding_MeV","source_line"),values.T))

def finite_size_per_baryon(a, z, density, sigma0, b):
    """Surface and unscreened continuous-charge Coulomb energies, MeV/A."""
    a,z,density=np.broadcast_arrays(np.asarray(a,float),np.asarray(z,float),np.asarray(density,float))
    if (not all(np.all(np.isfinite(x)) for x in (a,z,density)) or np.any(a<=0) or
            np.any(z<=0) or np.any(z>=a) or np.any(density<=0)):
        raise ValueError("finite A>Z>0 and positive density required")
    sigma=np.array([surface_sigma(float(y),sigma0,b) for y in (z/a).flat]).reshape(a.shape)
    radius=(3*a/(4*math.pi*density))**(1/3)
    return 4*math.pi*radius**2*sigma/a, 3*E2*z*z/(5*radius*a)

@dataclass
class DensityGrid:
    a: np.ndarray
    z: np.ndarray
    density: np.ndarray
    bulk_binding_per_baryon: np.ndarray

    def __post_init__(self):
        self.a=np.asarray(self.a,float);self.z=np.asarray(self.z,float)
        self.density=np.asarray(self.density,float)
        self.bulk_binding_per_baryon=np.asarray(self.bulk_binding_per_baryon,float)
        if (self.density.ndim!=1 or len(self.density)<3 or
                np.any(np.diff(self.density)<=0) or self.density[0]<=0 or
                self.a.shape!=self.z.shape or self.a.ndim!=1 or len(self.a)==0 or
                not np.all(np.isfinite(self.a)) or not np.all(np.isfinite(self.z)) or
                np.any(self.a<=0) or np.any(self.z<=0) or np.any(self.z>=self.a) or
                self.bulk_binding_per_baryon.shape!=(len(self.a),len(self.density)) or
                not np.all(np.isfinite(self.bulk_binding_per_baryon)) or
                not np.allclose(np.diff(self.density),np.diff(self.density)[0])):
            raise ValueError("finite uniform increasing density grid and matching arrays required")
        self.radius=(3*self.a[:,None]/(4*math.pi*self.density))**(1/3)
        self.surface_area_per_baryon=4*math.pi*self.radius**2/self.a[:,None]
        self.coulomb_per_baryon=3*E2*self.z[:,None]**2/(5*self.radius*self.a[:,None])

    def minima(self, sigma0: float, b: float):
        sigma=np.array([surface_sigma(float(z/a),sigma0,b) for a,z in zip(self.a,self.z)])
        energy=self.bulk_binding_per_baryon+self.surface_area_per_baryon*sigma[:,None]+self.coulomb_per_baryon
        index=np.argmin(energy,axis=1)
        row=np.arange(len(self.a))
        best=energy[row,index].copy()
        density=self.density[index].copy()
        interior=(index>0)&(index<len(self.density)-1)
        r=row[interior];i=index[interior]
        left,center,right=energy[r,i-1],energy[r,i],energy[r,i+1]
        curvature=left-2*center+right
        good=curvature>0
        shift=np.zeros(len(r))
        shift[good]=(left[good]-right[good])/(2*curvature[good])
        # A discrete minimum has a local quadratic vertex inside +/-half a cell.
        shift=np.clip(shift,-.5,.5)
        best[r]=center+0.5*(right-left)*shift+0.5*curvature*shift**2
        density[r]+=shift*(self.density[1]-self.density[0])
        return -best*self.a,density,~interior

    def coarsen(self):
        return DensityGrid(self.a,self.z,self.density[::2],self.bulk_binding_per_baryon[:,::2])

def prepare_grid(model_id: str, data: dict, density: np.ndarray) -> DensityGrid:
    a,z=data["A"],data["Z"]
    fractions,inverse=np.unique(z/a,return_inverse=True)
    bulk=np.array([[baryons(model_id,float(n),float(y))["epsilon_baryons"]/n-MASS
                    for n in density] for y in fractions])
    return DensityGrid(a,z,density,bulk[inverse])

def fit_surface(grid: DensityGrid, observed_binding, *, weighting="mass_2percent", training=None):
    observed_binding=np.asarray(observed_binding,float)
    if observed_binding.shape!=grid.a.shape or not np.all(np.isfinite(observed_binding)):
        raise ValueError("finite binding vector matching isotope count required")
    if weighting=="mass_2percent":
        scale=.02*((grid.a-grid.z)*M_NEUTRON+grid.z*M_PROTON-observed_binding)
    elif weighting=="binding_2percent":
        scale=.02*observed_binding
    else:
        raise ValueError("unknown candidate weighting")
    if np.any(scale<=0):
        raise ValueError("positive uncertainty required")
    if training is None: training=np.ones(len(grid.a),dtype=bool)
    training=np.asarray(training,bool)
    if training.shape!=grid.a.shape or sum(training)<3:
        raise ValueError("at least three training nuclei required")
    def residual(log_parameters):
        sigma0,b=np.exp(log_parameters[0]),np.exp(log_parameters[1])-16
        prediction,_,_=grid.minima(sigma0,b)
        return ((observed_binding-prediction)/scale)[training]
    low=np.log([.01,.01]);high=np.log([10.,1e6])
    trials=[]
    for initial in ([1.,20.],[1.,100.],[1.,1000.]):
        opt=least_squares(residual,np.log(initial),bounds=(low,high),
                          xtol=1e-11,ftol=1e-11,gtol=1e-13,max_nfev=250)
        trials.append(opt)
    winner=min(trials,key=lambda x:float(x.fun@x.fun))
    sigma0,b=float(np.exp(winner.x[0])),float(np.exp(winner.x[1])-16)
    prediction,density,edge=grid.minima(sigma0,b)
    residuals=observed_binding-prediction
    singular=np.linalg.svd(winner.jac,compute_uv=False)
    mask_bounds=(abs(winner.x-low)<1e-5)|(abs(winner.x-high)<1e-5)
    def metric(mask):
        if not np.any(mask):return None
        r=residuals[mask]
        return {"n":int(sum(mask)),"rms_mass_residual_MeV":float(np.sqrt(np.mean(r*r))),
                "rms_binding_per_baryon_residual_MeV":float(np.sqrt(np.mean((r/grid.a[mask])**2))),
                "max_absolute_mass_residual_MeV":float(np.max(abs(r))),
                "mean_scaled_square":float(np.mean((r/scale[mask])**2))}
    result={"sigma0_MeV_fm2":sigma0,"b":b,"weighting":weighting,
        "training":metric(training),"heldout":metric(~training),"all":metric(np.ones(len(grid.a),bool)),
        "density_min_fm3":float(min(density)),"density_max_fm3":float(max(density)),
        "density_boundary_hits":int(sum(edge)),"parameter_bounds_active":mask_bounds.tolist(),
        "optimizer_success":bool(winner.success),"optimizer_message":winner.message,
        "nfev":int(winner.nfev),
        "jacobian_condition_number":float(singular[0]/singular[-1]) if singular[-1]>np.finfo(float).eps*singular[0] else None,
        "locally_identifiable_two_parameters":bool(singular[-1]>np.finfo(float).eps*singular[0]),
        "multistart_scaled_sumsquares":[float(x.fun@x.fun) for x in trials],
        "parameter_search_bounds":{"sigma0":[.01,10.],"b_plus_16":[.01,1e6]},
        "prediction_binding_MeV":prediction,"optimized_density_fm3":density}
    return result
