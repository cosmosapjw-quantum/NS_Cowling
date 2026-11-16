"""Candidate spherical GDFM CLD cell; not source-calibrated.

The legacy cell uses pure neutron gas; two_phase_cell includes variational
proton-bearing gas. Both use uniform electrons, no muons and a sharp interface. The fitted
surface parameters are the separate AME2408 vacuum-droplet experiment. This
module does not supply an outer crust or certify the global minimum.
"""
import json
import hashlib
from dataclasses import dataclass, asdict
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, least_squares
from functools import lru_cache
from .cold_leptons import ELECTRON_MASS, MUON_MASS, free_fermion
from .gdfm_uniform import MASS, baryons
from .native import project_root
from .scurto_interface import spherical_phi, surface_sigma, surface_sigma_prime

E2 = 1.43996448


@dataclass(frozen=True)
class SurfaceParameters:
    """Immutable p=3 surface prescription, with explicit model and input identity."""
    model_id: str
    sigma0_MeV_fm2: float
    b: float
    label: str
    source_sha256: str
    p: int = 3

    def __post_init__(self):
        if self.model_id not in {f'SPG_M{i}' for i in range(1,6)}:
            raise ValueError('unknown surface model')
        if not (math.isfinite(self.sigma0_MeV_fm2) and 0 < self.sigma0_MeV_fm2 <= 1e6
                and math.isfinite(self.b) and -16 < self.b <= 1e100 and self.p == 3):
            raise ValueError('surface parameters outside supported p=3 domain')
        if not self.label or len(self.source_sha256)!=64 or any(c not in '0123456789abcdef' for c in self.source_sha256):
            raise ValueError('surface label and SHA256 provenance required')


@lru_cache(maxsize=25)
def _surface_from_content(model, raw):
    values=json.loads(raw)['models'][model]['all_mass']
    return SurfaceParameters(model,float(values['sigma0_MeV_fm2']),float(values['b']),
                             'AME2408_candidate',hashlib.sha256(raw).hexdigest())


def default_surface_parameters(model):
    """Read the current candidate bytes; cache only by exact content identity."""
    path=project_root()/'data/research_generator_20261116/candidate_surface/candidate_surface_fit.json'
    return _surface_from_content(model,path.read_bytes())


def surface_parameters(model):
    """Backward-compatible numeric pair for the current default prescription."""
    p=default_surface_parameters(model)
    return p.sigma0_MeV_fm2,p.b


def _resolve_surface(model,surface):
    p=default_surface_parameters(model) if surface is None else surface
    if not isinstance(p,SurfaceParameters) or p.model_id!=model:
        raise ValueError('surface prescription model does not match bulk model')
    return p


def cell(model, nb, x, *, surface=None):
    """Energy and gradient for x=(volume fraction, radius, cluster density, yp).

    The gas density is constrained by baryon conservation, n_e=f*n_i*y.
    Gradient is in MeV/fm^3 per corresponding coordinate.
    """
    f, radius, ni, y = map(float, x)
    if not all(map(math.isfinite, (nb, f, radius, ni, y))) or not (nb > 0 and 0 < f < 1 and radius > 0 and ni > 0 and 1e-8 <= y <= 1-1e-8):
        raise ValueError('invalid cell parameters')
    ng = (nb-f*ni)/(1-f)
    if ng < 0 and ng >= -64*np.finfo(float).eps*max(nb, ni):
        ng = 0.0
    if not 0 <= ng < ni:
        raise ValueError('pure neutron gas density outside [0, cluster density)')
    surface = _resolve_surface(model,surface)
    sigma0,b = surface.sigma0_MeV_fm2,surface.b
    cluster = baryons(model, ni, y)
    gas = (baryons(model, ng, 0.0) if ng > 0 else
           {'epsilon_baryons': 0.0, 'mu_n': MASS, 'mu_p': MASS})
    ne = f*ni*y
    ee, _, mue = free_fermion(ne, ELECTRON_MASS)
    phi = spherical_phi(f)
    dphi = (1-f**(-2/3))/5
    sigma = surface_sigma(y, sigma0, b)
    dsigma = surface_sigma_prime(y, sigma0, b)
    coul = 2*math.pi*E2*f*phi*radius**2*(ni*y)**2
    surf = 3*f*sigma/radius
    eps = f*cluster['epsilon_baryons']+(1-f)*gas['epsilon_baryons']+ee+coul+surf
    mun, mup, mug = cluster['mu_n'], cluster['mu_p'], gas['mu_n']
    grad = np.array([
        cluster['epsilon_baryons']-gas['epsilon_baryons']+mug*(ng-ni)+mue*ni*y
        + 2*math.pi*E2*radius**2*(ni*y)**2*(phi+f*dphi)+3*sigma/radius,
        (2*coul-surf)/radius,
        f*(mun+(mup-mun)*y-mug+mue*y)+2*coul/ni,
        f*ni*(mup-mun+mue)+2*coul/y+3*f*dsigma/radius,
    ])
    return eps, grad, {'n_gas':ng, 'n_e':ne, 'mu_n_gas':mug,
        'epsilon_coulomb':coul, 'epsilon_surface':surf, 'epsilon_electron':ee,
        'n_proton_cluster':ni*y, 'sigma0':sigma0, 'b':b,
        'mu_e':mue, 'muon_insertion_margin_MeV':MUON_MASS-mue,
        'gas_proton_insertion_margin_MeV':gas['mu_p']-mug+mue
            -4*math.pi*E2*f*phi*radius**2*ni*y/(1-f)}


def solve(model, nb, starts, *, surface=None):
    """Search explicit finite bounds; expose all attempts and unscaled residuals."""
    surface=_resolve_surface(model,surface)
    lower_ni = max(0.06, nb*1.001)
    if not 0 < nb < 0.06 or lower_ni >= 0.25:
        raise ValueError('candidate search supports 0 < n_B < 0.06 fm^-3')
    bounds = [(1e-5, .99999), (2., 40.), (lower_ni, .25), (.005, .45)]
    def unpack(z):
        t, radius, ni, y = z
        return np.array([t*nb/ni, radius, ni, y])
    def objective(z):
        x = unpack(z)
        e, g, _ = cell(model, nb, x, surface=surface)
        t, _, ni, _ = z
        return (e-MASS*nb)/nb, np.array([g[0]*nb/ni, g[1], g[2]-g[0]*x[0]/ni, g[3]])/nb
    attempts=[]
    for start in starts:
        result = minimize(objective, np.asarray(start,float), method='L-BFGS-B', jac=True,
                          bounds=bounds, options={'maxiter':1000,'ftol':1e-14,'gtol':1e-8,'maxls':50})
        x=unpack(result.x); energy, gradient, state=cell(model,nb,x,surface=surface)
        attempts.append({'success':bool(result.success), 'message':str(result.message),
            'coordinates':[float(v) for v in x], 'search_coordinates':[float(v) for v in result.x],
            'energy':float(energy), 'gradient':[float(v) for v in gradient],
            'transformed_gradient':[float(v) for v in objective(result.x)[1]], 'state':state})
    return min(attempts,key=lambda a:a['energy']), attempts


# Reduced variational solver. The legacy four-coordinate solve above is kept
# to reproduce the initial follow-up; new sequences use this analytically
# eliminated-radius formulation and expose finite-domain stationary candidates.
def reduced_cell(model, nb, z, *, surface=None):
    """Return E/A relative to MASS and its envelope gradient in (t,n_i,y).

    t=f*n_i/nb is the cluster baryon fraction. The exact spherical virial
    minimum eliminates R; t=1 is the physical no-free-neutron boundary.
    """
    t, ni, y = map(float, z)
    if not (math.isfinite(nb) and nb > 0 and 0 < t <= 1 and ni > nb and 1e-8 <= y < .5):
        raise ValueError('invalid reduced CLD coordinates')
    f = t*nb/ni
    surface = _resolve_surface(model,surface)
    sigma0,b = surface.sigma0_MeV_fm2,surface.b
    sigma = surface_sigma(y, sigma0, b)
    radius = (3*sigma/(4*math.pi*E2*spherical_phi(f)*(ni*y)**2))**(1/3)
    eps, gradient, state = cell(model, nb, [f, radius, ni, y], surface=surface)
    reduced_gradient = np.array([gradient[0]/ni,
        (gradient[2]-gradient[0]*f/ni)/nb, gradient[3]/nb])
    return (eps-MASS*nb)/nb, reduced_gradient, {
        'coordinates': [f, radius, ni, y], 'energy': eps,
        'surface_parameters':asdict(surface),
        'gradient': gradient.tolist(), 'state': state}


def _reduced_hessian(model, nb, z, *, surface=None):
    z = np.asarray(z,float)
    h = np.array([2e-5,2e-6,2e-6])
    columns=[]
    for j in range(3):
        a=z.copy();b=z.copy();a[j]+=h[j];b[j]-=h[j]
        if j==0 and a[j]>1:
            a[j]=z[j]; b[j]=z[j]-h[j]
        columns.append((reduced_cell(model,nb,a,surface=surface)[1]-reduced_cell(model,nb,b,surface=surface)[1])/(a[j]-b[j]))
    matrix=np.asarray(columns).T
    return (matrix+matrix.T)/2


def solve_reduced(model, nb, starts=None, *, surface=None):
    """Finite multistart search plus analytic-gradient stationarity polishing.

    The old 1e-4 residual limit is not relaxed: this solver requires 1e-5
    projected derivatives in the declared reduced coordinates. Source-table
    agreement does not enter the acceptance decision.
    """
    surface=_resolve_surface(model,surface)
    if not math.isfinite(nb) or not 1e-7 <= nb <= .12:
        raise ValueError('computational search supports 1e-7 <= n_B <= .12 fm^-3')
    lo=np.array([1e-8,max(.025,nb*(1+1e-6)),.0001])
    hi=np.array([1.,.25,.499999])
    if starts is None:
        starts=[[.9,.16,.3],[.5,.15,.14],[.2,.14,.08]]
    bounds=list(zip(lo,hi));attempts=[]
    def objective(z):
        e,g,_=reduced_cell(model,nb,z,surface=surface)
        return e,g
    for initial in starts:
        seed=np.clip(np.asarray(initial,float),lo+1e-10,hi-1e-10)
        result=minimize(objective,seed,method='L-BFGS-B',jac=True,bounds=bounds,
            options={'maxiter':800,'ftol':1e-14,'gtol':1e-9,'maxls':50})
        z=result.x.copy()
        before={'success':bool(result.success),'message':str(result.message),
            'coordinates':z.tolist(),'energy_per_baryon':float(result.fun),
            'gradient':objective(z)[1].tolist()}
        # Physical t=1 boundary is a two-variable stationary problem.
        boundary=z[0]>1-1e-6 and objective(z)[1][0]<0
        if boundary:
            z[0]=1.
            fit=least_squares(lambda q:objective([1.,q[0],q[1]])[1][1:],z[1:],
                bounds=(lo[1:]+1e-12,hi[1:]-1e-12),jac='3-point',diff_step=1e-5,
                xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=120)
            z=np.array([1.,*fit.x])
        else:
            fit=least_squares(lambda q:objective(q)[1],np.clip(z,lo+1e-11,hi-1e-11),
                bounds=(lo+1e-12,hi-1e-12),jac='3-point',diff_step=1e-5,
                xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=120)
            z=fit.x
        energy,gradient,detail=reduced_cell(model,nb,z,surface=surface)
        projected=gradient.copy()
        if z[0]>=1-1e-8: projected[0]=max(gradient[0],0.)
        active=[j for j in (1,2) if min(z[j]-lo[j],hi[j]-z[j])<1e-6*(hi[j]-lo[j])]
        if z[0]-lo[0]<1e-6*(hi[0]-lo[0]): active.append(0)
        violations=[]
        if not fit.success:violations.append('STATIONARITY_SOLVER_FAILED')
        if max(abs(projected))>1e-5:violations.append('STATIONARITY_RESIDUAL')
        if active:violations.append('SEARCH_BOUND_ACTIVE')
        state=detail['state']
        if state['muon_insertion_margin_MeV']<0:violations.append('MUON_INSERTION_UNSTABLE')
        if state['gas_proton_insertion_margin_MeV']<0:violations.append('GAS_PROTON_INSERTION_UNSTABLE')
        # A saddle root must not be promoted merely because its gradient is zero.
        mineig=None
        if not active and max(abs(projected))<1e-4:
            hessian=_reduced_hessian(model,nb,z,surface=surface)
            if z[0]>=1-1e-8:hessian=hessian[1:,1:];scale=np.array([.1,.1])
            else:scale=np.array([1.,.1,.1])
            mineig=float(np.linalg.eigvalsh(hessian*scale[:,None]*scale[None,:])[0])
            if mineig < -1e-6:violations.append('NEGATIVE_LOCAL_CURVATURE')
        mu=state['mu_n_gas']+z[0]*gradient[0]
        f,r,ni,y=detail['coordinates'];a=4*math.pi*r**3*ni/3
        row={**detail,'n_B':float(nb),'search_coordinates':z.tolist(),
            'energy_per_baryon':float(energy),'reduced_gradient':gradient.tolist(),
            'max_projected_gradient':float(max(abs(projected))),
            'physical_boundary':'NEUTRON_GAS_ABSENT' if z[0]>=1-1e-8 else 'INTERIOR_GAS',
            'minimum_scaled_hessian_eigenvalue':mineig,'violations':violations,
            'status':'STATIONARY_CANDIDATE' if not violations else 'HOLD',
            'mu_B_MeV':float(mu),'pressure_MeV_fm3':float(nb*mu-detail['energy']),
            'composition':{'Yn':float(1-z[0]),'Yp':0.,'Ye':float(z[0]*y),'Ymu':0.,
                'Aav':float(a),'Zav':float(a*y),'Yh':float(z[0]),
                'cluster_number_fraction':float(z[0]/a)},
            'initial_search':before,'polish_success':bool(fit.success),
            'polish_message':str(fit.message),'search_bounds':[lo.tolist(),hi.tolist()]}
        attempts.append(row)
    best=min(attempts,key=lambda item:item['energy_per_baryon'])
    return best,attempts



def two_phase_cell(model, nb, x, *, surface=None):
    """Scurto Eqs.21-25 with proton-bearing gas; x=(f,R,ni,yi,yg)."""
    f,r,ni,yi,yg=map(float,x)
    if not (0<f<1 and r>0 and ni>0 and 1e-8<=yi<1 and 0<=yg<1):
        raise ValueError('invalid two-phase cell')
    ng=(nb-f*ni)/(1-f)
    if ng<0 and ng>=-64*np.finfo(float).eps*max(nb,ni):ng=0.
    if not 0<=ng<ni:raise ValueError('invalid constrained gas density')
    cluster=baryons(model,ni,yi)
    gas=baryons(model,ng,yg) if ng>0 else {'epsilon_baryons':0.,'mu_n':MASS,'mu_p':MASS}
    ne=f*ni*yi+(1-f)*ng*yg
    ee,pe,mue=free_fermion(ne,ELECTRON_MASS)
    contrast=ni*yi-ng*yg
    if contrast<=0:raise ValueError('positive cluster proton-density contrast required')
    phi=spherical_phi(f);dphi=(1-f**(-2/3))/5
    surface=_resolve_surface(model,surface)
    sigma0,b=surface.sigma0_MeV_fm2,surface.b
    sigma=surface_sigma(yi,sigma0,b);ds=surface_sigma_prime(yi,sigma0,b)
    coul=2*math.pi*E2*f*phi*r*r*contrast**2;surf=3*f*sigma/r
    mi=cluster['mu_n']+(cluster['mu_p']-cluster['mu_n'])*yi
    mg=gas['mu_n']+(gas['mu_p']-gas['mu_n'])*yg
    gradient=np.array([
        cluster['epsilon_baryons']-gas['epsilon_baryons']+mg*(ng-ni)+mue*ni*(yi-yg)
            +2*math.pi*E2*r*r*(contrast**2*(phi+f*dphi)+2*f*phi*contrast*yg*(ni-ng)/(1-f))+3*sigma/r,
        (2*coul-surf)/r,
        f*(mi-mg+mue*(yi-yg))+2*coul/contrast*(yi+f*yg/(1-f)),
        f*ni*(cluster['mu_p']-cluster['mu_n']+mue)+2*coul*ni/contrast+3*f*ds/r,
        (1-f)*ng*(gas['mu_p']-gas['mu_n']+mue)-2*coul*ng/contrast])
    eps=f*cluster['epsilon_baryons']+(1-f)*gas['epsilon_baryons']+ee+coul+surf
    partial_nb=mg+mue*yg-2*coul*yg/(contrast*(1-f))
    state={'n_gas':ng,'gas_proton_fraction':yg,'n_e':ne,'mu_n_gas':gas['mu_n'],
        'mu_p_gas':gas['mu_p'],'epsilon_coulomb':coul,'epsilon_surface':surf,
        'epsilon_electron':ee,'n_proton_cluster':ni*yi,'sigma0':sigma0,'b':b,'mu_e':mue,
        'muon_insertion_margin_MeV':MUON_MASS-mue,
        'gas_proton_insertion_margin_MeV':gas['mu_p']-gas['mu_n']+mue-2*coul/(contrast*(1-f)),
        'partial_energy_nb_MeV':partial_nb}
    return eps,gradient,state


def reduced_two_phase_cell(model,nb,z, *, surface=None):
    """Reduced variables t,n_i,y_i,v with y_g=v*y_i and analytic optimal R."""
    t,ni,yi,v=map(float,z)
    if not (0<t<=1 and ni>nb>0 and 1e-8<=yi<.5 and 0<=v<1):
        raise ValueError('invalid reduced two-phase coordinates')
    f=t*nb/ni;ng=nb*(1-t)/(1-f);yg=v*yi
    contrast=ni*yi-ng*yg
    surface=_resolve_surface(model,surface)
    sigma0,b=surface.sigma0_MeV_fm2,surface.b
    r=(3*surface_sigma(yi,sigma0,b)/(4*math.pi*E2*spherical_phi(f)*contrast**2))**(1/3)
    eps,g,state=two_phase_cell(model,nb,[f,r,ni,yi,yg],surface=surface)
    grad=np.array([g[0]/ni,(g[2]-g[0]*f/ni)/nb,(g[3]+v*g[4])/nb,yi*g[4]/nb])
    return (eps-MASS*nb)/nb,grad,{'energy':eps,'coordinates':[f,r,ni,yi,yg],
        'surface_parameters':asdict(surface),
        'gradient':g.tolist(),'state':state}


def solve_two_phase(model,nb,starts=None,*,force_proton_search=False,surface=None):
    """Complementary proton-gas solve; retain the restricted KKT branch too."""
    surface=_resolve_surface(model,surface)
    restricted,restricted_attempts=solve_reduced(model,nb,starts,surface=surface)
    if restricted['state']['gas_proton_insertion_margin_MeV']>=0 and not force_proton_search:
        restricted['gas_protons']='KKT_ZERO';restricted['state']['gas_proton_fraction']=0.
        return restricted,restricted_attempts
    lo=np.array([1e-8,max(.025,nb*(1+1e-6)),.0001,0.])
    hi=np.array([1.,.25,.499999,.999999])
    base=np.array(restricted['search_coordinates'])
    def objective(z):return reduced_two_phase_cell(model,nb,z,surface=surface)[:2]
    attempts=[]
    for v0 in (.001,.05,.3):
        seed=np.clip(np.r_[base,v0],lo+1e-10,hi-1e-10)
        result=minimize(objective,seed,method='L-BFGS-B',jac=True,bounds=list(zip(lo,hi)),
            options={'maxiter':1000,'ftol':1e-14,'gtol':1e-9,'maxls':50})
        if result.x[3]<1e-7 and objective(result.x)[1][3]>=0 and restricted['status']=='STATIONARY_CANDIDATE':
            boundary=dict(restricted)
            boundary['gas_protons']='KKT_ZERO_FROM_POSITIVE_START'
            boundary['initial_search']={'success':bool(result.success),'message':str(result.message),
                'coordinates':result.x.tolist(),'energy_per_baryon':float(result.fun),'initial_gas_ratio':v0}
            attempts.append(boundary)
            continue
        fit=least_squares(lambda q:objective(q)[1],np.clip(result.x,lo+1e-10,hi-1e-10),
            bounds=(lo+1e-12,hi-1e-12),jac='3-point',diff_step=1e-5,
            xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=160)
        z=fit.x;e,g,detail=reduced_two_phase_cell(model,nb,z,surface=surface)
        projected=g.copy()
        if z[3]<=1e-8:projected[3]=min(g[3],0.)
        if z[0]>=1-1e-8:projected[0]=max(g[0],0.)
        active=[j for j in (1,2) if min(z[j]-lo[j],hi[j]-z[j])<1e-6*(hi[j]-lo[j])]
        if z[0]-lo[0]<1e-6*(hi[0]-lo[0]):active.append(0)
        if hi[3]-z[3]<1e-6:active.append(3)
        violations=[]
        if not fit.success:violations.append('STATIONARITY_SOLVER_FAILED')
        if max(abs(projected))>1e-5:violations.append('STATIONARITY_RESIDUAL')
        if active:violations.append('SEARCH_BOUND_ACTIVE')
        state=detail['state']
        if state['muon_insertion_margin_MeV']<0:violations.append('MUON_INSERTION_UNSTABLE')
        mineig=None
        if not active and max(abs(projected))<1e-4:
            steps=np.array([2e-6,2e-7,2e-7,min(2e-6,z[3]*.1,(1-z[3])*.1)])
            matrix=[]
            for j,h in enumerate(steps):
                a=z.copy();b=z.copy();a[j]+=h;b[j]-=h
                if j==0 and a[j]>1:a[j]=z[j]
                matrix.append((objective(a)[1]-objective(b)[1])/(a[j]-b[j]))
            hessian=np.asarray(matrix).T;hessian=(hessian+hessian.T)/2
            scale=np.array([1.,.1,.1,1.])
            mineig=float(np.linalg.eigvalsh(hessian*scale[:,None]*scale[None,:])[0])
            if mineig< -1e-6:violations.append('NEGATIVE_LOCAL_CURVATURE')
        mu=state['partial_energy_nb_MeV']+z[0]*g[0]
        f,r,ni,yi,yg=detail['coordinates'];a=4*math.pi*r**3*ni/3
        attempts.append({**detail,'n_B':float(nb),'search_coordinates':z.tolist(),
            'energy_per_baryon':float(e),'reduced_gradient':g.tolist(),
            'max_projected_gradient':float(max(abs(projected))),
            'physical_boundary':'PROTON_BEARING_GAS','gas_protons':'VARIATIONAL_POSITIVE',
            'minimum_scaled_hessian_eigenvalue':mineig,'violations':violations,
            'status':'STATIONARY_CANDIDATE' if not violations else 'HOLD',
            'mu_B_MeV':float(mu),'pressure_MeV_fm3':float(nb*mu-detail['energy']),
            'composition':{'Yn':float((1-z[0])*(1-yg)),'Yp':float((1-z[0])*yg),
                'Ye':float(z[0]*yi+(1-z[0])*yg),'Ymu':0.,'Aav':float(a),'Zav':float(a*yi),
                'Yh':float(z[0]),'cluster_number_fraction':float(z[0]/a)},
            'initial_search':{'success':bool(result.success),'message':str(result.message),
                'coordinates':result.x.tolist(),'energy_per_baryon':float(result.fun)},
            'polish_success':bool(fit.success),'polish_message':str(fit.message),
            'search_bounds':[lo.tolist(),hi.tolist()]})
    candidates=attempts+[restricted]
    best=min(candidates,key=lambda row:row['energy_per_baryon'])
    return best,restricted_attempts+attempts
