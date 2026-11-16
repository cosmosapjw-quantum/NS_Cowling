"""Candidate spherical GDFM CLD cell; not source-calibrated.

Pure neutron gas, uniform electrons, no muons, sharp interface. The fitted
surface parameters are the separate AME2408 vacuum-droplet experiment. This
module does not supply an outer crust or certify the global minimum.
"""
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from .cold_leptons import ELECTRON_MASS, MUON_MASS, free_fermion
from .gdfm_uniform import MASS, baryons
from .native import project_root
from .scurto_interface import spherical_phi, surface_sigma, surface_sigma_prime

E2 = 1.43996448


def surface_parameters(model):
    path = project_root() / 'data/research_generator_20261116/candidate_surface/candidate_surface_fit.json'
    values = json.loads(path.read_text())['models'][model]['all_mass']
    return values['sigma0_MeV_fm2'], values['b']


def cell(model, nb, x):
    """Energy and gradient for x=(volume fraction, radius, cluster density, yp).

    The gas density is constrained by baryon conservation, n_e=f*n_i*y.
    Gradient is in MeV/fm^3 per corresponding coordinate.
    """
    f, radius, ni, y = map(float, x)
    if not all(map(math.isfinite, (nb, f, radius, ni, y))) or not (nb > 0 and 0 < f < 1 and radius > 0 and ni > 0 and 1e-8 <= y <= 1-1e-8):
        raise ValueError('invalid cell parameters')
    ng = (nb-f*ni)/(1-f)
    if not 0 < ng < ni:
        raise ValueError('pure neutron gas density outside (0, cluster density)')
    sigma0, b = surface_parameters(model)
    cluster = baryons(model, ni, y)
    gas = baryons(model, ng, 0.0)
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


def solve(model, nb, starts):
    """Search explicit finite bounds; expose all attempts and unscaled residuals."""
    lower_ni = max(0.06, nb*1.001)
    if not 0 < nb < 0.06 or lower_ni >= 0.25:
        raise ValueError('candidate search supports 0 < n_B < 0.06 fm^-3')
    bounds = [(1e-5, .99999), (2., 40.), (lower_ni, .25), (.005, .45)]
    def unpack(z):
        t, radius, ni, y = z
        return np.array([t*nb/ni, radius, ni, y])
    def objective(z):
        x = unpack(z)
        e, g, _ = cell(model, nb, x)
        t, _, ni, _ = z
        return (e-MASS*nb)/nb, np.array([g[0]*nb/ni, g[1], g[2]-g[0]*x[0]/ni, g[3]])/nb
    attempts=[]
    for start in starts:
        result = minimize(objective, np.asarray(start,float), method='L-BFGS-B', jac=True,
                          bounds=bounds, options={'maxiter':1000,'ftol':1e-14,'gtol':1e-8,'maxls':50})
        x=unpack(result.x); energy, gradient, state=cell(model,nb,x)
        attempts.append({'success':bool(result.success), 'message':str(result.message),
            'coordinates':[float(v) for v in x], 'search_coordinates':[float(v) for v in result.x],
            'energy':float(energy), 'gradient':[float(v) for v in gradient],
            'transformed_gradient':[float(v) for v in objective(result.x)[1]], 'state':state})
    return min(attempts,key=lambda a:a['energy']), attempts
