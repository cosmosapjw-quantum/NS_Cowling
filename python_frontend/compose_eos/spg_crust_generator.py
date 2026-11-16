"""Composition-aware standalone SPG spherical CLD sequences.

The candidate surface fit is explicit. Independent variational checks govern
implementation acceptance; external CompOSE differences are diagnostics.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import brentq, least_squares
from .candidate_gdfm_cld import MASS, reduced_cell, solve_two_phase
from .gdfm_uniform import beta_state
from .eos_data import EOSData


def neutron_drip(model):
    """Locate the constrained t=1 endpoint where the neutron insertion cost vanishes."""
    def vacuum(n):
        fit=least_squares(lambda q:reduced_cell(model,n,[1.,q[0],q[1]])[1][1:], [.16,.33],
            bounds=([.05,.01],[.25,.499]),jac='3-point',diff_step=1e-5,
            xtol=1e-12,ftol=1e-12,gtol=1e-12)
        e,g,detail=reduced_cell(model,n,[1.,*fit.x])
        if not fit.success or max(abs(g[1:]))>1e-5:raise RuntimeError('vacuum endpoint stationarity failed')
        return fit.x,g,detail
    n=brentq(lambda x:vacuum(x)[1][0],1e-6,.002,xtol=2e-14)
    z,g,detail=vacuum(n)
    return {'n_B':float(n),'mu_B_MeV':float(MASS+g[0]),
        'pressure_MeV_fm3':float(n*(MASS+g[0])-detail['energy']),
        'stationarity_gradient':g.tolist(),'coordinates':detail['coordinates'],
        'energy':detail['energy'],'definition':'neutron KKT boundary in the declared common938.9MeV rest-energy convention'}


def homogeneous_endpoint(model, bracket):
    """Two locally neutral phases of the same bulk: identify the CLD support end.

    Only the cluster-side endpoint is exported as crust. The uniform phase is
    used to determine a thermodynamic endpoint, not to construct a matched EOS.
    """
    def state(n):return solve_two_phase(model,float(n))[0]
    def delta(n):return state(n)['energy']-beta_state(model,n)['epsilon']
    cross=brentq(delta,*bracket,xtol=2e-10)
    def residual(v):
        c=state(v[0]);u=beta_state(model,float(v[1]))
        return [c['mu_B_MeV']-u['mu_n'],10*(c['pressure_MeV_fm3']-u['pressure'])]
    fit=least_squares(residual,[cross*.995,cross*1.005],
        bounds=([bracket[0]*.98,bracket[0]*.98],[bracket[1]*1.02,bracket[1]*1.02]),
        diff_step=2e-5,xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=80)
    low=state(fit.x[0]);high=beta_state(model,float(fit.x[1]))
    errors=residual(fit.x)
    if not fit.success or max(abs(v) for v in errors)>1e-6 or fit.x[1]<=fit.x[0]:
        raise RuntimeError('CLD/uniform endpoint coexistence did not converge')
    if low['status']!='STATIONARY_CANDIDATE':raise RuntimeError('CLD endpoint has unresolved variational conditions')
    pressure=(low['pressure_MeV_fm3']+high['pressure'])/2
    chemical=(low['mu_B_MeV']+high['mu_n'])/2
    return {'n_cluster_fm3':float(fit.x[0]),'n_uniform_fm3':float(fit.x[1]),
        'pressure_MeV_fm3':float(pressure),'mu_B_MeV':float(chemical),
        'fixed_density_energy_crossing':float(cross),'residuals':errors,
        'cluster_state':low,'uniform_state':high,
        'construction':'local-neutral same-functional common tangent; no core model splicing'}


def sequence(model,densities):
    """Use continuation plus three fixed seeds; preserve every computed attempt."""
    densities=np.asarray(densities,float)
    if densities.ndim!=1 or len(densities)<3 or np.any(np.diff(densities)<=0):
        raise ValueError('strictly increasing density sequence of at least3points required')
    rows=[];previous=None
    for n in densities:
        starts=[[.9,.16,.3],[.5,.15,.14],[.2,.14,.08]]
        if previous is not None:starts.insert(0,previous)
        best,attempts=solve_two_phase(model,float(n),starts)
        previous=best['search_coordinates'][:3]
        rows.append({'n_B':float(n),'best':best,'attempts':attempts})
    return rows


def to_eos_data(rows):
    """Export full free/cluster composition using the native Yh baryon convention."""
    if any(row['best']['status']!='STATIONARY_CANDIDATE' for row in rows):
        raise ValueError('unresolved variational rows cannot be exported as stationary EOS')
    n=np.array([row['n_B'] for row in rows]);states=[row['best'] for row in rows]
    def col(key):return np.array([state[key] for state in states])
    def comp(key):return np.array([state['composition'][key] for state in states])
    e=col('energy')/(n*MASS)-1;zero=np.zeros(len(n))
    return EOSData(nb=n,mn=MASS,mp=MASS,Q1=col('pressure_MeV_fm3')/n,
        Q2=zero.copy(),Q3=col('mu_B_MeV')/MASS-1,
        Q4=-np.array([s['state']['mu_e'] for s in states])/MASS,Q5=zero.copy(),
        Q6=e,Q7=e.copy(),q_extra=zero.copy(),Iphase=np.full(len(n),7,dtype=int),
        Yn=comp('Yn'),Yp=comp('Yp'),Ye=comp('Ye'),Ymu=comp('Ymu'),
        Aav=comp('Aav'),Zav=comp('Zav'),Yh=comp('Yh'),
        composition_file_present=True,composition_available=True,beta_equilibrated=True)
