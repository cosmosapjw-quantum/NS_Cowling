"""Generate the five declared cold core models on caller-supplied density grids.

The computational states depend on published functionals/energy samples, not
CompOSE target thermodynamics. Source comparison and density support are separate
from implementation checks. Units: fm^-3 and MeV; pressure/energy MeV fm^-3.
"""
from __future__ import annotations
import numpy as np
from . import apr_uniform, apr_gibbs, bbb_uniform, gm1_uniform, skyrme_uniform
from .eos_data import EOSData

MODEL_CONTRACTS = {
    'APR_APR': {'default_interval_fm3':[.085,1.34], 'species':'npemu',
        'representation':'Published21parameterAPRfunctional; bulk two-charge Gibbs',
        'paper_implementation':'AUTHOR_FUNCTIONAL_WITH_EXPLICIT_BULK_GIBBS_ASSUMPTION',
        'source':apr_uniform.SOURCE},
    'BBB_BHF2': {'default_interval_fm3':[.08,1.5], 'species':'npemu',
        'representation':'PublishedParis+UrbanaSNM/PNMenergies; parabolic isospin; PCHIP energy derivative',
        'paper_implementation':'PUBLISHED_ENERGY_TABLE_RECONSTRUCTION; original BHF kernel is not re-solved',
        'microscopic_support_max_fm3':.8,'source':bbb_uniform.SOURCE},
    'GM_GM1': {'default_interval_fm3':[.08,1.0], 'species':'npe',
        'representation':'Original rounded nonlinear RMF couplings; common nucleon rest energy939MeV',
        'paper_implementation':'PUBLISHED_PARAMETER_FUNCTIONAL','source':gm1_uniform.SOURCE},
    'RG_KDE0v1': {'default_interval_fm3':[.08,1.8], 'species':'npemu',
        'representation':'Published full asymmetric Skyrme functional',
        'paper_implementation':'PUBLISHED_PARAMETER_FUNCTIONAL','source':skyrme_uniform.parameter_provenance('RG_KDE0v1')},
    'RG_SLY4': {'default_interval_fm3':[.08,1.4], 'species':'npemu',
        'representation':'Published full asymmetric Skyrme functional',
        'paper_implementation':'PUBLISHED_PARAMETER_FUNCTIONAL','source':skyrme_uniform.parameter_provenance('RG_SLY4')},
}

def state_solver(model_id):
    """Return a state evaluator and its phase context without consulting target tables."""
    if model_id not in MODEL_CONTRACTS: raise ValueError('unknown core model')
    context={}
    if model_id=='APR_APR':
        context=apr_gibbs.gibbs_coexistence()
        evaluate=lambda n:apr_gibbs.gibbs_equilibrium_state(float(n),context)
    elif model_id=='BBB_BHF2':evaluate=lambda n:bbb_uniform.beta_state(float(n))
    elif model_id=='GM_GM1':evaluate=lambda n:gm1_uniform.beta_state(float(n),include_muons=False)
    else:
        def evaluate(n):
            raw=skyrme_uniform.beta_equilibrium(model_id,float(n))
            return {**raw,'n_b':raw['n_b_fm3'],'epsilon':raw['energy_density_mev_fm3'],
                'pressure':raw['pressure_mev_fm3'],'mu_n':raw['mu_n_mev'],
                'mu_p':raw['mu_p_mev'],'mu_e':raw['mu_e_mev']}
    return evaluate,context

def generate_core(model_id, densities):
    """Return an EOSData table, solved states and exact declared model/phase context."""
    grid=np.asarray(densities,dtype=float)
    if (grid.ndim!=1 or grid.size<3 or np.any(~np.isfinite(grid)) or
        np.any(grid<=0) or np.any(np.diff(grid)<=0)):
        raise ValueError('at least three finite positive increasing densities required')
    solve,context=state_solver(model_id)
    if model_id=='APR_APR':
        edges=[context[k]['n_b'] for k in ('left','right')]
        grid=np.unique(np.concatenate((grid,[n for n in edges if grid[0]<n<grid[-1]])))
        mn,mp=apr_uniform.M_N,apr_uniform.M_P
    elif model_id.startswith('RG_'):mn,mp=skyrme_uniform.MN,skyrme_uniform.MP
    else:mn=mp=939.0
    rows=[solve(n) for n in grid]
    col=lambda k:np.asarray([r[k] for r in rows],dtype=float)
    phase=np.full(grid.size,6,dtype=int)
    if model_id=='APR_APR':
        phase=np.array([6 if r['branch']=='low' else 7 if r['branch']=='high' else 8 for r in rows])
    zero=np.zeros(grid.size)
    table=EOSData(nb=grid,mn=mn,mp=mp,Q1=col('pressure')/grid,Q2=zero.copy(),
        Q3=col('mu_n')/mn-1,Q4=-col('mu_e')/mn,Q5=zero.copy(),
        Q6=col('epsilon')/(grid*mn)-1,Q7=col('epsilon')/(grid*mn)-1,
        q_extra=zero.copy(),Iphase=phase,Yn=col('Y_n'),Yp=col('Y_p'),
        Ye=col('Y_e'),Ymu=col('Y_mu'),Aav=np.ones(grid.size),Zav=zero.copy(),
        Yh=zero.copy(),composition_file_present=True,composition_available=True,beta_equilibrated=True)
    closure=max(float(np.max(abs(table.Yn+table.Yp-1))),float(np.max(abs(table.Yp-table.Ye-table.Ymu))))
    euler=float(np.max(abs(col('epsilon')+col('pressure')-grid*col('mu_n'))/(grid*mn)))
    beta=float(np.max(abs(col('mu_n')-col('mu_p')-col('mu_e'))))
    if closure>1e-10 or euler>1e-10 or beta>1e-7:
        raise RuntimeError('generated core violates its charge/baryon/Euler/beta contract')
    diagnostics={'fraction_closure_max_abs':closure,'euler_scaled_max_abs':euler,
        'beta_max_abs_MeV':beta,'composition_phase_codes':{'6':'uniformlow/default','7':'APRhigh','8':'APRbulkGibbs'},
        'acceptance':'Independent functional/constraints; CompOSE differences have no unconditional rejection threshold',
        'requested_interval_fm3':[float(grid[0]),float(grid[-1])],
        'physical_support':'Determined separately by density-resolved derivatives and model assumptions'}
    return table,rows,{'model_contract':MODEL_CONTRACTS[model_id], 'phase_context':context,'implementation_checks':diagnostics}
