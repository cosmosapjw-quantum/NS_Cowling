"""Cold barotropic analytic EOS, independent SI oracles (not a stellar calibration)."""
import os
from pathlib import Path
import pytest
from tools.runtime import CaseRunner

@pytest.fixture(scope='module')
def poly(tmp_path_factory):
    runner=CaseRunner(Path(__file__).resolve().parents[1],tmp_path_factory.mktemp('poly'),
                      workspace=os.environ.get('NS_WORKSPACE','solutions'))
    art=runner.build_fpm()
    return lambda pars={}:runner.run_artifact(art,'poly_state',pars)

def test_polytrope_thermodynamics(poly):
    r=poly()
    assert r['status']=='ok'
    assert r['energy_j_m3']==pytest.approx(1.51e34,rel=2e-14)
    assert r['mu_j']==pytest.approx(1.52e-10,rel=2e-14)
    assert r['chi']==pytest.approx(2e32/(1.52e34),rel=2e-14)

def test_pressure_inverse_first_success(poly):
    r=poly({'query':'pressure','value':1e32})
    assert r['status']=='ok'
    assert r['pressure_pa']==pytest.approx(1e32,rel=2e-14)

def test_failure_does_not_export_old_state(poly):
    assert poly()['status']=='ok'
    r=poly({'n_b_m3':-1.})
    assert r['status']=='out_of_domain'
    assert r['result_valid'] is False and 'energy_j_m3' not in r

import math
import json
from decimal import Decimal, localcontext

BASE=dict(n_anchor_m3=1e44,p_anchor_pa=1e32,gamma=2.0,c_j=1.5e-10,n_min_m3=1e40,n_max_m3=1e46)

def analytic(n,g=2.,c=1.5e-10,na=1e44,pa=1e32):
    p=pa*(n/na)**g
    e=c*n+p/(g-1)
    mu=c+(g/(g-1))*(p/n)
    return p,e,mu,g*p/(e+p)

def state_close(r,n,g=2.,c=1.5e-10,na=1e44,pa=1e32,rtol=5e-13):
    assert r['status']=='ok',r
    for k,v in zip(('pressure_pa','energy_j_m3','mu_j','chi'),analytic(n,g,c,na,pa)):
        assert r[k]==pytest.approx(v,rel=rtol)
    assert r['n_b_m3']==pytest.approx(n,rel=rtol)

@pytest.mark.parametrize('n',[1e40,1e42,1e44,1e45,1e46])
def test_gamma_two_decimal_oracle(poly,n):
    r=poly({'n_b_m3':n})
    with localcontext() as ctx:
        ctx.prec=60
        d=Decimal(str(n));p=Decimal('1e32')*(d/Decimal('1e44'))**2
        energy=Decimal('1.5e-10')*d+p
        assert r['energy_j_m3']==pytest.approx(float(energy),rel=5e-14)
        assert r['pressure_pa']==pytest.approx(float(p),rel=5e-14)
    state_close(r,n)

@pytest.mark.parametrize('g',[1.2,1.5,5./3.,2.,2.5,3.])
def test_general_gamma_on_causal_interval(poly,g):
    n=1.7e44
    r=poly({'n_b_m3':n,'gamma':g,'n_min_m3':5e43,'n_max_m3':2e44})
    state_close(r,n,g)

@pytest.mark.parametrize('g',[1.5,2.,2.5])
@pytest.mark.parametrize('n',[5e43,1e44,2e44])
def test_analytic_inverses_with_endpoint_binding(poly,g,n):
    cfg={'gamma':g,'n_min_m3':5e43,'n_max_m3':2e44}
    s=poly(dict(cfg,n_b_m3=n));assert s['status']=='ok'
    for by,key in [('pressure','pressure_pa'),('mu','mu_j')]:
        r=poly(dict(cfg,query=by,value=s[key]))
        state_close(r,n,g,rtol=2e-12)
        if by=='mu':
            expected=n*n/(g*s['pressure_pa'])
            assert r['dn_dmu_m3_per_j']==pytest.approx(expected,rel=3e-12)

@pytest.mark.parametrize('pars',[
    {'gamma':1.},{'gamma':0.5},{'p_anchor_pa':0.},{'p_anchor_pa':-1.},
    {'n_anchor_m3':0.},{'c_j':0.},{'c_j':-1.}, {'n_min_m3':0.},
    {'n_min_m3':1e46},{'n_max_m3':1e40},{'n_min_m3':2e44},
    {'constructor':'k','k_si':0.},{'constructor':'k','k_si':-1.},
])
def test_invalid_model_parameters_rejected(poly,pars):
    r=poly(pars)
    assert r['status']=='invalid_input'
    assert r['result_valid'] is False and r['model_ready'] is False
    assert 'energy_j_m3' not in r

@pytest.mark.parametrize('special',['nan_n','inf_pressure','nan_mu','nan_gamma','inf_c'])
def test_nonfinite_direct_fortran_input(poly,special):
    r=poly({'special_input':special})
    assert r['status']=='nonfinite_input'
    assert r['result_valid'] is False and r['state_poisoned'] is True
    assert 'energy_j_m3' not in r

@pytest.mark.parametrize('n',[0.,-1.,math.nextafter(1e40,-math.inf),math.nextafter(1e46,math.inf)])
def test_density_domain_is_closed_and_not_extrapolated(poly,n):
    r=poly({'n_b_m3':n})
    assert r['status']=='out_of_domain' and not r['result_valid']
    assert 'pressure_pa' not in r

@pytest.mark.parametrize('query,key',[('pressure','p_min_pa'),('mu','mu_min_j')])
def test_inverse_domain_one_ulp_is_not_silently_clipped(poly,query,key):
    s=poly(); value=math.nextafter(s[key],-math.inf)
    r=poly({'query':query,'value':value})
    assert r['status']=='out_of_domain' and not r['result_valid']

@pytest.mark.parametrize('query,key',[('pressure','p_max_pa'),('mu','mu_max_j')])
def test_inverse_above_upper_bound_rejected(poly,query,key):
    s=poly(); r=poly({'query':query,'value':math.nextafter(s[key],math.inf)})
    assert r['status']=='out_of_domain'

def test_chemical_potential_c_is_not_vacuum_extension(poly):
    r=poly({'query':'mu','value':BASE['c_j']})
    assert r['status']=='out_of_domain'

@pytest.mark.parametrize('side',[0,1,-1])
def test_known_one_sided_labels_on_smooth_branch_agree(poly,side):
    r=poly({'query':'pressure','value':1e32,'side':side})
    state_close(r,1e44)
    assert r['side']==side

def test_unknown_side_is_rejected(poly):
    r=poly({'query':'pressure','value':1e32,'side':42})
    assert r['status']=='invalid_input' and not r['result_valid']

def test_complete_interval_admission_not_just_one_query(poly):
    # Query itself is causal but upper end is not. The entire supplied branch is rejected.
    r=poly({'gamma':3.})
    assert r['status']=='out_of_domain' and r['category']=='physics'
    assert r['model_ready'] is False
    assert 'causal' in r['detail']
    s=poly({'gamma':3.,'n_max_m3':2e44});assert s['status']=='ok'

def test_finite_k_constructor_agrees_with_anchor(poly):
    r=poly({'constructor':'k','k_si':1e-56,'n_b_m3':1.4e44})
    state_close(r,1.4e44,rtol=3e-13)

def test_anchor_avoids_unrepresentable_n_power(poly):
    p={'gamma':3.,'c_j':1e32,'n_anchor_m3':1e120,'p_anchor_pa':1e150,
       'n_min_m3':5e119,'n_max_m3':2e120,'n_b_m3':1.5e120}
    r=poly(p);state_close(r,1.5e120,3.,1e32,1e120,1e150,rtol=5e-13)
    # n**Gamma overflows real64, but the EOS values and K=1e-210 are finite.
    s=poly(dict(p,constructor='k',k_si=1e-210));state_close(s,1.5e120,3.,1e32,1e120,1e150,rtol=1e-12)

@pytest.mark.parametrize('pars',[
    {'gamma':1.5,'n_anchor_m3':1.,'p_anchor_pa':1e308,'c_j':1.,'n_min_m3':0.5,'n_max_m3':2.},
    {'gamma':2.,'n_anchor_m3':1.,'p_anchor_pa':1.,'c_j':1.,'n_min_m3':1e-200,'n_max_m3':2.},
])
def test_true_representation_range_failure_is_not_zero(poly,pars):
    r=poly(pars)
    assert r['status']=='range_error' and not r['result_valid']
    assert r['state_poisoned'] is True

def test_unresolvable_mu_interval_rejected(poly):
    r=poly({'n_anchor_m3':1.,'p_anchor_pa':1e-20,'c_j':1.,'gamma':2.,
            'n_min_m3':0.5,'n_max_m3':2.,'n_b_m3':1.})
    assert r['status']=='stagnation' and not r['model_ready']

def test_energy_offset_changes_energy_and_mu_not_pressure(poly):
    n=1.3e44;dc=2e-11
    a=poly({'n_b_m3':n});b=poly({'n_b_m3':n,'c_j':BASE['c_j']+dc})
    assert b['pressure_pa']==a['pressure_pa']
    assert (b['energy_j_m3']-a['energy_j_m3'])/n==pytest.approx(dc,rel=2e-13)
    assert b['mu_j']-a['mu_j']==pytest.approx(dc,rel=2e-13)
    assert b['chi']<a['chi']

def test_independent_energy_pressure_derivatives(poly):
    cfg={'gamma':1.6,'n_min_m3':5e43,'n_max_m3':2e44}
    n=1.1e44;s=poly(dict(cfg,n_b_m3=n));errs=[]
    for h in [0.04,0.02,0.01]:
        a=poly(dict(cfg,n_b_m3=n*(1-h)));b=poly(dict(cfg,n_b_m3=n*(1+h)))
        de=(b['energy_j_m3']-a['energy_j_m3'])/(2*n*h)
        dp=(b['pressure_pa']-a['pressure_pa'])/(2*n*h)
        errs.append(abs(de/s['mu_j']-1))
        assert dp/de==pytest.approx(s['chi'],rel=7e-5)
    assert errs[-1]<errs[0]/12
    assert errs[-1]<1e-7

def test_first_law_integral_independent_pressure_path(poly):
    # Simpson on dimensionless u=n/n_anchor; no Python stellar/EOS solver.
    cfg={'gamma':2.4,'n_min_m3':5e43,'n_max_m3':2e44};na=BASE['n_anchor_m3']
    states=[poly(dict(cfg,n_b_m3=na*(1+i/64))) for i in range(65)]
    values=[s['pressure_pa']/(na*(1+i/64)**2) for i,s in enumerate(states)]
    target=states[-1]['energy_j_m3']/states[-1]['n_b_m3']-states[0]['energy_j_m3']/na
    def simp(vals,h):return h*(vals[0]+vals[-1]+4*sum(vals[1:-1:2])+2*sum(vals[2:-1:2]))/3
    e32=abs(simp(values[::2],1/32)-target);e64=abs(simp(values,1/64)-target)
    assert e64<e32/10
    assert e64/abs(target)<3e-10

def test_independent_pressure_and_density_mu_derivatives(poly):
    cfg={'n_min_m3':5e43,'n_max_m3':2e44,'gamma':1.7}
    s=poly(dict(cfg,n_b_m3=1.2e44));mu=s['mu_j'];h=(mu-BASE['c_j'])*1e-4
    a=poly(dict(cfg,query='mu',value=mu-h));b=poly(dict(cfg,query='mu',value=mu+h))
    c=poly(dict(cfg,query='mu',value=mu))
    assert (b['pressure_pa']-a['pressure_pa'])/(2*h)==pytest.approx(s['n_b_m3'],rel=2e-8)
    assert (b['n_b_m3']-a['n_b_m3'])/(2*h)==pytest.approx(c['dn_dmu_m3_per_j'],rel=2e-8)

def test_failed_call_clears_state_and_failed_init_clears_model(poly):
    r=poly({'action':'lifecycle'})
    assert r['status']=='out_of_domain'
    assert not r['result_valid'] and r['state_poisoned']
    assert r['failed_reinit_cleared'] and r['uninitialized_evaluation_rejected']

def test_interleaved_models_do_not_share_mutable_state(poly):
    r=poly({'action':'interleaved'})
    assert r['status']=='ok'
    assert r['independent_models'] and r['first_model_unchanged']

def test_run_contract_names_eos_and_not_tov(poly):
    r=poly();receipt=json.loads(Path(r['_evidence']['run_receipt']).read_text())
    assert receipt['contract']['units']=='SI'
    assert receipt['contract']['model']=='cold_polytrope'
    assert receipt['contract']['eos']=='analytic_anchor_polytrope_v1'
    assert receipt['contract']['closure']=='barotropic_equilibrium_derivative'
    assert receipt['contract']['surface']=='not_applicable'


def test_failed_reinit_does_not_emit_stale_domain(poly):
    r=poly({'action':'lifecycle'})
    assert r['model_ready'] is False
    assert 'p_min_pa' not in r and 'mu_max_j' not in r

@pytest.mark.parametrize('pars',[{'query':'bad'},{'action':'bad'},{'constructor':'bad'},{'special_input':'bad'}])
def test_unknown_request_options_fail_without_state(poly,pars):
    r=poly(pars)
    assert r['status']=='invalid_input' and not r['result_valid']
    assert 'energy_j_m3' not in r

def test_mu_derivative_range_failure_clears_a_previously_valid_state(poly):
    cfg=dict(n_anchor_m3=1e200,p_anchor_pa=1e-8,gamma=2.,c_j=1e-200,
             n_min_m3=5e199,n_max_m3=2e200,n_b_m3=1e200)
    s=poly(cfg);assert s['status']=='ok'
    r=poly(dict(cfg,query='mu',value=s['mu_j']))
    assert r['status']=='range_error'
    assert not r['result_valid'] and r['state_poisoned']
    assert 'dn_dmu_m3_per_j' not in r and 'energy_j_m3' not in r

def test_gibbs_relation_across_the_interval(poly):
    for n in (1e42,1e44,1e46):
        s=poly({'n_b_m3':n})
        assert s['mu_j']==pytest.approx((s['energy_j_m3']+s['pressure_pa'])/s['n_b_m3'],rel=3e-14)

def test_interleaved_diagnostic_does_not_serialize_an_uncomputed_mu_derivative(poly):
    r=poly({'action':'interleaved','query':'mu'})
    assert r['status']=='ok' and r['independent_models']
    assert 'dn_dmu_m3_per_j' not in r
