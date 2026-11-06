"""SI TOV and explicit positive-pressure surface contract."""
import math,os
from pathlib import Path
import pytest
from tools.runtime import CaseRunner

@pytest.fixture(scope='module')
def star(tmp_path_factory):
    runner=CaseRunner(Path(__file__).resolve().parents[1],tmp_path_factory.mktemp('tov'),
                      workspace=os.environ.get('NS_WORKSPACE','solutions'))
    art=runner.build_fpm()
    return lambda pars={},case='tov_baseline':runner.run_artifact(art,case,pars)

def test_tov_invariants_and_lapse(star):
    r=star()
    assert r['status']=='ok',r
    assert r['mass_non_decreasing'] and r['pressure_non_increasing']
    assert r['surface_event_closed'] and r['lapse_contract_closed']

def test_radial_si_rhs_first(star):
    r=star({},'tov_rhs')
    assert r['status']=='ok',r
    assert r['dm_dr_kg_m']>0 and r['dp_dr_pa_m']<0

G=6.67430e-11;C=299792458.

def rows(r):
    return [{k:r[f'{tag}_{i:05d}'] for k,tag in [('r','r'),('m','m'),('p','p'),('epsilon','eps'),
            ('n','n'),('mu','mu'),('chi','chi'),('nu','nu'),('lambda','lambda'),('dnu','dnu'),('q','q')]}
            for i in range(1,r['profile_n']+1)]

@pytest.mark.parametrize('radius,mass,p',[(1e3,1e27,1e34),(1e4,1e30,1e30),(2.,1e19,1e32)])
def test_si_rhs_independent_component_values(star,radius,mass,p):
    r=star({'radius_m':radius,'mass_kg':mass,'pressure_pa':p},'tov_rhs')
    assert r['status']=='ok'
    n=1e44*(p/1e32)**0.5; eps=1.5e-10*n+p
    grad=G*(mass+4*math.pi*radius**3*p/C**2)/(C**2*radius**2*(1-2*G*mass/(radius*C**2)))
    assert r['dm_dr_kg_m']==pytest.approx(4*math.pi*radius**2*eps/C**2,rel=2e-13)
    assert r['dp_dr_pa_m']==pytest.approx(-(eps+p)*grad,rel=2e-13)
    assert r['dnu_dr_per_m']==pytest.approx(grad,rel=2e-13)

@pytest.mark.parametrize('pars',[
    {'radius_m':0.},{'radius_m':-1.},{'mass_kg':-1.},{'pressure_pa':0.},
    {'radius_m':1e3,'mass_kg':1e32},{'special_input':'nan_r'},
    {'special_input':'fail_callback'},{'special_input':'wrong_pressure'},
    {'special_input':'invalid_state'},{'special_input':'nan_energy'},
])
def test_bad_radial_rhs_rejected_without_numeric_output(star,pars):
    r=star(pars,'tov_rhs');assert r['status']!='ok'
    assert 'dm_dr_kg_m' not in r

@pytest.mark.parametrize('pars',[
    {'pc_pa':0.},{'pc_pa':-1.},{'pc_pa':1e40},{'p_cut_pa':0.},{'p_cut_pa':1e10},
    {'p_cut_pa':1e35},{'r0_m':0.},{'r0_m':1e-9},{'r0_m':1e5},{'r_limit_m':100.},
    {'rtol':-1.},{'atol_r_m':0.},{'atol_mass_kg':0.},{'atol_nu':0.},
    {'dq_initial':0.},{'dq_min':1.,'dq_max':0.1},{'max_steps':1},
    {'special_input':'nan_pc'},{'special_input':'sharp_event'},
    {'surface_policy':'vacuum'},{'eos_kind':'unknown'},
    {'use_stops':True,'stop1_pa':1e20}, {'use_stops':True,'stop2_pa':1e32},
])
def test_tov_failure_never_exports_stale_background(star,pars):
    r=star(pars);assert r['status']!='ok',r
    assert r['result_valid'] is False and 'radius_m' not in r and 'mass_kg' not in r

def test_pressure_surface_and_profile_values(star):
    r=star({'emit_profile':True,'use_stops':True})
    assert r['status']=='ok',r
    a=rows(r);assert all(a[i+1]['r']>a[i]['r'] for i in range(len(a)-1))
    assert a[-1]['p']==r['p_cut_pa'] and r['surface_residual_pa']==0.
    assert r['stops_hit']==2 and r['coordinate']=='log_pressure'
    assert not r['exact_vacuum_surface']
    for p in [1e31,1e28]:assert min(abs(x['p']/p-1) for x in a)<4e-14
    for s in a[::max(1,len(a)//30)]:
        n=1e44*math.sqrt(s['p']/1e32)
        assert s['epsilon']==pytest.approx(1.5e-10*n+s['p'],rel=2e-13)
        assert math.exp(-2*s['lambda'])==pytest.approx(1-2*G*s['m']/(s['r']*C**2),rel=3e-14)
        assert s['mu']*math.exp(s['nu'])==pytest.approx(a[0]['mu']*math.exp(a[0]['nu']),rel=2e-8)

def test_reinitialized_failed_profile_is_invalid(star):
    r=star({'special_input':'lifecycle'});assert r['status']!='ok'
    assert not r['result_valid'] and 'radius_m' not in r

def test_global_lapse_offset_removed(star):
    a=star();b=star({'nu_start':17.})
    assert b['radius_m']==pytest.approx(a['radius_m'],rel=2e-10)
    assert b['mass_kg']==pytest.approx(a['mass_kg'],rel=2e-10)
    assert b['nu_surface']==pytest.approx(a['nu_surface'],abs=8e-15)
    assert b['lapse_shift']-a['lapse_shift']==pytest.approx(-17.,abs=2e-14)

def uniform_radius(pc,p_cut,eps):
    s=(eps+pc)/(eps+3*pc)
    u=s*(eps+3*p_cut)/(eps+p_cut)
    return math.sqrt(3*C**4*(1-u*u)/(8*math.pi*G*eps))

@pytest.mark.parametrize('pc',[1e33,1e34,5e34])
def test_uniform_density_exact_schwarzschild_interior(star,pc):
    eps=1e35;cut=1e25
    r=star({'eos_kind':'uniform_benchmark','pc_pa':pc,'p_cut_pa':cut,'rtol':1e-10,'emit_profile':True})
    assert r['status']=='ok',r
    R=uniform_radius(pc,cut,eps);M=4*math.pi*eps*R**3/(3*C**2)
    assert r['radius_m']==pytest.approx(R,rel=3e-8)
    assert r['mass_kg']==pytest.approx(M,rel=8e-8)
    assert r['cowling_background_admitted'] is False
    for s in rows(r)[::30]:
        rr=uniform_radius(pc,s['p'],eps)
        # Center series truncation and analytic cancellation dominate the very first samples.
        assert s['r']==pytest.approx(rr,rel=2e-5,abs=2e-3)
        assert s['m']==pytest.approx(4*math.pi*eps*s['r']**3/(3*C**2),rel=8e-8)

def test_three_level_ode_refinement_against_exact_solution(star):
    exact=uniform_radius(1e34,1e25,1e35)
    out=[star({'eos_kind':'uniform_benchmark','p_cut_pa':1e25,'rtol':t,'atol_r_m':t,
               'atol_mass_kg':1e30*t*1e-3,'atol_nu':t*1e-3}) for t in [1e-5,1e-7,1e-9]]
    err=[abs(x['radius_m']/exact-1) for x in out]
    assert err[-1]<2e-8 and err[-1]<err[0]/10,err
    assert out[-1]['accepted']>out[0]['accepted']

def test_center_and_cutoff_are_separate_refinements(star):
    a=[star({'r0_m':x,'rtol':1e-10}) for x in [10.,1.,0.1]]
    assert abs(a[1]['radius_m']-a[2]['radius_m'])<abs(a[0]['radius_m']-a[1]['radius_m'])+2e-6
    b=[star({'p_cut_pa':x,'rtol':1e-10}) for x in [1e26,1e24,1e22]]
    assert b[0]['radius_m']<b[1]['radius_m']<b[2]['radius_m']
    assert b[0]['mass_kg']<b[1]['mass_kg']<b[2]['mass_kg']

def test_weak_field_polytrope_has_lane_emden_limit(star):
    # Gamma=2 has n=1 Lane-Emden R=pi*a in Newtonian limit, a^2=K_mass/(2*pi*G).
    km=1e-56/(1.5e-10/C**2)**2
    Rn=math.pi*math.sqrt(km/(2*math.pi*G))
    a=star({'pc_pa':1e27,'p_cut_pa':1e14,'rtol':1e-10})
    b=star({'pc_pa':1e26,'p_cut_pa':1e14,'rtol':1e-10})
    assert a['status']==b['status']=='ok',(a,b)
    assert abs(b['radius_m']/Rn-1)<abs(a['radius_m']/Rn-1)
    assert abs(b['radius_m']/Rn-1)<3e-4

def test_stellar_run_metadata_is_not_support_only(star):
    import json
    r=star();meta=json.loads(Path(r['_evidence']['run_receipt']).read_text())
    c=meta['contract']
    assert c['model']=='static_spherical_TOV'
    assert c['surface']=='finite_cutoff_schwarzschild_time'
    assert c['integration_coordinate']=='log_pressure'
    assert c['units']=='SI' and c['interface']=='smooth_only'

def test_center_series_coefficient_independent(star):
    r=star({'r0_m':10.,'emit_profile':True});a=rows(r)[0]
    pc=1e34;ec=1.5e-10*1e45+pc;r0=10.
    alpha=2*math.pi*G*(ec+pc)*(ec+3*pc)/(3*C**4)
    assert a['m']==pytest.approx(4*math.pi*ec*r0**3/(3*C**2),rel=2e-14)
    assert a['p']==pytest.approx(pc-alpha*r0*r0,rel=2e-14)
