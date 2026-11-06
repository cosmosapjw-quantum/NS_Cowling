"""The integral check must not construct a tautology using the eigenfrequency."""
from pathlib import Path
import json, math
import numpy as np
import pytest
from tools.variational_check import rayleigh_ratio, hermite, profile_integrals, VariationalError
from tools.native_eos_reference import NativeEOS
from tools.precision_study import ROOT, run_one

def test_eigenvalue_dependent_robin_surface_mass():
    # y=sin(k*x), y(0)=0, y'(1)=k² B y(1), k=pi/4.
    k=math.pi/4;B=1/k
    numerator=k*k*(.5+math.sin(2*k)/(4*k))
    bulk=.5-math.sin(2*k)/(4*k);surface=B*math.sin(k)**2
    assert rayleigh_ratio(numerator,bulk,surface)==pytest.approx(k*k,rel=1e-14)
    assert abs(rayleigh_ratio(numerator,bulk,0.)-k*k)>.1

def test_invalid_integral_data_are_rejected():
    for n,d,b in [(float('nan'),1.,0.),(1.,0.,0.),(1.,1.,-1.),(-1.,1.,0.)]:
        with pytest.raises(VariationalError):rayleigh_ratio(n,d,b)

def test_hermite_derivative_does_not_require_frequency():
    s=np.linspace(0,1,21);v,ds=hermite(1.,4.,0.,6.,1.,s)
    assert np.max(abs(v-(1+3*s*s)))<1e-14
    assert np.max(abs(ds-6*s))<1e-14

@pytest.fixture(scope='module')
def shape_case(tmp_path_factory):
    out=tmp_path_factory.mktemp('variational')/'star'
    data=run_one(out,dict(alpha=0.,rtol=1e-11,mode_rtol=1e-11,step_cap=.05),workspace='calculation')
    a,b=[np.genfromtxt(out/f'result_mode1_{x}.csv',delimiter=',',names=True) for x in ['inner','outer']]
    eos=NativeEOS(ROOT/'data/normalized/compose_69.csv')
    return data,a,b,eos

def test_integral_is_frequency_independent_and_scale_invariant(shape_case):
    meta,a,b,eos=shape_case
    p=profile_integrals(meta,a,b,eos,8)
    changed={**meta,'f_hz':meta['f_hz']*1.1,'p1_hz':-999.}
    q=profile_integrals(changed,a,b,eos,8)
    assert p['rayleigh_hz']==q['rayleigh_hz']
    a2=a.copy();b2=b.copy()
    for c in ['w_over_r3','v_over_r2']:a2[c]*=-3;b2[c]*=7
    z=profile_integrals(meta,a2,b2,eos,8)
    assert z['rayleigh_hz']==pytest.approx(p['rayleigh_hz'],abs=1e-9)
    assert abs(p['rayleigh_hz']-meta['f_hz'])<2e-6

def test_endpoint_duplicates_do_not_hide_physical_jumps(shape_case):
    meta,a,b,eos=shape_case
    z=np.insert(a,10,a[10]);good=profile_integrals(meta,z,b,eos,8)
    assert math.isfinite(good['rayleigh_hz'])
    z['v_over_r2'][11]*=1.01
    with pytest.raises(VariationalError,match='conflicting'):profile_integrals(meta,z,b,eos,8)

def test_conflicting_matching_geometries_are_rejected(shape_case):
    meta,a,b,eos=shape_case;b=b.copy();b['r_m'][-1]*=1.1
    with pytest.raises(VariationalError):profile_integrals(meta,a,b,eos,8)

def test_wrong_shape_is_not_hidden_by_saved_eigenfrequency(shape_case):
    meta,a,b,eos=shape_case;a=a.copy()
    q=profile_integrals(meta,a,b,eos,8)
    a['v_over_r2']*=1+.001*np.sin(np.pi*a['r_m']/a['r_m'][-1])**2
    bad=profile_integrals(meta,a,b,eos,8)
    assert abs(bad['rayleigh_hz']-q['rayleigh_hz'])>1e-4

def test_integral_gate_never_accepts_shared_residual_drift():
    from tools.variational_check import assess_integral_ladder
    budgets={x:{'empirical_allowance_hz':1e-6} for x in ['f','p1']}
    rows=[{side:{'modes':{m:{'difference_hz':v,'quadrature_change_hz':0.,'shooting_hz':2500.} for m in ['f','p1']}} for side in ['reference','variant']} for v in [1e-3,1e-4,1e-5,1e-8]]
    assert assess_integral_ladder(rows,budgets)['status']=='integral_consistent'
    for row in rows[-2:]:row['variant']['modes']['p1']['difference_hz']=.01
    assert assess_integral_ladder(rows,budgets)['status']=='integral_consistency_unresolved'
