import numpy as np
import pytest
from compose_eos.candidate_gdfm_cld import (
    MASS, cell, two_phase_cell, reduced_cell, reduced_two_phase_cell, solve_reduced, solve_two_phase)
from compose_eos.spg_crust_generator import neutron_drip, to_eos_data

@pytest.mark.parametrize('model',['SPG_M1','SPG_M3','SPG_M5'])
def test_two_phase_energy_full_gradient(model):
    x=np.array([.3,10.,.12,.12,.003]);nb=.055
    eps,g,s=two_phase_cell(model,nb,x)
    for j,h in enumerate([2e-6,2e-4,2e-6,2e-6,2e-7]):
        a=x.copy();b=x.copy();a[j]+=h;b[j]-=h
        numeric=(two_phase_cell(model,nb,a)[0]-two_phase_cell(model,nb,b)[0])/(2*h)
        assert g[j]==pytest.approx(numeric,abs=3e-6,rel=3e-6)
    h=2e-6
    numeric=(two_phase_cell(model,nb+h,x)[0]-two_phase_cell(model,nb-h,x)[0])/(2*h)
    assert s['partial_energy_nb_MeV']==pytest.approx(numeric,abs=3e-6,rel=3e-8)
    assert x[1]*g[1]==pytest.approx(2*s['epsilon_coulomb']-s['epsilon_surface'],abs=1e-14)


def test_two_phase_zero_gas_protons_recovers_original_functional():
    x=[.12,8.,.15,.14]
    a=cell('SPG_M1',.03,x);b=two_phase_cell('SPG_M1',.03,[*x,0.])
    assert a[0]==pytest.approx(b[0],abs=1e-14)
    np.testing.assert_allclose(a[1],b[1][:4],atol=1e-13,rtol=1e-13)


def test_reduced_positive_gas_gradient_and_virial():
    z=np.array([.62,.085,.052,.005]);nb=.078
    e,g,s=reduced_two_phase_cell('SPG_M3',nb,z)
    for j,h in enumerate([2e-6,2e-7,2e-7,2e-6]):
        a=z.copy();b=z.copy();a[j]+=h;b[j]-=h
        num=(reduced_two_phase_cell('SPG_M3',nb,a)[0]-reduced_two_phase_cell('SPG_M3',nb,b)[0])/(2*h)
        assert g[j]==pytest.approx(num,abs=2e-6,rel=2e-5)
    assert s['state']['epsilon_surface']==pytest.approx(2*s['state']['epsilon_coulomb'],rel=1e-13)


def test_old_m3_abnormal_line_search_is_resolved_by_stationary_polish():
    best,attempts=solve_reduced('SPG_M3',.005011872336)
    assert best['status']=='STATIONARY_CANDIDATE'
    assert best['max_projected_gradient']<1e-8
    assert best['minimum_scaled_hessian_eigenvalue']>0


def test_positive_proton_gas_restores_insertion_stationarity():
    best,_=solve_two_phase('SPG_M1',.0745)
    assert best['status']=='STATIONARY_CANDIDATE'
    assert best['state']['gas_proton_fraction']>0
    assert abs(best['state']['gas_proton_insertion_margin_MeV'])<1e-7
    assert best['max_projected_gradient']<1e-8


def test_no_neutron_gas_boundary_has_envelope_pressure():
    n=1e-5;h=1e-8
    best,_=solve_reduced('SPG_M1',n)
    assert best['physical_boundary']=='NEUTRON_GAS_ABSENT'
    plus=solve_reduced('SPG_M1',n+h)[0]['energy']
    minus=solve_reduced('SPG_M1',n-h)[0]['energy']
    mu=(plus-minus)/(2*h)
    assert best['mu_B_MeV']==pytest.approx(mu,abs=3e-6)
    assert best['mu_B_MeV']<MASS


def test_drip_condition_and_composition_export_contract():
    d=neutron_drip('SPG_M1')
    assert abs(d['mu_B_MeV']-MASS)<1e-6
    rows=[{'n_B':n,'best':solve_two_phase('SPG_M1',n)[0]} for n in [.001,.02,.0745]]
    data=to_eos_data(rows)
    np.testing.assert_allclose(data.Yn+data.Yp+data.Yh,1,atol=2e-15)
    np.testing.assert_allclose(data.Yp+data.Zav/data.Aav*data.Yh,data.Ye,atol=2e-15)
    assert data.C3_applicable.all() and data.C4_applicable.all()
    rows[1]['best']['status']='HOLD'
    with pytest.raises(ValueError):to_eos_data(rows)
