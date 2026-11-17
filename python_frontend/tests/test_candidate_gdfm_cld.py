import numpy as np
import pytest
from compose_eos.candidate_gdfm_cld import cell

@pytest.mark.parametrize('model', ['SPG_M1', 'SPG_M3', 'SPG_M5'])
def test_full_energy_gradient_matches_independent_differences(model):
    x = np.array([0.12, 8.0, 0.15, 0.14])
    energy, grad, state = cell(model, .03, x)
    assert np.isfinite(energy)
    for i, step in enumerate([1e-6, 1e-4, 1e-6, 1e-6]):
        plus=x.copy(); minus=x.copy(); plus[i]+=step; minus[i]-=step
        numeric=(cell(model,.03,plus)[0]-cell(model,.03,minus)[0])/(2*step)
        assert grad[i] == pytest.approx(numeric, abs=2e-6, rel=2e-6)
    assert state['n_gas']*(1-x[0])+x[0]*x[2] == pytest.approx(.03, abs=1e-15)


def test_invalid_gas_is_rejected():
    with pytest.raises(ValueError):
        cell('SPG_M1', .01, [0.2, 8., 0.15, 0.1])


def test_gas_proton_insertion_margin_against_energy_difference():
    from compose_eos.gdfm_uniform import baryons
    from compose_eos.cold_leptons import free_fermion, ELECTRON_MASS
    from compose_eos.scurto_interface import spherical_phi
    from compose_eos.candidate_gdfm_cld import E2
    import math
    f,r,ni,y = x = np.array([.12,8.,.15,.14])
    _,_,s = cell('SPG_M1',.03,x)
    ng=s['n_gas']; p=ni*y
    def varying(pg):
        eg=baryons('SPG_M1',ng,pg/ng)['epsilon_baryons']
        ee=free_fermion(f*p+(1-f)*pg,ELECTRON_MASS)[0]
        ec=2*math.pi*E2*f*spherical_phi(f)*r*r*(p-pg)**2
        return (1-f)*eg+ee+ec
    h=1e-8
    numerical=(varying(h)-varying(0))/(h*(1-f))
    assert s['gas_proton_insertion_margin_MeV'] == pytest.approx(numerical, abs=0.02)
