import numpy as np
import pytest
from compose_eos.bsk22_extended_representation import ExtendedAnchoredBSk22,extended_source_state
from compose_eos import bsk22_analytic as source

@pytest.fixture
def representation():
    return ExtendedAnchoredBSk22(outer_density=.0002694320064010902,
        outer_pressure=.0004992059167689819,outer_mu_b=939.5654205200001,
        outer_epsilon=.2526497904790187)

def test_original_guard_and_lower_formula_unchanged():
    with pytest.raises(ValueError):source.state(.003)
    for n in (.0003,.001,.002):
        assert extended_source_state(n)==source.state(n)

@pytest.mark.parametrize('n',[.0021,.003,.006,.0099])
def test_extended_first_law_and_positive_composition(representation,n):
    s=representation.state(n);h=n*2e-5
    left=representation.state(n-h);right=representation.state(n+h)
    derivative=(right['epsilon']-left['epsilon'])/(2*h)
    assert derivative==pytest.approx(s['mu_b_euler'],abs=2e-7)
    assert n*s['mu_b_euler']-s['epsilon']==pytest.approx(s['pressure'],rel=1e-10)
    assert s['Y_n_free']+s['Y_p_free']+s['Y_h']==pytest.approx(1,abs=1e-15)
    assert s['Y_e']==pytest.approx(s['Z_cluster']/s['A_cluster']*s['Y_h'],abs=1e-15)
    assert s['Z_cluster']==40
    assert s['Y_p_free']==0
    assert s['pressure_density_derivative']>0

def test_extension_continuity_without_new_anchor(representation):
    edge=representation.state(.002)
    right=representation.state(np.nextafter(.002,1))
    for key in ('epsilon','pressure','mu_b_euler','Y_e','Y_n_free'):
        assert right[key]==pytest.approx(edge[key],abs=1e-12)

def test_source_composition_chemical_residual_remains_visible(representation):
    s=representation.state(.006)
    assert abs(s['source_composition_mu_n_residual_MeV'])>1e-4
    assert s['original_fit_thermodynamic_relative_residual']!=0.
    assert s['pressure_fit_thermodynamic_relative_residual']==0.

@pytest.mark.parametrize('n',[.010001,float('inf'),float('nan'),-.001])
def test_extension_refuses_outside_declared_support(representation,n):
    with pytest.raises(ValueError):representation.state(n)

def test_carrier_preserves_composition(representation):
    data,rows=representation.eos_data(np.geomspace(representation.n0,.01,17))
    np.testing.assert_allclose(data.eps,[r['epsilon'] for r in rows],rtol=2e-15)
    np.testing.assert_allclose(data.Ye,data.Yp+data.Zav/data.Aav*data.Yh,atol=1e-15)
    assert data.Yn_present.all() and data.Yh_present.all()
