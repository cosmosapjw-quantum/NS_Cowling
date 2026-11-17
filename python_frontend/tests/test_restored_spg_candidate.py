import numpy as np
import pytest
from compose_eos.restored_spg_candidate import RestoredSPGCandidate,MN,MP
from compose_eos.candidate_gdfm_cld import reduced_cell,default_surface_parameters,MASS

def test_species_linear_gradient_is_analytic_and_changes_beta_functional():
    c=RestoredSPGCandidate('SPG_M1')
    n=.003;z=np.array([.3,.14,.16])
    e,g,_=c.cell(n,z)
    original,oldg,_=reduced_cell('SPG_M1',n,z,surface=c.surface)
    assert e-original==pytest.approx(MN-MASS+(MP-MN)*z[0]*z[2],abs=1e-12)
    assert g[2]-oldg[2]==pytest.approx((MP-MN)*z[0],abs=1e-12)
    for j,h in enumerate([1e-5,1e-6,1e-6]):
        a=z.copy();b=z.copy();a[j]+=h;b[j]-=h
        fd=(c.cell(n,a)[0]-c.cell(n,b)[0])/(2*h)
        assert fd==pytest.approx(g[j],rel=2e-5,abs=2e-6)

def test_actual_reoptimized_state_has_first_law_charge_and_stable_gas():
    c=RestoredSPGCandidate('SPG_M3');n=.003;h=n*1e-4
    s=c.state(n)
    derivative=(c.state(n+h)['epsilon']-c.state(n-h)['epsilon'])/(2*h)
    assert derivative==pytest.approx(s['mu_B'],abs=1e-5)
    assert s['epsilon']+s['pressure']==pytest.approx(n*s['mu_B'],rel=1e-14)
    assert s['gas_proton_insertion_margin_MeV']>=0
    assert s['projected_gradient_max']<1e-6
    assert s['Y_e']==pytest.approx(s['coordinates'][0]*s['coordinates'][2])
    assert s['surface_parameters']['source_sha256']==c.surface.source_sha256

@pytest.mark.parametrize('n',[0.,.0301,float('nan'),float('inf')])
def test_candidate_declared_support_refused(n):
    with pytest.raises(ValueError):RestoredSPGCandidate('SPG_M1').state(n)

def test_surface_model_mismatch_refused():
    with pytest.raises(ValueError):
        RestoredSPGCandidate('SPG_M1',surface=default_surface_parameters('SPG_M2'))

@pytest.mark.parametrize('mass',[0.,float('nan'),float('inf')])
def test_invalid_rest_energy_refused(mass):
    with pytest.raises(ValueError):RestoredSPGCandidate('SPG_M1',mn=mass)


@pytest.mark.parametrize("field",["surface","mn","mp","difference"])
def test_cached_configuration_cannot_be_mutated(field):
    from dataclasses import FrozenInstanceError
    c=RestoredSPGCandidate('SPG_M1')
    c.state(.002)
    replacement=default_surface_parameters('SPG_M2') if field=="surface" else 1.
    with pytest.raises(FrozenInstanceError):
        setattr(c,field,replacement)

def test_returned_nested_state_cannot_poison_cache():
    from copy import deepcopy
    c=RestoredSPGCandidate('SPG_M1')
    result=c.state(.002);expected=deepcopy(result)
    result["epsilon"]=-1234.
    result["surface_parameters"]["label"]="POISON"
    result["all_three_attempts"][0]["coordinates"][0]=-9.
    assert c.state(.002)==expected

def test_separate_rest_energy_configurations_have_separate_cache_entries():
    a=RestoredSPGCandidate('SPG_M1')
    b=RestoredSPGCandidate('SPG_M1',mn=MN+.1,mp=MP+.1)
    sa=a.state(.002);sb=b.state(.002)
    assert sb["mu_B"]-sa["mu_B"]==pytest.approx(.1,abs=1e-7)
    assert sb["epsilon"]-sa["epsilon"]==pytest.approx(.002*.1,abs=1e-10)


def test_stationary_saddle_is_refused_by_curvature_check():
    from compose_eos.restored_spg_candidate import _local_minimum_check
    z=np.array([.5,.14,.2]);lo=np.array([1e-7,.035,.001]);hi=np.array([1.,.25,.499])
    gradient=lambda q:np.array([q[0]-.5,-(q[1]-.14),q[2]-.2])
    with pytest.raises(RuntimeError,match="negative local curvature"):
        _local_minimum_check(gradient,z,lo,hi)

def test_artificial_bound_is_refused_even_with_zero_gradient():
    from compose_eos.restored_spg_candidate import _local_minimum_check
    lo=np.array([1e-7,.035,.001]);hi=np.array([1.,.25,.499])
    z=np.array([.5,lo[1]+1e-12,.2])
    with pytest.raises(RuntimeError,match="artificial search bound"):
        _local_minimum_check(lambda q:q-z,z,lo,hi)

def test_physical_zero_gas_boundary_has_two_dimensional_tangent_check():
    from compose_eos.restored_spg_candidate import _local_minimum_check
    z=np.array([1.,.14,.2]);lo=np.array([1e-7,.035,.001]);hi=np.array([1.,.25,.499])
    result=_local_minimum_check(lambda q:np.array([-1.,q[1]-.14,q[2]-.2]),z,lo,hi)
    assert result["curvature_free_coordinates"]==[1,2]
    assert result["minimum_scaled_hessian_eigenvalue"]>0


def test_marginal_zero_gas_boundary_inward_saddle_is_refused():
    from compose_eos.restored_spg_candidate import _local_minimum_check
    z=np.array([1.,.14,.2]);lo=np.array([1e-7,.035,.001]);hi=np.array([1.,.25,.499])
    # Gradient vanishes at t=1 but inward t displacement lowers E.
    gradient=lambda q:np.array([-(q[0]-1.),q[1]-.14,q[2]-.2])
    with pytest.raises(RuntimeError,match="negative local curvature"):
        _local_minimum_check(gradient,z,lo,hi)

def test_marginal_zero_gas_boundary_checks_inward_t_curvature():
    from compose_eos.restored_spg_candidate import _local_minimum_check
    z=np.array([1.,.14,.2]);lo=np.array([1e-7,.035,.001]);hi=np.array([1.,.25,.499])
    result=_local_minimum_check(lambda q:q-z,z,lo,hi)
    assert result["curvature_free_coordinates"]==[0,1,2]
    assert result["minimum_scaled_hessian_eigenvalue"]>0
