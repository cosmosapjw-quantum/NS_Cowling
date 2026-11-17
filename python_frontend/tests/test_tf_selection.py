import pytest
from compose_eos.tf_selection import select_geometry_candidates, SUCCESS

FAMILIES=[(3,"droplet"),(2,"droplet")]

def state(dimension=3,**changes):
    value=dict(dimension=dimension,seed="droplet",n_b=.0625,epsilon=58.75,
        energy_per_baryon=940.,numerical_status="STATIONARY",muon_threshold_crossed=False,
        baryon_constraint_error=0.,charge_constraint_error=0.,
        scalar_pde_residual=0.,vector_pde_residual=0.,electron_pde_residual=0.)
    value.update(changes);return value

def test_failed_family_is_unresolved():
    failure=dict(dimension=2,seed="droplet",n_b=.0625,exception="electron line search")
    out=select_geometry_candidates([state(),failure],FAMILIES)
    assert out["status"]=="HOLD_UNRESOLVED_GEOMETRY"
    assert out["selected"] is None and out["unresolved_families"]==[(2,"droplet")]

def test_retry_resolves_family_without_erasing_exception():
    failure=dict(dimension=2,seed="droplet",n_b=.0625,exception="electron line search")
    out=select_geometry_candidates([state(),failure,state(2)],FAMILIES)
    assert out["status"]==SUCCESS and len(out["attempt_diagnostics"])==3
    assert not out["attempt_diagnostics"][1]["stationary_supported"]

def test_lower_feasible_unstationary_state_is_not_hidden():
    low=state(2,energy_per_baryon=939.,epsilon=58.6875,numerical_status="REQUIRES_REFINEMENT")
    out=select_geometry_candidates([state(),state(2),low],FAMILIES)
    assert out["status"]=="LOWER_UNCONVERGED_ATTEMPT"
    assert out["lower_unconverged_attempt_indices"]==[2]

def test_infeasible_trial_energy_does_not_block_resolved_geometry():
    low=state(2,energy_per_baryon=900.,baryon_constraint_error=.01,numerical_status="REQUIRES_REFINEMENT")
    out=select_geometry_candidates([state(),state(2),low],FAMILIES)
    assert out["status"]==SUCCESS
    assert out["attempt_diagnostics"][2]["feasibility"]=="INFEASIBLE"

def test_unknown_lower_energy_remains_unresolved():
    low=dict(dimension=2,seed="droplet",n_b=.0625,energy_per_baryon=939.)
    out=select_geometry_candidates([state(),state(2),low],FAMILIES)
    assert out["status"]=="LOWER_UNCONVERGED_ATTEMPT"

@pytest.mark.parametrize("changes",[
    {"baryon_constraint_error":1e-5},
    {"charge_constraint_error":1e-5},
    {"scalar_pde_residual":1e-5},
    {"electron_pde_residual":1e-5},
    {"epsilon":float("nan")},
    {"muon_threshold_crossed":True},
    {"n_b":.063},
])
def test_stationarity_label_does_not_override_physical_support(changes):
    out=select_geometry_candidates([state(),state(2,**changes)],FAMILIES,n_b=.0625)
    assert out["selected"] is None
    assert out["unresolved_families"]==[(2,"droplet")]
