import pytest
from compose_cowling.mass_sequence import extract_stable_branch, find_mass_bracket, paired_sensitivity


def point(pc, mass, radius=12000.0, index=0):
    return {"scan_index": index, "pc_pa": pc, "mass_kg": mass, "radius_m": radius, "status": "ok"}


def test_first_turning_point_truncates_the_branch():
    data = [point(1,.8,index=0), point(2,1.2,index=1), point(3,1.5,index=2),
            point(4,1.49,index=3), point(5,1.45,index=4), point(6,1.4,index=5)]
    branch = extract_stable_branch(data)
    assert [p["pc_pa"] for p in branch.points] == [1,2,3]
    assert branch.turning_point_found
    assert branch.maximum_mass_kg == pytest.approx(1.5)
    assert branch.status == "first_turning_point_resolved"


def test_failed_tov_sample_is_a_hard_sequence_barrier():
    data = [point(1,.8,index=0), point(2,1.0,index=1),
            {"scan_index":2,"pc_pa":3,"status":"failed"},
            point(4,1.3,index=3), point(5,1.4,index=4)]
    branch = extract_stable_branch(data)
    assert [p["pc_pa"] for p in branch.points] == [1,2]
    assert branch.status == "sequence_gap_truncated"
    assert branch.maximum_mass_kg is None
    assert branch.branch_endpoint_mass_kg == pytest.approx(1.0)


def test_mass_bracket_uses_adjacent_points_and_log_interpolation():
    branch = extract_stable_branch([point(1,.5,10000,0),point(4,1.5,16000,1),point(16,2,18000,2)])
    result = find_mass_bracket(branch, 1.0)
    assert result["pc_lo_pa"] == 1
    assert result["pc_hi_pa"] == 4
    assert result["pc_seed_pa"] == pytest.approx(2.0)
    assert result["radius_guess_m"] == pytest.approx((10000*16000)**0.5)
    assert find_mass_bracket(branch, 2.1) is None


def test_six_method_pair_requires_six_valid_modes_and_pair_mean():
    methods = ["css","cubic_F","cubic_P","maxwell","quintic_F","quintic_P"]
    rows = [{"crust":"C","core":"K","method":m,"mass_msun":1.0,
             "f_hz":2000+i,"p1_hz":6000+10*i,"valid":True} for i,m in enumerate(methods)]
    paired, deltas = paired_sensitivity(rows, methods)
    mean_f = sum(r["f_hz"] for r in rows)/6
    assert paired[0]["pair_complete"]
    assert paired[0]["f_mean_hz"] == pytest.approx(mean_f)
    assert deltas[0]["delta_f_percent"] == pytest.approx(100*(2000/mean_f-1))
    incomplete, _ = paired_sensitivity(rows[:-1], methods)
    assert not incomplete[0]["pair_complete"]
    assert incomplete[0]["n_valid_methods"] == 5
    assert incomplete[0]["f_range_percent"] is None
