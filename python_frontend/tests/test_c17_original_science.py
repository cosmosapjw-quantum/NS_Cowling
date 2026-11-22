"""C17 original-source extended branch scientific receipt tests.

Read the actual bounded M2/M3 source outputs (not fabricated mock data).
Verify physical admissions, explicit block propagation, original EOS identity,
phase-selective overlap gating and deterministic separate scan payloads.
"""
from __future__ import annotations
import gzip,hashlib,json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[2]
D=ROOT/"data/cowling_c17_extended_branch_20261122"
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"

@pytest.fixture(scope="module",params=("M2","M3"))
def original(request):
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SHA
    model=request.param
    obj=json.loads((D/"REAL_P2_INITIAL"/f"{model}.json").read_text())
    compressed=(D/"REAL_P2_INITIAL"/f"{model}_all_scan_points.json.gz").read_bytes()
    assert hashlib.sha256(compressed).hexdigest()==obj["original_scan_payload_gzip_sha256"]
    raw=gzip.decompress(compressed)
    assert hashlib.sha256(raw).hexdigest()==obj["original_scan_payload_uncompressed_sha256"]
    points=json.loads(raw)
    return model,obj,points

def test_c17_source_locked_20_density_controls_per_model(original):
    model,data,scans=original
    assert data["original_source_archive_sha256"]==SHA
    assert len(data["densities"])==20
    assert len(scans)==20
    assert data["Maxwell_small_core_onset"]["claim_scope"].startswith("small-core")
    for row in data["densities"]:
        assert row["n_c_fm3"]>data["Maxwell_small_core_onset"]["high_density_boundary_fm3"]
        if "scan_point_failures" not in row:
            assert row.get("status")=="TOV_OR_GEOMETRY_REJECTED"
            continue
        assert row["scan_point_failures"]==0
        key=f"{row['density_factor']:.5f}"
        orig=scans[key]
        assert row["native_scan_points_sha256"]==hashlib.sha256(
            json.dumps(orig,sort_keys=True,separators=(",",":")).encode()).hexdigest()
        assert row["native_scan_rhs_evaluations_total"]==sum(r["native_rhs_evaluations"] for r in orig)
        assert row["scan_sign_change_candidate_count"]==len(row["candidate_eigenvalues"])
        for c in row["candidate_eigenvalues"]:
            assert c["status"] in ("PHYSICAL_ROOT_ADMITTED","ROOT_OR_INTERFACE_REJECTED")
            if c["status"]=="PHYSICAL_ROOT_ADMITTED":
                assert c["source_eos_rhs_callback_calls"]==0
                assert c["relative_join_mismatch"]<1e-5
                assert c["phasewise_nodal_backend"]=="fortran_phasewise"
        for name,nodes in (("f",0),("p1",1)):
            b=row["local_branches"][name]
            if b["frequency_hz"] is None:
                assert b["status"] not in ("OVERLAP_AND_NODE_CONTINUATION_ADMITTED","LOCAL_NODE_CONTINUATION_ADMITTED")
            else:
                assert b["W_nodes"]==b["V_nodes"]==nodes
                assert b["relative_join_mismatch"]<1e-5
                if b["status"]=="OVERLAP_AND_NODE_CONTINUATION_ADMITTED":
                    overlap=b["cross_density_shape_overlap"]
                    assert overlap["high_phase_overlap"]>=.95
                    assert overlap["low_phase_overlap"]>=.95

def test_c17_explicit_maxwell_onset_threshold_classification(original):
    model,data,scans=original
    t=data["Maxwell_small_core_onset"]
    if model=="M2":
        assert t["exceeds_Lindblom_small_core_critical"] is True
        assert t["jump_minus_threshold"]>0
        assert .057<t["factor_of_reference_central_density"]<.058
    else:
        assert t["exceeds_Lindblom_small_core_critical"] is False
        assert t["jump_minus_threshold"]<0
        assert .488<t["factor_of_reference_central_density"]<.490


def test_c17_blocked_branch_does_not_get_fabricated_frequency(original):
    model,data,scans=original
    rows={r["density_factor"]:r for r in data["densities"]}
    if model=="M2":
        assert rows[.55]["local_branches"]["p1"]["frequency_hz"] is None
        assert rows[.6]["local_branches"]["p1"]["frequency_hz"] is not None
    else:
        assert rows[1.05]["local_branches"]["p1"]["status"]=="PHASE_SHAPE_OVERLAP_BELOW_THRESHOLD"
        diag=rows[1.05]["local_branches"]["p1"]["candidate_overlap_diagnostics"]
        assert diag[0]["diagnostic"]["high_phase_overlap"]<.95
        assert rows[1.10]["local_branches"]["p1"]["frequency_hz"] is None


def test_c17_one_bounded_midpoint_does_not_change_acceptance_gate():
    path=D/"BOUNDED_DENSITY_MIDPOINTS.json"
    data=json.loads(path.read_text())
    assert data["original_source_sha256"]==SHA
    assert data["max_midpoints_per_edge"]==1
    assert data["physical_acceptance_change"] is False
    assert len(data["results"])==4
    success={(r["model"],r["mode"],r["factors"][0],r["factors"][-1])
             for r in data["results"] if r["status"]=="TWO_HOP_FORTRAN_CONTINUATION_ADMITTED"}
    assert ("M2","p1",.6,.55) in success
    assert ("M3","p1",1.02,1.05) in success
