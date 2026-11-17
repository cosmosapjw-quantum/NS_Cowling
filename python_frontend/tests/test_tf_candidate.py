"""Native spatial TF registry has the frozen F homogeneous physics limit."""
from pathlib import Path
import numpy as np
import pytest
from compose_eos import tf_candidate, rmf_candidate
from compose_eos.gppva_rmf import beta_state

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/new_model_maxwell_20261117/sources/PARAMETER_INPUTS.json"

@pytest.fixture(scope="module")
def drivers(tmp_path_factory):
    where = tmp_path_factory.mktemp("tf_candidate")
    tf = tf_candidate.build_driver(ROOT, where / "tf_build")
    bulk = rmf_candidate.build_driver(ROOT, where / "bulk_build")
    return tf, bulk, where

@pytest.mark.parametrize("model", ["TM1_CANONICAL", "DDMEX_ORIGINAL_TABLE_RHO0", "NL3_PROVIDER_PRECISION"])
def test_new_registry_uniform_limit_matches_frozen_parameter_driver(drivers, model):
    tf, bulk, where = drivers
    params = rmf_candidate.load_parameters(SOURCE, model)
    states = rmf_candidate.evaluate_states(bulk, params, [0.04], work_dir=where / (model+"_bulk"))
    s = tf_candidate.run_state(tf, model, 0.04, 16, 3, "uniform", 30.,
                               work_dir=where / (model+"_tf"), evaluation_only=True,
                               initial_nn=np.full(16,.04*(1-states["Y_p"][0])),
                               initial_np=np.full(16,.04*states["Y_p"][0]))
    assert s["status"] == "PROFILE_EVALUATED", s
    assert s["geometry"] == "uniform"
    assert not s["accepted_stationary"]
    assert s["muon_threshold_crossed"] == 0
    assert abs(s["epsilon"]-states["epsilon"][0]) < 2e-8
    assert abs(s["mu_B"]-states["mu_B"][0]) < 1e-4
    assert abs(s["pressure"]-states["pressure"][0]) < 5e-6
    assert abs(s["Y_p"]-states["Y_p"][0]) < 2e-6

def test_existing_ddme2_uniform_limit_preserved(drivers):
    tf, _, where = drivers
    b = beta_state("DDME2", 0.01, include_muons=False)
    s = tf_candidate.run_state(tf, "DDME2", 0.01, 16, 3, "uniform", 40.,
                               work_dir=where / "DDME2_tf")
    assert s["accepted_stationary"], s
    assert abs(s["epsilon"]-b["epsilon"]) < 2e-8
    assert abs(s["pressure"]-b["pressure"]) < 5e-6

def test_exported_tf_transport_has_no_invented_abundance_masks(drivers):
    from compose_eos.api import diagnose_maxwell_phase
    tf, _, where = drivers
    rows=[tf_candidate.run_state(tf,"DDME2",n,16,3,"uniform",40.,
          work_dir=where/f"transport_{i}") for i,n in enumerate([.008,.01,.012])]
    data,_=tf_candidate.export_compose(where/"transport_export",rows,mn=939.,mp=939.,model="DDME2")
    assert not np.any(data.Yn_present)
    assert not np.any(data.Yp_present)
    assert not np.any(data.C3_applicable)
    diagnostic=diagnose_maxwell_phase(data)
    assert diagnostic["global_validation_category"] != "COMPOSITION_UNSUPPORTED"
