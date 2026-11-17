"""Native parameter-injected uniform RMF physical checks (no fitted source offsets)."""
import os
from pathlib import Path
import numpy as np
import pytest
from compose_eos.rmf_candidate import build_driver, load_parameters, evaluate_states

ROOT = Path(os.environ.get("RMF_CANDIDATE_ROOT", Path(__file__).resolve().parents[2]))
SOURCE = Path(os.environ.get("RMF_CANDIDATE_SOURCE", ROOT / "data/new_model_maxwell_20261117/sources/PARAMETER_INPUTS.json"))

@pytest.fixture(scope="module")
def executable(tmp_path_factory):
    return build_driver(ROOT, tmp_path_factory.mktemp("rmf-driver"))

@pytest.mark.parametrize("model", [
    "NL3_PROVIDER_PRECISION", "PK1_CANONICAL", "TM1_CANONICAL",
    "DDMEX_PROVIDER_ROUNDED_N0", "DDMEX_ORIGINAL_TABLE_RHO0"])
def test_physical_beta_and_first_law(executable, model, tmp_path):
    p = load_parameters(SOURCE, model)
    density=np.array([.015,.05,.1,.2,.4,.8,1.])
    s=evaluate_states(executable,p,density,work_dir=tmp_path/model)
    assert np.all(s["effective_mass_n"]>0)
    assert np.all(s["effective_mass_p"]>0)
    assert np.all(s["scalar_curvature"]>0)
    assert np.all(s["scalar_root_count"]>=1)
    assert np.max(np.abs(s["beta_residual_mev"]))<2e-8
    assert np.max(np.abs(s["charge_residual_fm3"]))<2e-13
    assert np.max(np.abs(s["euler_residual_mevfm3"]))<2e-8
    assert np.max(np.abs(s["scalar_residual_fm3"]))<2e-11
    for n in [.1,.4]:
        h=n*2e-4
        v=evaluate_states(executable,p,[n-h,n,n+h],work_dir=tmp_path/model/f"fd-{n}")
        dedn=(v["epsilon"][2]-v["epsilon"][0])/(2*h)
        dPdn=(v["pressure"][2]-v["pressure"][0])/(2*h)
        dmudn=(v["mu_B"][2]-v["mu_B"][0])/(2*h)
        assert abs(dedn-v["mu_B"][1])/abs(v["mu_B"][1])<2e-6
        assert abs(dPdn-n*dmudn)/max(1,abs(dPdn))<2e-6

def test_reject_invalid_density(executable, tmp_path):
    p=load_parameters(SOURCE,"NL3_PROVIDER_PRECISION")
    for density in [[0.],[-.1],[float("nan")],[float("inf")]]:
        with pytest.raises(ValueError,match="positive"):
            evaluate_states(executable,p,density,work_dir=tmp_path/"invalid")
