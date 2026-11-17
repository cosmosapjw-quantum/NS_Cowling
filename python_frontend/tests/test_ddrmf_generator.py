"""Composition serialization of independently solved homogeneous DD-RMF matter."""
import numpy as np
import pytest

from compose_eos.api import read_compose
from compose_eos.ddrmf_generator import generate_ddrmf_core


@pytest.mark.parametrize("model", ["GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW"])
def test_generated_composition_roundtrips_and_covers_muon_threshold(tmp_path, model):
    target = tmp_path / model
    receipt = generate_ddrmf_core(model, [0.001, 0.03, 0.15, 0.5], target)
    data = read_compose(target / "compose")
    assert receipt["status"] == "GENERATED_AND_NATIVE_READBACK_VERIFIED"
    assert data.composition_available and data.composition_file_present
    np.testing.assert_allclose(data.Yn + data.Yp, 1.0, atol=2e-14, rtol=0)
    np.testing.assert_allclose(data.Yp, data.Ye + data.Ymu, atol=2e-11, rtol=0)
    assert data.Ymu[0] == 0 and data.Ymu[-1] > 0
    assert data.Ymu_present.all()
    assert (data.Yh == 0).all()
    assert receipt["composition"]["C3"]["status"] == "CHECKED"
    assert receipt["composition"]["C4"]["status"] == "CHECKED"
    assert receipt["max_abs_beta_residual_mev"] < 2e-8
    np.testing.assert_allclose(data.eps + data.p, data.nb * data.mub, rtol=2e-10)
    with pytest.raises(FileExistsError):
        generate_ddrmf_core(model, [0.01, 0.1, 0.2], target)


@pytest.mark.parametrize("density", [[0.1, 0.1, 0.2], [0.2, 0.1, 0.3],
                                    [0.1, float("nan"), 0.3], [0, 0.1, 0.2]])
def test_invalid_grids_rejected_before_writing(tmp_path, density):
    target = tmp_path / "invalid"
    with pytest.raises(ValueError):
        generate_ddrmf_core("GPPVA_DD2", density, target)
    assert not target.exists()


def test_unknown_model_rejected_before_writing(tmp_path):
    with pytest.raises(ValueError, match="supported"):
        generate_ddrmf_core("APR_APR", [0.1, 0.2, 0.3], tmp_path / "unknown")
    assert not (tmp_path / "unknown").exists()


def test_failed_native_solve_retains_failed_receipt(tmp_path, monkeypatch):
    import json
    import compose_eos.ddrmf_generator as generator

    def fail(model, density):
        raise RuntimeError("test native failure")
    monkeypatch.setattr(generator, "evaluate_ddrmf", fail)
    target = tmp_path / "failed"
    with pytest.raises(RuntimeError, match="test native failure"):
        generator.generate_ddrmf_core("GPPVA_DD2", [0.01, 0.1, 0.2], target)
    receipt = json.loads((target / "generation_receipt.json").read_text())
    assert receipt["status"] == "FAILED"
    assert receipt["failure_type"] == "RuntimeError"
    assert not (target / "compose").exists()
