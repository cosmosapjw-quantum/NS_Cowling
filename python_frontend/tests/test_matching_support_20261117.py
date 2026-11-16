"""Declared support stays a precondition at the public matching boundary."""
import json

import numpy as np
import pytest

from compose_eos import EOSData
from compose_eos import api


def _table(tmp_path, phases):
    n = np.array([.003, .01, .02, .03, .04, .075])
    z = np.zeros(len(n))
    data = EOSData(nb=n, mn=939., mp=938., q_extra=z.copy(), Iphase=np.asarray(phases),
                   **{f"Q{i}": z.copy() for i in range(1, 8)},
                   **{key: z.copy() for key in ("Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh")},
                   beta_equilibrated=True, source_dir=str(tmp_path))
    (tmp_path / "eos.nb").write_text("1\n6\n" + "\n".join(map(str, n)) + "\n")
    return data


def _support(tmp_path):
    (tmp_path / "SUPPORT.json").write_text(json.dumps({
        "scope": "computed_native_samples", "cross_phase_interpolation": "UNAPPROVED",
        "cross_segment_interpolation": "UNAPPROVED"}))


def _must_not_match(*args):
    raise RuntimeError("native matching entered")


def test_declared_unapproved_cross_phase_input_rejected_before_native(tmp_path, monkeypatch):
    data = _table(tmp_path, [1, 1, 1, 9, 9, 9])
    _support(tmp_path)
    monkeypatch.setattr(api, "run_native", _must_not_match)
    with pytest.raises(ValueError, match="support.*phase|phase.*support"):
        api.match_crust_core(data, data, method="css", n_bridge=8)


def test_in_memory_cross_segment_concatenation_rejected(tmp_path, monkeypatch):
    data = _table(tmp_path, [1] * 6)
    _support(tmp_path)
    data.nb[0] = .002  # outside this segment, the historical lower-crust endpoint
    monkeypatch.setattr(api, "run_native", _must_not_match)
    with pytest.raises(ValueError, match="support.*source|source.*support"):
        api.match_crust_core(data, data, method="css", n_bridge=8)


@pytest.mark.parametrize("declared", [False, True])
def test_constant_noncanonical_phase_is_not_inferred_to_be_invalid(tmp_path, monkeypatch, declared):
    data = _table(tmp_path, [9] * 6)
    if declared:
        _support(tmp_path)
    monkeypatch.setattr(api, "run_native", _must_not_match)
    with pytest.raises(RuntimeError, match="native matching entered"):
        api.match_crust_core(data, data, method="css", n_bridge=8)


def test_sharp_maxwell_api_also_obeys_declared_support(tmp_path, monkeypatch):
    data = _table(tmp_path, [1, 1, 1, 9, 9, 9])
    _support(tmp_path)
    monkeypatch.setattr(api, "run_native", _must_not_match)
    with pytest.raises(ValueError, match="support.*phase|phase.*support"):
        api.construct_maxwell_transition(data, data, assume_local_charge_neutrality=True)


@pytest.mark.parametrize("change", ["edit", "delete"])
def test_sealed_support_cannot_be_removed_or_rewritten(tmp_path, monkeypatch, change):
    import hashlib
    data = _table(tmp_path, [9] * 6)
    _support(tmp_path)
    path = tmp_path / "SUPPORT.json"
    (tmp_path / "RESULT.json").write_text(json.dumps({"output_sha256": {
        "SUPPORT.json": hashlib.sha256(path.read_bytes()).hexdigest()}}))
    if change == "edit":
        path.write_text('{}')
    else:
        path.unlink()
    monkeypatch.setattr(api, "run_native", _must_not_match)
    with pytest.raises(ValueError, match="support.*sha256|support.*missing"):
        api.match_crust_core(data, data, method="css", n_bridge=8)
