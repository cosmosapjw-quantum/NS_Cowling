"""Generator transport closure: support evidence, finite rows and sparse masks.

Artificial fixtures below test transport, not scientific EOS validity. The one
historical TF readback uses the retained native cohort; no campaign is rerun.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from compose_eos.api import read_compose
from compose_eos.native_generation import verify_native_generation
from test_native_generation_receipts import core_request, evidence

ROOT = Path(__file__).resolve().parents[2]


def test_native_support_is_hashed_when_present(tmp_path):
    request = core_request()
    evidence(tmp_path, request)
    support = tmp_path / "SUPPORT.json"
    support.write_text(json.dumps({"scope": "computed_native_samples",
                                  "cross_phase_interpolation": "UNAPPROVED"}))
    result = verify_native_generation(tmp_path, request)
    assert result["output_sha256"]["SUPPORT.json"] == hashlib.sha256(support.read_bytes()).hexdigest()


@pytest.mark.parametrize("mutation", ["change", "delete"])
def test_sealed_support_cannot_change_or_disappear(tmp_path, mutation):
    request = core_request()
    record = evidence(tmp_path, request)
    support = tmp_path / "SUPPORT.json"
    support.write_text(json.dumps({"scope": "computed_native_samples",
                                  "cross_phase_interpolation": "UNAPPROVED"}))
    record["output_sha256"] = verify_native_generation(tmp_path, request)["output_sha256"]
    # Explicitly seal it so this catches ignored expected hashes independently
    # of whether the unsealed readback above already includes SUPPORT.json.
    record["output_sha256"]["SUPPORT.json"] = hashlib.sha256(support.read_bytes()).hexdigest()
    (tmp_path / "RESULT.json").write_text(json.dumps(record))
    if mutation == "change":
        support.write_text(json.dumps({"cross_phase_interpolation": "APPROVED",
                                      "global_minimum_proven": True}))
    else:
        support.unlink()
    with pytest.raises(ValueError, match="checksum|sealed.*missing"):
        verify_native_generation(tmp_path, request)


@pytest.mark.parametrize("column", range(3, 10))
@pytest.mark.parametrize("value", ["NaN", "Infinity"])
def test_nonfinite_native_thermodynamics_rejected(tmp_path, column, value):
    request = core_request()
    evidence(tmp_path, request)
    path = tmp_path / "eos.thermo"
    lines = path.read_text().splitlines()
    row = lines[2].split()
    row[column] = value
    lines[2] = " ".join(row)
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="readback failed"):
        verify_native_generation(tmp_path, request)


def test_historical_tf_nonuniform_composition_masks_remain_unavailable():
    path = ROOT / "data/research_native_20261116/public_cohort/GPPVA_DD2/attempt_001/output/generated"
    table = read_compose(path)
    nonuniform = table.Iphase != 6
    assert np.any(nonuniform)
    for name in ("Yn_present", "Yp_present", "Abar_present", "Zbar_present", "Yh_present"):
        assert not np.any(getattr(table, name)[nonuniform]), name
    assert np.all(table.Ye_present[nonuniform])
    assert np.all(table.Ymu_present[nonuniform])
    assert not np.any(table.C3_applicable[nonuniform])
    assert not np.any(table.C4_applicable[nonuniform])
