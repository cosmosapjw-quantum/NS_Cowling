"""T15 difference classification without relative-error singularities at the null."""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from tools.runtime import CaseRunner

ROOT = Path(__file__).resolve().parents[1]


def _run(tmp_path, **kw):
    defaults = dict(f_ref_hz=2500.0, f_cmp_hz=2500.0, p1_ref_hz=7200.0, p1_cmp_hz=7200.0, resolution_hz=1e-5)
    defaults.update(kw)
    runner = CaseRunner(ROOT, tmp_path / "runs", workspace=os.environ.get("NS_WORKSPACE", "solutions"))
    return runner.run_case("mode_difference", defaults)


def test_exact_null_uses_absolute_resolution_not_relative_error(tmp_path):
    r = _run(tmp_path)
    assert r["status"] == "ok", r
    assert r["delta_f_hz"] == 0.0
    assert r["delta_p1_hz"] == 0.0
    assert r["f_interpretation"] == "below_resolution"
    assert r["p1_interpretation"] == "below_resolution"
    assert r["f_relative_defined"] is True
    assert r["p1_relative_defined"] is True
    assert r["delta_f_relative"] == 0.0
    assert r["delta_p1_relative"] == 0.0


def test_pair_resolution_is_separate_from_mode_identity(tmp_path):
    r = _run(tmp_path, f_cmp_hz=2500.000001, p1_cmp_hz=7200.25)
    assert r["status"] == "ok", r
    assert r["f_interpretation"] == "below_resolution"
    assert r["p1_interpretation"] == "resolved"
    assert r["delta_p1_hz"] == pytest.approx(0.25, rel=0, abs=1e-12)
    assert r["delta_p1_relative"] == pytest.approx(0.25 / 7200.0, rel=2e-14)
    assert r["identity_compared"] is False
