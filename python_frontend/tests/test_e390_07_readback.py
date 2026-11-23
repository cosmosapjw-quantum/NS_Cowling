"""E390-07 readback regression: the verify action does NOT emit row count;
the orchestrator must read it from the 'read' action on the written CompOSE
directory. This test pins that contract."""
import sys, tempfile
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
sys.path.insert(0, str(ROOT / "scripts" / "m7"))

from compose_eos.native import run_native


def _needs_backend():
    import os
    return os.environ.get("COMPOSE_EOS_FORTRAN_EXE") is None


@pytest.mark.skipif(_needs_backend(), reason="requires the native Fortran backend")
def test_verify_action_does_not_emit_row_count():
    """The native 'verify' action emits overall/failed_checks but NOT 'n'.
    The orchestrator must NOT read row count from verify."""
    with tempfile.TemporaryDirectory(prefix="vtest_") as td:
        out_json = Path(td) / "verify.json"
        r = run_native("verify", ROOT / "eos_generator" / "eos_data" / "GPPVA_DD2_matched", out_json)
        assert r.get("status") == "ok"
        # verify emits overall (PASS/FAIL) but NOT "n"
        assert "overall" in r
        assert "n" not in r, "verify must not emit row count; orchestrator reads it from 'read'"


@pytest.mark.skipif(_needs_backend(), reason="requires the native Fortran backend")
def test_read_action_emits_row_count():
    """The native 'read' action emits 'n' (row count) -- the source the
    orchestrator must use for readback verification."""
    with tempfile.TemporaryDirectory(prefix="rtest_") as td:
        internal = Path(td) / "internal"
        internal.mkdir()
        r = run_native("read", ROOT / "eos_generator" / "eos_data" / "GPPVA_DD2_matched", internal)
        assert r.get("status") == "ok"
        assert "n" in r and int(r["n"]) > 0


@pytest.mark.skipif(_needs_backend(), reason="requires the native Fortran backend")
def test_native_pipeline_readback_passes_for_composition_absent_core(tmp_path):
    """A real causal GM1 thermodynamic table roundtrips without composition."""
    import shutil
    from compose_eos import read_compose
    from m7_batch_orchestrator import _native_pipeline
    source = ROOT / "eos_generator/eos_data/GM_GM1"
    assert run_native("verify", source, tmp_path / "source_verify.json")["overall"] == "PASS"
    core = tmp_path / "causal_core_without_composition"
    core.mkdir()
    # Preserve the physical thermodynamic bytes; omit only the composition sidecar.
    for name in ("eos.nb", "eos.t", "eos.yq", "eos.thermo"):
        shutil.copyfile(source / name, core / name)
    assert not read_compose(core).composition_file_present
    r = _native_pipeline("css", str(ROOT / "eos_generator/eos_data/SPG_M2"),
                         str(core), str(tmp_path / "pipeline"))
    assert r["readback_status"] == "PASS"
    assert r["combined_n"] == r["readback_n"]
    assert r["combined_n"] > 4
    assert r["c1_c11_overall"] == "PASS"
    assert r["c1_c11_failed_checks"] == ""


@pytest.mark.skipif(_needs_backend(), reason="requires the native Fortran backend")
def test_historical_apr_acausal_tail_is_rejected_before_readback(tmp_path):
    """The original APR fixture has an acausal retained tail; preserve rejection."""
    import json
    from m7_batch_orchestrator import BackendUnavailable, _native_pipeline
    core = ROOT / "eos_generator/eos_data/APR_APR"
    report_path = tmp_path / "apr_verify.json"
    summary = run_native("verify", core, report_path)
    assert summary["overall"] == "FAIL" and "C10" in summary["failed_checks"]
    report = json.loads(report_path.read_text())
    assert report["acausal_point_count"] > 0
    assert report["max_cs2"] > 1
    with pytest.raises(BackendUnavailable, match="exhaustive native endpoint search"):
        _native_pipeline("css", str(ROOT / "eos_generator/eos_data/SPG_M2"),
                         str(core), str(tmp_path / "rejected_pipeline"))
    assert not (tmp_path / "rejected_pipeline/internal/eos_data.csv").exists()


@pytest.mark.skipif(_needs_backend(), reason="requires the native Fortran backend")
def test_native_pipeline_readback_passes_for_composition_present_pair():
    """End-to-end: both composition-present (SPG_M2 + GM_GM1) must readback PASS."""
    from m7_batch_orchestrator import _native_pipeline
    with tempfile.TemporaryDirectory(prefix="e39007_") as td:
        r = _native_pipeline("css", str(ROOT / "eos_generator/eos_data/SPG_M2"),
                             str(ROOT / "eos_generator/eos_data/GM_GM1"), td)
        assert r["readback_status"] == "PASS"
        assert r["combined_n"] == r["readback_n"]
        assert r["combined_n"] > 4
        assert r["c1_c11_overall"] in ("PASS", "FAIL")
