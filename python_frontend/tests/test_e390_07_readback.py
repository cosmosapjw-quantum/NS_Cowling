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
def test_native_pipeline_readback_passes_for_composition_absent_core():
    """End-to-end: a composition-absent core (APR_APR, no eos.compo) spliced
    with a composition-present crust (SPG_M2) via match-core must readback
    PASS with combined_n == readback_n."""
    from m7_batch_orchestrator import _native_pipeline
    with tempfile.TemporaryDirectory(prefix="e39007_") as td:
        r = _native_pipeline("css", str(ROOT / "eos_generator/eos_data/SPG_M2"),
                             str(ROOT / "eos_generator/eos_data/APR_APR"), td)
        assert r["readback_status"] == "PASS"
        assert r["combined_n"] == r["readback_n"]
        assert r["combined_n"] > 4
        # C1-C11 overall is recorded (PASS or FAIL) but does not block readback
        assert r["c1_c11_overall"] in ("PASS", "FAIL")
        assert isinstance(r["c1_c11_failed_checks"], str)


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
