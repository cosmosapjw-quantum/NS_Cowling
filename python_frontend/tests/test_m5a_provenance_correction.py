"""§6 regression: M5A provenance correction. All 7 audit findings fixed.

Verifies the corrected M5A receipt truthfully reports the verification
backend, uses the real BBB_BHF2 proton mass (not the invented M_N*0.9938),
quantifies the first-law reset, records the mu_t mismatch, explicitly
classifies Delta_eps=0 as CONTINUOUS_UNCLASSIFIED, and declares the output as a
barotrope (not a CompOSE-generated EOS). The misnamed GPPVA_DD2 file is gone.
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "maxwell_audit"))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_ROOT", str(ROOT))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_EXE",
                      str(ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"))

RECEIPT = ROOT / "download" / "m5a_real_css" / "E390-05A_M5A_REAL_CSS_RECEIPT_20261114.json"
OLD_MISNAMED = ROOT / "download" / "m5a_real_css" / "GPPVA_DD2_m5a_real_css.json"
NEW_NAMED = ROOT / "download" / "m5a_real_css" / "BBB_BHF2_m5a_real_css.json"

M_N = 939.56536


@pytest.fixture(scope="module")
def rep():
    return json.loads(RECEIPT.read_text())


def test_m5a_receipt_exists(rep):
    assert RECEIPT.exists()


def test_m5a_finding1_python_check_not_native(rep):
    """Finding 1: the consistency check field is renamed to
    python_consistency_check_pass (NOT the misleading native_c1_c11_overall_pass).
    The native_c1_c11_verification and native_fortran_compose_readback fields
    are TRUE only because the native Fortran read/write/verify pathway was
    actually executed (scripts/maxwell_audit/m5a_native_readback_verify.py);
    they carry concrete evidence in native_readback_evidence. The §6 spec
    allows these to be TRUE when actual native calls succeed."""
    sample = rep["results"]["nB0.4_cs20.3333_dEps0.0"]
    # the misleading field name is gone
    assert "python_consistency_check_pass" in sample
    assert "native_c1_c11_overall_pass" not in sample
    # the native verification fields are TRUE with evidence (not BLOCKED)
    assert sample["native_c1_c11_verification"] is True
    assert sample["native_fortran_compose_readback"] is True
    # the truthful verification-backend summary
    vbt = rep["verification_backend_truthful"]
    assert vbt["python_consistency_check"] is True
    assert vbt["native_c1_c11_verification"] is True
    assert vbt["native_fortran_compose_readback"] is True
    # concrete evidence must be present
    ev = rep["native_readback_evidence"]
    assert ev["native_write_executed"] is True
    assert ev["native_read_executed"] is True
    assert ev["native_verify_executed"] is True
    assert ev["readback_n_match"] is True
    assert ev["c1_c11_overall_pass"] is True
    # failed_checks is empty (either "" or []) -- no checks failed
    assert ev["c1_c11_failed_checks"] in ("", [], None)


def test_m5a_finding2_filename_corrected(rep):
    """Finding 2: the misnamed GPPVA_DD2 file is removed; the correctly-named
    BBB_BHF2 file exists."""
    assert not OLD_MISNAMED.exists(), "misnamed GPPVA_DD2_m5a_real_css.json must be removed"
    assert NEW_NAMED.exists(), "BBB_BHF2_m5a_real_css.json must exist"
    assert rep["output_barotrope_file"] == "BBB_BHF2_m5a_real_css.json"


def test_m5a_finding3_real_proton_mass(rep):
    """Finding 3: the EOSData uses the real BBB_BHF2 mp = 938.272 MeV, not the
    invented M_N*0.9938 = 933.74 MeV."""
    real_mp = rep["hadronic_source_real_mp_mev"]
    invented = rep["hadronic_source_invented_mp_was_M_N_times_0_9938_mev"]
    assert real_mp == pytest.approx(938.272, rel=1e-6)
    assert invented == pytest.approx(M_N * 0.9938, rel=1e-9)
    assert abs(invented - real_mp) > 4.0  # the 4.53 MeV misidentification


def test_m5a_finding4_first_law_reset_quantified(rep):
    """Finding 4: the first-law reset eps := mu*n - P is quantified. The max
    relative change vs the source energy density is recorded (1.77e-8)."""
    reset = rep["first_law_reset_max_relative_change"]
    assert reset > 0, "first-law reset must be quantified (nonzero)"
    assert reset < 1e-6, "first-law reset must be small (< 1e-6)"


def test_m5a_finding5_mu_t_mismatch_recorded(rep):
    """Finding 5: the mu_t mismatch (interpolated source mu_B vs
    (eps_l+P_t)/n_l) is recorded per transition density, along with P and eps
    interpolation residuals."""
    for key in ["nB0.3_cs20.3333_dEps0.0", "nB0.4_cs20.3333_dEps0.0", "nB0.5_cs20.3333_dEps0.0"]:
        r = rep["results"][key]
        assert "mu_t_mismatch_relative" in r
        assert "P_interp_residual_relative" in r
        assert "eps_interp_residual_relative" in r
        assert r["mu_t_mismatch_relative"] >= 0


def test_m5a_finding6_delta_eps_classification_explicit(rep):
    """Finding 6: Delta_eps=0 is explicitly CONTINUOUS_UNCLASSIFIED (zero jump
    alone is NOT a crossover/2nd-order inference); Delta_eps>0 is
    SHARP_FIRST_ORDER."""
    cont = rep["results"]["nB0.4_cs20.3333_dEps0.0"]
    sharp = rep["results"]["nB0.4_cs20.3333_dEps100.0"]
    assert cont["delta_eps_classification_explicit"] == "CONTINUOUS_UNCLASSIFIED"
    assert sharp["delta_eps_classification_explicit"] == "SHARP_FIRST_ORDER"


def test_m5a_finding7_barotrope_not_compose(rep):
    """Finding 7: the saved JSON is declared a barotrope (n_B, mu_B, P, eps),
    NOT a CompOSE-generated EOS. The native Fortran read-back may be run as a
    verification step (and recorded in native_readback_evidence), but that
    does NOT promote the barotrope to a CompOSE-generated multi-phase EOS --
    the headline artifact is still the (n_B, mu_B, P, eps) barotrope."""
    assert rep["output_is_compose_generated_eos"] is False
    assert rep["output_is_barotrope_only"] is True
    # the native read-back is a verification backend, not a CompOSE-EOS claim:
    # it is TRUE only because the native I/O pathway was actually executed
    # (with evidence), and output_is_compose_generated_eos stays False.
    assert rep["verification_backend_truthful"]["native_fortran_compose_readback"] is True
    assert rep["verification_backend_truthful"]["native_c1_c11_verification"] is True
    # but the artifact is still declared a barotrope, not a CompOSE EOS
    assert rep["output_is_compose_generated_eos"] is False
