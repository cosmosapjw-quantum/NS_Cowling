import numpy as np
import pytest

from mode_identity.diagnostics import (
    ModeProfile, PhaseProfile, diagnose_mode, phase_profile_overlap, track_spectra,
)


def make_profile(*, scale=1.0, sign=1.0, mode="f", freq=1800.0, model="M0", conversion="slow"):
    rc, re = np.linspace(0.2, 1.0, 101), np.linspace(1.0, 2.0, 121)
    if mode == "f":
        wc, we = rc**3, 1.0 + 0.1 * (re - 1.0)
    elif mode == "p1":
        wc, we = rc**3 * (rc - 0.65), 1.0 + 0.2 * np.sin(np.pi * (re - 1.0))
    elif mode == "core":
        wc = rc**3 * (1.0 + 0.1 * np.sin(2.0 * np.pi * rc))
        we = 0.03 * (1.0 + 0.1 * (re - 1.0))
    else:
        raise ValueError(mode)
    core = PhaseProfile(rc, sign * scale * wc, sign * scale * 0.2 * rc**2,
                        np.full(rc.size, 1e35), np.full(rc.size, -0.2), np.full(rc.size, 0.1))
    env = PhaseProfile(re, sign * scale * we, sign * scale * 0.15 * re,
                       np.full(re.size, 0.8e35), np.full(re.size, -0.15), np.full(re.size, 0.08))
    return ModeProfile(core, env, freq, ell=2, conversion=conversion, model_id=model)


def test_norm_and_trapping_are_positive_and_normalized():
    d = diagnose_mode(make_profile())
    assert d["norm"] > 0 and d["norm_core"] > 0 and d["norm_envelope"] > 0
    assert d["core_trapping_fraction"] + d["envelope_trapping_fraction"] == pytest.approx(1.0)


def test_interface_polarization_has_physical_radius_factor():
    d = diagnose_mode(make_profile())
    assert d["interface_Pi_inner"]["Pi"] == pytest.approx(np.sqrt(6.0) * 0.2)


def test_norm_is_quadratic_and_trapping_is_scale_invariant():
    a, b = diagnose_mode(make_profile()), diagnose_mode(make_profile(scale=1e-9))
    assert b["norm"] / a["norm"] == pytest.approx(1e-18, rel=2e-8)
    assert b["core_trapping_fraction"] == pytest.approx(a["core_trapping_fraction"], rel=1e-12)


def test_overlap_is_one_under_global_sign_and_scale_change():
    a, b = make_profile(), make_profile(scale=7.0, sign=-1.0, freq=1900.0)
    ov = phase_profile_overlap(a, b)
    assert ov["overall"] == pytest.approx(1.0, abs=1e-12)
    assert ov["core"] == pytest.approx(1.0, abs=1e-12)
    assert ov["envelope"] == pytest.approx(1.0, abs=1e-12)


def test_global_overlap_retains_relative_phase_flip_across_interface():
    a = make_profile()
    e = a.envelope
    env = PhaseProfile(e.r_m, -e.W_m3, -e.V_m2, e.enthalpy_j_m3, e.nu, e.lambda_)
    b = ModeProfile(a.core, env, a.frequency_hz, a.ell, a.conversion, "flip")
    ov = phase_profile_overlap(a, b)
    assert ov["core"] == pytest.approx(1.0, abs=1e-12)
    assert ov["envelope"] == pytest.approx(1.0, abs=1e-12)
    assert ov["overall"] < 0.95


def test_overlap_symmetry_and_phase_resolution():
    a, b = make_profile(), make_profile(mode="core", freq=1900.0)
    ab, ba = phase_profile_overlap(a, b), phase_profile_overlap(b, a)
    assert ab["overall"] == pytest.approx(ba["overall"], rel=1e-12, abs=1e-12)
    assert ab["core"] != pytest.approx(ab["envelope"], abs=1e-4)


def test_complex_profiles_do_not_get_real_sign_reversal_label():
    p = make_profile()
    c = PhaseProfile(p.core.r_m, p.core.W_m3.astype(complex)*(1+0.2j),
                     p.core.V_m2.astype(complex)*(1+0.2j), p.core.enthalpy_j_m3, p.core.nu, p.core.lambda_)
    e = PhaseProfile(p.envelope.r_m, p.envelope.W_m3.astype(complex)*(1-0.3j),
                     p.envelope.V_m2.astype(complex)*(1-0.3j), p.envelope.enthalpy_j_m3, p.envelope.nu, p.envelope.lambda_)
    assert diagnose_mode(ModeProfile(c, e, p.frequency_hz))["interface_W_sign_reversal"] is None


def test_track_spectra_keeps_branch_id_separate_from_character():
    prev = [make_profile(mode="f", freq=1800.0), make_profile(mode="p1", freq=4900.0)]
    curr = [make_profile(mode="p1", freq=4950.0), make_profile(mode="f", freq=1820.0)]
    out = track_spectra(prev, curr, ["B_f", "B_p1"], min_overlap=0.5)
    by_id = {m["spectral_branch_id"]: m["current_index"] for m in out["matches"]}
    assert by_id["B_f"] == 1 and by_id["B_p1"] == 0
    assert all(m["physical_character"] == "unresolved" for m in out["matches"])


def test_ambiguous_duplicate_profiles_are_not_forced_into_ids():
    p = [make_profile(freq=1800.0), make_profile(freq=1810.0)]
    c = [make_profile(freq=1805.0), make_profile(freq=1815.0)]
    out = track_spectra(p, c, ["B0", "B1"], ambiguity_margin=0.05)
    assert all(m["spectral_branch_id"] is None for m in out["matches"])


def test_slow_and_rapid_continuations_are_separate():
    with pytest.raises(ValueError):
        phase_profile_overlap(make_profile(), make_profile(conversion="rapid"))


def test_nonmatching_interface_geometry_is_rejected():
    p = make_profile()
    e = p.envelope
    bad = PhaseProfile(e.r_m + 0.1, e.W_m3, e.V_m2, e.enthalpy_j_m3, e.nu, e.lambda_)
    with pytest.raises(ValueError):
        diagnose_mode(ModeProfile(p.core, bad, p.frequency_hz))
