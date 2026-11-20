"""C7 physics checks: source barotrope closure, never frozen-composition inference."""
from __future__ import annotations

import io
import json
import zipfile
from pathlib import Path

import pytest

from compose_eos import stellar_sharp
from compose_eos.whole_crust import WholeCrustBarotrope

ROOT = Path(__file__).resolve().parents[2]
ZIP = ROOT / "data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
H_POINTS = (0.005, 0.025, 0.055, 0.09, 0.14, 0.18)


def _eos(model, method):
    with zipfile.ZipFile(ZIP) as bundle:
        with zipfile.ZipFile(io.BytesIO(bundle.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inner:
            blob = inner.read(f"FRESH_5x2_WHOLE_STARS/{model}/{method}/REPRESENTATION.json")
    return stellar_sharp.SharpContractEnthalpyEOS(
        WholeCrustBarotrope.from_contract(json.loads(blob)))


@pytest.mark.parametrize("model,join", (("SPG_M2", "maxwell"), ("SPG_M3", "css")))
def test_real_source_equilibrium_derivatives_and_latent_heat(model, join):
    eos = _eos(model, join)
    result = stellar_sharp.audit_equilibrium_cowling_source(eos, H_POINTS)
    assert result["status"] == "PASS_SCOPED_EQUILIBRIUM_BAROTROPIC"
    assert result["perturbation_closure"] == "GAMMA1_EQUALS_GAMMA_EQ"
    assert result["frozen_composition_gamma1"] == "NOT_DERIVABLE_FROM_EQUILIBRIUM_BAROTROPE"
    assert result["composition_g_mode_admitted"] is False
    assert len(result["smooth_samples"]) == len(H_POINTS)
    assert len(result["interfaces"]) == 1
    assert result["max_rel_dPdh_error"] < 5e-5
    assert result["max_rel_cs2_error"] < 5e-5
    assert result["interfaces"][0]["delta_epsilon_MeV_fm3"] > 0
    assert result["interfaces"][0]["max_rel_coexistence_error"] < 1e-9
    assert min(row["cs2_over_c2"] for row in result["smooth_samples"]) > 0


def test_cs2_corruption_fails_physics_audit_not_a_fast_mode_pass(monkeypatch):
    eos = _eos("SPG_M3", "css")
    original = eos.state_at_enthalpy

    def poisoned(h, side="high"):
        state = original(h, side)
        if abs(h - 0.055) < 1e-13:
            state["cs2_over_c2"] *= 1.1
        return state

    monkeypatch.setattr(eos, "state_at_enthalpy", poisoned)
    result = stellar_sharp.audit_equilibrium_cowling_source(eos, (0.055,))
    assert result["status"] == "FAIL_SOURCE_THERMODYNAMIC_INCONSISTENCY"
    assert result["max_rel_cs2_error"] > 0.09


def test_sampling_exact_knot_is_rejected_not_smoothed(monkeypatch):
    eos = _eos("SPG_M2", "maxwell")
    with pytest.raises(ValueError, match="smooth-interval"):
        stellar_sharp.audit_equilibrium_cowling_source(eos, (eos.transitions[0]["h"],))


def test_adiabatic_nonidentifiability_counterexample():
    # epsilon(n,x)=mn + A*n^2/2 + K*(x-a*n/n0)^2/2.
    # Every K has the SAME equilibrium P(epsilon); frozen derivative differs.
    n, m, a, A, n0 = 0.16, 900.0, 0.03, 500.0, 0.16
    eq = A*n/(m+A*n)
    frozen_zero = (A*n + 0*n*(a/n0)**2)/(m+A*n)
    frozen_nonzero = (A*n + 100*n*(a/n0)**2)/(m+A*n)
    assert eq == frozen_zero
    assert frozen_nonzero > eq

@pytest.mark.slow
@pytest.mark.parametrize("model,join,nc,conversion,bracket", (
    ("SPG_M2", "maxwell", 0.30385128926615257, "rapid", (5110.0, 5170.0)),
    ("SPG_M3", "css", 0.192510767556007, "slow", (3710.0, 3760.0)),
))
def test_real_native_mode_reports_exact_thermodynamic_closure(
        model, join, nc, conversion, bracket):
    from compose_eos.stellar_sharp_native import (
        solve_sharp_background_native, NativeSharpCowlingDVODE)
    from compose_eos.stellar_cowling import CowlingSettings

    eos = _eos(model, join)
    bg = solve_sharp_background_native(
        eos, nc, native_samples=2049, surface_h=3e-8,
        rtol=3e-10, atol=3e-11, max_step_h=.004)
    settings = CowlingSettings(surface_h=3e-8, rtol=1e-9,
                               atol=1e-12, max_step_h=.002)
    mode = NativeSharpCowlingDVODE(
        bg, settings, samples=4097, coefficient_backend="hybrid_cs2").solve_mode(
            bracket, conversion=conversion, profile_points=350)
    assert mode["relative_join_mismatch"] < 1e-5
    assert mode["perturbation_closure"] == "GAMMA1_EQUALS_GAMMA_EQ"
    assert mode["frozen_composition_gamma1"] == "NOT_DERIVABLE_FROM_EQUILIBRIUM_BAROTROPE"
    assert mode["bulk_composition_g_modes_supported"] is False
