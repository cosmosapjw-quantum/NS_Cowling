"""Adversarial native Maxwell root isolation regressions."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import EOSData, construct_maxwell_transition, read_compose, write_compose


def _dyadic_quadratic_phase(rest_mass: float, stiffness: float) -> EOSData:
    """P=K n² on exactly representable knots and CompOSE normalization.

    Both phases have mu knots 1025.2578125 + j/64 MeV, j=0,...,128.
    Their unique coexistence point in the overlap is mu=1026 MeV,
    exactly halfway between two common knots. Binary arithmetic therefore
    encounters an exactly zero pressure difference during bisection.

    The spacing also resolves the independent verifier's log-density cutoff;
    a five-knot version exercises the same root but masks every finite-
    difference stencil and therefore cannot certify C10/C11 through the API.
    """
    mu = np.arange(129, dtype=float) / 64.0 + 1025.2578125
    nb = (mu - rest_mass) / (2.0 * stiffness)
    zeros = np.zeros_like(nb)
    q7 = (rest_mass + stiffness * nb) / 1024.0 - 1.0
    return EOSData(
        nb=nb, mn=1024.0, mp=1024.0,
        Q1=stiffness * nb, Q2=zeros.copy(), Q3=mu / 1024.0 - 1.0,
        Q4=zeros.copy(), Q5=zeros.copy(), Q6=q7.copy(), Q7=q7.copy(),
        q_extra=zeros.copy(), Iphase=np.zeros(nb.size, dtype=int),
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        Il=1, iT=1, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=True,
    )


def test_native_maxwell_accepts_exact_zero_bisection_midpoint():
    low = _dyadic_quadratic_phase(rest_mass=1024.0, stiffness=1.0)
    high = _dyadic_quadratic_phase(rest_mass=1025.0, stiffness=0.25)
    out, report = construct_maxwell_transition(
        low, high, assume_local_charge_neutrality=True, return_report=True
    )
    assert report["status"] == "ok"
    assert report["pressure_crossings"] == 1
    assert report["mu_transition"] == 1026.0
    assert report["pressure_transition"] == pytest.approx(1.0, abs=1e-13)
    assert report["n_low"] == pytest.approx(1.0, abs=1e-13)
    assert report["n_high"] == pytest.approx(2.0, abs=1e-13)
    assert report["latent_energy_density"] == pytest.approx(1026.0, abs=1e-10)
    assert np.all(np.diff(out.nb) > 0.0)


def test_native_maxwell_phase_only_sidecar_preserves_readback(tmp_path):
    low = _dyadic_quadratic_phase(rest_mass=1024.0, stiffness=1.0)
    high = _dyadic_quadratic_phase(rest_mass=1025.0, stiffness=0.25)
    out = construct_maxwell_transition(low, high, assume_local_charge_neutrality=True)
    assert out.composition_file_present is True
    assert out.composition_available is False
    assert out.nT == out.nYq == 1
    target = tmp_path / "maxwell"
    write_compose(target, out)
    restored = read_compose(target)
    assert restored.composition_file_present is True
    assert restored.composition_available is False
    for field in ("nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "q_extra"):
        np.testing.assert_allclose(getattr(restored, field), getattr(out, field), rtol=5e-13, atol=1e-14)
    for field in ("Iphase", "row_status", "Yn_present", "Yp_present", "Ye_present",
                  "Ymu_present", "Abar_present", "Zbar_present", "Yh_present",
                  "C3_applicable", "C4_applicable"):
        np.testing.assert_array_equal(getattr(restored, field), getattr(out, field))
    assert np.all(restored.row_status == 1)  # phase-only row; abundances absent
    for field in ("Yn_present", "Yp_present", "Ye_present", "Ymu_present",
                  "Abar_present", "Zbar_present", "Yh_present", "C3_applicable", "C4_applicable"):
        assert not np.any(getattr(restored, field))
    for line, phase in zip((target / "eos.compo").read_text().splitlines(), out.Iphase):
        assert line.split()[3:] == [str(int(phase)), "0", "0"]


def test_native_maxwell_pipeline_survivor_keeps_phase_only_sidecar(tmp_path):
    """One admitted synthetic EOS pair exercises the real census I/O path."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from match_core_crust import run_pipeline

    low_path, high_path = tmp_path / "low", tmp_path / "high"
    write_compose(low_path, _dyadic_quadratic_phase(rest_mass=1024.0, stiffness=1.0))
    write_compose(high_path, _dyadic_quadratic_phase(rest_mass=1025.0, stiffness=0.25))
    result = run_pipeline("maxwell", low_path, high_path, tmp_path / "result")
    assert result["survived"] is True, result
    assert result["readback_validation"]["status"] == "PASS"
    assert result["admission"]["status"] == "ADMITTED"
    assert (tmp_path / "result/compose/eos.compo").is_file()
