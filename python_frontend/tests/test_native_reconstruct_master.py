from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import EOSData, NativeExecutionError, reconstruct_master


HBARC = 197.3269804
ME = 0.51099895
MMU = 105.6583745


def master_source(*, composition=False, outer_crust=False, dissolution=False,
                  nadd=0, q_extra=None, n=96):
    nb = np.geomspace(1.0e-7, 1.0e-2, n) if outer_crust else np.geomspace(0.05, 0.8, n)
    x = np.log(nb)
    mn = 939.0
    # Linear F(x) is reproduced exactly by the local cubic and PCHIP paths.
    free = 1000.0 + 30.0 * x
    q1 = np.full(n, 30.0)
    mu = free + q1
    zeros = np.zeros(n)
    phases = np.zeros(n, dtype=int)
    yn = np.zeros(n)
    yp = np.zeros(n)
    ye = np.zeros(n)
    ymu = np.zeros(n)
    abar = np.zeros(n)
    zbar = np.zeros(n)
    yh = np.zeros(n)

    if composition:
        yp[:] = 0.05
        yh[:] = 0.90
        abar[:] = 100.0
        zbar[:] = 40.0
        yn[:] = 0.05
        ye[:] = 0.41
        if dissolution:
            k = int(0.65 * n)
            # The last nuclear row has an inconsistent electron fraction.
            ye[k] = 0.10
            # Rows above k are homogeneous; the dissolution repair interpolates
            # the charge fraction from the previous crust row to the next core row.
            yp[k + 1:] = 0.15
            yh[k + 1:] = 0.0
            abar[k + 1:] = 0.0
            zbar[k + 1:] = 0.0
            yn[k + 1:] = 0.85
        if outer_crust:
            split = 20
            yp[:split] = 0.0
            yh[:split] = 1.0
            abar[:split] = 56.0
            zbar[:split] = 26.0
            yn[:split] = 0.0
            ye[:split] = 26.0 / 56.0
            yp[split:] = 0.10
            yh[split:] = 0.0
            abar[split:] = 0.0
            zbar[split:] = 0.0
            yn[split:] = 0.90
            ye[split:] = 0.10

    extra = zeros.copy() if q_extra is None else np.full(n, float(q_extra))
    return EOSData(
        nb=nb, mn=mn, mp=mn,
        Q1=q1, Q2=zeros.copy(), Q3=mu / mn - 1.0,
        Q4=np.full(n, -0.01), Q5=zeros.copy(),
        Q6=free / mn - 1.0, Q7=free / mn - 1.0,
        q_extra=extra, Iphase=phases,
        Yn=yn, Yp=yp, Ye=ye, Ymu=ymu,
        Aav=abar, Zav=zbar, Yh=yh,
        Il=1, iT=1, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=nadd,
        composition_file_present=composition, composition_available=composition,
        beta_equilibrated=True,
    )


def fermi_mu(nb, fraction, mass):
    kf = (3.0 * np.pi**2 * nb * fraction)**(1.0 / 3.0)
    return np.sqrt(mass**2 + (HBARC * kf)**2)


def test_native_master_reconstruction_obeys_exact_first_law_without_inventing_composition():
    source = master_source(composition=False)
    q4_before = source.Q4.copy()
    source_nb = source.nb.copy()

    out, report = reconstruct_master(
        source, smooth=True, smooth_sigma=0.05, return_report=True
    )

    assert report["status"] == "ok"
    assert report["composition_reconstructed"] is False
    assert report["composition_preserved"] is False
    assert report["composition_provenance"] == "source_composition_absent_and_not_invented"
    assert out.composition_available is False
    np.testing.assert_array_equal(out.Q4, q4_before)
    np.testing.assert_array_equal(out.Q2, np.zeros(out.n))
    np.testing.assert_array_equal(out.Q5, np.zeros(out.n))
    np.testing.assert_array_equal(out.Q6, out.Q7)
    np.testing.assert_allclose(out.F, source.F, rtol=3e-13, atol=3e-13)
    np.testing.assert_allclose(out.Q1, 30.0, rtol=1e-11, atol=3e-12)
    np.testing.assert_allclose(out.mub, out.F + out.Q1, rtol=0.0, atol=5e-13)
    np.testing.assert_allclose(out.p, out.nb * out.Q1, rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(out.eps, out.nb * out.F, rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(out.eps + out.p, out.nb * out.mub, rtol=2e-14)
    assert np.all(np.diff(out.nb) > 0.0)
    np.testing.assert_array_equal(source.nb, source_nb)


def test_native_master_unsmoothed_pchip_path_preserves_free_energy_knots():
    source = master_source(composition=False)
    source.Q7[:] = source.Q6

    out = reconstruct_master(source, smooth=False)

    np.testing.assert_array_equal(out.F, source.F)
    np.testing.assert_allclose(out.Q1, 30.0, rtol=2e-12, atol=2e-12)
    np.testing.assert_allclose(out.mub, out.F + out.Q1, rtol=2e-13)


def test_native_master_reconstruction_solves_electron_muon_and_charge_constraints():
    source = master_source(composition=True)
    out = reconstruct_master(source, smooth=True)

    target_charge = source.Yp + (source.Zav / source.Aav) * source.Yh
    np.testing.assert_allclose(out.Yn + out.Yp + out.Yh, 1.0, atol=2e-15)
    np.testing.assert_allclose(out.Ye + out.Ymu, target_charge, rtol=2e-12, atol=2e-13)
    assert np.all(out.Ye >= 0.0)
    assert np.all(out.Ymu >= 0.0)
    assert np.any(out.Ymu > 0.0)
    mu_e = np.array([fermi_mu(n, y, ME) for n, y in zip(out.nb, out.Ye)])
    mu_mu = np.array([fermi_mu(n, y, MMU) for n, y in zip(out.nb, out.Ymu)])
    active = out.Ymu > 1.0e-10
    np.testing.assert_allclose(mu_e[active], mu_mu[active], rtol=2e-11, atol=2e-10)
    np.testing.assert_allclose(out.Q4 * out.mn, -mu_e, rtol=2e-12, atol=2e-12)
    np.testing.assert_array_equal(out.Q5, np.zeros(out.n))
    assert out.composition_available is True
    assert out.composition_file_present is True


def test_native_master_dissolves_an_inconsistent_last_nuclear_row():
    source = master_source(composition=True, dissolution=True)
    k = int(0.65 * source.n)

    out = reconstruct_master(source, smooth=False, fix_dissolution=True, dissolution_tol=0.015)

    assert out.Aav[k] == 0.0
    assert out.Zav[k] == 0.0
    assert out.Yh[k] == 0.0
    assert out.Yp[k] == pytest.approx(0.5 * (0.4 * source.Yh[k - 1] + source.Yp[k + 1]), rel=2e-12)
    np.testing.assert_allclose(out.Ye[k] + out.Ymu[k], out.Yp[k], rtol=2e-12, atol=2e-13)


def test_native_master_outer_crust_uses_relativistic_electron_and_lattice_pressure():
    source = master_source(composition=True, outer_crust=True, n=96)
    out, report = reconstruct_master(source, smooth=False, return_report=True)

    nb = out.nb[0]
    ye = out.Ye[0]
    ne = nb * ye
    kf = (3.0 * np.pi**2 * ne)**(1.0 / 3.0)
    xe = HBARC * kf / ME
    if xe < 0.05:
        bracket = xe**5 * (1.6 + xe**2 * (-4.0 / 7.0 + xe**2 / 3.0))
    else:
        bracket = xe * (2.0 * xe**2 - 3.0) * np.sqrt(1.0 + xe**2) + 3.0 * np.arcsinh(xe)
    pe = ME**4 / (24.0 * np.pi**2 * HBARC**3) * bracket
    eps_lat = -0.9 * (4.0 / (9.0 * np.pi))**(1.0 / 3.0) * (1.0 / 137.035999084) * 26.0**(2.0 / 3.0) * HBARC * kf * ne
    expected_q1 = (pe + eps_lat / 3.0) / nb
    assert out.Q1[0] == pytest.approx(expected_q1, rel=3e-11, abs=1e-13)
    assert report["outer_crust_correction_applied"] is True
    np.testing.assert_allclose(out.eps + out.p, out.nb * out.mub, rtol=2e-13, atol=2e-13)


def test_native_master_preserves_extra_thermo_column_unless_explicitly_declared_cs2():
    source = master_source(composition=False, nadd=1, q_extra=0.37)
    out = reconstruct_master(source, smooth=False)
    np.testing.assert_array_equal(out.q_extra, source.q_extra)

    rebuilt, report = reconstruct_master(source, smooth=False, recompute_extra_cs2=True, return_report=True)
    expected = 30.0 / (rebuilt.F + 30.0)
    np.testing.assert_allclose(rebuilt.q_extra, expected, rtol=1e-10, atol=2e-12)
    assert report["q_extra_recomputed_as_cs2"] is True


def test_native_master_rejects_phase_pressure_plateau_instead_of_smoothing_across_it():
    source = master_source(composition=False)
    source.Iphase[source.n // 2:] = 1
    source.Q1[source.n // 2] = source.nb[source.n // 2 - 1] * source.Q1[source.n // 2 - 1] / source.nb[source.n // 2]
    with pytest.raises(NativeExecutionError, match="pressure plateau across an Iphase boundary"):
        reconstruct_master(source, smooth=True)


def test_native_master_rejects_invalid_smoothing_and_temperature_controls():
    source = master_source(composition=False)
    with pytest.raises(ValueError, match="smooth_sigma"):
        reconstruct_master(source, smooth_sigma=0.0)
    source.temperature_mev = 1.0
    with pytest.raises(NativeExecutionError, match="cold T=0"):
        reconstruct_master(source)


def test_native_master_rejects_nonzero_entropy_even_with_cold_metadata():
    source = master_source(composition=False)
    source.Q2[4] = 1.e-3
    with pytest.raises(NativeExecutionError, match="zero source entropy Q2"):
        reconstruct_master(source)


def test_native_master_rejects_trapped_lepton_potential_for_neutrino_transparent_rebuild():
    source = master_source(composition=False)
    source.Q5[4] = 1.e-3
    with pytest.raises(NativeExecutionError, match="neutrino-transparent beta equilibrium"):
        reconstruct_master(source)

\n\ndef test_native_master_rejects_cold_source_with_inconsistent_free_and_internal_energy():
    source = master_source(composition=False)
    source.Q7[5] += 1.0e-4
    with pytest.raises(NativeExecutionError, match="cold-source Q6=Q7"):
        reconstruct_master(source, smooth=False)


def test_native_master_rejects_source_baryon_chemical_potential_euler_mismatch():
    source = master_source(composition=False)
    source.Q3[5] += 1.0e-3
    with pytest.raises(NativeExecutionError, match="Euler chemical potential identity"):
        reconstruct_master(source, smooth=False)


def test_native_master_preserves_compose_il_zero_charge_potential_convention():
    source = master_source(composition=False)
    source.Il = 0
    source.charge_fraction_grid = 0.25
    source.Q4[:] = -1.0e-2
    mu_source = source.F + source.Q1 - source.charge_fraction_grid * source.Q4 * source.mn
    source.Q3[:] = mu_source / source.mn - 1.0

    out = reconstruct_master(source, smooth=False)

    expected_q3 = (out.F + out.Q1) / out.mn - 1.0 - source.charge_fraction_grid * out.Q4
    np.testing.assert_allclose(out.Q3, expected_q3, rtol=0.0, atol=2e-14)
    np.testing.assert_allclose(out.mub, mu_source, rtol=2e-14, atol=2e-12)
