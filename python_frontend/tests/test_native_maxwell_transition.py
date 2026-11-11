from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import EOSData, NativeExecutionError, construct_maxwell_transition, check_consistency


def quadratic_phase(*, rest_mass: float, stiffness: float, n_min: float,
                    n_max: float, count: int = 81) -> EOSData:
    """Exact cold barotrope P=K n^2, eps=m n+K n^2."""
    nb = np.linspace(n_min, n_max, count)
    pressure = stiffness * nb**2
    epsilon = rest_mass * nb + pressure
    mu = (epsilon + pressure) / nb
    zeros = np.zeros_like(nb)
    q7 = epsilon / (nb * 939.0) - 1.0
    return EOSData(
        nb=nb, mn=939.0, mp=939.0,
        Q1=pressure / nb, Q2=zeros.copy(), Q3=mu / 939.0 - 1.0,
        Q4=zeros.copy(), Q5=zeros.copy(), Q6=q7.copy(), Q7=q7.copy(),
        q_extra=zeros.copy(), Iphase=np.zeros(nb.size, dtype=int),
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        Il=1, iT=1, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=True,
    )


def acausal_cubic_phase() -> EOSData:
    """P=K n^3 and eps=m n+K n^3/2; c_s^2 exceeds unity at high n."""
    nb = np.linspace(0.5, 2.0, 81)
    stiffness = 1000.0
    rest_mass = 940.0
    pressure = stiffness * nb**3
    epsilon = rest_mass * nb + 0.5 * stiffness * nb**3
    mu = (epsilon + pressure) / nb
    zeros = np.zeros_like(nb)
    q7 = epsilon / (nb * 939.0) - 1.0
    return EOSData(
        nb=nb, mn=939.0, mp=939.0,
        Q1=pressure / nb, Q2=zeros.copy(), Q3=mu / 939.0 - 1.0,
        Q4=zeros.copy(), Q5=zeros.copy(), Q6=q7.copy(), Q7=q7.copy(),
        q_extra=zeros.copy(), Iphase=np.zeros(nb.size, dtype=int),
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        Il=1, iT=1, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=True,
    )


def test_compose_phase_only_sidecar_roundtrips_Iphase_without_fabricating_composition(tmp_path):
    from compose_eos import read_compose, write_compose

    data = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    data.Iphase[:] = np.where(data.nb < 0.85, 4, 3)
    data.composition_file_present = True
    data.composition_available = False

    target = tmp_path / "phase_only_compose"
    write_compose(target, data)
    restored = read_compose(target)

    assert restored.composition_file_present is True
    assert restored.composition_available is False
    np.testing.assert_array_equal(restored.Iphase, data.Iphase)
    for field in (restored.Yn, restored.Yp, restored.Ye, restored.Ymu,
                  restored.Aav, restored.Zav, restored.Yh):
        assert np.all(field == 0.0)
    rows = (target / "eos.compo").read_text(encoding="utf-8").splitlines()
    assert len(rows) == data.n
    for i, row in enumerate(rows):
        fields = row.split()
        assert len(fields) == 6
        assert int(fields[3]) == int(data.Iphase[i])
        assert fields[4:] == ["0", "0"]


def test_native_maxwell_transition_recovers_exact_analytic_coexistence():
    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = quadratic_phase(rest_mass=940.0, stiffness=0.25, n_min=1.3, n_max=2.8)

    out, report = construct_maxwell_transition(low, high, assume_local_charge_neutrality=True, return_report=True)

    assert report["status"] == "ok"
    assert report["pressure_crossings"] == 1
    assert report["is_exact_maxwell"] is True
    assert report["tov_interface_support"] == "not_verified"
    assert report["mu_transition"] == pytest.approx(941.0, rel=2e-10)
    assert report["pressure_transition"] == pytest.approx(1.0, rel=2e-10)
    assert report["n_low"] == pytest.approx(1.0, rel=2e-10)
    assert report["n_high"] == pytest.approx(2.0, rel=2e-10)
    assert report["epsilon_low"] == pytest.approx(940.0, rel=2e-10)
    assert report["epsilon_high"] == pytest.approx(1881.0, rel=2e-10)
    assert report["density_jump"] == pytest.approx(1.0, rel=2e-10)
    assert report["latent_energy_density"] == pytest.approx(941.0, rel=2e-10)
    assert report["maxwell_identity_residual"] < 1e-11

    assert np.all(np.diff(out.nb) > 0.0)
    assert np.all(np.diff(out.p) >= -1e-12)
    assert np.all(np.diff(out.eps) > 0.0)
    assert out.composition_file_present is False
    assert out.composition_available is False
    assert np.all(out.Q2 == 0.0)
    assert np.all(out.Q5 == 0.0)
    assert np.all(out.Q6 == out.Q7)
    assert np.count_nonzero(out.Iphase == 1) == report["low_phase_rows"]
    assert np.count_nonzero(out.Iphase == 2) == report["high_phase_rows"]
    boundary = np.flatnonzero(np.diff(out.Iphase) != 0)
    assert boundary.size == 1
    i = int(boundary[0])
    assert out.p[i] == pytest.approx(out.p[i + 1], rel=2e-11)
    assert out.mub[i] == pytest.approx(out.mub[i + 1], rel=2e-11)
    assert out.nb[i + 1] > out.nb[i]
    np.testing.assert_allclose(out.eps[i + 1] - out.eps[i],
                               out.mub[i] * (out.nb[i + 1] - out.nb[i]),
                               rtol=2e-11, atol=2e-11)

    verification = check_consistency(out)
    assert verification["overall"]["pass_"] is True
    assert verification["C6"]["pass_"] is True
    assert verification["C10"]["pass_"] is True
    assert verification["C11"]["pass_"] is True


def test_native_maxwell_rejects_composition_instead_of_silently_dropping_it():
    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = quadratic_phase(rest_mass=940.0, stiffness=0.25, n_min=1.3, n_max=2.8)
    low.composition_file_present = True
    with pytest.raises(NativeExecutionError, match="composition-bearing inputs are rejected"):
        construct_maxwell_transition(low, high, assume_local_charge_neutrality=True)


def test_native_maxwell_rejects_nonoverlapping_mu_domains():
    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = quadratic_phase(rest_mass=950.0, stiffness=0.25, n_min=1.3, n_max=2.8)
    with pytest.raises(NativeExecutionError, match="no overlapping baryon-chemical-potential interval"):
        construct_maxwell_transition(low, high, assume_local_charge_neutrality=True)


def test_native_maxwell_rejects_acausal_phase_before_finding_a_transition():
    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = acausal_cubic_phase()
    with pytest.raises(NativeExecutionError, match="violates the causal bound c_s\\^2/c\\^2<=1"):
        construct_maxwell_transition(low, high, assume_local_charge_neutrality=True)


def test_native_maxwell_rejects_provider_first_law_mismatch():
    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = quadratic_phase(rest_mass=940.0, stiffness=0.25, n_min=1.3, n_max=2.8)
    high.Q3[10] += 1.e-2
    with pytest.raises(NativeExecutionError, match="provider Q3 chemical potential disagrees with the Euler value"):
        construct_maxwell_transition(low, high, assume_local_charge_neutrality=True)


def test_maxwell_requires_explicit_local_neutrality_assumption():
    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = quadratic_phase(rest_mass=940.0, stiffness=0.25, n_min=1.3, n_max=2.8)
    with pytest.raises(ValueError, match="confirm that both input phases are locally charge neutral"):
        construct_maxwell_transition(low, high)


def test_maxwell_phase_barotrope_cannot_be_written_as_complete_compose_without_composition(tmp_path):
    from compose_eos import NativeFeatureUnavailable, write_compose

    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = quadratic_phase(rest_mass=940.0, stiffness=0.25, n_min=1.3, n_max=2.8)
    out = construct_maxwell_transition(
        low, high, assume_local_charge_neutrality=True
    )
    with pytest.raises(NativeFeatureUnavailable, match="phase composition"):
        write_compose(tmp_path / "bad-compose-export", out)


def test_native_cli_writer_rejects_phase_barotrope_without_composition(tmp_path):
    import subprocess
    from compose_eos import executable

    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = quadratic_phase(rest_mass=940.0, stiffness=0.25, n_min=1.3, n_max=2.8)
    out = construct_maxwell_transition(
        low, high, assume_local_charge_neutrality=True
    )
    internal = tmp_path / "internal"
    destination = tmp_path / "compose"
    internal.mkdir()
    destination.mkdir()
    out.write_internal_directory(internal)
    result = subprocess.run(
        [str(executable()), "write", str(internal), str(destination)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0
    assert "cannot export a phase-discontinuous barotrope" in result.stdout
    assert list(destination.iterdir()) == []


def test_native_maxwell_rejects_mismatched_charge_fraction_metadata():
    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = quadratic_phase(rest_mass=940.0, stiffness=0.25, n_min=1.3, n_max=2.8)
    high.charge_fraction_grid = 0.1
    with pytest.raises(NativeExecutionError, match="same CompOSE lepton-content and charge-fraction conventions"):
        construct_maxwell_transition(low, high, assume_local_charge_neutrality=True)


def test_native_maxwell_rejects_provider_specific_nadd():
    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = quadratic_phase(rest_mass=940.0, stiffness=0.25, n_min=1.3, n_max=2.8)
    high.q_extra_count = 1
    with pytest.raises(NativeExecutionError, match="provider-specific Nadd"):
        construct_maxwell_transition(low, high, assume_local_charge_neutrality=True)

def test_native_maxwell_respects_il0_compose_q3_charge_potential_convention():
    """I_l=0 Q3 is Euler mu minus Yq*muq, including across the jump."""
    yq = 0.1
    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = quadratic_phase(rest_mass=940.0, stiffness=0.25, n_min=1.3, n_max=2.8)
    for phase in (low, high):
        phase.Il = 0
        phase.charge_fraction_grid = yq
        phase.Q4[:] = 0.2
        phase.Q3[:] -= yq * phase.Q4

    out, report = construct_maxwell_transition(
        low, high, assume_local_charge_neutrality=True, return_report=True
    )

    assert report["status"] == "ok"
    mu_euler = (out.eps + out.p) / out.nb
    expected_q3 = (mu_euler / out.mn) - 1.0 - yq * out.Q4
    np.testing.assert_allclose(out.Q3, expected_q3, rtol=2e-13, atol=2e-14)
    np.testing.assert_allclose(out.mub, mu_euler - yq * out.Q4 * out.mn,
                               rtol=2e-13, atol=2e-12)
    verification = check_consistency(out)
    assert verification["C2"]["pass_"] is True
    assert verification["overall"]["pass_"] is True

