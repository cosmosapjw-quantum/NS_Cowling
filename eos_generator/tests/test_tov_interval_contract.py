"""Counterexamples to pressure/density mapping and interval-domain admission."""
from dataclasses import replace
from types import SimpleNamespace
import numpy as np
import pytest

from compose_eos.tov import (
    BarotropicEOS, TOVConfig, _certify_rhs_causality,
    find_causal_domain, solve_tov_star, MEV_FM3_TO_PA, C,
)


def _table():
    eps = np.geomspace(200.0, 2000.0, 12)
    p = 0.9 * (eps - 150.0)
    p[6] *= 0.97
    h = eps + p
    nb = 0.16 * np.exp(np.cumsum(np.r_[0.0, np.diff(eps) / ((h[1:] + h[:-1]) / 2.0)]))
    return BarotropicEOS(nb, p, eps)


def test_causal_density_limit_maps_to_the_admitted_pressure():
    eos = _table()
    domain = find_causal_domain(eos)
    actual_pc = float(eos.pressure_of_nb(domain.nb_max))
    assert actual_pc == pytest.approx(domain.pressure_max, rel=2e-13)
    assert float(eos.cs2_of_pressure(actual_pc)) <= 1.0 + 1e-8


def test_tov_checks_pressure_limit_even_with_stale_density_mapping():
    eos = _table()
    domain = find_causal_domain(eos)
    # Reproduce the formerly allowed density from the independent n(epsilon)
    # interpolation, and a caller's old cached density limit.
    old_nb = float(eos.nb_of_energy(domain.eps_max))
    stale = replace(domain, nb_max=old_nb)
    assert float(eos.pressure_of_nb(old_nb)) > domain.pressure_max * 1.0001
    with pytest.raises(ValueError, match="causal domain"):
        solve_tov_star(eos, old_nb, causal_domain=stale)


def test_central_expansion_uses_the_rhs_energy_density(monkeypatch):
    import compose_eos.tov as module
    eos = _table()
    nc = 0.35
    seen = {}
    def record(rhs, span, state, **kwargs):
        seen["r0"], seen["m0"] = span[0], state[0]
        return SimpleNamespace(t_events=[], success=False, nfev=0, message="IC capture")
    monkeypatch.setattr(module, "solve_ivp", record)
    solve_tov_star(eos, nc)
    central_rho = seen["m0"] / (4.0 * np.pi * seen["r0"]**3 / 3.0)
    expected = float(eos.energy_of_pressure(eos.pressure_of_nb(nc))) * MEV_FM3_TO_PA / C**2
    assert central_rho == pytest.approx(expected, rel=2e-14)


def test_first_unresolved_interval_is_not_overridden_by_later_violation():
    eos = _table()
    x, _, safe, undecided = _certify_rhs_causality(eos, 1e-8, min_width=0.5)
    assert undecided
    assert x == safe == pytest.approx(np.log(eos.pressure[0]))


@pytest.mark.parametrize("field", ["causal_tolerance", "max_step_m", "max_radius_m", "initial_radius_m"])
@pytest.mark.parametrize("value", [np.nan, np.inf])
def test_config_refuses_nonfinite_controls(field, value):
    with pytest.raises(ValueError, match="finite"):
        TOVConfig(**{field: value})


@pytest.mark.parametrize("value", [np.nan, np.inf])
def test_domain_refuses_nonfinite_tolerance(value):
    with pytest.raises(ValueError, match="finite"):
        find_causal_domain(_table(), tolerance=value)


def _phase_fixture():
    from compose_eos import EOSData
    from compose_eos.phase_transitions import MaxwellEOS
    low = _table()
    pt = low.pressure[-1]
    et = low.energy_density[-1]
    ph = np.linspace(pt, 3000.0, 10)
    eh = et + 200.0 + 2.0 * (ph - pt)
    n0 = low.nb[-1] * (eh[0] + pt) / (et + pt)
    nh = n0 * ((eh + ph) / (eh[0] + pt))**(2.0/3.0)
    n = np.r_[low.nb, nh]
    p = np.r_[low.pressure, ph]
    eps = np.r_[low.energy_density, eh]
    mu = (eps + p) / n
    z = np.zeros_like(n)
    data = EOSData(nb=n, mn=939.0, mp=939.0, Q1=p/n, Q2=z,
        Q3=mu/939.0-1.0, Q4=z, Q5=z, Q6=eps/n/939.0-1.0,
        Q7=eps/n/939.0-1.0, q_extra=z, Iphase=np.r_[np.ones(low.nb.size,int),np.full(nh.size,2)],
        Yn=z, Yp=z, Ye=z, Ymu=z, Aav=z, Zav=z, Yh=z, composition_available=False,
        beta_equilibrated=True)
    return MaxwellEOS(data, {"pressure_MeV_fm3":pt, "low_index":low.nb.size-1,
        "high_index":low.nb.size, "delta_eps_MeV_fm3":200.0,
        "seidov_small_core_stable":True})


def test_phase_tov_refuses_acausal_shell_with_causal_center():
    from compose_eos.phase_transitions import solve_tov_with_phase_transition
    eos = _phase_fixture()
    # Central high-phase speed is about 0.5; the low-phase shell passes through
    # the acausal pocket near P=530 MeV fm^-3.
    with pytest.raises(ValueError, match="causal"):
        solve_tov_with_phase_transition(eos, 2500.0)


def test_phase_tov_allows_causal_prefix_below_unused_acausal_shell():
    from compose_eos.phase_transitions import solve_tov_with_phase_transition
    result = solve_tov_with_phase_transition(_phase_fixture(), 200.0)
    assert result["status"] == "SUCCESS"
    assert result["phase_interfaces"] == []


def _stars(masses):
    from compose_eos.tov import TOVStarResult
    return [TOVStarResult("SUCCESS" if m is not None else "SOLVER_FAILURE",
        0.1*(i+1), m, 14.0-i if m is not None else None, m is not None,
        1e-10, 10, "test") for i,m in enumerate(masses)]


def test_sequence_summary_stops_at_first_turnover_before_second_branch():
    from compose_eos.tov import stable_branch_summary
    summary = stable_branch_summary(_stars([0.5, 1.3, 1.2, 1.6, 1.8]))
    assert summary["Mmax_sampled_solar"] == pytest.approx(1.3)
    assert summary["R14_km_linear_on_stable_branch"] is None
    assert summary["radial_stability_validated"] is False
    assert summary["max_on_sampled_segment_solar"] == pytest.approx(1.3)


def test_sequence_summary_never_interpolates_across_failed_star():
    from compose_eos.tov import stable_branch_summary
    summary = stable_branch_summary(_stars([0.5, 1.1, None, 1.6, 1.8]))
    assert summary["R14_km_linear_on_stable_branch"] is None
    assert summary["Mmax_sampled_solar"] == pytest.approx(1.1)
    assert summary["Mmax_peak_classification"] == "GAP_LIMITED"


def test_sequence_summary_stops_at_unresolved_plateau():
    from compose_eos.tov import stable_branch_summary
    summary = stable_branch_summary(_stars([0.5, 1.2, 1.2, 1.6]))
    assert summary["Mmax_sampled_solar"] == pytest.approx(1.2)
    assert summary["Mmax_peak_classification"] == "PLATEAU_LIMITED"
    assert summary["R14_km_linear_on_stable_branch"] is None
