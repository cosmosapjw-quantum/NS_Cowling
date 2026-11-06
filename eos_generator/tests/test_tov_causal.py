import numpy as np
import pytest

from compose_eos.tov import BarotropicEOS, TOVConfig, find_causal_domain, solve_tov_sequence, solve_tov_star


def test_pchip_loglog_power_law_and_derivative_are_exact():
    eps = np.geomspace(1.0, 10.0, 32)
    nb = np.geomspace(0.05, 0.5, 32)
    p = 0.03 * eps ** 1.7
    eos = BarotropicEOS(nb=nb, pressure=p, energy_density=eps)
    q = np.geomspace(1.1, 9.0, 17)
    expected_p = 0.03 * q ** 1.7
    expected_cs2 = 1.7 * expected_p / q
    assert np.allclose(eos.pressure_of_energy(q), expected_p, rtol=2e-12)
    assert np.allclose(eos.cs2(q), expected_cs2, rtol=2e-10)


def test_causal_cutoff_interpolates_first_cs2_equals_one_crossing():
    eps = np.geomspace(1.0, 5.0, 64)
    nb = eps ** 0.6
    p = 0.2 * eps**2
    eos = BarotropicEOS(nb=nb, pressure=p, energy_density=eps)
    domain = find_causal_domain(eos)
    assert domain.status == "ACASUAL_CUTOFF"
    assert domain.eps_max == pytest.approx(2.5, rel=2e-4)
    assert domain.cs2_at_limit == pytest.approx(1.0, rel=2e-4)


def test_barotrope_rejects_nonmonotone_pressure_and_out_of_range_query():
    with pytest.raises(ValueError, match="strictly increasing"):
        BarotropicEOS([1, 2, 3, 4], [1, 0.9, 2, 3], [2, 3, 4, 5])
    eos = BarotropicEOS([1, 2, 3, 4], [1, 2, 3, 4], [2, 3, 4, 5])
    with pytest.raises(ValueError, match="outside tabulated range"):
        eos.pressure_of_energy(6.0)


def test_tov_star_terminates_at_surface_and_returns_physical_values():
    nb = np.geomspace(0.02, 1.1, 180)
    # Smooth, monotone, causal test EOS in CompOSE units (MeV fm^-3).
    eps = 930.0 * nb + 180.0 * nb**2
    p = 90.0 * nb**2
    eos = BarotropicEOS(nb=nb, pressure=p, energy_density=eps)
    cfg = TOVConfig(rtol=2e-8, max_radius_m=3.0e6, max_step_m=2.0e4)
    result = solve_tov_star(eos, central_nb=0.55, config=cfg)
    assert result.status == "SUCCESS"
    assert result.surface_event
    assert 0.1 < result.mass_solar < 4.0
    assert 4.0 < result.radius_km < 30.0
    assert result.surface_pressure_mev_fm3 == pytest.approx(p[0])


def test_tov_star_refuses_central_density_outside_causal_domain():
    eps = np.geomspace(1.0, 5.0, 64)
    nb = eps ** 0.6
    p = 0.2 * eps**2
    eos = BarotropicEOS(nb=nb, pressure=p, energy_density=eps)
    domain = find_causal_domain(eos)
    with pytest.raises(ValueError, match="causal domain"):
        solve_tov_star(eos, central_nb=nb[-1], causal_domain=domain)


def test_first_law_nb_diagnostic_matches_analytic_polytrope():
    eps = np.geomspace(1.0, 20.0, 96)
    k = 0.02
    p = k * eps**2
    nb = 0.1 * eps / (1.0 + k * eps)
    eos = BarotropicEOS(nb=nb, pressure=p, energy_density=eps)
    from compose_eos.tov import first_law_nb_consistency
    result = first_law_nb_consistency(eos, oversampling=24)
    assert result["max_abs_log_ratio"] < 2e-6


def test_causal_domain_reports_full_table_when_no_crossing():
    eps = np.geomspace(1.0, 10.0, 48)
    p = 0.02 * eps**1.2
    nb = eps**0.8
    eos = BarotropicEOS(nb, p, eps)
    domain = find_causal_domain(eos)
    assert domain.status == "FULL_TABLE_CAUSAL"
    assert domain.nb_max == pytest.approx(nb[-1])


def test_barotrope_rejects_nan_and_nonpositive_entries():
    with pytest.raises(ValueError, match="finite"):
        BarotropicEOS([1, 2, 3, 4], [1, 2, np.nan, 4], [2, 3, 4, 5])
    with pytest.raises(ValueError, match="positive"):
        BarotropicEOS([1, 2, 3, 4], [1, 2, 0, 4], [2, 3, 4, 5])


def test_barotrope_rejects_mismatched_shapes():
    with pytest.raises(ValueError, match="same-length"):
        BarotropicEOS([1, 2, 3, 4], [1, 2, 3], [2, 3, 4, 5])


def test_tov_config_rejects_invalid_tolerances_and_steps():
    with pytest.raises(ValueError, match="rtol"):
        TOVConfig(rtol=1.0)
    with pytest.raises(ValueError, match="positive"):
        TOVConfig(max_step_m=0.0)


def test_sequence_rejects_too_few_points():
    eos = BarotropicEOS(np.geomspace(0.1, 1.0, 16),
                        np.geomspace(0.01, 10.0, 16),
                        np.geomspace(100.0, 1000.0, 16))
    with pytest.raises(ValueError, match="at least 4"):
        solve_tov_sequence(eos, n_points=3)


def test_sequence_summary_marks_r14_absent_when_not_bracketed():
    from compose_eos.tov import TOVStarResult, stable_branch_summary
    rows = [
        TOVStarResult("SUCCESS", 0.1, 0.2, 14.0, True, 1e-10, 100, "ok"),
        TOVStarResult("SUCCESS", 0.2, 0.7, 13.0, True, 1e-10, 100, "ok"),
        TOVStarResult("SUCCESS", 0.3, 1.1, 12.0, True, 1e-10, 100, "ok"),
        TOVStarResult("SUCCESS", 0.4, 1.0, 11.0, True, 1e-10, 100, "ok"),
    ]
    summary = stable_branch_summary(rows)
    assert summary["R14_km_linear_on_stable_branch"] is None
    assert summary["R14_in_sampled_stable_branch"] is False


def test_tov_star_reports_no_surface_event_not_a_fabricated_star():
    nb = np.geomspace(0.02, 1.1, 80)
    eps = 930.0 * nb + 180.0 * nb**2
    p = 90.0 * nb**2
    eos = BarotropicEOS(nb, p, eps)
    cfg = TOVConfig(rtol=1e-8, max_radius_m=100.0, max_step_m=10.0)
    result = solve_tov_star(eos, 0.55, config=cfg)
    assert result.status == "NO_SURFACE_EVENT"
    assert result.mass_solar is None
    assert result.radius_km is None


def test_plotting_handles_no_gap_and_multiple_gap_indices(tmp_path):
    from compose_eos import EOSData, make_figures
    nb = np.array([0.01, 0.02, 0.03, 0.2, 0.21, 0.22, 2.0, 2.1])
    f_per_baryon = 930.0 + 100.0 * nb
    pressure = 100.0 * nb**2
    eps = 930.0 * nb + 100.0 * nb**2
    mn = 939.0
    data = EOSData(
        nb=nb, mn=mn, mp=mn,
        Q1=pressure / nb, Q2=np.zeros_like(nb),
        Q3=(930.0 + 200.0 * nb) / mn - 1.0,
        Q4=np.zeros_like(nb), Q5=np.zeros_like(nb),
        Q6=f_per_baryon / mn - 1.0,
        Q7=eps / nb / mn - 1.0,
        Iphase=np.zeros_like(nb, dtype=int),
        Yn=np.zeros_like(nb), Yp=np.zeros_like(nb), Ye=np.zeros_like(nb),
        Ymu=np.zeros_like(nb), Aav=np.zeros_like(nb), Zav=np.zeros_like(nb),
        Yh=np.zeros_like(nb), composition_available=False,
    )
    out = make_figures(data, str(tmp_path), formats=("png",))
    assert len(out) == 5
    assert all(__import__("pathlib").Path(p).is_file() for p in out)


def test_sequence_peak_at_upper_boundary_is_not_claimed_as_global_mmax():
    from compose_eos.tov import TOVStarResult, stable_branch_summary
    rows = [
        TOVStarResult("SUCCESS", 0.1, 0.4, 14.0, True, 1e-10, 100, "ok"),
        TOVStarResult("SUCCESS", 0.2, 0.9, 12.0, True, 1e-10, 100, "ok"),
        TOVStarResult("SUCCESS", 0.3, 1.3, 11.0, True, 1e-10, 100, "ok"),
    ]
    summary = stable_branch_summary(rows)
    assert summary["Mmax_peak_classification"] == "UPPER_DOMAIN_LIMITED"
    assert summary["Mmax_is_global_maximum_claimed"] is False


def test_causality_scan_detects_between_knot_supersonic_derivative():
    increments = np.array([0.22919693, 0.16743684, 0.16919822, 0.1665051,
                           0.17724121, 1.65216718, 0.63981576, 0.1990403,
                           0.65375929])
    logeps = np.linspace(0.0, 5.0, len(increments) + 1)
    logp = np.r_[0.0, np.cumsum(increments)]
    eps = np.exp(logeps)
    pressure = np.exp(logp)
    nb = eps.copy()
    eos = BarotropicEOS(nb, pressure, eps)
    assert float(np.max(eos.cs2(eps))) < 1.0
    domain = find_causal_domain(eos, samples_per_interval=16)
    assert domain.status == "ACASUAL_CUTOFF"
    assert domain.eps_max < eps[-1]
    assert domain.cs2_at_limit == pytest.approx(1.0, abs=1e-8)
