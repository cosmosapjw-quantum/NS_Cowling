"""Generated crust+CSS+core EOS -> phase-segmented TOV / Cowling adapter.

The native CSS matcher produces one table with a first-order crust->CSS jump
(two rows with equal P and mu_B) and a continuous CSS->core join.  The adapter
labels the rows below/above the jump as two phases for the phase-segmented
reference TOV/Cowling path.  The CLD inner crusts carry no outer crust, so the
stellar surface is the crust table's minimum pressure (truncation radius).
"""
import json
from pathlib import Path

import numpy as np
import pytest

from compose_eos.native import run_native
from compose_cowling.css_adapter import css_match_maxwell_eos, phase_table_from_css_match
from compose_cowling.direct_barotrope import DirectBarotrope
from compose_cowling.maxwell import (
    MEV_FM3_TO_J_M3, shoot_maxwell_cowling, solve_maxwell_tov)

ROOT = Path(__file__).resolve().parents[2]
CRUST = ROOT / "download/reconstruction_v3/SPG_M1"
CORE = ROOT / "eos_generator/eos_data/GM_GM1"


@pytest.fixture(scope="module")
def matched(tmp_path_factory):
    if not (CRUST / "eos.thermo").is_file():
        pytest.skip("v3 SPG_M1 crust not present")
    tmp = tmp_path_factory.mktemp("css")
    dirs = {k: tmp / k for k in ("crust", "core", "matched", "compose")}
    for d in dirs.values():
        d.mkdir()
    run_native("read", str(CRUST), str(dirs["crust"]))
    run_native("read", str(CORE), str(dirs["core"]))
    report = run_native("match-core", str(dirs["crust"]), str(dirs["core"]),
                        str(dirs["matched"]), "css", "8", "true")
    run_native("write", str(dirs["matched"]), str(dirs["compose"]))
    return dirs["compose"], report


def test_phase_table_splits_exactly_at_the_crust_css_jump(matched):
    compose, report = matched
    table = phase_table_from_css_match(compose, report)
    lo, hi = report["css_transition_low_index"], report["css_transition_high_index"]
    np.testing.assert_array_equal(table.Iphase, report["output_phase_ids"])
    assert table.Iphase[lo] != table.Iphase[hi]
    assert abs(table.p[lo] - table.p[hi]) <= 1e-10 * table.p[lo] and table.nb[hi] > table.nb[lo]


def test_adapter_rejects_report_that_is_not_a_css_match(matched):
    compose, report = matched
    bad = dict(report, construction_kind="smooth_near_pressure_plateau")
    with pytest.raises(ValueError, match="css"):
        phase_table_from_css_match(compose, bad)
    bad = dict(report, css_transition_high_index=report["css_transition_low_index"] + 2)
    with pytest.raises(ValueError, match="index"):
        phase_table_from_css_match(compose, bad)


def test_maxwell_eos_from_generated_css_table_keeps_interface_thermodynamics(matched):
    compose, report = matched
    eos, meta = css_match_maxwell_eos(compose, report)
    t = eos.interface_thermo
    assert t["maxwell_identity_relative_residual"] < 1e-10
    assert abs(t["pressure_mev_fm3"] - report["css_transition_pressure"]) <= 1e-12 * t["pressure_mev_fm3"]
    assert abs(t["chemical_potential_mev"] - report["css_transition_mu"]) <= 2e-8 * t["chemical_potential_mev"]
    assert meta["low_branch_first_law_anchor"] == "interface"
    assert np.isfinite(meta["low_branch_surface_density_drift_relative"])
    assert meta["surface_definition"].startswith("truncation")
    assert meta["pressure_cutoff_mev_fm3"] == pytest.approx(eos.low.pressure[0])


def test_phase_segmented_tov_and_cowling_run_on_generated_eos(matched):
    compose, report = matched
    eos, _ = css_match_maxwell_eos(compose, report)
    star = solve_maxwell_tov(eos, 60.0 * MEV_FM3_TO_J_M3)
    assert 0.0 < star.mass_solar < 3.0
    assert 0.0 < star.interface.pressure_pa < star.p_c_pa
    assert star.inner.r_end_m < star.radius_m
    slow = shoot_maxwell_cowling(star, 1500.0, conversion="slow", max_step_m=200.0)
    rapid = shoot_maxwell_cowling(star, 1500.0, conversion="rapid", max_step_m=200.0)
    assert np.isfinite(slow.residual) and np.isfinite(rapid.residual)
    assert slow.residual != rapid.residual


def test_interface_anchor_reproduces_last_knot_exactly():
    p = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    e = np.array([100.0, 150.0, 210.0, 280.0, 360.0])
    n = np.array([0.10, 0.14, 0.19, 0.24, 0.30])
    mu = (e + p) / n
    last = DirectBarotrope(p, e, n, mu, anchor="last")
    first = DirectBarotrope(p, e, n, mu)
    assert last.evaluate(p[-1]).number_density == pytest.approx(n[-1], rel=1e-14)
    assert first.evaluate(p[0]).number_density == pytest.approx(n[0], rel=1e-14)
    with pytest.raises(ValueError, match="anchor"):
        DirectBarotrope(p, e, n, mu, anchor="middle")


@pytest.mark.parametrize("rel, ok", [(2.2e-16, True), (5e-11, True), (1e-6, False)])
def test_adapter_uses_declared_pressure_gate_not_bit_equality(matched, monkeypatch, rel, ok):
    """CompOSE round trip re-forms P = n_B*Q1, so the two jump rows may differ by an ulp."""
    import compose_cowling.css_adapter as adapter
    compose, report = matched
    real = adapter.read_compose(str(compose))
    hi = report["css_transition_high_index"]

    def fake_read(_path):
        data = __import__("copy").deepcopy(real)
        data.Q1[hi] *= 1.0 + rel
        return data

    monkeypatch.setattr(adapter, "read_compose", fake_read)
    if ok:
        table = adapter.phase_table_from_css_match(compose, report)
        assert table.Iphase[hi] == report["bridge_phase_id"]
    else:
        with pytest.raises(ValueError, match="first-order jump"):
            adapter.phase_table_from_css_match(compose, report)


def test_css_provider_preserves_exact_constant_sound_speed_and_energy(matched):
    compose, report = matched
    eos, _ = css_match_maxwell_eos(compose, report)
    join = report["css_core_bridge_index"] - report["css_transition_high_index"]
    pressures = np.geomspace(eos.high.pressure[0], eos.high.pressure[join], 101)[:-1]
    speed = report["css_single_sound_speed"]
    eps0 = report["css_transition_epsilon_high"]
    pt = report["css_transition_pressure"]
    for pressure in pressures:
        state = eos.high.evaluate(pressure)
        assert state.sound_speed_sq == pytest.approx(speed, rel=2e-13)
        assert state.energy == pytest.approx(eps0 + (pressure - pt) / speed, rel=2e-13)


def test_css_core_join_preserves_one_sided_derivatives_and_continuous_state(matched):
    compose, report = matched
    eos, _ = css_match_maxwell_eos(compose, report)
    join = report["css_core_bridge_index"] - report["css_transition_high_index"]
    pressure = eos.high.pressure[join]
    left = eos.high.evaluate(np.nextafter(pressure, 0.))
    right = eos.high.evaluate(pressure)
    assert left.sound_speed_sq == pytest.approx(report["css_single_sound_speed"], rel=2e-13)
    assert right.sound_speed_sq != pytest.approx(left.sound_speed_sq, rel=1e-5)
    assert right.energy == pytest.approx(left.energy, rel=2e-12)
    assert right.number_density == pytest.approx(left.number_density, rel=2e-12)
    assert right.chemical_potential == pytest.approx(left.chemical_potential, rel=2e-12)


@pytest.mark.parametrize("field,value", [("css_single_sound_speed", .75),
    ("css_core_bridge_index", 1), ("css_core_bridge_index", 4.5)])
def test_css_provider_rejects_report_that_disagrees_with_bridge(matched, field, value):
    compose, report = matched
    with pytest.raises(ValueError, match="CSS|css"):
        css_match_maxwell_eos(compose, dict(report, **{field: value}))


def test_css_provider_rejects_overflowing_analytic_bridge(matched):
    compose, report = matched
    with pytest.raises(ValueError, match="CSS|css"):
        css_match_maxwell_eos(compose, dict(report, css_single_sound_speed=np.nextafter(0., 1.)))


@pytest.mark.parametrize("side", ["crust", "core"])
def test_stellar_intake_keeps_boundaries_and_refuses_undeclared_interpolation(matched, monkeypatch, side):
    import compose_cowling.css_adapter as adapter
    compose, report = matched
    data = adapter.read_compose(compose)
    index = (report["css_transition_low_index"] // 2 if side == "crust"
             else report["css_core_bridge_index"] + 2)
    data.Iphase[index] = 123456
    receipt = dict(report, output_phase_ids=data.Iphase.tolist())
    monkeypatch.setattr(adapter, "read_compose", lambda _: data.copy())
    # The intake result and source keep every phase; no provider is constructed.
    np.testing.assert_array_equal(adapter.phase_table_from_css_match(compose, receipt).Iphase, data.Iphase)
    with pytest.raises(adapter.UnsupportedSourcePhaseBoundary) as exc:
        adapter.css_match_maxwell_eos(compose, receipt)
    assert (index - 1, index) in exc.value.boundaries
    assert data.Iphase[index] == 123456


def test_stellar_intake_rejects_erased_phase_sidecar(matched, monkeypatch):
    import compose_cowling.css_adapter as adapter
    compose, report = matched
    data = adapter.read_compose(compose)
    data.Iphase[:] = 0
    monkeypatch.setattr(adapter, "read_compose", lambda _: data)
    with pytest.raises(ValueError, match="phase layout"):
        adapter.phase_table_from_css_match(compose, report)
