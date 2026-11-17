"""Analytic controls for physical decoding and phase/missingness boundaries."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location(
    "reference_comparison", Path(__file__).resolve().parents[2] / "scripts/compare_generated_reference.py")
module = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = module
SPEC.loader.exec_module(module)


def analytic_table(path, *, mass, nb_origin=1, t_origin=0, y_origin=1,
                   phase=None, missing_particle=(), yq=0., densities=None, pressure_scale=1.):
    path.mkdir()
    n = np.array([.1, .2, .3, .4] if densities is None else densities)
    if phase is None:
        phase = [3] * len(n)
    p = pressure_scale * 2 * n ** 2
    eps = 900 * n + 2 * n ** 2
    mu = 900 + 4 * n
    for name, origin, values in (("nb", nb_origin, n), ("t", t_origin, [0.]), ("yq", y_origin, [yq])):
        (path / f"eos.{name}").write_text(f"{origin}\n{origin + len(values) - 1}\n" +
                                         "\n".join(f"{v:.17g}" for v in values) + "\n")
    rows = [f"{mass:.17g} {mass - 1:.17g} 1"]
    compo = []
    for i, density in enumerate(n):
        q = [p[i] / density, 0., mu[i] / mass - 1., 0., 0., eps[i] / density / mass - 1., eps[i] / density / mass - 1.]
        rows.append(f"{t_origin} {nb_origin+i} {y_origin} " + " ".join(f"{v:.17g}" for v in q) + " 0")
        particle = "0 0" if i in missing_particle else "1 11 0.1 0"
        compo.append(f"{t_origin} {nb_origin+i} {y_origin} {phase[i]} {particle}")
    (path / "eos.thermo").write_text("\n".join(rows) + "\n")
    (path / "eos.compo").write_text("\n".join(compo) + "\n")
    return module.read_cold_table(path)


def compare(g, r, **kwargs):
    return module.compare_tables(g, r, generated_phases={3}, reference_phases={3},
                                 density_min=.1, density_max=.4, **kwargs)


def test_reference_mass_reencoding_preserves_physics_with_arbitrary_origins(tmp_path):
    g = analytic_table(tmp_path / "g", mass=939., nb_origin=7, t_origin=-3, y_origin=5)
    r = analytic_table(tmp_path / "r", mass=1000., nb_origin=-2, t_origin=9, y_origin=-1)
    result, records = compare(g, r)
    assert len(records) == 4
    assert g.metadata["axis_origins"] == {"t": -3, "nb": 7, "yq": 5}
    for field in module.THERMO_FIELDS:
        assert result["metrics"][field]["max_abs_difference"] < 3e-13
    assert (tmp_path / "g/eos.thermo").read_bytes() != (tmp_path / "r/eos.thermo").read_bytes()


def test_missing_composition_not_zero_or_interpolated_across_gap(tmp_path):
    g = analytic_table(tmp_path / "g", mass=939., missing_particle=(1,))
    r = analytic_table(tmp_path / "r", mass=939.)
    result, rows = compare(g, r, composition_fields=("y_p_free",))
    assert np.isnan(g.fields["y_p_free"][1])
    assert np.isnan(rows[1]["generated_y_p_free"])
    assert result["metrics"]["y_p_free"]["compared_rows"] == 3
    pchip, _ = module._interpolate_segment(g.n, g.fields["y_p_free"], np.array([.15, .25]))
    assert np.all(np.isnan(pchip))


def test_same_phase_id_does_not_bridge_a_different_phase(tmp_path):
    g = analytic_table(tmp_path / "g", mass=939., phase=(3, 4, 3, 3))
    r = analytic_table(tmp_path / "r", mass=939.)
    result, rows = compare(g, r)
    assert [r["n_b"] for r in rows] == [.3, .4]
    assert len(result["overlaps"]) == 1


def test_q3_mu_is_not_silently_replaced_by_euler_enthalpy(tmp_path):
    table = analytic_table(tmp_path / "table", mass=939.)
    thermo = tmp_path / "table/eos.thermo"
    rows = thermo.read_text().splitlines()
    tokens = rows[1].split()
    tokens[5] = str(float(tokens[5]) + 1. / 939.)
    rows[1] = " ".join(tokens)
    thermo.write_text("\n".join(rows) + "\n")
    changed = module.read_cold_table(thermo.parent)
    assert changed.fields["mu_b"][0] - table.fields["mu_b"][0] == pytest.approx(1.)
    assert changed.fields["enthalpy_b"][0] == table.fields["enthalpy_b"][0]
    assert changed.metadata["max_abs_enthalpy_minus_mu_b_MeV"] == pytest.approx(1.)


def test_different_physical_charge_slices_rejected(tmp_path):
    g = analytic_table(tmp_path / "g", mass=939.)
    r = analytic_table(tmp_path / "r", mass=939., yq=.1)
    with pytest.raises(ValueError, match="charge-fraction"):
        compare(g, r)


def test_consistent_scalar_temperature_alias_recorded_and_mixed_alias_rejected(tmp_path):
    table = analytic_table(tmp_path / "table", mass=939., t_origin=1)
    (table.path / "eos.t").write_text("0\n0\n0\n")
    changed = module.read_cold_table(table.path)
    assert changed.metadata["axis_origins"]["t"] == 0
    assert changed.metadata["actual_scalar_row_labels"]["t"] == 1
    thermo = table.path / "eos.thermo"
    rows = thermo.read_text().splitlines()
    tokens = rows[1].split()
    tokens[0] = "0"
    rows[1] = " ".join(tokens)
    thermo.write_text("\n".join(rows) + "\n")
    with pytest.raises(ValueError, match="scalar T row labels"):
        module.read_cold_table(table.path)


def test_relative_difference_at_zero_not_infinite():
    summary = module._summary(np.array([.1, .2]), np.array([1., 2.]),
                              np.array([0., 2.]), np.array([1., 2.]))
    assert summary["relative_defined_rows"] == 1
    assert summary["max_abs_difference"] == 1.
    assert summary["max_abs_relative_difference"] == 0.


def test_derivative_diagnostic_requires_interior_phase_support(tmp_path):
    g = analytic_table(tmp_path / "g", mass=939.)
    r = analytic_table(tmp_path / "r", mass=939.)
    result, rows = compare(g, r)
    assert result["first_law_diagnostic"]["reference"]["compared_rows"] == 0
    assert all(np.isnan(row["reference_pressure_from_energy_cubic"]) for row in rows)


def test_pressure_difference_decomposes_without_adjusting_energy(tmp_path):
    density = np.linspace(.1, .4, 9)
    g = analytic_table(tmp_path / "g", mass=939., densities=density)
    r = analytic_table(tmp_path / "r", mass=939., densities=density, pressure_scale=1.2)
    result, _ = compare(g, r)
    assert result["metrics"]["pressure"]["max_abs_relative_difference"] == pytest.approx(1 / 6)
    d = result["first_law_diagnostic"]
    assert d["reference"]["max_abs_relative_difference"] == pytest.approx(1 / 6, abs=1e-10)
    residual = d["pressure_difference_after_reference_energy_derivative_diagnostic"]
    assert residual["compared_rows"] == 5
    assert residual["max_abs_relative_difference"] < 1e-10
