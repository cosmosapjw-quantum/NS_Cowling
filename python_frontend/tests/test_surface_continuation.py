import csv
import math
from pathlib import Path

import numpy as np
import pytest

from compose_cowling.surface_continuation import (
    NORMALIZED_FIELDS,
    ThermodynamicSurfaceContinuation,
    build_extended_projected_rows,
    write_extended_projected_csv,
)


def raw_rows():
    # The first three real projected knots of GPPVA_FSU2/BBB_BHF2/css.
    return [
        {
            "source_row": "1", "i_t": "0", "i_nb": "1", "i_yq": "1",
            "n_b_fm3": "4.6796200000000002e-10", "n_b_m3": "4.6796199999999998e35",
            "pressure_mev_fm3": "9.0000000794219998e-12", "pressure_pa": "1.4419589833248072e21",
            "energy_mev_fm3": "4.3541347345308269e-07", "energy_j_m3": "6.9760929329530834e25",
            "mu_b_mev": "930.46545115449999", "mu_b_j": "1.4907700045840082e-10",
            "q1_p_over_nb_mev": "0.019232330999999998", "q2_entropy_per_baryon": "0",
            "q3_mu_b_over_mn_minus_1": "-0.0096852007571416321", "q4_mu_q_over_mn": "-0.00066813571529932881",
            "q5_mu_l_over_mn": "0", "q6_f_over_nb_mn_minus_1": "-0.009705670149089074",
            "q7_e_over_nb_mn_minus_1": "-0.009705670149089074",
        },
        {
            "source_row": "2", "i_t": "0", "i_nb": "2", "i_yq": "1",
            "n_b_fm3": "5.8496400000000005e-10", "n_b_m3": "5.8496400000000002e35",
            "pressure_mev_fm3": "1.2858829966793143e-11", "pressure_pa": "2.0602116913374969e21",
            "energy_mev_fm3": "5.4428022096665206e-07", "energy_j_m3": "8.7203305238112679e25",
            "mu_b_mev": "930.47278088330006", "mu_b_j": "1.4907817481042254e-10",
            "q1_p_over_nb_mev": "0.021982258680522462", "q2_entropy_per_baryon": "0",
            "q3_mu_b_over_mn_minus_1": "-0.0096773995659247269", "q4_mu_q_over_mn": "-0.00068613959596614742",
            "q5_mu_l_over_mn": "0", "q6_f_over_nb_mn_minus_1": "-0.0097007878884802778",
            "q7_e_over_nb_mn_minus_1": "-0.0097007957662511579",
        },
        {
            "source_row": "3", "i_t": "0", "i_nb": "3", "i_yq": "1",
            "n_b_fm3": "7.3384599999999987e-10", "n_b_m3": "7.3384599999999981e35",
            "pressure_mev_fm3": "1.8427710474894465e-11", "pressure_pa": "2.9524447140992956e21",
            "energy_mev_fm3": "6.8281149820965563e-07", "energy_j_m3": "1.0939846278580431e26",
            "mu_b_mev": "930.48122619750006", "mu_b_j": "1.4907952789893034e-10",
            "q1_p_over_nb_mev": "0.025111141131646787", "q2_entropy_per_baryon": "0",
            "q3_mu_b_over_mn_minus_1": "-0.0096684110326845651", "q4_mu_q_over_mn": "-0.00070680905180909562",
            "q5_mu_l_over_mn": "0", "q6_f_over_nb_mn_minus_1": "-0.0096951185017650765",
            "q7_e_over_nb_mn_minus_1": "-0.0096951373712688138",
        },
    ]


def test_local_tangent_extension_matches_cutoff_value_slope_and_chemical_potential():
    model = ThermodynamicSurfaceContinuation.from_projected_rows(raw_rows(), alpha=0.0, family="power_law")
    state = model.state_at_energy(model.epsilon_cut_mev_fm3)
    assert state["pressure_mev_fm3"] == pytest.approx(model.pressure_cut_mev_fm3, rel=2e-14)
    assert state["number_density_fm3"] == pytest.approx(model.number_density_cut_fm3, rel=2e-14)
    assert state["chemical_potential_mev"] == pytest.approx(model.chemical_potential_cut_mev, rel=2e-14)
    assert state["cs2"] == pytest.approx(model.cs2_cut, rel=2e-14)
    assert model.gamma_s == pytest.approx(
        model.epsilon_cut_mev_fm3 / model.pressure_cut_mev_fm3 * model.cs2_cut,
        rel=2e-14,
    )
    gamma_eq = (model.epsilon_cut_mev_fm3 + model.pressure_cut_mev_fm3) / model.pressure_cut_mev_fm3 * model.cs2_cut
    assert gamma_eq == pytest.approx(model.gamma_s + model.cs2_cut, rel=2e-14)


@pytest.mark.parametrize("family,alpha", [("exp_quadratic", 0.25), ("poly_quadratic", -0.5)])
def test_first_law_and_chemical_potential_identity_hold_in_the_tail(family, alpha):
    model = ThermodynamicSurfaceContinuation.from_projected_rows(raw_rows(), alpha=alpha, family=family)
    for x in (1e-9, 1e-5, 0.03, 0.4, 0.9):
        eps = model.epsilon_cut_mev_fm3 * x
        h = 2e-5
        sm = model.state_at_energy(eps * math.exp(-h))
        s0 = model.state_at_energy(eps)
        sp = model.state_at_energy(eps * math.exp(h))
        dlnn_dlnE = (math.log(sp["number_density_fm3"]) - math.log(sm["number_density_fm3"])) / (2*h)
        assert dlnn_dlnE == pytest.approx(
            eps / (eps + s0["pressure_mev_fm3"]), rel=2e-8, abs=2e-10
        )
        if x >= 0.03:  # At x<<1, subtracting mu_b~m_n in float64 loses derivative precision.
            h_mu = 3e-4
            sm_mu = model.state_at_energy(eps * math.exp(-h_mu))
            sp_mu = model.state_at_energy(eps * math.exp(h_mu))
            dmu_de = (sp_mu["chemical_potential_mev"] - sm_mu["chemical_potential_mev"]) / (2*h_mu*eps)
            dp_de = (sp_mu["pressure_mev_fm3"] - sm_mu["pressure_mev_fm3"]) / (2*h_mu*eps)
            assert s0["number_density_fm3"] * dmu_de == pytest.approx(dp_de, rel=2e-6)
        assert s0["chemical_potential_mev"] * s0["number_density_fm3"] == pytest.approx(
            eps + s0["pressure_mev_fm3"], rel=2e-14
        )


def test_curvature_family_preserves_cutoff_tangent_and_changes_interior():
    rows = raw_rows()
    baseline = ThermodynamicSurfaceContinuation.from_projected_rows(rows, alpha=0.0)
    alternatives = (
        ThermodynamicSurfaceContinuation.from_projected_rows(rows, alpha=0.5, family="exp_quadratic"),
        ThermodynamicSurfaceContinuation.from_projected_rows(rows, alpha=-0.5, family="exp_quadratic"),
        ThermodynamicSurfaceContinuation.from_projected_rows(rows, alpha=0.5, family="poly_quadratic"),
        ThermodynamicSurfaceContinuation.from_projected_rows(rows, alpha=-0.5, family="poly_quadratic"),
    )
    eps = baseline.epsilon_cut_mev_fm3 * 0.2
    p0 = baseline.state_at_energy(eps)["pressure_mev_fm3"]
    for model in alternatives:
        assert model.state_at_energy(baseline.epsilon_cut_mev_fm3)["cs2"] == pytest.approx(baseline.cs2_cut)
        assert abs(model.state_at_energy(eps)["pressure_mev_fm3"] / p0 - 1.0) > 1e-3
    for model in alternatives:
        for x in np.geomspace(1e-12, 1.0, 100):
            state = model.state_at_energy(model.epsilon_cut_mev_fm3 * x)
            assert 0.0 < state["cs2"] < 1.0
            assert state["number_density_fm3"] > 0.0
            assert state["chemical_potential_mev"] > 0.0


def test_extended_rows_are_positive_ordered_and_leave_source_rows_unchanged():
    original = raw_rows()
    extended = build_extended_projected_rows(
        original, alpha=0.0, tail_points=64, min_energy_ratio=1e-10
    )
    assert len(extended) == len(original) + 64
    assert extended[-len(original):] == original
    for key in ("pressure_pa", "energy_j_m3", "n_b_m3", "mu_b_j"):
        vals = [float(row[key]) for row in extended]
        assert all(math.isfinite(v) and v > 0.0 for v in vals)
        assert all(b > a for a, b in zip(vals, vals[1:]))


def test_invalid_tail_domain_and_unphysical_curvature_are_rejected():
    with pytest.raises(ValueError):
        build_extended_projected_rows(raw_rows(), alpha=0.0, tail_points=8, min_energy_ratio=1e-10)
    with pytest.raises(ValueError):
        build_extended_projected_rows(raw_rows(), alpha=4.0, family="exp_quadratic", tail_points=32, min_energy_ratio=1e-10)


def test_tangent_polytrope_number_density_matches_closed_form_first_law():
    model = ThermodynamicSurfaceContinuation.from_projected_rows(raw_rows(), alpha=0.0, family="power_law")
    beta = model.pressure_cut_mev_fm3 / model.epsilon_cut_mev_fm3
    for x in (1e-12, 1e-8, 0.01, 0.3, 0.95):
        state = model.state_at_energy(model.epsilon_cut_mev_fm3 * x)
        expected = model.number_density_cut_fm3 * x * (
            (1.0 + beta) / (1.0 + beta * x ** (model.gamma_s - 1.0))
        ) ** (1.0 / (model.gamma_s - 1.0))
        assert state["number_density_fm3"] == pytest.approx(expected, rel=3e-13)


def test_writer_emits_native_normalized_header_and_hash_receipts(tmp_path):
    source = tmp_path / "source.csv"
    target = tmp_path / "extended.csv"
    with source.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=NORMALIZED_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(raw_rows())
    receipt = write_extended_projected_csv(
        source, target, alpha=0.0, family="power_law", tail_points=32, min_energy_ratio=1e-8
    )
    with target.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        assert tuple(reader.fieldnames) == NORMALIZED_FIELDS
    assert receipt["source_row_count"] == 3
    assert receipt["output_row_count"] == 35
    assert receipt["tail_row_count"] == 32
    assert receipt["source_sha256"]
    assert receipt["output_sha256"]
    assert rows[-3:] == raw_rows()
    assert float(rows[0]["pressure_pa"]) > 0.0
