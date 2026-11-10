from __future__ import annotations

import pytest

from scripts.r16_minimal_hermite_pressure_repair import minimal_admissible_delta_p
from scripts.r16_minimal_repair_tov_pilot import repair_normalized_rows


def _normalized_rows():
    factor = 1.602176634e32
    mu_j = [1.490e-10, 1.491e-10, 1.492e-10]
    n_m3 = [4.0e36, 4.2e36, 4.5e36]
    mu_mev = [x / 1.602176634e-13 for x in mu_j]
    n_fm3 = [x / 1.0e45 for x in n_m3]
    dp = []
    for i in range(2):
        repair = minimal_admissible_delta_p(
            mu0=mu_j[i], n0=n_m3[i], mu1=mu_j[i+1], n1=n_m3[i+1],
            delta_p_source=1.0e30,
            interior_margin_fraction=1.0e-8,
        )
        dp.append(
            repair["delta_p_lower"] * 0.8
            if i == 0 else 0.5 * (repair["delta_p_lower"] + repair["delta_p_upper"])
        )
    p_pa = [1.0e20, 1.0e20 + dp[0], 1.0e20 + dp[0] + dp[1]]
    rows = []
    for i in range(3):
        eps = mu_j[i] * n_m3[i] - p_pa[i]
        p_mev = p_pa[i] / factor
        eps_mev = eps / factor
        q3 = 0.0
        rows.append({
            "source_row": str(i + 1),
            "mu_b_j": f"{mu_j[i]:.17g}",
            "n_b_m3": f"{n_m3[i]:.17g}",
            "pressure_pa": f"{p_pa[i]:.17g}",
            "energy_j_m3": f"{eps:.17g}",
            "mu_b_mev": f"{mu_mev[i]:.17g}",
            "n_b_fm3": f"{n_fm3[i]:.17g}",
            "pressure_mev_fm3": f"{p_mev:.17g}",
            "energy_mev_fm3": f"{eps_mev:.17g}",
            "q1_p_over_nb_mev": f"{p_mev/n_fm3[i]:.17g}",
            "q3_mu_b_over_mn_minus_1": f"{q3:.17g}",
            "q7_e_over_nb_mn_minus_1": f"{eps_mev/(n_fm3[i]*mu_mev[i])-1.0:.17g}",
        })
    return rows


def test_repaired_table_preserves_knots_and_restores_fortran_cell_admission():
    source = _normalized_rows()
    repaired, audit = repair_normalized_rows(source)

    assert audit["row_count"] == 3
    assert audit["cell_count"] == 2
    assert audit["adjusted_cells"] == 1
    assert audit["preserved_cells"] == 1
    assert [r["source_row"] for r in repaired] == ["1", "2", "3"]
    assert [r["mu_b_j"] for r in repaired] == [r["mu_b_j"] for r in source]
    assert [r["n_b_m3"] for r in repaired] == [r["n_b_m3"] for r in source]
    assert float(repaired[1]["pressure_pa"]) > float(source[1]["pressure_pa"])
    assert float(repaired[2]["pressure_pa"]) > float(source[2]["pressure_pa"])
    for row in repaired:
        lhs = float(row["energy_j_m3"]) + float(row["pressure_pa"])
        rhs = float(row["mu_b_j"]) * float(row["n_b_m3"])
        assert lhs == pytest.approx(rhs, rel=2e-15)
