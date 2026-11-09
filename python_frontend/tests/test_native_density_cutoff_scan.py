from __future__ import annotations

import math

import pytest

from compose_cowling.native_table_cutoff_scan import select_native_density_knots


def _rows(count: int = 217):
    rows = []
    for i in range(count):
        fraction = i / (count - 1)
        n_b = 10.0 ** (-9.33 + 5.76 * fraction)
        rows.append({
            "source_row": str(i + 1),
            "pressure_pa": f"{10.0 ** (21.0 + 5.0 * fraction):.17e}",
            "energy_mev_fm3": f"{max(1.0e-16, n_b * 930.0):.17e}",
            "n_b_fm3": f"{n_b:.17e}",
            "mu_b_mev": f"{930.0 + 0.1 * fraction:.17e}",
        })
    # A first inner-crust-like row lies above the R15 selector boundary.
    rows.append({
        "source_row": str(count + 1),
        "pressure_pa": f"{10.0 ** 27.0:.17e}",
        "energy_mev_fm3": f"{3.0:.17e}",
        "n_b_fm3": f"{3.0e-3:.17e}",
        "mu_b_mev": f"{941.0:.17e}",
    })
    return rows


def test_log_density_selector_uses_distinct_existing_rows_and_preserves_fields():
    rows = _rows()
    original = [dict(row) for row in rows]

    selected = select_native_density_knots(
        rows, target_count=8, maximum_density_fm3=1.0e-3
    )

    assert len(selected) == 8
    assert selected[0].row_index == 0
    assert selected[-1].row_index == 216
    assert len({k.row_index for k in selected}) == 8
    assert [k.row_index for k in selected] == sorted(k.row_index for k in selected)
    assert all(k.n_b_fm3 <= 1.0e-3 for k in selected)
    assert all(k.row_index < 217 for k in selected)
    for knot in selected:
        source = original[knot.row_index]
        assert knot.source_row == source["source_row"]
        assert knot.pressure_pa == float(source["pressure_pa"])
        assert knot.energy_mev_fm3 == float(source["energy_mev_fm3"])
        assert knot.n_b_fm3 == float(source["n_b_fm3"])
        assert knot.mu_b_mev == float(source["mu_b_mev"])
    assert rows == original


def test_log_density_targets_cover_the_entire_eligible_density_segment():
    rows = _rows()
    selected = select_native_density_knots(rows, target_count=8)

    assert selected[0].n_b_fm3 == float(rows[0]["n_b_fm3"])
    assert selected[-1].n_b_fm3 == float(rows[216]["n_b_fm3"])
    assert selected[0].actual_offset_dex == pytest.approx(0.0, abs=1e-14)
    assert selected[-1].actual_offset_dex == pytest.approx(
        math.log10(float(rows[216]["n_b_fm3"]) / float(rows[0]["n_b_fm3"])),
        rel=0.0,
        abs=1e-14,
    )
    assert all(
        b.target_log10_n_b_fm3 > a.target_log10_n_b_fm3
        for a, b in zip(selected, selected[1:])
    )


def test_log_density_selector_rejects_non_monotone_pressure_or_density():
    rows = _rows()
    rows[20]["pressure_pa"], rows[21]["pressure_pa"] = (
        rows[21]["pressure_pa"], rows[20]["pressure_pa"]
    )
    with pytest.raises(ValueError, match="pressure must be strictly increasing"):
        select_native_density_knots(rows)

    rows = _rows()
    rows[20]["n_b_fm3"], rows[21]["n_b_fm3"] = (
        rows[21]["n_b_fm3"], rows[20]["n_b_fm3"]
    )
    with pytest.raises(ValueError, match="baryon density must be strictly increasing"):
        select_native_density_knots(rows)


def test_log_density_selector_rejects_too_many_cutoffs_and_bad_domain():
    rows = _rows()
    with pytest.raises(ValueError, match="native rows"):
        select_native_density_knots(rows, target_count=218)

    with pytest.raises(ValueError, match="finite and positive"):
        select_native_density_knots(rows, maximum_density_fm3=float("nan"))
