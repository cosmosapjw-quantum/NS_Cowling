
import pytest

from compose_cowling.native_table_cutoff_scan import (
    load_native_eos_rows,
    select_native_pressure_knots,
)


def _rows():
    return [
        {
            "source_row": str(i + 1),
            "pressure_pa": f"{10.0 ** (21.0 + 0.25 * i):.17e}",
            "energy_mev_fm3": f"{1.0e-7 * (i + 1):.17e}",
            "n_b_fm3": f"{1.0e-10 * (i + 1):.17e}",
            "mu_b_mev": f"{930.0 + 0.1 * i:.17e}",
        }
        for i in range(16)
    ]


def test_cutoff_selector_returns_existing_knots_without_interpolation(tmp_path):
    path = tmp_path / "native.csv"
    rows = _rows()
    header = list(rows[0])
    path.write_text(
        ",".join(header) + "\n"
        + "".join(",".join(row[key] for key in header) + "\n" for row in rows),
        encoding="utf-8",
    )

    original_bytes = path.read_bytes()
    loaded = load_native_eos_rows(path)
    selected = select_native_pressure_knots(loaded, (0.0, 0.5, 1.0, 2.0))

    assert path.read_bytes() == original_bytes
    assert [k.row_index for k in selected] == [0, 2, 4, 8]
    assert [k.source_row for k in selected] == ["1", "3", "5", "9"]
    for knot in selected:
        assert knot.pressure_pa == float(rows[knot.row_index]["pressure_pa"])
        assert knot.energy_mev_fm3 == float(rows[knot.row_index]["energy_mev_fm3"])
        assert knot.n_b_fm3 == float(rows[knot.row_index]["n_b_fm3"])
        assert knot.mu_b_mev == float(rows[knot.row_index]["mu_b_mev"])


def test_cutoff_selector_rejects_non_monotone_native_order():
    rows = _rows()
    rows[4]["pressure_pa"], rows[5]["pressure_pa"] = rows[5]["pressure_pa"], rows[4]["pressure_pa"]
    with pytest.raises(ValueError, match="strictly increasing in source order"):
        select_native_pressure_knots(rows, (0.0, 0.5, 1.0))


def test_cutoff_selector_rejects_targets_that_collapse_to_the_same_native_knot():
    rows = _rows()
    with pytest.raises(ValueError, match="distinct native knots"):
        select_native_pressure_knots(rows, (0.0, 0.01, 0.02))


def test_native_table_loader_requires_physical_columns_and_finite_positive_values(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text(
        "source_row,pressure_pa\n1,1.0\n2,2.0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="missing required EOS columns"):
        load_native_eos_rows(path)

    good = _rows()
    header = list(good[0])
    good[2]["pressure_pa"] = "nan"
    path.write_text(
        ",".join(header) + "\n"
        + "".join(",".join(row[key] for key in header) + "\n" for row in good),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="finite and positive"):
        load_native_eos_rows(path)


def test_cutoff_selector_rejects_target_outside_native_table_domain():
    with pytest.raises(ValueError, match="outside the native EOS table"):
        select_native_pressure_knots(_rows(), (0.0, 10.0))
