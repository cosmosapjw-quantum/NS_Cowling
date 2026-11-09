from __future__ import annotations

import pytest

from scripts.run_native_outer_crust_r15 import PAIRS, METHODS, select_cohort_pairs


def _full_reference_rows():
    cores = ("BBB_BHF2", "GM_GM1", "RG_SLY4")
    rows = []
    for crust_index in range(13):
        crust = f"CRUST_{crust_index:02d}"
        for core in cores:
            for method in METHODS:
                rows.append({
                    "crust": crust,
                    "core": core,
                    "matching_method": method,
                })
    return rows


def test_default_cohort_remains_the_frozen_three_pair_pilot():
    assert select_cohort_pairs(_full_reference_rows()) == PAIRS
    assert len(PAIRS) == 3


def test_full_cohort_resolves_39_complete_pairs():
    pairs = select_cohort_pairs(_full_reference_rows(), full_cohort=True)
    assert len(pairs) == 39
    assert len({crust for crust, _ in pairs}) == 13
    assert {core for _, core in pairs} == {"BBB_BHF2", "GM_GM1", "RG_SLY4"}
    assert pairs == tuple(sorted(pairs))


def test_full_cohort_rejects_incomplete_six_method_pair():
    rows = _full_reference_rows()
    rows = [
        row for row in rows
        if not (
            row["crust"] == "CRUST_00"
            and row["core"] == "BBB_BHF2"
            and row["matching_method"] == "quintic_P"
        )
    ]
    with pytest.raises(RuntimeError, match="not six-method complete"):
        select_cohort_pairs(rows, full_cohort=True)
