from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import build_crust_only, read_compose, write_compose, NativeExecutionError


EXPECTED = {
    "SPG_M1": 365, "SPG_M2": 377, "SPG_M3": 366, "SPG_M4": 375, "SPG_M5": 363,
    "GPPVA_DD2": 283, "GPPVA_DDME2": 288, "GPPVA_FSU2": 269, "GPPVA_FSU2H": 300,
    "GPPVA_FSU2R": 296, "GPPVA_NL3wrL55": 298, "GPPVA_TM1e": 303,
}


def test_e390_03b_builds_all_twelve_admitted_crust_only_tables(tmp_path: Path) -> None:
    source_root = ROOT / "eos_generator" / "eos_data"
    output_root = tmp_path / "crust-only-attempt"
    receipts = build_crust_only(output_root, source_root=source_root)

    assert len(receipts) == 12
    manifest = json.loads((output_root / "crust_only_manifest.v1.json").read_text())
    coverage = json.loads((output_root / "crust_coverage_report.v1.json").read_text())
    assert manifest["status"] == "ALL_ADMITTED_CRUST_COMPONENTS_CONSTRUCTED_NATIVE_READBACK_VERIFIED"
    assert manifest["total_output_rows"] == 3883
    assert coverage["admitted_model_count"] == 12
    assert coverage["blocked"]["GPPVA_TW"].startswith("BLOCKED")

    source_fields = (
        "nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "q_extra",
        "Iphase", "Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh",
        "row_status", "Yn_present", "Yp_present", "Ye_present",
        "Ymu_present", "Abar_present", "Zbar_present", "Yh_present",
        "C3_applicable", "C4_applicable",
    )
    for receipt in receipts:
        model = receipt["model_id"]
        n = EXPECTED[model]
        selection = receipt["selection"]
        assert selection["first_row"] == 1
        assert selection["outer_crust_rows"] == [1, 217]
        assert selection["inner_crust_rows"] == [218, n]
        assert selection["last_row"] == n
        assert selection["output_point_count"] == n
        assert selection["source_grid_has_no_rows_between_drip_and_provider_splice"] is True
        assert selection["no_interpolation_or_extrapolation"] is True
        assert selection["no_homogeneous_core_tail_imported"] is True
        assert receipt["checks"]["native_readback"] == "PASS"
        assert receipt["checks"]["BSk22_neutron_drip_boundary"] == "PASS"
        assert receipt["checks"]["C3_C4_masks_exact_preservation"] == "PASS"
        assert receipt["checks"]["source_row_status_lossless_preservation"] == "PASS"

        source = read_compose(source_root / model)
        output_dir = Path(receipt["output"]["directory"])
        crust = read_compose(output_dir)
        assert crust.n == source.n == n
        assert crust.inb_start == source.inb_start
        assert np.all(np.diff(crust.nb) > 0.0)
        assert np.all(crust.Iphase[:217] == 0)
        assert np.all(crust.Iphase[217:] == 7)
        assert np.isclose(crust.nb[216], 2.69021e-4, rtol=0.0, atol=2e-9)
        assert np.searchsorted(crust.nb, 2.0e-3, side="left") + 1 == 218
        for name in source_fields:
            np.testing.assert_array_equal(getattr(crust, name), getattr(source, name))
        for name in ("mn", "mp", "Il", "temperature_mev", "charge_fraction_grid",
                     "q_extra_count", "composition_file_present", "nT", "nYq",
                     "t_first", "t_last", "yq_first", "yq_last",
                     "selected_iT", "selected_iYq"):
            assert getattr(crust, name) == getattr(source, name)
        assert crust.composition_file_present is True
        if model.startswith("SPG_"):
            assert crust.q_extra_count == 0
            assert receipt["construction_provenance"]["inner_crust_method"] == "CLD"
        else:
            assert crust.q_extra_count == 1
            assert set(np.unique(crust.q_extra).tolist()).issubset({0.0, 1.0, 2.0, 3.0})
            assert receipt["construction_provenance"]["inner_crust_method"] == "Thomas-Fermi"
        saved = json.loads((output_dir / "E390-03B_RECEIPT.json").read_text())
        assert saved["model_id"] == model
        assert set(saved["output"]["file_sha256"]) == {
            "eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo"
        }


def test_e390_03b_rejects_gppva_tw_provider_conflict(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="GPPVA_TW remains blocked"):
        build_crust_only(
            tmp_path / "tw-attempt",
            source_root=ROOT / "eos_generator" / "eos_data",
            model_ids=["GPPVA_TW"],
        )


def test_e390_03b_rejects_modified_source_grid_before_native_extraction(tmp_path: Path) -> None:
    source_root = ROOT / "eos_generator" / "eos_data"
    damaged_root = tmp_path / "damaged-source"
    damaged_root.mkdir()
    shutil.copytree(source_root / "SPG_M1", damaged_root / "SPG_M1")
    nb_file = damaged_root / "SPG_M1" / "eos.nb"
    nb_file.write_bytes(nb_file.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="Git blob identity mismatch"):
        build_crust_only(
            tmp_path / "damaged-output", source_root=damaged_root, model_ids=["SPG_M1"]
        )


def test_e390_03b_never_overwrites_an_existing_attempt(tmp_path: Path) -> None:
    output_root = tmp_path / "existing"
    output_root.mkdir()
    with pytest.raises(FileExistsError, match="fresh attempt path"):
        build_crust_only(
            output_root,
            source_root=ROOT / "eos_generator" / "eos_data",
            model_ids=["SPG_M1"],
        )


def test_e390_03b_preserves_source_fail_status_without_relaxing_default_writer(tmp_path: Path) -> None:
    source_root = ROOT / "eos_generator" / "eos_data"
    source = read_compose(source_root / "SPG_M1")
    assert np.any(source.row_status == 3), "fixture must retain source rows marked FAIL"

    # Normal writes remain fail-closed; only the dedicated crust projection path
    # opts into preserving finite source values without physical repair.
    with pytest.raises(NativeExecutionError, match="row-level composition status is FAIL"):
        write_compose(tmp_path / "strict-output", source)
