from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import build_core_only, read_compose


def test_e390_03a_builds_all_five_admitted_core_only_tables(tmp_path: Path) -> None:
    source_root = ROOT / "eos_generator" / "eos_data"
    output_root = tmp_path / "core-only-attempt"
    receipts = build_core_only(output_root, source_root=source_root)

    assert len(receipts) == 5
    manifest = json.loads((output_root / "core_only_manifest.v1.json").read_text())
    assert manifest["status"] == "ALL_REQUESTED_COMPONENTS_CONSTRUCTED_AND_NATIVE_READBACK_VERIFIED"

    expected = {
        "APR_APR": (70, 172, 103),
        "BBB_BHF2": (70, 84, 15),
        "GM_GM1": (49, 349, 301),
        "RG_KDE0v1": (967, 1240, 274),
        "RG_SLY4": (144, 197, 54),
    }
    for receipt in receipts:
        model = receipt["model_id"]
        first, last, n = expected[model]
        selection = receipt["selection"]
        assert selection["first_row"] == first
        assert selection["last_row"] == last
        assert selection["output_point_count"] == n
        assert selection["no_interpolation_or_extrapolation"] is True
        assert receipt["checks"]["native_readback"] == "PASS"
        assert receipt["checks"]["Q1_Q7_and_Qextra_exact_row_preservation"] == "PASS"
        assert receipt["checks"]["Iphase_and_row_composition_presence_exact_preservation"] == "PASS"
        for filename, digest in receipt["output"]["file_sha256"].items():
            assert len(digest) == 64
            assert (Path(receipt["output"]["directory"]) / filename).is_file()

        output = read_compose(receipt["output"]["directory"])
        assert output.n == n
        assert np.all(np.diff(output.nb) > 0.0)
        assert output.inb_start >= 1
        saved = json.loads(
            (Path(receipt["output"]["directory"]) / "E390-03A_RECEIPT.json").read_text()
        )
        assert saved["model_id"] == model


def test_e390_03a_rejects_modified_source_grid_before_native_extraction(tmp_path: Path) -> None:
    source_root = ROOT / "eos_generator" / "eos_data"
    damaged_root = tmp_path / "damaged-source"
    damaged_root.mkdir()
    damaged_model = damaged_root / "APR_APR"
    shutil.copytree(source_root / "APR_APR", damaged_model)
    nb_file = damaged_model / "eos.nb"
    nb_file.write_bytes(nb_file.read_bytes() + b"\n")

    with pytest.raises(ValueError, match="Git blob identity mismatch"):
        build_core_only(
            tmp_path / "damaged-output",
            source_root=damaged_root,
            model_ids=["APR_APR"],
        )


def test_e390_03a_never_overwrites_an_existing_attempt(tmp_path: Path) -> None:
    output_root = tmp_path / "existing"
    output_root.mkdir()
    with pytest.raises(FileExistsError, match="fresh attempt path"):
        build_core_only(
            output_root,
            source_root=ROOT / "eos_generator" / "eos_data",
            model_ids=["APR_APR"],
        )
