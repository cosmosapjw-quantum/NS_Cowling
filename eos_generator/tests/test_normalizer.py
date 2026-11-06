"""Regression tests for the portable CompOSE composition normalizer."""
from pathlib import Path

import pytest

from scripts.normalize_compose_composition import (
    _parse_compo_row,
    normalize_compo_file,
    normalize_tree,
)


def test_empty_compose_row_is_preserved_and_not_invented_as_neutrons():
    raw = "1 2 3 1 0 0"
    normalized, record = _parse_compo_row(raw, Path("eos.compo"), 1)
    assert normalized == raw
    assert record is None


def test_nuclear_group_is_written_as_compose_number_abundance():
    # Yn=0.1, Yp=0.2, Yh=0.7; Z/A=26/56 and Ye=0.525.
    raw = "1 1 1 1 4 10 0.1 11 0.2 0 0.525 1 0 1 0 56 26 0.7"
    normalized, record = _parse_compo_row(raw, Path("eos.compo"), 1)
    fields = normalized.split()
    assert fields[4] == "4"
    nquad_at = 5 + 2 * int(fields[4])
    assert fields[nquad_at] == "1"
    Abar, Zbar, YI = map(float, fields[nquad_at + 2:nquad_at + 5])
    assert Abar == 56.0
    assert Zbar == 26.0
    assert YI == pytest.approx(0.7 / 56.0, rel=1e-14)
    assert record["normalization"]["baryon_residual"] < 1e-14
    assert record["normalization"]["charge_residual"] < 1e-14


def test_reconstructs_missing_free_neutron_from_charge_closure():
    # The provider has no free-neutron species. With Yp=0, charge neutrality
    # fixes Yh=0.7 and baryon closure then fixes Yn=0.3.
    raw = "1 1 1 1 2 11 0 0 0.325 1 0 56 26 0.7"
    normalized, record = _parse_compo_row(raw, Path("eos.compo"), 1)
    fields = normalized.split()
    npairs = int(fields[4])
    pairs = {int(fields[5 + 2*i]): float(fields[6 + 2*i]) for i in range(npairs)}
    assert pairs[10] == pytest.approx(0.3)
    assert pairs[11] == pytest.approx(0.0)
    assert record["adjustments"]["Yn"] == pytest.approx(0.3)


def test_ambiguous_or_multiple_groups_are_rejected_not_guessed():
    raw = "1 1 1 1 0 2 0 56 26 0.5 1 4 2 0.1"
    with pytest.raises(ValueError, match="Nquad=2"):
        _parse_compo_row(raw, Path("eos.compo"), 1)


def test_normalize_file_records_provenance(tmp_path):
    source = tmp_path / "source.compo"
    target = tmp_path / "target.compo"
    source.write_text(
        "1 1 1 1 4 10 0.1 11 0.2 0 0.525 1 0 1 0 56 26 0.7\n",
        encoding="utf-8",
    )
    result = normalize_compo_file(source, target)
    assert result["rows"] == 1
    assert result["nuclear_group_rows"] == 1
    assert result["rows_reconstructed"] == 0
    assert len(result["source_sha256"]) == 64
    assert len(result["normalized_sha256"]) == 64


def test_normalize_tree_does_not_modify_input_and_writes_manifest(tmp_path):
    source = tmp_path / "raw"
    model = source / "MODEL"
    model.mkdir(parents=True)
    (model / "eos.compo").write_text(
        "1 1 1 1 4 10 0.1 11 0.2 0 0.525 1 0 1 0 56 26 0.7\n",
        encoding="utf-8",
    )
    input_hash = (model / "eos.compo").read_bytes()
    target = tmp_path / "normalized"
    manifest = normalize_tree(source, target)
    assert (model / "eos.compo").read_bytes() == input_hash
    assert (target / "MODEL" / "eos.compo").is_file()
    assert (target / "COMPOSITION_NORMALIZATION_MANIFEST.json").is_file()
    assert manifest["models"][0]["status"] == "normalized"


def test_normalize_tree_rejects_output_inside_input(tmp_path):
    source = tmp_path / "raw"
    (source / "MODEL").mkdir(parents=True)
    with pytest.raises(ValueError, match="inside the input tree"):
        normalize_tree(source, source / "normalized")
