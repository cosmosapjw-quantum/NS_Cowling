from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import sharp_maxwell_source_probe as probe


def test_git_blob_identity_matches_git_hash_object(tmp_path: Path) -> None:
    source = tmp_path / "source.dat"
    source.write_bytes(b"immutable\x00source\n")
    expected = subprocess.run(
        ["git", "hash-object", str(source)], capture_output=True, text=True, check=True
    ).stdout.strip()
    assert probe._fingerprint(source)["git_blob_sha1"] == expected


def test_raw_nb_and_thermo_are_read_rowwise(tmp_path: Path) -> None:
    nb = tmp_path / "eos.nb"
    nb.write_text("1\n2\n0.1\n0.2\n", encoding="utf-8")
    assert probe._read_raw_density(nb) == (1, [0.1, 0.2])

    thermo = tmp_path / "eos.thermo"
    thermo.write_text(
        "938.9 938.9 1\n"
        "0 1 1 1 2 3 4 5 6 7 0\n"
        "0 2 1 8 9 10 11 12 13 14 1 2\n",
        encoding="utf-8",
    )
    nadd, q_rows, extras, header = probe._read_raw_thermo(thermo)
    assert nadd == [0, 1]
    assert q_rows == [[1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14]]
    assert extras == [0.0, 2.0]
    assert header == (938.9, 938.9, 1)


def test_raw_composition_preserves_phase_and_compose_group_fraction(tmp_path: Path) -> None:
    compo = tmp_path / "eos.compo"
    compo.write_text(
        "0 1 1 0 2 10 0.2 0 0.1 1 1 56 26 0.5\n",
        encoding="utf-8",
    )
    phase, fields = probe._read_raw_composition(compo)
    assert phase == [0]
    assert fields["Yn"] == [0.2]
    assert fields["Ye"] == [0.1]
    assert fields["Yh"] == [0.5]


def test_nuclear_group_fraction_is_not_multiplied_by_average_mass(tmp_path: Path) -> None:
    # The selected CompOSE crusts include A_av=56, Y_I=1 rows.
    compo = tmp_path / "eos.compo"
    compo.write_text("0 1 1 0 2 10 0 11 0 1 1 56 26 1\n", encoding="utf-8")
    _, fields = probe._read_raw_composition(compo)
    assert fields["Aav"] == [56.0]
    assert fields["Yh"] == [1.0]


def test_histograms_and_phase_boundary_use_compose_knot_indices() -> None:
    assert probe._histogram([0, 7, 7, 0]) == {"0": 2, "7": 2}
    assert probe._phase_boundaries([0, 0, 7, 7]) == [
        {"left_knot": 2, "right_knot": 3, "from": 0, "to": 7}
    ]


def test_extra_column_is_not_silently_coerced() -> None:
    with pytest.raises(ValueError, match="non-integer-valued"):
        probe._extra_histogram([0.0, 1.25])
