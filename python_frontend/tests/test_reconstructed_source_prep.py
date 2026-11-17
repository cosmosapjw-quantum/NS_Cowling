from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from prepare_eos390_reconstructed_inputs import (  # noqa: E402
    IntakeError,
    collect_source_members,
    expected_files,
    extract_source_members,
)


def _inventory():
    models = {
        "MODEL_A": [
            {"name": "eos.nb", "size_bytes": 3},
            {"name": "eos.thermo", "size_bytes": 4},
        ]
    }
    models.update({f"EMPTY_{index:02d}": [] for index in range(17)})
    return {
        "counts": {"files": 2},
        "models": models,
        "replacement_archive": {
            "filename": "source.zip",
            "size_bytes": 123,
            "sha256_from_delivery_manifest": "0" * 64,
            "dropbox_path": "/source.zip",
            "dropbox_revision": "rev-test",
        },
        "supersedes": {"lost_filename": "lost.tar.zst", "status": "LOST_USER_CONFIRMED"},
    }


def _zip(path: Path, members: dict[str, bytes]):
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, content in members.items():
            zf.writestr(name, content)


def test_collects_exact_source_subtree_and_ignores_other_eos_tree(tmp_path):
    path = tmp_path / "source.zip"
    _zip(
        path,
        {
            "release/code_package/eos_data/MODEL_A/eos.nb": b"abc",
            "release/code_package/eos_data/MODEL_A/eos.thermo": b"1234",
            "release/code_package/eos_data_normalized/MODEL_A/eos.nb": b"not-source",
        },
    )
    with zipfile.ZipFile(path) as zf:
        actual = collect_source_members(zf, expected_files(_inventory()))
    assert set(actual) == {"MODEL_A/eos.nb", "MODEL_A/eos.thermo"}


def test_duplicate_source_member_fails_closed(tmp_path):
    path = tmp_path / "duplicate.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("a/code_package/eos_data/MODEL_A/eos.nb", b"abc")
        zf.writestr("b/code_package/eos_data/MODEL_A/eos.nb", b"abc")
        zf.writestr("a/code_package/eos_data/MODEL_A/eos.thermo", b"1234")
    with zipfile.ZipFile(path) as zf:
        with pytest.raises(IntakeError, match="duplicate source member"):
            collect_source_members(zf, expected_files(_inventory()))


def test_missing_and_unlisted_members_fail_closed(tmp_path):
    missing = tmp_path / "missing.zip"
    _zip(missing, {"x/code_package/eos_data/MODEL_A/eos.nb": b"abc"})
    with zipfile.ZipFile(missing) as zf:
        with pytest.raises(IntakeError, match="missing inventoried"):
            collect_source_members(zf, expected_files(_inventory()))

    extra = tmp_path / "extra.zip"
    _zip(
        extra,
        {
            "x/code_package/eos_data/MODEL_A/eos.nb": b"abc",
            "x/code_package/eos_data/MODEL_A/eos.thermo": b"1234",
            "x/code_package/eos_data/MODEL_A/unlisted": b"x",
        },
    )
    with zipfile.ZipFile(extra) as zf:
        with pytest.raises(IntakeError, match="absent from checked-in inventory"):
            collect_source_members(zf, expected_files(_inventory()))


def test_size_mismatch_and_path_traversal_fail_closed(tmp_path):
    wrong_size = tmp_path / "wrong-size.zip"
    _zip(
        wrong_size,
        {
            "x/code_package/eos_data/MODEL_A/eos.nb": b"abcd",
            "x/code_package/eos_data/MODEL_A/eos.thermo": b"1234",
        },
    )
    with zipfile.ZipFile(wrong_size) as zf:
        with pytest.raises(IntakeError, match="size mismatch"):
            collect_source_members(zf, expected_files(_inventory()))

    traversal = tmp_path / "traversal.zip"
    _zip(
        traversal,
        {
            "../escape.txt": b"bad",
            "x/code_package/eos_data/MODEL_A/eos.nb": b"abc",
            "x/code_package/eos_data/MODEL_A/eos.thermo": b"1234",
        },
    )
    with zipfile.ZipFile(traversal) as zf:
        with pytest.raises(IntakeError, match="unsafe ZIP member"):
            collect_source_members(zf, expected_files(_inventory()))


def test_extract_writes_fresh_hash_manifest_and_never_normalizes(tmp_path):
    path = tmp_path / "source.zip"
    source_a = b"abc"
    source_b = b"1234"
    _zip(
        path,
        {
            "x/code_package/eos_data/MODEL_A/eos.nb": source_a,
            "x/code_package/eos_data/MODEL_A/eos.thermo": source_b,
        },
    )
    with zipfile.ZipFile(path) as zf:
        members = collect_source_members(zf, expected_files(_inventory()))
        records = extract_source_members(
            zf, members, _inventory(), tmp_path / "out", archive_sha256="f" * 64
        )
    assert (tmp_path / "out/MODEL_A/eos.nb").read_bytes() == source_a
    assert (tmp_path / "out/MODEL_A/eos.thermo").read_bytes() == source_b
    sums = (tmp_path / "out/SHA256SUMS.txt").read_text()
    assert len(sums.splitlines()) == 2
    assert all(len(row["sha256"]) == 64 for row in records)
    receipt = json.loads((tmp_path / "out/source_manifest.json").read_text())
    assert receipt["status"] == "RECONSTRUCTED_FROM_SURVIVING_V3_2_ZIP"
    assert receipt["identity_note"].startswith("This is a new reconstruction")
