"""Lossless archive transport and hostile archive rejection; no EOS numerics."""
import gzip
import hashlib
import importlib
import io
import json
from pathlib import Path
import sys
import tarfile

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
package = importlib.import_module("package_eos_reaudit_20261123")
extract = importlib.import_module("extract_eos_reaudit_candidate_20261123")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def tar_bytes(files):
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode="w", format=tarfile.USTAR_FORMAT) as tf:
        for name, value in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(value)
            info.mode = 0o644
            tf.addfile(info, io.BytesIO(value))
    return stream.getvalue()


def campaign_fixture(path, public=3, private=1):
    path.mkdir()
    (path / "artifacts").mkdir()
    records, contents = [], {}
    for i in range(public + private):
        cell = f"cell_{i:06d}"
        r = dict(cell_id=cell, status="SURVIVED", private=i >= public,
                 crust_source_id="test:crust", core_source_id="test:core",
                 method="maxwell" if i == 1 else "quintic_P",
                 source_coverage="FINITE_PHASE_COMPONENT_ONLY")
        data = {"PROVENANCE.json": json.dumps(r).encode(),
                "compose/eos.thermo": b"939 938 1\n1 1 1 2.0 0 0 0 0 0 0\n",
                "internal/eos_data.csv": b"n,P,epsilon\n0.1,2.0,94.0\n",
                "opaque/extra.bin": b"\x00\xff full original candidate payload \n"}
        raw = tar_bytes(data)
        compressed = gzip.compress(raw, mtime=0)
        archive = path / "artifacts" / (cell + ".tar.gz")
        if not r["private"]:
            archive.write_bytes(compressed)
        else:
            r["native_report"] = {"secret": "PRIVATE_NUMERICAL_PAYLOAD_MUST_NOT_PUBLISH"}
            r["failure_trace"] = {"path": "PRIVATE_TRACE_MUST_NOT_PUBLISH"}
        r["artifact"] = dict(path=str(archive), sha256=sha(compressed), bytes=len(compressed))
        records.append(r)
        contents[cell] = (data, raw, compressed)
    settings = dict(pilot=False, full_requested_cells=len(records), catalog_sha256="c"*64,
                    binary_sha256="d"*64, runner_sha256="e"*64)
    summary = dict(status="COMPLETED", input_and_binary_hashes_unchanged=True,
                   requested_cells=len(records), recorded_cells=len(records),
                   survivors=len(records), statuses={"SURVIVED": len(records)}, settings=settings)
    for name, value in [("SUMMARY.json", summary), ("SURVIVORS.json", records),
                        ("SETTINGS.json", settings), ("INPUT_IDENTITIES.json", {})]:
        (path / name).write_text(json.dumps(value))
    (path / "RESULTS.jsonl").write_text("".join(json.dumps(r)+"\n" for r in records))
    return records, contents


def test_full_tar_roundtrip_preserves_all_files_and_excludes_private(tmp_path):
    campaign = tmp_path / "campaign"
    records, contents = campaign_fixture(campaign)
    before = {str(p): sha(p.read_bytes()) for p in campaign.rglob("*") if p.is_file()}
    delivery = tmp_path / "delivery"
    index = package.package_campaign(campaign, delivery, shard_cells=128, workers=1, xz_preset=1)
    assert index["status"] == "COMPLETED_LOSSLESS_PUBLIC_DELIVERY"
    assert set(index["cells"]) == {r["cell_id"] for r in records if not r["private"]}
    assert set(index["private_cells"]) == {records[-1]["cell_id"]}
    text = (delivery / "DELIVERY_INDEX.json").read_text()
    assert "PRIVATE_NUMERICAL_PAYLOAD" not in text and "PRIVATE_TRACE" not in text
    assert "artifact" not in index["private_cells"][records[-1]["cell_id"]]
    for cell, record in index["cells"].items():
        expected, raw, compressed = contents[cell]
        assert record["original_gzip_sha256"] == sha(compressed)
        assert record["raw_tar_sha256"] == sha(raw)
        with tarfile.open(delivery / record["shard"], "r:xz") as tf:
            assert tf.extractfile(record["member"]).read() == raw
        target = tmp_path / (cell + "_extracted")
        result = extract.extract_candidate(delivery / "DELIVERY_INDEX.json", cell, target)
        assert result["status"] == "EXTRACTED_HASH_VERIFIED"
        actual = {str(p.relative_to(target)): p.read_bytes() for p in target.rglob("*") if p.is_file()}
        assert actual == expected
    assert before == {str(p): sha(p.read_bytes()) for p in campaign.rglob("*") if p.is_file()}


@pytest.mark.parametrize("change", ["incomplete", "identity_changed", "count_mismatch"])
def test_packaging_requires_completed_consistent_campaign(tmp_path, change):
    campaign = tmp_path / "campaign"
    campaign_fixture(campaign)
    p = campaign / "SUMMARY.json"
    summary = json.loads(p.read_text())
    if change == "incomplete": summary["status"] = "RUNNING"
    if change == "identity_changed": summary["input_and_binary_hashes_unchanged"] = False
    if change == "count_mismatch": summary["recorded_cells"] -= 1
    p.write_text(json.dumps(summary))
    with pytest.raises(ValueError):
        package.package_campaign(campaign, tmp_path / "rejected", workers=1, xz_preset=1)
    assert not (tmp_path / "rejected").exists()


def test_original_gzip_hash_mismatch_blocks_publication(tmp_path):
    campaign = tmp_path / "campaign"
    campaign_fixture(campaign)
    (campaign / "artifacts/cell_000000.tar.gz").write_bytes(b"changed source")
    with pytest.raises(ValueError, match="hash"):
        package.package_campaign(campaign, tmp_path / "rejected", workers=1, xz_preset=1)
    assert not (tmp_path / "rejected").exists()


def test_summary_change_during_packaging_blocks_publication(tmp_path, monkeypatch):
    campaign = tmp_path / "campaign"
    campaign_fixture(campaign)
    original = package._build_shard
    def mutate(*args, **kwargs):
        result = original(*args, **kwargs)
        p = campaign / "SUMMARY.json"
        p.write_bytes(p.read_bytes() + b"\n")
        return result
    monkeypatch.setattr(package, "_build_shard", mutate)
    with pytest.raises(ValueError, match="changed"):
        package.package_campaign(campaign, tmp_path / "rejected", workers=1, xz_preset=1)
    assert not (tmp_path / "rejected").exists()


@pytest.mark.parametrize("which", ["shard", "member"])
def test_extractor_rejects_hash_mismatch_before_output(tmp_path, which):
    campaign = tmp_path / "campaign"
    campaign_fixture(campaign, public=1, private=0)
    delivery = tmp_path / "delivery"
    index = package.package_campaign(campaign, delivery, workers=1, xz_preset=1)
    entry = index["cells"]["cell_000000"]
    if which == "shard":
        p = delivery / entry["shard"]
        p.write_bytes(p.read_bytes() + b"corrupt")
    else:
        entry["raw_tar_sha256"] = "0"*64
        (delivery / "DELIVERY_INDEX.json").write_text(json.dumps(index))
    with pytest.raises(ValueError, match="hash"):
        extract.extract_candidate(delivery / "DELIVERY_INDEX.json", "cell_000000", tmp_path / "rejected")
    assert not (tmp_path / "rejected").exists()


def unsafe_delivery(path, name, kind="file"):
    path.mkdir()
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as tf:
        info = tarfile.TarInfo(name)
        if kind == "symlink":
            info.type = tarfile.SYMTYPE; info.linkname = "../escaped"
        elif kind == "hardlink":
            info.type = tarfile.LNKTYPE; info.linkname = "../escaped"
        elif kind == "fifo":
            info.type = tarfile.FIFOTYPE
        else:
            info.size = 1
        tf.addfile(info, io.BytesIO(b"x") if kind in ("file", "duplicate") else None)
        if kind == "duplicate":
            tf.addfile(info, io.BytesIO(b"x"))
    raw = buffer.getvalue()
    shard = path / "shard.tar.xz"
    with tarfile.open(shard, "w:xz", preset=1) as tf:
        member = tarfile.TarInfo("cell_000000.tar"); member.size = len(raw)
        tf.addfile(member, io.BytesIO(raw))
    index = dict(schema="ns_cowling.eos_delivery.v1", status="COMPLETED_LOSSLESS_PUBLIC_DELIVERY",
        cells={"cell_000000": dict(private=False, shard=shard.name, member="cell_000000.tar",
             raw_tar_bytes=len(raw), raw_tar_sha256=sha(raw))}, private_cells={},
        shards=[dict(path=shard.name, sha256=sha(shard.read_bytes()), bytes=shard.stat().st_size, cells=1)])
    p = path / "DELIVERY_INDEX.json"; p.write_text(json.dumps(index))
    return p


@pytest.mark.parametrize("name,kind", [
    ("../escaped", "file"), ("/absolute", "file"), ("C:/absolute", "file"),
    ("nested\\escaped", "file"), ("link", "symlink"), ("hard", "hardlink"),
    ("pipe", "fifo"), ("duplicate", "duplicate"),
])
def test_safe_extraction_rejects_paths_links_special_files_and_duplicates(tmp_path, name, kind):
    index = unsafe_delivery(tmp_path / "delivery", name, kind)
    with pytest.raises(ValueError):
        extract.extract_candidate(index, "cell_000000", tmp_path / "rejected")
    assert not (tmp_path / "rejected").exists()
    assert not (tmp_path / "escaped").exists()


def test_existing_output_is_preserved(tmp_path):
    campaign = tmp_path / "campaign"
    campaign_fixture(campaign, public=1, private=0)
    delivery = tmp_path / "delivery"
    package.package_campaign(campaign, delivery, workers=1, xz_preset=1)
    out = tmp_path / "existing"; out.mkdir(); (out / "keep").write_text("immutable")
    with pytest.raises(FileExistsError):
        extract.extract_candidate(delivery / "DELIVERY_INDEX.json", "cell_000000", out)
    with pytest.raises(FileExistsError):
        package.package_campaign(campaign, delivery, workers=1, xz_preset=1)
    assert (out / "keep").read_text() == "immutable"


def test_sha256_supports_declared_python310_without_file_digest(tmp_path, monkeypatch):
    path = tmp_path / "reference"
    path.write_bytes(b"abc")
    monkeypatch.delattr(hashlib, "file_digest", raising=False)
    assert extract.file_sha256(path) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
