#!/usr/bin/env python3
"""Rebuild a checksum-locked EOS390 source cohort from the surviving v3.2 ZIP.

This does not recreate the lost tar.zst byte-for-byte. It extracts only the
original provider files listed in the checked-in inventory, hashes them, and
creates a new, explicitly versioned tar.zst. No normalization is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import stat
import subprocess
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any


DEFAULT_INVENTORY = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "eos390_reconstructed_source_inventory.v1.json"
)
SCRIPT_VERSION = "1.0.0"


class IntakeError(RuntimeError):
    """Raised when source identity or ZIP contents do not meet the inventory."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expected_files(inventory: dict[str, Any]) -> dict[str, int]:
    result: dict[str, int] = {}
    models = inventory.get("models")
    if not isinstance(models, dict) or len(models) != 18:
        raise IntakeError("inventory must declare exactly 18 EOS models")
    for model, files in models.items():
        if not isinstance(model, str) or not model or "/" in model or "\\" in model:
            raise IntakeError(f"invalid model directory name: {model!r}")
        if not isinstance(files, list):
            raise IntakeError(f"inventory model {model} must contain a file list")
        for record in files:
            name = record.get("name")
            size = record.get("size_bytes")
            if (
                not isinstance(name, str)
                or not name
                or "/" in name
                or "\\" in name
                or name in {".", ".."}
                or not isinstance(size, int)
                or size < 0
            ):
                raise IntakeError(f"invalid file record under {model}: {record!r}")
            key = f"{model}/{name}"
            if key in result:
                raise IntakeError(f"duplicate inventory path: {key}")
            result[key] = size
    if len(result) != inventory.get("counts", {}).get("files"):
        raise IntakeError("inventory file count does not match its counts.files field")
    return result


def _validate_member_name(name: str) -> PurePosixPath:
    if "\\" in name:
        raise IntakeError(f"ZIP member uses a backslash path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts:
        raise IntakeError(f"unsafe ZIP member path: {name!r}")
    return path


def collect_source_members(
    archive: zipfile.ZipFile, expected: dict[str, int]
) -> dict[str, zipfile.ZipInfo]:
    """Find exactly one code_package/eos_data/<model>/<file> tree."""
    found: dict[str, zipfile.ZipInfo] = {}
    known_models = {key.split("/", 1)[0] for key in expected}

    for info in archive.infolist():
        path = _validate_member_name(info.filename)
        if info.is_dir():
            continue
        parts = path.parts
        roots = [
            index
            for index in range(len(parts) - 1)
            if parts[index] == "code_package" and parts[index + 1] == "eos_data"
        ]
        if not roots:
            continue
        if len(roots) != 1:
            raise IntakeError(f"ambiguous EOS source root in ZIP member: {info.filename}")
        tail = parts[roots[0] + 2 :]
        if not tail:
            continue
        model = tail[0]
        if model not in known_models:
            raise IntakeError(f"unexpected directory under code_package/eos_data: {model}")
        if len(tail) != 2:
            raise IntakeError(
                f"unexpected nested source path (inventory expects model/file): {info.filename}"
            )
        key = f"{model}/{tail[1]}"
        if key not in expected:
            raise IntakeError(f"source member absent from checked-in inventory: {key}")
        if stat.S_ISLNK(info.external_attr >> 16):
            raise IntakeError(f"symbolic link is not accepted as a source file: {info.filename}")
        if info.file_size != expected[key]:
            raise IntakeError(
                f"size mismatch for {key}: ZIP={info.file_size}, inventory={expected[key]}"
            )
        if key in found:
            raise IntakeError(f"duplicate source member in ZIP: {key}")
        found[key] = info

    missing = sorted(set(expected) - set(found))
    if missing:
        raise IntakeError("ZIP is missing inventoried source files: " + ", ".join(missing))
    return found


def extract_source_members(
    archive: zipfile.ZipFile,
    members: dict[str, zipfile.ZipInfo],
    inventory: dict[str, Any],
    output_root: Path,
    archive_sha256: str,
) -> list[dict[str, Any]]:
    output_root.mkdir(parents=True, exist_ok=True)
    if any(output_root.iterdir()):
        raise IntakeError(f"refusing to write into non-empty output directory: {output_root}")

    records: list[dict[str, Any]] = []
    for key in sorted(members):
        info = members[key]
        target = output_root / key
        target.parent.mkdir(parents=True, exist_ok=True)
        digest = hashlib.sha256()
        byte_count = 0
        with archive.open(info, "r") as source, target.open("xb") as destination:
            for chunk in iter(lambda: source.read(1024 * 1024), b""):
                destination.write(chunk)
                digest.update(chunk)
                byte_count += len(chunk)
        if byte_count != info.file_size:
            raise IntakeError(f"streamed byte count mismatch for {key}")
        records.append(
            {
                "path": key,
                "source_zip_member": info.filename,
                "size_bytes": byte_count,
                "sha256": digest.hexdigest(),
            }
        )

    sums = "".join(f"{row['sha256']}  {row['path']}\n" for row in records)
    (output_root / "SHA256SUMS.txt").write_text(sums, encoding="utf-8", newline="\n")
    source_manifest = {
        "schema": "ns-cowling.eos390-reconstructed-source-manifest",
        "schema_version": "1.0.0",
        "status": "RECONSTRUCTED_FROM_SURVIVING_V3_2_ZIP",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "preparation_script_version": SCRIPT_VERSION,
        "source_archive": {
            "filename": inventory["replacement_archive"]["filename"],
            "size_bytes": inventory["replacement_archive"]["size_bytes"],
            "sha256": archive_sha256,
            "dropbox_path": inventory["replacement_archive"]["dropbox_path"],
            "dropbox_revision": inventory["replacement_archive"]["dropbox_revision"],
            "manifest_sha256_claim": inventory["replacement_archive"][
                "sha256_from_delivery_manifest"
            ],
        },
        "lost_archive": inventory["supersedes"],
        "identity_note": (
            "This is a new reconstruction from the surviving ZIP source subtree. "
            "It is not byte-identical to the lost tar.zst."
        ),
        "models": sorted(inventory["models"]),
        "file_count": len(records),
        "files": records,
    }
    (output_root / "source_manifest.json").write_text(
        json.dumps(source_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return records


def create_tar_zst(output_root: Path, output_archive: Path, model_names: list[str]) -> None:
    if output_archive.exists():
        raise IntakeError(f"refusing to overwrite existing archive: {output_archive}")
    if output_root.resolve() in output_archive.resolve().parents:
        raise IntakeError("output tar.zst must be outside the extracted source directory")
    output_archive.parent.mkdir(parents=True, exist_ok=True)
    staged = output_archive.with_name(output_archive.name + ".partial")
    members = ["SHA256SUMS.txt", "source_manifest.json", *sorted(model_names)]
    env = os.environ.copy()
    env["ZSTD_CLEVEL"] = "19"
    try:
        subprocess.run(
            [
                "tar",
                "--zstd",
                "--sort=name",
                "--mtime=@0",
                "--owner=0",
                "--group=0",
                "--numeric-owner",
                "-cf",
                str(staged),
                "-C",
                str(output_root),
                *members,
            ],
            check=True,
            env=env,
        )
        subprocess.run(
            ["tar", "--zstd", "-tf", str(staged)],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        os.replace(staged, output_archive)
    except (OSError, subprocess.CalledProcessError) as exc:
        staged.unlink(missing_ok=True)
        raise IntakeError(
            "could not create or verify tar.zst; GNU tar with zstd support is required"
        ) from exc


def prepare(
    archive_path: Path,
    inventory_path: Path,
    output_root: Path,
    output_archive: Path,
) -> dict[str, Any]:
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    expected = expected_files(inventory)
    archive_meta = inventory["replacement_archive"]
    actual_size = archive_path.stat().st_size
    if actual_size != archive_meta["size_bytes"]:
        raise IntakeError(
            f"archive size mismatch: got {actual_size}, expected {archive_meta['size_bytes']}"
        )
    actual_sha = sha256_file(archive_path)
    expected_sha = archive_meta["sha256_from_delivery_manifest"]
    if actual_sha != expected_sha:
        raise IntakeError(
            f"archive SHA-256 mismatch: got {actual_sha}, expected manifest {expected_sha}"
        )

    with zipfile.ZipFile(archive_path, "r") as zf:
        if len(zf.infolist()) != archive_meta["zip_members_reported"]:
            raise IntakeError(
                f"ZIP member count mismatch: got {len(zf.infolist())}, "
                f"expected {archive_meta['zip_members_reported']}"
            )
        corrupt_member = zf.testzip()
        if corrupt_member is not None:
            raise IntakeError(f"ZIP CRC failure in member: {corrupt_member}")
        members = collect_source_members(zf, expected)
        if len(members) != inventory["counts"]["files"]:
            raise IntakeError("selected source member count does not match inventory")
        records = extract_source_members(zf, members, inventory, output_root, actual_sha)

    create_tar_zst(output_root, output_archive, sorted(inventory["models"]))
    return {
        "status": "RECONSTRUCTED",
        "archive": str(output_archive),
        "archive_sha256": sha256_file(output_archive),
        "source_archive_sha256": actual_sha,
        "models": len(inventory["models"]),
        "files": len(records),
        "source_bytes": sum(row["size_bytes"] for row in records),
        "output_root": str(output_root),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--output-archive", type=Path)
    args = parser.parse_args()
    output_archive = args.output_archive or (
        args.output_root.parent / "NS_Cowling_EOS390_reconstructed_inputs_v2.tar.zst"
    )
    try:
        result = prepare(args.archive, args.inventory, args.output_root, output_archive)
    except (OSError, ValueError, KeyError, zipfile.BadZipFile, IntakeError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
