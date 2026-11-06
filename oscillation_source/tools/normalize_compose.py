#!/usr/bin/env python3
"""Strict CompOSE cold-table intake and SI normalization for T13.

The parser intentionally preserves source row order and duplicate index triples.
It never repairs a table by sorting, deduplicating, clipping, or extrapolating.
"""
from __future__ import annotations

from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any
from uuid import uuid4

MEV_FM3_TO_J_M3 = 1.602176634e32  # exact from SI eV and fm definitions
FM3_INV_TO_M3_INV = 1.0e45
MEV_TO_J = 1.602176634e-13


class ComposeFormatError(ValueError):
    """Raised when a CompOSE source file violates the declared input contract."""


def _tokens(path: Path) -> list[str]:
    try:
        return Path(path).read_text(encoding="ascii").split()
    except (OSError, UnicodeError) as exc:
        raise ComposeFormatError(f"cannot read {path}: {exc}") from exc


def parse_nb(path: Path) -> dict[str, Any]:
    tok = _tokens(path)
    if len(tok) < 3:
        raise ComposeFormatError("eos.nb is too short")
    try:
        grid_code = int(tok[0])
        max_index = int(tok[1])
        values = [float(x) for x in tok[2:]]
    except ValueError as exc:
        raise ComposeFormatError("invalid eos.nb numeric token") from exc
    if max_index < 0 or len(values) != max_index + 1:
        raise ComposeFormatError(
            f"eos.nb grid count mismatch: max_index={max_index}, values={len(values)}"
        )
    if any(not (x > 0.0) for x in values):
        raise ComposeFormatError("eos.nb requires positive baryon densities")
    # Grid files define the index mapping. Do not silently sort a malformed grid.
    if any(b <= a for a, b in zip(values, values[1:])):
        raise ComposeFormatError("eos.nb is not strictly increasing in source order")
    return {"grid_code": grid_code, "max_index": max_index, "n_b_fm3": values}


def parse_thermo(path: Path) -> dict[str, Any]:
    lines = [line for line in Path(path).read_text(encoding="ascii").splitlines() if line.strip()]
    if len(lines) < 2:
        raise ComposeFormatError("eos.thermo is too short")
    head = lines[0].split()
    if len(head) != 3:
        raise ComposeFormatError("eos.thermo header must have neutron mass, proton mass, lepton flag")
    try:
        m_n_mev = float(head[0])
        m_p_mev = float(head[1])
        lepton_flag = int(head[2])
    except ValueError as exc:
        raise ComposeFormatError("invalid eos.thermo header") from exc
    if m_n_mev <= 0 or m_p_mev <= 0 or lepton_flag not in (0, 1):
        raise ComposeFormatError("invalid eos.thermo mass/lepton convention")

    rows: list[dict[str, Any]] = []
    for source_row, line in enumerate(lines[1:], start=2):
        fields = line.split()
        if len(fields) < 11:
            raise ComposeFormatError(f"thermo row {source_row} has fewer than 11 required entries")
        try:
            i_t, i_nb, i_yq = map(int, fields[:3])
            q = [float(x) for x in fields[3:10]]
            n_add = int(fields[10])
            additional = [float(x) for x in fields[11:]]
        except ValueError as exc:
            raise ComposeFormatError(f"invalid thermo row {source_row}") from exc
        if n_add < 0 or len(additional) != n_add:
            raise ComposeFormatError(
                f"thermo row {source_row} additional-quantity count mismatch"
            )
        rows.append(
            {
                "source_row": source_row,
                "i_t": i_t,
                "i_nb": i_nb,
                "i_yq": i_yq,
                "q1_p_over_nb_mev": q[0],
                "q2_entropy_per_baryon": q[1],
                "q3_mu_b_over_mn_minus_1": q[2],
                "q4_mu_q_over_mn": q[3],
                "q5_mu_l_over_mn": q[4],
                "q6_f_over_nb_mn_minus_1": q[5],
                "q7_e_over_nb_mn_minus_1": q[6],
                "additional": additional,
            }
        )
    return {
        "m_n_mev": m_n_mev,
        "m_p_mev": m_p_mev,
        "leptons_included": bool(lepton_flag),
        "rows": rows,
    }


def normalize_cold_table(nb_path: Path, thermo_path: Path) -> dict[str, Any]:
    nb = parse_nb(nb_path)
    thermo = parse_thermo(thermo_path)
    n_grid = nb["n_b_fm3"]
    normalized: list[dict[str, Any]] = []
    triples: list[tuple[int, int, int]] = []

    for row in thermo["rows"]:
        i_nb = row["i_nb"]
        if i_nb < 0 or i_nb >= len(n_grid):
            raise ComposeFormatError(
                f"density index {i_nb} outside eos.nb range 0..{len(n_grid)-1}"
            )
        n_fm3 = n_grid[i_nb]
        q1 = row["q1_p_over_nb_mev"]
        q3 = row["q3_mu_b_over_mn_minus_1"]
        q7 = row["q7_e_over_nb_mn_minus_1"]
        pressure_mev_fm3 = n_fm3 * q1
        energy_mev_fm3 = n_fm3 * thermo["m_n_mev"] * (1.0 + q7)
        mu_b_mev = thermo["m_n_mev"] * (1.0 + q3)
        triple = (row["i_t"], i_nb, row["i_yq"])
        triples.append(triple)
        normalized.append(
            {
                **row,
                "n_b_fm3": n_fm3,
                "n_b_m3": n_fm3 * FM3_INV_TO_M3_INV,
                "pressure_mev_fm3": pressure_mev_fm3,
                "pressure_pa": pressure_mev_fm3 * MEV_FM3_TO_J_M3,
                "energy_mev_fm3": energy_mev_fm3,
                "energy_j_m3": energy_mev_fm3 * MEV_FM3_TO_J_M3,
                "mu_b_mev": mu_b_mev,
                "mu_b_j": mu_b_mev * MEV_TO_J,
            }
        )

    counts = Counter(triples)
    duplicates = [list(t) for t, count in counts.items() if count > 1]
    mu_values = [r["mu_b_mev"] for r in normalized]
    mu_nonmono = any(b < a for a, b in zip(mu_values, mu_values[1:]))
    return {
        "source_rows_preserved": True,
        "normalized_units": "SI",
        "energy_includes_rest_mass": True,
        "m_n_mev": thermo["m_n_mev"],
        "m_p_mev": thermo["m_p_mev"],
        "leptons_included": thermo["leptons_included"],
        "duplicate_index_triples": duplicates,
        "mu_b_nonmonotone_in_source_order": mu_nonmono,
        "rows": normalized,
    }


NORMALIZED_FIELDS = [
    "source_row", "i_t", "i_nb", "i_yq", "n_b_fm3", "n_b_m3",
    "pressure_mev_fm3", "pressure_pa", "energy_mev_fm3", "energy_j_m3",
    "mu_b_mev", "mu_b_j", "q1_p_over_nb_mev", "q2_entropy_per_baryon",
    "q3_mu_b_over_mn_minus_1", "q4_mu_q_over_mn", "q5_mu_l_over_mn",
    "q6_f_over_nb_mn_minus_1", "q7_e_over_nb_mn_minus_1",
]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _load_metadata(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ComposeFormatError(f"cannot read metadata {path}: {exc}") from exc
    if not isinstance(value, dict) or not isinstance(value.get("member_sha256"), dict):
        raise ComposeFormatError("metadata lacks member_sha256 authority")
    return value


def _verify_raw_members(metadata: dict[str, Any], raw_dir: Path) -> dict[str, str]:
    raw_dir = Path(raw_dir)
    verified: dict[str, str] = {}
    for name, expected in metadata["member_sha256"].items():
        path = raw_dir / name
        if not path.is_file() or path.is_symlink():
            raise ComposeFormatError(f"missing raw member: {name}")
        actual = _sha256(path)
        if actual != expected:
            raise ComposeFormatError(
                f"raw member hash mismatch for {name}: expected {expected}, got {actual}"
            )
        verified[name] = actual
    return verified


def _samples_match(rows: list[dict[str, Any]], samples: list[dict[str, Any]]) -> bool:
    for sample in samples:
        try:
            row = rows[int(sample["zero_based_row"])]
        except (KeyError, TypeError, ValueError, IndexError):
            return False
        for key in ("n_b_fm3", "pressure_pa", "energy_j_m3", "mu_b_mev"):
            if key not in sample or not math.isclose(
                float(row[key]), float(sample[key]), rel_tol=5e-14, abs_tol=0.0
            ):
                return False
    return bool(samples)


def _write_csv_atomic(path: Path, rows: list[dict[str, Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_symlink():
        raise ComposeFormatError("normalized CSV destination is a symlink")
    tmp = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        with tmp.open("x", encoding="ascii", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=NORMALIZED_FIELDS, extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row[key] for key in NORMALIZED_FIELDS})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def _write_json_atomic(path: Path, value: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_symlink():
        raise ComposeFormatError("receipt destination is a symlink")
    raw = (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode()
    tmp = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        with tmp.open("xb") as handle:
            handle.write(raw)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def write_normalized_dataset(
    metadata_path: Path, raw_dir: Path, out_csv: Path, out_receipt: Path
) -> dict[str, Any]:
    """Verify pinned raw members, normalize eos.nb/eos.thermo, and write evidence.

    This is an intake transform only. It does not interpolate, repair, infer hidden
    phase boundaries, or admit the table for stellar/mode calculations.
    """
    metadata = _load_metadata(metadata_path)
    verified = _verify_raw_members(metadata, raw_dir)
    result = normalize_cold_table(Path(raw_dir) / "eos.nb", Path(raw_dir) / "eos.thermo")
    rows = result["rows"]
    expected_count = metadata.get("nb_grid_points")
    if expected_count is not None and len(rows) != int(expected_count):
        raise ComposeFormatError(
            f"row count {len(rows)} does not match metadata nb_grid_points={expected_count}"
        )
    samples = metadata.get("independent_sample_rows", [])
    sample_ok = _samples_match(rows, samples)
    if not sample_ok:
        raise ComposeFormatError("independent sample-row recovery failed")

    _write_csv_atomic(out_csv, rows)
    receipt: dict[str, Any] = {
        "schema_version": 1,
        "task": "T13",
        "status": "ok",
        "compose_id": metadata.get("compose_id"),
        "model_name": metadata.get("name"),
        "raw_sha256_verified": bool(metadata.get("raw_sha256_verified")),
        "verified_member_sha256": verified,
        "observed_archive_sha256": metadata.get("observed_archive_sha256"),
        "official_checksum_declared_sha256": metadata.get("official_checksum_declared_sha256"),
        "official_checksum_matches_observed": metadata.get("official_checksum_matches_observed"),
        "checksum_discrepancy_preserved": metadata.get("checksum_discrepancy_preserved"),
        "distribution_status": metadata.get("distribution_status"),
        "row_count": len(rows),
        "normalized_units": result["normalized_units"],
        "energy_includes_rest_mass": result["energy_includes_rest_mass"],
        "source_rows_preserved": result["source_rows_preserved"],
        "duplicate_index_triples": result["duplicate_index_triples"],
        "mu_b_nonmonotone_in_source_order": result["mu_b_nonmonotone_in_source_order"],
        "segment_count": 1,
        "segment_policy": "one source-order segment; no physical jump inferred from spacing alone",
        "physical_jumps": [],
        "jump_detection_scope": "explicit source rows only; derivative/physics admission deferred to T14",
        "sample_rows_independently_recovered": sample_ok,
        "independent_sample_rows": samples,
        "normalized_csv_sha256": _sha256(out_csv),
        "scientific_admission": "T13_INTAKE_ONLY_NOT_T14_TABLE_EOS_ADMISSION",
    }
    _write_json_atomic(out_receipt, receipt)
    return receipt


def main(argv=None) -> int:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--raw-dir', type=Path, required=True)
    parser.add_argument('--out-csv', type=Path, required=True)
    parser.add_argument('--out-receipt', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        receipt = write_normalized_dataset(
            args.metadata, args.raw_dir, args.out_csv, args.out_receipt
        )
    except ComposeFormatError as exc:
        parser.exit(2, f'{exc}\n')
    print(json.dumps({'status': receipt['status'], 'compose_id': receipt['compose_id'],
                      'row_count': receipt['row_count'],
                      'normalized_csv_sha256': receipt['normalized_csv_sha256']},
                     sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
