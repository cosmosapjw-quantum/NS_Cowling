from __future__ import annotations

import hashlib
import math

import pytest

from scripts.r15_native_density_preflight import audit_native_table


def _write_table(path, *, count=40, maximum=1.0e-3):
    rows = [
        "source_row,pressure_pa,energy_mev_fm3,n_b_fm3,mu_b_mev"
    ]
    for i in range(count):
        fraction = i / (count - 1)
        n_b = 10.0 ** (-8.0 + 5.7 * fraction)
        p = 10.0 ** (20.0 + 5.5 * fraction)
        eps = max(1.0e-15, n_b * 930.0)
        rows.append(f"{i + 1},{p:.17e},{eps:.17e},{n_b:.17e},{930.0 + fraction:.17e}")
    rows.append(f"{count + 1},{1.0e28:.17e},{3.0:.17e},{3.0e-3:.17e},{941.0:.17e}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_audit_selects_table_local_native_rows_without_fixed_index_assumption(tmp_path):
    path = tmp_path / "variable_mesh.csv"
    sha = _write_table(path, count=41)
    record, knots = audit_native_table(
        path,
        name="synthetic_variable_mesh",
        expected_sha256=sha,
        central_pressure_pa=1.0e35,
    )
    assert len(knots) == 8
    assert knots[0]["source_row_one_based"] == 1
    assert knots[-1]["source_row_one_based"] == record["last_eligible_source_row"]
    assert len(set(record["selected_row_index_signature_zero_based"])) == 8
    assert record["unchanged"] is True
    assert record["boundary_to_first_ineligible"]["density_jump_ratio"] > 1.0


def test_audit_rejects_projected_table_hash_mismatch(tmp_path):
    path = tmp_path / "table.csv"
    sha = _write_table(path)
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        audit_native_table(
            path,
            name="hash_mismatch",
            expected_sha256="0" * 64,
            central_pressure_pa=1.0e35,
        )


def test_audit_rejects_source_row_reordering(tmp_path):
    path = tmp_path / "reordered.csv"
    _write_table(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    lines[2], lines[3] = lines[3], lines[2]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="source_row order mismatch"):
        audit_native_table(
            path,
            name="reordered",
            expected_sha256=sha,
            central_pressure_pa=1.0e35,
        )


def test_audit_rejects_cutoff_pressure_above_central_pressure(tmp_path):
    path = tmp_path / "pressure_domain.csv"
    sha = _write_table(path)
    with pytest.raises(ValueError, match="not below central pressure"):
        audit_native_table(
            path,
            name="bad_pc",
            expected_sha256=sha,
            central_pressure_pa=1.0,
        )
