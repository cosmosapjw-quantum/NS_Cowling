"""E390-10B phase 5: lossless multidim round-trip (read -> write -> read).

Phase 5 extends the writer to emit the full nT x nYq x n tensor and verifies
that a multidim table round-trips bit-for-bit through the internal directory:
  read_compose(dir) -> write_internal_directory(tmp) -> from_internal_directory(tmp)
preserves every (T,Yq,density) cell in q_full and q_extra_full.

Contract verified here:
  1. A 3T x 2Yq x n Nadd=1 table round-trips: every q_full cell and every
     q_extra_full cell matches the original read to machine precision.
  2. A 3T x 1Yq x n Nadd=0 table (multidim T only) round-trips.
  3. A 1T x 3Yq x n Nadd=2 table (multidim Yq only) round-trips.
  4. The cold slice (1T x 1Yq) still round-trips (no tensor allocated).
"""
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import read_compose, write_compose  # noqa: E402


def _write_multi_t_table(d: Path, *, n: int, nT: int, nYq: int = 1, nadd: int = 0):
    d.mkdir(parents=True, exist_ok=True)
    nb = np.linspace(0.05, 0.05 + 0.001 * (n - 1), n)
    (d / "eos.nb").write_text(f"1\n{n}\n" + "\n".join(f"{x:.10e}" for x in nb) + "\n")
    t_vals = [0.0, 10.0, 30.0, 50.0, 100.0][:nT]
    (d / "eos.t").write_text(f"1\n{nT}\n" + "\n".join(f"{t:.10e}" for t in t_vals) + "\n")
    yq_vals = [0.0, 0.1, 0.3, 0.5][:nYq]
    (d / "eos.yq").write_text(f"1\n{nYq}\n" + "\n".join(f"{y:.10e}" for y in yq_vals) + "\n")
    lines = ["939.0 939.0 1"]
    for it in range(1, nT + 1):
        for iy in range(1, nYq + 1):
            for i in range(n):
                q1 = 1.0e-3 + 1.0e-5 * (i + 10 * (iy - 1) + 100 * (it - 1))
                row = f"{it} {i+1} {iy} {q1:.10e} 0 0 0 0 0 0"
                if nadd == 0:
                    row += " 0"
                else:
                    row += f" {nadd}"
                    for slot in range(1, nadd + 1):
                        row += f" {100*slot + i + 10*(iy-1) + 100*(it-1):.10e}"
                lines.append(row)
    (d / "eos.thermo").write_text("\n".join(lines) + "\n")


# --- 1. 3T x 2Yq x n Nadd=1 round-trips bit-for-bit ----------------------

def test_3t_2yq_nadd1_round_trips_bit_for_bit(tmp_path):
    _write_multi_t_table(tmp_path / "src", n=4, nT=3, nYq=2, nadd=1)
    eos1 = read_compose(tmp_path / "src")
    # write internal directory
    outdir = tmp_path / "roundtrip"
    eos1.write_internal_directory(outdir)
    # read back
    from compose_eos import EOSData
    eos2 = EOSData.from_internal_directory(outdir)
    assert eos2.nT == 3 and eos2.nYq == 2
    assert eos2.q_full.shape == (3, 2, 7, 4)
    assert eos2.q_extra_full.shape == (1, 3, 2, 4)
    # every q_full cell matches
    np.testing.assert_allclose(eos2.q_full, eos1.q_full, rtol=0, atol=0)
    # every q_extra_full cell matches
    np.testing.assert_allclose(eos2.q_extra_full, eos1.q_extra_full, rtol=0, atol=0)
    # 1D surface matches the first (T,Yq) cell
    np.testing.assert_allclose(eos2.Q1, eos1.Q1, rtol=0, atol=0)
    np.testing.assert_allclose(eos2.q_extra, eos1.q_extra, rtol=0, atol=0)


# --- 2. 3T x 1Yq x n Nadd=0 (multidim T only) round-trips ---------------

def test_3t_1yq_nadd0_round_trips(tmp_path):
    _write_multi_t_table(tmp_path / "src", n=3, nT=3, nYq=1, nadd=0)
    eos1 = read_compose(tmp_path / "src")
    outdir = tmp_path / "roundtrip"
    eos1.write_internal_directory(outdir)
    from compose_eos import EOSData
    eos2 = EOSData.from_internal_directory(outdir)
    assert eos2.nT == 3 and eos2.nYq == 1
    assert eos2.q_full.shape == (3, 1, 7, 3)
    np.testing.assert_allclose(eos2.q_full, eos1.q_full, rtol=0, atol=0)


# --- 3. 1T x 3Yq x n Nadd=2 (multidim Yq only) round-trips -------------

def test_1t_3yq_nadd2_round_trips(tmp_path):
    _write_multi_t_table(tmp_path / "src", n=4, nT=1, nYq=3, nadd=2)
    eos1 = read_compose(tmp_path / "src")
    outdir = tmp_path / "roundtrip"
    eos1.write_internal_directory(outdir)
    from compose_eos import EOSData
    eos2 = EOSData.from_internal_directory(outdir)
    assert eos2.nT == 1 and eos2.nYq == 3
    assert eos2.q_full.shape == (1, 3, 7, 4)
    assert eos2.q_extra_full.shape == (2, 1, 3, 4)
    np.testing.assert_allclose(eos2.q_full, eos1.q_full, rtol=0, atol=0)
    np.testing.assert_allclose(eos2.q_extra_full, eos1.q_extra_full, rtol=0, atol=0)


# --- 4. cold slice (1T x 1Yq) still round-trips (no tensor) -----------

def test_cold_slice_round_trips_after_phase5(tmp_path):
    _write_multi_t_table(tmp_path / "src", n=4, nT=1, nYq=1, nadd=0)
    eos1 = read_compose(tmp_path / "src")
    assert eos1.nT == 1 and eos1.nYq == 1
    assert eos1.q_full.size == 0  # cold slice: no tensor
    outdir = tmp_path / "roundtrip"
    eos1.write_internal_directory(outdir)
    from compose_eos import EOSData
    eos2 = EOSData.from_internal_directory(outdir)
    assert eos2.nT == 1 and eos2.nYq == 1
    assert eos2.q_full.size == 0
    np.testing.assert_allclose(eos2.Q1, eos1.Q1, rtol=0, atol=0)
