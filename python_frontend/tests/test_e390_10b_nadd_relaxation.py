"""E390-10B phase 2: Nadd relaxation (reader accepts arbitrary Nadd>=0).

Phase 2 removes the hard nadd<=1 bound on the eos.thermo reader. The reader
now accepts any non-negative Nadd, stores every extra column in the
q_extra_full(nadd_max,nT,nYq,n) tensor, and keeps the 1D q_extra(n) surface
mirroring slot 1 for backward compatibility with M3-M6 downstream code.

Contract verified here:
  1. A synthetic Nadd=2 cold-slice table reads: nadd=2, nadd_max=2, and the
     q_extra_full tensor carries both extra columns (slot 1 == 1D q_extra).
  2. A synthetic Nadd=3 table reads and round-trips through write_compose
     (CompOSE format) bit-for-bit in both extra columns.
  3. Cold-slice Nadd in {0,1} is unaffected (covered by the existing suite:
     read_compose on RG_SLY4 [Nadd=1] and reconstruct_master [Nadd=0/1]).
  4. The internal CSV round-trips Nadd>=2 via the q_extra_2, q_extra_3, ...
     optional columns (absent columns default to 0).
"""
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import EOSData, read_compose, write_compose  # noqa: E402


def _write_nadd_table(d: Path, *, n: int, nadd: int, extra_vals=None):
    """Write a minimal CompOSE cold-slice dir with an arbitrary Nadd>=0.

    Each thermo row has width 11 + nadd (iT inb iYq Q1..Q7 Nadd [extra_1..]).
    extra_vals is a dict {slot_index (1-based): np.array(n)} of extra column
    values; slots default to distinct signature values 100*slot + k.
    """
    d.mkdir(parents=True, exist_ok=True)
    nb = np.linspace(0.05, 0.05 + 0.001 * (n - 1), n)
    (d / "eos.nb").write_text(f"1\n{n}\n" + "\n".join(f"{x:.10e}" for x in nb) + "\n")
    (d / "eos.t").write_text(f"1\n1\n0.0000000000e+00\n")
    (d / "eos.yq").write_text(f"1\n1\n0.0000000000e+00\n")
    lines = ["939.0 939.0 1"]
    for i in range(n):
        q1 = 1.0e-3 + 1.0e-5 * i
        # row = iT inb iYq (3) + Q1..Q7 (7) + Nadd (1) + extras (nadd)
        row = f"1 {i+1} 1 {q1:.10e} 0 0 0 0 0 0 {nadd}"
        for slot in range(1, nadd + 1):
            if extra_vals is not None and slot in extra_vals:
                v = float(extra_vals[slot][i])
            else:
                v = 100.0 * slot + i  # distinct signature per slot
            row += f" {v:.10e}"
        lines.append(row)
    (d / "eos.thermo").write_text("\n".join(lines) + "\n")


# --- 1. a synthetic Nadd=2 table reads into q_extra_full ---------------------

def test_nadd2_table_reads_into_q_extra_full(tmp_path):
    _write_nadd_table(tmp_path / "nadd2", n=6, nadd=2)
    eos = read_compose(tmp_path / "nadd2")
    assert eos.q_extra_count == 2, "nadd must be 2"
    assert eos.nadd_max == 2, "nadd_max must be 2"
    # slot 1 mirrors the 1D q_extra surface
    expected_slot1 = np.array([100, 101, 102, 103, 104, 105], dtype=float)
    expected_slot2 = np.array([200, 201, 202, 203, 204, 205], dtype=float)
    assert np.allclose(eos.q_extra, expected_slot1)
    # q_extra_full carries both slots (cold slice: nT=nYq=1)
    assert eos.q_extra_full.shape == (2, 1, 1, 6)
    assert np.allclose(eos.q_extra_full[0, 0, 0, :], expected_slot1)
    assert np.allclose(eos.q_extra_full[1, 0, 0, :], expected_slot2)
    # q_full stays empty (cold slice, no multidim reader)
    assert eos.q_full.size == 0


# --- 2. Nadd=3 reads and round-trips through write_compose ------------------

def test_nadd3_round_trips_through_compose_format(tmp_path):
    src = tmp_path / "nadd3_src"
    dst = tmp_path / "nadd3_dst"
    n = 5
    extra = {1: np.array([10, 11, 12, 13, 14], dtype=float),
             2: np.array([20, 21, 22, 23, 24], dtype=float),
             3: np.array([30, 31, 32, 33, 34], dtype=float)}
    _write_nadd_table(src, n=n, nadd=3, extra_vals=extra)
    eos = read_compose(src)
    assert eos.q_extra_count == 3 and eos.nadd_max == 3
    assert eos.q_extra_full.shape == (3, 1, 1, n)
    for slot in (1, 2, 3):
        assert np.allclose(eos.q_extra_full[slot - 1, 0, 0, :], extra[slot])
    # write back as CompOSE and re-read; every extra column must round-trip
    write_compose(dst, eos)
    back = read_compose(dst)
    assert back.q_extra_count == 3
    assert back.nadd_max == 3
    assert back.q_extra_full.shape == (3, 1, 1, n)
    for slot in (1, 2, 3):
        assert np.allclose(back.q_extra_full[slot - 1, 0, 0, :], extra[slot], atol=1e-9)


# --- 3. internal CSV round-trips Nadd>=2 via q_extra_N columns ----------------

def test_nadd2_internal_directory_round_trip(tmp_path):
    _write_nadd_table(tmp_path / "nadd2int", n=4, nadd=2)
    eos = read_compose(tmp_path / "nadd2int")
    out_dir = tmp_path / "internal"
    eos.write_internal_directory(out_dir)
    back = EOSData.from_internal_directory(out_dir)
    assert back.nadd_max == 2
    assert back.q_extra_full.shape == (2, 1, 1, 4)
    assert np.allclose(back.q_extra_full[0, 0, 0, :], eos.q_extra_full[0, 0, 0, :])
    assert np.allclose(back.q_extra_full[1, 0, 0, :], eos.q_extra_full[1, 0, 0, :])


# --- 4. Nadd=0 table still works (backward compatible) ----------------------

def test_nadd0_table_still_reads(tmp_path):
    _write_nadd_table(tmp_path / "nadd0", n=5, nadd=0)
    eos = read_compose(tmp_path / "nadd0")
    assert eos.q_extra_count == 0
    assert eos.nadd_max == 0
    # no q_extra_full tensor allocated for Nadd=0
    assert eos.q_extra_full.size == 0
