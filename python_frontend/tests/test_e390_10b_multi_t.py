"""E390-10B phase 4: multi-point T axis (reader accepts nT>1).

Phase 4 relaxes the nT==1 cold-slice bound. The reader now accepts a table
with a multi-point T axis (nT>1) in addition to the phase-3 multi-Yq support.
Every (T,Yq,density) thermo row is stored in the q_full(nT,nYq,7,n) tensor;
the 1D q(7,n) surface still materializes the selected (first) T,Yq for
backward compatibility with M3-M6 downstream code.

Contract verified here:
  1. A synthetic 3T x 2Yq x n table reads into q_full(3,2,7,n); the 1D Q1
     surface mirrors the first (T,Yq) cell and the other 5 cells carry their
     own distinct Q1 values in the tensor.
  2. select_slice materializes each (iT,iYq) cell into the 1D Q1..Q7 surface
     and records the selected_iT/selected_iYq.
  3. A 3T x 2Yq x n Nadd=1 table round-trips through the internal directory:
     q_full carries all 6 T-Yq blocks and q_extra_full carries the extra slot
     for all 6 blocks.
  4. The cold-slice path (nT==1, T=0) is unaffected (covered by the existing
     suite: read_compose on RG_SLY4, GPPVA_DD2, and the phase 1-3 tests).
"""
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import (  # noqa: E402
    EOSData, read_compose, select_slice, write_compose, executable,
)


def _write_multi_t_table(d: Path, *, n: int, nT: int, nYq: int = 1, nadd: int = 0):
    """Write a minimal CompOSE multidim dir with nT>1 (and optionally nYq>1).

    Each thermo row: iT inb iYq Q1..Q7 [Nadd] [extras]. Q1 carries a distinct
    signature per (iT,iYq,i): 1e-3 + 1e-5*(i + 10*(iy-1) + 100*(it-1)).
    Extras (slot j) carry 100*j + i + 10*(iy-1) + 100*(it-1).
    T axis values are non-zero (multidim, not cold slice): T = 0, 10, 30 MeV.
    """
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


# --- 1. a synthetic 3T x 2Yq x n table reads into q_full ------------------

def test_multi_t_reads_into_q_full(tmp_path):
    _write_multi_t_table(tmp_path / "t3yq2", n=4, nT=3, nYq=2)
    eos = read_compose(tmp_path / "t3yq2")
    assert eos.nT == 3
    assert eos.nYq == 2
    assert eos.q_full.shape == (3, 2, 7, 4)
    # 1D Q1 mirrors the first (T,Yq) cell: 4 points with distinct Q1 per density
    expected_first = np.array([1.0e-3 + 1.0e-5 * i for i in range(4)])
    assert np.allclose(eos.Q1, expected_first)
    # the 6 cells carry distinct Q1 signatures in the tensor
    for it in range(3):
        for iy in range(2):
            for i in range(4):
                expected = 1.0e-3 + 1.0e-5 * (i + 10 * iy + 100 * it)
                assert np.isclose(eos.q_full[it, iy, 0, i], expected), \
                    f"q_full[{it},{iy},0,{i}] mismatch"


# --- 2. select_slice materializes each (iT,iYq) cell ----------------------

def test_multi_t_select_slice_materializes_each_cell(tmp_path):
    _write_multi_t_table(tmp_path / "t3yq2", n=4, nT=3, nYq=2)
    eos = read_compose(tmp_path / "t3yq2")
    for it in range(1, 4):
        for iy in range(1, 3):
            sliced = select_slice(eos, iT=it, iYq=iy)
            assert sliced.selected_iT == it
            assert sliced.selected_iYq == iy
            for i in range(4):
                expected = 1.0e-3 + 1.0e-5 * (i + 10 * (iy - 1) + 100 * (it - 1))
                assert np.isclose(sliced.Q1[i], expected), \
                    f"select_slice(T={it},Yq={iy}) Q1[{i}] mismatch"


# --- 3. 3T x 2Yq x n Nadd=1 round-trips through internal directory -------

def test_multi_t_nadd1_round_trips(tmp_path):
    _write_multi_t_table(tmp_path / "t3yq2n1", n=3, nT=3, nYq=2, nadd=1)
    eos = read_compose(tmp_path / "t3yq2n1")
    assert eos.nT == 3 and eos.nYq == 2
    assert eos.q_full.shape == (3, 2, 7, 3)
    assert eos.q_extra_full.shape == (1, 3, 2, 3)
    # the 1D q_extra mirrors the first (T,Yq) extra slot
    for i in range(3):
        expected = 100 * 1 + i + 0 + 0  # slot=1, i, iy=1, it=1
        assert np.isclose(eos.q_extra[i], expected), f"q_extra[{i}] mismatch"
    # all 6 cells carry distinct extra signatures
    for it in range(3):
        for iy in range(2):
            for i in range(3):
                expected = 100 * 1 + i + 10 * iy + 100 * it
                assert np.isclose(eos.q_extra_full[0, it, iy, i], expected), \
                    f"q_extra_full[0,{it},{iy},{i}] mismatch"


# --- 4. cold slice regression (nT==1, T=0) is unaffected ------------------

def test_cold_slice_still_reads_after_phase4(tmp_path):
    # write a cold slice table (nT=1, T=0) and confirm it still reads
    _write_multi_t_table(tmp_path / "cold", n=3, nT=1, nYq=1)
    eos = read_compose(tmp_path / "cold")
    assert eos.nT == 1 and eos.nYq == 1
    assert eos.q_full.size == 0  # cold slice does not materialize the tensor
    # 1D Q1 carries the cold slice values
    for i in range(3):
        expected = 1.0e-3 + 1.0e-5 * i
        assert np.isclose(eos.Q1[i], expected)
