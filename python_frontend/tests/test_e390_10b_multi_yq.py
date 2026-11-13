"""E390-10B phase 3: multi-point Yq axis (reader accepts nYq>1 with nT==1, T=0).

Phase 3 relaxes the nYq==1 cold-slice bound. The reader now accepts a table
with a multi-point Yq axis (nYq>1) as long as nT==1 and T=0. Every (Yq,density)
thermo row is stored in the q_full(nT,nYq,7,n) tensor; the 1D q(7,n) surface
still materializes the selected (first) Yq for backward compatibility with
M3-M6 downstream code.

Contract verified here:
  1. A synthetic 1T x 3Yq x n table reads into q_full(1,3,7,n); the 1D Q1
     surface mirrors the first Yq and the other two Yq blocks carry their own
     distinct Q1 values in the tensor.
  2. select_slice materializes each (iT=1, iYq) cell into the 1D Q1..Q7 surface
     and records the selected_iYq.
  3. A 1T x 2Yq x n Nadd=2 table round-trips through the internal directory:
     q_full carries both Yq blocks and q_extra_full carries both extra slots
     for both Yq blocks.
  4. The internal-directory round-trip preserves nYq and the per-Yq Q1 values.
"""
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import EOSData, read_compose, select_slice  # noqa: E402

INPUTS = ROOT / "inputs" / "raw_compose"


def _write_multi_yq_table(d: Path, *, n: int, nYq: int, nadd: int = 0):
    """Write a minimal CompOSE cold-slice dir with nT==1, T=0 and nYq>1.

    Each thermo row: iT inb iYq Q1..Q7 [Nadd] [extras]. Q1 carries a distinct
    signature per (iy, i): 1e-3 + 1e-5*(i + 10*(iy-1)). Extras (slot j) carry
    100*j + i + 10*(iy-1).
    """
    d.mkdir(parents=True, exist_ok=True)
    nb = np.linspace(0.05, 0.05 + 0.001 * (n - 1), n)
    (d / "eos.nb").write_text(f"1\n{n}\n" + "\n".join(f"{x:.10e}" for x in nb) + "\n")
    (d / "eos.t").write_text("1\n1\n0.0000000000e+00\n")
    yq_vals = "\n".join(f"{0.1 * j:.10e}" for j in range(nYq))
    (d / "eos.yq").write_text(f"1\n{nYq}\n{yq_vals}\n")
    lines = ["939.0 939.0 1"]
    for iy in range(1, nYq + 1):
        for i in range(n):
            q1 = 1.0e-3 + 1.0e-5 * (i + 10 * (iy - 1))
            row = f"1 {i+1} {iy} {q1:.10e} 0 0 0 0 0 0"
            if nadd == 0:
                row += " 0"
            else:
                row += f" {nadd}"
                for slot in range(1, nadd + 1):
                    row += f" {100*slot + i + 10*(iy-1):.10e}"
            lines.append(row)
    (d / "eos.thermo").write_text("\n".join(lines) + "\n")


# --- 1. a synthetic 1T x 3Yq x n table reads into q_full --------------------

def test_multi_yq_reads_into_q_full(tmp_path):
    if not (INPUTS / "RG_SLY4").exists():
        pytest.skip("inputs/raw_compose not staged (E390-00)")
    _write_multi_yq_table(tmp_path / "yq3", n=4, nYq=3)
    eos = read_compose(tmp_path / "yq3")
    assert eos.nYq == 3
    assert eos.nT == 1
    assert eos.q_full.shape == (1, 3, 7, 4)
    # 1D Q1 mirrors the first Yq
    assert np.allclose(eos.Q1, [1.0e-3 + 1.0e-5 * i for i in range(4)])
    # each Yq block carries its own Q1 in the tensor
    for iy in range(3):
        expected = np.array([1.0e-3 + 1.0e-5 * (i + 10 * iy) for i in range(4)])
        assert np.allclose(eos.q_full[0, iy, 0, :], expected)


# --- 2. select_slice materializes each Yq cell -----------------------------

def test_select_slice_materializes_each_yq(tmp_path):
    if not (INPUTS / "RG_SLY4").exists():
        pytest.skip("inputs/raw_compose not staged (E390-00)")
    _write_multi_yq_table(tmp_path / "yq3s", n=5, nYq=3)
    eos = read_compose(tmp_path / "yq3s")
    for iyq in (1, 2, 3):
        cell = select_slice(eos, iT=1, iYq=iyq)
        expected = np.array([1.0e-3 + 1.0e-5 * (i + 10 * (iyq - 1)) for i in range(5)])
        assert np.allclose(cell.Q1, expected), f"Yq={iyq} Q1 mismatch"
        assert cell.iYq == iyq


# --- 3. 1T x 2Yq x n Nadd=2 round-trips through the internal directory ----

def test_multi_yq_nadd2_internal_directory_round_trip(tmp_path):
    if not (INPUTS / "RG_SLY4").exists():
        pytest.skip("inputs/raw_compose not staged (E390-00)")
    _write_multi_yq_table(tmp_path / "yq2n2", n=4, nYq=2, nadd=2)
    eos = read_compose(tmp_path / "yq2n2")
    assert eos.nYq == 2
    assert eos.q_full.shape == (1, 2, 7, 4)
    assert eos.nadd_max == 2
    assert eos.q_extra_full.shape == (2, 1, 2, 4)
    out_dir = tmp_path / "internal"
    eos.write_internal_directory(out_dir)
    back = EOSData.from_internal_directory(out_dir)
    assert back.nYq == 2
    assert back.q_full.shape == (1, 2, 7, 4)
    assert back.nadd_max == 2
    assert back.q_extra_full.shape == (2, 1, 2, 4)
    # every (Yq, extra slot) round-trips
    for iy in range(2):
        exp_q1 = np.array([1.0e-3 + 1.0e-5 * (i + 10 * iy) for i in range(4)])
        assert np.allclose(back.q_full[0, iy, 0, :], exp_q1)
        for slot in (1, 2):
            exp_v = np.array([100 * slot + i + 10 * iy for i in range(4)], dtype=float)
            assert np.allclose(back.q_extra_full[slot - 1, 0, iy, :], exp_v)
