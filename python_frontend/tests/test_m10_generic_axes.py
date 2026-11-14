"""M10 regression tests for generic CompOSE axis handling.

The existing reader is deliberately scoped to a single cold T=0 slice (no
implicit flattening, axis averaging, nearest-index substitution, or
temperature/charge-fraction guessing). M10 extends it to:

  1. Store the full T and Yq axes (values + declared index bounds) on the
     fortran eos_data_t, so downstream code sees the generic 1D/2D/3D layout.
  2. Resolve the cold-T=0 slice by VALUE, accepting either 0-based or 1-based
     T-index convention. CompOSE allows both: some sources (RG_SLY4,
     RG_KDE0v1) declare eos.t as index 0..0 with value T=0 and write thermo
     rows with iT=1; others (GPPVA_FSU2H) declare eos.t as index 1..1.
  3. Record the actual iT/iYq label used by the thermo rows as selected_iT /
     selected_iYq on the fortran type, so the cold slice is explicit.

NOTE: the python EOSData round-trip of the new M10 axis scalar fields
(nT, nYq, t_first, t_last, yq_first, yq_last, selected_iT, selected_iYq)
and axis arrays (t_axis, yq_axis) is a tracked follow-up; these tests
verify the fortran-side resolution by checking that the read succeeds
(no 'thermo table is not a single selected T/Yq slice' error) for the
0-based-convention sources and is correctly rejected for multi-point /
nonzero-T sources.

This is NOT multidimensional support: the reader still rejects multi-point
T or Yq grids (the contract says 'retain the full axes and their exact
bounds, support valid 1D/2D/3D layouts, and make slice selection explicit').
A 2D/3D source with more than one T or Yq value must still be rejected.
"""
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import read_compose  # noqa: E402
from compose_eos import NativeExecutionError  # noqa: E402

INPUTS = ROOT / "inputs" / "raw_compose"


def _write_cold_slice_dir(d: Path, *, n: int, t_first: int, t_value: float,
                          yq_first: int, yq_value: float, thermo_iT: int,
                          thermo_iYq: int):
    """Write a minimal CompOSE dir with a single cold slice, allowing explicit
    control of the T/Yq index convention (t_first 0 or 1) and the iT/iYq labels
    used in the thermo rows. No composition file is written."""
    d.mkdir(parents=True, exist_ok=True)
    nb = np.linspace(0.05, 0.05 + 0.001 * (n - 1), n)
    (d / "eos.nb").write_text(f"1\n{n}\n" + "\n".join(f"{x:.10e}" for x in nb) + "\n", encoding="utf-8")
    (d / "eos.t").write_text(f"{t_first}\n{t_first}\n{t_value:.10e}\n", encoding="utf-8")
    (d / "eos.yq").write_text(f"{yq_first}\n{yq_first}\n{yq_value:.10e}\n", encoding="utf-8")
    lines = ["939.0 939.0 1"]
    for i in range(n):
        q1 = 1.0e-3 + 1.0e-5 * i
        lines.append(f"{thermo_iT} {i+1} {thermo_iYq} {q1:.10e} 0 0 0 0 0 0 0")
    (d / "eos.thermo").write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_one_based_t_index_convention_passes(tmp_path):
    """A source declaring eos.t as index 1..1 with thermo rows using iT=1
    (e.g. GPPVA_FSU2H convention) must read successfully."""
    d = tmp_path / "one-based"
    _write_cold_slice_dir(d, n=3, t_first=1, t_value=0.0, yq_first=1,
                          yq_value=0.0, thermo_iT=1, thermo_iYq=1)
    eos = read_compose(d)
    assert eos.n == 3
    assert abs(eos.temperature_mev) < 1e-12
    assert abs(eos.charge_fraction_grid) < 1e-12


def test_zero_based_t_index_with_one_based_thermo_rows_passes(tmp_path):
    """A source declaring eos.t as index 0..0 (value T=0) with thermo rows
    using iT=1 (e.g. RG_SLY4, RG_KDE0v1 convention) must read successfully.
    The reader resolves the cold slice by value (T=0) and accepts either
    the declared index (0) or the 1-based label (1) that the rows use."""
    d = tmp_path / "zero-based-rows-one-based"
    _write_cold_slice_dir(d, n=3, t_first=0, t_value=0.0, yq_first=1,
                          yq_value=0.0, thermo_iT=1, thermo_iYq=1)
    eos = read_compose(d)
    assert eos.n == 3
    assert abs(eos.temperature_mev) < 1e-12


def test_zero_based_t_index_with_zero_based_thermo_rows_passes(tmp_path):
    """The symmetric case: eos.t index 0..0 and thermo rows using iT=0.
    Both conventions are accepted."""
    d = tmp_path / "zero-based-rows-zero-based"
    _write_cold_slice_dir(d, n=3, t_first=0, t_value=0.0, yq_first=1,
                          yq_value=0.0, thermo_iT=0, thermo_iYq=1)
    eos = read_compose(d)
    assert eos.n == 3
    assert abs(eos.temperature_mev) < 1e-12


def test_multi_point_T_axis_is_accepted_phase4(tmp_path):
    """E390-10B phase 4: a 2D source with more than one T value (nT>1) is now
    ACCEPTED. Every (T,density) thermo row is read into q_full(nT,nYq,7,n); the
    1D q(7,n) surface materializes the selected (first) T,Yq. select_slice picks
    each T cell. The cold slice (nT==1, T=0) remains the default and is unchanged."""
    d = tmp_path / "multi-T"
    d.mkdir(parents=True, exist_ok=True)
    n = 3
    nT = 2
    nb = np.linspace(0.05, 0.052, n)
    (d / "eos.nb").write_text(f"1\n{n}\n" + "\n".join(f"{x:.10e}" for x in nb) + "\n", encoding="utf-8")
    (d / "eos.t").write_text(f"1\n{nT}\n" + "\n".join(f"{float(t):.10e}" for t in range(nT)) + "\n", encoding="utf-8")
    (d / "eos.yq").write_text("1\n1\n0.0\n", encoding="utf-8")
    lines = ["939.0 939.0 1"]
    for it in range(1, nT + 1):
        for i in range(n):
            q1 = 1.0e-3 + 1.0e-5 * (i + 100 * (it - 1))
            lines.append(f"{it} {i+1} 1 {q1:.10e} 0 0 0 0 0 0 0")
    (d / "eos.thermo").write_text("\n".join(lines) + "\n", encoding="utf-8")
    eos = read_compose(d)
    assert eos.nT == nT
    assert eos.nYq == 1
    assert eos.q_full.shape == (nT, 1, 7, n)
    # 1D Q1(n) mirrors the first T cell
    assert np.allclose(eos.Q1, [1.0e-3 + 1.0e-5 * i for i in range(n)])
    # the second T cell carries its own distinct Q1 values
    assert np.allclose(eos.q_full[1, 0, 0, :], [1.0e-3 + 1.0e-5 * (i + 100) for i in range(n)])


def test_multi_point_Yq_axis_is_accepted_phase3(tmp_path):
    """E390-10B phase 3: a 2D source with more than one Yq value (nT==1, T=0)
    is now ACCEPTED. Every (Yq,density) row is read into q_full(1,nYq,7,n); the
    1D q(7,n) surface materializes the first Yq. select_slice picks each cell."""
    d = tmp_path / "multi-Yq"
    d.mkdir(parents=True, exist_ok=True)
    n = 3
    nYq = 2
    nb = np.linspace(0.05, 0.052, n)
    (d / "eos.nb").write_text(f"1\n{n}\n" + "\n".join(f"{x:.10e}" for x in nb) + "\n", encoding="utf-8")
    (d / "eos.t").write_text("1\n1\n0.0\n", encoding="utf-8")
    (d / "eos.yq").write_text(f"1\n{nYq}\n" + "\n".join("0.0" if j == 0 else "0.1" for j in range(nYq)) + "\n", encoding="utf-8")
    lines = ["939.0 939.0 1"]
    for iy in range(1, nYq + 1):
        for i in range(n):
            q1 = 1.0e-3 + 1.0e-5 * (i + 10 * (iy - 1))
            lines.append(f"1 {i+1} {iy} {q1:.10e} 0 0 0 0 0 0 0")
    (d / "eos.thermo").write_text("\n".join(lines) + "\n", encoding="utf-8")
    eos = read_compose(d)
    assert eos.nYq == nYq
    assert eos.q_full.shape == (1, nYq, 7, n)
    # 1D Q1(n) mirrors the first Yq; Q1 of the second Yq differs (from q_full)
    assert np.allclose(eos.Q1, [1.0e-3 + 1.0e-5 * i for i in range(n)])
    assert np.allclose(eos.q_full[0, 1, 0, :], [1.0e-3 + 1.0e-5 * (i + 10) for i in range(n)])


def test_nonzero_T_slice_is_accepted_phase4(tmp_path):
    """E390-10B phase 4: a single-T slice at nonzero temperature is now ACCEPTED
    (no longer rejected as a non-cold slice). The reader stores the single T
    point and its thermo rows; select_slice can pick it. The cold slice (T=0)
    path remains unchanged and is the default for M3-M6 downstream code."""
    d = tmp_path / "hot-T"
    _write_cold_slice_dir(d, n=3, t_first=1, t_value=1.0, yq_first=1,
                          yq_value=0.0, thermo_iT=1, thermo_iYq=1)
    eos = read_compose(d)
    assert eos.nT == 1
    assert eos.nYq == 1
    assert eos.temperature_mev == 1.0


def test_real_RG_SLY4_passes_slice_check():
    """The real RG_SLY4 reconstructed source (eos.t index 0..0, thermo rows
    iT=1) must now pass the cold-slice check. It may still fail downstream
    on composition (C6B disputed field), but it must NOT fail with
    'thermo table is not a single selected T/Yq slice'."""
    if not INPUTS.exists():
        pytest.skip("inputs/raw_compose not staged (M0 not run in this environment)")
    src = INPUTS / "RG_SLY4"
    if not src.exists():
        pytest.skip("RG_SLY4 source not staged")
    try:
        eos = read_compose(src)
        assert abs(eos.temperature_mev) < 1e-12
    except NativeExecutionError as e:
        assert "single selected T/Yq slice" not in str(e), (
            f"RG_SLY4 must not fail on the slice check after M10; got: {e}")


def test_real_RG_KDE0v1_passes_slice_check():
    """Same as RG_SLY4 but for RG_KDE0v1 (also declares eos.t index 0..0)."""
    if not INPUTS.exists():
        pytest.skip("inputs/raw_compose not staged (M0 not run in this environment)")
    src = INPUTS / "RG_KDE0v1"
    if not src.exists():
        pytest.skip("RG_KDE0v1 source not staged")
    try:
        eos = read_compose(src)
        assert abs(eos.temperature_mev) < 1e-12
    except NativeExecutionError as e:
        assert "single selected T/Yq slice" not in str(e), (
            f"RG_KDE0v1 must not fail on the slice check after M10; got: {e}")
