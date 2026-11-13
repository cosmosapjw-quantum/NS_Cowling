"""EOS-10 regression: M10 axis metadata must round-trip losslessly through the
internal-directory (write_internal_directory -> from_internal_directory) and
through copy().

The Fortran reader already stores the full T/Yq axes and their declared index
bounds on eos_data_t (compose_eos_io.f90: nT, nYq, t_first, t_last, yq_first,
yq_last, selected_iT, selected_iYq, t_axis, yq_axis), and the Python EOSData
dataclass has the corresponding fields. But the internal-directory round-trip
(write_internal_directory / from_internal_directory) and copy() previously
dropped these fields, so a read -> write -> read round-trip lost the generic
CompOSE axis layout. This is the EOS-10 acceptance criterion "T/nB/Yq axes
and index bounds preserved" + "lossless round-trip".

This test does NOT exercise the native Fortran read (that is covered by
test_m10_generic_axes.py); it exercises the Python internal-directory
serialization path directly, which is the lossless-round-trip gate.
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import EOSData  # noqa: E402


def _make_eos_with_m10_axes(n: int = 4) -> EOSData:
    """Build a minimal EOSData populated with non-trivial M10 axis metadata
    (the kind a native read of a 1-based cold slice would produce)."""
    nb = np.linspace(0.05, 0.05 + 0.001 * (n - 1), n)
    q = 1.0e-3 + 1.0e-5 * np.arange(n)
    eos = EOSData(
        nb=nb, mn=939.0, mp=939.0,
        Q1=q, Q2=np.zeros(n), Q3=np.zeros(n), Q4=np.zeros(n),
        Q5=np.zeros(n), Q6=np.zeros(n), Q7=np.zeros(n),
        q_extra=np.zeros(n), Iphase=np.zeros(n, dtype=int),
        Yn=np.zeros(n), Yp=np.zeros(n), Ye=np.zeros(n), Ymu=np.zeros(n),
        Aav=np.zeros(n), Zav=np.zeros(n), Yh=np.zeros(n),
        Il=1, iT=1, iYq=1, inb_start=1,
        temperature_mev=0.0, charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=False,
        # M10 generic axes: a single cold-T=0 slice declared 1..1, plus a
        # 3-point Yq axis that the reader resolved to the cold Yq=0 point.
        nT=1, nYq=3,
        t_first=1, t_last=1, yq_first=1, yq_last=3,
        selected_iT=1, selected_iYq=2,
        t_axis=np.array([0.0]),
        yq_axis=np.array([0.0, 0.05, 0.1]),
    )
    return eos


def test_m10_axis_metadata_round_trips_through_internal_directory(tmp_path):
    """write_internal_directory -> from_internal_directory must preserve the
    M10 axis scalar fields (nT, nYq, t_first, t_last, yq_first, yq_last,
    selected_iT, selected_iYq) and the axis arrays (t_axis, yq_axis)."""
    eos = _make_eos_with_m10_axes()
    eos.write_internal_directory(tmp_path)
    again = EOSData.from_internal_directory(tmp_path)
    assert again.nT == eos.nT, f"nT lost: {again.nT} != {eos.nT}"
    assert again.nYq == eos.nYq, f"nYq lost: {again.nYq} != {eos.nYq}"
    assert again.t_first == eos.t_first, f"t_first lost: {again.t_first} != {eos.t_first}"
    assert again.t_last == eos.t_last, f"t_last lost: {again.t_last} != {eos.t_last}"
    assert again.yq_first == eos.yq_first, f"yq_first lost"
    assert again.yq_last == eos.yq_last, f"yq_last lost"
    assert again.selected_iT == eos.selected_iT, f"selected_iT lost"
    assert again.selected_iYq == eos.selected_iYq, f"selected_iYq lost"
    np.testing.assert_array_equal(again.t_axis, eos.t_axis, err_msg="t_axis lost")
    np.testing.assert_array_equal(again.yq_axis, eos.yq_axis, err_msg="yq_axis lost")


def test_m10_axis_metadata_round_trips_through_copy():
    """copy() must preserve the M10 axis scalar fields and arrays."""
    eos = _make_eos_with_m10_axes()
    cp = eos.copy()
    assert cp.nT == eos.nT
    assert cp.nYq == eos.nYq
    assert cp.t_first == eos.t_first
    assert cp.t_last == eos.t_last
    assert cp.yq_first == eos.yq_first
    assert cp.yq_last == eos.yq_last
    assert cp.selected_iT == eos.selected_iT
    assert cp.selected_iYq == eos.selected_iYq
    np.testing.assert_array_equal(cp.t_axis, eos.t_axis)
    np.testing.assert_array_equal(cp.yq_axis, eos.yq_axis)


def test_m10_double_round_trip_is_stable(tmp_path):
    """A second round-trip must reproduce the first exactly (idempotent)."""
    eos = _make_eos_with_m10_axes()
    eos.write_internal_directory(tmp_path / "a")
    first = EOSData.from_internal_directory(tmp_path / "a")
    first.write_internal_directory(tmp_path / "b")
    second = EOSData.from_internal_directory(tmp_path / "b")
    for fld in ("nT", "nYq", "t_first", "t_last", "yq_first", "yq_last",
                "selected_iT", "selected_iYq"):
        assert getattr(second, fld) == getattr(first, fld), f"{fld} drifted on 2nd round-trip"
    np.testing.assert_array_equal(second.t_axis, first.t_axis)
    np.testing.assert_array_equal(second.yq_axis, first.yq_axis)
