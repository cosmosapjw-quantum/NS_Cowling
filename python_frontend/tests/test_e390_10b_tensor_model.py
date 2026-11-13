"""E390-10B phase 1: full multidim tensor data model extension (no reader change).

Phase 1 adds the q_full / q_extra_full / nadd_max fields to EOSData and a
select_slice bridge that materializes a (iT, iYq) cell into the 1D Q/q_extra
surface, WITHOUT changing the cold-slice reader. The contract:

  1. Cold-slice readers leave q_full / q_extra_full empty by default; the
     existing M3-M6 production path is unaffected.
  2. q_full is opt-in: only a multidim reader (phase 2+) populates it.
  3. select_slice copies the chosen cell into Q1..Q7 / q_extra and records
     selected_iT / selected_iYq.
  4. nadd_max round-trips through the internal namelist.
  5. copy() propagates the tensor fields (deep copy, no shared mutation).
  6. select_slice on an empty q_full raises a clear ValueError (cold-slice
     EOSData has no tensor to slice).
"""
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import EOSData, select_slice, read_compose  # noqa: E402

INPUTS = ROOT / "inputs" / "raw_compose"


def _cold_slice_eos(n: int = 8) -> EOSData:
    """A minimal directly-constructed cold-slice EOSData (no q_full tensor)."""
    nb = np.linspace(0.05, 0.12, n)
    q = np.zeros((7, n))
    q[0] = 1.0e-3  # Q1 (pressure-ish, positive for validate)
    return EOSData(
        nb=nb, mn=939.0, mp=939.0,
        Q1=q[0], Q2=q[1], Q3=q[2], Q4=q[3], Q5=q[4], Q6=q[5], Q7=q[6],
        q_extra=np.zeros(n), Iphase=np.zeros(n, dtype=np.int64),
        Yn=np.zeros(n), Yp=np.zeros(n), Ye=np.zeros(n), Ymu=np.zeros(n),
        Aav=np.zeros(n), Zav=np.zeros(n), Yh=np.zeros(n),
        composition_available=False, composition_file_present=False,
    )


def _multidim_eos(nT: int = 2, nYq: int = 3, n: int = 5, nadd_max: int = 1) -> EOSData:
    """A synthetic multidim EOSData with a populated q_full / q_extra_full tensor.

    q_full has shape (nT, nYq, 7, n); q_extra_full has shape (nadd_max, nT, nYq, n),
    matching the Fortran eos_data_t storage order so the select_slice copy is a
    contiguous slice. Each cell (iT, iYq) is filled with a distinct signature so
    select_slice can verify the correct cell was copied:
    Q1[iT,iYq,k] = 100*iT + 10*iYq + k (iT, iYq are 1-based cell labels).
    """
    base = _cold_slice_eos(n)
    q_full = np.zeros((nT, nYq, 7, n))
    for iT in range(nT):
        for iYq in range(nYq):
            for k in range(n):
                q_full[iT, iYq, 0, k] = 100 * (iT + 1) + 10 * (iYq + 1) + k
                # Q2..Q7 leave 0 for simplicity; Q1 carries the signature
    q_extra_full = np.zeros((nadd_max, nT, nYq, n))
    for iT in range(nT):
        for iYq in range(nYq):
            for k in range(n):
                q_extra_full[0, iT, iYq, k] = 1000 * (iT + 1) + 100 * (iYq + 1) + k
    base.q_full = q_full
    base.q_extra_full = q_extra_full
    base.nT = nT
    base.nYq = nYq
    base.nadd_max = nadd_max
    base.q_extra_count = nadd_max
    # mirror the first extra slot into the 1D q_extra surface so write/read
    # round-trips the slot-1 values (the writer serializes q_extra then reads
    # it back into q_extra_full[0]).
    base.q_extra = q_extra_full[0, 0, 0, :].copy()
    return base


# --- 1. cold-slice path is unaffected: q_full/q_extra_full empty by default ---

def test_cold_slice_eosdata_has_empty_q_full_by_default():
    eos = _cold_slice_eos()
    assert eos.q_full.size == 0
    assert eos.q_extra_full.size == 0
    assert eos.nadd_max == 0


def test_cold_slice_read_compose_leaves_q_full_empty():
    """read_compose on a real cold-slice source must NOT populate q_full.

    Phase 2 update: a cold-slice Nadd=1 source now materializes q_extra_full
    (nadd_max=1, shape (1,1,1,n)) mirroring the 1D q_extra surface — this is
    the Nadd-relaxation contract. q_full stays empty because the cold slice
    has nT=nYq=1 and the multidim tensor is only materialized by a phase-3/4
    multidim reader. nadd_max reflects the source's declared Nadd.
    """
    if not (INPUTS / "RG_SLY4").exists():
        pytest.skip("inputs/raw_compose not staged (E390-00)")
    eos = read_compose(INPUTS / "RG_SLY4")
    assert eos.q_full.size == 0, "cold-slice reader must leave q_full empty (multidim opt-in)"
    # RG_SLY4 is Nadd=1: q_extra_full is populated with nadd_max=1, mirroring q_extra
    assert eos.nadd_max == 1
    assert eos.q_extra_full.ndim == 4
    assert eos.q_extra_full.shape == (1, 1, 1, eos.n)
    assert np.allclose(eos.q_extra_full[0, 0, 0, :], eos.q_extra)


# --- 2. select_slice materializes the chosen cell into the 1D surface -------

def test_select_slice_copies_chosen_cell_into_Q_surface():
    eos = _multidim_eos(nT=2, nYq=3, n=5, nadd_max=1)
    out = select_slice(eos, iT=2, iYq=3)
    # Q1 signature for cell (2,3): 100*2 + 10*3 + k = 230 + k
    expected_q1 = np.array([230, 231, 232, 233, 234], dtype=float)
    assert np.allclose(out.Q1, expected_q1)
    # q_extra signature: 1000*2 + 100*3 + k = 2300 + k
    expected_qe = np.array([2300, 2301, 2302, 2303, 2304], dtype=float)
    assert np.allclose(out.q_extra, expected_qe)
    assert out.selected_iT == 2
    assert out.selected_iYq == 3


def test_select_slice_records_selected_indices_and_preserves_tensor():
    eos = _multidim_eos(nT=2, nYq=3, n=5, nadd_max=1)
    out = select_slice(eos, iT=1, iYq=2)
    assert out.selected_iT == 1
    assert out.selected_iYq == 2
    # the tensor is preserved on the view so another cell can be selected
    assert out.q_full.shape == (2, 3, 7, 5)
    assert out.q_extra_full.shape == (1, 2, 3, 5)


def test_select_slice_does_not_mutate_source_or_other_cells():
    eos = _multidim_eos(nT=2, nYq=2, n=3, nadd_max=1)
    out_a = select_slice(eos, iT=1, iYq=1)
    out_b = select_slice(eos, iT=2, iYq=2)
    # source tensor is unchanged. NumPy is 0-based: q_full[1,1] is cell (iT=2,iYq=2)
    # with signature 100*2 + 10*2 + k = 220 + k. q_full axis order is (nT,nYq,7,n).
    assert np.allclose(eos.q_full[1, 1, 0, :], [220, 221, 222])
    # the two slices differ (deep copy, not a shared reference)
    assert not np.allclose(out_a.Q1, out_b.Q1)
    assert np.allclose(out_a.Q1, [110, 111, 112])
    assert np.allclose(out_b.Q1, [220, 221, 222])


# --- 3. select_slice on empty q_full raises a clear error --------------------

def test_select_slice_raises_on_empty_q_full():
    eos = _cold_slice_eos()
    with pytest.raises(ValueError, match="q_full is not materialized"):
        select_slice(eos, iT=1, iYq=1)


def test_select_slice_raises_on_out_of_range_indices():
    eos = _multidim_eos(nT=2, nYq=3, n=5, nadd_max=1)
    with pytest.raises(ValueError, match="iT/iYq out of range"):
        select_slice(eos, iT=3, iYq=1)  # iT=3 > nT=2
    with pytest.raises(ValueError, match="iT/iYq out of range"):
        select_slice(eos, iT=1, iYq=0)  # iYq=0 < 1


def test_select_slice_raises_on_inconsistent_tensor_shape():
    eos = _multidim_eos(nT=2, nYq=3, n=5, nadd_max=1)
    # corrupt the tensor shape so it disagrees with data.nT/nYq
    eos.q_full = np.zeros((4, 4, 7, 5))  # wrong axis lengths
    with pytest.raises(ValueError, match="disagree with data.nT"):
        select_slice(eos, iT=1, iYq=1)


# --- 4. nadd_max round-trips through the internal namelist -------------------

def test_nadd_max_round_trips_through_internal_directory(tmp_path):
    eos = _multidim_eos(nT=2, nYq=2, n=4, nadd_max=1)
    # write the internal directory (1D surface + q_extra_N columns for slots 2+;
    # the nadd_max=1 case serializes only the q_extra slot which mirrors q_extra_full[0])
    eos.write_internal_directory(tmp_path)
    back = EOSData.from_internal_directory(tmp_path)
    assert back.nadd_max == 1
    assert back.q_full.size == 0  # q_full is in-memory only (multidim reader opt-in)
    # phase 2: q_extra_full IS serialized for nadd_max>=1; slot 1 mirrors q_extra
    assert back.q_extra_full.ndim == 4
    assert back.q_extra_full.shape[0] == 1
    assert np.allclose(back.q_extra_full[0, 0, 0, :], eos.q_extra_full[0, 0, 0, :])


def test_nadd_max_defaults_to_zero_for_legacy_namelist(tmp_path):
    """A legacy namelist without nadd_max must read back as 0 (backward compat)."""
    eos = _cold_slice_eos(n=4)
    eos.write_internal_directory(tmp_path)
    # strip the nadd_max line to simulate a legacy namelist
    nml = (tmp_path / "eos_data.nml").read_text()
    nml = nml.replace(" nadd_max=0\n", "")
    (tmp_path / "eos_data.nml").write_text(nml)
    back = EOSData.from_internal_directory(tmp_path)
    assert back.nadd_max == 0


# --- 5. copy() propagates the tensor fields (deep copy) ----------------------

def test_copy_propagates_q_full_and_q_extra_full():
    eos = _multidim_eos(nT=2, nYq=2, n=3, nadd_max=1)
    cp = eos.copy()
    assert cp.q_full.shape == (2, 2, 7, 3)
    assert cp.q_extra_full.shape == (1, 2, 2, 3)
    assert cp.nadd_max == 1
    assert cp.nT == 2 and cp.nYq == 2
    # deep copy: mutating the copy must not affect the source
    cp.q_full[0, 0, 0, 0] = -999.0
    assert eos.q_full[0, 0, 0, 0] != -999.0


def test_copy_of_cold_slice_keeps_empty_tensor():
    eos = _cold_slice_eos()
    cp = eos.copy()
    assert cp.q_full.size == 0
    assert cp.q_extra_full.size == 0
    assert cp.nadd_max == 0
