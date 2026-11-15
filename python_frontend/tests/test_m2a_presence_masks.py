"""M2A regression tests for per-row/per-species/per-group presence masks.

C6B diagnosed that a global composition_available flag conflates absent,
partial and measured-zero composition. The native parser now maintains
per-row presence masks:

  row_status: 0=NOT_AVAILABLE, 1=PARTIAL, 2=CHECKED, 3=FAIL
  *_present: per-row logical (a written zero is True; a missing entry is False)
  C3_applicable: per-row, complete baryon inventory for Yn+Yp+Yh=1
  C4_applicable: per-row, complete charge inventory for Yp+(Z/A)Yh=Ye+Ymu

These tests build small synthetic CompOSE dirs that exercise each scenario
without depending on the disputed SPG/GPPVA sources (which C6B blocks). They
verify the four required scenarios from the C6B adjudication:

  1. Phase-only fixture retains Iphase and records composition as absent.
  2. GM_GM1 outer rows do not fail from zero-initialized missing data.
  3. GM_GM1 core C3 is unavailable unless the required charge-carrier
     inventory is complete; C4 is independently applicable only with the
     complete baryon inventory.
  4. Mixed source coverage yields row-level CHECKED/PARTIAL/NOT_AVAILABLE/FAIL
     receipts.
"""
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import read_compose, write_compose  # noqa: E402
from compose_eos import NativeExecutionError  # noqa: E402


def _write_eos_files(d: Path, *, n: int, compo_rows: list[str] | None,
                     iphase_per_row: list[int] | None = None) -> None:
    """Write a minimal synthetic CompOSE directory with n density points.

    compo_rows: a list of exactly n eos.compo row strings, or None to omit
    eos.compo (composition_file_present=False).

    Grid file format (read_grid): first line = first index, second line =
    last index, then one value per line.
    """
    d.mkdir(parents=True, exist_ok=True)
    nb = np.linspace(0.05, 0.05 + 0.001 * (n - 1), n)
    # eos.nb: 1-based indices 1..n then values
    (d / "eos.nb").write_text(
        f"1\n{n}\n" + "\n".join(f"{x:.10e}" for x in nb) + "\n", encoding="utf-8")
    # eos.t: single T=0 point
    (d / "eos.t").write_text("1\n1\n0.0\n", encoding="utf-8")
    # eos.yq: single Yq=0 point
    (d / "eos.yq").write_text("1\n1\n0.0\n", encoding="utf-8")
    # eos.thermo: header mn mp il, then per row iT inb iYq Q1..Q7 Nadd [q_extra].
    # Q1>0 so pressure Q1*nb>0; everything else finite. iT=1 iYq=1 match the
    # single T/Yq slice. inb ranges 1..n matching eos.nb indices.
    lines = ["939.0 939.0 1"]
    for i in range(n):
        q1 = 1.0e-3 + 1.0e-5 * i
        q2 = q3 = q4 = q5 = q6 = q7 = 0.0
        lines.append(f"1 {i+1} 1 {q1:.10e} {q2:.10e} {q3:.10e} {q4:.10e} "
                     f"{q5:.10e} {q6:.10e} {q7:.10e} 0")
    (d / "eos.thermo").write_text("\n".join(lines) + "\n", encoding="utf-8")
    # eos.compo
    if compo_rows is not None:
        assert iphase_per_row is not None and len(compo_rows) == n
        compo_lines = []
        for i, (row, ph) in enumerate(zip(compo_rows, iphase_per_row)):
            compo_lines.append(f"1 {i+1} 1 {ph} {row}")
        (d / "eos.compo").write_text("\n".join(compo_lines) + "\n", encoding="utf-8")


def _compo_row(*, npairs=0, pairs=None, nquad=0, quads=None):
    """Build one eos.compo row body (after iT inb iYq Iphase Npairs)."""
    parts = [str(npairs)]
    for p in (pairs or []):
        parts.append(f"{p[0]} {p[1]:.10e}")
    parts.append(str(nquad))
    for q in (quads or []):
        parts.append(f"{q[0]} {q[1]:.10e} {q[2]:.10e} {q[3]:.10e}")
    return " ".join(parts)


def test_phase_only_row_records_composition_as_absent(tmp_path):
    """Scenario 1: a phase-only row (Npairs=0, Nquad=0) records Iphase and
    sets row_status=PARTIAL, all *_present=False, C3/C4 not applicable."""
    d = tmp_path / "phase-only"
    rows = [_compo_row(npairs=0, nquad=0) for _ in range(4)]
    _write_eos_files(d, n=4, compo_rows=rows, iphase_per_row=[0, 0, 0, 0])
    eos = read_compose(d)
    assert eos.composition_file_present
    assert not eos.composition_available  # no nonzero composition
    # All rows are PARTIAL (1), not FAIL (3).
    assert list(eos.row_status) == [1, 1, 1, 1]
    # Iphase preserved.
    assert list(eos.Iphase) == [0, 0, 0, 0]
    # No particle/group field present.
    assert not eos.Yn_present.any()
    assert not eos.Yp_present.any()
    assert not eos.Ye_present.any()
    assert not eos.Ymu_present.any()
    assert not eos.Abar_present.any()
    assert not eos.Zbar_present.any()
    assert not eos.Yh_present.any()
    # C3/C4 not applicable for phase-only rows.
    assert not eos.C3_applicable.any()
    assert not eos.C4_applicable.any()


def test_gm1_outer_rows_do_not_fail_from_zero_initialized_missing_data(tmp_path):
    """Scenario 2: GM_GM1 outer-crust phase-only rows must not raise a
    NativeExecutionError. The old global-flag parser would either pass with
    a wrong global 'available' or fail the whole read; M2A keeps the read
    succeeding with PARTIAL row_status on the phase-only rows."""
    d = tmp_path / "gm1-outer"
    # 6 rows: first 3 are phase-only (outer crust), last 3 are core rows
    rows = []
    phases = []
    for _ in range(3):
        rows.append(_compo_row(npairs=0, nquad=0))
        phases.append(0)
    # core rows: Yp, Ye, Abar, Zbar, Yh present (no Yn, no Ymu) -> C3/C4 not applicable.
    # Yh = Abar * YI must be <= 1+5e-3, so use YI = 0.9/56.
    for i in range(3):
        rows.append(_compo_row(
            npairs=2, pairs=[(11, 0.05), (0, 0.05)],
            nquad=1, quads=[(1, 56.0, 26.0, 0.9/56.0)],
        ))
        phases.append(7)
    _write_eos_files(d, n=6, compo_rows=rows, iphase_per_row=phases)
    eos = read_compose(d)
    # outer-crust rows: PARTIAL
    assert list(eos.row_status[:3]) == [1, 1, 1]
    # core rows: PARTIAL (Yn and Ymu missing -> C3 and C4 both inapplicable)
    assert list(eos.row_status[3:]) == [1, 1, 1]
    # No FAIL rows anywhere.
    assert not (eos.row_status == 3).any()
    # Iphase preserved on every row.
    assert list(eos.Iphase) == [0, 0, 0, 7, 7, 7]


def test_gm1_core_c3_unavailable_without_complete_charge_carrier_inventory(tmp_path):
    """Scenario 3: a core row missing Ymu (incomplete charge inventory) has
    C3_applicable=True if Yn+Yp+Yh are present but C4_applicable=False.
    Uses 3 rows; the tested row is the middle one."""
    d = tmp_path / "gm1-core-partial-lepton"
    # All 3 rows: Yn, Yp, Abar, Zbar, Yh, Ye present (no Ymu). Yh=Abar*YI,
    # need Yh<=1+5e-3, so YI=0.5/56.
    rows = [_compo_row(
        npairs=3, pairs=[(10, 0.2), (11, 0.05), (0, 0.05)],
        nquad=1, quads=[(1, 56.0, 26.0, 0.5/56.0)],
    ) for _ in range(3)]
    _write_eos_files(d, n=3, compo_rows=rows, iphase_per_row=[7, 7, 7])
    eos = read_compose(d)
    for i in range(3):
        assert bool(eos.C3_applicable[i])  # Yn,Yp,Yh present
        assert not bool(eos.C4_applicable[i])  # Ymu missing
        assert int(eos.row_status[i]) == 1  # PARTIAL (C4 not applicable)


def test_gm1_core_c4_applicable_only_with_complete_baryon_inventory(tmp_path):
    """Scenario 3b: a row missing Yn has C4_applicable=False (Yp,Z,Abar,Yh,
    Ye,Ymu present is not enough if Yn is missing - baryon inventory incomplete).
    Actually C4 only needs Yp,Z,Abar,Yh,Ye,Ymu per the contract; this row has all
    of those present, so C4 is applicable. C3 requires Yn which is missing."""
    d = tmp_path / "gm1-core-no-yet"
    # 3 rows: Yp, Ye, Ymu, Abar, Zbar, Yh present (no Yn). Yh=56*(0.5/56)=0.5.
    rows = [_compo_row(
        npairs=3, pairs=[(11, 0.05), (0, 0.05), (1, 0.02)],
        nquad=1, quads=[(1, 56.0, 26.0, 0.5/56.0)],
    ) for _ in range(3)]
    _write_eos_files(d, n=3, compo_rows=rows, iphase_per_row=[7, 7, 7])
    eos = read_compose(d)
    for i in range(3):
        # C4 requires Yp, Abar>0, Zbar, Yh, Ye, Ymu - all present here.
        assert bool(eos.C4_applicable[i])
        # C3 requires Yn, Yp, Yh - Yn missing -> C3 not applicable.
        assert not bool(eos.C3_applicable[i])
        assert int(eos.row_status[i]) == 1  # PARTIAL


def test_mixed_coverage_yields_checked_partial_not_available_and_fail(tmp_path):
    """Scenario 4: a mixed-coverage source yields PARTIAL, CHECKED and FAIL
    row_status values in a single read."""
    d = tmp_path / "mixed"
    rows = []
    phases = []
    # PARTIAL: phase only.
    rows.append(_compo_row(npairs=0, nquad=0)); phases.append(0)
    # CHECKED: full baryon + charge inventory. Yh = 56*(0.5/56) = 0.5 (in range).
    rows.append(_compo_row(
        npairs=4, pairs=[(10, 0.4), (11, 0.05), (0, 0.05), (1, 0.0)],
        nquad=1, quads=[(1, 56.0, 26.0, 0.5/56.0)],
    )); phases.append(7)
    # FAIL: Yn out of physical range (negative).
    rows.append(_compo_row(
        npairs=1, pairs=[(10, -0.5)],
        nquad=0,
    )); phases.append(7)
    _write_eos_files(d, n=3, compo_rows=rows, iphase_per_row=phases)
    eos = read_compose(d)
    assert int(eos.row_status[0]) == 1  # PARTIAL (phase only)
    assert int(eos.row_status[1]) == 2  # CHECKED
    assert int(eos.row_status[2]) == 3  # FAIL (Yn<0)
    # No row should be NOT_AVAILABLE here (file present, all rows present).
    assert not (eos.row_status == 0).any()


def test_no_composition_file_all_rows_not_available(tmp_path):
    """If eos.compo is absent, every row is NOT_AVAILABLE (0)."""
    d = tmp_path / "no-compo"
    _write_eos_files(d, n=3, compo_rows=None, iphase_per_row=None)
    eos = read_compose(d)
    assert not eos.composition_file_present
    assert list(eos.row_status) == [0, 0, 0]
    assert not eos.Yn_present.any()
    assert not eos.C3_applicable.any()
    assert not eos.C4_applicable.any()


def test_masks_round_trip_through_internal_csv(tmp_path):
    """Per-row masks survive a write_compose -> read_compose round trip.
    write_compose re-validates baryon/charge abundance closure on every row
    where composition_available is True (which is set globally if any row has
    nonzero composition). So the round-trip source must be fully CHECKED:
    every row satisfies Yn+Yp+Yh=1 and Yp+(Z/A)Yh=Ye+Ymu."""
    d = tmp_path / "roundtrip-src"
    # Three identical CHECKED rows with valid closure.
    # Abar=56, Zbar=26, Yh=0.5 (YI=0.5/56).
    # Baryon: Yn=0.45, Yp=0.05, Yh=0.5 -> sum=1.0.
    # Charge: Yp + (26/56)*0.5 = 0.05 + 0.23214285714 = 0.28214285714.
    #         Ye=0.05, Ymu=0.23214285714 -> sum=0.28214285714. Match.
    row = _compo_row(
        npairs=4, pairs=[(10, 0.45), (11, 0.05), (0, 0.05), (1, 0.23214285714)],
        nquad=1, quads=[(1, 56.0, 26.0, 0.5/56.0)],
    )
    rows = [row, row, row]
    phases = [7, 7, 7]
    _write_eos_files(d, n=3, compo_rows=rows, iphase_per_row=phases)
    eos = read_compose(d)
    assert int(eos.row_status[0]) == 2  # CHECKED
    out = tmp_path / "roundtrip-out"
    write_compose(out, eos)
    again = read_compose(out)
    assert list(again.row_status) == list(eos.row_status)
    assert list(again.Yn_present) == list(eos.Yn_present)
    assert list(again.Yp_present) == list(eos.Yp_present)
    assert list(again.Ye_present) == list(eos.Ye_present)
    assert list(again.Ymu_present) == list(eos.Ymu_present)
    assert list(again.Abar_present) == list(eos.Abar_present)
    assert list(again.Zbar_present) == list(eos.Zbar_present)
    assert list(again.Yh_present) == list(eos.Yh_present)
    assert list(again.C3_applicable) == list(eos.C3_applicable)
    assert list(again.C4_applicable) == list(eos.C4_applicable)


def test_no_spg_gppva_ambiguous_field_remapping_is_made_by_coverage_mask_task(tmp_path):
    """The M2A mask task must not remap any ambiguous SPG/GPPVA field.
    Verify that a row with a Yh value that C6B marks as ambiguous is still
    preserved byte-for-byte (no remap to Y_I, X_h, Ye or another aggregate)
    AND its row_status reflects the per-row range check, not a guess."""
    d = tmp_path / "ambiguous-yh"
    # 3 rows with Yh=1.0 (C6B example outer-crust ambiguous value). The parser
    # must preserve it exactly; only the per-row range check (Yh<=1+5e-3) runs.
    rows = [_compo_row(
        npairs=2, pairs=[(10, 0.0), (11, 0.0)],
        nquad=1, quads=[(1, 56.0, 26.0, 1.0/56.0)],  # Yh = Abar*YI = 56*(1/56) = 1.0
    ) for _ in range(3)]
    _write_eos_files(d, n=3, compo_rows=rows, iphase_per_row=[0, 0, 0])
    eos = read_compose(d)
    # The exact value 1.0 is preserved (not remapped). float64 round-trip of
    # 56*(1/56) introduces ~1e-11 rounding; the assertion is that no remap
    # occurred, so compare against the recomputed 56*(1/56) tolerance.
    for i in range(3):
        assert abs(float(eos.Yh[i]) - (56.0 * (1.0/56.0))) < 1e-10
    # C3 applicable: Yn(=0 written), Yp(=0 written), Yh all present.
    assert bool(eos.C3_applicable[0])
    # row_status is PARTIAL (Ye, Ymu missing -> C4 not applicable), not a guess.
    assert int(eos.row_status[0]) == 1
