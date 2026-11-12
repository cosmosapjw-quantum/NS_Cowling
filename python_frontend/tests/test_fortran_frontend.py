import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import (read_compose, write_compose, check_consistency,
                         NativeFeatureUnavailable, NativeExecutionError)
SAMPLE = ROOT / "eos_generator" / "eos_data_normalized" / "GPPVA_FSU2H"


def _append_duplicate_first_record(path: Path, *, skip_header: bool = False) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    records = [line for line in lines if line.strip() and not line.lstrip().startswith("#")]
    if not records:
        raise AssertionError(f"fixture has no records: {path}")
    first = records[1] if skip_header else records[0]
    lines.append(first)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_fortran_reader_verifier_and_roundtrip(tmp_path):
    eos = read_compose(SAMPLE)
    assert eos.n == 300
    assert eos.composition_available
    assert np.all(np.diff(eos.nb) > 0)
    assert np.array_equal(eos.p, eos.nb * eos.Q1)
    report = check_consistency(eos)
    assert report["overall"]["pass_"]
    out = tmp_path / "written"
    write_compose(out, eos)
    again = read_compose(out)
    for name in ("nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7",
                 "Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh"):
        assert np.array_equal(getattr(eos, name), getattr(again, name)), name


def test_alternate_eos_roundtrip(tmp_path):
    source = ROOT / "eos_generator" / "eos_data_normalized" / "BBB_BHF2"
    eos = read_compose(source)
    path = tmp_path / "roundtrip"
    write_compose(path, eos)
    again = read_compose(path)
    assert np.array_equal(eos.Q1, again.Q1)
    assert np.array_equal(eos.nb, again.nb)
    assert again.composition_available == eos.composition_available
    for name in ("Q2","Q3","Q4","Q5","Q6","Q7"):
        assert np.array_equal(getattr(eos,name),getattr(again,name))


def test_native_density_gap_is_available_and_remaining_generator_operations_fail_explicitly():
    from compose_eos import fill_density_gap, reconstruct_master, match_crust_core, NativeExecutionError
    eos = read_compose(SAMPLE)
    filled = fill_density_gap(eos, n_fill=3, gap_phase=7)
    assert filled.n == eos.n + 3
    # Reconstruction is now implemented natively. This particular production
    # table is deliberately refused by the v1 interpolant's physical gate;
    # it must not fall back to the legacy Python smoother or emit an unchecked EOS.
    with pytest.raises(NativeExecutionError, match="master reconstruction .*failure"):
        reconstruct_master(eos)
    with pytest.raises(NativeFeatureUnavailable):
        match_crust_core(eos)


def test_normalize_compose_and_verify_csv(tmp_path):
    from compose_eos import normalize_compose
    norm = tmp_path / "normalized.csv"
    res = normalize_compose(SAMPLE, norm)
    assert res.is_file()
    assert res.stat().st_size > 0


def test_differential_corpus_18_eos():
    import json
    diff_file = ROOT / "fortran_eos" / "differential_corpus_v4" / "fortran_compose_differential.json"
    assert diff_file.is_file()
    data = json.loads(diff_file.read_text())
    assert data["n_tables"] == 18
    assert data["n_tested"] == 18
    assert data["n_failed"] == 0
    assert data["roundtrip_exact"] == 18
    assert len(data["tables"]) == 18


def test_tov_solver_star_integration():
    from compose_cowling.solver import solve_tov_star
    source = ROOT / "eos_generator" / "eos_data_normalized" / "BBB_BHF2"
    eos = read_compose(source)
    star = solve_tov_star(eos, pc_pa=1.5e34)
    assert star["status"] == "ok"
    assert star["backend"] == "Fortran"
    assert star["mass_nominal_solar"] > 1.0
    assert star["radius_m"] > 9000.0


def test_combined_spectrum_234_modes_integrity():
    import csv
    spec_file = ROOT / "combined_spectrum_records.csv"
    assert spec_file.is_file()
    with spec_file.open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 234
    assert len(set(r["name"] for r in rows)) == 234
    for r in rows:
        assert r["mode_status"] == "VALID"
        f_w = int(r["f_inner_W_nodes"]) + int(r["f_outer_W_nodes"])
        f_v = int(r["f_inner_V_nodes"]) + int(r["f_outer_V_nodes"])
        assert f_w == 0 and f_v == 0
        p1_w = int(r["p1_inner_W_nodes"]) + int(r["p1_outer_W_nodes"])
        p1_v = int(r["p1_inner_V_nodes"]) + int(r["p1_outer_V_nodes"])
        assert p1_w == 1 and p1_v == 1



def test_internal_metadata_parses_fortran_l1_logical_values(tmp_path):
    from compose_eos.eos_data import _read_metadata

    metadata = tmp_path / "eos_data.nml"
    metadata.write_text(
        "&metadata\n"
        " composition_file_present=F\n"
        " composition_available=F\n"
        " beta_equilibrated=T\n"
        "/\n",
        encoding="utf-8",
    )
    parsed = _read_metadata(metadata)
    assert parsed["composition_file_present"] is False
    assert parsed["composition_available"] is False
    assert parsed["beta_equilibrated"] is True


def test_duplicate_thermo_density_key_is_rejected(tmp_path):
    import shutil

    source = tmp_path / "duplicate-thermo"
    shutil.copytree(SAMPLE, source)
    _append_duplicate_first_record(source / "eos.thermo", skip_header=True)
    with pytest.raises(NativeExecutionError, match="duplicate eos.thermo row key"):
        read_compose(source)


def test_duplicate_composition_density_key_is_rejected(tmp_path):
    import shutil

    source = tmp_path / "duplicate-compo"
    shutil.copytree(SAMPLE, source)
    _append_duplicate_first_record(source / "eos.compo")
    with pytest.raises(NativeExecutionError, match="duplicate eos.compo row key"):
        read_compose(source)


def test_row_level_composition_presence_roundtrip_preserves_partial_and_explicit_zero(tmp_path):
    eos = read_compose(SAMPLE)
    assert eos.composition_file_present
    assert eos.composition_available

    # Row 0 is an explicit phase-only row. Preserve Iphase but do not invent
    # particle/group abundances from numerical zero placeholders.
    for name in ("Yn_present", "Yp_present", "Ye_present", "Ymu_present",
                 "Abar_present", "Zbar_present", "Yh_present"):
        getattr(eos, name)[0] = False
    for name in ("Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh"):
        getattr(eos, name)[0] = 0.0
    eos.C3_applicable[0] = False
    eos.C4_applicable[0] = False
    eos.row_status[0] = 1

    # A zero muon fraction is nevertheless explicitly present in the selected row.
    zero_muon = int(np.argmin(np.abs(eos.Ymu)))
    eos.Ymu[zero_muon] = 0.0
    eos.Ymu_present[zero_muon] = True

    out = tmp_path / "partial-roundtrip"
    write_compose(out, eos)
    again = read_compose(out)

    assert again.Iphase[0] == eos.Iphase[0]
    assert not again.Yn_present[0] and not again.Yp_present[0]
    assert not again.Ye_present[0] and not again.Ymu_present[0]
    assert not again.Abar_present[0] and not again.Zbar_present[0] and not again.Yh_present[0]
    assert again.row_status[0] == 1
    assert again.Ymu_present[zero_muon]
    assert again.Ymu[zero_muon] == 0.0


def test_python_frontend_refuses_missing_row_presence_masks(tmp_path):
    eos = read_compose(SAMPLE)
    eos.Yn_present = np.zeros(0, dtype=bool)
    with pytest.raises(ValueError, match="missing row-level presence mask Yn_present"):
        write_compose(tmp_path / "lossy-export", eos)
