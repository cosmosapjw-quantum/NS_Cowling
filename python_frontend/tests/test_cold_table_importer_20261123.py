"""Native raw-table transport: analytic thermodynamics and hostile input."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import pytest

ROOT = Path(__file__).resolve().parents[2]
LAYOUTS = {"gmrs_2022": (11, 7, 8), "gogny_2020": (9, 4, 5),
           "bsk_2018": (28, 9, 7), "gogny_full_2020": (8, 1, 2)}


@pytest.fixture(scope="module")
def importer(tmp_path_factory):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran required")
    build = tmp_path_factory.mktemp("native_raw_import")
    binary = build / "importer"
    subprocess.run([
        "gfortran", "-std=f2018", "-O0", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
        "-ffree-line-length-none", "-J" + str(build), "-I" + str(build),
        str(ROOT / "fortran_eos/src/compose_eos_types.f90"),
        str(ROOT / "fortran_eos/src/compose_eos_io.f90"),
        str(ROOT / "fortran_eos/app/compose_eos_import_cold_table.f90"),
        "-o", str(binary),
    ], check=True, capture_output=True, text=True)
    return binary


def source(path, layout, alteration=None):
    columns, energy, pressure = LAYOUTS[layout]
    lines = []
    for i, n in enumerate((.01, .02, .03)):
        row = ["0"] * columns
        row[0] = format(n, ".17g")
        row[energy] = str(i + 1)
        row[pressure] = (".01", ".04", ".09")[i]
        lines.append(row)
    if alteration:
        alteration(lines)
    path.write_text("# analytic cold barotrope: epsilon=m*n+100*n**2\n\n" +
                    "\n".join(" ".join(row) for row in lines) + "\n")


def run_import(binary, layout, raw, out, mass="939", proton="938", bounds=None):
    args = [str(binary), layout, str(raw), str(out), mass, proton]
    if bounds is not None:
        args.extend(bounds)
    return subprocess.run(args, capture_output=True, text=True)


@pytest.mark.parametrize("layout", LAYOUTS)
def test_native_projection_preserves_known_analytic_table(importer, tmp_path, layout):
    raw = tmp_path / "source.dat"
    source(raw, layout)
    digest = hashlib.sha256(raw.read_bytes()).hexdigest()
    out = tmp_path / "projection"
    result = run_import(importer, layout, raw, out)
    assert result.returncode == 0, result.stdout + result.stderr
    with (out / "normalized_eos.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 3
    for i, row in enumerate(rows):
        n = (.01, .02, .03)[i]
        assert float(row["n_b_fm3"]) == n
        assert float(row["pressure_mev_fm3"]) == pytest.approx(100*n*n, rel=2e-15)
        assert float(row["energy_mev_fm3"]) == pytest.approx(939*n+100*n*n, rel=2e-15)
        assert float(row["mu_b_mev"]) == pytest.approx(939+200*n, rel=2e-15)
    report = json.loads((out / "IMPORT.json").read_text())
    assert report["status"] == "NATIVE_RAW_PROJECTION_READBACK_VERIFIED"
    assert report["rows"] == 3
    assert report["source_columns_one_based"] == [1, LAYOUTS[layout][1]+1, LAYOUTS[layout][2]+1]
    assert report["physical_acceptance_asserted"] is False
    assert report["composition_available"] is False
    assert not (out / "eos.compo").exists()
    assert hashlib.sha256(raw.read_bytes()).hexdigest() == digest


@pytest.mark.parametrize("alteration", [
    lambda rows: rows[1].__setitem__(0, "NaN"),
    lambda rows: rows[1].__setitem__(7, "Infinity"),
    lambda rows: rows[1].__setitem__(8, "bad"),
    lambda rows: rows[1].__setitem__(0, ".01"),
    lambda rows: rows[1].__setitem__(0, "0"),
    lambda rows: rows[1].append("9"),
    lambda rows: rows[1].pop(),
    lambda rows: rows.__delitem__(slice(1, None)),
])
def test_malformed_rows_fail_before_output_creation(importer, tmp_path, alteration):
    raw = tmp_path / "source.dat"
    source(raw, "gmrs_2022", alteration)
    out = tmp_path / "rejected"
    result = run_import(importer, "gmrs_2022", raw, out)
    assert result.returncode != 0
    assert "ERROR:" in result.stdout
    assert not out.exists()


@pytest.mark.parametrize("mass", ["NaN", "Infinity", "-939", "0"])
def test_mass_normalization_must_be_explicit_finite_positive(importer, tmp_path, mass):
    raw = tmp_path / "source.dat"
    source(raw, "gmrs_2022")
    out = tmp_path / "rejected"
    result = run_import(importer, "gmrs_2022", raw, out, mass)
    assert result.returncode != 0
    assert not out.exists()


def test_existing_output_is_never_reused_or_overwritten(importer, tmp_path):
    raw = tmp_path / "source.dat"
    source(raw, "gmrs_2022")
    out = tmp_path / "existing"
    out.mkdir()
    marker = out / "keep"
    marker.write_text("original evidence")
    result = run_import(importer, "gmrs_2022", raw, out)
    assert result.returncode != 0
    assert marker.read_text() == "original evidence"


def test_unknown_layout_is_not_guessed(importer, tmp_path):
    raw = tmp_path / "source.dat"
    source(raw, "gmrs_2022")
    out = tmp_path / "rejected"
    result = run_import(importer, "guess_columns", raw, out)
    assert result.returncode != 0
    assert not out.exists()


@pytest.mark.parametrize("mass,proton", [
    ("939 extra", "938"), ("939,938", "938"), ("939 938", "938"),
    ("1*939", "938"), ("939", "938 extra"),
])
def test_mass_arguments_are_single_numeric_tokens(importer, tmp_path, mass, proton):
    raw = tmp_path / "source.dat"
    source(raw, "gmrs_2022")
    out = tmp_path / "rejected"
    result = run_import(importer, "gmrs_2022", raw, out, mass, proton)
    assert result.returncode != 0
    assert not out.exists()


TOTAL_ENERGY_LAYOUTS = {
    "user_core_energy_density": (7, 0, 1, 2),
    "user_crust_energy_density": (5, 1, 2, 3),
    "n_eps_pressure": (3, 0, 1, 2),
}


def total_energy_source(path, layout):
    columns, density, energy, pressure = TOTAL_ENERGY_LAYOUTS[layout]
    lines = []
    for n, e, p in zip((".01", ".02", ".03", ".04", ".05"),
                       ("9.4", "18.82", "28.26", "37.72", "47.2"),
                       (".01", ".04", ".09", ".16", ".25")):
        row = ["0"] * columns
        row[density], row[energy], row[pressure] = n, e, p
        lines.append(" ".join(row))
    path.write_text("# total epsilon=939*n+100*n^2; P=100*n^2\n" +
                    "\n".join(lines) + "\n")


@pytest.mark.parametrize("layout", TOTAL_ENERGY_LAYOUTS)
@pytest.mark.parametrize("header_mass", ["939", "900"])
def test_total_energy_is_not_shifted_by_header_normalization(importer, tmp_path, layout, header_mass):
    raw = tmp_path / "source.dat"
    total_energy_source(raw, layout)
    digest = hashlib.sha256(raw.read_bytes()).hexdigest()
    out = tmp_path / "projection"
    result = run_import(importer, layout, raw, out, mass=header_mass)
    assert result.returncode == 0, result.stdout + result.stderr
    rows = list(csv.DictReader((out / "normalized_eos.csv").open()))
    assert len(rows) == 5
    for row, n, e, p in zip(rows, (.01, .02, .03, .04, .05),
                          (9.4, 18.82, 28.26, 37.72, 47.2), (.01, .04, .09, .16, .25)):
        assert float(row["n_b_fm3"]) == n
        assert float(row["energy_mev_fm3"]) == pytest.approx(e, rel=2e-15)
        assert float(row["pressure_mev_fm3"]) == pytest.approx(p, rel=2e-15)
        assert float(row["mu_b_mev"]) == pytest.approx(939+200*n, rel=2e-15)
    report = json.loads((out / "IMPORT.json").read_text())
    assert report["energy_representation"] == "total_energy_density_mev_fm3"
    assert report["rest_mass_offset_mev"] == 0
    assert report["neutron_header_mass_mev"] == float(header_mass)
    assert report["raw_rows"] == report["selected_rows"] == report["rows"] == 5
    assert report["excluded_rows"] == 0
    assert report["density_selection_applied"] is False
    assert report["requested_density_domain_fm3"] is None
    assert report["source_columns_one_based"] == [x+1 for x in TOTAL_ENERGY_LAYOUTS[layout][1:]]
    evidence = list(csv.DictReader((out / "SOURCE_VALUES.csv").open()))
    assert float(evidence[0]["energy_mev_fm3"]) == 9.4
    assert hashlib.sha256(raw.read_bytes()).hexdigest() == digest


@pytest.mark.parametrize("layout", TOTAL_ENERGY_LAYOUTS)
@pytest.mark.parametrize("bounds,expected_n,excluded_low,excluded_high", [
    ((".02", ".04"), [.02, .03, .04], 1, 1),
    (("0", ".03"), [.01, .02, .03], 0, 2),
])
def test_explicit_density_selection_preserves_rows_and_receipt(
        importer, tmp_path, layout, bounds, expected_n, excluded_low, excluded_high):
    raw = tmp_path / "source.dat"
    total_energy_source(raw, layout)
    out = tmp_path / "projection"
    result = run_import(importer, layout, raw, out, bounds=bounds)
    assert result.returncode == 0, result.stdout + result.stderr
    rows = list(csv.DictReader((out / "normalized_eos.csv").open()))
    assert [float(row["n_b_fm3"]) for row in rows] == expected_n
    report = json.loads((out / "IMPORT.json").read_text())
    assert report["raw_rows"] == 5
    assert report["selected_rows"] == report["rows"] == 3
    assert report["excluded_rows"] == 2
    assert report["excluded_lower_rows"] == excluded_low
    assert report["excluded_upper_rows"] == excluded_high
    assert report["source_row_start_zero_based"] == excluded_low
    assert report["source_row_stop_exclusive"] == 5-excluded_high
    assert report["density_selection_applied"] is True
    assert report["requested_density_domain_fm3"] == list(map(float, bounds))
    assert report["raw_density_domain_fm3"] == [.01, .05]
    assert report["selected_density_domain_fm3"] == [expected_n[0], expected_n[-1]]
    assert report["physical_acceptance_asserted"] is False
    assert len(list(csv.DictReader((out / "SOURCE_VALUES.csv").open()))) == 5


@pytest.mark.parametrize("bounds", [
    (".05", ".02"), ("NaN", ".03"), ("-.01", ".04"),
    (".02", ".02"), (".015", ".035"), (".01 extra", ".05"),
])
def test_invalid_or_insufficient_density_window_is_rejected(importer, tmp_path, bounds):
    raw = tmp_path / "source.dat"
    total_energy_source(raw, "user_core_energy_density")
    out = tmp_path / "rejected"
    result = run_import(importer, "user_core_energy_density", raw, out, bounds=bounds)
    assert result.returncode != 0
    assert not out.exists()


def test_selection_does_not_hide_invalid_excluded_rows(importer, tmp_path):
    raw = tmp_path / "source.dat"
    total_energy_source(raw, "user_core_energy_density")
    raw.write_text(raw.read_text().replace(".01 9.4", "NaN 9.4"))
    out = tmp_path / "rejected"
    result = run_import(importer, "user_core_energy_density", raw, out, bounds=(".02", ".05"))
    assert result.returncode != 0
    assert not out.exists()


@pytest.mark.parametrize("header_mass", ["939", "900"])
def test_secondary_bbp_cgs_units_and_explicit_fm3_selection(importer, tmp_path, header_mass):
    raw = tmp_path / "source.dat"
    raw.write_text("# rho[g/cm^3] P[dyn/cm^2] n[cm^-3]\n"
                   "1e12 1e30 1e36\n2e12 4e30 2e36\n3e12 9e30 3e36\n"
                   "4e12 16e30 4e36\n5e12 25e30 5e36\n")
    digest = hashlib.sha256(raw.read_bytes()).hexdigest()
    out = tmp_path / "projection"
    result = run_import(importer, "bbp_cgs_2008", raw, out, mass=header_mass,
                        bounds=(".002", ".004"))
    assert result.returncode == 0, result.stdout + result.stderr
    rows = list(csv.DictReader((out / "normalized_eos.csv").open()))
    # Independent decimal unit oracle: 1 g/cm3 -> 5.609588603804452e-13 MeV/fm3,
    # 1 dyn/cm2 -> 6.241509074460763e-34 MeV/fm3; neither is a rest-mass addition.
    for row, n, epsilon, pressure in zip(
            rows, (.002, .003, .004),
            (1.1219177207608904, 1.6828765811413356, 2.2438354415217808),
            (.002496603629784305, .005617358167014687, .00998641451913722)):
        assert float(row["n_b_fm3"]) == pytest.approx(n, rel=2e-15)
        assert float(row["energy_mev_fm3"]) == pytest.approx(epsilon, rel=2e-15)
        assert float(row["pressure_mev_fm3"]) == pytest.approx(pressure, rel=2e-15)
        assert float(row["mu_b_mev"]) == pytest.approx((epsilon+pressure)/n, rel=2e-15)
    receipt = json.loads((out / "IMPORT.json").read_text())
    assert receipt["raw_rows"] == 5 and receipt["selected_rows"] == 3
    assert receipt["excluded_lower_rows"] == receipt["excluded_upper_rows"] == 1
    assert receipt["source_columns_one_based"] == [3, 1, 2]
    assert receipt["source_units"] == ["cm^-3", "g cm^-3", "dyn cm^-2"]
    assert receipt["unit_conversion"] == "CGS_TO_MEV_FM3_EXACT_SI_CONSTANTS"
    assert receipt["rest_mass_offset_mev"] == 0
    assert receipt["density_selection_applied"] is True
    assert receipt["physical_acceptance_asserted"] is False
    assert len(list(csv.DictReader((out / "SOURCE_VALUES.csv").open()))) == 5
    assert hashlib.sha256(raw.read_bytes()).hexdigest() == digest


def test_gogny_full_explicit_split_does_not_interpolate_or_relabel(importer, tmp_path):
    raw = tmp_path / "source.dat"
    raw.write_text("# crust\n.01 1 .01 1 2 0 0 0\n.02 2 .04 1 2 0 0 0\n"
                   "#Core\n.03 3 .09 1 2 0 0 0\n.04 4 .16 1 2 0 0 0\n.05 5 .25 1 2 0 0 0\n")
    out = tmp_path / "projection"
    result = run_import(importer, "gogny_full_2020", raw, out, bounds=(".03", ".05"))
    assert result.returncode == 0, result.stdout + result.stderr
    receipt = json.loads((out / "IMPORT.json").read_text())
    assert receipt["raw_rows"] == 5 and receipt["selected_rows"] == 3
    assert receipt["source_row_start_zero_based"] == 2
    assert receipt["source_row_stop_exclusive"] == 5
    assert receipt["rest_mass_offset_mev"] == 939
    assert receipt["source_columns_one_based"] == [1, 2, 3]
    assert receipt["phase_labels_reconstructed"] is False
    rows = list(csv.DictReader((out / "normalized_eos.csv").open()))
    assert [float(row["n_b_fm3"]) for row in rows] == [.03, .04, .05]
    assert [float(row["energy_mev_fm3"]) for row in rows] == pytest.approx([28.26, 37.72, 47.2], rel=2e-15)


@pytest.mark.parametrize("bad", ["NaN 1e30 1e36", "1e12 1e30 0", "1e12 1e30 1e36 0"])
def test_secondary_bbp_cgs_malformed_rows_rejected_even_when_excluded(importer, tmp_path, bad):
    raw = tmp_path / "source.dat"
    raw.write_text(bad + "\n2e12 4e30 2e36\n3e12 9e30 3e36\n4e12 16e30 4e36\n")
    out = tmp_path / "rejected"
    result = run_import(importer, "bbp_cgs_2008", raw, out, bounds=(".002", ".004"))
    assert result.returncode != 0
    assert not out.exists()
