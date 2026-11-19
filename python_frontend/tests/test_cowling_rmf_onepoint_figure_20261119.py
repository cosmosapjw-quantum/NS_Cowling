"""Figure source identity, CSV exact deltas and independent F90 parity checks."""
import csv,hashlib,json
from pathlib import Path
from PIL import Image
import pytest

ROOT=Path(__file__).resolve().parents[2]
D=ROOT/"data/cowling_original_rmf_onepoint_20261119"
FIG=D/"figures"

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def test_plot_figure_bound_to_exact_three_original_source_knots():
    outcome=read(FIG/"FIGURE_RECEIPT.json")
    src=D/"RESULT.json"
    assert outcome["original_source_sha256"]==sha(src)
    assert outcome["no_input_interpolation"] is True
    assert outcome["input_source_has_original_meson_fields"] is False
    assert outcome["source_model_names"]==["FSU2","TM1e","DDME2"]
    for typ,row in outcome["figures"].items():
        p=ROOT/row["relative_path"]
        assert p.is_file()
        assert p.stat().st_size==row["bytes"]
        assert sha(p)==row["sha256"]
    with Image.open(FIG/"COMPOSE_vs_LOCAL_RMF_exact_onepoint.png") as img:
        assert img.format=="PNG"
        assert img.width>=1500 and img.height>=950
    svg=(FIG/"COMPOSE_vs_LOCAL_RMF_exact_onepoint.svg").read_text()
    assert "CompOSE" in svg and "Chemical-potential" in svg


def test_onepoint_csv_uses_original_source_hash_and_correct_difference_sign():
    z=read(D/"RESULT.json")["models"]
    with (FIG/"COMPOSE_vs_LOCAL_RMF_exact_onepoint.csv").open() as f:
        rows=list(csv.DictReader(f))
    assert [r["model"] for r in rows]==["FSU2","TM1e","DDME2"]
    for row in rows:
        name=row["model"];s=z[name]["source"];local=z[name]["local_beta_equilibrium"]
        assert float(row["n_B_original_fm3"])==s["original_source_n_B_fm3"]
        assert float(row["P_CompOSE_MeV_fm3"])==s["original_pressure_MeV_fm3"]
        assert float(row["P_local_MeV_fm3"])==local["pressure"]
        assert float(row["DeltaP_relative_percent"])==pytest.approx(
            100*(local["pressure"]/s["original_pressure_MeV_fm3"]-1),abs=1e-12)
        assert float(row["Deltaepsilon_MeV_fm3"])==pytest.approx(
            local["epsilon"]-s["original_epsilon_MeV_fm3"],abs=1e-12)
        assert float(row["Deltamu_local_minus_Q3_MeV"])==pytest.approx(
            local["muB"]-s["original_mu_B_Q3_MeV"],abs=1e-12)
        assert float(row["Deltamu_local_minus_Euler_MeV"])==pytest.approx(
            local["muB"]-s["original_mu_B_Euler_MeV"],abs=1e-12)
        assert row["source_has_original_meson_field_data"]=="False"


def test_five_native_RMF_original_onepoint_independent_fortran_parities():
    z=read(D/"FORTRAN_COMPONENT_PARITY.json")
    assert z["status"]=="ALL_FIVE_LOCAL_STATES_INDEPENDENT_NATIVE_F90_COMPONENT_PARITY_VERIFIED"
    assert len(z["independent_cases"])==5
    assert all(x["compiled_with_bounds_check"] for x in [z["build"]])
    for case in z["independent_cases"].values():
        assert case["status"]=="INDEPENDENT_FORTRAN_PYTHON_TERM_PARITY_PASSED"
        assert case["local_fields_not_original_provider"]
        assert case["maximum_component_delta_MeV_fm3"]<1e-8
        assert case["maximum_total_and_field_delta"]<1e-8
        assert max(abs(x) for x in case["Fortran_EOM_fm3"].values())<1e-7
        assert len(case["termwise_fortran_minus_python"])>=18
    for source,expected in z["build"]["source_sha256"].items():
        assert sha(ROOT/source)==expected



def test_full_density_compose_curve_preserves_213_exact_source_points():
    source=ROOT/"data/cowling_rmf_reconciliation_20261119/EXACT_DENSITY_AUDIT.json"
    receipt=read(FIG/"DENSITY_TREND_RECEIPT.json")
    assert receipt["source_SHA256"]==sha(source)
    assert receipt["no_source_interpolation_or_physics_fitting"] is True
    assert receipt["model_count"]=={"FSU2":75,"TM1e":63,"DDME2":75}
    assert receipt["total_rows"]==213
    old=read(source)["models"]
    csvfile=FIG/"COMPOSE_vs_LOCAL_RMF_original_density_trend.csv"
    with csvfile.open() as stream:
        rows=list(csv.DictReader(stream))
    assert len(rows)==213
    for name in ("FSU2","TM1e","DDME2"):
        actual=[x for x in rows if x["model"]==name]
        original=old[name]["original_beta_equilibrium_fraction_comparison"]["homogeneous_rows"]
        assert len(actual)==len(original)
        for a,b in zip(actual,original):
            assert float(a["n_B_fm3"])==b["density_fm3"]
            assert float(a["P_CompOSE_MeV_fm3"])==b["pressure_provider"]
            assert float(a["P_local_MeV_fm3"])==b["pressure_model"]
            assert float(a["Relative_deltaP_percent"])==pytest.approx(
                100*(b["pressure_model"]/b["pressure_provider"]-1),abs=1e-12)
            assert float(a["Delta_mu_Q3_minus_Euler_CompOSE_MeV"])==pytest.approx(
                b["mu_Q3_provider"]-b["mu_Euler_provider"],abs=1e-12)
            assert a["original_meson_components_available"]=="False"
    for filename,meta in receipt["figures"].items():
        p=FIG/filename
        assert sha(p)==meta["sha256"]
        assert p.stat().st_size==meta["size_bytes"]


def test_full_density_plot_is_readable_and_includes_each_model():
    with Image.open(FIG/"COMPOSE_vs_LOCAL_RMF_original_density_trend.png") as im:
        assert im.width>=1500 and im.height>=950
        assert im.format=="PNG"
    svg=(FIG/"COMPOSE_vs_LOCAL_RMF_original_density_trend.svg").read_text()
    assert all(x in svg for x in ("FSU2","TM1e","DDME2"))



def test_onepoint_nested_evidence_report_zip_source_hashes_and_crc():
    import io,zipfile
    zpath=D/"artifacts/COWLING_RMF_ONEPOINT_COMPOSE_AUTHOR_COMPARE_20261119.zip"
    receipt=read(D/"artifacts/REPORT_PACKAGE_RECEIPT.json")
    assert sha(zpath)==receipt["sha256"]
    assert zpath.stat().st_size==receipt["size_bytes"]
    assert receipt["original_source_binaries_included"] is False
    assert receipt["all_CRC_SHA_tests_passed"] is True
    with zipfile.ZipFile(zpath) as z:
        assert z.testzip() is None
        index=json.loads(z.read("INTERNAL_SOURCE_IDENTITY_SHA256.json"))
        assert index["original_shen_author_archive_not_included"] is True
        assert index["original_gppva_compose_archives_not_included"] is True
        assert index["original_Pais_meson_fields_or_generator_obtained"] is False
        assert index["exact_original_density_trend_rows"]==213
        assert len(index["members"])==27
        for name,record in index["members"].items():
            assert name in z.namelist()
            assert hashlib.sha256(z.read(name)).hexdigest()==record["sha256"]
            assert len(z.read(name))==record["size_bytes"]
        assert all("shen_EOS4_original_eos4.t00.gz" not in x for x in z.namelist())
