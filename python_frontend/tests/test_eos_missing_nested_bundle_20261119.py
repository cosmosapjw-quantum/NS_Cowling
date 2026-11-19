"""Independent content, source tier and numerical acceptance tests for nested EOS."""
import hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"data/eos_missing_data_recalc_20261119"
MASTER=DATA/"artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
NAMES=[*(f"SPG_M{i}" for i in range(1,6)),
       *(f"GPPVA_{m}" for m in
         ("DD2","DDME2","TW","FSU2","FSU2R","FSU2H","NL3wrL55","TM1e"))]
def sha(raw):return hashlib.sha256(raw).hexdigest()
def load(path):return json.loads(Path(path).read_text())
def inner_zip(name):
    with zipfile.ZipFile(MASTER) as z:
        return io.BytesIO(z.read(name))


def test_all_13_energy_primary_native_recomputed_are_not_overwritten_original():
    from compose_eos import read_compose
    report=load(DATA/"CORRECTED13_RESULT.json")
    assert sorted(report["passed"])==sorted(NAMES)
    assert report["blocked"]==[]
    for name in NAMES:
        rec=report["model_outputs"][name]
        assert rec["readback_verified"]
        assert rec["source_original_knots_preserved"]
        assert rec["changed_fields"]==["Q1","Q3"]
        a=read_compose(ROOT/rec["source"])
        b=read_compose(ROOT/rec["dest"])
        np.testing.assert_array_equal(a.nb,b.nb)
        np.testing.assert_array_equal(a.Q7,b.Q7)
        np.testing.assert_array_equal(a.Iphase,b.Iphase)
        assert np.any(a.Q1!=b.Q1)
        assert a.n==b.n


def test_six_physically_recomputed_native_crust_candidates():
    spg=load(DATA/"reconstructed_spg_fresh17/SUMMARY.json")
    assert spg["status"]=="COMPLETED_SCOPED"
    assert spg["source_calibrated"] is False
    for n in range(1,6):
        r=spg["models"][f"SPG_M{n}"]
        assert r["rows"]==17
        assert r["status"]=="GENERATED_STATIONARY_NATIVE_COMPOSITION_READBACK_VERIFIED"
    F=load(DATA/"reconstructed_gppva_fsu2_fresh52/RECALC_REPORT.json")
    assert F["status"]=="FINISHED_NATIVE_FORTRAN_REGENERATED_AND_RELOADED"
    assert F["density_count"]==52
    assert F["byte_identical_to_prior"]


def test_pk1_is_newly_native_fortran_reproduced_all_4097_densities():
    receipt=load(DATA/"pk1_core_fresh/RECALC_REPORT.json")
    assert receipt["status"]=="DIRECT_FORTRAN_PK1_REEVALUATED"
    assert receipt["density_grid_count"]==4097
    assert receipt["validation"]["physical_state_checks"]
    assert receipt["bytes_equal_to_previous"]
    assert sha((DATA/"pk1_core_fresh/table/compose/eos.thermo").read_bytes())==receipt["fresh_sha256"]


def test_bsk22_analytic_surface_anchor_and_selected_outer_are_distinct():
    z=load(DATA/"bsk22_outer_and_surface/RECALC_REPORT.json")
    assert z["source_outer_rows"]==217
    assert z["independent_surface_rows"]==256
    assert z["surface_status"]=="FRESH_ANALYTIC_FIRSTLAW_EVALUATION_FROM_SOURCE_ANCHOR"
    assert z["exact_anchor_energy_relative_error"]<1e-13
    assert z["outer_source_role"].find("NOT REGENERATED ETFSI")>=0
    lower=load(DATA/"bsk22_lower_inner_fresh/RECEIPT.json")
    assert lower["status"]=="PUBLISHED_ANALYTIC_REPRESENTATION_VERIFIED_SCOPED"
    assert lower["versions"]["source_fit"]["points"]==129
    assert lower["versions"]["energy_derivative"]["points"]==129


def test_all_10_native_tov_stars_and_one_slow_p1_node_uncertainty_preserved():
    stars=load(DATA/"spg_pk1_complete_stars/RESULT.json")["models"]
    assert sorted(stars)==sorted(f"SPG_M{i}" for i in range(1,6))
    cases=[]
    for name,rows in stars.items():
        for method in ("maxwell","css"):
            case=rows[method]
            assert case["status"]=="NATIVE_1P4_EXECUTED"
            assert case["outer_original_217_knots_unchanged"]
            solver=case["native_star"]
            assert abs(solver["native_TOV"]["M_Msun"]-1.4)<3e-6
            for conversion in ("slow","rapid"):
                assert "f" in solver["native_cowling"][conversion]["accepted_modes"]
            cases.append((name,method,solver["status"]))
    complete=[x for x in cases if x[2]=="FIXED_1P4_NATIVE_TOV_F_P1_BOTH_CONVERSIONS_VERIFIED"]
    partial=[x for x in cases if x[2]=="FIXED_1P4_NATIVE_TOV_PARTIAL_MODE_IDENTITY"]
    assert len(complete)==9
    assert partial==[("SPG_M3","css","FIXED_1P4_NATIVE_TOV_PARTIAL_MODE_IDENTITY")]


def test_nuda_original_licence_and_extention_not_miscounted_as_24_models():
    z=load(DATA/"external_original/NUDA_SOURCE_AUDIT.json")
    assert z["nominal_independent_source_models"]==23
    assert z["total"]==24
    assert z["additional_nonindependent_extension"]=="2020-MVCD-D1MS-full.dat"
    assert z["model_groups"]=={"BSk_2018":4,"GMRS_2022":15,
                                "Gogny_MVCD":4,"Negele_Vautherin_1973":1}
    assert z["license"]=="CC BY-NC-ND 4.0"
    for name,record in z["original_datasets"].items():
        assert record["raw_bytes_unchanged"]
        assert not record["physical_model_recomputed"]
        assert sha((DATA/"external_original"/name).read_bytes())==record["sha256"]


def test_nested_zip_contains_all_seven_verified_independent_members():
    assert MASTER.is_file()
    with zipfile.ZipFile(MASTER) as outer:
        assert outer.testzip() is None
        manifest=json.loads(outer.read("BUNDLE_MANIFEST.json"))
        assert len(manifest["inner_archives"])==7
        for entry in manifest["inner_archives"]:
            blob=outer.read(entry["name"])
            assert sha(blob)==entry["sha256"]
            with zipfile.ZipFile(io.BytesIO(blob)) as inner:
                assert inner.testzip() is None
                files=json.loads(inner.read("INTERNAL_PROVENANCE_AND_SHA256.json"))
                assert files["count"]==entry["n_files"]
                assert set(files["file_sha256"]).issubset(set(inner.namelist()))
                for path,h in files["file_sha256"].items():
                    assert sha(inner.read(path))==h
        assert manifest["results"]["original_13_recalculated_count"]==13
        assert manifest["results"]["strict_slow_rapid_f_p1_complete"]==9
        assert manifest["results"]["external_original_named_models"]==23
        assert "README_KO.md" in outer.namelist()


def test_outer_and_recomputed_13_are_never_mislabeled_original():
    with zipfile.ZipFile(inner_zip("01_CORRECTED_13_SPG_GPPVA.zip")) as z:
        contents=z.namelist()
        for name in NAMES:
            assert f"ORIGINAL_FULL_SOURCE_UNMODIFIED/{name}/eos.thermo" in contents
            assert f"FRESH_ENERGY_PRIMARY_CORRECTED_SELECTED_PHASE/{name}/eos.thermo" in contents
            orig=z.read(f"ORIGINAL_FULL_SOURCE_UNMODIFIED/{name}/eos.thermo")
            updated=z.read(f"FRESH_ENERGY_PRIMARY_CORRECTED_SELECTED_PHASE/{name}/eos.thermo")
            assert sha(orig)!=sha(updated)


def test_external_zip_bytes_preserved_unmodified_and_contains_license():
    with zipfile.ZipFile(inner_zip("06_NUDA_23_ORIGINAL_MODELS_PLUS_D1MS_FULL.zip")) as z:
        audit=json.loads(z.read("NUDA_SOURCE_AUDIT.json"))
        assert len(audit["original_datasets"])==24
        assert "UPSTREAM_SOURCE_LICENSE_CC_BY_NC_ND_4_0.txt" in z.namelist()
        for name,row in audit["original_datasets"].items():
            arc="UPSTREAM_UNMODIFIED_NUCLEARDATAPY_1_0_2/"+name
            assert sha(z.read(arc))==row["sha256"]
