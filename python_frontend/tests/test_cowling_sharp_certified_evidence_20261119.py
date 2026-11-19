"""Pinned source-bearing evidence for sharp native/DOP853 precision upgrade."""
from pathlib import Path
import hashlib
import json
import math
import zipfile
import io

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"data/cowling_certified_modes_20261119"
ARCHIVE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"

def sha(data):
    return hashlib.sha256(data).hexdigest()

def get():
    return json.loads((DATA/"RESULT.json").read_text())


def test_source_derived_original_inputs_are_pinned_and_not_relabeled():
    data=get()
    assert data["original_parent_commit"]=="1e7823ecd62a2914d40c0948d3811410b1ea16da"
    assert data["immutable_derived_source_zip_sha256"]==sha(ARCHIVE.read_bytes())
    assert data["sources_are_source_derived_not_original_Pais_or_SPG"] is True
    assert data["status"]=="BOUNDED_REFERENCE_GATED_RESULTS_RECORDED"
    with zipfile.ZipFile(ARCHIVE) as top:
        inner=io.BytesIO(top.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))
        with zipfile.ZipFile(inner) as sub:
            receipts=json.loads(sub.read("INTERNAL_PROVENANCE_AND_SHA256.json"))
            for m,join in (("SPG_M3","css"),("SPG_M2","maxwell")):
                key=f"FRESH_5x2_WHOLE_STARS/{m}/{join}/REPRESENTATION.json"
                stored=sha(sub.read(key))
                assert stored==receipts["file_sha256"][key]
                expected=next(z["source_contract_sha256"] for n,z in
                    data["modes"].items() if f"{m}_{join}".lower() in n.lower())
                assert stored==expected


def test_m3_low_mass_native_failure_is_preserved_and_reference_recovers_mode():
    row=get()["modes"]["SPG_M3_CSS_0p7_slow_p1"]
    assert row["status"]=="REFERENCE_ONLY_NATIVE_NOT_CERTIFIED"
    assert row["selected_backend"]=="PYTHON_DOP853_DIRECT_EOS"
    assert abs(row["frequency_hz"]-3736.2089)<.03
    assert row["node_gate"]=="STRICT_P1"
    assert row["expected_mode_identification"]=="STRICT_CONFIRMED"
    assert not row["native_convergence"]["native_converged"]
    assert len(row["native_attempts"])==2
    assert all(a["status"]=="REJECTED" for a in row["native_attempts"])
    assert row["python_reference"]["relative_join_mismatch"]<1e-6


def test_m3_medium_mass_frequency_exists_but_mode_name_is_conditional():
    row=get()["modes"]["SPG_M3_CSS_1p4_slow_W1V2"]
    assert row["status"]=="REFERENCE_ONLY_NATIVE_NOT_CERTIFIED"
    assert row["node_gate"]=="W1_V2_CLASSIFICATION_UNRESOLVED"
    assert row["expected_mode_identification"]=="CONDITIONAL_NUMERIC_FAMILY_NOT_STRICT_P1"
    assert abs(row["frequency_hz"]-5329.2027)<.03
    assert row["native_convergence"]["reason"]=="COEFFICIENT_SAMPLING_NOT_CONVERGED"
    assert row["native_convergence"]["comparison"]["difference_hz"]>.5


def test_m2_medium_mass_has_native_source_coefficient_and_python_oracle_agreement():
    row=get()["modes"]["SPG_M2_MAXWELL_1p4_slow_p1"]
    assert row["status"]=="CERTIFIED_NATIVE_REFERENCE_AGREEMENT"
    assert row["node_gate"]=="STRICT_P1"
    assert row["native_convergence"]["native_converged"]
    assert row["native_convergence"]["comparison"]["eigenfunction_overlap"]>.999
    assert row["native_convergence"]["comparison"]["difference_hz"]<.05
    assert row["native_reference_difference_hz"]<.05
    assert abs(row["frequency_hz"]-5148.001)<.02
    assert row["python_reference"]["relative_join_mismatch"]<1e-6
