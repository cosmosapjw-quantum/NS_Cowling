"""Immutable numerical evidence: source-callback Fortran vs PCHIP and DOP853.

No physical EOS input, Maxwell transfer or slow/rapid conversion rule was
modified. This is a bounded numerical comparison, not an exact error bound.
"""
import hashlib,json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"data/cowling_direct_source_callback_20261119"
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def load(name):
    return json.loads((DATA/name).read_text())


def test_reconstructed_eos_binary_identity_preserved_and_reference_locked():
    measured=load("RESULT.json")
    assert measured["original_derived_eos_zip_sha256"]==sha(SOURCE)
    assert measured["original_derived_eos_zip_sha256"]=="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
    assert measured["phasewise_sharp_interfaces_unchanged"]
    assert measured["no_new_density_interpolation_between_phases"]
    assert len(measured["cases"])==3


def test_old_native_PCHIP_output_bitwise_numeric_unchanged_by_source_callback():
    original=load("RESULT_BEFORE_RESTART.json")["cases"]
    latest=load("RESULT.json")["cases"]
    for case in latest:
        for key in ("pchip_8193","pchip_16385"):
            x=original[case]["variants"][key]
            y=latest[case]["variants"][key]
            assert x["status"]==y["status"]
            if y["status"]=="NATIVE_ROOT_ADMITTED":
                assert y["frequency_hz"]==x["frequency_hz"]
                assert y["delta_against_DOP853_hz"]==x["delta_against_DOP853_hz"]


def test_direct_source_accurate_on_physical_p1_and_retains_W1V2():
    z=load("RESULT.json")["cases"]
    for name,record in z.items():
        local=record["variants"]["direct_source_129"]
        assert local["status"]=="NATIVE_ROOT_ADMITTED"
        assert abs(local["delta_against_DOP853_hz"])<.005
        assert local["relative_join_mismatch"]<1e-5
        assert local["source_callback_rhs_calls"]>1000
        assert len(local["tabulation_points_per_phase"])==2
        assert all(k==33 for k in local["tabulation_points_per_phase"])
    assert z["SPG_M3_CSS_0p7_slow_p1"]["variants"]["pchip_8193"]["status"]=="NATIVE_ROOT_REJECTED"
    assert z["SPG_M3_CSS_0p7_slow_p1"]["variants"]["pchip_16385"]["status"]=="NATIVE_ROOT_REJECTED"
    problematic=z["SPG_M3_CSS_1p4_slow_W1V2"]["variants"]
    assert (problematic["direct_source_129"]["nodes_W"],
            problematic["direct_source_129"]["nodes_V"])==(1,2)
    assert abs(problematic["pchip_8193"]["delta_against_DOP853_hz"])>.05
    assert abs(problematic["pchip_16385"]["delta_against_DOP853_hz"])>.03


def test_rapid_conversion_direct_source_against_python_DOP853():
    r=load("RAPID.json")
    assert r["junction_conditions_unchanged"] is True
    assert r["conversion"]=="rapid"
    native=r["modes"]["direct_source"]
    old=r["modes"]["pchip"]
    oracle=r["oracle"]
    assert native["status"]==old["status"]=="ADMITTED"
    assert (native["W_nodes"],native["V_nodes"])==(1,1)
    assert abs(native["frequency_hz"]-oracle["frequency_hz"])<.001
    assert abs(old["frequency_hz"]-oracle["frequency_hz"])>.01
    assert native["traction_residual"]<1e-9
    assert native["kinematic_residual"]<1e-9
    assert native["join"]<1e-5
    assert native["callback_rhs_count"]>1000


def test_strong_knot_restarts_repair_observed_dvode_admission_sensitivity():
    old=load("TOLERANCE_BEFORE_RESTART.json")
    new=load("TOLERANCE.json")
    assert old["original_contract_sha256"]==new["original_contract_sha256"]
    assert old["original_archive_sha256"]==new["original_archive_sha256"]
    assert len(old["trials"])==len(new["trials"])==3
    assert sum(x["status"]=="ADMITTED" for x in old["trials"])==1
    assert all(x["status"]=="ADMITTED" for x in new["trials"])
    values=[x["frequency_hz"] for x in new["trials"]]
    assert max(values)-min(values)<.005
    assert all(x["join"]<1e-5 for x in new["trials"])
    assert all((x["W_nodes"],x["V_nodes"])==(1,1) for x in new["trials"])
    assert all(x["source_rhs_evaluations"]>1000 for x in new["trials"])
