"""Physical source-locked evidence for optional hybrid Cowling cs² backend."""
from __future__ import annotations
import hashlib,json
from pathlib import Path
import pytest

from scripts.summarize_hybrid_accuracy_20261120 import choose_verified_backend

ROOT=Path(__file__).resolve().parents[2]
D=ROOT/"data/cowling_hybrid_geometry_20261120"
EOS=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"


def load(relative):
    return json.loads((D/relative).read_text())


def test_hybrid_benchmark_uses_exact_frozen_source_archive():
    source=load("RESULT.json")
    actual=hashlib.sha256(EOS.read_bytes()).hexdigest()
    assert actual==source["frozen_EOS_bundle_sha256"]
    assert actual=="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
    assert set(source["cases"])=={
       "M3_CSS_0p7_SLOW_P1","M2_MAXWELL_1p4_RAPID_P1"}
    assert source["status"]=="BOUNDED_BENCHMARK_RECORDED_NOT_GENERAL_PERFORMANCE_PROOF"


def test_exact_cs2_recovers_low_mass_mode_while_32769_pchip_still_rejects():
    z=load("RESULT.json")["cases"]["M3_CSS_0p7_SLOW_P1"]
    a=z["backend_variants"]
    assert all(a[f"PCHIP_{n}"]["status"]=="NATIVE_ROOT_REJECTED"
           for n in (8193,16385,32769))
    for n in (2049,4097):
        row=a[f"HYBRID_CS2_{n}"]
        assert row["status"]=="MATCHED_REFERENCE_PASS"
        assert row["nodes"]==[1,1]
        assert abs(row["reference_delta_hz"])<.005
        assert row["source_restart_count"]>30
        assert row["source_callback_calls"]>1000
        assert row["join_relative"]<1e-5
    direct=a["DIRECT_SOURCE_129"]
    candidate=a["HYBRID_CS2_2049"]
    assert direct["status"]=="MATCHED_REFERENCE_PASS"
    assert candidate["total_seconds"]<direct["total_seconds"]


def test_rapid_mode_accuracy_cost_is_not_misrepresented_as_global_speedup():
    z=load("RESULT.json")["cases"]["M2_MAXWELL_1p4_RAPID_P1"]
    a=z["backend_variants"]
    assert a["PCHIP_8193"]["status"]=="MATCHED_REFERENCE_PASS"
    assert abs(a["PCHIP_8193"]["reference_delta_hz"])>.005
    assert abs(a["PCHIP_16385"]["reference_delta_hz"])<.005
    assert abs(a["PCHIP_32769"]["reference_delta_hz"])<.005
    assert a["HYBRID_CS2_4097"]["status"]=="MATCHED_REFERENCE_PASS"
    assert abs(a["HYBRID_CS2_4097"]["reference_delta_hz"])<.001
    assert a["HYBRID_CS2_4097"]["total_seconds"]<a["DIRECT_SOURCE_129"]["total_seconds"]


def test_empirical_selector_retains_nonadmitted_accuracy_scopes():
    raw=load("RESULT.json")["cases"]
    report=load("DECISIONS.json")
    assert report["input_benchmark_sha256"]==hashlib.sha256(
        (D/"RESULT.json").read_bytes()).hexdigest()
    m3=raw["M3_CSS_0p7_SLOW_P1"]
    m2=raw["M2_MAXWELL_1p4_RAPID_P1"]
    assert choose_verified_backend(m3,.005)["selected"]=="HYBRID_CS2_2049"
    assert choose_verified_backend(m3,.001)["selected"] is None
    assert choose_verified_backend(m2,.05)["selected"]=="PCHIP_8193"
    assert choose_verified_backend(m2,.005)["selected"]=="PCHIP_16385"
    assert choose_verified_backend(m2,.001)["selected"]=="HYBRID_CS2_4097"
    for bad in (-1,0,float("nan"),float("inf"),True):
        with pytest.raises(ValueError):
            choose_verified_backend(m2,bad)
