"""Reproducible bounded scientific receipts for COWLING-EOS-SOURCE-CALIBRATION-04.

These checks validate persisted source identities and observed numerical
inequalities. They supplement, not replace, actual Fortran integration tests.
"""
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"data/cowling_eos_source_calibration_20261119"


def load(path):
    return json.loads((DATA/path).read_text())


def test_independent_official_archive_sha_identical_to_zenodo():
    z=load("zenodo/ZENODO_INDEPENDENT_ARCHIVE_READBACK.json")
    assert set(z)=={"FSU2","TM1e","DDME2"}
    for name,x in z.items():
        assert x["same_bytes_as_compose_download"] is True
        assert x["compose_checksum_file_mismatch"] is True
        assert len(x["observed_zip_sha256"])==64
        assert len(x["zenodo_md5"])==32


def test_same_named_functional_not_mistaken_for_identical_provider_eos():
    z=load("source_audit/RESULT.json")["models"]
    assert set(z)=={"FSU2","TM1e","DDME2"}
    for name in ("FSU2","TM1e"):
        x=z[name]["original_homogeneous_core_same_density_comparison"]
        assert x["classification"]=="SAME_NAMED_PARAMETER_FAMILY_BUT_SOURCE_NOT_NUMERICALLY_EQUIVALENT"
        assert x["pressure_abs_rel_max"]>.01
        assert x["n_evaluated"]>40
    dd=z["DDME2"]["original_homogeneous_core_same_density_comparison"]
    assert dd["pressure_abs_rel_max"]<1e-4
    assert dd["classification"]=="HIGH_NUMERIC_AGREEMENT_NOT_BIT_IDENTITY"


def test_break_waypoints_do_not_create_spurious_frequency_improvement():
    experiment=load("sharp_waypoints/RESULT.json")
    count=0
    for case in experiment["cases"].values():
        for ns,receipt in case["mode_results"].items():
            trials=receipt["modes"]
            a=trials["baseline_continuous_DVODE"]
            b=trials["physical_break_waypoints"]
            for conversion in ("slow","rapid"):
                for label in ("f","p1"):
                    native=a[conversion][label]
                    forced=b[conversion][label]
                    assert native["status"]==forced["status"]=="CHECKED"
                    assert abs(native["frequency_hz"]-forced["frequency_hz"])<1e-5
                    assert forced["source_break_waypoint_count"]>0
                    count+=1
    assert count==4*3*2*2


def test_error_component_is_cowling_eos_coefficients_not_native_tov():
    x=load("resolution_separation/RESULT.json")
    source=x["cases"]["same_FSU2_maxwell"]
    py=source["python_background"]
    slow_p1=lambda n: py["coefficient_sampling"][str(n)]["modes"]["slow"]["p1"]["f_Hz"]
    ref=4897.917834807901 # independent Python DOP853 mode as checked in 03
    assert abs(slow_p1(32769)-ref)<.003
    assert abs(slow_p1(4097)-ref)>.08
    for _,run in source["native_background"].items():
        f=run["fine_cowling_samples_32769"]["modes"]["slow"]["p1"]["f_Hz"]
        assert abs(f-slow_p1(32769))<.001
    assert py["coefficient_sampling"]["4097"]["modes"]["slow"]["p1"]["status"]=="CHECKED"


def test_tm1e_strict_classification_unresolved_but_conditional_continued():
    x=load("tm1e_mode_identity/RESULT.json")
    slow=x["cases"]["slow"]
    assert len(slow)==7
    unresolved=[v for v in slow if v["continuation"] is not None
        and v["continuation"]["node_pair_safe_tracker"]["status"]=="UNRESOLVED"]
    assert len(unresolved)==2
    assert [round(v["nc_fm3"],2) for v in unresolved]==[.35,.38]
    for r in unresolved:
        c=r["continuation"]["conditional_shape_continuity_tracker"]
        assert c["status"]=="CONTINUED_WITH_V_NODE_CHANGE"
        assert c["pairs"][0]["absolute_overlap"]>=.995
    assert all(row["continuation"] is None or
          row["continuation"]["conditional_shape_continuity_tracker"]["status"]
          in ("CONTINUED","CONTINUED_WITH_V_NODE_CHANGE") for row in slow)


def test_tm1e_v_node_extra_zero_is_on_both_sides_and_reproduced():
    x=load("tm1e_mode_identity/NODE_LOCALIZATION.json")
    row=x["cases"]["0.36680854236573207"]
    interface=row["physical_interface_r_km"]
    for ns in ("160","450","900","1800"):
        r=row["results"][ns]
        assert (r["nodes_W"],r["nodes_V"])==(1,2)
        inner=r["phases"]["0"]["V_crossings_km"]
        outer=r["phases"]["1"]["V_crossings_km"]
        assert len(inner)==len(outer)==1
        assert inner[0]<interface<outer[0]
        assert 0<outer[0]-interface<.01
    assert all(r["nodes_V"]==1 for r in x["cases"]["0.38"]["results"].values())


def test_phase_local_repair_physically_differs_from_cross_seam_bridge():
    x=load("phase_local_fixed_mass/RESULT.json")["cases"]
    for name in ("FSU2","DDME2"):
        r=x[name]
        assert r["status"]=="FIXED_1P4_NATIVE_TOV_F_P1_BOTH_CONVERSIONS_VERIFIED"
        delta=r["phasewise_minus_cross_seam_firstlaw"]
        assert delta["R_km"]<-.02
        assert delta["mode_frequency_hz"]["slow"]["p1"]>40.
    tm=x["TM1e"]
    assert tm["status"]=="FIXED_1P4_NATIVE_TOV_PARTIAL_MODE_IDENTITY"
    assert tm["native_cowling"]["slow"]["status"]=="PARTIAL_NODE_RESOLVED_SPECTRUM"
    assert tm["native_cowling"]["rapid"]["status"]=="F_P1_COMPLETE"
