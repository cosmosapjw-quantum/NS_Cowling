"""Origin-isolated, source-sealed and full-crust gates for exhaustive EOS census."""
import importlib,json
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[2]

def mod():
    return importlib.import_module("all_origin_maxwell_css_20261118")

def test_full_origin_catalog_and_no_cross_tier_pairing():
    m=mod()
    catalog=m.catalog()
    src=catalog["sources"]
    assert len(src)==95
    pairs=m.enumerate_pairs(src)
    assert len(pairs)==635
    assert all(a["tier"]==b["tier"] for a,b in pairs)
    assert set(x["tier"] for x in src)=={"original","corrected","reconstructed"}
    assert len(set(z["source_id"] for z in src))==len(src)

def test_public_original_outer_crust_and_native_missing_outer():
    m=mod()
    d={x["source_id"]:x for x in m.catalog()["sources"]}
    for name in ("SPG_M1","GPPVA_DDME2"):
        assert d["public_original:"+name+":crust"]["outer_present"]
        assert d["public_corrected:"+name+":crust"]["outer_present"]
        assert not d["native_reconstructed:"+name+":crust"]["outer_present"]

def test_independent_external_crust_included_as_not_ready_without_thermo():
    m=mod()
    d={x["source_id"]:x for x in m.catalog()["sources"]}
    assert d["external_nuda:1973-Negele-Vautherin:crust"]["status"]=="BLOCKED_MISSING_PRESSURE"
    assert not d["external_nuda:1973-Negele-Vautherin:crust"]["full_outer_inner"]
    assert d["external_nuda:2020-MVCD-D1S:crust"]["core_independent_model"]
    assert any("BPS" in x["source_id"] for x in d.values())
    assert any("BBP" in x["source_id"] for x in d.values())

def test_cannot_count_css_only_or_maxwell_only_as_both():
    m=mod()
    assert m.classify_both({"css":"PASS","maxwell":"NO_CROSSING","whole":"PASS"})=="NOT_BOTH"
    assert m.classify_both({"css":"PASS","maxwell":"PASS","whole":"OUTER_ABSENT"})=="FINITE_ONLY_NOT_WHOLE"
    assert m.classify_both({"css":"PASS","maxwell":"PASS","whole":"PASS"})=="BOTH_WHOLE_PASS"

def test_fast_scalar_precheck_preserves_overlapping_possible_roots():
    m=mod()
    from scipy.interpolate import CubicHermiteSpline
    x=np.array([1.,2.,3.])
    a=CubicHermiteSpline(x,np.array([0.,1.,2.]),np.ones(3))
    b=CubicHermiteSpline(x,np.array([.5,.5,.5]),np.zeros(3))
    roots=m.bracket_roots(a,b,1.,3.)
    assert len(roots)==1 and abs(roots[0]-1.5)<1e-8

def test_frozen_campaign_seven_dual_finite_and_five_conditional_stars():
    from collections import Counter
    root=ROOT/"data/all_origin_maxwell_css_20261118"
    campaign=json.loads((root/"CAMPAIGN.json").read_text())
    stars=json.loads((root/"whole_stars/CAMPAIGN.json").read_text())
    c=Counter(x["both_status"] for x in campaign["cases"])
    assert len(campaign["cases"])==635
    assert c["FINITE_ONLY_NOT_WHOLE"]==7
    assert c.get("BOTH_WHOLE_PASS",0)==0
    assert stars["both_conditional_whole"]==5
    for crust in ("SPG_M1","SPG_M2","SPG_M3","SPG_M4","SPG_M5"):
        r=json.loads((root/"whole_stars"/crust/"RESULT.json").read_text())
        assert all(v["status"]=="CONDITIONAL_WHOLE_NS_VERIFIED_TOV" for v in r["methods"].values())
        assert not r["outer_reference_sha256"]==""


def test_original_and_corrected_classes_have_no_dual_finite_result():
    data=json.loads((ROOT/"data/all_origin_maxwell_css_20261118/CAMPAIGN.json").read_text())
    for row in data["cases"]:
        if row["tier"]!="reconstructed":
            assert row["both_status"]!="BOTH_WHOLE_PASS"
        assert row["tier"] in ("original","corrected","reconstructed")
        if row["both_status"]=="FINITE_ONLY_NOT_WHOLE":
            assert row["css"]["status"]=="PASS"
            assert row["maxwell"]["status"].startswith("PASS")
            assert row["crust_source_id"].split(":")[0] in ("native_reconstructed","public_original")
            assert row["same_functional"] is False

def test_matching_scope_is_only_maxwell_and_css():
    import compare_hermite_as_interpolation_only_20261118 as m
    assert m.MATCHING_METHODS==("maxwell","css")
    assert m.INTERPOLATION_METHODS==("cubic","quintic")
    assert m.EXPECTED_VARIANTS==("cubic+maxwell","quintic+maxwell","cubic+css","quintic+css")
    assert all("cubic_F" not in v and "quintic_P" not in v for v in m.EXPECTED_VARIANTS)


def test_hermite_resample_preserves_source_knots_and_uses_no_core_provider_mixing():
    import compare_hermite_as_interpolation_only_20261118 as m
    from compose_eos import read_compose
    d=read_compose(ROOT/"data/maxwell_three_families_20261117/prepared_v2/native_generated/SPG_M1/compose")
    before=(d.nb.copy(),d.p.copy(),d.eps.copy())
    transformed,receipt=m.resample_eos(d,"cubic")
    assert receipt["original_rows"]==len(d.nb)
    assert receipt["inserted_rows"]>0
    assert np.all(np.isfinite(transformed.nb))
    assert np.all(np.diff(transformed.nb)>0)
    # Every source value is retained as an original knot without interpolated
    # pressure values replacing it.
    for n,p,e in zip(*before):
        i=np.where(np.isclose(transformed.nb,n,rtol=0,atol=1e-14))[0]
        assert len(i)==1
        assert abs(transformed.p[i[0]]-p)<1e-9
        assert abs(transformed.eps[i[0]]-e)<1e-8
    np.testing.assert_array_equal(d.nb,before[0])
    np.testing.assert_array_equal(d.p,before[1])
    np.testing.assert_array_equal(d.eps,before[2])

def test_native_maxwell_on_same_cubic_quintic_interpolated_inputs():
    from pathlib import Path
    root=ROOT/"data/all_origin_maxwell_css_20261118/HERMITED_INPUT_2X2"
    receipt=json.loads((root/"NATIVE_MAXWELL/CAMPAIGN.json").read_text())
    assert receipt["matching"]=="maxwell"
    assert receipt["interpolation_orders"]==["cubic","quintic"]
    assert receipt["jobs"]==14
    assert receipt["native_maxwell_accepted"]==14
    assert all(x["status"]=="PASS_NATIVE_MAXWELL_MATCHED" for x in receipt["results"])
    for x in receipt["results"]:
        p=root/"NATIVE_MAXWELL"/x["tier"]/(x["crust"]+"__"+x["core"])/x["order"]/"RESULT.json"
        z=json.loads(p.read_text())
        assert z["status"]=="PASS_NATIVE_MAXWELL_MATCHED"
        assert 0<z["n_low"]<z["n_high"]
        assert z["transition_P"]>0 and z["latent_energy"]>0
        assert z["fortran_maxwell_identity"]<1e-8
        assert z["serialized_readback_abs"]<1e-9
