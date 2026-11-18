"""Scientific evidence and source-chain guards for COWLING-MAXWELL-NATIVE-03."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from compose_eos import read_compose
from run_bbp_bps_matched_core_20261119 import read_bbp_secondary,make_bridge
from compose_eos.whole_crust import WholeCrustBarotrope

BASE=ROOT/"data/cowling_maxwell_native_20261119"
def load(path):return json.loads(path.read_text())
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def test_same_rmf_original_fsu2_native_two_matching_methods():
    d=load(BASE/"same_rmf/SURVEY.json")
    f=next(x for x in d["cases"] if x.get("label")=="GPPVA_FSU2"
           and x.get("tier")=="original")
    assert f["dual_pass_count"]==1
    x=next(row for row in f["matches"] if row.get("maxwell",{}).get("status")=="PASS_NATIVE_FINITE")
    assert x["css"]["status"]=="PASS_NATIVE_FINITE"
    pt=x["maxwell"]["P_t"]
    for interp in ("cubic","quintic"):
        z=x["interpolation"][interp]
        assert z["status"]=="POSITIVE_JUMP_ROOT"
        assert z["roots"][0]["n_high"]>z["roots"][0]["n_low"]
        assert abs(z["roots"][0]["P_MeVfm3"]-pt)<2e-5
    assert f["original_source_Q3_minus_Euler_mu_max_abs_MeV"]<0.002


@pytest.mark.parametrize("model",[ "FSU2","TM1e","DDME2"])
def test_official_full_eos_css_source_projection_and_modes(model):
    d=load(BASE/"provider_original_unified/ORIGINAL_CONTROLS.json")
    a=d["models"][model]
    assert a["independent_thermo_matches_zip_member"] is True
    assert a["provider_zip_checksum_matches"] is False
    assert a["first_law_corrected_split_methods"]["css"]["status"]=="PASS_NATIVE_FINITE"
    assert a["raw_original_split_methods"]["maxwell"]["status"]=="NATIVE_REJECTED"
    if model=="DDME2":
        projection=a["source_projection"]
        assert projection["retained_source_rows"]==413
        assert projection["original_rows"]==415
        assert projection["rejected_tail_source_rows"]==[414,415]
        prefix=read_compose(ROOT/projection["projection_path"])
        assert prefix.n==413
        assert np.all(np.diff(prefix.nb)>0)
        # Source table is *never* replaced in place by the trimmed projection.
        lines=(BASE/"provider_original_unified/DDME2/raw/eos.nb").read_text().splitlines()
        assert int(lines[1])==415
    sharp=load(BASE/"official_unified_css_sharp/RESULT.json")["cases"][model]
    assert sharp["status"]=="CSS_THREE_LAYER_NATIVE_STELLAR_CHECKED"
    assert len(sharp["native_fortran"]["native_tov"]["interfaces"])==1
    assert abs(sharp["native_fortran"]["mass_delta_msun"])<1e-5
    assert abs(sharp["native_fortran"]["radius_delta_km"])<1e-4
    for conv in ("slow","rapid"):
        c=sharp["native_fortran"]["modes"][conv]
        assert c["status"]=="NATIVE_MODES_VERIFIED"
        assert {m["label"] for m in c["modes"]}=={"f","p1"}
        assert all(m["join_relative"]<1e-5 and m["interface_max_traction"]<1e-9
                   for m in c["modes"])


def test_bbp_table_origin_and_low_density_bridge_precision():
    n,p,e,receipt=read_bbp_secondary()
    assert len(n)==35 and receipt["source_kind"].startswith("SECONDARY_")
    assert receipt["original_primary_BBP_bytes"]=="NOT_IN_HAND"
    d=load(BASE/"bbp_bps/RESULT.json")
    assert d["source"]["source_pdf_sha256"]==digest(
        BASE/"bbp_bps/source_reference/diener_2008_thesis.pdf")
    bps=read_compose(ROOT/"data/maxwell_three_families_20261117/intake/references/APR_APR")
    joined,bridge=make_bridge(n[:26],p[:26],e[:26],bps)
    assert bridge["status"]=="ROUNDED_FIRST_N_WITHIN_PUBLISHED_PRECISION_AND_CAUSAL"
    assert abs(bridge["delta_n_fm3"])<5e-7
    assert joined[2].validate()["status"]=="PASS"
    assert d["source_reconstruction"]["max_abs_delta_energy_per_baryon_MeV"]>5


def test_bbp_bps_css_native_four_mode_scope():
    d=load(BASE/"bbp_bps/RESULT.json")
    rows={c["core"]:c for c in d["combinations"]
          if c["branch"]=="pressure_primary_first_law_BBP"}
    assert set(rows)=={"GPPVA_FSU2","GPPVA_TM1e","GPPVA_DDME2","GM_GM1","PK1_CANONICAL"}
    for name,c in rows.items():
        assert c["css"]["status"]=="PASS_NATIVE_FINITE"
        assert c["maxwell"]["status"]=="REJECTED"
        x=c["css_three_layer"]
        assert x["stellar"]["native_tov"]["mass_msun"]>0
        assert x["stellar"]["native_tov"]["R_km"]>0
        if name in {"GPPVA_FSU2","GPPVA_TM1e","PK1_CANONICAL"}:
            assert x["status"]=="THREE_LAYER_NATIVE_TOV_F_P1_BOTH_CONVERSIONS_VERIFIED"
            for conv in ("slow","rapid"):
                mode=x["stellar"]["modes"][conv]
                assert mode["status"]=="NATIVE_MODES_VERIFIED"
                assert {a["label"] for a in mode["modes"]}=={"f","p1"}
        else:
            assert x["status"]=="THREE_LAYER_NATIVE_TOV_PARTIAL_MODE_IDENTITY"
            assert x["stellar"]["modes"]["slow"]["status"]=="PARTIAL_MODE_SET"


def test_native_phasewise_controls_preserve_slow_rapid_difference():
    d=load(BASE/"three_layer_original_fsu2/RESULT.json")
    for name in ("maxwell","css"):
        x=d["methods"][name]
        assert x["status"]=="THREE_LAYER_NATIVE_TOV_COWLING_CHECKED"
        assert x["admitted"]["join_index_k"]==0
        w=WholeCrustBarotrope.from_contract(load(
           BASE/"three_layer_original_fsu2"/name/"REPRESENTATION.json"))
        assert len(w.sharp_interfaces)==1
        star=x["admitted"]["stellar"]
        assert abs(star["mass_delta_msun"])<1e-5
        for conv in ("slow","rapid"):
            assert star["modes"][conv]["status"]=="NATIVE_MODES_VERIFIED"
            assert all(t["join_relative"]<1e-5 for t in star["modes"][conv]["modes"])
        f_slow=star["modes"]["slow"]["modes"][0]["frequency_hz"]
        f_rapid=star["modes"]["rapid"]["modes"][0]["frequency_hz"]
        assert abs(f_slow-f_rapid)>0.001
