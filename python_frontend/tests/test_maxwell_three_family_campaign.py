"""Family/source-domain contracts of the production 195-cell campaign."""
from pathlib import Path
import sys
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT / "python_frontend")]
import maxwell_three_family_campaign as campaign
from compose_eos.eos_data import EOSData


def table():
    n = np.array([.001, .002, .003, .004, .005, .006])
    z = np.zeros(len(n))
    return EOSData(nb=n, mn=939., mp=938., Q1=n.copy(), Q2=z.copy(), Q3=n.copy(),
                   Q4=z.copy(), Q5=z.copy(), Q6=n.copy(), Q7=n.copy(),
                   q_extra=np.array([1.,1.,1.,2.,2.,3.]), Iphase=np.full(len(n),7),
                   Yn=z.copy(), Yp=z.copy(), Ye=z.copy(), Ymu=z.copy(), Aav=z.copy(),
                   Zav=z.copy(), Yh=z.copy(), q_extra_count=1, beta_equilibrated=True)


def intervals(valid):
    return [{"left_row_zero_based": i, "right_row_zero_based": i+1, "valid": v}
            for i,v in enumerate(valid)]


def test_family_guard_rejects_every_cross_family_direction():
    for low in campaign.FAMILIES:
        for high in campaign.FAMILIES:
            a,b = {"family":low,"role":"crust"},{"family":high,"role":"core"}
            if low == high:
                campaign.check_same_family(a,b)
            else:
                with pytest.raises(ValueError,match="cross-family"):
                    campaign.check_same_family(a,b)


def test_components_never_bridge_bad_intervals_or_admit_two_knots():
    assert campaign.valid_components(intervals([True,True,False,True,True,False,True]),8) == [(0,3),(3,6)]
    assert campaign.valid_components(intervals([False,True]),3) == []
    assert campaign.valid_components(intervals([True,True]),3) == [(0,3)]
    with pytest.raises(ValueError,match="cover every"):
        campaign.valid_components(intervals([True]),3)
    wrong=intervals([True,True]);wrong[0]["right_row_zero_based"]=2
    with pytest.raises(ValueError,match="indexing"):
        campaign.valid_components(wrong,3)


def test_public_core_boundary_uses_provider_cut_even_if_phase_labels_are_identical():
    data=table();data.Iphase[:]=3
    entry={"homogeneous_core_selection":{"start_zero_based":3,"stop_exclusive":6,"evidence":"provider"}}
    selected,receipt=campaign.role_slice(data,entry,"core","public_original")
    assert np.array_equal(selected.nb,data.nb[3:])
    assert receipt["excluded_lower_rows"] == 3
    assert np.array_equal(data.nb,table().nb)


def test_documented_auxiliary_geometry_is_preserved_without_thermo_repair(tmp_path,monkeypatch):
    data=table();before=data.copy()
    monkeypatch.setattr(campaign,"ROOT",tmp_path)
    entry={"model":"GPPVA_DD2","family":"public_original","path":"raw/DD2",
           "table_sha256":{},"extra_column_primary_evidence":{"page":2},
           "extra_column_semantics":"provider geometry"}
    out,receipt=campaign.auxiliary_projection(data,entry,{"source_start_zero_based":217},tmp_path/"projection")
    for name in ("nb","Q1","Q2","Q3","Q4","Q5","Q6","Q7"):
        assert np.array_equal(getattr(out,name),getattr(before,name))
    assert out.q_extra_count == 0
    assert np.array_equal(out.Iphase,[1,1,1,2,2,3])
    assert np.array_equal(data.q_extra,before.q_extra)
    assert data.q_extra_count == 1
    payload=campaign.load(tmp_path/receipt["payload_path"])
    assert payload["rows"][0]["source_row_zero_based"] == 217
    assert [r["payload"][0] for r in payload["rows"]] == before.q_extra.tolist()


def test_unknown_auxiliary_not_given_invented_physical_meaning(tmp_path,monkeypatch):
    data=table();data.q_extra=np.linspace(-.01,1.2,data.n)
    monkeypatch.setattr(campaign,"ROOT",tmp_path)
    entry={"model":"RG_KDE0v1","family":"public_original","path":"raw/KDE","table_sha256":{}}
    out,receipt=campaign.auxiliary_projection(data,entry,{"source_start_zero_based":0},tmp_path/"projection")
    payload=campaign.load(tmp_path/receipt["payload_path"])
    assert payload["semantics"] == "UNKNOWN_PROVIDER_AUXILIARY"
    assert np.array_equal(out.Iphase,data.Iphase)
    assert payload["rows"][-1]["payload"] == [1.2]


def test_kernel_attempt_means_native_command_was_invoked(monkeypatch):
    state={"kernel_invoked":False}
    def rejected(action,*args):
        raise RuntimeError("a concrete native command failed")
    monkeypatch.setattr(campaign.eos_api,"run_native",rejected)
    with pytest.raises(ValueError):
        with campaign.observe_kernel(state):
            raise ValueError("preflight failed before native command")
    assert not state["kernel_invoked"]
    with pytest.raises(RuntimeError):
        with campaign.observe_kernel(state):
            campaign.eos_api.run_native("maxwell-transition","low","high","out","true")
    assert state["kernel_invoked"]
    assert campaign.eos_api.run_native is rejected
