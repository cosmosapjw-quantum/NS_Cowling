"""Public generation contract checks; no target-table agreement gate."""
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from compose_eos import generation_dispatch as dispatch
from compose_eos.cli import main


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def native(path):
    path.mkdir(parents=True, exist_ok=True)
    for name in dispatch.NATIVE_FILES:
        (path / name).write_text("native-output-fixture\n")


def core_evidence(output, points=3):
    model = "GM_GM1"
    directory = output / model
    native(directory)
    rows = [{"n_b": n} for n in (.1, .2, .4)]
    dump(output / "SUMMARY.json", {"status": "COMPLETE", "models": {
        model: {"points": points, "native_readback": "PASS"}}})
    dump(directory / "RECEIPT.json", {"rows": rows, "native_readback": "PASS",
        "output_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in directory.glob("eos.*")}})
    return directory


def test_catalog_is_exact_and_scope_is_explicit():
    catalog = dispatch.generator_catalog()
    assert set(catalog) == set(dispatch.MODEL_IDS)
    assert len(catalog) == 22
    assert all(catalog[x]["region"] == "core_only" for x in dispatch.CORE_IDS)
    assert all(catalog[x]["region"] == "inner_crust_only" for x in dispatch.SPG_IDS)
    assert all("masks unavailable" in catalog[x]["composition"] for x in dispatch.TF_IDS)
    assert "kernel is not re-solved" in catalog["BBB_BHF2"]["paper_implementation"]


@pytest.mark.parametrize("model,options", [
    ("GM_GM1", {"points": 2}), ("GM_GM1", {"n_min": float("nan")}),
    ("GM_GM1", {"n_min": .4, "n_max": .1}),
    ("GM_GM1", {"densities": [.1,.2,.3]}),
    ("SPG_M1", {"points": 16}), ("SPG_M1", {"points": 25.0}),
    ("SPG_M1", {"points": 18}), ("SPG_M1", {"variant": "invented"}),
    ("SPG_M1", {"n_min": .001}),
    ("GPPVA_DD2", {}), ("GPPVA_DD2", {"densities":[.001,.002,.003], "resolutions":[0]}),
    ("GPPVA_DD2", {"densities":[.001,.003,.002], "resolutions":[16]}),
    ("GPPVA_DD2", {"densities":[.001,.002,.003], "resolutions":[32,16]}),
    ("GPPVA_DD2", {"densities":[.001,.002,.003], "resolutions":[16], "points":17}),
    ("GM_GM1", {"workers":0}), ("unknown", {}),
])
def test_invalid_requests_fail_before_creating_output(model, options):
    with pytest.raises(ValueError):
        dispatch.generation_request(model, **options)


def test_spg_candidate_defaults_and_explicit_tf_grid():
    assert dispatch.generation_request("SPG_M1")["points"] == 129
    assert dispatch.generation_request("SPG_M1")["variant"] == "AME2408_candidate"
    request = dispatch.generation_request("GPPVA_TW", densities=[.001,.002,.003],
                                          resolutions=[16,32])
    assert request["resolutions"] == [16,32]


@pytest.fixture
def checkout(tmp_path):
    root = tmp_path / "checkout"
    (root / ".git").mkdir(parents=True)
    for script in dispatch.SCRIPTS.values():
        p = root / "scripts" / script
        p.parent.mkdir(exist_ok=True)
        p.write_text("# source runner")
    (root / "python_frontend/compose_eos").mkdir(parents=True)
    (root / "python_frontend/compose_eos/api.py").write_text("# api")
    return root


def test_requires_source_checkout(tmp_path):
    with pytest.raises(ValueError, match="source checkout"):
        dispatch.repository_root(tmp_path)


def test_nonzero_child_exit_and_failure_evidence_survive(tmp_path, checkout, monkeypatch):
    out = tmp_path / "result"
    def child(command, **kwargs):
        dump(Path(command[2]) / "FAILURE.json", {"reason":"actual child failure"})
        return SimpleNamespace(returncode=7)
    monkeypatch.setattr(dispatch.subprocess, "run", child)
    assert dispatch.run_generation("GM_GM1", out, root=checkout, points=3) == 7
    assert (out / "generated/FAILURE.json").exists()
    receipt = json.loads((out / "DISPATCH.json").read_text())
    assert receipt["child_returncode"] == 7 and receipt["status"] == "CHILD_FAILED"


def test_zero_exit_without_outputs_is_failure(tmp_path, checkout, monkeypatch):
    monkeypatch.setattr(dispatch.subprocess, "run", lambda *a, **kw: SimpleNamespace(returncode=0))
    out = tmp_path / "result"
    assert dispatch.run_generation("GM_GM1", out, root=checkout, points=3) == 1
    assert json.loads((out / "DISPATCH.json").read_text())["status"] == "EVIDENCE_OR_EXECUTION_FAILED"


def test_existing_output_is_preserved(tmp_path, checkout):
    (tmp_path / "result").mkdir()
    sentinel = tmp_path / "result/sentinel"
    sentinel.write_text("keep")
    with pytest.raises(FileExistsError):
        dispatch.run_generation("GM_GM1", sentinel.parent, root=checkout, points=3)
    assert sentinel.read_text() == "keep"


def test_core_evidence_requires_requested_rows_and_native_hashes(tmp_path):
    request = dispatch.generation_request("GM_GM1", points=3, n_min=.1, n_max=.4)
    directory = core_evidence(tmp_path)
    assert dispatch.verify_generation(tmp_path, request)["rows"] == 3
    (directory / "eos.compo").write_text("corrupt")
    with pytest.raises(ValueError, match="checksum"):
        dispatch.verify_generation(tmp_path, request)


def test_core_missing_requested_grid_row_is_failure(tmp_path):
    request = dispatch.generation_request("GM_GM1", points=3, n_min=.1, n_max=.4)
    directory = core_evidence(tmp_path)
    record = json.loads((directory / "RECEIPT.json").read_text())
    record["rows"][1]["n_b"] = .25
    dump(directory / "RECEIPT.json", record)
    with pytest.raises(ValueError, match="missing"):
        dispatch.verify_generation(tmp_path, request)


def test_spg_hold_cannot_be_promoted_by_zero_exit(tmp_path):
    request = dispatch.generation_request("SPG_M1", points=17)
    dump(tmp_path / "SUMMARY.json", {"status":"PARTIAL_WITH_HOLD","cases":[]})
    dump(tmp_path / "SPG_M1/AME2408_candidate/RECEIPT.json", {"status":"HOLD"})
    with pytest.raises(ValueError):
        dispatch.verify_generation(tmp_path, request)


def tf_evidence(output):
    model = "GPPVA_DD2"
    densities = [.001,.002,.003]
    directory = output / model
    native(directory / "compose")
    dump(output / "CAMPAIGN.json", {"source_unchanged_during_execution":True,
        "models":[{"model":model,"rows":3}], "densities_fm3":densities, "resolutions":[16]})
    dump(directory / "RECEIPT.json", {"model_id":model,"stages":[{
        "points":16,"requested_points":3,"selected_stationary_points":3}],
        "export":{"status":"EXPORTED_NATIVE_READBACK_VERIFIED","row_count":3}})
    rows = [{"n_b":n,"model_id":model,"points":16,"numerical_status":"STATIONARY", "muon_threshold_crossed":False,
             "baryon_constraint_error":0.,"charge_constraint_error":0.,
             "scalar_pde_residual":0.,"vector_pde_residual":0.,"electron_pde_residual":0.,
             "epsilon":1.,"energy_per_baryon":939.,"pressure":.01,"mu_B":939.,"mu_Q_lepton":10.}
            for n in densities]
    dump(directory / "composition.json", {"row_count":3,"rows":rows})
    with gzip.open(directory / "selected_states.json.gz","wt") as stream:
        json.dump(rows,stream)
    seal_tf_export(directory)
    return directory


def test_tf_success_needs_actual_profile_density_rows(tmp_path):
    request = dispatch.generation_request("GPPVA_DD2", densities=[.001,.002,.003],resolutions=[16])
    directory = tf_evidence(tmp_path)
    assert dispatch.verify_generation(tmp_path, request)["rows"] == 3
    with gzip.open(directory / "selected_states.json.gz","wt") as stream:
        json.dump([{"n_b":.001},{"n_b":.003}],stream)
    seal_tf_export(directory)
    with pytest.raises(ValueError, match="profile/composition"):
        dispatch.verify_generation(tmp_path, request)


def test_tf_partial_segments_are_not_complete(tmp_path):
    request = dispatch.generation_request("GPPVA_DD2", densities=[.001,.002,.003],resolutions=[16])
    directory = tf_evidence(tmp_path)
    record=json.loads((directory/"RECEIPT.json").read_text())
    record["export"]["status"]="PARTIAL_SEGMENTS_NATIVE_READBACK_VERIFIED"
    dump(directory/"RECEIPT.json",record)
    with pytest.raises(ValueError,match="all-row"):
        dispatch.verify_generation(tmp_path,request)


def test_cli_catalog_and_dispatch(tmp_path,monkeypatch,capsys):
    assert main(["generators"]) == 0
    assert len(json.loads(capsys.readouterr().out)) == 22
    calls=[]
    monkeypatch.setattr(dispatch,"run_generation",lambda *a,**k: calls.append((a,k)) or 7)
    assert main(["generate","SPG_M1",str(tmp_path/"new"),"--points","17"]) == 7
    assert calls[0][0][0] == "SPG_M1"
    assert calls[0][1]["points"] == 17


def seal_tf_export(directory):
    required = ["compose/" + name for name in dispatch.NATIVE_FILES]
    required += ["composition.json", "selected_states.json.gz"]
    if (directory / "INTERVAL_SUPPORT.json").exists():
        required.append("INTERVAL_SUPPORT.json")
    record = json.loads((directory / "RECEIPT.json").read_text())
    record["export"]["output_sha256"] = {
        name: hashlib.sha256((directory/name).read_bytes()).hexdigest() for name in required}
    dump(directory / "RECEIPT.json", record)


def test_tf_native_corruption_is_rejected(tmp_path):
    request = dispatch.generation_request("GPPVA_DD2", densities=[.001,.002,.003],resolutions=[16])
    directory = tf_evidence(tmp_path)
    (directory / "compose/eos.thermo").write_text("corrupt table")
    with pytest.raises(ValueError,match="checksum"):
        dispatch.verify_generation(tmp_path,request)


@pytest.mark.parametrize("omit", ["compose/eos.thermo", "composition.json", "selected_states.json.gz", "INTERVAL_SUPPORT.json"])
def test_tf_required_hash_coverage(tmp_path, omit):
    request = dispatch.generation_request("GPPVA_DD2", densities=[.001,.002,.003],resolutions=[16])
    directory = tf_evidence(tmp_path)
    dump(directory / "INTERVAL_SUPPORT.json", {"interpolation":"NOT_ADMITTED"})
    seal_tf_export(directory)
    record = json.loads((directory/"RECEIPT.json").read_text())
    del record["export"]["output_sha256"][omit]
    dump(directory/"RECEIPT.json",record)
    with pytest.raises(ValueError,match="checksums"):
        dispatch.verify_generation(tmp_path,request)


@pytest.mark.parametrize("key,value", [
    ("numerical_status","REQUIRES_REFINEMENT"),
    ("model_id","GPPVA_TW"),
    ("points",32),
    ("muon_threshold_crossed",True),
    ("baryon_constraint_error",1e-3),
    ("charge_constraint_error",1e-3),
    ("scalar_pde_residual",1e-3),
    ("electron_pde_residual",1e-3),
    ("epsilon",float("nan")),
    ("energy_per_baryon",float("inf")),
    ("pressure",float("nan")),
    ("mu_B",float("inf")),
    ("mu_Q_lepton",float("nan")),
])
def test_tf_selected_row_requires_physical_evidence_even_after_rehash(tmp_path,key,value):
    request = dispatch.generation_request("GPPVA_DD2", densities=[.001,.002,.003],resolutions=[16])
    directory = tf_evidence(tmp_path)
    path = directory / "selected_states.json.gz"
    with gzip.open(path,"rt") as stream:
        rows=json.load(stream)
    rows[0][key]=value
    with gzip.open(path,"wt") as stream:
        json.dump(rows,stream)
    seal_tf_export(directory)
    with pytest.raises(ValueError,match="selected TF row"):
        dispatch.verify_generation(tmp_path,request)


def test_tf_all_requested_stage_meshes_must_match(tmp_path):
    request = dispatch.generation_request("GPPVA_DD2", densities=[.001,.002,.003],resolutions=[8,16])
    directory = tf_evidence(tmp_path)
    campaign=json.loads((tmp_path/"CAMPAIGN.json").read_text())
    campaign["resolutions"]=[8,16]
    dump(tmp_path/"CAMPAIGN.json",campaign)
    record=json.loads((directory/"RECEIPT.json").read_text())
    record["stages"].insert(0,{"points":12,"requested_points":3,"selected_stationary_points":3})
    dump(directory/"RECEIPT.json",record)
    with pytest.raises(ValueError,match="all-row"):
        dispatch.verify_generation(tmp_path,request)
