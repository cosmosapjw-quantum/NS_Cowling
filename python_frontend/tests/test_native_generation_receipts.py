"""Independent transport verification of native generated tables.

Fixtures are deliberately artificial tables, not physical EOS evidence. They
exercise the real Fortran reader and diagnostics instead of trusting receipts.
"""
import json
import shutil
from pathlib import Path
import pytest
from compose_eos.generation_dispatch import generation_request
from compose_eos.native_generation import verify_native_generation, NATIVE_FILES

def write_table(path, densities):
    path.mkdir(parents=True, exist_ok=True)
    (path/"eos.nb").write_text("1\n"+str(len(densities))+"\n"+"\n".join(format(x,".17g") for x in densities)+"\n")
    (path/"eos.t").write_text("1\n1\n0\n")
    (path/"eos.yq").write_text("1\n1\n0\n")
    (path/"eos.thermo").write_text("939 939 1\n"+"".join(
        f"1 {i} 1 .1 0 .01 -.001 0 .02 .02 0\n" for i in range(1,len(densities)+1)))
    (path/"eos.compo").write_text("".join(
        f"1 {i} 1 6 4 10 .95 11 .05 0 .05 1 0 1 1 1 0 0\n" for i in range(1,len(densities)+1)))

def evidence(path, request, densities=None):
    densities = list(densities or request.get("densities") or [.1,.2,.4])
    write_table(path,densities)
    record={"backend":"fortran","model":request["model"],"family":request["family"],
            "status":"REQUESTED_ROWS_NATIVE_VERIFIED","native_readback":True,
            "rows":len(densities),"meshes":request.get("resolutions",[]),"held_nodes":0,
            "lower_crust_segments":bool(request.get("include_lower_crust"))}
    (path/"RESULT.json").write_text(json.dumps(record))
    for mesh in request.get("resolutions",[]):
        write_table(path/f"mesh_{mesh}",densities)
    if request.get("include_lower_crust"):
        write_table(path/"outer",[1e-5,1e-4,2.69432e-4])
        write_table(path/"lower_source_fit",[2.69403e-4,.001,.002])
        (path/"crust_segments.json").write_text(json.dumps({"cross_segment_interpolation":"UNAPPROVED"}))
    return record

def core_request(model="GM_GM1"):
    return generation_request(model,points=3,n_min=.1,n_max=.4)

def tf_request(**kwargs):
    return generation_request("GPPVA_DD2",densities=[.03,.04,.05],resolutions=[16,24],**kwargs)

def test_actual_core_readback_and_native_diagnostics(tmp_path):
    request=core_request();evidence(tmp_path,request)
    result=verify_native_generation(tmp_path,request)
    assert result["native_frontend_readback"]=={"generated":3}
    assert set(result["output_sha256"])==set(NATIVE_FILES)
    assert result["native_consistency"]["generated"]
    # Artificial source physical failures are diagnostics, not a receipt gate.
    assert result["rows"]==3

def test_fake_receipt_cannot_bless_malformed_tables(tmp_path):
    request=core_request();evidence(tmp_path,request)
    for name in NATIVE_FILES:(tmp_path/name).write_text("not a native table\n")
    with pytest.raises(ValueError,match="readback failed"):
        verify_native_generation(tmp_path,request)

def test_same_count_wrong_core_density_rejected(tmp_path):
    request=core_request();evidence(tmp_path,request,[.1,.25,.4])
    with pytest.raises(ValueError,match="density grid row missing"):
        verify_native_generation(tmp_path,request)

def test_apr_allows_inserted_in_window_rows(tmp_path):
    request=core_request("APR_APR");evidence(tmp_path,request,[.1,.15,.2,.3,.4])
    assert verify_native_generation(tmp_path,request)["rows"]==5

def test_apr_cannot_hide_missing_requested_row(tmp_path):
    request=core_request("APR_APR");evidence(tmp_path,request,[.1,.15,.25,.3,.4])
    with pytest.raises(ValueError,match="density grid row missing"):
        verify_native_generation(tmp_path,request)

def test_family_identity_is_not_inferred_from_model(tmp_path):
    request=core_request();record=evidence(tmp_path,request)
    record["family"]="tf";(tmp_path/"RESULT.json").write_text(json.dumps(record))
    with pytest.raises(ValueError,match="identity"):
        verify_native_generation(tmp_path,request)

def test_tf_reads_every_requested_mesh_and_hashes_them(tmp_path):
    request=tf_request();evidence(tmp_path,request)
    result=verify_native_generation(tmp_path,request)
    assert result["native_frontend_readback"]=={"generated":3,"mesh_16":3,"mesh_24":3}
    assert len(result["output_sha256"])==15

def test_tf_same_count_wrong_mesh_density_rejected(tmp_path):
    request=tf_request();evidence(tmp_path,request)
    write_table(tmp_path/"mesh_16",[.03,.045,.05])
    with pytest.raises(ValueError,match="mesh_16 density grid"):
        verify_native_generation(tmp_path,request)

def test_tf_missing_mesh_not_rescued_by_receipt(tmp_path):
    request=tf_request();evidence(tmp_path,request)
    shutil.rmtree(tmp_path/"mesh_16")
    with pytest.raises(ValueError,match="output missing"):
        verify_native_generation(tmp_path,request)

def test_main_tf_output_must_equal_finest_mesh(tmp_path):
    request=tf_request();evidence(tmp_path,request)
    p=tmp_path/"mesh_24/eos.thermo";p.write_text(p.read_text().replace(".1 0 .01",".2 0 .01"))
    with pytest.raises(ValueError,match="highest requested mesh"):
        verify_native_generation(tmp_path,request)

def test_lower_segments_are_read_back_without_forcing_a_seam(tmp_path):
    request=tf_request(include_lower_crust=True);evidence(tmp_path,request)
    result=verify_native_generation(tmp_path,request)
    assert result["native_frontend_readback"]["outer"]==3
    assert result["native_frontend_readback"]["lower_source_fit"]==3
    assert "outer/eos.thermo" in result["output_sha256"]

def test_lower_malformed_composition_cannot_pass(tmp_path):
    request=tf_request(include_lower_crust=True);evidence(tmp_path,request)
    (tmp_path/"lower_source_fit/eos.compo").write_text("invalid\n")
    with pytest.raises(ValueError,match="readback failed"):
        verify_native_generation(tmp_path,request)


def test_native_physical_diagnostics_are_hashed_without_becoming_pass_gates(tmp_path):
    request=core_request();record=evidence(tmp_path,request)
    (tmp_path/"sound_speed.csv").write_text("n_b,cs2_over_c2_fd\n.4,1.2\n")
    result=verify_native_generation(tmp_path,request)
    assert "sound_speed.csv" in result["output_sha256"]
    record["output_sha256"]=result["output_sha256"]
    (tmp_path/"RESULT.json").write_text(json.dumps(record))
    (tmp_path/"sound_speed.csv").write_text("modified evidence")
    with pytest.raises(ValueError,match="checksum mismatch"):
        verify_native_generation(tmp_path,request)


def test_native_attachment_evidence_is_hashed(tmp_path):
    request=tf_request(include_lower_crust=True);evidence(tmp_path,request)
    (tmp_path/"attachment_diagnostics.json").write_text(json.dumps({
        "equilibrium_established":False,"cross_interval_supported":False}))
    result=verify_native_generation(tmp_path,request)
    assert "attachment_diagnostics.json" in result["output_sha256"]


@pytest.mark.parametrize("name", ["sound_speed.csv", "attachment_diagnostics.json"])
def test_sealed_native_diagnostic_cannot_disappear(tmp_path, name):
    request=core_request();record=evidence(tmp_path,request)
    (tmp_path/name).write_text("native evidence")
    record["output_sha256"]=verify_native_generation(tmp_path,request)["output_sha256"]
    (tmp_path/"RESULT.json").write_text(json.dumps(record))
    (tmp_path/name).unlink()
    with pytest.raises(ValueError,match="sealed native diagnostic missing"):
        verify_native_generation(tmp_path,request)
