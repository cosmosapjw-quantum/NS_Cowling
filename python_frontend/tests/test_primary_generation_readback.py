"""Real public native generation, optional abundances and diagnostic seals."""
import json
import os
from pathlib import Path
import pytest
from compose_eos.generation_dispatch import generation_request, run_generation
from compose_eos.native_generation import verify_native_generation
from compose_eos.primary_model_contracts import PRIMARY_IDS

@pytest.mark.parametrize("model", PRIMARY_IDS)
def test_public_primary_generation_readback_and_seal(model, tmp_path, monkeypatch):
    # CI supplies the compiled executable; normal source-checkout execution
    # uses the same source-signature-checked build as other native tests.
    from compose_eos.native import executable
    exe=executable()
    monkeypatch.setenv("COMPOSE_EOS_NATIVE_BUILD_DIR",str(exe.parent))
    monkeypatch.setenv("COMPOSE_EOS_FORTRAN_EXE",str(exe))
    output=tmp_path/"run"
    assert run_generation(model,output)==0
    request=generation_request(model)
    evidence=verify_native_generation(output/"generated",request)
    assert evidence["rows"]>=26
    compo=output/"generated/eos.compo"
    assert compo.exists()==(model=="BPS1971_sequence")
    for name in request["contract"]["diagnostics"]:
        assert name in evidence["output_sha256"]
    # A source-consistency failure is preserved, not turned into a generation
    # failure or a fabricated microphysical-equivalence receipt.
    if model=="BBP1971_table5":
        diagnostics=json.loads((output/"generated/primary_source_diagnostics.json").read_text())
        assert diagnostics
        assert evidence["rows"]==26
    diagnostic=output/"generated"/request["contract"]["diagnostics"][0]
    diagnostic.write_bytes(diagnostic.read_bytes()+b"\n")
    with pytest.raises(ValueError,match="checksum mismatch"):
        verify_native_generation(output/"generated",request)

def test_bps_small_grid_fails_before_output():
    with pytest.raises(ValueError):
        generation_request("BPS1971_sequence",points=15)
