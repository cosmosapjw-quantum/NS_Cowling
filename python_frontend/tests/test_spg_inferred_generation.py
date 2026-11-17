"""Public path and provenance regression for the versioned SPG candidate."""
import json
from pathlib import Path
import pytest
from compose_eos.generation_dispatch import generation_request, generator_catalog
from compose_eos.native_generation import surface_input, write_request
from compose_eos.native_model_contracts import SPG_INFERRED_VARIANT, SPG_INFERRED_SURFACE
ROOT=Path(__file__).resolve().parents[2]

@pytest.mark.parametrize("model", [f"SPG_M{i}" for i in range(1,6)])
def test_explicit_inferred_request(model,tmp_path):
    request=generation_request(model,variant=SPG_INFERRED_VARIANT)
    assert request["points"]==129
    assert generation_request(model)["variant"]=="AME2408_candidate"
    surface=surface_input(ROOT,model,SPG_INFERRED_VARIANT)
    source=json.loads((ROOT/SPG_INFERRED_SURFACE["path"]).read_text())
    assert surface["sigma0"]==source["parameters"][model]["sigma0"]
    assert surface["surface_b"]==source["parameters"][model]["b"]
    assert surface["original_mass_calibration"]=="UNRESOLVED"
    assert "NOT_ORIGINAL_AUTHOR_CALIBRATION" in surface["source_equivalence"]
    text=write_request(ROOT,tmp_path/"generated",request).read_text()
    assert "fixed_density_source" in text
    assert format(surface["sigma0"],".17e") in text
    assert SPG_INFERRED_VARIANT in generator_catalog()[model]["variants"]

def test_inferred_input_is_immutable(tmp_path):
    path=tmp_path/SPG_INFERRED_SURFACE["path"]
    path.parent.mkdir(parents=True)
    path.write_bytes((ROOT/SPG_INFERRED_SURFACE["path"]).read_bytes()+b" ")
    with pytest.raises(ValueError,match="checksum mismatch"):
        surface_input(tmp_path,"SPG_M1",SPG_INFERRED_VARIANT)
