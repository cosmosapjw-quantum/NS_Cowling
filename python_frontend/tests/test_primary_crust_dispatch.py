"""Primary-literature crust models use the production native generator."""
from pathlib import Path
import pytest
from compose_eos.generation_dispatch import generation_request, generator_catalog, _command

ROOT = Path(__file__).resolve().parents[2]
MODELS = ("BPS1971_sequence", "BBP1971_table5", "BBP1971_pressure_integral", "NV1974_author_fit")

@pytest.mark.parametrize("model", MODELS)
def test_primary_models_have_explicit_native_contract(model, tmp_path):
    request = generation_request(model)
    assert request["family"] == "primary_crust"
    assert request["contract"]["source_equivalence"] != "FULL_ORIGINAL_MICROPHYSICS"
    assert request["primary_source_path"]
    command = _command(ROOT, tmp_path / "generated", request)
    assert Path(command[0]).name == "run-generator.sh"
    content = Path(command[1]).read_text()
    assert "family='primary_crust'" in content
    assert "primary_source_path=" in content
    assert "sigma0=" not in content

@pytest.mark.parametrize("options", [
    {"include_lower_crust": True}, {"n_min": 1e-4}, {"n_max": .08},
    {"densities": [.001,.002,.003]}, {"resolutions": [16]},
    {"variant": "modern_masses"}, {"points": 2}, {"points": True},
])
def test_primary_models_reject_unrelated_physical_options(options):
    with pytest.raises(ValueError):
        generation_request("NV1974_author_fit", **options)

def test_raw_primary_table_is_not_resampled():
    assert generation_request("BBP1971_table5")["points"] == 26
    with pytest.raises(ValueError):
        generation_request("BBP1971_table5", points=257)
