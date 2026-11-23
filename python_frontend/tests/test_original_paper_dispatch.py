"""Original-paper models have explicit native provenance and fixed source domains."""
import pytest
from compose_eos.generation_dispatch import generation_request, generator_catalog
from compose_eos.primary_model_contracts import PRIMARY_IDS

MODELS = (
    "BBP1971_original_published", "BBP1971_original_variational", "BBP1971_variational_maxwell",
    "NV1973_author_fit", "BPS_Myers1966_magic14", "BPS_Myers1966_magic20",
)

@pytest.mark.parametrize("model", MODELS)
def test_original_paper_models_have_native_source_contract(model):
    request = generation_request(model)
    assert request["family"] == "primary_crust"
    assert request["points"] >= 17
    contract = request["contract"]
    assert contract["source_path"].startswith("data/eos_original_papers_20261123/")
    assert contract["doi"]
    assert len(contract["source_document_sha256"]) == 64
    assert contract["source_equivalence"]
    assert contract["diagnostics"]

@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("options", [{"include_lower_crust": True}, {"n_min": 0.001}, {"points": True}])
def test_original_models_reject_changes_to_fixed_source_contract(model, options):
    with pytest.raises(ValueError):
        generation_request(model, **options)

def test_historical_primary_campaign_model_list_is_preserved():
    assert PRIMARY_IDS == ("BPS1971_sequence", "BBP1971_table5", "BBP1971_pressure_integral", "NV1974_author_fit")

def test_original_nv_fit_declares_same_physical_function_as_1974():
    contract = generator_catalog()["NV1973_author_fit"]
    assert contract["semantic_equivalent_model"] == "NV1974_author_fit"
    assert contract["full_hartree_fock_rerun"] is False

def test_bbp_original_closures_are_explicitly_distinct():
    catalog = generator_catalog()
    left = catalog["BBP1971_original_published"]
    right = catalog["BBP1971_original_variational"]
    assert left["source_path"] == right["source_path"]
    assert left["derivative_closure"] != right["derivative_closure"]

def test_myers_magic_choices_do_not_change_source_constants():
    catalog = generator_catalog()
    left = catalog["BPS_Myers1966_magic14"]
    right = catalog["BPS_Myers1966_magic20"]
    assert left["source_path"] == right["source_path"]
    assert (left["magic_choice"], right["magic_choice"]) == (14,20)
