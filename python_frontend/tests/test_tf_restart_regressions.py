"""Actual failed cell profiles guard numerical restarts without weaker physics."""
import gzip,importlib.util,json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[2]

@pytest.mark.parametrize("case",[0,1])
def test_energy_floor_and_nearly_absent_species_counterexamples(case):
    path=ROOT/"data/research_completion_20261116/gppva_tf/numerical_counterexamples.json.gz"
    old=json.loads(gzip.decompress(path.read_bytes()))[case]
    assert old["numerical_status"]=="REQUIRES_REFINEMENT"
    if case==0:
        assert not old["optimizer_success"]
        assert old["max_chemical_KKT_residual_MeV"]<.02
    else:
        assert old["max_chemical_KKT_residual_MeV"]>.02
    spec=importlib.util.spec_from_file_location("tf_polish",ROOT/"scripts/tf_polish.py")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    final=module.polish_profile(old)
    assert final["numerical_status"]=="STATIONARY"
    assert final["max_chemical_KKT_residual_MeV"]<.02
    assert final["radius_stationarity_MeV"]<2e-4
    assert final["baryon_constraint_error"]<2e-10
    assert final["charge_constraint_error"]<2e-10
    assert final["epsilon"]<=old["epsilon"]+1e-10
    assert final["restart_attempts"]
