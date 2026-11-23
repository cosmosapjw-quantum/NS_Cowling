"""Public all-method routing and complete native search receipts."""
from pathlib import Path
from compose_eos import read_compose, match_crust_core

ROOT = Path(__file__).resolve().parents[2]
METHOD = "first_law_constrained_hermite_pressure_bridge"


def test_public_first_law_bridge_has_complete_method_identity_and_search_receipt():
    crust = read_compose(ROOT / "eos_generator/eos_data_normalized/GPPVA_FSU2H")
    core = read_compose(ROOT / "eos_generator/eos_data_normalized/GM_GM1")
    matched, report = match_crust_core(
        crust, core, method=METHOD, n_bridge=40, drop_composition=True)
    assert matched.n > 40
    assert report["method"] == report["native_report"]["method"] == METHOD
    search = report["window_search"]
    assert search["pairs_total"] >= search["pairs_examined"] >= 1
    assert search["candidates_constructed"] >= 1
    assert search["candidates_rejected"] >= 0
    assert search["search_exhausted"] is False
    assert report["native_report"]["first_law_relative_residual"] <= 5e-9
    assert report["consistency"]["overall"]["pass_"]
