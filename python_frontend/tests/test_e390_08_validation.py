"""E390-08 validation framework regression tests."""
import sys, json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))


def _needs_backend():
    import os
    return os.environ.get("COMPOSE_EOS_FORTRAN_EXE") is None


def _sample_batch_exists():
    return (ROOT / "download" / "m7_batch_v2_20261114" / "M7_BATCH_MANIFEST_v2.json").is_file()


@pytest.mark.skipif(_needs_backend() or not _sample_batch_exists(),
                    reason="requires the native backend + the sample batch")
def test_e390_08_validation_produces_all_four_receipts():
    from compose_eos.e390_m8_validation import (
        c1_c11_report, differential_comparison, edge_case_suite)
    batch = ROOT / "download" / "m7_batch_v2_20261114"
    c = c1_c11_report(batch)
    assert c["total_tasks"] == 390
    assert "summary" in c and "PASS" in c["summary"]
    d = differential_comparison(batch)
    assert d["n_pairs"] >= 1
    # each pair must have method entries with transition observables
    for p in d["pairs"]:
        assert len(p["methods"]) >= 1
        for m, info in p["methods"].items():
            assert "P_t_mean" in info and "mu_t_mean" in info
    e = edge_case_suite(batch)
    assert "category_counts" in e
    # the sample batch has composition-absent-core tasks (APR_APR)
    assert e["category_counts"].get("composition_absent_core", 0) > 0


@pytest.mark.skipif(_needs_backend() or not _sample_batch_exists(),
                    reason="requires the native backend + the sample batch")
def test_e390_08_differential_methods_agree_on_transition_pressure():
    """The 6 smooth methods must produce the same P_t for a given pair
    (structural consistency; the method differential is in the residual)."""
    from compose_eos.e390_m8_validation import differential_comparison
    batch = ROOT / "download" / "m7_batch_v2_20261114"
    d = differential_comparison(batch)
    for p in d["pairs"]:
        pts = [info["P_t_mean"] for info in p["methods"].values()
               if info.get("P_t_mean") is not None]
        if len(pts) >= 2:
            spread = max(pts) - min(pts)
            # methods must agree on P_t to within 1e-6 (structural consistency)
            assert spread < 1e-6, (
                f"methods disagree on P_t for {p['crust']}/{p['core']}: spread {spread}")


@pytest.mark.skipif(_needs_backend(), reason="requires the native backend")
def test_e390_08_edge_case_identifies_composition_absent_core():
    """APR_APR and BBB_BHF2 cores have no eos.compo -> composition-absent."""
    from compose_eos.eos_data import EOSData
    from compose_eos.api import read_compose
    for core in ("APR_APR", "BBB_BHF2"):
        d = read_compose(ROOT / "eos_generator" / "eos_data" / core)
        # composition-absent = no eos.compo sidecar
        assert not d.composition_file_present, f"{core} unexpectedly has eos.compo"
