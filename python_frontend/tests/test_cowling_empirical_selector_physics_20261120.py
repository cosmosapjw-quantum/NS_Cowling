"""Regression: a fast run is not admissible without its physical mode evidence."""
from __future__ import annotations

import pytest
from scripts.summarize_hybrid_accuracy_20261120 import choose_verified_backend


def _candidate(seconds=1.0):
    return {
        "status": "MATCHED_REFERENCE_PASS",
        "backend": "hybrid_cs2",
        "samples": 2049,
        "frequency_hz": 3000.0005,
        "reference_delta_hz": 0.0005,
        "total_seconds": seconds,
        "join_relative": 1e-7,
        "nodes": [1, 1],
    }


def _case():
    return {
        "reference_hz": 3000.0,
        "backend_variants": {
            "FAST": _candidate(0.1),
            "VALID_P1": _candidate(1.0),
        },
    }


@pytest.mark.parametrize("mutations", [
    {"nodes": [0, 0]},
    {"join_relative": -1e-7},
    {"frequency_hz": 3001.0},
])
def test_selector_rejects_false_pass_physical_evidence(mutations):
    case = _case()
    case["backend_variants"]["FAST"].update(mutations)
    chosen = choose_verified_backend(case, 0.005)
    assert chosen["selected"] == "VALID_P1"


def test_selector_returns_no_admitted_candidate_if_only_wrong_mode_exists():
    case = _case()
    del case["backend_variants"]["VALID_P1"]
    case["backend_variants"]["FAST"]["nodes"] = [0, 0]
    chosen = choose_verified_backend(case, 0.005)
    assert chosen["selected"] is None
    assert chosen["status"] == "NO_VERIFIED_CANDIDATE_AT_REQUESTED_FREQUENCY_TOLERANCE"


def test_selector_still_admits_consistent_observed_p1():
    case = _case()
    assert choose_verified_backend(case, 0.005)["selected"] == "FAST"
