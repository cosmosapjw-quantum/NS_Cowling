"""Missing support and finite construction restrictions are not model no-go results."""
import pytest
from scripts.match_report_admission import classify_failure

@pytest.mark.parametrize("message,native_category,expected", [
    ("the two phases have no overlapping baryon-chemical-potential interval", "NO_COMMON_MU_SUPPORT", "input_contract"),
    ("", "NO_COMMON_MU_SUPPORT", "input_contract"),
    ("matching source support forbids cross-phase interpolation", None, "input_contract"),
    ("no convex matching window", None, "implementation_scope"),
    ("unsupported canonical profile: exact polynomial failed interval construction", None, "implementation_scope"),
])
def test_source_support_is_distinct_from_physical_impossibility(message, native_category, expected):
    assert classify_failure("matching", message, native_category) == expected
