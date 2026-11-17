"""Native topology controls built from exact cubic P(mu), MeV/fm^3."""
import numpy as np
import pytest
from matching_theory_soft_barotrope import soft_branch
from compose_eos import NativeExecutionError, construct_maxwell_transition


def phase(mu, perturbation="none"):
    x = mu - 990.0
    p = 8.0 + 0.12 * x + 0.001 * x**2
    n = 0.12 + 0.002 * x
    if perturbation == "three":
        p = p + 1e-6 * (x + 20.0) * x * (x - 20.0)
        n = n + 1e-6 * (3.0 * x**2 - 400.0)
    elif perturbation == "reverse":
        p = p - 0.01 * x
        n = n - 0.01
    out = soft_branch(n)
    out.Q1 = p / n
    out.Q3 = mu / out.mn - 1.0
    out.Q6 = (mu - p / n) / out.mn - 1.0
    out.Q7 = out.Q6.copy()
    return out


@pytest.mark.parametrize("case,category", [
    ("three", "MULTIPLE_COEXISTENCE_CROSSINGS"),
    ("reverse", "NONPOSITIVE_DENSITY_JUMP"),
    ("gap", "NO_COMMON_MU_SUPPORT"),
])
def test_maxwell_preserves_topology_and_support_rejections(case, category):
    mu = np.array([950.0, 960.0, 1030.0, 1040.0])
    low = phase(mu)
    high = phase(mu, case)
    if case == "gap":
        low = phase(np.array([950.0, 960.0, 970.0]))
        high = phase(np.array([1010.0, 1020.0, 1030.0]))
    with pytest.raises(NativeExecutionError) as failure:
        construct_maxwell_transition(
            low, high, assume_local_charge_neutrality=True
        )
    assert failure.value.category == category
