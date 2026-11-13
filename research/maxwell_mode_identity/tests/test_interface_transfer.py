import numpy as np
import pytest

from mode_identity.interface_transfer import (
    apply_interface_map,
    interface_map,
    interface_residuals,
)


@pytest.mark.parametrize("conversion", ["slow", "rapid"])
def test_transfer_preserves_pressure_perturbation_and_expected_kinematic_quantity(conversion):
    alpha, A, nu_prime = 1.7, 0.83, 0.21
    yin = np.array([2.4, -0.7])
    yout = apply_interface_map(conversion, alpha, A, nu_prime, yin)
    residuals = interface_residuals(conversion, alpha, A, nu_prime, yin, yout)
    assert residuals["q_scaling_residual"] < 1e-14
    if conversion == "slow":
        assert residuals["slow_radial_displacement_residual"] < 1e-14
        assert yout[0] == pytest.approx(yin[0])
    else:
        assert residuals["rapid_isobar_displacement_residual"] < 1e-14
        assert yout[1] == pytest.approx(yin[1])


@pytest.mark.parametrize("conversion", ["slow", "rapid"])
def test_map_determinant_and_zero_jump_limit(conversion):
    alpha, A, nu_prime = 2.3, 0.61, 0.17
    obj = interface_map(conversion, alpha, A, nu_prime)
    assert obj.determinant == pytest.approx(alpha, rel=1e-14)
    np.testing.assert_allclose(interface_map(conversion, 1.0, A, nu_prime).matrix,
                               np.eye(2), rtol=0.0, atol=0.0)


def test_slow_and_rapid_are_distinct_for_nonzero_jump():
    yin = np.array([1.2, 0.9])
    slow = apply_interface_map("slow", 1.8, 0.7, 0.2, yin)
    rapid = apply_interface_map("rapid", 1.8, 0.7, 0.2, yin)
    assert not np.allclose(slow, rapid)


@pytest.mark.parametrize("args", [(0.0, 1.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 0.0)])
def test_invalid_coefficients_rejected(args):
    with pytest.raises(ValueError):
        interface_map("slow", *args)


def test_invalid_conversion_and_state_rejected():
    with pytest.raises(ValueError):
        interface_map("intermediate", 1.0, 1.0, 1.0)
    with pytest.raises(ValueError):
        apply_interface_map("slow", 1.0, 1.0, 1.0, np.array([1.0]))
