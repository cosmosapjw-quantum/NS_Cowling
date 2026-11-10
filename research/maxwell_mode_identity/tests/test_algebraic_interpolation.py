import numpy as np
import pytest
from research.maxwell_mode_identity.finite_conversion_interpolation import algebraic_interpolation_map, interpolation_constraint_residual

def test_theta_zero_is_exact_slow_map():
    k,A,g=3.7,0.42,0.013
    got=algebraic_interpolation_map(0.0,k,A,g).matrix
    np.testing.assert_allclose(got,[[1,0],[(k-1)*g/A,k]],rtol=0,atol=2e-15)

def test_theta_one_is_exact_rapid_map():
    k,A,g=3.7,0.42,0.013
    got=algebraic_interpolation_map(1.0,k,A,g).matrix
    np.testing.assert_allclose(got,[[k,(k-1)*A/g],[0,1]],rtol=0,atol=2e-13)

@pytest.mark.parametrize("theta",[0.0,0.1,0.37,0.5,0.9,1.0])
def test_determinant_and_weighted_concomitant_jump(theta):
    k,A,g=4.2,0.71,0.021
    m=algebraic_interpolation_map(theta,k,A,g)
    assert m.determinant==pytest.approx(k,rel=2e-14,abs=2e-14)
    assert (1.0/k)*m.determinant==pytest.approx(1.0,rel=2e-14,abs=2e-14)

def test_interpolated_constraint_and_pressure_continuity():
    theta,k,A,g=0.43,2.8,0.63,0.017
    yi=np.array([1.3,-0.27]); yo=algebraic_interpolation_map(theta,k,A,g).matrix@yi
    residual=interpolation_constraint_residual(theta,k,A,g,yi,yo)
    assert residual["pressure_residual"]<1e-13
    assert residual["interpolated_kinematic_residual"]<1e-13

def test_zero_jump_is_identity_for_all_theta():
    for theta in np.linspace(0.0,1.0,11):
        np.testing.assert_allclose(algebraic_interpolation_map(float(theta),1.0,0.7,0.02).matrix,np.eye(2),rtol=0,atol=2e-15)

def test_theta_outside_unit_interval_is_rejected():
    with pytest.raises(ValueError): algebraic_interpolation_map(1.01,2.0,0.7,0.02)

def test_nonpositive_or_nonfinite_parameters_are_rejected():
    for p in [(-0.1,2.0,0.7,0.02),(0.5,0.0,0.7,0.02),(0.5,2.0,0.0,0.02),(0.5,2.0,0.7,float("nan"))]:
        with pytest.raises(ValueError): algebraic_interpolation_map(*p)
