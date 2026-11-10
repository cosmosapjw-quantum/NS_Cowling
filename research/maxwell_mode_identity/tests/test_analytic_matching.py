import numpy as np
import pytest

from research.maxwell_mode_identity.analytic_matching import (
    matching_determinant,
    basis_change_determinant,
    gauge_transform_system,
    gauge_transform_interface,
    interface_matching_form,
    implicit_interface_determinant,
    transfer_matching_determinant,
    smooth_match_radius_factor,
    dimensionless_interpolation_transfer,
    transfer_condition_number_diagnostic,
)


def test_complex_bilinear_matching_determinant_and_basis_covariance():
    pc = np.array([[1+2j], [0.3-0.4j]], dtype=complex)
    ps = np.array([[0.7-0.2j], [1.1+0.9j]], dtype=complex)
    d0 = matching_determinant(pc, ps)
    ac = np.array([[2.0+0.5j]])
    ass = np.array([[0.4-0.3j]])
    d1, factor = basis_change_determinant(pc, ps, ac, ass)
    assert d1 == pytest.approx(d0 * np.linalg.det(ac) * np.linalg.det(ass))
    assert factor == pytest.approx(np.linalg.det(ac) * np.linalg.det(ass))
    # The bilinear determinant must not silently use a Hermitian product.
    hermitian_det = np.linalg.det(np.column_stack([pc[:, 0], np.conjugate(ps[:, 0])]))
    assert not np.isclose(d0, hermitian_det)


def test_matching_determinant_shape_and_rank_contract():
    pc = np.eye(3, dtype=complex)[:, :1]
    ps = np.eye(3, dtype=complex)[:, 1:]
    assert matching_determinant(pc, ps) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        matching_determinant(np.ones((2, 1)), np.ones((3, 1)))
    with pytest.raises(ValueError):
        matching_determinant(np.ones((2, 2)), np.ones((2, 1)))


def test_gauge_transformed_matrix_has_correct_trace_cocycle():
    M = np.array([[0.2+0.1j, 0.3], [-0.4j, -0.7+0.2j]], dtype=complex)
    G = np.array([[2.0+0.2j, 0.1], [0.3j, 1.4-0.1j]], dtype=complex)
    Gp = np.array([[0.02j, 0.03], [-0.04, -0.01j]], dtype=complex)
    Mz = gauge_transform_system(M, G, Gp)
    expected = np.trace(M) - np.trace(np.linalg.solve(G, Gp))
    assert np.trace(Mz) == pytest.approx(expected)


def test_gauge_interface_and_matching_determinant_cocycle():
    Ti = np.array([[1.2, 0.4], [-0.2, 0.9]], dtype=complex)
    Gi = np.array([[2.0, 0.1j], [0.2, 1.3]], dtype=complex)
    Go = np.array([[1.1+0.2j, 0.0], [0.3, 0.8-0.1j]], dtype=complex)
    Tz = gauge_transform_interface(Ti, Gi, Go)
    assert np.linalg.det(Tz) == pytest.approx(np.linalg.det(Ti) * np.linalg.det(Gi) / np.linalg.det(Go))
    yi = np.array([0.3+0.2j, 1.0-0.4j])
    ys = np.array([1.1-0.1j, 0.4+0.5j])
    dy = matching_determinant((Ti @ yi).reshape(2, 1), ys.reshape(2, 1))
    dz = matching_determinant((Tz @ np.linalg.solve(Gi, yi)).reshape(2, 1), np.linalg.solve(Go, ys).reshape(2, 1))
    assert dy == pytest.approx(np.linalg.det(Go) * dz)


def test_smooth_match_radius_determinant_factor_matches_matrix_exponential():
    M = np.array([[0.03+0.02j, 0.01], [-0.04j, -0.02+0.01j]], dtype=complex)
    Y0 = np.array([[1.0+0.2j, 0.4-0.1j], [0.1+0.3j, 1.2-0.2j]])
    dr = 0.7
    U = _expm2(M * dr)
    Y1 = U @ Y0
    predicted = np.exp(np.trace(M) * dr) * np.linalg.det(Y0)
    assert np.linalg.det(Y1) == pytest.approx(predicted, rel=2e-13, abs=2e-14)
    assert smooth_match_radius_factor(np.trace(M), dr) == pytest.approx(np.exp(np.trace(M) * dr))


def _expm2(A):
    from scipy.linalg import expm
    return expm(A)


@pytest.mark.parametrize("theta", [0.0, 0.2, 0.5, 0.8, 1.0])
def test_implicit_interface_blocks_reproduce_la10_transfer(theta):
    kappa, A, nu_p = 3.4, 0.73, 0.019
    form = interface_matching_form(theta, kappa, A, nu_p)
    assert form.kind == "transfer"
    T = form.T
    Ci, Co = form.C_inner, form.C_outer
    np.testing.assert_allclose(Co @ T + Ci, 0.0, atol=2e-12, rtol=0.0)
    assert np.linalg.det(T) == pytest.approx(kappa, rel=3e-14)


def test_implicit_determinant_is_analytic_unit_times_transfer_determinant():
    theta, kappa, A, nu_p = 0.37, 2.6, 0.42+0.3j, 0.017
    ui = np.array([1.0+0.2j, 0.4-0.1j])
    vo = np.array([0.2+0.5j, 1.1-0.3j])
    form = interface_matching_form(theta, kappa, A, nu_p)
    e_imp = implicit_interface_determinant(form, ui, vo)
    e_tr = transfer_matching_determinant(form.T, ui, vo)
    assert e_imp == pytest.approx(A * nu_p * e_tr)


def test_zero_frequency_intermediate_transfer_pole_but_implicit_form_is_finite():
    theta, kappa, nu_p = 0.35, 2.7, 0.02
    f0 = interface_matching_form(theta, kappa, 0.0, nu_p)
    assert f0.kind == "implicit"
    assert np.all(np.isfinite(f0.C_inner)) and np.all(np.isfinite(f0.C_outer))
    assert np.linalg.matrix_rank(np.column_stack([f0.C_inner, f0.C_outer])) == 2
    assert np.linalg.det(f0.C_outer) == pytest.approx(0.0, abs=0.0)
    # At small nonzero A the graph transfer exists but its lower-left entry
    # diverges as 1/A; this is a chart singularity, not a determinant zero.
    f1 = interface_matching_form(theta, kappa, 1.0e-10, nu_p)
    assert f1.kind == "transfer"
    assert abs(f1.T[1, 0]) > 1.0e8
    assert np.linalg.det(f1.T) == pytest.approx(kappa, rel=1e-8)


def test_rapid_endpoint_has_regular_zero_frequency_extension():
    kappa, nu_p = 3.2, 0.015
    form = interface_matching_form(1.0, kappa, 0.0, nu_p)
    assert form.kind == "transfer"
    np.testing.assert_allclose(form.T, [[kappa, 0.0], [0.0, 1.0]], atol=0.0, rtol=0.0)
    assert np.linalg.det(form.T) == pytest.approx(kappa)


def test_zero_jump_is_identity_even_at_zero_frequency():
    for theta in [0.0, 0.3, 0.7, 1.0]:
        form = interface_matching_form(theta, 1.0, 0.0, 0.02)
        assert form.kind == "transfer"
        np.testing.assert_array_equal(form.T, np.eye(2))


def test_intermediate_zero_frequency_rank_drop_for_zero_jump_is_not_hidden():
    # kappa=1 is the no-interface case, so the helper takes the exact identity
    # branch before constructing the rank-deficient implicit blocks.
    form = interface_matching_form(0.4, 1.0, 0.0, 0.02)
    np.testing.assert_array_equal(form.T, np.eye(2))


def test_singular_gauge_and_invalid_interface_data_are_rejected():
    with pytest.raises(ValueError):
        gauge_transform_system(np.eye(2), np.array([[1, 0], [0, 0]]), np.eye(2))
    with pytest.raises(ValueError):
        interface_matching_form(0.4, 2.0, 0.0, 0.0)
    with pytest.raises(ValueError):
        matching_determinant(np.ones((2, 1)), np.ones((2, 0)))


def test_norm_based_complex_normalization_fails_cauchy_riemann():
    # f(z)=sqrt(1+|z|^2) is a typical Hermitian norm factor. Its Wirtinger
    # derivative dbar(f)=(x+i*y)/(2*f) is nonzero away from z=0.
    x, y = 0.7, -0.9
    f = np.sqrt(1.0 + x*x + y*y)
    dbar = complex(x, y) / (2.0*f)
    assert abs(dbar) > 0.0


def test_basis_multiplier_with_zero_inside_contour_changes_argument_count():
    # A basis factor nonzero on |lambda|=1 but vanishing at lambda=0.2
    # is not an admissible holomorphic unit on the enclosed disk.
    lam0, basis_zero = 0.1+0.0j, 0.2+0.0j
    d = lambda z: z-lam0
    unit = lambda z: z-basis_zero
    circle = np.exp(1j*np.linspace(0.0, 2*np.pi, 101))
    assert np.min(np.abs(unit(circle))) >= 0.8 - 1e-12
    assert d(basis_zero) != 0.0 and unit(basis_zero) == 0.0
    # Product gains an interior zero although the basis factor is nonzero on contour.
    assert d(basis_zero)*unit(basis_zero) == 0.0


def test_dimensionless_transfer_condition_scales_like_A_inverse_squared():
    theta, kappa, nu_p, L = 0.4, 2.8, 0.02, 1000.0
    x1, x2 = 2.0e-3, 1.0e-3
    A1, A2 = x1*nu_p*L, x2*nu_p*L
    f1 = interface_matching_form(theta, kappa, A1, nu_p)
    f2 = interface_matching_form(theta, kappa, A2, nu_p)
    T1 = dimensionless_interpolation_transfer(f1, L)
    T2 = dimensionless_interpolation_transfer(f2, L)
    np.testing.assert_allclose(T1, [[1+theta*(kappa-1), theta*(kappa-1)*x1],
                                    [(1-theta)*(kappa-1)/x1, 1+(1-theta)*(kappa-1)]])
    c1 = transfer_condition_number_diagnostic(f1, L)
    c2 = transfer_condition_number_diagnostic(f2, L)
    assert c2/c1 == pytest.approx(4.0, rel=2e-2)
    # This is a real-frequency conditioning diagnostic, not a holomorphic scalar.
    with pytest.raises(ValueError):
        transfer_condition_number_diagnostic(interface_matching_form(theta,kappa,1j*A2,nu_p),L)


def test_rapid_regularized_implicit_determinant_unit():
    theta, kappa, A, nu_p = 1.0, 2.3, 0.19+0.06j, 0.014
    ui = np.array([0.8+0.2j, -0.3+0.4j])
    vo = np.array([1.2-0.1j, 0.6+0.7j])
    form = interface_matching_form(theta, kappa, A, nu_p)
    e_imp = implicit_interface_determinant(form, ui, vo)
    e_transfer = transfer_matching_determinant(form.T, ui, vo)
    assert e_imp == pytest.approx(-nu_p*e_transfer)


def test_general_three_state_basis_cocycle():
    pc = np.array([[1+0.2j, 0.3], [0.1j, 1.4-0.1j], [0.4, -0.2j]])
    ps = np.array([[0.2-0.3j], [1.1+0.2j], [0.7]])
    ac = np.array([[1.2+0.1j, 0.2], [-0.1j, 0.8-0.2j]])
    ass = np.array([[1.7-0.4j]])
    d0 = matching_determinant(pc, ps)
    d1, factor = basis_change_determinant(pc, ps, ac, ass)
    assert d1 == pytest.approx(d0*factor)
    assert factor == pytest.approx(np.linalg.det(ac)*np.linalg.det(ass))


def test_implicit_matching_scalar_has_finite_A_to_zero_limit():
    theta, kappa, nu_p = 0.37, 2.6, 0.017
    ui = np.array([1.0+0.2j, 0.4-0.1j])
    vo = np.array([0.2+0.5j, 1.1-0.3j])
    e0 = implicit_interface_determinant(interface_matching_form(theta,kappa,0.0,nu_p),ui,vo)
    e1 = implicit_interface_determinant(interface_matching_form(theta,kappa,1e-10,nu_p),ui,vo)
    assert e1 == pytest.approx(e0, rel=3e-7, abs=1e-13)
